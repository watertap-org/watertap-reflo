from pyomo.environ import (
    Var,
    Constraint,
    Param,
    value,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

import idaes.logger as idaeslog
from idaes.core import UnitModelBlockData, useDefault, declare_process_block_class
from idaes.core.util.scaling import calculate_scaling_factors, set_scaling_factor
from idaes.core.util.config import is_physical_parameter_block

from watertap.costing.util import register_costing_parameter_block
from watertap.core import InitializationMixin
from watertap_contrib.reflo.core import SolarEnergyBaseData
from watertap.core.solvers import get_solver

from watertap_contrib.reflo.costing.util import (
    make_capital_cost_var,
    make_variable_operating_cost_var,
    make_fixed_operating_cost_var,
)

solver = get_solver()

"""
Set of dummy treatment and energy generation units used to test the costing package.
"""

__author__ = "Kurban Sitterley"

############################################################################
############################################################################


@declare_process_block_class("DummyTreatmentUnit")
class DummyTreatmentUnitData(InitializationMixin, UnitModelBlockData):
    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(default=False, domain=In([False])),
    )

    CONFIG.declare(
        "has_holdup",
        ConfigValue(default=False, domain=In([False])),
    )

    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
        ),
    )

    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
        ),
    )

    def build(self):
        super().build()

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = False

        self.properties = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Unit properties", **tmp_dict
        )

        self.chemical_dose = Param(
            initialize=2.2e-2,
            mutable=True,
            units=pyunits.kg / pyunits.m**3,
            doc="Chemical dose",
        )

        self.design_var_a = Var(
            initialize=42,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Test treatment unit design variable",
        )

        self.design_var_b = Var(
            initialize=1.23e-4,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Test treatment unit design variable",
        )

        self.capital_var = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Test treatment unit capital variable",
        )

        self.fixed_operating_var = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Test treatment unit fixed operating variable",
        )

        self.variable_operating_var = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Test treatment unit variable operating variable",
        )

        self.electricity_consumption = Var(
            initialize=1e4,
            units=pyunits.kilowatt,
            bounds=(0, None),
            doc="Constant energy consumption",
        )

        self.heat_consumption = Var(
            initialize=2e4,
            units=pyunits.kilowatt,
            bounds=(0, None),
            doc="Constant heat consumption",
        )

        dimensionless_flow_vol = pyunits.convert(
            self.properties[0].flow_vol_phase["Liq"] * pyunits.m**-3 * pyunits.s,
            to_units=pyunits.dimensionless,
        )
        dimensionless_flow_mass = pyunits.convert(
            self.properties[0].flow_mass_phase_comp["Liq", "TDS"]
            * pyunits.s
            * pyunits.kg**-1,
            to_units=pyunits.dimensionless,
        )
        dimensionless_conc = pyunits.convert(
            self.properties[0].conc_mass_phase_comp["Liq", "TDS"]
            * pyunits.m**3
            * pyunits.kg**-1,
            to_units=pyunits.dimensionless,
        )

        @self.Constraint(doc="Capital variable calculation")
        def eq_capital_var(b):
            return b.capital_var == b.design_var_a * dimensionless_flow_vol

        @self.Constraint(doc="Fixed operating variable calculation")
        def eq_fixed_operating_var(b):
            return b.fixed_operating_var == b.design_var_b * dimensionless_conc

        @self.Constraint(doc="Variable operating variable calculation")
        def eq_variable_operating_var(b):
            return (
                b.variable_operating_var
                == (b.design_var_a * b.design_var_b) * dimensionless_flow_mass
            )

    def initialize_build(self):
        solve_log = idaeslog.getSolveLogger(self.name, tag="unit")

        opt = get_solver()

        flags = self.properties.initialize(hold_state=True)
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
            if not check_optimal_termination(res):
                res = opt.solve(self, tee=slc.tee)  # just try again!

        self.properties.release_state(flags)

    def calculate_scaling_factors(self):

        set_scaling_factor(self.design_var_a, 1 / value(self.design_var_a))
        set_scaling_factor(self.design_var_b, 1 / value(self.design_var_b))
        set_scaling_factor(self.capital_var, 1 / value(self.capital_var))
        set_scaling_factor(
            self.fixed_operating_var, 1 / value(self.fixed_operating_var)
        )
        set_scaling_factor(
            self.variable_operating_var, 1 / value(self.variable_operating_var)
        )
        set_scaling_factor(
            self.electricity_consumption, 1 / value(self.electricity_consumption)
        )
        set_scaling_factor(self.heat_consumption, 1 / value(self.heat_consumption))

    @property
    def default_costing_method(self):
        return cost_dummy_treatment_unit


def build_dummy_treatment_unit_param_block(blk):

    blk.capital_cost_param = Var(
        initialize=1e4,
        units=pyunits.USD_2000,
        bounds=(0, None),
        doc="Capital cost param",
    )

    blk.fixed_operating_cost_param = Var(
        initialize=3,
        units=pyunits.USD_2011 / pyunits.m**3,
        bounds=(0, None),
        doc="Operating cost param",
    )

    blk.variable_operating_cost_param = Var(
        initialize=4.2e-1,
        units=pyunits.USD_2020 / pyunits.m**3,
        bounds=(0, None),
        doc="Operating cost param",
    )


def build_chem_cost_param_block(blk):

    blk.cost = Param(
        mutable=True,
        initialize=0.68,
        units=pyunits.USD_2016 / pyunits.kg,
        doc="Chemical unit cost",
    )

    blk.purity = Param(
        mutable=True,
        initialize=0.56,
        units=pyunits.dimensionless,
        doc="Chemical purity",
    )

    blk.parent_block().register_flow_type("chemical", blk.cost / blk.purity)


@register_costing_parameter_block(
    build_rule=build_chem_cost_param_block,
    parameter_block_name="test_chemical",
)
@register_costing_parameter_block(
    build_rule=build_dummy_treatment_unit_param_block,
    parameter_block_name="dummy_treatment_unit",
)
def cost_dummy_treatment_unit(blk):

    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)
    make_variable_operating_cost_var(blk)

    blk.costing_package.add_cost_factor(blk, "TIC")

    blk.capital_cost_constraint = Constraint(
        expr=blk.capital_cost
        == pyunits.convert(
            blk.costing_package.dummy_treatment_unit.capital_cost_param
            * blk.unit_model.capital_var,
            to_units=blk.costing_package.base_currency,
        )
    )

    blk.fixed_operating_cost_constraint = Constraint(
        expr=blk.fixed_operating_cost
        == pyunits.convert(
            blk.costing_package.dummy_treatment_unit.fixed_operating_cost_param
            * blk.unit_model.properties[0].flow_vol_phase["Liq"]
            * blk.unit_model.fixed_operating_var,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    blk.variable_operating_cost_constraint = Constraint(
        expr=blk.variable_operating_cost
        == pyunits.convert(
            blk.costing_package.dummy_treatment_unit.variable_operating_cost_param
            * blk.unit_model.properties[0].flow_vol_phase["Liq"]
            * blk.unit_model.variable_operating_var,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    blk.costing_package.cost_flow(
        blk.unit_model.electricity_consumption,
        "electricity",
    )

    blk.costing_package.cost_flow(
        blk.unit_model.heat_consumption,
        "heat",
    )

    blk.costing_package.cost_flow(
        pyunits.convert(
            blk.unit_model.chemical_dose
            * blk.unit_model.properties[0].flow_vol_phase["Liq"],
            to_units=pyunits.kg / blk.costing_package.base_period,
        ),
        "chemical",
    )


############################################################################
############################################################################


@declare_process_block_class("DummyTreatmentNoHeatUnit")
class DummyTreatmentNoHeatUnitData(DummyTreatmentUnitData):
    CONFIG = DummyTreatmentUnitData.CONFIG()

    def build(self):
        super().build()
        self.del_component(self.heat_consumption)
        self.del_component(self.fixed_operating_var)
        self.del_component(self.variable_operating_var)

    def calculate_scaling_factors(self):

        set_scaling_factor(self.design_var_a, 1 / value(self.design_var_a))
        set_scaling_factor(self.design_var_b, 1 / value(self.design_var_b))
        set_scaling_factor(self.capital_var, 1 / value(self.capital_var))
        set_scaling_factor(
            self.electricity_consumption, 1 / value(self.electricity_consumption)
        )

    @property
    def default_costing_method(self):
        return cost_dummy_treatment_no_heat_unit


@register_costing_parameter_block(
    build_rule=build_dummy_treatment_unit_param_block,
    parameter_block_name="dummy_treatment_no_heat_unit",
)
def cost_dummy_treatment_no_heat_unit(blk):

    make_capital_cost_var(blk)

    blk.costing_package.add_cost_factor(blk, None)

    blk.capital_cost_constraint = Constraint(
        expr=blk.capital_cost
        == pyunits.convert(
            blk.costing_package.dummy_treatment_no_heat_unit.capital_cost_param
            * blk.unit_model.capital_var,
            to_units=blk.costing_package.base_currency,
        )
    )

    blk.costing_package.cost_flow(
        blk.unit_model.electricity_consumption,
        "electricity",
    )


############################################################################
############################################################################


@declare_process_block_class("DummyElectricityUnit")
class DummyElectricityUnitData(SolarEnergyBaseData):
    """
    Test unit for electricity generation.
    Generates zero heat.
    """

    CONFIG = SolarEnergyBaseData.CONFIG()
    CONFIG.solar_model_type = "physical"

    def build(self):
        super().build()

        self.heat.fix(0)

    @property
    def default_costing_method(self):
        return cost_dummy_electricity_unit


def build_dummy_electricity_unit_param_block(blk):

    blk.capital_per_watt = Var(
        initialize=0.3,
        units=pyunits.USD_2019 / pyunits.watt,
        bounds=(0, None),
        doc="Capital cost per watt",
    )

    blk.fixed_operating_per_watt = Var(
        initialize=0.042,
        units=pyunits.USD_2019 / (pyunits.watt * pyunits.year),
        bounds=(0, None),
        doc="Operating cost per watt",
    )


@register_costing_parameter_block(
    build_rule=build_dummy_electricity_unit_param_block,
    parameter_block_name="dummy_electricity_unit",
)
def cost_dummy_electricity_unit(blk):

    blk.costing_package.has_electricity_generation = True

    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.costing_package.add_cost_factor(blk, None)
    blk.capital_cost_constraint = Constraint(
        expr=blk.capital_cost
        == pyunits.convert(
            blk.costing_package.dummy_electricity_unit.capital_per_watt
            * blk.unit_model.electricity,
            to_units=blk.costing_package.base_currency,
        )
    )

    blk.fixed_operating_cost_constraint = Constraint(
        expr=blk.fixed_operating_cost
        == pyunits.convert(
            blk.costing_package.dummy_electricity_unit.fixed_operating_per_watt
            * blk.unit_model.electricity,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    # Generating is negative by convention
    blk.costing_package.cost_flow(-1 * blk.unit_model.electricity, "electricity")


############################################################################
############################################################################


@declare_process_block_class("DummyHeatUnit")
class DummyHeatUnitData(SolarEnergyBaseData):
    CONFIG = SolarEnergyBaseData.CONFIG()
    CONFIG.solar_model_type = "physical"

    def build(self):
        super().build()

        self.electricity.fix(20)

    @property
    def default_costing_method(self):
        return cost_dummy_heat_unit


def build_dummy_heat_unit_param_block(blk):

    blk.capital_per_watt = Var(
        initialize=0.15,
        units=pyunits.USD_2019 / pyunits.watt,
        bounds=(0, None),
        doc="Cost per watt",
    )

    blk.fixed_operating_per_watt = Var(
        initialize=0.019,
        units=pyunits.USD_2019 / (pyunits.watt * pyunits.year),
        bounds=(0, None),
        doc="Cost per watt",
    )


@register_costing_parameter_block(
    build_rule=build_dummy_heat_unit_param_block,
    parameter_block_name="dummy_heat_unit",
)
def cost_dummy_heat_unit(blk):

    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    blk.costing_package.add_cost_factor(blk, None)
    blk.capital_cost_constraint = Constraint(
        expr=blk.capital_cost
        == pyunits.convert(
            blk.costing_package.dummy_heat_unit.capital_per_watt * blk.unit_model.heat,
            to_units=blk.costing_package.base_currency,
        )
    )

    blk.fixed_operating_cost_constraint = Constraint(
        expr=blk.fixed_operating_cost
        == pyunits.convert(
            blk.costing_package.dummy_heat_unit.fixed_operating_per_watt
            * blk.unit_model.heat,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )

    # Generating is negative by convention
    blk.costing_package.cost_flow(-1 * blk.unit_model.heat, "heat")
    # Heat generatig units could require energy
    blk.costing_package.cost_flow(blk.unit_model.electricity, "electricity")
