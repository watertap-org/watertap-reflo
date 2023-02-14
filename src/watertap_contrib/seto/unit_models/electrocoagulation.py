###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

from copy import deepcopy

# Import Pyomo libraries
from pyomo.environ import (
    Set,
    Var,
    check_optimal_termination,
    Param,
    Suffix,
    log,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.constants import Constants
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import StrEnum

from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog

from watertap.core import InitializationMixin

__author__ = "Kurban Sitterley"

_log = idaeslog.getLogger(__name__)


class ElectrodeMaterial(StrEnum):
    aluminum = "aluminum"
    iron = "iron"


@declare_process_block_class("Electrocoagulation")
class ElectrocoagulationData(InitializationMixin, UnitModelBlockData):
    """
    Zero order electrocoagulation model
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False.""",
        ),
    )

    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
    **default** - False.""",
        ),
    )

    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )

    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
    and used when constructing these,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}""",
        ),
    )

    CONFIG.declare(
        "electrode_material",
        ConfigValue(
            default="aluminum",
            domain=In(ElectrodeMaterial),
            description="Electrode material",
        ),
    )

    def build(self):
        super().build()

        if "TDS" not in self.config.property_package.component_list:
            raise ConfigurationError("TDS must be in feed stream")

        if self.config.electrode_material == ElectrodeMaterial.aluminum:
            self.mw_electrode_material = Param(
                initialize=0.027,
                units=pyunits.kg / pyunits.mol,
                doc="Molecular weight of electrode material",
            )
            self.valence_electrode_material = Param(
                initialize=3,
                units=pyunits.dimensionless,
                doc="Number of valence electrons of electrode material",
            )
            self.density_electrode_material = Param(
                initialize=2710,
                units=pyunits.kg / pyunits.m**3,
                doc="Density of electrode material",
            )
            self.metal_dose_to_toc_ratio = Param(
                initialize=1,
                units=pyunits.kg / pyunits.kg,
                doc="Coagulant dose to inlet TOC ratio",
            )

        if self.config.electrode_material == ElectrodeMaterial.iron:
            self.mw_electrode_material = Param(
                initialize=0.056,
                units=pyunits.kg / pyunits.mol,
                doc="Molecular weight of electrode material",
            )
            self.valence_electrode_material = Param(
                initialize=2,
                units=pyunits.dimensionless,
                doc="Number of valence electrons of electrode material",
            )
            self.density_electrode_material = Param(
                initialize=7860,
                units=pyunits.kg / pyunits.m**3,
                doc="Density of electrode material",
            )
            self.metal__dose_to_toc_ratio = Param(
                initialize=2,
                units=pyunits.kg / pyunits.kg,
                doc="Coagulant dose to inlet TOC ratio",
            )

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True  # inlet block is an inlet
        self.properties_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of inlet", **tmp_dict
        )

        # Add outlet and waste block
        tmp_dict["defined_state"] = False
        self.properties_out = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of outlet",
            **tmp_dict,
        )

        self.properties_waste = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of waste", **tmp_dict
        )

        self.add_port(name="inlet", block=self.properties_in)
        self.add_port(name="outlet", block=self.properties_out)
        self.add_port(name="waste", block=self.properties_waste)

        self.number_cells_system = Param(
            initialize=3,
            mutable=True,
            doc="Number total cells for system - electrochemical cell > flotation cell > sedimentation cell",
        )

        self.number_redundant_cells = Param(
            initialize=2,
            mutable=True,
            doc="Number redundant electrochemical cells",
        )

        self.current_per_reactor = Param(
            initialize=3000,
            units=pyunits.ampere,
            mutable=True,
            doc="Current required per reactor",
        )

        self.tds_to_cond_conversion = Param(
            initialize=5e3,
            mutable=True,
            units=(pyunits.mg * pyunits.m) / (pyunits.liter * pyunits.S),
            doc="Conersion factor for mg/L TDS to S/m",
        )

        self.removal_efficiency = Param(
            self.config.property_package.component_set,
            initialize=0.7,
            mutable=True,
            doc="Removal efficiency",
        )

        self.vol_recovery = Param(
            initialize=0.99, mutable=True, doc="Volumetric recovery"
        )

        self.electrode_width = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m,
            doc="Electrode width",
        )

        self.electrode_height = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m,
            doc="Electrode height",
        )

        self.electrode_thick = Var(
            initialize=0.001,
            bounds=(0, 0.1),
            units=pyunits.m,
            doc="Electrode thickness",
        )

        self.electrode_mass = Var(
            initialize=10,
            bounds=(0, None),
            units=pyunits.kg,
            doc="Electrode mass",
        )

        self.electrode_area_total = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Total electrode area",
        )

        self.electrode_area_per = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**2,
            doc="Electrode area",
        )

        self.electrode_volume_per = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Electrode volume",
        )

        self.electrode_gap = Var(
            initialize=0.005,
            bounds=(0.001, 0.2),
            units=pyunits.m,
            doc="Electrode gap",
        )

        self.electrolysis_time = Var(
            initialize=30,
            bounds=(2, 200),
            units=pyunits.minute,
            doc="Electrolysis time",
        )

        self.number_electrode_pairs = Var(
            initialize=5,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of electrode pairs",
        )

        self.number_cells = Var(
            initialize=5,
            bounds=(0, None),
            units=pyunits.dimensionless,
            doc="Number of cells",
        )

        self.applied_current = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.ampere,
            doc="Applied current",
        )

        self.current_efficiency = Var(
            initialize=1,
            bounds=(0.9, 2.5),
            units=pyunits.kg / pyunits.kg,
            doc="Current efficiency",
        )

        self.cell_voltage = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.volt,
            doc="Cell voltage",
        )

        self.potential_balance = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.volt,
            doc="Potential balance for voltage", #???
        )

        self.reactor_volume = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.m**3,
            doc="Reactor volume total (electrochemical + flotation + sedimentation)",
        )

        self.metal_loading = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.kg / pyunits.liter,
            doc="Metal loading",
        )

        # self.metal_loading_rate = Var(
        #     initialize=1,
        #     units=pyunits.kg,
        #     doc="Metal loading",
        # )

        self.ohmic_resistance = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.ohm,
            doc="Ohmic resistance of solution",
        )

        self.charge_loading_rate = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.coulomb / pyunits.liter,
            doc="Charge loading rate",
        )

        self.current_density = Var(
            initialize=1,
            bounds=(1, 2000),
            units=pyunits.ampere / pyunits.m**2,
            doc="Current density",
        )

        @self.Constraint(doc="Charge loading rate equation")
        def eq_charge_loading_rate(b):
            return b.charge_loading_rate == (
                b.applied_current
                * pyunits.convert(b.electrolysis_time, to_units=pyunits.second)
            ) / pyunits.convert(b.reactor_volume, to_units=pyunits.liter)

        @self.Constraint(doc="Metal loading rate equation")
        def eq_metal_loading_rate(b):
            return b.metal_loading == (
                b.current_efficiency * b.charge_loading_rate * b.mw_electrode_material
            ) / (b.valence_electrode_material * Constants.faraday_constant)

        @self.Constraint(doc="Total current required")
        def eq_applied_current(b):
            flow_in = pyunits.convert(
                b.properties_in[0].flow_vol, to_units=pyunits.liter / pyunits.second
            )
            return b.applied_current == (
                flow_in
                * b.metal_loading
                * b.valence_electrode_material
                * Constants.faraday_constant
            ) / (b.current_efficiency * b.mw_electrode_material)

        @self.Constraint(doc="Total electrode area required")
        def eq_electrode_area_total(b):
            return b.electrode_area_total == b.applied_current / b.current_density

        # @self.Constraint(doc="Number of cells")
        # def eq_number_cells(b):
        #     return b.number_cells == b.applied_current / b.current_per_reactor

        @self.Constraint(doc="Cell voltage")
        def eq_cell_voltage(b):
            return b.cell_voltage == b.potential_balance + b.applied_current * b.ohmic_resistance

        @self.Constraint(doc="Area per electrode")
        def eq_electrode_area_per(b):
            return b.electrode_area_per == b.electrode_area_total / (
                b.number_electrode_pairs * 2
            )

        @self.Constraint(doc="Electrode width")
        def eq_electrode_width(b):
            return b.electrode_width == (2 * b.electrode_area_per) ** 0.5

        @self.Constraint(doc="Electrode height")
        def eq_electrode_height(b):
            return b.electrode_height == b.electrode_area_per / b.electrode_width

        @self.Constraint(doc="Electrode volume")
        def eq_electrode_volume_per(b):
            return (
                b.electrode_volume_per
                == b.electrode_width * b.electrode_height * b.electrode_thick
            )

        @self.Constraint(doc="Total reactor volume")
        def eq_reactor_volume(b):
            flow_vol = b.properties_in[0].flow_vol
            return b.reactor_volume == pyunits.convert(
                flow_vol * b.electrolysis_time,
                to_units=pyunits.m**3,
            ) / b.number_cells

        @self.Expression(doc="Conductivity")
        def conductivity(b):
            tds = pyunits.convert(
                b.properties_in[0].conc_mass_comp["TDS"],
                to_units=pyunits.mg / pyunits.L,
            )
            return tds / b.tds_to_cond_conversion

        @self.Constraint(doc="Ohmic resistance")
        def eq_ohmic_resistance(b):
            return b.ohmic_resistance == b.electrode_gap / (
                b.conductivity * b.electrode_area_per * b.number_cells
            )

        @self.Constraint(doc="Electrode mass")
        def eq_electrode_mass(b):
            return (
                b.electrode_mass
                == b.electrode_volume_per * b.density_electrode_material
            )

        @self.Constraint(doc="Effluent flow")
        def eq_effluent_flow(b):
            prop_in = b.properties_in[0]
            prop_out = b.properties_out[0]
            return prop_out.flow_vol == prop_in.flow_vol * b.vol_recovery
 
        @self.Constraint(
            self.config.property_package.component_set, doc="Component removal"
        )
        def eq_component_removal(b, j):
            prop_in = b.properties_in[0]
            prop_waste = b.properties_waste[0]
            return (
                b.removal_efficiency[j] * prop_in.flow_mass_comp[j]
                == prop_waste.flow_mass_comp[j]
            )

        @self.Constraint(doc="Flow balance")
        def eq_flow_balance(b):
            prop_in = b.properties_in[0]
            prop_out = b.properties_out[0]
            prop_waste = b.properties_waste[0]
            return prop_in.flow_vol == prop_out.flow_vol + prop_waste.flow_vol

        @self.Constraint(
            self.config.property_package.component_set, doc="Component mass balance"
        )
        def eq_mass_balance(b, j):
            prop_in = b.properties_in[0]
            prop_out = b.properties_out[0]
            prop_waste = b.properties_waste[0]
            return prop_in.flow_mass_comp[j] == prop_out.flow_mass_comp[j] + prop_waste.flow_mass_comp[j]

    def initialize_build(
        blk,
        state_args=None,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for initialization routines

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        flags = blk.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info("Initialization Step 1a Complete.")
        # ---------------------------------------------------------------------
        # Initialize other state blocks
        # Set state_args from inlet state
        if state_args is None:
            blk.state_args = state_args = {}
            state_dict = blk.properties_in[
                blk.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        blk.flags = flags
        state_args_out = deepcopy(state_args)

        for j in blk.properties_out.component_list:
            if j == "H2O":
                state_args_out["flow_vol"] = state_args["flow_vol"] * blk.vol_recovery
            else:
                state_args_out["conc_mass_comp"][j] = state_args["conc_mass_comp"][j] * (1 - blk.removal_efficiency[j])

        blk.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        init_log.info("Initialization Step 1b Complete.")

        state_args_waste = deepcopy(state_args)
        for j in blk.properties_waste.component_list:
            if j == "H2O":
                state_args_waste["flow_vol"] = state_args["flow_vol"] * (1 - blk.vol_recovery)
            else:
                state_args_waste["conc_mass_comp"][j] = state_args["conc_mass_comp"][j] * blk.removal_efficiency[j]
        
        blk.properties_waste.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_waste,
        )

        init_log.info("Initialization Step 1c Complete.")

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        # blk.properties_in.release_state(flags, outlvl=outlvl)
        # init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {blk.name} failed to initialize")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        iscale.set_scaling_factor(self.electrode_width, 10)

        iscale.set_scaling_factor(self.electrode_height, 10)

        iscale.set_scaling_factor(self.electrode_thick, 1e3)

        iscale.set_scaling_factor(self.electrode_mass, 1)

        iscale.set_scaling_factor(self.electrode_area_total, 1e-2)

        iscale.set_scaling_factor(self.electrode_area_per, 10)

        iscale.set_scaling_factor(self.electrode_volume_per, 1)

        iscale.set_scaling_factor(self.electrode_gap, 10)

        iscale.set_scaling_factor(self.electrolysis_time, 0.1)

        iscale.set_scaling_factor(self.number_electrode_pairs, 0.1)

        iscale.set_scaling_factor(self.number_cells, 1)

        iscale.set_scaling_factor(self.applied_current, 1e-4)

        iscale.set_scaling_factor(self.current_efficiency, 1)

        iscale.set_scaling_factor(self.reactor_volume, 0.1)

        iscale.set_scaling_factor(self.metal_loading, 1e6)

        iscale.set_scaling_factor(self.ohmic_resistance, 1e5)

        iscale.set_scaling_factor(self.charge_loading_rate, 1e-2)

        iscale.set_scaling_factor(self.current_density, 1e-2)


        # transforming constraints
        sf = iscale.get_scaling_factor(self.metal_loading)
        iscale.constraint_scaling_transform(self.eq_metal_loading_rate, sf)

        sf = iscale.get_scaling_factor(self.charge_loading_rate)
        iscale.constraint_scaling_transform(self.eq_charge_loading_rate, sf)

        sf = iscale.get_scaling_factor(self.ohmic_resistance)
        iscale.constraint_scaling_transform(self.eq_ohmic_resistance, sf)

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.inlet,
                "Liquid Outlet": self.outlet,
                "Waste Outlet": self.waste,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):

        # TODO
        var_dict = {}

        return {"vars": var_dict}
