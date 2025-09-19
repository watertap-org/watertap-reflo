#################################################################################
# WaterTAP Copyright (c) 2020-2024, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################


from pyomo.environ import (
    Constraint,
    NonNegativeReals,
    Var,
    Param,
    log,
    exp,
    value,
    units as pyunits,
)
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import StrEnum

from watertap.costing.util import register_costing_parameter_block
from ..util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
    make_variable_operating_cost_var,
)

__author__ = "Kurban Sitterley"

blm_costing_params_dict = {
    2500: {
        "logging_testing": {"intercept": 251.48, "slope": 4.5193},
        "drilling": {"intercept": 338.05, "slope": 17.611},
        "piping": {"base": 116.65, "exponent": 0.3786},
        "casing": {"intercept": 384.36, "slope": 20.545},
        "grouting": {"base": 145.01, "exp_coeff": 0.0708},
    },
    5000: {
        "logging_testing": {"intercept": 306.27, "slope": 5.5412},
        "drilling": {"intercept": 616.91, "slope": 31.783},
        "piping": {"base": 213.65, "exponent": 0.3748},
        "casing": {"intercept": 698.25, "slope": 37.591},
        "grouting": {"base": 268.83, "exp_coeff": 0.0708},
    },
    7500: {
        "logging_testing": {"intercept": 361.17, "slope": 6.5739},
        "drilling": {"intercept": 889.24, "slope": 46.01},
        "piping": {"base": 309.88, "exponent": 0.3733},
        "casing": {"intercept": 956.01, "slope": 51.407},
        "grouting": {"base": 368.65, "exp_coeff": 0.0709},
    },
    10000: {
        "logging_testing": {"intercept": 415.55, "slope": 7.5694},
        "drilling": {"intercept": 1166.5, "slope": 60.407},
        "piping": {"base": 407.38, "exponent": 0.3741},
        "casing": {"intercept": 1228.1, "slope": 65.52},
        "grouting": {"base": 465.73, "exp_coeff": 0.0708},
    },
}


class DeepWellInjectionCostMethod(StrEnum):
    as_capex = "as_capex"
    as_opex = "as_opex"
    blm = "blm"


def cost_deep_well_injection(blk, cost_method=DeepWellInjectionCostMethod.blm):

    if cost_method == DeepWellInjectionCostMethod.blm:
        cost_deep_well_injection_blm(blk)
    elif cost_method == DeepWellInjectionCostMethod.as_capex:
        cost_deep_well_injection_as_capex(blk)
    elif cost_method == DeepWellInjectionCostMethod.as_opex:
        cost_deep_well_injection_as_opex(blk)
    else:
        raise ConfigurationError(
            f"{blk.unit_model.name} received invalid argument for cost_method:"
            f" {cost_method}. Argument must be a member of the DeepWellInjectionCostMethod Enum."
        )


def build_deep_well_injection_cost_as_opex_param_block(blk):
    costing = blk.parent_block()
    blk.dwi_lcow = Var(
        initialize=0.0587,
        bounds=(0, None),
        units=costing.base_currency / pyunits.m**3,
        doc="Assumed LCOW of deep well injection",
    )


def build_deep_well_injection_cost_as_capex_param_block(blk):
    costing = blk.parent_block()
    blk.dwi_lcow = Var(
        initialize=0.0587,
        bounds=(0, None),
        units=costing.base_currency / pyunits.m**3,
        doc="Assumed LCOW of deep well injection",
    )


@register_costing_parameter_block(
    build_rule=build_deep_well_injection_cost_as_capex_param_block,
    parameter_block_name="deep_well_injection",
)
def cost_deep_well_injection_as_capex(blk):
    """
    Provied LCOW is used to calculate the equivalent capital costs.
    All of the LCOW is attributable to CAPEX.
    """
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)

    blk.base_period_flow = pyunits.convert(
        blk.unit_model.properties[0].flow_vol_phase["Liq"]
        * blk.costing_package.utilization_factor,
        to_units=pyunits.m**3 / blk.costing_package.base_period,
    )

    blk.capital_cost_constraint = Constraint(
        expr=blk.capital_cost
        == pyunits.convert(
            (blk.costing_package.deep_well_injection.dwi_lcow * blk.base_period_flow)
            / blk.costing_package.capital_recovery_factor,
            to_units=blk.costing_package.base_currency,
        )
    )

    blk.fixed_operating_cost_constraint = Constraint(
        expr=blk.fixed_operating_cost
        == -blk.costing_package.maintenance_labor_chemical_factor * blk.capital_cost
    )


@register_costing_parameter_block(
    build_rule=build_deep_well_injection_cost_as_opex_param_block,
    parameter_block_name="deep_well_injection",
)
def cost_deep_well_injection_as_opex(blk):
    """
    Provied LCOW is used to calculate the equivalent operating costs.
    All of the LCOW is attributable to OPEX.
    """
    make_capital_cost_var(blk)
    make_variable_operating_cost_var(blk)
    blk.costing_package.add_cost_factor(blk, None)

    blk.capital_cost_constraint = Constraint(
        expr=blk.capital_cost == 0 * blk.costing_package.base_currency
    )

    blk.base_period_flow = pyunits.convert(
        blk.unit_model.properties[0].flow_vol_phase["Liq"]
        * blk.costing_package.utilization_factor,
        to_units=pyunits.m**3 / blk.costing_package.base_period,
    )

    blk.variable_operating_cost_constraint = Constraint(
        expr=blk.variable_operating_cost
        == pyunits.convert(
            blk.costing_package.deep_well_injection.dwi_lcow * blk.base_period_flow,
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )


def build_deep_well_injection_cost_blm_param_block(blk):
    """
    Costing DWI unit model from BLM reference

    Membrane Concentrate Disposal: Practices and Regulation
    Michael C. Mickley, P.E., Ph.D.
    Prepared for US Dept. of Interior, Bureau of Land Management
    Agreement No. 98-FC-81-0054
    April 2006

    Data from plots in this document was extracted and fit using Excel.

    """

    blk.packing_capital_cost_slope = Param(
        initialize=21.153,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Packing capital cost equation - base",
    )

    blk.packing_capital_cost_intercept = Param(
        initialize=37.31,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Packing capital cost equation - intercept",
    )

    blk.monitoring_well_capital_cost_base = Param(
        initialize=31.337,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Monitoring well capital cost equation - base",
    )

    blk.monitoring_well_capital_cost_exponent = Param(
        initialize=0.3742,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Monitoring well capital cost equation - exponent",
    )

    blk.mobilization_capital_cost_base = Param(
        initialize=31.337,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Mobilization capital cost equation - base",
    )

    blk.mobilization_capital_cost_exponent = Param(
        initialize=0.3742,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Mobilization capital cost equation - exponent",
    )

    # All default values are for 2,500 ft depth

    blk.drilling_capital_cost_intercept = Param(
        initialize=338.05,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Drilling capital cost equation - intercept",
    )

    blk.drilling_capital_cost_slope = Param(
        initialize=17.611,
        mutable=True,
        units=pyunits.kUSD_2001 * pyunits.inches**-1,
        doc="Drilling capital cost equation - slope",
    )

    blk.piping_capital_cost_base = Param(
        initialize=116.65,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Injection piping capital cost equation - base",
    )

    blk.piping_capital_cost_exponent = Param(
        initialize=0.3786,
        mutable=True,
        units=pyunits.dimensionless,
        doc="Injection piping capital cost equation - exponent",
    )

    blk.casing_capital_cost_intercept = Param(
        initialize=20.545,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Casing capital cost equation - intercept",
    )

    blk.casing_capital_cost_slope = Param(
        initialize=384.36,
        mutable=True,
        units=pyunits.kUSD_2001 * pyunits.inches**-1,
        doc="Casing capital cost equation - slope",
    )

    blk.grouting_capital_cost_base = Param(
        initialize=145.01,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Grouting capital cost equation - base",
    )

    blk.grouting_capital_cost_exp_coeff = Param(
        initialize=0.0708,
        mutable=True,
        units=pyunits.inches**-1,
        doc="Grouting capital cost equation - exponential coefficient",
    )

    blk.logging_testing_capital_cost_intercept = Param(
        initialize=251.48,
        mutable=True,
        units=pyunits.kUSD_2001,
        doc="Logging + testing capital cost equation - intercept",
    )

    blk.logging_testing_capital_cost_slope = Param(
        initialize=4.5193,
        mutable=True,
        units=pyunits.kUSD_2001 * pyunits.inches**-1,
        doc="Logging + testing capital cost equation - slope",
    )


@register_costing_parameter_block(
    build_rule=build_deep_well_injection_cost_blm_param_block,
    parameter_block_name="deep_well_injection",
)
def cost_deep_well_injection_blm(blk):
    """
    Capital and operating costs for deep well injection
    based on BLM reference
    """
    dwi_params = blk.costing_package.deep_well_injection
    dwi = blk.unit_model
    pipe_diam_dimensionless = pyunits.convert(
        dwi.pipe_diameter * pyunits.inches**-1, to_units=pyunits.dimensionless
    )
    injection_well_depth_dimensionless = pyunits.convert(
        dwi.injection_well_depth * pyunits.ft**-1, to_units=pyunits.dimensionless
    )
    monitoring_well_depth_dimensionless = pyunits.convert(
        dwi.monitoring_well_depth * pyunits.ft**-1, to_units=pyunits.dimensionless
    )
    make_capital_cost_var(blk)

    costing_params_depth_dict = blm_costing_params_dict[
        int(value(injection_well_depth_dimensionless))
    ]
    # change value of costing parameters based on
    # user provided injection_well_depth
    for cv, d in costing_params_depth_dict.items():
        for p, v in d.items():
            cvp = getattr(dwi_params, f"{cv}_capital_cost_{p}")
            cvp.set_value(v)

    blk.drilling_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Drilling capital costs",
    )

    blk.piping_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Injection piping capital costs",
    )

    blk.packing_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Packing capital costs",
    )

    blk.casing_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Casing capital costs",
    )

    blk.grouting_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Grouting capital costs",
    )

    blk.monitoring_well_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Monitoring well capital costs",
    )

    blk.mobilization_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Mobilization capital costs",
    )

    blk.logging_testing_capital_cost = Var(
        initialize=1e5,
        domain=NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Logging + testing capital costs",
    )

    capital_cost_expr = 0

    blk.drilling_capital_cost_constraint = Constraint(
        expr=blk.drilling_capital_cost
        == pyunits.convert(
            dwi_params.drilling_capital_cost_slope * dwi.pipe_diameter
            + dwi_params.drilling_capital_cost_intercept,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.drilling_capital_cost

    blk.piping_capital_cost_constraint = Constraint(
        expr=blk.piping_capital_cost
        == pyunits.convert(
            dwi_params.piping_capital_cost_base
            * pipe_diam_dimensionless**dwi_params.piping_capital_cost_exponent,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.piping_capital_cost

    blk.packing_capital_cost_constraint = Constraint(
        expr=blk.packing_capital_cost
        == pyunits.convert(
            dwi_params.packing_capital_cost_slope * log(pipe_diam_dimensionless)
            + dwi_params.packing_capital_cost_intercept,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.packing_capital_cost

    blk.casing_capital_cost_constraint = Constraint(
        expr=blk.casing_capital_cost
        == pyunits.convert(
            dwi_params.casing_capital_cost_slope * dwi.pipe_diameter
            + dwi_params.casing_capital_cost_intercept,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.casing_capital_cost

    blk.grouting_capital_cost_constraint = Constraint(
        expr=blk.grouting_capital_cost
        == pyunits.convert(
            dwi_params.grouting_capital_cost_base
            * exp(dwi_params.grouting_capital_cost_exp_coeff * dwi.pipe_diameter),
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.grouting_capital_cost

    blk.monitoring_well_capital_cost_constraint = Constraint(
        expr=blk.monitoring_well_capital_cost
        == pyunits.convert(
            dwi_params.monitoring_well_capital_cost_base
            * monitoring_well_depth_dimensionless
            ** dwi_params.monitoring_well_capital_cost_exponent,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.monitoring_well_capital_cost

    blk.mobilization_capital_cost_constraint = Constraint(
        expr=blk.mobilization_capital_cost
        == pyunits.convert(
            dwi_params.mobilization_capital_cost_base
            * injection_well_depth_dimensionless
            ** dwi_params.mobilization_capital_cost_exponent,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.mobilization_capital_cost

    blk.logging_testing_capital_cost_constraint = Constraint(
        expr=blk.logging_testing_capital_cost
        == pyunits.convert(
            dwi_params.logging_testing_capital_cost_slope * dwi.pipe_diameter
            + dwi_params.logging_testing_capital_cost_intercept,
            to_units=blk.costing_package.base_currency,
        )
    )

    capital_cost_expr += blk.logging_testing_capital_cost

    blk.costing_package.add_cost_factor(blk, None)
    blk.capital_cost_constraint = Constraint(expr=blk.capital_cost == capital_cost_expr)

    blk.pumping_power_required = Var(
        initialize=1e3,
        bounds=(0, None),
        units=pyunits.kilowatt,
        doc="Pumping power required",
    )

    # calculate power required based on user provided required pumping pressure
    # TODO: add capital cost for pump?
    blk.pumping_power_required_constraint = Constraint(
        expr=blk.pumping_power_required
        == pyunits.convert(
            dwi.properties[0].flow_vol_phase["Liq"] * dwi.injection_pressure,
            to_units=pyunits.kilowatt,
        )
    )
    blk.costing_package.cost_flow(blk.pumping_power_required, "electricity")
