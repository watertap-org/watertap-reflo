import pyomo.environ as pyo
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.reflo.costing.util import (
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)
from idaes.core.util.constants import Constants

# Costing equations from:



def build_air_stripping_cost_param_block(blk):

    costing = blk.parent_block()

    blk.pressure_ambient = pyo.Param(
        initialize=101325,
        units=pyo.units.Pa,
        mutable=True,
        doc="Ambient pressure",
    )

    blk.blower_efficiency = pyo.Var(
        initialize=0.4,
        units=pyo.units.dimensionless,
        bounds=(0, 1),
        doc="Blower efficiency",
    )

    blk.pump_efficiency = pyo.Var(
        initialize=0.85,
        units=pyo.units.dimensionless,
        bounds=(0, 1),
        doc="Pump efficiency",
    )

    blk.power_blower_denom_coeff = pyo.Var(
        initialize=0.283,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Blower power equation denominator coefficient",
    )

    blk.power_blower_exponent = pyo.Var(
        initialize=0.283,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Blower power equation exponent",
    )

    blk.capital_cost_tower_A_param = pyo.Var(
        initialize=45.2,
        units=pyo.units.USD_1991/pyo.units.feet,
        bounds=(0, None),
        doc="Tower capital cost A parameter",
    )

    blk.capital_cost_tower_B_param = pyo.Var(
        initialize=3.5,
        units=pyo.units.USD_1991/(pyo.units.inch*pyo.units.ft),
        bounds=(0, None),
        doc="Tower capital cost B parameter",
    )

    blk.capital_cost_tower_C_param = pyo.Var(
        initialize=-0.0077,
        units=pyo.units.USD_1991/(pyo.units.inch**2*pyo.units.ft),
        bounds=(0, None),
        doc="Tower capital cost C parameter",
    )

    blk.capital_cost_port_A_param = pyo.Var(
        initialize=-31.6,
        units=pyo.units.USD_1991,
        bounds=(0, None),
        doc="Port capital cost A parameter",
    )

    blk.capital_cost_port_B_param = pyo.Var(
        initialize=72.8,
        units=pyo.units.USD_1991/pyo.units.inch,
        bounds=(0, None),
        doc="Port capital cost B parameter",
    )

    blk.capital_cost_port_C_param = pyo.Var(
        initialize=-2.8,
        units=pyo.units.USD_1991/pyo.units.inch**2,
        bounds=(0, None),
        doc="Port capital cost C parameter",
    )

    blk.capital_cost_port_D_param = pyo.Var(
        initialize=0.11,
        units=pyo.units.USD_1991/pyo.units.inch**3,
        bounds=(0, None),
        doc="Port capital cost D parameter",
    )

    blk.capital_cost_pipe_A_param = pyo.Var(
        initialize=133.8,
        units=pyo.units.USD_1991,
        bounds=(0, None),
        doc="Liquid inlet/outlet piping capital cost A parameter",
    )

    blk.capital_cost_pipe_B_param = pyo.Var(
        initialize=42,
        units=pyo.units.USD_1991/pyo.units.inch,
        bounds=(0, None),
        doc="Liquid inlet/outlet piping capital cost B parameter",
    )

    blk.capital_cost_pipe_C_param = pyo.Var(
        initialize=4.8,
        units=pyo.units.USD_1991/pyo.units.inch**2,
        bounds=(0, None),
        doc="Liquid inlet/outlet piping capital cost C parameter",
    )

    blk.capital_cost_pipe_air_param = pyo.Var(
        initialize=1.05,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Air inlet piping capital cost parameter",
    )

    blk.capital_cost_tray_rings_A_param = pyo.Var(
        initialize=70.4,
        units=pyo.units.USD_1991,
        bounds=(0, None),
        doc="Tray rings capital cost A parameter",
    )

    blk.capital_cost_tray_rings_B_param = pyo.Var(
        initialize=4.45,
        units=pyo.units.USD_1991/pyo.units.inch,
        bounds=(0, None),
        doc="Tray rings capital cost B parameter",
    )

    blk.capital_cost_tray_rings_C_param = pyo.Var(
        initialize=1.73e-2,
        units=pyo.units.USD_1991/pyo.units.inch**2,
        bounds=(0, None),
        doc="Tray rings capital cost C parameter",
    )

    blk.capital_cost_tray_A_param = pyo.Var(
        initialize=658.1,
        units=pyo.units.USD_1991,
        bounds=(0, None),
        doc="Tray capital cost A parameter",
    )

    blk.capital_cost_tray_B_param = pyo.Var(
        initialize=-6.5,
        units=pyo.units.USD_1991/pyo.units.inch,
        bounds=(0, None),
        doc="Tray capital cost B parameter",
    )

    blk.capital_cost_tray_C_param = pyo.Var(
        initialize=0.22,
        units=pyo.units.USD_1991/pyo.units.inch**2,
        bounds=(0, None),
        doc="Tray capital cost C parameter",
    )

    blk.capital_cost_plate_A_param = pyo.Var(
        initialize=20.6,
        units=pyo.units.USD_1991,
        bounds=(0, None),
        doc="Plate capital cost A parameter",
    )

    blk.capital_cost_plate_B_param = pyo.Var(
        initialize=1.1,
        units=pyo.units.USD_1991/pyo.units.inch,
        bounds=(0, None),
        doc="Plate capital cost B parameter",
    )

    blk.capital_cost_plate_C_param = pyo.Var(
        initialize=9.7e-2,
        units=pyo.units.USD_1991/pyo.units.inch**2,
        bounds=(0, None),
        doc="Plate capital cost C parameter",
    )

    blk.capital_cost_mister_A_param = pyo.Var(
        initialize=46.4,
        units=pyo.units.USD_1991,
        bounds=(0, None),
        doc="Mister capital cost A parameter",
    )

    blk.capital_cost_mister_B_param = pyo.Var(
        initialize=9.3,
        units=pyo.units.USD_1991/pyo.units.inch,
        bounds=(0, None),
        doc="Mister capital cost B parameter",
    )

    blk.capital_cost_mister_C_param = pyo.Var(
        initialize=0.14,
        units=pyo.units.USD_1991/pyo.units.inch**2,
        bounds=(0, None),
        doc="Mister capital cost C parameter",
    )

    blk.capital_cost_pump_base_param = pyo.Var(
        initialize=9.84e3,
        units=pyo.units.USD_2005,
        bounds=(0, None),
        doc="Water pump capital cost base",
    )

    blk.capital_cost_pump_denom_param = pyo.Var(
        initialize=4,
        units=pyo.units.kilowatt,
        bounds=(0, None),
        doc="Water pump capital cost denominator parameter",
    )

    blk.capital_cost_pump_exponent = pyo.Var(
        initialize=0.55,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Water pump capital cost exponent",
    )

    blk.capital_cost_blower_intercept = pyo.Var(
        initialize=4450,
        units=pyo.units.USD_2013,
        bounds=(0, None),
        doc="Blower capital cost intercept",
    )

    blk.capital_cost_blower_base = pyo.Var(
        initialize=57,
        units=(pyo.units.USD_2013*pyo.units.hour)/pyo.units.m**3,
        bounds=(0, None),
        doc="Blower capital cost base",
    )

    blk.capital_cost_pump_exponent = pyo.Var(
        initialize=0.8,
        units=pyo.units.dimensionless,
        bounds=(0, None),
        doc="Blower capital cost exponent",
    )

    blk.capital_cost_packing = pyo.Var(
        initialize=5500,
        units=pyo.units.USD_2013/pyo.units.m**3,
        bounds=(0, None),
        doc="Capital cost of packing material per m3",
    )


@register_costing_parameter_block(
    build_rule=build_air_stripping_cost_param_block,
    parameter_block_name="air_stripping",
)
def cost_air_stripping(blk):

    ax = blk.costing_package.air_stripping
    make_capital_cost_var(blk)
    make_fixed_operating_cost_var(blk)

    base_currency = blk.config.flowsheet_costing_block.base_currency

    blk.tower_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Aluminum tower cost",
    )

    blk.port_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Aluminum access port system cost",
    )

    blk.piping_liq_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Liquid inlet/outlet piping cost",
    )

    blk.piping_air_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Air inlet piping cost",
    )

    blk.tray_ring_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Tray ring cost",
    )

    blk.tower_internals_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Tray ring cost",
    )

    blk.packing_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Packing material cost",
    )

    blk.mist_eliminator_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Mist eliminator cost",
    )

    blk.pump_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Water pump cost",
    )

    blk.blower_cost = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=base_currency,
        doc="Air blower cost",
    )

    blk.blower_power = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=pyo.units.kilowatt,
        doc="Air blower power requirement",
    )

    blk.pump_power = pyo.Var(
        initialize=100,
        bounds=(0, None),
        units=pyo.units.kilowatt,
        doc="Water pump power requirement",
    )


    # blk.electricity_flow = pyo.Expression(
    #     expr=1
    # )
    # blk.costing_package.cost_flow(blk.electricity_flow, "electricity")
