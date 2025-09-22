import logging
import warnings

from pyomo.environ import (
    ConcreteModel,
    Block,
    Var,
    Constraint,
    Objective,
    NonNegativeReals,
    TransformationFactory,
    assert_optimal_termination,
    value,
    units as pyunits,
)
from pyomo.network import Arc

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.models.unit_models import (
    Product,
    Feed,
    Mixer,
    MixingType,
    MomentumMixingType,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import (
    FlowsheetBlock,
    UnitModelCostingBlock,
    MaterialFlowBasis,
    MaterialBalanceType,
)
from idaes.core.util.initialization import propagate_state

from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.property_models.NaCl_prop_pack import NaClParameterBlock
from watertap.property_models.multicomp_aq_sol_prop_pack import MCASParameterBlock
from watertap.property_models.seawater_prop_pack import SeawaterParameterBlock
from watertap.property_models.unit_specific.cryst_prop_pack import (
    NaClParameterBlock as CrystParameterBlock,
)
from watertap.property_models.water_prop_pack import (
    WaterParameterBlock as SteamParameterBlock,
)
from watertap.unit_models.pressure_changer import Pump

from watertap_contrib.reflo.core import REFLODatabase
from watertap_contrib.reflo.costing import (
    TreatmentCosting,
    EnergyCosting,
    REFLOSystemCosting,
)
from watertap_contrib.reflo.flowsheets.KBHDP.components import *
from watertap_contrib.reflo.flowsheets.KBHDP.utils import solve


# mute some of the logs
for x in (
    "idaes",
    "watertap",
    "idaes.core.util.scaling",
    "idaes.init.blocks",
    "idaes.core.base.costing_base",
):
    try:
        idaeslog.getLogger(x).setLevel(idaeslog.CRITICAL)
    except Exception:
        logging.getLogger(x).setLevel(logging.CRITICAL)
logging.getLogger().setLevel(logging.CRITICAL)
warnings.filterwarnings("ignore")


def build_zld_ro(ro_recovery=0.5, Qin=4):

    m = build_ro_treatment()
    add_ro_connections(m)
    set_ro_operating_conditions(m.fs.treatment, Qin)
    apply_ro_scaling(m)
    init_ro_treatment(m.fs.treatment)
    add_ro_recovery_constraint(m, m.fs.treatment.RO, ro_recovery)
    optimize_zld_ro(m, water_recovery=ro_recovery)
    results = solve(m)
    assert_optimal_termination(results)

    return m


def add_ro_recovery_constraint(m, blk, ro_recovery):
    m.fs.treatment.water_recovery = Var(
        initialize=ro_recovery,
        bounds=(0, 0.99),
        domain=NonNegativeReals,
        units=pyunits.dimensionless,
        doc="RO Water Recovery",
    )

    blk.eq_water_recovery = Constraint(
        expr=blk.feed.properties[0].flow_vol * m.fs.treatment.water_recovery
        == blk.product.properties[0].flow_vol
    )


def build_ro_treatment():
    m = ConcreteModel()
    m.db = REFLODatabase()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.MCAS_properties = MCASParameterBlock(
        solute_list=[
            "Alkalinity_2-",
            "Ca_2+",
            "Cl_-",
            "Mg_2+",
            "K_+",
            "SiO2",
            "Na_+",
            "SO4_2-",
        ],
        material_flow_basis=MaterialFlowBasis.mass,
    )

    m.fs.RO_properties = NaClParameterBlock()
    m.fs.UF_properties = WaterParameterBlock(solute_list=["tds", "tss"])
    m.fs.cryst_properties = CrystParameterBlock()
    m.fs.properties_md = SeawaterParameterBlock()
    m.fs.vapor_properties = SteamParameterBlock()

    treatment = m.fs.treatment = Block()

    treatment.feed = Feed(property_package=m.fs.MCAS_properties)
    treatment.RO_waste = Product(property_package=m.fs.RO_properties)
    treatment.sludge = Product(property_package=m.fs.UF_properties)
    treatment.UF_waste = Product(property_package=m.fs.UF_properties)

    treatment.EC = FlowsheetBlock(dynamic=False)
    treatment.UF = FlowsheetBlock(dynamic=False)
    treatment.pump = Pump(property_package=m.fs.RO_properties)
    treatment.RO = FlowsheetBlock(dynamic=False)

    treatment.MCAS_to_TDS_translator = TranslatorMCAStoZO(
        inlet_property_package=m.fs.MCAS_properties,
        outlet_property_package=m.fs.UF_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=False,
    )

    treatment.TDS_to_NaCl_translator = TranslatorZOtoNaCl(
        inlet_property_package=m.fs.UF_properties,
        outlet_property_package=m.fs.RO_properties,
        has_phase_equilibrium=False,
        outlet_state_defined=True,
    )

    treatment.RO_product_to_translator = TranslatorNaCltoNaCl(
        inlet_property_package=m.fs.RO_properties,
        outlet_property_package=m.fs.cryst_properties,
    )

    build_EC(treatment.EC, prop_package=m.fs.UF_properties)
    build_UF(treatment.UF, prop_package=m.fs.UF_properties)
    build_ro(treatment.RO, prop_package=m.fs.RO_properties, number_of_stages=1)

    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-1, index=("Liq", "H2O")
    )
    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp", 10**-1, index=("Liq", "NaCl")
    )

    treatment.EC.feed.properties[0].conc_mass_comp
    treatment.RO.disposal.properties[0].conc_mass_phase_comp

    m.fs.UF_properties.set_default_scaling("flow_mass_comp", 1e-2, index=("H2O"))
    m.fs.UF_properties.set_default_scaling("flow_mass_comp", 1, index=("tds"))
    m.fs.UF_properties.set_default_scaling("flow_mass_comp", 1e5, index=("tss"))

    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 1, index=("Liq", "H2O")
    )
    m.fs.RO_properties.set_default_scaling(
        "flow_mass_phase_comp", 1e2, index=("Liq", "NaCl")
    )

    return m


def add_ro_connections(m):
    treatment = m.fs.treatment

    treatment.feed_to_translator = Arc(
        source=treatment.feed.outlet,
        destination=treatment.MCAS_to_TDS_translator.inlet,
    )

    treatment.translator_to_EC = Arc(
        source=treatment.MCAS_to_TDS_translator.outlet,
        destination=treatment.EC.feed.inlet,
    )

    treatment.EC_to_UF = Arc(
        source=treatment.EC.product.outlet,
        destination=treatment.UF.feed.inlet,
    )

    treatment.EC_to_sludge = Arc(
        source=treatment.EC.disposal.outlet,
        destination=treatment.sludge.inlet,
    )

    treatment.UF_to_translator3 = Arc(
        source=treatment.UF.product.outlet,
        destination=treatment.TDS_to_NaCl_translator.inlet,
    )

    treatment.UF_to_waste = Arc(
        source=treatment.UF.disposal.outlet,
        destination=treatment.UF_waste.inlet,
    )

    treatment.translator_to_pump = Arc(
        source=treatment.TDS_to_NaCl_translator.outlet,
        destination=treatment.pump.inlet,
    )

    treatment.pump_to_ro = Arc(
        source=treatment.pump.outlet,
        destination=treatment.RO.feed.inlet,
    )

    treatment.ro_to_product = Arc(
        source=treatment.RO.product.outlet,
        destination=treatment.RO_product_to_translator.inlet,
    )

    treatment.ro_to_md = Arc(
        source=treatment.RO.disposal.outlet,
        destination=treatment.RO_waste.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def set_inlet_conditions(m, Qin=4):

    print(f'\n{"=======> SETTING OPERATING CONDITIONS <=======":^60}\n')

    Qin = Qin * pyunits.Mgallon / pyunits.day
    rho = 1000 * pyunits.kg / pyunits.m**3
    print('\n=======> SETTING FEED CONDITIONS <======="\n')

    inlet_dict = {
        "Ca_2+": 0.61 * pyunits.kg / pyunits.m**3,
        "Mg_2+": 0.161 * pyunits.kg / pyunits.m**3,
        "Alkalinity_2-": 0.0821 * pyunits.kg / pyunits.m**3,
        "SiO2": 0.13 * pyunits.kg / pyunits.m**3,
        "Cl_-": 5.5 * pyunits.kg / pyunits.m**3,
        "Na_+": 5.5 * pyunits.kg / pyunits.m**3,
        "K_+": 0.016 * pyunits.kg / pyunits.m**3,
        "SO4_2-": 0.23 * pyunits.kg / pyunits.m**3,
    }

    # initialize feed
    m.fs.treatment.feed.pressure[0].fix(101325)
    m.fs.treatment.feed.temperature[0].fix(293)
    m.fs.treatment.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(Qin * rho)
    m.fs.treatment.feed.properties[0].flow_vol_phase

    for solute, solute_conc in inlet_dict.items():
        m.fs.treatment.feed.properties[0].flow_mass_phase_comp["Liq", solute].fix(
            Qin * solute_conc
        )
        m.fs.MCAS_properties.set_default_scaling(
            "flow_mass_phase_comp",
            1
            / value(
                m.fs.treatment.feed.properties[0].flow_mass_phase_comp["Liq", solute]
            ),
            index=("Liq", solute),
        )

    m.fs.MCAS_properties.set_default_scaling(
        "flow_mass_phase_comp",
        1 / value(m.fs.treatment.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"]),
        index=("Liq", "H2O"),
    )


def set_ro_operating_conditions(blk, Qin=4, RO_pressure=20e5):
    m = blk.model()
    # Set inlet conditions and operating conditions for each unit
    set_inlet_conditions(m, Qin=Qin)
    set_EC_operating_conditions(blk.EC)
    set_UF_op_conditions(blk.UF)
    blk.pump.efficiency_pump.fix(0.8)
    blk.pump.control_volume.properties_out[0].pressure.fix(RO_pressure)
    set_ro_system_operating_conditions(blk.RO, mem_area=10000)


def apply_ro_scaling(m):

    add_ec_scaling(m.fs.treatment.EC)
    add_UF_scaling(m.fs.treatment.UF)
    add_ro_scaling(m.fs.treatment.RO)

    iscale.set_scaling_factor(
        m.fs.treatment.UF_waste.properties[0.0].flow_mass_comp["tds"], 1e3
    )
    iscale.calculate_scaling_factors(m)


def init_ro_treatment(blk):

    blk.feed.initialize()
    propagate_state(blk.feed_to_translator)

    blk.MCAS_to_TDS_translator.initialize()
    propagate_state(blk.translator_to_EC)

    init_EC(blk.EC)
    propagate_state(blk.EC_to_UF)
    propagate_state(blk.EC_to_sludge)
    blk.sludge.initialize()

    init_UF(blk.UF)
    propagate_state(blk.UF_to_translator3)
    propagate_state(blk.UF_to_waste)
    blk.UF_waste.initialize()

    blk.TDS_to_NaCl_translator.initialize()
    propagate_state(blk.translator_to_pump)

    blk.pump.initialize()

    propagate_state(blk.pump_to_ro)

    init_ro_system(blk.RO)
    propagate_state(blk.ro_to_product)
    propagate_state(blk.ro_to_md)

    blk.RO.disposal.properties[0].flow_vol_phase
    blk.RO.disposal.initialize()


def optimize_zld_ro(m, water_recovery=0.5, fixed_pressure=None, ro_mem_area=None):
    treatment = m.fs.treatment
    print("\n\nDOF before optimization: ", degrees_of_freedom(m))

    if water_recovery is not None:
        print(f"\n------- Fixed Recovery at {100*water_recovery}% -------")
        m.fs.treatment.water_recovery.fix(water_recovery)
    else:
        lower_bound = 0.5
        upper_bound = 0.8
        print(f"\n------- Unfixed Recovery -------")
        print(f"Lower Bound: {lower_bound}")
        print(f"Upper Bound: {upper_bound}")
        m.fs.treatment.water_recovery.unfix()
        m.fs.treatment.water_recovery.setlb(lower_bound)
        m.fs.treatment.water_recovery.setub(upper_bound)

    if fixed_pressure is not None:
        print(f"\n------- Fixed RO Pump Pressure at {fixed_pressure} -------\n")
        treatment.pump.control_volume.properties_out[0].pressure.fix(fixed_pressure)
    else:
        lower_bound = 100 * pyunits.psi
        upper_bound = 900 * pyunits.psi
        print(f"------- Unfixed RO Pump Pressure -------")
        print(f"Lower Bound: {value(lower_bound)} {pyunits.get_units(lower_bound)}")
        print(f"Upper Bound: {value(upper_bound)} {pyunits.get_units(upper_bound)}")
        treatment.pump.control_volume.properties_out[0].pressure.unfix()
        treatment.pump.control_volume.properties_out[0].pressure.setlb(lower_bound)
        treatment.pump.control_volume.properties_out[0].pressure.setub(upper_bound)

    if ro_mem_area is not None:
        print(f"\n------- Fixed RO Membrane Area at {ro_mem_area} -------\n")
        treatment.RO.total_membrane_area.fix(ro_mem_area)
        for _, stage in treatment.RO.stage.items():
            stage.module.area.unfix()
            stage.module.area.setub(1e6)
    else:
        lower_bound = 1e3
        upper_bound = 2e5
        print(f"\n------- Unfixed RO Membrane Area -------")
        print(f"Lower Bound: {lower_bound} m2")
        print(f"Upper Bound: {upper_bound} m2")
        print("\n")
        for _, stage in treatment.RO.stage.items():
            stage.module.area.unfix()
            stage.module.area.setub(1e6)

    for _, stage in treatment.RO.stage.items():
        stage.module.width.setub(5000)
        stage.module.feed_side.velocity[0, 0].unfix()
        stage.module.feed_side.velocity[0, 1].setlb(0.0)
        stage.module.feed_side.K.setlb(1e-6)
        stage.module.feed_side.friction_factor_darcy.setub(50)
        stage.module.flux_mass_phase_comp.setub(1)
        # stage.module.flux_mass_phase_comp.setlb(1e-5)
        stage.module.feed_side.cp_modulus.setub(10)
        stage.module.rejection_phase_comp.setlb(1e-4)
        stage.module.feed_side.N_Re.setlb(1)
        stage.module.recovery_mass_phase_comp.setlb(1e-7)


def add_zld_md(m, md_water_recovery=0.5):

    treat = m.fs.treatment

    Q_md = value(
        pyunits.convert(
            m.fs.treatment.RO.disposal.properties[0].flow_vol_phase["Liq"],
            to_units=pyunits.m**3 / pyunits.s,
        )
    )
    Cin = value(
        pyunits.convert(
            m.fs.treatment.RO.disposal.properties[0].conc_mass_phase_comp[
                "Liq", "NaCl"
            ],
            to_units=pyunits.g / pyunits.L,
        )
    )

    m.inlet_flow_rate = pyunits.convert(
        Q_md * pyunits.m**3 / pyunits.s, to_units=pyunits.m**3 / pyunits.s
    )
    m.inlet_flow_rate = m.fs.treatment.RO.disposal.properties[0].flow_vol_phase["Liq"]
    m.inlet_salinity = pyunits.convert(
        Cin * pyunits.g / pyunits.liter, to_units=pyunits.kg / pyunits.m**3
    )
    m.water_recovery = md_water_recovery

    m.fs.treatment.sw_to_nacl_product = TranslatorSWtoNaCl(
        inlet_property_package=m.fs.properties_md,
        outlet_property_package=m.fs.cryst_properties,
    )

    m.fs.treatment.sw_to_nacl_disposal = TranslatorSWtoNaCl(
        inlet_property_package=m.fs.properties_md,
        outlet_property_package=m.fs.cryst_properties,
    )

    m.fs.treatment.md = FlowsheetBlock(dynamic=False)
    build_md(m.fs.treatment.md, prop_package=m.fs.properties_md)

    treat.md.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        treat.RO.disposal.properties[0].flow_mass_phase_comp["Liq", "H2O"].value
    )
    treat.md.feed.properties[0].flow_mass_phase_comp["Liq", "TDS"].fix(
        treat.RO.disposal.properties[0].flow_mass_phase_comp["Liq", "NaCl"].value
    )
    treat.md.feed.properties[0].temperature.fix(298.15)
    treat.md.feed.properties[0].pressure.fix(101325)
    treat.md.feed.properties[0].conc_mass_phase_comp
    treat.md.feed.properties[0].flow_vol_phase
    treat.md.feed.initialize()

    init_md(treat.md)

    treat.md_to_product = Arc(
        source=treat.md.permeate.outlet, destination=treat.sw_to_nacl_product.inlet
    )

    treat.md_disposal_to_nacl_translator = Arc(
        source=treat.md.concentrate.outlet,
        destination=treat.sw_to_nacl_disposal.inlet,
    )

    propagate_state(treat.md_to_product)
    treat.sw_to_nacl_product.initialize()

    propagate_state(treat.md_disposal_to_nacl_translator)
    treat.sw_to_nacl_disposal.initialize()

    TransformationFactory("network.expand_arcs").apply_to(m)

    results = solve(m.fs.treatment.md)
    assert_optimal_termination(results)


def add_zld_mec(m):

    print("\n--------- Adding MEC---------\n")
    treat = m.fs.treatment
    treat.mec = FlowsheetBlock(dynamic=False)

    build_MEC(
        treat.mec,
        prop_package=m.fs.cryst_properties,
        vapor_prop_package=m.fs.vapor_properties,
    )

    total_feed_H2O_mass = value(
        treat.sw_to_nacl_disposal.outlet.flow_mass_phase_comp[0, "Liq", "H2O"]
    )
    total_feed_NaCl_mass = value(
        treat.sw_to_nacl_disposal.outlet.flow_mass_phase_comp[0, "Liq", "NaCl"]
    )

    treat.mec.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].set_value(
        total_feed_H2O_mass
    )
    treat.mec.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].set_value(
        total_feed_NaCl_mass
    )
    treat.mec.unit.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].unfix()
    treat.mec.unit.inlet.flow_mass_phase_comp[0, "Sol", "NaCl"].set_value(0)
    # treat.mec.unit.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].unfix()
    # treat.mec.unit.inlet.flow_mass_phase_comp[0, "Vap", "H2O"].set_value(0)
    treat.mec.unit.inlet.temperature[0].set_value(273.15 + 30.51)
    treat.mec.unit.inlet.pressure[0].set_value(101325)

    set_MEC_op_conditions(
        treat.mec,
        feed_H2O=total_feed_H2O_mass,
        feed_NaCl=total_feed_NaCl_mass,
    )
    scale_MEC(
        treat.mec,
        cryst_prop_pack=m.fs.cryst_properties,
        vapor_prop_pack=m.fs.vapor_properties,
    )

    iscale.calculate_scaling_factors(treat.mec)

    init_MEC(treat.mec)
    report_MEC(treat.mec)

    treat.mec.feed.properties[0].flow_mass_phase_comp["Liq", "H2O"].fix(
        total_feed_H2O_mass
    )
    treat.mec.feed.properties[0].flow_mass_phase_comp["Liq", "NaCl"].fix(
        total_feed_NaCl_mass
    )
    treat.mec.feed.properties[0].flow_mass_phase_comp["Sol", "NaCl"].fix(0)
    # feed steam flow is dictated by unit.heating_steam[0] so we don't fix the inlet steam flow here
    # m.fs.feed.properties[0].flow_mass_phase_comp["Vap", "H2O"].fix(0)
    treat.mec.feed.properties[0].temperature.fix(
        value(treat.sw_to_nacl_disposal.outlet.temperature[0])
    )
    treat.mec.feed.properties[0].pressure.fix(
        value(treat.sw_to_nacl_disposal.outlet.pressure[0])
    )
    treat.mec.feed.properties[0].conc_mass_phase_comp[...]
    treat.mec.feed.properties[0].flow_vol_phase
    treat.mec.feed.initialize()

    treat.mec.product.properties[0].conc_mass_phase_comp[...]
    treat.mec.product.properties[0].flow_vol_phase
    treat.mec.product.initialize()

    results = solve(treat.mec)
    assert_optimal_termination(results)

    # rescale_MEC(
    #     treat.mec,
    #     flow_mass_phase_water_total=total_feed_H2O_mass,
    #     flow_mass_phase_salt_total=total_feed_NaCl_mass,
    #     cryst_prop_pack=m.fs.cryst_properties,
    #     vapor_prop_pack=m.fs.vapor_properties,
    # )


def add_product_stream(m):
    treat = m.fs.treatment

    treat.product_NaCl_mixer = Mixer(
        property_package=m.fs.cryst_properties,
        num_inlets=3,
        inlet_list=["ro_product", "md_product", "cryst_product"],
        material_balance_type=MaterialBalanceType.componentPhase,
        energy_mixing_type=MixingType.extensive,
        momentum_mixing_type=MomentumMixingType.none,
    )

    treat.product = Product(property_package=m.fs.cryst_properties)

    treat.md_translator_to_product_NaCl_mixer = Arc(
        source=treat.sw_to_nacl_product.outlet,
        destination=treat.product_NaCl_mixer.md_product,
    )

    treat.cryst_to_product_NaCl_mixer = Arc(
        source=treat.mec.unit.outlet,
        destination=treat.product_NaCl_mixer.cryst_product,
    )

    treat.ro_to_product_NaCl_mixer = Arc(
        source=treat.RO_product_to_translator.outlet,
        destination=treat.product_NaCl_mixer.ro_product,
    )

    treat.product_NaCl_mixer_to_product = Arc(
        source=treat.product_NaCl_mixer.outlet,
        destination=treat.product.inlet,
    )

    TransformationFactory("network.expand_arcs").apply_to(m)

    treat.sw_to_nacl_product.initialize()

    propagate_state(treat.md_translator_to_product_NaCl_mixer)
    propagate_state(treat.ro_to_product_NaCl_mixer)
    propagate_state(treat.cryst_to_product_NaCl_mixer)

    treat.product_NaCl_mixer.initialize()
    treat.product_NaCl_mixer.outlet.pressure[0].fix()

    propagate_state(treat.product_NaCl_mixer_to_product)

    treat.product.properties[0].flow_vol
    treat.product.properties[0].flow_vol_phase
    treat.product.initialize()


def add_treatment_costing(blk, heat_price=0, electricity_price=0):

    # Add treatment costing
    blk.costing = TreatmentCosting()

    blk.pump.costing = UnitModelCostingBlock(
        flowsheet_costing_block=blk.costing,
    )
    blk.costing.heat_cost.fix(heat_price)
    blk.costing.electricity_cost.fix(electricity_price)

    add_EC_costing(blk.EC, costing_block=blk.costing)
    add_UF_costing(blk.UF, costing_block=blk.costing)
    add_ro_costing(blk.RO, costing_block=blk.costing)

    iscale.constraint_scaling_transform(
        blk.EC.unit.costing.capital_cost_floc_constraint, 1e-6
    )
    iscale.constraint_scaling_transform(
        blk.EC.unit.costing.capital_cost_reactor_constraint, 1e-3
    )
    iscale.constraint_scaling_transform(
        blk.EC.unit.costing.capital_cost_power_supply_constraint, 1e-6
    )
    iscale.constraint_scaling_transform(
        blk.EC.unit.costing.capital_cost_electrodes_constraint, 1e-3
    )

    blk.md.unit.add_costing_module(blk.costing)

    add_MEC_costing(blk.mec, costing_block=blk.costing)

    iscale.constraint_scaling_transform(
        blk.mec.unit.costing.total_capital_cost_effect_1_constraint, 1e-7
    )
    iscale.constraint_scaling_transform(
        blk.mec.unit.costing.total_capital_cost_effect_2_constraint, 1e-7
    )
    iscale.constraint_scaling_transform(
        blk.mec.unit.costing.total_capital_cost_effect_3_constraint, 1e-7
    )
    iscale.constraint_scaling_transform(
        blk.mec.unit.costing.total_capital_cost_effect_4_constraint, 1e-7
    )

    blk.costing.cost_process()
    blk.costing.add_LCOW(blk.product.properties[0].flow_vol)

    blk.costing.initialize()

    print("\n--------- Treatment Costing Initialization Complete ---------\n")


def add_system_costing(m, heat_price=0.0166, electricity_price=0.04989):

    m.fs.costing = REFLOSystemCosting()
    m.fs.costing.heat_cost_buy.fix(heat_price)
    m.fs.costing.electricity_cost_buy.set_value(electricity_price)
    m.fs.costing.cost_process()

    m.fs.energy.costing.add_LCOE()
    m.fs.energy.costing.add_LCOH()
    m.fs.costing.add_LCOT(m.fs.treatment.product.properties[0].flow_vol)

    m.fs.energy.costing.initialize()
    m.fs.treatment.costing.initialize()
    m.fs.costing.initialize_build()

    print("\n--------- System Costing Initialization Complete ---------\n")


def build_energy(m):

    m.fs.energy = Block()
    m.fs.energy.costing = EnergyCosting()
    m.fs.energy.costing.heat_cost.set_value(0)
    m.fs.energy.costing.electricity_cost.fix(0)

    m.fs.energy.cst = FlowsheetBlock()

    build_CST(m.fs.energy.cst)
    set_CST_op_conditions(m.fs.energy.cst)
    # init_CST(m.fs.energy.cst)
    add_CST_costing(m.fs.energy.cst, costing_block=m.fs.energy.costing)
    add_CST_costing_scaling(m.fs.energy.cst)

    results = solve(m)
    assert_optimal_termination(results)

    m.fs.energy.pv = FlowsheetBlock()
    build_pv(m.fs.energy.pv)
    set_pv_op_conditions(m.fs.energy.pv)
    # init_pv(m.fs.energy.pv)
    add_pv_costing(m.fs.energy.pv, costing_block=m.fs.energy.costing)

    m.fs.energy.costing.cost_process()
    m.fs.energy.costing.initialize()

    iscale.constraint_scaling_transform(
        m.fs.energy.cst.unit.costing.direct_cost_constraint, 1e-5
    )

    results = solve(m)
    assert_optimal_termination(results)
    print("\n--------- Energy Costing Initialization Complete ---------\n")


def report_all_treatment_results(m, w=30):
    print(f"{'*' * (3 * w)}")
    print(f"{'*' * (3 * w)}")
    title = "KBHDP ZLD Stream Table"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    flow_in = pyunits.convert(
        m.fs.treatment.feed.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.Mgallons / pyunits.day,
    )
    flow_out = pyunits.convert(
        m.fs.treatment.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.Mgallons / pyunits.day,
    )
    inlet_salinity = pyunits.convert(
        m.fs.treatment.EC.feed.properties[0].conc_mass_comp["tds"],
        to_units=pyunits.g / pyunits.L,
    )
    sys_recovery = flow_out / flow_in
    # ------
    ro_feed_flow = pyunits.convert(
        m.fs.treatment.RO.feed.properties[0].flow_vol,
        to_units=pyunits.Mgallons / pyunits.day,
    )
    ro_product_flow = pyunits.convert(
        m.fs.treatment.RO.product.properties[0].flow_vol,
        to_units=pyunits.Mgallons / pyunits.day,
    )
    ro_brine_flow = pyunits.convert(
        m.fs.treatment.RO.disposal.properties[0].flow_vol,
        to_units=pyunits.Mgallons / pyunits.day,
    )
    ro_brine_conc = pyunits.convert(
        m.fs.treatment.RO.disposal.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
        to_units=pyunits.g / pyunits.L,
    )
    # ------
    md_feed_flow = pyunits.convert(
        m.inlet_flow_rate, to_units=pyunits.Mgallons / pyunits.day
    )
    md_product_flow = pyunits.convert(
        m.fs.treatment.md.permeate.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.Mgallons / pyunits.day,
    )
    md_brine_flow = pyunits.convert(
        m.fs.treatment.md.concentrate.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.Mgallons / pyunits.day,
    )
    md_feed_conc = pyunits.convert(
        m.fs.treatment.md.feed.properties[0].conc_mass_phase_comp["Liq", "TDS"],
        to_units=pyunits.g / pyunits.L,
    )
    md_brine_conc = pyunits.convert(
        m.fs.treatment.md.concentrate.properties[0].conc_mass_phase_comp["Liq", "TDS"],
        to_units=pyunits.g / pyunits.L,
    )
    # ------
    mec_feed_flow = pyunits.convert(
        m.fs.treatment.mec.feed.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.Mgallons / pyunits.day,
    )
    mec_product_flow = pyunits.convert(
        m.fs.treatment.mec.product.properties[0].flow_vol_phase["Liq"],
        to_units=pyunits.Mgallons / pyunits.day,
    )
    mec_brine_flow = value(mec_feed_flow) - value(mec_product_flow)
    mec_feed_conc = pyunits.convert(
        m.fs.treatment.mec.feed.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
        to_units=pyunits.g / pyunits.L,
    )
    mec_brine_conc = pyunits.convert(
        m.fs.treatment.mec.product.properties[0].conc_mass_phase_comp["Liq", "NaCl"],
        to_units=pyunits.g / pyunits.L,
    )
    mec_solids_flow = m.fs.treatment.mec.product.properties[0].flow_mass_phase_comp[
        "Sol", "NaCl"
    ]

    print(f'{"Feed Flow":<{w}s}{f"{value(flow_in):<{w},.2f}"}{"MGD":<{w}s}')
    print(f'{"Product Flow":<{w}s}{f"{value(flow_out):<{w},.2f}"}{"MGD":<{w}s}')
    print(
        f'{"System Recovery":<{w}s}{f"{value(sys_recovery)*100:<{w},.1f}"}{"%":<{w}s}'
    )
    print(f'{"Inlet Salinity":<{w}s}{f"{value(inlet_salinity):<{w},.1f}"}{"g/L":<{w}s}')
    print(f"{'-' * (1 * w)}")
    print(f'{"RO Feed Flow":<{w}s}{f"{value(ro_feed_flow):<{w},.2f}"}{"MGD":<{w}s}')
    print(
        f'{"RO Product Flow":<{w}s}{f"{value(ro_product_flow):<{w},.2f}"}{"MGD":<{w}s}'
    )
    print(f'{"RO Brine Flow":<{w}s}{f"{value(ro_brine_flow):<{w},.2f}"}{"MGD":<{w}s}')
    print(
        f'{"RO Brine Concentration":<{w}s}{f"{value(ro_brine_conc):<{w},.1f}"}{"g/L":<{w}s}'
    )
    print(f"{'-' * (1 * w)}")
    print(f'{"MD Feed Flow":<{w}s}{f"{value(md_feed_flow):<{w},.2f}"}{"MGD":<{w}s}')
    print(
        f'{"MD Product Flow":<{w}s}{f"{value(md_product_flow):<{w},.2f}"}{"MGD":<{w}s}'
    )
    print(f'{"MD Brine Flow":<{w}s}{f"{value(md_brine_flow):<{w},.2f}"}{"MGD":<{w}s}')
    print(f'{"MD Feed Salinity":<{w}s}{f"{value(md_feed_conc):<{w},.1f}"}{"g/L":<{w}s}')
    print(
        f'{"MD Brine Salinity":<{w}s}{f"{value(md_brine_conc):<{w},.1f}"}{"g/L":<{w}s}'
    )
    print(f"{'-' * (1 * w)}")
    print(f'{"MEC Feed Flow":<{w}s}{f"{value(mec_feed_flow):<{w},.2f}"}{"MGD":<{w}s}')
    print(
        f'{"MEC Product Flow":<{w}s}{f"{value(mec_product_flow):<{w},.2f}"}{"MGD":<{w}s}'
    )
    print(f'{"MEC Brine Flow":<{w}s}{f"{mec_brine_flow:<{w},.2f}"}{"MGD":<{w}s}')
    print(
        f'{"MEC Feed Salinity":<{w}s}{f"{value(mec_feed_conc):<{w},.1f}"}{"g/L":<{w}s}'
    )
    print(
        f'{"MEC Brine Salinity":<{w}s}{f"{value(mec_brine_conc):<{w},.1f}"}{"g/L":<{w}s}'
    )
    print(
        f'{"MEC Solids Flow":<{w}s}{f"{value(mec_solids_flow):<{w},.1f}"}{"kg/s":<{w}s}'
    )

    title = "KBHDP ZLD Unit Results"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")

    report_EC(m.fs.treatment.EC, w=w)
    report_UF(m.fs.treatment.UF, w=w)
    report_RO(m.fs.treatment.RO, w=w)
    report_MD(m.fs.treatment.md, w=w)
    report_MEC(m.fs.treatment.mec, w=w)


def report_all_energy_results(m, w=30):
    print(f"{'*' * (3 * w)}")
    print(f"{'*' * (3 * w)}")
    title = "KBHDP ZLD Energy Results"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")

    report_CST(m.fs.energy.cst, w=w)
    report_PV(m.fs.energy.pv, w=w)


def display_costing_results(m, w=30):

    print(f"{'*' * (3 * w)}")
    print(f"{'*' * (3 * w)}")

    title = "KBHDP ZLD Unit Costing Results"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")

    print_EC_costing_breakdown(m.fs.treatment.EC, w=w)
    print_UF_costing_breakdown(m.fs.treatment.UF, w=w)
    print_RO_costing_breakdown(m.fs.treatment.RO, w=w)
    print_MD_costing_breakdown(m.fs.treatment.md, w=w)
    print_MEC_costing_breakdown(m.fs.treatment.mec, w=w)
    print_CST_costing_breakdown(m.fs.energy.cst, w=w)
    print_PV_costing_breakdown(m.fs.energy.pv, w=w)

    title = "KBHDP ZLD System Costing Results"
    side = int(((3 * w) - len(title)) / 2) - 1
    header = "=" * side + f" {title} " + "=" * side
    print(f"\n{header}\n")
    print(f'{"Parameter":<{w}s}{"Value":<{w}s}{"Units":<{w}s}')
    print(f"{'-' * (3 * w)}")
    print(f'{"LCOT":<{w}s}{f"{value(m.fs.costing.LCOT):<{w},.2f}"}{"$/m3":<{w}s}')
    print(
        f'{"LCOW":<{w}s}{f"{value(m.fs.treatment.costing.LCOW):<{w},.2f}"}{"$/m3":<{w}s}'
    )
    print(
        f'{"LCOE":<{w}s}{f"{value(m.fs.energy.costing.LCOE):<{w},.2f}"}{"$/kWh":<{w}s}'
    )
    print(
        f'{"LCOH":<{w}s}{f"{value(m.fs.energy.costing.LCOH):<{w},.2f}"}{"$/kWh":<{w}s}'
    )
    print(f"{'-' * (1 * w)}")
    print(
        f'{"Grid Frac. Electricity":<{w}s}{f"{value(m.fs.costing.frac_elec_from_grid):<{w},.3f}"}{"-":<{w}s}'
    )
    print(
        f'{"Grid Frac. Heat":<{w}s}{f"{value(m.fs.costing.frac_heat_from_grid):<{w},.3f}"}{"-":<{w}s}'
    )
    print(f"{'-' * (1 * w)}")
    print(
        f'{"Total Treatment CAPEX":<{w}s}{f"{value(m.fs.treatment.costing.total_capital_cost):<{w},.0f}"}{"$":<{w}s}'
    )
    print(
        f'{"Total Treatment OPEX":<{w}s}{f"{value(m.fs.treatment.costing.total_operating_cost):<{w},.0f}"}{"$/year":<{w}s}'
    )
    print(f"{'-' * (1 * w)}")
    print(
        f'{"Total Energy CAPEX":<{w}s}{f"{value(m.fs.energy.costing.total_capital_cost):<{w},.0f}"}{"$":<{w}s}'
    )
    print(
        f'{"Total Energy OPEX":<{w}s}{f"{value(m.fs.energy.costing.total_operating_cost):<{w},.0f}"}{"$/year":<{w}s}'
    )
    print(f"{'-' * (1 * w)}")
    print(
        f'{"Total CAPEX":<{w}s}{f"{value(m.fs.costing.total_capital_cost):<{w},.0f}"}{"$":<{w}s}'
    )
    print(
        f'{"Total OPEX":<{w}s}{f"{value(m.fs.costing.total_operating_cost):<{w},.0f}"}{"$/year":<{w}s}'
    )
    print(f"{'.' * (3 * w)}")
    for flow in m.fs.treatment.costing.aggregate_flow_costs:
        if flow == "electricity":
            print(
                f'{f"{flow.upper()} Flow":<{w}s}{f"{m.fs.costing.aggregate_flow_electricity_purchased():<{w},.0f}"}{pyunits.get_units(m.fs.costing.aggregate_flow_electricity_purchased)}'
            )
        elif flow == "heat":
            print(
                f'{f"{flow.upper()} Flow":<{w}s}{f"{m.fs.costing.aggregate_flow_heat_purchased():<{w},.0f}"}{pyunits.get_units(m.fs.costing.aggregate_flow_heat_purchased)}'
            )
        else:
            print(
                f'{f"{flow.upper()} Cost":<{w}s}{f"${m.fs.treatment.costing.aggregate_flow_costs[flow]():<{w},.0f}"}{pyunits.get_units(m.fs.treatment.costing.aggregate_flow_costs[flow])}'
            )
    print("\n\n")


def report_all_results(m, w=30):
    report_all_treatment_results(m, w=w)
    report_all_energy_results(m, w=w)
    display_costing_results(m, w=w)


def main(Qin=4, ro_recovery=0.8, md_water_recovery=0.7):

    m = build_zld_ro(ro_recovery=ro_recovery, Qin=Qin)

    build_energy(m)

    add_zld_md(m, md_water_recovery=md_water_recovery)
    results = solve(m)
    assert_optimal_termination(results)

    add_zld_mec(m)
    add_product_stream(m)
    results = solve(m.fs.treatment.mec)
    assert_optimal_termination(results)
    results = solve(m)
    assert_optimal_termination(results)

    add_treatment_costing(m.fs.treatment, heat_price=0, electricity_price=0)
    results = solve(m)
    assert_optimal_termination(results)

    # build_energy(m)

    add_system_costing(m)
    m.fs.energy.cst.unit.system_capacity.unfix()
    m.fs.energy.pv.unit.system_capacity.unfix()

    m.fs.obj = Objective(expr=m.fs.costing.LCOT)
    # m.fs.energy.cst.unit.system_capacity.setlb(25)
    iscale.calculate_scaling_factors(m)

    results = solve(m)
    assert_optimal_termination(results)

    # m.fs.energy.cst.unit.system_capacity.fix()
    # m.fs.energy.pv.unit.system_capacity.fix()

    # results = solve(m)
    # assert_optimal_termination(results)

    m.fs.costing.frac_elec_from_grid.fix(0.9)
    m.fs.costing.frac_heat_from_grid.fix(0.9)

    results = solve(m)
    assert_optimal_termination(results)

    m.fs.costing.frac_elec_from_grid.fix(0.5)
    m.fs.costing.frac_heat_from_grid.fix(0.5)
    results = solve(m)
    assert_optimal_termination(results)

    report_all_results(m, w=25)

    m.fs.treatment.sw_to_nacl_product.properties_in[0].temperature.fix()
    m.fs.treatment.sw_to_nacl_product.properties_in[0].pressure.fix()

    print(f"DOF FINAL = {degrees_of_freedom(m)}")

    return m


if __name__ == "__main__":

    m = main(ro_recovery=0.8, md_water_recovery=0.7)
