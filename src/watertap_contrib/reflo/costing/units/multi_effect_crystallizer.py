#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pyomo.environ as pyo
from idaes.core.util.misc import StrEnum
from idaes.core.util.constants import Constants
from idaes.core.util.exceptions import ConfigurationError
from watertap.costing.util import register_costing_parameter_block
from watertap_contrib.reflo.costing.util import make_capital_cost_var


class MultiEffectCrystallizerCostType(StrEnum):
    mass_basis = "mass_basis"
    volume_basis = "volume_basis"


def build_recovered_nacl_cost_param_block(blk):

    blk.cost = pyo.Param(
        mutable=True,
        initialize=0,
        doc="Recovered cost (sale price) of NaCl",
        units=pyo.units.USD_2020 / pyo.units.kg,
    )

    costing = blk.parent_block()
    costing.register_flow_type("NaCl_recovered", blk.cost)


def build_steam_cost_param_block(blk):

    blk.cost = pyo.Param(
        mutable=True,
        initialize=0.004,
        doc="Steam cost",
        units=pyo.units.USD_2018 * pyo.units.meter**-3,
    )

    costing = blk.parent_block()
    costing.register_flow_type("steam", blk.cost)


def build_multi_effect_crystallizer_cost_param_block(blk):

    # Crystallizer operating cost information from literature
    blk.fob_unit_cost = pyo.Var(
        initialize=675000,
        doc="Forced circulation crystallizer reference free-on-board cost (Woods, 2007)",
        units=pyo.units.USD_2007,
    )

    blk.ref_capacity = pyo.Var(
        initialize=1,
        doc="Forced circulation crystallizer reference crystal capacity (Woods, 2007)",
        units=pyo.units.kg / pyo.units.s,
    )

    blk.ref_exponent = pyo.Var(
        initialize=0.53,
        doc="Forced circulation crystallizer cost exponent factor (Woods, 2007)",
        units=pyo.units.dimensionless,
    )

    blk.iec_percent = pyo.Var(
        initialize=1.43,
        doc="Forced circulation crystallizer installed equipment cost (Diab and Gerogiorgis, 2017)",
        units=pyo.units.dimensionless,
    )

    blk.volume_cost = pyo.Var(
        initialize=16320,
        doc="Forced circulation crystallizer cost per volume (Yusuf et al., 2019)",
        units=pyo.units.USD_2007,  ## TODO: Needs confirmation, but data is from Perry apparently
    )

    blk.vol_basis_exponent = pyo.Var(
        initialize=0.47,
        doc="Forced circulation crystallizer volume-based cost exponent (Yusuf et al., 2019)",
        units=pyo.units.dimensionless,
    )

    blk.heat_exchanger_capital_factor = pyo.Var(
        initialize=420,
        units=pyo.units.USD_2018 / (pyo.units.meter**2),
        doc="Heat exchanger cost per area",  # TODO: need reference and costing year
    )

    blk.heat_exchanger_endplates_capital_factor = pyo.Var(
        initialize=1020,
        units=pyo.units.USD_2018,
        doc="Heat exchanger endplates cost per area",  # TODO: need reference and costing year
    )

    blk.heat_exchanger_endplates_capital_basis = pyo.Var(
        initialize=10,
        units=pyo.units.meter**2,
        doc="Heat exchanger endplates cost per area basis",  # TODO: need reference and costing year
    )

    blk.heat_exchanger_endplates_capital_exponent = pyo.Var(
        initialize=0.6,
        units=pyo.units.dimensionless,
        doc="Heat exchanger endplates cost per area basis",  # TODO: need reference and costing year
    )


@register_costing_parameter_block(
    build_rule=build_recovered_nacl_cost_param_block,
    parameter_block_name="nacl_recovered",
)
@register_costing_parameter_block(
    build_rule=build_steam_cost_param_block,
    parameter_block_name="steam",
)
@register_costing_parameter_block(
    build_rule=build_multi_effect_crystallizer_cost_param_block,
    parameter_block_name="multi_effect_crystallizer",
)
def cost_multi_effect_crystallizer(
    blk, cost_type=MultiEffectCrystallizerCostType.mass_basis
):
    """
    Function for costing the forced circulation crystallizer by the mass flow of produced crystals.
    The operating cost model assumes that heat is supplied via condensation of saturated steam (see Dutta et al.)

    Args:
        cost_type: Option for crystallizer cost function type - volume or mass basis
    """
    global costing_package

    costing_package = blk.costing_package
    make_capital_cost_var(blk)
    costing_package.add_cost_factor(blk, "TIC")

    if cost_type == MultiEffectCrystallizerCostType.mass_basis:
        effect_costing_method = cost_crystallizer_effect_by_crystal_mass
    elif cost_type == MultiEffectCrystallizerCostType.volume_basis:
        effect_costing_method = cost_crystallizer_effect_by_volume
    else:
        raise ConfigurationError(
            f"{blk.unit_model.name} received invalid argument for cost_type:"
            f" {cost_type}. Argument must be a member of the MultiEffectCrystallizerCostType Enum."
        )

    total_capex_expr = 0

    for effect_number, eff in blk.unit_model.effects.items():

        effect_capex_expr = 0

        effect_capex_var = pyo.Var(
            initialize=1e5,
            units=costing_package.base_currency,
            doc=f"Capital cost effect {effect_number}",
        )

        # add capital of crystallizer
        effect_capex_expr += effect_costing_method(eff.effect)
        # add capital of heat exchangers
        effect_capex_expr += cost_crystallizer_heat_exchanger(eff.effect)

        effect_capex_constr = pyo.Constraint(
            expr=effect_capex_var == effect_capex_expr,
            doc=f"Constraint for capital cost of effect {effect_number}.",
        )
        blk.add_component(f"capital_cost_effect_{effect_number}", effect_capex_var)
        blk.add_component(
            f"capital_cost_effect_{effect_number}_constraint", effect_capex_constr
        )
        total_capex_expr += effect_capex_expr
        _cost_effect_flows(eff.effect, effect_number)

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost == blk.cost_factor * total_capex_expr
    )


def cost_crystallizer_heat_exchanger(effect):

    capital_cost_hx_effect = 0

    capital_cost_hx = pyo.units.convert(
        costing_package.multi_effect_crystallizer.heat_exchanger_capital_factor
        * effect.heat_exchanger_area,
        to_units=costing_package.base_currency,
    )

    capital_cost_hx_effect += capital_cost_hx

    dimensionless_hx_area = pyo.units.convert(
        (
            effect.heat_exchanger_area
            / costing_package.multi_effect_crystallizer.heat_exchanger_endplates_capital_basis
        ),
        to_units=pyo.units.dimensionless,
    )

    capital_cost_hx_endplates = pyo.units.convert(
        costing_package.multi_effect_crystallizer.heat_exchanger_endplates_capital_factor
        * (dimensionless_hx_area)
        ** costing_package.multi_effect_crystallizer.heat_exchanger_endplates_capital_exponent,
        to_units=costing_package.base_currency,
    )

    capital_cost_hx_effect += capital_cost_hx_endplates

    return capital_cost_hx_effect


def cost_crystallizer_effect_by_crystal_mass(effect):
    """
    Mass-based capital cost for forced circulation crystallizer
    """

    capital_cost_effect = pyo.units.convert(
        (
            costing_package.multi_effect_crystallizer.iec_percent
            * costing_package.multi_effect_crystallizer.fob_unit_cost
            * (
                sum(
                    effect.properties_solids[0].flow_mass_phase_comp["Sol", j]
                    for j in effect.config.property_package.solute_set
                )
                / costing_package.multi_effect_crystallizer.ref_capacity
            )
            ** costing_package.multi_effect_crystallizer.ref_exponent
        ),
        to_units=costing_package.base_currency,
    )

    return capital_cost_effect


def cost_crystallizer_effect_by_volume(effect):
    """
    Volume-based capital cost for forced circulation crystallizer
    """

    capital_cost_effect = pyo.units.convert(
        (
            costing_package.multi_effect_crystallizer.volume_cost
            * (
                (
                    pyo.units.convert(
                        effect.volume_suspension
                        * (effect.height_crystallizer / effect.height_slurry),
                        to_units=(pyo.units.ft) ** 3,
                    )
                )
                / pyo.units.ft**3
            )
            ** costing_package.multi_effect_crystallizer.vol_basis_exponent
        ),
        to_units=costing_package.base_currency,
    )

    return capital_cost_effect


def _cost_effect_flows(effect, effect_number):

    costing_package.cost_flow(
        pyo.units.convert(
            (
                effect.magma_circulation_flow_vol
                * effect.dens_mass_slurry
                * Constants.acceleration_gravity
                * effect.pump_head_height
                / effect.efficiency_pump
            ),
            to_units=pyo.units.kW,
        ),
        "electricity",
    )

    costing_package.cost_flow(
        effect.properties_solids[0].flow_mass_phase_comp["Sol", "NaCl"],
        "NaCl_recovered",
    )

    if effect_number == 1:
        effect.steam_flow = pyo.units.convert(
            (effect.work_mechanical[0] / _compute_steam_specific_energy(effect)),
            to_units=pyo.units.m**3 / pyo.units.s,
        )
        costing_package.cost_flow(
            effect.steam_flow,
            "steam",
        )
        # costing_package.cost_flow(effect.heating_steam[0].flow_vol_phase["Vap"],
        #     "steam",
        # )


def _compute_steam_specific_energy(effect):
    """
    Function for computing saturated steam specific energy for thermal heating estimation.

    Args:
        pressure_sat:   Steam gauge pressure in bar

    Out:
        Steam thermal capacity (latent heat of condensation * density) in kJ/m3
    """
    # TODO: add specific volume to property package?

    # pressure_sat = effect.steam_pressure
    # pressure_sat = effect.heating_steam[0].pressure_sat
    # 1. Compute saturation temperature of steam: computed from El-Dessouky expression
    # tsat_constants = [
    #     42.6776 * pyo.units.K,
    #     -3892.7 * pyo.units.K,
    #     1000 * pyo.units.kPa,
    #     -9.48654 * pyo.units.dimensionless,
    # ]
    # psat = (
    #     pyo.units.convert(pressure_sat, to_units=pyo.units.kPa)
    #     # + 101.325 * pyo.units.kPa
    # )
    # temperature_sat = tsat_constants[0] + tsat_constants[1] / (
    #     pyo.log(psat / tsat_constants[2]) + tsat_constants[3]
    # )

    # 2. Compute latent heat of condensation/vaporization: computed from Sharqawy expression
    # t = temperature_sat - 273.15 * pyo.units.K
    t = effect.heating_steam[0].temperature - 273.15 * pyo.units.K
    enth_mass_units = pyo.units.J / pyo.units.kg
    t_inv_units = pyo.units.K**-1
    dh_constants = [
        2.501e6 * enth_mass_units,
        -2.369e3 * enth_mass_units * t_inv_units**1,
        2.678e-1 * enth_mass_units * t_inv_units**2,
        -8.103e-3 * enth_mass_units * t_inv_units**3,
        -2.079e-5 * enth_mass_units * t_inv_units**4,
    ]
    dh_vap = (
        dh_constants[0]
        + dh_constants[1] * t
        + dh_constants[2] * t**2
        + dh_constants[3] * t**3
        + dh_constants[4] * t**4
    )
    # dh_vap = pyo.units.convert(dh_vap, to_units=pyo.units.kJ / pyo.units.kg)
    effect.dh_vap = dh_vap

    # 3. Compute specific volume: computed from Affandi expression (Eq 5)
    t_critical = 647.096 * pyo.units.K
    # t_red = temperature_sat / t_critical  # Reduced temperature
    t_red = effect.heating_steam[0].temperature / t_critical  # Reduced temperature

    sp_vol_constants = [
        -7.75883 * pyo.units.dimensionless,
        3.23753 * pyo.units.dimensionless,
        2.05755 * pyo.units.dimensionless,
        -0.06052 * pyo.units.dimensionless,
        0.00529 * pyo.units.dimensionless,
    ]
    log_sp_vol = (
        sp_vol_constants[0]
        + sp_vol_constants[1] * (pyo.log(1 / t_red)) ** 0.4
        + sp_vol_constants[2] / (t_red**2)
        + sp_vol_constants[3] / (t_red**4)
        + sp_vol_constants[4] / (t_red**5)
    )
    sp_vol = pyo.exp(log_sp_vol) * pyo.units.m**3 / pyo.units.kg

    # 4. Return specific energy: density * latent heat
    # return dh_vap / sp_vol
    return pyo.units.convert(
        effect.heating_steam[0].dh_vap_mass / sp_vol,
        to_units=pyo.units.kJ / pyo.units.m**3,
    )
