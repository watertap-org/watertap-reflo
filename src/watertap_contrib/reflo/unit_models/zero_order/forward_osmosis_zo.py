#################################################################################
# WaterTAP Copyright (c) 2020-2025, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from copy import deepcopy
from pyomo.environ import (
    Var,
    Suffix,
    check_optimal_termination,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.config import is_physical_parameter_block
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.tables import (
    create_stream_table_dataframe,
)
import idaes.logger as idaeslog
from watertap.core.solvers import get_solver
from watertap_contrib.reflo.costing.units.forward_osmosis_zo import (
    cost_forward_osmosis,
)

_log = idaeslog.getLogger(__name__)
__author__ = "Zhuoran Zhang"


@declare_process_block_class("ForwardOsmosisZO")
class ForwardOsmosisZOData(UnitModelBlockData):
    """
    Forward Osmosis - ZO model
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
    **default** = False. """,
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
        "property_package_water",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for feed, distillate, brine and cooling water properties",
            doc="""Property parameter object used to define water property calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_draw_solution",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for heating and motive steam properties",
            doc="""Property parameter object used to define steasm property calculations,
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

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        """
        Specify system configurations
        """
        self.recovery_ratio = Var(
            initialize=0.3,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Recovery ratio of FO",
        )

        self.nanofiltration_recovery_ratio = Var(
            initialize=0.8,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Recovery ratio of nanofiltration",
        )

        self.regeneration_temp = Var(
            initialize=90 + 273.15,
            bounds=(0, None),
            units=pyunits.K,
            doc="Regeneration temperature of draw solution in the separator",
        )

        self.separator_temp_loss = Var(
            initialize=1,
            bounds=(0, None),
            units=pyunits.K,
            doc="Temperature loss of draw solution in the separator",
        )

        self.dp_brine = Var(
            initialize=0,
            bounds=(0, None),
            units=pyunits.Pa,
            doc="Desired pressure difference of draw solution over brine osmotic pressure",
        )

        self.heat_mixing = Var(
            initialize=75.6,
            bounds=(0, None),
            units=pyunits.MJ / pyunits.m**3,
            doc="Heat of mixing in membrane (per m3 of separated water)",
        )

        """
        Add intermediate variables
        """
        self.heat_transfer_to_brine = Var(
            initialize=50,
            bounds=(0, None),
            units=pyunits.MJ / pyunits.m**3,
            doc="Heat of mixing transferred to brine (per m3 of product water)",
        )

        self.heat_transfer_to_weak_draw = Var(
            initialize=55,
            bounds=(0, None),
            units=pyunits.MJ / pyunits.m**3,
            doc="Heat of mixing transferred to weak draw solution (per m3 of product water)",
        )

        self.delta_temp_membrane = Var(
            initialize=5,
            bounds=(0, None),
            units=pyunits.K,
            doc="Temperature difference between membrane and outlet flow due to heat of mixing",
        )

        self.membrane_temp = Var(
            initialize=20 + 273.15,
            bounds=(0, None),
            units=pyunits.K,
            doc="Membrane temperature",
        )

        """
        Add block for feed water
        """
        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package_water
        tmp_dict["defined_state"] = True

        self.feed_props = self.config.property_package_water.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of feed water",
            **tmp_dict,
        )

        """
        Add block for brine
        """
        tmp_dict["defined_state"] = False

        self.brine_props = self.config.property_package_water.state_block_class(
            self.flowsheet().config.time, doc="Material properties of brine", **tmp_dict
        )

        """
        Add block for strong draw solution entering the membrane module
        """
        tmp_dict["parameters"] = self.config.property_package_draw_solution
        self.strong_draw_props = (
            self.config.property_package_draw_solution.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of strong draw solution",
                **tmp_dict,
            )
        )

        """
        Add block for weak draw solution leaving the membrane module
        """
        self.weak_draw_props = (
            self.config.property_package_draw_solution.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of weak draw solution",
                **tmp_dict,
            )
        )

        """
        Add block for regenerated draw solution from the separator
        """
        self.reg_draw_props = (
            self.config.property_package_draw_solution.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of regenerated draw solution",
                **tmp_dict,
            )
        )

        """
        Add block for the product water (contaminated with drawsolution) from the separator
        """
        self.product_props = (
            self.config.property_package_draw_solution.state_block_class(
                self.flowsheet().config.time,
                doc="Material properties of product water",
                **tmp_dict,
            )
        )

        # Add ports
        self.add_port(name="feed", block=self.feed_props)
        self.add_port(name="brine", block=self.brine_props)
        self.add_port(name="strong_draw", block=self.strong_draw_props)
        self.add_port(name="weak_draw", block=self.weak_draw_props)
        self.add_port(name="reg_draw", block=self.reg_draw_props)
        self.add_port(name="product", block=self.product_props)

        """
        Mass balances
        """

        @self.Constraint(doc="Brine volumetric flow rate")
        def eq_brine_vol_flow(b):
            return b.brine_props[0].flow_vol_phase["Liq"] == b.feed_props[
                0
            ].flow_vol_phase["Liq"] * (1 - b.recovery_ratio)

        @self.Constraint(doc="Brine salinity")
        def eq_brine_salinity(b):
            return b.brine_props[0].conc_mass_phase_comp["Liq", "TDS"] == b.feed_props[
                0
            ].conc_mass_phase_comp["Liq", "TDS"] / (1 - b.recovery_ratio)

        @self.Constraint(doc="Product water flow rate")
        def eq_product_water_mass_flow(b):
            return (
                b.product_props[0].flow_vol_phase["Liq"]
                == b.feed_props[0].flow_vol_phase["Liq"]
                * b.recovery_ratio
                / b.nanofiltration_recovery_ratio
            )

        @self.Constraint(doc="Draw solution mass remains same in the system")
        def eq_draw_sol_mass_balance(b):
            return (
                b.strong_draw_props[0].flow_mass_phase_comp["Liq", "DrawSolution"]
                == b.weak_draw_props[0].flow_mass_phase_comp["Liq", "DrawSolution"]
            )

        @self.Constraint(doc="Draw solution dilution")
        def eq_water_mass_in_weak_draw(b):
            return (
                b.weak_draw_props[0].flow_mass_phase_comp["Liq", "H2O"]
                == b.strong_draw_props[0].flow_mass_phase_comp["Liq", "H2O"]
                + b.feed_props[0].flow_mass_phase_comp["Liq", "H2O"] * b.recovery_ratio
            )

        """
        Energy balances
        """

        @self.Constraint(doc="Brine temperature")
        def eq_brine_temp(b):
            return (
                b.brine_props[0].temperature == b.membrane_temp + b.delta_temp_membrane
            )

        @self.Constraint(doc="Weak draw solution temperature")
        def eq_weak_draw_temp(b):
            return (
                b.weak_draw_props[0].temperature
                == b.membrane_temp + b.delta_temp_membrane
            )

        @self.Constraint(doc="Heat of mixing in membrane")
        def eq_heat_mixing(b):
            return (
                b.heat_mixing == b.heat_transfer_to_weak_draw + b.heat_transfer_to_brine
            )

        @self.Constraint(
            doc="Brine and weak solution approaching same temperature at the outlet of FO membrane"
        )
        def eq_temp_dif_membrane1(b):
            return b.delta_temp_membrane == pyunits.convert(
                b.heat_transfer_to_weak_draw
                * b.product_props[0].flow_vol_phase["Liq"]
                / b.weak_draw_props[0].dens_mass_phase["Liq"]
                / b.weak_draw_props[0].flow_vol_phase["Liq"]
                / b.weak_draw_props[0].cp_mass_phase["Liq"],
                to_units=pyunits.K,
            )

        @self.Constraint(
            doc="Brine and weak solution approaching same temperature at the outlet of FO membrane"
        )
        def eq_temp_dif_membrane2(b):
            return b.delta_temp_membrane == pyunits.convert(
                b.heat_transfer_to_brine
                * b.product_props[0].flow_vol_phase["Liq"]
                / b.brine_props[0].dens_mass_phase["Liq"]
                / b.brine_props[0].flow_vol_phase["Liq"]
                / b.brine_props[0].cp_mass_phase["Liq"],
                to_units=pyunits.K,
            )

        @self.Constraint(doc="Calculate membrane temperature")
        def eq_memb_temp(b):
            return b.membrane_temp == (
                b.strong_draw_props[0].dens_mass_phase["Liq"]
                * b.strong_draw_props[0].flow_vol_phase["Liq"]
                * b.strong_draw_props[0].cp_mass_phase["Liq"]
                * b.strong_draw_props[0].temperature
                + b.feed_props[0].dens_mass_phase["Liq"]
                * b.feed_props[0].flow_vol_phase["Liq"]
                * b.feed_props[0].cp_mass_phase["Liq"]
                * b.feed_props[0].temperature
            ) / (
                b.strong_draw_props[0].dens_mass_phase["Liq"]
                * b.strong_draw_props[0].flow_vol_phase["Liq"]
                * b.strong_draw_props[0].cp_mass_phase["Liq"]
                + b.feed_props[0].dens_mass_phase["Liq"]
                * b.feed_props[0].flow_vol_phase["Liq"]
                * b.feed_props[0].cp_mass_phase["Liq"]
            )

        """
        System configuration
        """

        @self.Constraint(doc="Required draw solution osmotic pressure")
        def eq_weak_draw_osm_pres(b):
            return (
                b.weak_draw_props[0].pressure_osm_phase["Liq"]
                == b.brine_props[0].pressure_osm_phase["Liq"] + b.dp_brine
            )

        @self.Constraint(
            self.config.property_package_draw_solution.component_list,
            doc="Draw solution regenerated at the same concentration as strong draw",
        )
        def eq_reg_draw_mass_flow(b, j):
            return (
                b.reg_draw_props[0].flow_mass_phase_comp["Liq", j]
                == b.strong_draw_props[0].flow_mass_phase_comp["Liq", j]
            )

        @self.Constraint(doc="Draw solution regenerated at designed temperature")
        def eq_reg_draw_temp(b):
            return (
                b.reg_draw_props[0].temperature
                == b.regeneration_temp - b.separator_temp_loss
            )

        @self.Constraint(doc="Product water separated at designed temperature")
        def eq_product_temp(b):
            return (
                b.product_props[0].temperature
                == b.regeneration_temp - b.separator_temp_loss
            )

        # Touch properties that need to be calculated
        self.reg_draw_props[0].flow_vol_phase["Liq"]
        self.reg_draw_props[0].mass_frac_phase_comp["Liq", "DrawSolution"]
        self.reg_draw_props[0].cp_mass_phase["Liq"]
        self.product_props[0].flow_vol_phase["Liq"]
        self.product_props[0].mass_frac_phase_comp["Liq", "DrawSolution"]
        self.product_props[0].cp_mass_phase["Liq"]

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
        flags = blk.feed_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args,
            hold_state=True,
        )
        init_log.info("FO Initialization Step 1 Complete.")
        # ---------------------------------------------------------------------
        # Initialize other state blocks
        if state_args is None:
            blk.state_args = state_args = {}
            state_dict = blk.feed_props[
                blk.flowsheet().config.time.first()
            ].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        # Initial guess of brine properties
        state_args_brine = deepcopy(state_args)
        for p, j in blk.brine_props.phase_component_set:
            if p == "Liq" and j == "H2O":
                state_args_brine["flow_mass_phase_comp"][(p, j)] = (
                    state_args["flow_mass_phase_comp"][(p, j)]
                    * (1 - blk.recovery_ratio)
                    * pyunits.kg
                    / pyunits.s
                )
            elif p == "Liq" and j == "TDS":
                state_args_brine["flow_mass_phase_comp"][(p, j)] = state_args[
                    "flow_mass_phase_comp"
                ][(p, j)]

        blk.brine_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_brine,
        )

        # Initialize draw solution blocks
        state_args_draw_solution = {}
        state_dict_draw_solution = blk.strong_draw_props[
            blk.flowsheet().config.time.first()
        ].define_port_members()

        for k in state_dict_draw_solution.keys():
            if state_dict_draw_solution[k].is_indexed():
                state_args_draw_solution[k] = {}
                for m in state_dict_draw_solution[k].keys():
                    state_args_draw_solution[k][m] = state_dict_draw_solution[k][
                        m
                    ].value
            else:
                state_args_draw_solution[k] = state_dict_draw_solution[k].value

        # Initial guess of weak draw properties
        state_args_weak_draw = deepcopy(state_args_draw_solution)
        for p, j in blk.weak_draw_props.phase_component_set:
            if p == "Liq" and j == "H2O":
                state_args_weak_draw["flow_mass_phase_comp"][(p, j)] = (
                    state_args["flow_mass_phase_comp"][(p, j)]
                    * blk.recovery_ratio
                    / blk.nanofiltration_recovery_ratio
                    * 2  # Approximated ratio of the weak draw flow to the feed flow
                    * pyunits.kg
                    / pyunits.s
                )
            elif p == "Liq" and j == "DrawSolution":
                state_args_weak_draw["flow_mass_phase_comp"][(p, j)] = (
                    state_args_weak_draw["flow_mass_phase_comp"][(p, "H2O")]
                )

        blk.weak_draw_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_weak_draw,
        )

        # Initial guess of strong draw properties
        state_args_strong_draw = deepcopy(state_args_draw_solution)
        for p, j in blk.strong_draw_props.phase_component_set:
            if p == "Liq" and j == "DrawSolution":
                state_args_strong_draw["flow_mass_phase_comp"][(p, j)] = (
                    state_args_weak_draw["flow_mass_phase_comp"][(p, "DrawSolution")]
                )
            elif p == "Liq" and j == "H2O":
                state_args_strong_draw["flow_mass_phase_comp"][(p, j)] = (
                    state_args_weak_draw["flow_mass_phase_comp"][(p, "DrawSolution")]
                    * 0.15  # typical draw : water in strong draw solution
                )

        blk.strong_draw_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_strong_draw,
        )

        # Initial guess of regenerated draw solution properties
        state_args_reg_draw = deepcopy(state_args_draw_solution)
        for p, j in blk.reg_draw_props.phase_component_set:
            if p == "Liq":
                state_args_reg_draw["flow_mass_phase_comp"][(p, j)] = (
                    state_args_strong_draw["flow_mass_phase_comp"][(p, j)]
                )
        state_args_reg_draw["temperature"] = 90 + 273.15

        blk.reg_draw_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_reg_draw,
        )

        # Initial guess of product properties
        state_args_product = deepcopy(state_args_draw_solution)
        for p, j in blk.product_props.phase_component_set:
            if p == "Liq" and j == "H2O":
                state_args_product["flow_mass_phase_comp"][(p, j)] = (
                    state_args["flow_mass_phase_comp"][(p, j)]
                    * (blk.recovery_ratio / blk.nanofiltration_recovery_ratio)
                    * pyunits.kg
                    / pyunits.s
                )
            elif p == "Liq" and j == "DrawSolution":
                state_args_product["flow_mass_phase_comp"][(p, j)] = (
                    state_args["flow_mass_phase_comp"][(p, "H2O")]
                    * (blk.recovery_ratio / blk.nanofiltration_recovery_ratio)
                    * 0.01  # Typical mass fraction of draw in product
                    * pyunits.kg
                    / pyunits.s
                )

        blk.product_props.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_product,
        )

        # Check degree of freedom
        assert degrees_of_freedom(blk) == 0

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info("FO Initialization Step 2 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.feed_props.release_state(flags, outlvl=outlvl)
        init_log.info("FO Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {blk.name} failed to initialize")

    def unfix_and_fix_freedom(self, mass_frac_strong_draw, mass_frac_product):

        # Unfix the mass flow rate of draw solution to apply the specified mass fraction
        self.strong_draw_props[0].flow_mass_phase_comp["Liq", "H2O"].unfix()
        self.strong_draw_props[0].mass_frac_phase_comp["Liq", "DrawSolution"].fix(
            mass_frac_strong_draw
        )

        self.product_props[0].flow_mass_phase_comp["Liq", "DrawSolution"].unfix()
        self.product_props[0].mass_frac_phase_comp["Liq", "DrawSolution"].fix(
            mass_frac_product
        )

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        if iscale.get_scaling_factor(self.recovery_ratio) is None:
            iscale.set_scaling_factor(self.recovery_ratio, 1e0)

        if iscale.get_scaling_factor(self.regeneration_temp) is None:
            iscale.set_scaling_factor(self.regeneration_temp, 1e-2)

        if iscale.get_scaling_factor(self.separator_temp_loss) is None:
            iscale.set_scaling_factor(self.separator_temp_loss, 1e0)

        if iscale.get_scaling_factor(self.heat_mixing) is None:
            iscale.set_scaling_factor(self.heat_mixing, 1e-2)

        if iscale.get_scaling_factor(self.heat_transfer_to_brine) is None:
            iscale.set_scaling_factor(self.heat_transfer_to_brine, 1e-2)

        if iscale.get_scaling_factor(self.heat_transfer_to_weak_draw) is None:
            iscale.set_scaling_factor(self.heat_transfer_to_weak_draw, 1e-2)

        if iscale.get_scaling_factor(self.dp_brine) is None:
            iscale.set_scaling_factor(self.dp_brine, 1e-6)

        if iscale.get_scaling_factor(self.delta_temp_membrane) is None:
            iscale.set_scaling_factor(self.delta_temp_membrane, 1e0)

        if iscale.get_scaling_factor(self.membrane_temp) is None:
            iscale.set_scaling_factor(self.membrane_temp, 1e-2)

        # Transforming constraints
        sf = iscale.get_scaling_factor(self.feed_props[0].flow_vol_phase["Liq"])
        iscale.constraint_scaling_transform(self.eq_brine_vol_flow, sf)

        sf = iscale.get_scaling_factor(
            self.feed_props[0].conc_mass_phase_comp["Liq", "TDS"]
        )
        iscale.constraint_scaling_transform(self.eq_brine_salinity, sf)

        sf = iscale.get_scaling_factor(self.feed_props[0].temperature)
        iscale.constraint_scaling_transform(self.eq_brine_temp, sf)

        sf = iscale.get_scaling_factor(self.weak_draw_props[0].temperature)
        iscale.constraint_scaling_transform(self.eq_weak_draw_temp, sf)

        sf = iscale.get_scaling_factor(
            self.strong_draw_props[0].flow_mass_phase_comp["Liq", "DrawSolution"]
        )
        iscale.constraint_scaling_transform(self.eq_draw_sol_mass_balance, sf)

        sf = iscale.get_scaling_factor(
            self.strong_draw_props[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        iscale.constraint_scaling_transform(self.eq_water_mass_in_weak_draw, sf)

        sf = iscale.get_scaling_factor(
            self.weak_draw_props[0].pressure_osm_phase["Liq"]
        )
        iscale.constraint_scaling_transform(self.eq_weak_draw_osm_pres, sf)

        sf = iscale.get_scaling_factor(self.heat_mixing)
        iscale.constraint_scaling_transform(self.eq_heat_mixing, sf)

        sf = iscale.get_scaling_factor(self.delta_temp_membrane)
        iscale.constraint_scaling_transform(self.eq_temp_dif_membrane1, sf)

        sf = iscale.get_scaling_factor(self.delta_temp_membrane)
        iscale.constraint_scaling_transform(self.eq_temp_dif_membrane2, sf)

        sf = iscale.get_scaling_factor(
            self.product_props[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        iscale.constraint_scaling_transform(self.eq_product_water_mass_flow, sf)

        sf = iscale.get_scaling_factor(self.reg_draw_props[0].temperature)
        iscale.constraint_scaling_transform(self.eq_reg_draw_temp, sf)

        for ind, c in self.eq_reg_draw_mass_flow.items():
            sf = iscale.get_scaling_factor(
                self.strong_draw_props[0].flow_mass_phase_comp["Liq", ind]
            )
            iscale.constraint_scaling_transform(c, sf)

        sf = iscale.get_scaling_factor(self.product_props[0].temperature)
        iscale.constraint_scaling_transform(self.eq_product_temp, sf)

        sf = iscale.get_scaling_factor(self.membrane_temp)
        iscale.constraint_scaling_transform(self.eq_memb_temp, sf)

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Water": self.feed_props,
                "Brine": self.brine_props,
                "Strong Draw": self.strong_draw_props,
                "Weak Draw": self.weak_draw_props,
                "Regenerated Draw": self.reg_draw_props,
                "Product Water": self.product_props,
            },
            time_point=time_point,
        )

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        var_dict["Strong draw solution volumetric flow rate (m3/s)"] = (
            self.strong_draw_props[0].flow_vol_phase["Liq"]
        )
        var_dict["Weak draw solution volumetric flow rate (m3/s)"] = (
            self.weak_draw_props[0].flow_vol_phase["Liq"]
        )
        var_dict["Mass fraction of weak draw solution"] = self.weak_draw_props[
            0
        ].mass_frac_phase_comp["Liq", "DrawSolution"]
        var_dict["Brine salinity (g/L)"] = self.brine_props[0].conc_mass_phase_comp[
            "Liq", "TDS"
        ]
        var_dict["Membrane temperature (K)"] = self.membrane_temp

        return {"vars": var_dict}

    @property
    def default_costing_method(self):
        return cost_forward_osmosis
