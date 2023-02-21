from copy import deepcopy
# check
# Import Pyomo libraries
from pyomo.environ import (
    Set,
    Var,
    Constraint,
    check_optimal_termination,
    Param,
    Suffix,
    value,
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

__author__ = "Kurban Sitterley, Abdiel Lugo"

_log = idaeslog.getLogger(__name__)


class SofteningTrainType(StrEnum):
    conventional = "conventional"
    solid_contact = "solid_contact"


class SofteningProcedureType(StrEnum):
    single_stage_lime = "single_stage_lime"
    excess_lime = "excess_lime"
    single_stage_lime_soda = "single_stage_lime_soda"
    excess_lime_soda = "excess_lime_soda"


@declare_process_block_class("ChemicalSoftening0D")
class ChemicalSoftening0DData(InitializationMixin, UnitModelBlockData):
    """
    Zero order chemical softening model
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
        "softening_train_type",
        ConfigValue(
            default=SofteningTrainType.conventional,
            domain=In(SofteningTrainType),
            description="Configuration of softening system",
        ),
    )

    CONFIG.declare(
        "softening_procedure_type",
        ConfigValue(
            default=SofteningProcedureType.single_stage_lime,
            domain=In(SofteningProcedureType),
            description="Procedure of softening system",
        ),
    )

    CONFIG.declare(
        "silica_removal",
        ConfigValue(
            default=False,
            domain=bool,
            description="Include silica removal in softening process -- assumes conditions are proper",
        ),
    )

    def build(self):
        super().build()

        ### REFERENCES ###

        #ion_set = self.config.property_package.ion_set

        # this creates blank scaling factors, which are populated later
        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["parameters"] = self.config.property_package
        tmp_dict["defined_state"] = True  # inlet block is an inlet
        self.properties_in = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of inlet", **tmp_dict
        )

        # Add outlet and waste block
        tmp_dict["defined_state"] = False  # outlet and waste block is not an inlet
        self.properties_out = self.config.property_package.state_block_class(
            self.flowsheet().config.time,
            doc="Material properties of outlet",
            **tmp_dict,
        )

        self.properties_waste = self.config.property_package.state_block_class(
            self.flowsheet().config.time, doc="Material properties of waste", **tmp_dict
        )

        # Add ports
        self.add_port(name="inlet", block=self.properties_in)
        self.add_port(name="outlet", block=self.properties_out)
        self.add_port(name="waste", block=self.properties_waste)

        prop_in = self.properties_in[0]
        prop_out = self.properties_out[0]
        prop_waste = self.properties_waste[0]
        comps = self.config.property_package.component_set

        # if "TSS" in comps:
        #     self.frac_TSS_removal = Param(
        #         initialize=0.85, doc="Default 85% removal of TSS"
        #     )
        

        removal_eff_dict = dict(
            zip(
                self.config.property_package.component_set,
                [
                    0.7 if j != "TDS" else 1e-3
                    for j in self.config.property_package.component_set
                ],
            )
        )

        self.removal_efficiency = Param(
            self.config.property_package.component_set,
            initialize=removal_eff_dict,
            # mutable=True,
            doc="Removal efficiency",
        )


        # self.lime_mw = Param(
        #     initialize=74, units=pyunits.g / pyunits.mol, doc="Molecular weight of lime"
        # )

        # self.ca_eff_target = Var(
        #     initialize=0.030, 
        #     units=pyunits.kg / pyunits.m**3
        # )

        # self.mg_eff_target = Var(
        #     initialize=0.020, 
        #     units=pyunits.kg / pyunits.m**3
        # )

        self.frac_vol_recovery = Param(
            initialize=0.99, units=pyunits.dimensionless, 
            doc="Fractional volumetric recovery of water"
        )

        # self.excess_lime = Var(initialize=100, units=pyunits.mg / pyunits.liter)

        self.retention_time_mixer = Var(
            initialize=1,
            units=pyunits.minutes,
            bounds=(0.1, 5),
            doc="Retention time for mixer",
        )

        self.retention_time_floc = Var(
            initialize=12,
            units=pyunits.minutes,
            bounds=(10, 45),
            doc="Retention time for flocculator",
        )

        self.retention_time_sed = Var(
            initialize=150,
            units=pyunits.minutes,
            bounds=(120, 240),
            doc="Retention time for sedimentation basin",
        )

        self.retention_time_recarb = Var(
            initialize=18,
            units=pyunits.minutes,
            bounds=(15, 30),
            doc="Retention time for recarbonation",
        )

        self.sedimentation_overflow = Var(
            initialize=90,
            units=pyunits.m / pyunits.day,
            bounds=(30, 100),
            doc="Sedimentation basin overflow rate",
        )

        self.volume_mixer = Var(
            initialize=10, units=pyunits.m**3, doc="Volume of mixer"
        )

        self.volume_floc = Var(
            initialize=10, units=pyunits.m**3, doc="Volume of flocculator"
        )

        self.volume_sed = Var(
            initialize=10, units=pyunits.m**3, doc="Volume of sedimentation basin"
        )

        self.volume_recarb = Var(
            initialize=10, units=pyunits.m**3, doc="Volume of recarbonator"
        )

        self.vel_gradient_mix = Var(
            initialize=500,
            units=pyunits.s**-1,
            bounds=(300, 1000),
            doc="Velocity gradient of rapid mixer",
        )

        self.vel_gradient_floc = Var(
            initialize=30,
            units=pyunits.s**-1,
            bounds=(20, 80),
            doc="Velocity gradient of flocculator",
        )

        self.vel_gradient_sed = Var(
            initialize=15,
            units=pyunits.s**-1,
            bounds=(10, 50),
            doc="Velocity gradient of sedimentation",  # ???
        )

        @self.Constraint(doc="Recovery")
        def eq_recovery(b):
            return prop_in.flow_vol == prop_out.flow_vol + prop_waste.flow_vol

        @self.Constraint(doc="Waste flow")
        def eq_waste_flow(b):
            return prop_out.flow_vol == prop_in.flow_vol * b.frac_vol_recovery

        @self.Constraint(comps, doc="Mass balance")
        def eq_mass_balance(b, j):
            return prop_in.flow_mass_comp[j] == prop_out.flow_mass_comp[j] + prop_waste.flow_mass_comp[j]

        # if "TSS" in comps:
        #     @self.Constraint(doc="TSS removal")
        #     def eq_effluent_tss(b):
        #         return prop_out.conc_mass_comp["TSS"] == prop_in.conc_mass_comp["TSS"] * (1 - b.frac_TSS_removal)

        @self.Constraint(comps,doc="Component Removal")
        def eq_component_removal(b,j):
                return prop_waste.conc_mass_comp[j] == prop_in.conc_mass_comp[j] * (b.removal_efficiency[j])

        # @self.Constraint(doc="Ca in effluent")
        # def eq_effluent_ca(b):
        #     return prop_out.conc_mass_caco3_comp["Ca_2+"] == b.ca_eff_target
        #     # return prop_out.conc_mass_comp["Ca_2+"] == b.ca_eff_target

        # if (
        #     self.config.softening_procedure_type
        #     is SofteningProcedureType.single_stage_lime
        # ):

            # if not self.config.silica_removal:

            #     @self.Constraint(doc="Mg in effluent")
            #     def eq_effluent_mg(b):
            #         return (
            #             prop_out.conc_mass_caco3_comp["Mg_2+"] == prop_in.conc_mass_caco3_comp["Mg_2+"]
            #             # prop_out.conc_mass_comp["Mg_2+"] == prop_in.conc_mass_comp["Mg_2+"]
            #         )

            # @self.Constraint(doc="Total hardness in effluent")
            # def eq_effluent_total_hardness(b):
            #     return prop_out.total_hardness == (
            #         prop_out.conc_mass_caco3_comp["Ca_2+"]
            #         + (prop_in.alkalinity - prop_in.conc_mass_caco3_comp["Ca_2+"])
            #     ) + (
            #         prop_in.total_hardness
            #         - (
            #             prop_in.conc_mass_caco3_comp["Ca_2+"]
            #             + (prop_in.alkalinity - prop_in.conc_mass_caco3_comp["Ca_2+"])
            #         )
            #     )

            # @self.Constraint(doc="Alkalinity in effluent")
            # def eq_effluent_alkalinity(b):
            #     return (
            #         prop_out.alkalinity
            #         == prop_in.alkalinity
            #         - prop_in.conc_mass_caco3_comp["Ca_2+"]
            #         # + b.ca_eff_target
            #         + prop_out.conc_mass_caco3_comp["Ca_2+"]
            #     )

            # if "TDS" in comps:
            #     if "SiO2" in comps:

            #         @self.Constraint(["TDS"], doc="TDS in effluent")
            #         def eq_effluent_tds(b, j):
            #             return prop_out.conc_mass_comp[
            #                  "TDS"
            #             ] == prop_in.conc_mass_comp["TDS"] - (
            #                 prop_in.conc_mass_comp[ "Ca_2+"]
            #                 - prop_out.conc_mass_comp["Ca_2+"]
            #             ) - (
            #                 prop_in.conc_mass_comp[ "SiO2"]
            #                 - prop_out.conc_mass_comp["SiO2"]
            #             )

            #     else:

            #         @self.Constraint(["TDS"], doc="TDS in effluent")
            #         def eq_effluent_tds(b, j):
            #             return prop_out.conc_mass_comp[
            #                 "TDS"
            #             ] == prop_in.conc_mass_comp[ "TDS"] - (
            #                 prop_in.conc_mass_comp[ "Ca_2+"]
            #                 - prop_out.conc_mass_comp[ "Ca_2+"]
            #             )



            # if "SiO2" in comps and self.config.silica_removal:
                
            #     self.mg_add_ratio = Param(
            #         initialize=9.1,
            #         units=pyunits.dimensionless,
            #         doc="Ratio of MgO: SiO2 to add for SiO2 removal"
            #     )

            #     self.mg_add = Var(
            #         initialize=10,
            #         units=pyunits.kg/pyunits.m**3,
            #         doc="Concentration MgO added for SiO2 removal"
            #     )

            #     @self.Constraint(doc="MgO addition")
            #     def eq_mg_add(b):
            #         return b.mg_add == b.mg_add_ratio * prop_in.conc_mass_comp["SiO2"]

            #     @self.Constraint(["Mg_2+"], doc="Mg in effluent")
            #     def eq_effluent_mg(b, j):
            #         return (
            #             prop_out.conc_mass_caco3_comp["Mg_2+"]
            #             == prop_in.conc_mass_caco3_comp["Mg_2+"] + b.mg_add
            #         )



        @self.Constraint(doc="Volume of mixer")
        def eq_volume_mixer(b):
            return (
                b.volume_mixer
                == b.properties_in[0].flow_vol* b.retention_time_mixer
            )

        @self.Constraint(doc="Volume of flocculator")
        def eq_volume_floc(b):
            return (
                b.volume_floc
                == b.properties_in[0].flow_vol * b.retention_time_floc
            )

        @self.Constraint(doc="Volume of sedimentation basin")
        def eq_volume_sed(b):
            return (
                b.volume_sed
                == b.properties_in[0].flow_vol* b.retention_time_sed
            )

        @self.Constraint(doc="Volume of recarbonation basin")
        def eq_volume_recarb(b):
            return (
                b.volume_recarb
                == b.properties_in[0].flow_vol* b.retention_time_recarb
            )

        # @self.Constraint(doc="Magnesium hydroxide acidity constant")
        # def eq_acidity_const_mg_oh2(b):
        #     return b.acidity_const_mg_oh2 == 4470.99 / b.properties_in[0].temperature + 0.01706 * b.properties_in[0].temperature - 6.0875

        # @self.Constraint()
        # def eq_excess_lime(b):
        #     return b.excess_lime == (10 ** -b.properties_out[0].pOH / 2) * 74000

        # @self.Constraint()
        # def eq_pH_modified(b):
        #     return b.properties_out[0].pH == log(10**-b.acidity_const_mg_oh2 / 10**-b.properties_out[0].pOH)

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

        state_args_out = deepcopy(state_args)

        # for j in blk.properties_out.component_list:
        #     if j == "H2O":
        #         state_args_out["flow_vol"] = value(
        #             state_args["flow_vol"] * blk.frac_vol_recovery
        #         )
        #     else:
        #         state_args_out["conc_mass_comp"][j] = value(
        #             state_args["conc_mass_comp"][j] * (1 - blk.removal_efficiency[j])
        #         )

        blk.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_out,
        )
        init_log.info("Initialization Step 1b Complete.")

        state_args_waste = deepcopy(state_args)

        # for j in blk.properties_waste.component_list:
        #     if j == "H2O":
        #         state_args_waste["flow_vol"] = value(
        #             state_args["flow_vol"] * (1 - blk.frac_vol_recovery)
        #         )
        #     else:
        #         state_args_waste["conc_mass_comp"][j] = value(
        #             state_args["conc_mass_comp"][j] * blk.removal_efficiency[j]
        #         )

        blk.properties_waste.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            state_args=state_args_waste,
        )

        blk.state_args_out = state_args_out
        blk.state_args_waste = state_args_waste

        init_log.info("Initialization Step 1c Complete.")

        # Solve unit
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info("Initialization Step 2 {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # Release Inlet state
        blk.properties_in.release_state(flags, outlvl=outlvl)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {blk.name} failed to initialize")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

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

        pass
