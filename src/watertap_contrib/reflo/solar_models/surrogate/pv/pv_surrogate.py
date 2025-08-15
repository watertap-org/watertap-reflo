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

import pandas as pd
from pyomo.environ import (
    Constraint,
    Expression,
    value,
    check_optimal_termination,
    units as pyunits,
)
from idaes.core import declare_process_block_class
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.logger as idaeslog

from watertap.core.solvers import get_solver
from watertap_contrib.reflo.core import SolarEnergyBaseData, SolarModelType, SolarSurrogateType
from watertap_contrib.reflo.costing.solar.pv import cost_pv

__author__ = "Zachary Binger, Matthew Boyd, Kurban Sitterley"


@declare_process_block_class("PVSurrogate")
class PVSurrogateData(SolarEnergyBaseData):
    """
    Surrogate model for PV
    """

    CONFIG = SolarEnergyBaseData.CONFIG()

    def build(self):
        super().build()
        self._tech_type = "PV"

        if not self.config.solar_model_type == SolarModelType.surrogate:
            err_msg = "The PV surrogate model can only be used with the surrogate model type."
            raise ConfigurationError(err_msg)

        # if self.config.scale_training_data:
        #     err_msg = "The PV surrogate model requires the input data not be scaled."
        #     err_msg += " Set the config argument 'scale_training_data' = False"
        #     raise ConfigurationError(err_msg)

        self.del_component(self.heat)

        if self.config.surrogate_model_file is not None:
            self.surrogate_model_file = self.config.surrogate_model_file
            self.load_surrogate()
        elif self.config.surrogate_model_type == SolarSurrogateType.rbf:
            self.create_rbf_surrogate()
        elif self.config.surrogate_model_type == SolarSurrogateType.polynomial:
            self.create_polynomial_surrogate()
        else:
            err_msg = "The PV surrogate model requires either an existing surrogate model file or a valid surrogate type."
            raise ConfigurationError(err_msg)

        if self.config.scale_training_data:
            self.annual_energy = Expression(
                expr=self.annual_energy_scaled / self.annual_energy_scaling)
            

        self.electricity_constraint = Constraint(
            expr=self.electricity
            == pyunits.convert(self.annual_energy, to_units=pyunits.kW)
            # expr=self.annual_energy
            # == pyunits.convert(self.electricity, to_units=pyunits.kWh / pyunits.year)
        )

    def calculate_scaling_factors(self):

        if iscale.get_scaling_factor(self.design_size) is None:
            sf = iscale.get_scaling_factor(self.design_size, default=1)
            iscale.set_scaling_factor(self.design_size, sf)

        if iscale.get_scaling_factor(self.annual_energy) is None:
            sf = iscale.get_scaling_factor(self.annual_energy, default=1, warning=True)
            iscale.set_scaling_factor(self.annual_energy, sf)

        if iscale.get_scaling_factor(self.electricity) is None:
            sf = iscale.get_scaling_factor(self.electricity, default=1, warning=True)
            iscale.set_scaling_factor(self.electricity, sf)

    def initialize(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        General wrapper for initialization routines

        Keyword Arguments:
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating which solver to use during
                     initialization (default = None)

        Returns: None
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")

        if solver is None:
            opt = get_solver(optarg)

        self.init_data = pd.DataFrame(
            {
                "design_size": [value(self.design_size)],
            }
        )
        self.init_output = self.surrogate.evaluate_surrogate(self.init_data)

        # Create solver
        res = opt.solve(self)

        init_log.info_high(f"Initialization Step 2 {idaeslog.condition(res)}")

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    @property
    def default_costing_method(self):
        return cost_pv
