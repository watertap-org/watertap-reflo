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
    Param,
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
from watertap_contrib.reflo.core import (
    SolarEnergyBaseData,
    SolarModelType,
    SolarSurrogateType,
)
from watertap_contrib.reflo.costing.solar.pv_battery import cost_pv_battery

__author__ = "Kurban Sitterley"


@declare_process_block_class("PVBatterySurrogate")
class PVBatterySurrogateData(SolarEnergyBaseData):
    """
    Surrogate model for PV+Battery System.
    """

    CONFIG = SolarEnergyBaseData.CONFIG()

    def build(self):
        super().build()

        self._tech_type = "PV_battery"

        if not self.config.solar_model_type == SolarModelType.surrogate:
            err_msg = "The PV+Battery surrogate model can only be used with the surrogate model type."
            raise ConfigurationError(err_msg)

        self.del_component(self.heat)

        self.single_inverter_capacity = Param(
            initialize=2507.2,
            units=pyunits.kW,
            doc="DC capacity of a single inverter",
        )

        self.DC_to_AC_ratio = Param(
            initialize=1.2,
            units=pyunits.dimensionless,
            doc="DC to AC ratio for PV system",
        )

        if self.config.surrogate_model_file is not None:
            self.surrogate_model_file = self.config.surrogate_model_file
            self.load_surrogate()
        elif self.config.surrogate_model_type == SolarSurrogateType.rbf:
            self.create_rbf_surrogate()
        elif self.config.surrogate_model_type == SolarSurrogateType.polynomial:
            self.create_polynomial_surrogate()
        else:
            err_msg = "The PV+Battery surrogate model requires either an existing surrogate model file or a valid surrogate type."
            raise ConfigurationError(err_msg)

        if "battery_power" not in self.input_labels:

            if "battery_power" in self.data.columns:
                self.battery_power = Param(
                    initialize=self.data["battery_power"].mean(),
                    units=pyunits.kW,
                    mutable=True,
                    doc="Battery power rating",
                )

            else:
                self.battery_power = Param(
                    initialize=100,
                    units=pyunits.kW,
                    mutable=True,
                    doc="Battery power rating",
                )

        if "hours_storage" not in self.input_labels:

            if "hours_storage" in self.data.columns:
                self.hours_storage = Param(
                    initialize=self.data["hours_storage"].mean(),
                    units=pyunits.hour,
                    mutable=True,
                    doc="Hours of storage",
                )

            else:
                self.hours_storage = Param(
                    initialize=12,
                    units=pyunits.hour,
                    mutable=True,
                    doc="Hours of storage",
                )

        if self.config.scale_training_data:
            self.electricity_annual = Expression(
                expr=self.electricity_annual_scaled / self.electricity_annual_scaling,
                doc="Annual electricity produced by PV system",
            )
            self.land_req = Expression(
                expr=self.land_req_scaled / self.land_req_scaling,
                doc="Land required for PV system",
            )

        self.inverter_capacity = Expression(
            expr=pyunits.convert(
                self.system_capacity / self.DC_to_AC_ratio, to_units=pyunits.kW
            ),
            doc="Inverter capacity for PV system",
        )

        self.electricity_constraint = Constraint(
            expr=self.electricity
            == pyunits.convert(self.electricity_annual, to_units=pyunits.kW)
        )

    def calculate_scaling_factors(self):

        if iscale.get_scaling_factor(self.system_capacity) is None:
            sf = iscale.get_scaling_factor(self.system_capacity, default=1)
            iscale.set_scaling_factor(self.system_capacity, sf)

        if iscale.get_scaling_factor(self.battery_power) is None:
            sf = iscale.get_scaling_factor(self.battery_power, default=1)
            iscale.set_scaling_factor(self.battery_power, sf)

        if iscale.get_scaling_factor(self.hours_storage) is None:
            sf = iscale.get_scaling_factor(self.hours_storage, default=1)
            iscale.set_scaling_factor(self.hours_storage, sf)

        if iscale.get_scaling_factor(self.land_req) is None:
            sf = iscale.get_scaling_factor(self.land_req, default=1)
            iscale.set_scaling_factor(self.land_req, sf)

        if iscale.get_scaling_factor(self.electricity_annual) is None:
            sf = iscale.get_scaling_factor(self.electricity_annual, default=1)
            iscale.set_scaling_factor(self.electricity_annual, sf)

        if iscale.get_scaling_factor(self.electricity) is None:
            sf = iscale.get_scaling_factor(self.electricity, default=1)
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
                "system_capacity": [value(self.system_capacity)],
                "battery_power": [value(self.battery_power)],
                "hours_storage": [value(self.hours_storage)],
            }
        )

        if self.config.surrogate_model_type == SolarSurrogateType.rbf:
            # BUG: evaluate_surrogate will raise error for polynomial surrogates
            # until IDAES dependency is updated
            self.init_output = self.surrogate.evaluate_surrogate(self.init_data)
            init_log.info_high(
                f"Initialization Step 1: Evaluate surrogate at design size {value(self.system_capacity)} kW, battery size {value(self.battery_power)*value(self.hours_storage)} kWh"
            )

        # Create solver
        res = opt.solve(self)

        init_log.info_high(f"Initialization Step 2: Solution {idaeslog.condition(res)}")

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    @property
    def default_costing_method(self):
        return cost_pv_battery
