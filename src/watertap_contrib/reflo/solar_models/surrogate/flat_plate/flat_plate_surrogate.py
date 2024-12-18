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
from idaes.core.util.exceptions import InitializationError
import idaes.logger as idaeslog

from watertap.core.solvers import get_solver
from watertap_contrib.reflo.core import SolarEnergyBaseData
from watertap_contrib.reflo.costing.solar.flat_plate import cost_flat_plate

__author__ = "Matthew Boyd, Kurban Sitterley"


@declare_process_block_class("FlatPlateSurrogate")
class FlatPlateSurrogateData(SolarEnergyBaseData):
    """
    Surrogate model for flat plate.
    """

    def build(self):
        super().build()

        self._tech_type = "flat_plate"
        self.add_surrogate_variables()
        self.get_surrogate_data()

        self.specific_heat_water = Param(
            initialize=4.181,  # defaults from SAM
            units=pyunits.kJ / (pyunits.kg * pyunits.K),
            doc="Specific heat of water",
        )

        self.dens_water = Param(
            initialize=1000,  # defaults from SAM
            units=pyunits.kg / pyunits.m**3,
            mutable=True,
            doc="Density of water",
        )

        self.temperature_cold = Param(
            initialize=293,  # defaults from SAM
            units=pyunits.degK,
            mutable=True,
            doc="Cold temperature",
        )

        self.factor_delta_T = Param(
            initialize=0.03,
            units=pyunits.degK,
            mutable=True,
            doc="Influent minus ambient temperature",
        )

        self.collector_area_per = Param(
            initialize=2.98,  # defaults from SAM
            units=pyunits.m**2,
            mutable=True,
            doc="Area for single collector",
        )

        self.FR_ta = Param(
            initialize=0.689,  # optical gain "a" in Hottel-Whillier-Bliss equation [hcoll = a - b*dT]; defaults from SAM
            units=pyunits.kilowatt / pyunits.m**2,
            mutable=True,
            doc="Product of collector heat removal factor (FR), cover transmittance (t), and shortwave absorptivity of absorber (a)",
        )

        self.FR_UL = Param(
            initialize=3.85,  # Thermal loss coeff "b" in Hottel-Whillier-Bliss equation [hcoll = a - b*dT]; defaults from SAM
            units=pyunits.kilowatt / (pyunits.m**2 * pyunits.degK),
            mutable=True,
            doc="Product of collector heat removal factor (FR) and overall heat loss coeff. of collector (UL)",
        )

        self.heat_annual = Expression(
            expr=self.heat_annual_scaled / self.heat_annual_scaling,
            doc="Annual heat generated by flat plate collector in kWh",
        )

        self.electricity_annual = Expression(
            expr=self.electricity_annual_scaled / self.electricity_annual_scaling,
            doc="Annual electricity consumed by flat plate collector in kWh",
        )

        if self.config.surrogate_model_file is not None:
            self.surrogate_file = self.config.surrogate_model_file
            self.load_surrogate()
        else:
            self.create_rbf_surrogate()

        self.heat_constraint = Constraint(
            expr=self.heat_annual
            == self.heat * pyunits.convert(1 * pyunits.year, to_units=pyunits.hour)
        )

        self.electricity_constraint = Constraint(
            expr=self.electricity_annual
            == self.electricity
            * pyunits.convert(1 * pyunits.year, to_units=pyunits.hour)
        )

        self.collector_area_total = Expression(
            expr=pyunits.convert(self.heat_load, to_units=pyunits.kilowatt)
            / (self.FR_ta - self.FR_UL * self.factor_delta_T)
        )

        self.number_collectors = Expression(
            expr=self.collector_area_total / self.collector_area_per
        )

        self.storage_volume = Expression(
            expr=pyunits.convert(
                (
                    (self.hours_storage * self.heat_load)
                    / (
                        self.specific_heat_water
                        * (
                            (self.temperature_hot + 273.15 * pyunits.degK)
                            - self.temperature_cold
                        )
                        * self.dens_water
                    )
                ),
                to_units=pyunits.m**3,
            )
        )

    def calculate_scaling_factors(self):

        if iscale.get_scaling_factor(self.hours_storage) is None:
            sf = iscale.get_scaling_factor(self.hours_storage, default=1)
            iscale.set_scaling_factor(self.hours_storage, sf)

        if iscale.get_scaling_factor(self.heat_load) is None:
            sf = iscale.get_scaling_factor(self.heat_load, default=1, warning=True)
            iscale.set_scaling_factor(self.heat_load, sf)

        if iscale.get_scaling_factor(self.temperature_hot) is None:
            sf = iscale.get_scaling_factor(
                self.temperature_hot, default=1, warning=True
            )
            iscale.set_scaling_factor(self.temperature_hot, sf)

        if iscale.get_scaling_factor(self.heat) is None:
            sf = iscale.get_scaling_factor(self.heat, default=1, warning=True)
            iscale.set_scaling_factor(self.heat, sf)

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
                "heat_load": [value(self.heat_load)],
                "hours_storage": [value(self.hours_storage)],
                "temperature_hot": [value(self.temperature_hot)],
            }
        )
        self.init_output = self.surrogate.evaluate_surrogate(self.init_data)
        self.heat_annual_scaled.set_value(self.init_output.heat_annual_scaled.values[0])
        self.electricity_annual_scaled.set_value(
            self.init_output.electricity_annual_scaled.values[0]
        )
        self.heat.set_value(value(self.heat_annual) / 8766)
        self.electricity.set_value(value(self.electricity_annual) / 8766)
        # Create solver
        res = opt.solve(self)

        init_log.info_high(f"Initialization Step 2 {idaeslog.condition(res)}")

        if not check_optimal_termination(res):
            raise InitializationError(f"Unit model {self.name} failed to initialize")

        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))

    @property
    def default_costing_method(self):
        return cost_flat_plate
