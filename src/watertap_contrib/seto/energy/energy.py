"""
This module contains common methods for determining the energy consumption or 
production for SETO modules.
"""

import idaes.logger as idaeslog

from pyomo.environ import Param, Var, units as pyunits

import PySAM.Pvsamv1 as pv
import PySAM.Grid as grid
import PySAM.Utilityrate5 as utilityrate
import PySAM.Singleowner as singleowner

import json
import numpy as np

# Some more inforation about this module
__author__ = "Kurban Sitterley"

# Set up logger
_log = idaeslog.getLogger(__name__)

def solar_energy(self):
    """
    Helper method for implementing electricity production from solar energy models.

    Two variables are added to the model:
        * electricity (indexed by time)

    """

    # Add electricity consumption to model
    self.electricity = Var(
        units=pyunits.kW,
        bounds=(None, None),
        doc="Electricity production of solar process",
    )

    self._fixed_perf_vars.append(self.electricity)