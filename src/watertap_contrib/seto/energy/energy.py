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
