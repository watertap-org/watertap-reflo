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
"""
This module contains the base class for interacting with WaterTAP data files
with zero-order model parameter data.
"""
import os
from copy import deepcopy
from watertap.core.wt_database import Database


class REFLODatabase(Database):
    """
    WaterTAP Database class.

    Used to instantiate an instance of a database for loading parameters
    associated with zero-order models in WaterTap.

    Args:
        dbpath - (optional) path to database folder containing yaml files

    Returns:
        an instance of a Database object linked to the provided database
    """

    def __init__(self, dbpath=None):
        self._cached_files = {}

        if dbpath is None:
            self._dbpath = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "..",
                "data",
                "technoeconomic",
            )
        else:
            self._dbpath = dbpath

            # Confirm valid path
            if not os.path.isdir(self._dbpath):
                raise OSError(
                    f"Could not find requested path {self._dbpath}. Please "
                    f"check that this path exists."
                )

        # Create placeholder _component_list attribute
        self._component_list = None
