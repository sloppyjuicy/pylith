# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file pylith/faults/KinSrcTimeHistory.py
#
# @brief Python object for time history data file slip time function.
#
# Factory: eq_kinematic_src

from .KinSrc import KinSrc
from .faults import KinSrcTimeHistory as ModuleKinSrc


class KinSrcTimeHistory(KinSrc, ModuleKinSrc):
    """
    Time history slip time function.

    The slip time function is given in a `TimeHistory` database.

    Implements `KinSrc`.
    """
    DOC_CONFIG = {
        "cfg": """
            [pylithapp.problem.interfaces.fault.eq_ruptures.rupture]
            origin_time = 10*year

            db_auxiliary_field = spatialdata.spatialdb.UniformDB
            db_auxiliary_field.description = Time history slip time function auxiliary field spatial database
            db_auxiliary_field.values = [initiation_time, final_slip_left_lateral, final_slip_opening]
            db_auxiliary_field.data = [0.0*s, -2.0*m, 0.0*m]

            time_history.description = Slip time function time history
            time_history.filename = slipfn.timedb
            """
    }


    import pythia.pyre.inventory

    from spatialdata.spatialdb.TimeHistory import TimeHistory
    dbTimeHistory = pythia.pyre.inventory.facility("time_history", factory=TimeHistory, family="temporal_database")
    dbTimeHistory.meta['tip'] = "Time history with normalized amplitude as a function of time."

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self, name="kinsrctimehistory"):
        """Constructor.
        """
        KinSrc.__init__(self, name)
        return

    def preinitialize(self, problem):
        """Do pre-initialization setup.
        """
        KinSrc.preinitialize(self, problem)

        ModuleKinSrc.setTimeHistoryDB(self, self.dbTimeHistory)
        return

    # PRIVATE METHODS ////////////////////////////////////////////////////

    def _createModuleObj(self):
        """Call constructor for module object for access to C++ object.
        """
        ModuleKinSrc.__init__(self)


# FACTORIES ////////////////////////////////////////////////////////////

def eq_kinematic_src():
    """Factory associated with KinSrcTimeHistory.
    """
    return KinSrcTimeHistory()


# End of file
