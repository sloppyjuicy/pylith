#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================
# @file tests/pytests/meshio/TestOutputTriggerTime.py
#
# @brief Unit testing of Python OutputTriggerTime object.

import unittest

from pylith.testing.UnitTestApp import TestComponent
from pylith.meshio.OutputTriggerTime import (OutputTriggerTime, output_trigger)


class TestOutputTriggerTime(TestComponent):
    """Unit testing of OutputTriggerTime object.
    """
    _class = OutputTriggerTime
    _factory = output_trigger


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestOutputTriggerTime))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file
