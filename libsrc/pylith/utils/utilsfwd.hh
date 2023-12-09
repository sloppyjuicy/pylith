// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file libsrc/utils/utilsfwd.hh
 *
 * @brief Forward declarations for PyLith utils objects.
 *
 * Including this header file eliminates the need to use separate
 * forward declarations.
 */

#if !defined(pylith_utils_utilsfwd_hh)
#define pylith_utils_utilsfwd_hh

namespace pylith {
    namespace utils {
        class PetscOptions;
        class PetscDefaults;

        class EventLogger;
        class GenericComponent;
        class PyreComponent;

        class PylithVersion;
        class PetscVersion;
        class DependenciesVersion;

        class TestArray;

    } // utils
} // pylith

#endif // pylith_utils_utilsfwd_hh

// End of file
