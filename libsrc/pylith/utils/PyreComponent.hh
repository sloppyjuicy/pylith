// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/**
 * @file libsrc/utils/PyreComponent.hh
 *
 * @brief C++ object for managing Pyre component and identifier information.
 */

#if !defined(pylith_utils_pyrecomponent_hh)
#define pylith_utils_pyrecomponent_hh

// Include directives ---------------------------------------------------
#include "utilsfwd.hh" // forward declarations

#include <string> // HASA std::string

// PyreComponent ----------------------------------------------------------
/** @brief C++ object for managing Pyre journals.
 *
 * Provides component name to facilitate consistent names for Pyre journals and
 * component identifier for use in error and debugging messages.
 */
class pylith::utils::PyreComponent {
    friend class TestPyreComponent; // unit testing

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /// Constructor
    PyreComponent(void);

    /// Destructor
    ~PyreComponent(void);

    /** Set name of component.
     *
     * @param value Name of component.
     */
    void setName(const char* value);

    /** Get name of component.
     *
     * @returns Name of component.
     */
    const char* getName(void) const;

    /** Set component identifier (identifies object in component hierarchy).
     *
     * @param value Component identifier.
     */
    void setIdentifier(const char* value);

    /** Get component identifier (identifies object in component hierarchy).
     *
     * @returns Component identifier.
     */
    const char* getIdentifier(void) const;

    // PRIVATE METHODS //////////////////////////////////////////////////////
private:

    std::string _name; ///< Name of component (used in journals).
    std::string _identifier; ///< Identifier for object in component hierarchy.

    // PRIVATE METHODS //////////////////////////////////////////////////////
private:

    PyreComponent(const PyreComponent&); ///< Not implemented
    const PyreComponent& operator=(const PyreComponent&); ///< Not implemented

}; // PyreComponent

#endif // pylith_utils_pyrecomponent_hh

// End of file
