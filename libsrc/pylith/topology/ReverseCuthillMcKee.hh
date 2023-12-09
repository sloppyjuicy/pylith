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
 * @file libsrc/topology/ReverseCuthillMcKee.hh
 *
 * @brief Interface to PETSc reverse Cuthill-McKee reordering.
 */

#if !defined(pylith_topology_reversecuthillmckee_hh)
#define pylith_topology_reversecuthillmckee_hh

// Include directives ---------------------------------------------------
#include "topologyfwd.hh" // forward declarations

// ReverseCuthillMcKee --------------------------------------------------
/// Interface to PETSc reverse Cuthill-McKee reordering.
class pylith::topology::ReverseCuthillMcKee { // ReverseCuthillMcKee
    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /** Reorder vertices and cells of mesh using PETSc routines
     * implementing reverse Cuthill-McKee algorithm.
     *
     * @param mesh PyLith finite-element mesh.
     */
    static
    void reorder(topology::Mesh* mesh);

}; // ReverseCuthillMcKee

#endif // pylith_topology_reversecuthillmckee_hh

// End of file
