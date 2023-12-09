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
 * @file tests/libtests/meshio/TestMeshIOPetsc.hh
 *
 * @brief C++ TestMeshIOPetsc object
 *
 * C++ unit testing for MeshIOPetsc.
 */

#if !defined(pylith_meshio_testmeshiopetsc_hh)
#define pylith_meshio_testmeshiopetsc_hh

#include "TestMeshIO.hh"

namespace pylith {
    namespace meshio {
        class TestMeshIOPetsc;
    } // meshio
} // pylith

class pylith::meshio::TestMeshIOPetsc : public TestMeshIO {
    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Constructor.
    TestMeshIOPetsc(TestMeshIO_Data* data);

    /// Destructor.
    ~TestMeshIOPetsc(void);

    /// Test filename()
    void testFilename(void);

    /// Test read().
    void testRead(void);

    // PROTECTED METHODS ////////////////////////////////////////////////
protected:

    MeshIOPetsc* _io; ///< Test subject.

}; // class TestMeshIOPetsc

#endif // pylith_meshio_testmeshiopetsc_hh

// End of file
