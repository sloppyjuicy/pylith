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
 * @file tests/libtests/meshio/TestMeshIOCubit.hh
 *
 * @brief C++ TestMeshIOCubit object
 *
 * C++ unit testing for MeshIOCubit.
 */

#if !defined(pylith_meshio_testmeshiocubit_hh)
#define pylith_meshio_testmeshiocubit_hh

#include "TestMeshIO.hh"

namespace pylith {
    namespace meshio {
        class TestMeshIOCubit;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestMeshIOCubit : public TestMeshIO {
    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Constructor.
    TestMeshIOCubit(TestMeshIO_Data* data);

    /// Destructor.
    ~TestMeshIOCubit(void);

    /// Test filename()
    void testFilename(void);

    /// Test read().
    void testRead(void);

    // PROTECTED METHODS ////////////////////////////////////////////////
protected:

    MeshIOCubit* _io; ///< Test subject.

}; // class TestMeshIOCubit

#endif // pylith_meshio_testmeshiocubit_hh

// End of file
