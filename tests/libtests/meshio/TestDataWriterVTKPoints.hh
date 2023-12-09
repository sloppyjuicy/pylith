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
 * @file tests/libtests/meshio/TestDataWriterVTKPoints.hh
 *
 * @brief C++ TestDataWriterVTKPoints object
 *
 * C++ unit testing for DataWriterVTKPoints.
 */

#if !defined(pylith_meshio_testdatawritervtkpoints_hh)
#define pylith_meshio_testdatawritervtkpoints_hh

#include "TestDataWriterVTK.hh" // ISA TestDataWriterVTK
#include "TestDataWriterPoints.hh" // ISA TestDataWriterPoints

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field

/// Namespace for pylith package
namespace pylith {
    namespace meshio {
        class TestDataWriterVTKPoints;
        class TestDataWriterVTKPoints_Data;
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterVTKPoints : public TestDataWriterVTK, public TestDataWriterPoints {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    TestDataWriterVTKPoints(TestDataWriterVTKPoints_Data* data);

    /// Destructor.
    ~TestDataWriterVTKPoints(void);

    /// Test openTimeStep() and closeTimeStep()
    void testTimeStep(void);

    /// Test writeVertexField.
    void testWriteVertexField(void);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
protected:

    /** Get test data.
     *
     * @returns Test data.
     */
    TestDataWriterPoints_Data* _getData(void);

    // PROTECTED MEMBDERS /////////////////////////////////////////////////
protected:

    TestDataWriterVTKPoints_Data* _data; ///< Data for testing.

}; // class TestDataWriterVTKPoints

// ------------------------------------------------------------------------------------------------
class pylith::meshio::TestDataWriterVTKPoints_Data : public TestDataWriterVTK_Data, public TestDataWriterPoints_Data {};

#endif // pylith_meshio_testdatawritervtkpoints_hh

// End of file
