// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "TestMeshIOPetsc.hh" // Implementation of class methods

#include "pylith/meshio/MeshIOPetsc.hh"
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END
#include "pylith/utils/journals.hh" // USES JournalingComponent

#include "catch2/catch_test_macros.hpp"

#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
// Constructor.
pylith::meshio::TestMeshIOPetsc::TestMeshIOPetsc(TestMeshIO_Data* data) :
    TestMeshIO(data) {
    _io = new MeshIOPetsc();assert(_io);
    _io->PyreComponent::setIdentifier("TestMeshIOPetsc");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor.
pylith::meshio::TestMeshIOPetsc::~TestMeshIOPetsc(void) {
    delete _io;_io = nullptr;
} // destructor


// ------------------------------------------------------------------------------------------------
// Test filename()
void
pylith::meshio::TestMeshIOPetsc::testFilename(void) {
    PYLITH_METHOD_BEGIN;
    assert(_io);

    const std::string& filename = "hi.txt";
    _io->setFilename(filename.c_str());
    CHECK(filename == std::string(_io->getFilename()));

    PYLITH_METHOD_END;
} // testFilename


// ------------------------------------------------------------------------------------------------
// Test read().
void
pylith::meshio::TestMeshIOPetsc::testRead(void) {
    PYLITH_METHOD_BEGIN;
    assert(_io);
    assert(_data);

    // Read mesh
    _io->setFilename(_data->filename.c_str());
    delete _mesh;_mesh = new topology::Mesh;assert(_mesh);
    _io->read(_mesh);

    pythia::journal::debug_t debug("TestMeshIOPetsc");
    if (debug.state()) {
        _mesh->view();
    } // if

    // Make sure mesh matches data
    _checkVals();

    PYLITH_METHOD_END;
} // testRead


// End of file
