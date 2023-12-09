// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

const char* pylith::meshio::DataWriterVTKDataFaultMeshTri3::_vertexFilename =
    "tri3_fault_vertex.vtk";

const char* pylith::meshio::DataWriterVTKDataFaultMeshTri3::_cellFilename =
    "tri3_fault_cell.vtk";

const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_time = 1.0;

const char* pylith::meshio::DataWriterVTKDataFaultMeshTri3::_timeFormat =
    "%3.1f";

const int pylith::meshio::DataWriterVTKDataFaultMeshTri3::_numVertices = 2;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataFaultMeshTri3::_vertexFields[4] = {
    { "pressure", topology::FieldBase::SCALAR, 1 },
    { "displacement", topology::FieldBase::VECTOR, 2 },
    { "stress", topology::FieldBase::TENSOR, 3 },
    { "other", topology::FieldBase::OTHER, 2 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_vertexFieldScalar[2*1] = {
    2.1, 3.2,
};
const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_vertexFieldVector[2*2] = {
    1.1, 2.2,
    3.3, 4.4,
};
const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_vertexFieldTensor[2*3] = {
    1.1, 1.2, 1.3,
    2.1, 2.2, 2.3,
};
const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_vertexFieldOther[2*2] = {
    1.2, 2.3,
    3.4, 4.5,
};

const int pylith::meshio::DataWriterVTKDataFaultMeshTri3::_numCells = 1;

const pylith::meshio::DataWriterData::FieldStruct
pylith::meshio::DataWriterVTKDataFaultMeshTri3::_cellFields[4] = {
    { "pressure", topology::FieldBase::SCALAR, 1 },
    { "traction", topology::FieldBase::VECTOR, 2 },
    { "stress", topology::FieldBase::TENSOR, 3 },
    { "other", topology::FieldBase::OTHER, 4 },
};
const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_cellFieldScalar[1*1] = {
    2.1,
};
const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_cellFieldVector[1*2] = {
    1.1, 2.2,
};
const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_cellFieldTensor[1*3] = {
    1.2, 2.3, 3.4,
};
const PylithScalar pylith::meshio::DataWriterVTKDataFaultMeshTri3::_cellFieldOther[1*4] = {
    1.2, 2.3, 3.4, 4.5,
};

pylith::meshio::DataWriterVTKDataFaultMeshTri3::DataWriterVTKDataFaultMeshTri3(void) { // constructor
    meshFilename = const_cast<char*>(_meshFilename);
    faultLabel = const_cast<char*>(_faultLabel);
    faultId = _faultId;

    timestepFilename = const_cast<char*>(_timestepFilename);
    vertexFilename = const_cast<char*>(_vertexFilename);
    cellFilename = const_cast<char*>(_cellFilename);

    time = _time;
    timeFormat = const_cast<char*>(_timeFormat);

    numVertices = _numVertices;
    assert(DataWriterData::numVertexFields == numVertexFields);
    vertexFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_vertexFields);
    vertexFields[0] = const_cast<PylithScalar*>(_vertexFieldScalar);
    vertexFields[1] = const_cast<PylithScalar*>(_vertexFieldVector);
    vertexFields[2] = const_cast<PylithScalar*>(_vertexFieldTensor);
    vertexFields[3] = const_cast<PylithScalar*>(_vertexFieldOther);

    numCells = _numCells;
    assert(DataWriterData::numCellFields == numCellFields);
    cellFieldsInfo = const_cast<DataWriterData::FieldStruct*>(_cellFields);
    cellFields[0] = const_cast<PylithScalar*>(_cellFieldScalar);
    cellFields[1] = const_cast<PylithScalar*>(_cellFieldVector);
    cellFields[2] = const_cast<PylithScalar*>(_cellFieldTensor);
    cellFields[3] = const_cast<PylithScalar*>(_cellFieldOther);
} // constructor


pylith::meshio::DataWriterVTKDataFaultMeshTri3::~DataWriterVTKDataFaultMeshTri3(void) {}


// End of file
