// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#if !defined(pylith_bc_forcedatadatatet4_hh)
#define pylith_bc_forcedatadatatet4_hh

#include "PointForceData.hh"

namespace pylith {
    namespace bc {
        class PointForceDataTet4;
    } // pylith
} // bc

class pylith::bc::PointForceDataTet4 : public PointForceData {
    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    PointForceDataTet4(void);

    /// Destructor
    ~PointForceDataTet4(void);

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    static const PylithScalar _tRef; ///< Reference time for rate of change of forces.
    static const PylithScalar _forceRate; ///< Rate of change of force.
    static const PylithScalar _tResidual; ///< Time for computing residual.

    static const int _numDOF; ///< Number of degrees of freedom at each point.
    static const int _numForceDOF; ///< Number of forces at points.
    static const int _numForcePts; ///< Number of points with forces.

    static const int _id; ///< Boundary condition identifier
    static const char* _label; ///< Label for boundary condition group

    static const int _forceDOF[]; ///< Degrees of freedom that are constrained at each point
    static const int _forcePoints[]; ///< Array of indices of points with forces.
    static const PylithScalar _forceInitial[]; ///< Forces at points.
    static const PylithScalar _residual[]; ///< Residual field.

    static const char* _meshFilename; ///< Filename for input mesh.
    static const char* _dbFilename; ///< Filename of simple spatial database.
};

#endif // pylith_bc_forcedatadatatet4_hh

// End of file
