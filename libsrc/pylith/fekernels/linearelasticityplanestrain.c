/* -*- C -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University of Chicago
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2015 University of California, Davis
 *
 * See COPYING for license information.
 *
 * ----------------------------------------------------------------------
 */

#include <portinfo>

#include "pylith/fekernels/linearelasticityplanestrain.h"

#include "pylith/fekernels/elasticity.h" /* USES Elasticity_f0_inertia, Elasticity_g0_bodyforce */

/* ====================================================================== 
 * Kernels for isotropic, linear elatsicity plane strain.
 * ====================================================================== 
 */


/* ---------------------------------------------------------------------- */
/** f0 function for isotropic linear elasticity plane strain.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_f0(const PylithInt dim,
							 const PylithInt numS,
							 const PylithInt numA,
							 const PylithInt sOff[],
							 const PylithInt sOff_x[],
							 const PylithScalar s[],
							 const PylithScalar s_t[],
							 const PylithScalar s_x[],
							 const PylithInt aOff[],
							 const PylithInt aOff_x[],
							 const PylithScalar a[],
							 const PylithScalar a_t[],
							 const PylithScalar a_x[],
							 const PylithReal t,
							 const PylithScalar x[],
							 PylithScalar f0[])
{ /* IsotropicLinearElasticityPlaneStrain_f0 */
  const PylithInt _dim = 2;

  const PylithInt _numS = 2; /* Number passed on to f0_inertia. */

  const PylithInt _numA = 1; /* Number passed on to f0_inertia. */
  const PylithInt i_density = 0;

  assert(_dim == dim);
  assert(2 == numS);
  assert(3 == numA || 4 == numA);
  assert(aOff);
  assert(aOff_x);

  pylith_fekernels_Elasticity_f0_inertia(_dim, _numS, _numA, sOff, sOff_x, s, s_t, s_x, &aOff[i_density], &aOff_x[i_density], a, a_t, a_x, t, x, f0);
} /* IsotropicLinearElasticityPlaneStrain_f0 */


/* ---------------------------------------------------------------------- */
/** g0 function for isotropic linear elasticity plane strain.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g0(const PylithInt dim,
							 const PylithInt numS,
							 const PylithInt numA,
							 const PylithInt sOff[],
							 const PylithInt sOff_x[],
							 const PylithScalar s[],
							 const PylithScalar s_t[],
							 const PylithScalar s_x[],
							 const PylithInt aOff[],
							 const PylithInt aOff_x[],
							 const PylithScalar a[],
							 const PylithScalar a_t[],
							 const PylithScalar a_x[],
							 const PylithReal t,
							 const PylithScalar x[],
							 PylithScalar g0[])
{ /* IsotropicLinearElasticityPlaneStrain_g0 */
  const PylithInt _dim = 2;

  const PylithInt _numS = 2;

  const PylithInt _numA = 1;
  const PylithInt i_bodyforce = 3;

  assert(_dim == dim);
  assert(2 == numS);
  assert(4 == numA);
  assert(aOff);
  assert(aOff_x);

  pylith_fekernels_Elasticity_g0_bodyforce(_dim, _numS, _numA, sOff, sOff_x, s, s_t, s_x, &aOff[i_bodyforce], &aOff_x[i_bodyforce], a, a_t, a_x, t, x, g0);
} /* IsotropicLinearElasticityPlaneStrain_g0 */


/* ---------------------------------------------------------------------- */
/** g1 function for isotropic linear elasticity plane strain WITHOUT initial stress and strain.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1(const PylithInt dim,
							 const PylithInt numS,
							 const PylithInt numA,
							 const PylithInt sOff[],
							 const PylithInt sOff_x[],
							 const PylithScalar s[],
							 const PylithScalar s_t[],
							 const PylithScalar s_x[],
							 const PylithInt aOff[],
							 const PylithInt aOff_x[],
							 const PylithScalar a[],
							 const PylithScalar a_t[],
							 const PylithScalar a_x[],
							 const PylithReal t,
							 const PylithScalar x[],
							 PylithScalar g1[])
{ /* IsotropicLinearElasticityPlaneStrain_g1 */
  const PylithInt _dim = 2;

  const PylithInt _numS = 1; /* Pass displacement field. */
  const PylithInt i_disp = 0;

  const PylithInt i_shearModulus = 1;
  const PylithInt i_bulkModulus = 2;

  const PylithInt numAMean = 1; /* Pass bulk modulus. */
  const PylithInt aOffMean[1] = { aOff[i_bulkModulus] };
  const PylithInt aOffMean_x[1] = { aOff_x[i_bulkModulus] };

  const PylithInt numADev = 1; /* Pass shear modulus. */
  const PylithInt aOffDev[1] = { aOff[i_shearModulus] };
  const PylithInt aOffDev_x[1] = { aOff_x[i_shearModulus] };

  assert(_dim == dim);
  assert(2 == numS);
  assert(3 == numA || 4 == numA);
  assert(sOff);
  assert(sOff_x);
  assert(aOff);
  assert(aOff_x);

  pylith_fekernels_IsotropicLinearElasticityPlaneStrain_meanStress(_dim, _numS, numAMean, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x, aOffMean, aOffMean_x, a, a_t, a_x, t, x, g1);
  pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress(_dim, _numS, numADev, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x, aOffDev, aOffDev_x, a, a_t, a_x, t, x, g1);
} /* IsotropicLinearElasticityPlaneStrain_g1 */


/* ---------------------------------------------------------------------- */
/** g1 function for isotropic linear elasticity plane strain with initial stress and strain.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1_initstate(const PylithInt dim,
								   const PylithInt numS,
								   const PylithInt numA,
								   const PylithInt sOff[],
								   const PylithInt sOff_x[],
								   const PylithScalar s[],
								   const PylithScalar s_t[],
								   const PylithScalar s_x[],
								   const PylithInt aOff[],
								   const PylithInt aOff_x[],
								   const PylithScalar a[],
								   const PylithScalar a_t[],
								   const PylithScalar a_x[],
								   const PylithReal t,
								   const PylithScalar x[],
								   PylithScalar g1[])
{ /* IsotropicLinearElasticityPlaneStrain_g1_initstate */
  const PylithInt _dim = 2;

  const PylithInt _numS = 1; /* Pass displacement field. */
  const PylithInt i_disp = 0;

  const PylithInt i_shearModulus = 1;
  const PylithInt i_bulkModulus = 2;
  const PylithInt i_istress = numA-2;
  const PylithInt i_istrain = numA-1;

  const PylithInt numAMean = 3; /* Pass shear modulus, initial stress, and initial strain. */
  const PylithInt aOffMean[3] = { aOff[i_bulkModulus], aOff[i_istress], aOff[i_istrain] };
  const PylithInt aOffMean_x[3] = { aOff_x[i_bulkModulus], aOff_x[i_istress], aOff_x[i_istrain] };

  const PylithInt numADev = 3; /* Pass bulk modulus, initial stress, and initial strain. */
  const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_istress], aOff[i_istrain] };
  const PylithInt aOffDev_x[3] = { aOff_x[i_shearModulus], aOff_x[i_istress], aOff_x[i_istrain] };

  assert(_dim == dim);
  assert(2 == numS);
  assert(5 == numA || 6 == numA);
  assert(sOff);
  assert(sOff_x);
  assert(aOff);
  assert(aOff_x);

  pylith_fekernels_IsotropicLinearElasticityPlaneStrain_meanStress_initstate(_dim, _numS, numAMean, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x, aOffMean, aOffMean_x, a, a_t, a_x, t, x, g1);
  pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress_initstate(_dim, _numS, numADev, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x, aOffDev, aOffDev_x, a, a_t, a_x, t, x, g1);
} /* IsotropicLinearElasticityPlaneStrain_g1_initstate */


/* ---------------------------------------------------------------------- */
/** Jf0 function for isotropic linear elasticity plane strain with implicit time stepping.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [2].
 * @param numA Number of registered subfields in auxiliary field [0].
 * @param sOff Offset of registered subfields in solution field [numS].
 * @param sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param a Auxiliary field with all subfields.
 * @param a_t Time derivative of auxiliary field.
 * @param a_x Gradient of auxiliary field.
 * @param t Time for residual evaluation.
 * @param utshift Coefficient for dF/ds_t term in Jacobian.
 * @param x Coordinates of point evaluation.
 * @param f0 Result [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0_uv_implicit(const PylithInt dim,
								      const PylithInt numS,
								      const PylithInt numA,
								      const PylithInt sOff[],
								      const PylithInt sOff_x[],
								      const PylithScalar s[],
								      const PylithScalar s_t[],
								      const PylithScalar s_x[],
								      const PylithInt aOff[],
								      const PylithInt aOff_x[],
								      const PylithScalar a[],
								      const PylithScalar a_t[],
								      const PylithScalar a_x[],
								      const PylithReal t,
								      const PylithReal utshift,
								      const PylithScalar x[],
								      PylithScalar Jf0[])
{ /* IsotropicLinearElasticityPlaneStrain_Jf0_uv_implicit */
  assert(0);
} /* IsotropicLinearElasticityPlaneStrain_Jf0_uv_implicit */


/** Jf0 function for isotropoc linear elasticity plane strain with explicit time stepping.
 *
 * @param dim Spatial dimension.
 * @param numS Number of registered subfields in solution field [2].
 * @param numA Number of registered subfields in auxiliary field [0].
 * @param sOff Offset of registered subfields in solution field [numS].
 * @param sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param s Solution field with all subfields.
 * @param s_t Time derivative of solution field.
 * @param s_x Gradient of solution field.
 * @param aOff Offset of registered subfields in auxiliary field [numA]
 * @param aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param a Auxiliary field with all subfields.
 * @param a_t Time derivative of auxiliary field.
 * @param a_x Gradient of auxiliary field.
 * @param t Time for residual evaluation.
 * @param utshift Coefficient for dF/ds_t term in Jacobian.
 * @param x Coordinates of point evaluation.
 * @param f0 Result [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0_uv_explicit(const PylithInt dim,
								       const PylithInt numS,
								       const PylithInt numA,
								       const PylithInt sOff[],
								       const PylithInt sOff_x[],
								       const PylithScalar s[],
								       const PylithScalar s_t[],
								       const PylithScalar s_x[],
								       const PylithInt aOff[],
								       const PylithInt aOff_x[],
								       const PylithScalar a[],
								       const PylithScalar a_t[],
								       const PylithScalar a_x[],
								       const PylithReal t,
								       const PylithReal utshift,
								       const PylithScalar x[],
								       PylithScalar Jf0[])
{ /* IsotropicLinearElasticityPlaneStrain_Jf0_uv_explicit */
  assert(0);
} /* IsotropicLinearElasticityPlaneStrain_Jf0_uv_explicit*/


/* ---------------------------------------------------------------------- */
/* g3_uu entry function for 2-D plane strain isotropic linear elasticity.
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jg3_uu(const PylithInt dim,
							     const PylithInt numS,
							     const PylithInt numA,
							     const PylithInt sOff[],
							     const PylithInt sOff_x[],
							     const PylithScalar s[],
							     const PylithScalar s_t[],
							     const PylithScalar s_x[],
							     const PylithInt aOff[],
							     const PylithInt aOff_x[],
							     const PylithScalar a[],
							     const PylithScalar a_t[],
							     const PylithScalar a_x[],
							     const PylithReal t,
							     const PylithReal utshift,
							     const PylithScalar x[],
							     PylithScalar Jg3[])
{ /* IsotropicLinearElasticityPlaneStrain_g3_uu */
  const PylithInt _dim = 2;

  const PylithInt _numS = 2;

  const PylithInt _numA = 3;
  const PylithInt i_shearModulus = 1;
  const PylithInt i_bulkModulus = 2;

  const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
  const PylithScalar shearModulus = a[aOff[i_shearModulus]];

  const PylithScalar lambda = bulkModulus - 2.0/3.0*shearModulus;
  const PylithScalar lambda2mu = bulkModulus + 4.0/3.0*shearModulus; 
   
  const PylithReal C1111 = lambda2mu;
  const PylithReal C2222 = lambda2mu;
  const PylithReal C1122 = lambda;
  const PylithReal C1212 = shearModulus;

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA >= numA);
  assert(aOff);

  /* j(f,g,df,dg) = C(f,df,g,dg)

     0: j0000 = C1111
     1: j0001 = C1112 = 0
     4: j0100 = C1121, symmetry C1112 = 0
     5: j0101 = C1122

     2: j0010 = C1211 = 0
     3: j0011 = C1212
     6: j0110 = C1221, symmetry C1212
     7: j0111 = C1222 = 0
  
     8: j1000 = C2111 = 0
     9: j1001 = C2112, symmetry C1212 
    12: j1100 = C2121, symmetry C1212
    13: j1101 = C2122, symmetry C1222 = 0

    10: j1010 = C2211, symmetry C1122
    11: j1011 = C2212, symmetry C1222 = 0
    14: j1110 = C2221, symmetry C1222 = 0
    15: j1111 = C2222
  */

  Jg3[ 0] -= C1111; /* j0000 */
  Jg3[ 3] -= C1212; /* j0011 */
  Jg3[ 5] -= C1122; /* j0101 */
  Jg3[ 6] -= C1212; /* j0110, C1221 */
  Jg3[ 9] -= C1212; /* j1001, C2112 */
  Jg3[10] -= C1122; /* j1010, C2211 */
  Jg3[12] -= C1212; /* j1100, C2121 */
  Jg3[15] -= C2222; /* j1111 */

} /* IsotropicLinearElasticityPlaneStrain_Jg3_uu */


/* ---------------------------------------------------------------------- */
/* Calculate mean stress contribution to f1 for 2-D plane strain
 * isotropic linear elasticity WITHOUT initial stress and strain.
 *
 * meanStress_ij = bulkModulus * strain_kk * delta_ij
 *
 * f1 = -meanStress
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_meanStress(const PylithInt dim,
								 const PylithInt numS,
								 const PylithInt numA,
								 const PylithInt sOff[],
								 const PylithInt sOff_x[],
								 const PylithScalar s[],
								 const PylithScalar s_t[],
								 const PylithScalar s_x[],
								 const PylithInt aOff[],
								 const PylithInt aOff_x[],
								 const PylithScalar a[],
								 const PylithScalar a_t[],
								 const PylithScalar a_x[],
								 const PylithReal t,
								 const PylithScalar x[],
								 PylithScalar f1[])
{ /* IsotropicLinearElasticityPlaneStrain_meanStress */
  const PylithInt _dim = 2;

  const PylithInt i_disp = 0;
  const PylithScalar* disp_x = &s_x[sOff[i_disp]];

  const PylithInt i_bulkModulus = 0;
  const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];

  PylithInt i;
  PylithScalar strainTrace = 0;

  assert(_dim == dim);
  assert(1 == numS);
  assert(1 == numA);
  assert(sOff);
  assert(aOff);
  assert(s_x);
  assert(a);
  assert(f1);

  for (i = 0; i < _dim; ++i) {
    strainTrace += disp_x[i*_dim+i];
  } /* for */
  for (i = 0; i < _dim; ++i) {
    f1[i*_dim+i] -= bulkModulus * strainTrace;
  } /* for */
} /* IsotropicLinearElasticityPlaneStrain_meanStress */


/* ---------------------------------------------------------------------- */
/* Calculate mean stress contribution to f1 for 2-D plane strain
 * isotropic linear elasticity WITH initial stress and initial strain.
 *
 * We compute the stress relative to a reference stress/strain state.
 *
 * meanStress_ij + meanInitialStress_ij = bulkModulus * (strain_kk + initialstrain_kk) * delta_ij
 *
 * f1 = -meanStress_ij
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_meanStress_initstate(const PylithInt dim,
									   const PylithInt numS,
									   const PylithInt numA,
									   const PylithInt sOff[],
									   const PylithInt sOff_x[],
									   const PylithScalar s[],
									   const PylithScalar s_t[],
									   const PylithScalar s_x[],
									   const PylithInt aOff[],
									   const PylithInt aOff_x[],
									   const PylithScalar a[],
									   const PylithScalar a_t[],
									   const PylithScalar a_x[],
									   const PylithReal t,
									   const PylithScalar x[],
									   PylithScalar f1[])
{ /* IsotropicLinearElasticityPlaneStrain_initState_meanStress */
  const PylithInt _dim = 2;

  const PylithInt i_disp = 0;
  const PylithScalar* disp_x = &s_x[sOff[i_disp]];

  const PylithInt i_bulkModulus = 0;
  const PylithInt i_istress = 1;
  const PylithInt i_istrain = 2;
  const PylithScalar bulkModulus = a[aOff[i_bulkModulus]];
  const PylithScalar* initialstress = &a[aOff[i_istress]];
  const PylithScalar* initialstrain = &a[aOff[i_istrain]];

  PylithInt i;
  PylithScalar strainTrace = 0;
  PylithScalar meanistress = 0;

  assert(_dim == dim);
  assert(1 == numS);
  assert(3 == numA);
  assert(sOff);
  assert(aOff);
  assert(s_x);
  assert(a);
  assert(f1);

  for (i=0; i < _dim; ++i) {
    strainTrace += disp_x[i*_dim+i] - initialstrain[i*_dim+i];
    meanistress += initialstress[i*_dim+i];
  } /* for */
  meanistress /= (PylithScalar)_dim;
  for (i = 0; i < _dim; ++i) {
    f1[i*_dim+i] -= bulkModulus * strainTrace - meanistress;
  } /* for */
} /* IsotropicLinearElasticityPlaneStrain_meanStress_initState */


/* ---------------------------------------------------------------------- */
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * elasticity WITHOUT initial stress and strain.
 *
 * devStress_ij = 2.0*shearModulus * strain_ij
 *
 * f1 = -devStress_ij
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress(const PylithInt dim,
								       const PylithInt numS,
								       const PylithInt numA,
								       const PylithInt sOff[],
								       const PylithInt sOff_x[],
								       const PylithScalar s[],
								       const PylithScalar s_t[],
								       const PylithScalar s_x[],
								       const PylithInt aOff[],
								       const PylithInt aOff_x[],
								       const PylithScalar a[],
								       const PylithScalar a_t[],
								       const PylithScalar a_x[],
								       const PylithReal t,
								       const PylithScalar x[],
								       PylithScalar f1[])
{ /* IsotropicLinearElasticityPlaneStrain_deviatoricStress */
  const PylithInt _dim = 2;

  const PylithInt i_disp = 0;
  const PylithScalar* disp_x = &s_x[sOff[i_disp]];

  const PylithInt i_shearModulus = 0;
  const PylithScalar shearModulus = a[aOff[i_shearModulus]];

  PylithInt i, j;

  assert(_dim == dim);
  assert(1 == numS);
  assert(1 == numA);
  assert(sOff);
  assert(aOff);
  assert(s_x);
  assert(a);
  assert(f1);

  for (i=0; i < _dim; ++i) {
    for (j=0; j < _dim; ++j) {
      f1[i*_dim+j] -= shearModulus * (disp_x[i*_dim+j] + disp_x[j*_dim+i]);
    } /* for */
  } /* for */
} /* IsotropicLinearElasticityPlaneStrain_deviatoricStress */


/* ---------------------------------------------------------------------- */
/* Calculate deviatoric stress for 2-D plane strain isotropic linear
 * elasticity WITH initial stress and initial strain.
 *
 * devStress_ij + (initialStress_ij - meanInitialStress_ij * delta_ij) = 2.0*shearModulus * (strain_ij + initialstrain_ij)
 *
 * f1 = -devStress_ij
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress_initstate(const PylithInt dim,
										 const PylithInt numS,
										 const PylithInt numA,
										 const PylithInt sOff[],
										 const PylithInt sOff_x[],
										 const PylithScalar s[],
										 const PylithScalar s_t[],
										 const PylithScalar s_x[],
										 const PylithInt aOff[],
										 const PylithInt aOff_x[],
										 const PylithScalar a[],
										 const PylithScalar a_t[],
										 const PylithScalar a_x[],
										 const PylithReal t,
										 const PylithScalar x[],
										 PylithScalar f1[])
{ /* deviatoricStress_IsotropicLinearElasticityPlaneStrain_initstate */
  const PylithInt _dim = 2;

  const PylithInt _numS = 1;
  const PylithInt i_disp = 0;
  const PylithScalar* disp_x = &s_x[sOff[i_disp]];

  const PylithInt _numA = 3;
  const PylithInt i_shearModulus = 0;
  const PylithInt i_istress = 1;
  const PylithInt i_istrain = 2;
  const PylithScalar shearModulus = a[aOff[i_shearModulus]];
  const PylithScalar* initialstress = &a[aOff[i_istress]];
  const PylithScalar* initialstrain = &a[aOff[i_istrain]];

  PylithInt i, j;
  PylithScalar meanistress = 0;

  assert(_dim == dim);
  assert(_numS == numS);
  assert(_numA == numA);
  assert(sOff);
  assert(aOff);
  assert(s_x);
  assert(a);
  assert(f1);

  for (i=0; i < _dim; ++i) {
    meanistress += initialstress[i*_dim+i];
  } /* for */
  meanistress /= (PylithScalar)_dim;
  for (i=0; i < _dim; ++i) {
    for (j=0; j < _dim; ++j) {
      f1[i*_dim+j] -= shearModulus * (disp_x[i*_dim+j] + disp_x[j*_dim+i] + initialstrain[i*_dim+j]) - initialstress[i*_dim+j];
    } /* for */
    f1[i*_dim+i] -= meanistress;
  } /* for */
} /* deviatoricStress_IsotropicLinearElasticityPlaneStrain_initstate */


/* ====================================================================== 
 * Kernels for incompressible elasticity.
 *
 * Solution fields:  [disp(dim), vel(dim), pres(1)]
 * Auxiliary fields: [density(1), shearModulus(1), bulkModulus(1), bodyforce(dim), initialstrain(dim*dim), initialstress(dim*dim)]
 *
 * ====================================================================== 
 */

/* ---------------------------------------------------------------------- */
/** f0 function for isotropic linear incompressible elasticity plane strain.
 * Solution fields:  [disp(dim), vel(dim)]
 * Auxiliary fields: [density(1), ...]
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_f0(
						const PylithInt dim,
						const PylithInt numS,
						const PylithInt numA,
						const PylithInt sOff[],
						const PylithInt sOff_x[],
						const PylithScalar s[],
						const PylithScalar s_t[],
						const PylithScalar s_x[],
						const PylithInt aOff[],
						const PylithInt aOff_x[],
						const PylithScalar a[],
						const PylithScalar a_t[],
						const PylithScalar a_x[],
						const PylithReal t,
						const PylithScalar x[],
						PylithScalar f0[])
{ /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_f0 */
  const PylithInt _dim = 2;

  const PylithInt _numS = 2; /* Passing on displacement and velocity. */
  const PylithInt i_disp = 0;
  const PylithInt i_vel = 1;

  const PylithInt _numA = 1; /* Passing on density. */
  const PylithInt i_density = 0;


  assert(_dim == dim);
  assert(3 == numS);
  assert(1 <= numA);
  assert(aOff);
  assert(aOff_x);

  const PylithInt _sOff[2] = {sOff[i_disp], sOff[i_vel]};
  const PylithInt _sOff_x[2] = {sOff_x[i_disp], sOff_x[i_vel]};

  pylith_fekernels_Elasticity_f0_inertia(_dim, _numS, _numA, _sOff, _sOff_x, s, s_t, s_x, &aOff[i_density], &aOff_x[i_density], a, a_t, a_x, t, x, f0);
} /* pylith_fekernels_IsotropicLinearElasticityIncompPlaneStrain_f0 */


/* ---------------------------------------------------------------------- */
/** g0 function for isotropic linear incompressible elasticity plane strain.
 *
 * Solution fields:  [disp(dim), vel(dim)]
 * Auxiliary fields: [density(1), shearModulus(1), bulkModulus(1), bodyforce(dim), ...]
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g0(
						const PylithInt dim,
						const PylithInt numS,
						const PylithInt numA,
						const PylithInt sOff[],
						const PylithInt sOff_x[],
						const PylithScalar s[],
						const PylithScalar s_t[],
						const PylithScalar s_x[],
						const PylithInt aOff[],
						const PylithInt aOff_x[],
						const PylithScalar a[],
						const PylithScalar a_t[],
						const PylithScalar a_x[],
						const PylithReal t,
						const PylithScalar x[],
						PylithScalar g0[])
{ /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g0 */
  const PylithInt _dim = 2;

  const PylithInt _numS = 2; /* Passing on displacement and velocity. */
  const PylithInt i_disp = 0;
  const PylithInt i_vel = 1;

  const PylithInt _numA = 1; /* Passing on body force. */
  const PylithInt i_bodyforce = 3;

  assert(_dim == dim);
  assert(3 == numS);
  assert(4 == numA || 6 == numA);
  assert(aOff);
  assert(aOff_x);

  const PylithInt _sOff[2] = {sOff[i_disp], sOff[i_vel]};
  const PylithInt _sOff_x[2] = {sOff_x[i_disp], sOff_x[i_vel]};

  pylith_fekernels_Elasticity_g0_bodyforce(_dim, _numS, _numA, _sOff, _sOff_x, s, s_t, s_x, &aOff[i_bodyforce], &aOff_x[i_bodyforce], a, a_t, a_x, t, x, g0);
} /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g0 */


/* ---------------------------------------------------------------------- */
/** g1 function for isotropic linear incompressible elasticity plane
 * strain WITHOUT initial stress and strain.
 *
 * Solution fields: [disp(dim), vel(dim), pres]
 * Auxiliary fields: [density(1), shearModulus(1), bulkModulus(1), ...]
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g1(
						const PylithInt dim,
						const PylithInt numS,
						const PylithInt numA,
						const PylithInt sOff[],
						const PylithInt sOff_x[],
						const PylithScalar s[],
						const PylithScalar s_t[],
						const PylithScalar s_x[],
						const PylithInt aOff[],
						const PylithInt aOff_x[],
						const PylithScalar a[],
						const PylithScalar a_t[],
						const PylithScalar a_x[],
						const PylithReal t,
						const PylithScalar x[],
						PylithScalar g1[])
{ /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g1 */
  const PylithInt _dim = 2;

  const PylithInt _numS = 1; /* Passing displacement field. */
  const PylithInt i_disp = 0;
  const PylithInt i_pres = 2;
  const PylithScalar pres = s[sOff[i_pres]];

  const PylithInt i_shearModulus = 1;

  const PylithInt numADev = 1; /* Passing shear modulus */
  const PylithInt aOffDev[1] = { aOff[i_shearModulus] };
  const PylithInt aOffDev_x[1] = { aOff_x[i_shearModulus] };

  PylithInt i;

  assert(_dim == dim);
  assert(3 == numS);
  assert(3 == numA || 4 == numA);
  assert(sOff);
  assert(sOff_x);
  assert(aOff);
  assert(aOff_x);

  pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress(_dim, _numS, numADev, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x, aOffDev, aOffDev_x, a, a_t, a_x, t, x, g1);

  for (i=0; i < _dim; ++i) {
    g1[i*_dim+i] += pres;
  } /* for */
} /* pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g1 */


/* ---------------------------------------------------------------------- */
/** g1 function for isotropic linear incompressible elasticity plane strain
 * with initial stress and strain.
 *
 * Solution fields:  [disp(dim), vel(dim), pres]
 * Auxiliary fields: [density(1), shearModulus(1), bulkModulus(1), ..., initialstress(dim*dim), initialstrain(dim*dim)]
 */
void
pylith_fekernels_IsotropicLinearIncompElasticityPlaneStrain_g1_initstate(
						const PylithInt dim,
						const PylithInt numS,
						const PylithInt numA,
						const PylithInt sOff[],
						const PylithInt sOff_x[],
						const PylithScalar s[],
						const PylithScalar s_t[],
						const PylithScalar s_x[],
						const PylithInt aOff[],
						const PylithInt aOff_x[],
						const PylithScalar a[],
						const PylithScalar a_t[],
						const PylithScalar a_x[],
						const PylithReal t,
						const PylithScalar x[],
						PylithScalar g1[])
{ /* IsotropicLinearIncompElasticityPlaneStrain_g1_initstate */
  const PylithInt _dim = 2;

  const PylithInt _numS = 1; /* Passing displacement field. */
  const PylithInt i_disp = 0;
  const PylithInt i_pres = 2;

  const PylithInt numADev = 3;
  const PylithInt i_shearModulus = 1;
  const PylithInt i_istress = numA-2;
  const PylithInt i_istrain = numA-1;

  const PylithScalar pres = s[sOff[i_pres]];

  const PylithInt aOffDev[3] = { aOff[i_shearModulus], aOff[i_istress], aOff[i_istrain] };
  const PylithInt aOffDev_x[3] = { aOff_x[i_shearModulus], aOff_x[i_istress], aOff_x[i_istrain] };

  PylithInt i;

  assert(_dim == dim);
  assert(3 == numS);
  assert(5 == numA || 6 == numA);
  assert(sOff);
  assert(sOff_x);
  assert(aOff);
  assert(aOff_x);

  pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress_initstate(_dim, _numS, numADev, &sOff[i_disp], &sOff_x[i_disp], s, s_t, s_x, aOffDev, aOffDev_x, a, a_t, a_x, t, x, g1);

  for (i=0; i < _dim; ++i) {
    g1[i*_dim+i] += pres;
  } /* for */
} /* IsotropicLinearIncompElasticityPlaneStrain_g1_initstate */


/* ---------------------------------------------------------------------- */
/* g3_uu entry function for 2-D plane strain incompressible isotropic linear
 * elasticity.
 *
 * Solution fields:  [disp(dim), vel(dim), pres]
 * Auxiliary fields: [density(1), shearModulus(1), bulkModulus[1], ...]
 */
void
pylith_fekernels_IsotropicIncompLinearElasticityPlaneStrain_Jg3_uu(
					const PylithInt dim,
					const PylithInt numS,
					const PylithInt numA,
					const PylithInt sOff[],
					const PylithInt sOff_x[],
					const PylithScalar s[],
					const PylithScalar s_t[],
					const PylithScalar s_x[],
					const PylithInt aOff[],
					const PylithInt aOff_x[],
					const PylithScalar a[],
					const PylithScalar a_t[],
					const PylithScalar a_x[],
					const PylithReal t,
					const PylithReal utshift,
					const PylithScalar x[],
					PylithScalar Jg3[])
{ /* pylith_fekernels_IsotropicIncompLinearElasticityPlaneStrain_Jg3_uu */
  const PylithInt _dim = 2;

  const PylithInt i_shearModulus = 1;
  const PylithScalar shearModulus = a[aOff[i_shearModulus]];

  const PylithReal C1111 = 5.0 * shearModulus/3.0;
  const PylithReal C2222 = C1111;
  const PylithReal C1122 = -shearModulus/3.0;
  const PylithReal C1212 = shearModulus;

  assert(_dim == dim);
  assert(3 == numS);
  assert(3 <= numA);
  assert(aOff);

  /* g(f,g,df,dg) = C^{\prime}(f,df,g,dg) - \frac{1}{6}C^{\prime}(f,df,h,h) \delta(g,dg)

     0: g0000 = C1111
     1: g0001 = C1112
     4: g0100 = C1121, symmetry C1112
     5: g0101 = C1122

     2: g0010 = C1211
     3: g0011 = C1212
     6: g0110 = C1221, symmetry C1212
     7: g0111 = C1222
  
     8: g1000 = C2111
     9: g1001 = C2112, symmetry C1212
    12: g1100 = C2121, symmetry C1212
    13: g1101 = C2122, symmetry C1222

    10: g1010 = C2211, symmetry C1122
    11: g1011 = C2212, symmetry C1222
    14: g1110 = C2221, symmetry C1222
    15: g1111 = C2222
  */

  Jg3[ 0] -= C1111; /* g0000 */
  Jg3[ 3] -= C1212; /* g0011 */
  Jg3[ 5] -= C1122; /* g0101 */
  Jg3[ 6] -= C1212; /* g0110, C1221 */
  Jg3[ 9] -= C1212; /* g1001, C2112 */
  Jg3[10] -= C1122; /* g1010, C2211 */
  Jg3[12] -= C1212; /* g1100, C2121 */
  Jg3[15] -= C2222; /* g1111 */
} /* pylith_fekernels_IsotropicIncompLinearElasticityPlaneStrain_Jg3_uu */
					      

/* ---------------------------------------------------------------------- */
/* g2_uv entry function for 2-D plane strain incompressible isotropic linear
 * elasticity.
 *
 * Solution fields = [disp(dim), vel(dim), pres]
 * Auxiliary fields = [density, shearModulus, bulkModulus, ...]
 */
void
pylith_fekernels_IsotropicIncompLinearElasticityPlaneStrain_Jg2_up(
					const PylithInt dim,
					const PylithInt numS,
					const PylithInt numA,
					const PylithInt sOff[],
					const PylithInt sOff_x[],
					const PylithScalar s[],
					const PylithScalar s_t[],
					const PylithScalar s_x[],
					const PylithInt aOff[],
					const PylithInt aOff_x[],
					const PylithScalar a[],
					const PylithScalar a_t[],
					const PylithScalar a_x[],
					const PylithReal t,
					const PylithReal utshift,
					const PylithScalar x[],
					PylithScalar Jg2[])
{ /* pylith_fekernels_IsotropicIncompLinearElasticityPlaneStrain_Jg2_up */
  const PylithInt _dim = 2;

  PylithInt i;

  assert(3 == numS);
  assert(3 <= numA);

  for (i=0; i < _dim; ++i) {
    Jg2[i*_dim+i] += 1.0;
  } /* for */
} /* pylith_fekernels_IsotropicIncompLinearElasticityPlaneStrain_Jg2_up */
