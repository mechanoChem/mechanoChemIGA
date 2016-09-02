#ifndef physicsHeaders_
#define physicsHeaders_

//extern "C" {
#include "petiga.h"
//}

#include "appCtx.h"

/**
 * The heart of the code. This computes the residual at the quadrature point level, including all of the associated stress, strain, etc.
*/

template <class T,unsigned int DIM, unsigned int DOF>
PetscErrorCode quadPtResidual(IGAPoint p,PetscReal dt2,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const T * U,
			PetscReal t0,const PetscScalar * U0,
			T *R,void *ctx);

/**
 * Compute the residual for the projection fields.
*/

template <int dim>
PetscErrorCode ProjectionResidual(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *ctx);

/**
 * Compute the jacobiam for the projection fields using Sacado automatic differentiation.
*/

PetscErrorCode ProjectionJacobian(IGAPoint p, const PetscScalar *U, PetscScalar *K, void *ctx);

/**
 * Function to project and output desired fields.
*/

template <int dim>
PetscErrorCode ProjectSolution(IGA iga, PetscInt step, Vec U, AppCtx *user);

/**
 * Function to define output.
*/

template <int dim>
PetscErrorCode OutputMonitor(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);

#endif
