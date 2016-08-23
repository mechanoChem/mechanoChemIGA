#ifndef physicsHeaders_
#define physicsHeaders_

//extern "C" {
#include "petiga.h"
//}

#include "appCtx.h"

template <class T,unsigned int DIM, unsigned int DOF>
PetscErrorCode quadPtResidual(IGAPoint p,PetscReal dt2,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const T * U,
			PetscReal t0,const PetscScalar * U0,
			T *R,void *ctx);

template <int dim>
PetscErrorCode ProjectionResidual(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *ctx);

PetscErrorCode ProjectionJacobian(IGAPoint p, const PetscScalar *U, PetscScalar *K, void *ctx);

template <int dim>
PetscErrorCode ProjectSolution(IGA iga, PetscInt step, Vec U, AppCtx *user);

template <int dim>
PetscErrorCode OutputMonitor(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);

#endif
