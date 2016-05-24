#ifndef physicsHeaders_
#define physicsHeaders_

//extern "C" {
#include "petiga.h"
//}
#include "../../applications/configurationalForces/3D/applicationHeaders.h"

typedef struct {
  PetscReal Ux, Uy;
#if DIM==3
  PetscReal Uz;
#endif
  PetscReal ux, uy;
#if DIM==3
  PetscReal uz;
#endif
} Field;

PetscErrorCode FormInitialCondition(IGA iga, Vec U, AppCtx *user);

template <class T>
PetscErrorCode Function(IGAPoint p,PetscReal dt2,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const T * U,
			PetscReal t0,const PetscScalar * U0,
			T *R,void *ctx);

PetscErrorCode E22Function(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *ctx);

PetscErrorCode E22Jacobian(IGAPoint p, const PetscScalar *U, PetscScalar *K, void *ctx);

PetscErrorCode ProjectSolution(IGA iga, PetscInt step, Vec U, AppCtx *user);

template <int dim>
PetscErrorCode OutputMonitor(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);

#endif
