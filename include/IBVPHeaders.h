#ifndef applicationHeaders_
#define applicationHeaders_

//extern "C" {
#include "petiga.h"
//}

#include "appCtx.h"

PetscErrorCode FormInitialCondition(IGA iga, Vec U, AppCtx *user);

template<unsigned int DIM,unsigned int DOF>
int setup(AppCtx& user,Vec *U,Vec *U0,TS &ts);

int boundaryConditions(AppCtx& user, double scale);

PetscErrorCode loadStep(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);

PetscErrorCode adaptiveTimeStep(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);

PetscErrorCode SNESConvergedTest(SNES snes, PetscInt it,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx);

int defineParameters(AppCtx& user);

#endif
