#ifndef applicationHeaders_
#define applicationHeaders_

//extern "C" {
#include "petiga.h"
//}

#include "appCtx.h"

template<unsigned int dim>
struct Field;

template<>
struct Field<2>{
  PetscReal Ux, Uy;
  PetscReal ux, uy;
};

template<>
struct Field<3>{
  PetscReal Ux, Uy, Uz;
  PetscReal ux, uy, uz;
};

PetscErrorCode FormInitialCondition(IGA iga, Vec U, AppCtx *user);

template<unsigned int DIM>
int setup(AppCtx& user,Vec *U,Vec *U0,TS &ts);

template<unsigned int dim>
int boundaryConditions(AppCtx& user, double scale);

template<int dim>
int timeStepSetup(AppCtx& user, TS& ts);

template <int dim>
PetscErrorCode loadStep(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);

template <int dim>
PetscErrorCode adaptiveTimeStep(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);

PetscErrorCode SNESConvergedTest(SNES snes, PetscInt it,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx);

PetscErrorCode SNESConvergedTest_Interactive(SNES snes, PetscInt it,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx);

int setConvergenceTest(AppCtx& user, TS& ts);

int defineParameters(AppCtx& user);

#endif
