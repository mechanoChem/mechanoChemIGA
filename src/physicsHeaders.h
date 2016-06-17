#ifndef physicsHeaders_
#define physicsHeaders_

//extern "C" {
#include "petiga.h"
//}
#include "../applications/configurationalForces/applicationHeaders.h"
//#include "../../applications/configurationalForces/3D/parameters.h"

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

template <class T,unsigned int DIM, unsigned int DOF>
PetscErrorCode Function(IGAPoint p,PetscReal dt2,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const T * U,
			PetscReal t0,const PetscScalar * U0,
			T *R,void *ctx);

template<unsigned int DIM, unsigned int DOF>
PetscErrorCode Residual(IGAPoint p,PetscReal dt,
                        PetscReal shift,const PetscScalar *V,
                        PetscReal t,const PetscScalar *U,
                        PetscReal t0,const PetscScalar *U0, 
			PetscScalar *R,void *ctx);

template<unsigned int DIM, unsigned int DOF>
PetscErrorCode Jacobian(IGAPoint p,PetscReal dt,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const PetscScalar *U,
			PetscReal t0,const PetscScalar *U0,
			PetscScalar *K,void *ctx);

template <int dim>
PetscErrorCode E22Function(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *ctx);

PetscErrorCode E22Jacobian(IGAPoint p, const PetscScalar *U, PetscScalar *K, void *ctx);

template <int dim>
PetscErrorCode ProjectSolution(IGA iga, PetscInt step, Vec U, AppCtx *user);

template <int dim>
PetscErrorCode OutputMonitor(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);

#endif
