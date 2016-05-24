#ifndef genericHeaders_
#define genericHeaders_

//extern "C" {
#include "petiga.h"
//}
#include "../applications/configurationalForces/3D/applicationHeaders.h"

PetscErrorCode Residual(IGAPoint p,PetscReal dt,
                        PetscReal shift,const PetscScalar *V,
                        PetscReal t,const PetscScalar *U,
                        PetscReal t0,const PetscScalar *U0, 
			PetscScalar *R,void *ctx);

template<unsigned int DOF>
PetscErrorCode Jacobian(IGAPoint p,PetscReal dt,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const PetscScalar *U,
			PetscReal t0,const PetscScalar *U0,
			PetscScalar *K,void *ctx);

enum fieldType{SCALAR,VECTOR,TENSOR};

template <class T, unsigned int dim, unsigned int dof>
void computeField(fieldType type, unsigned int index, IGAPoint p, const T* U, T* _value=0, T* _grad=0, T* _hess=0);

template<unsigned int DOF>
int init(AppCtx& user, PetscInt N, PetscInt p);

#endif
