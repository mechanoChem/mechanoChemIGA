#ifndef genericHeaders_
#define genericHeaders_

//extern "C" {
#include "petiga.h"
//}
#include "appCtx.h"
#include "physicsHeaders.h"

//For a compile time "pow" function:
template <int base,int exp>
struct power{
	static const int value = base*power<base,exp-1>::value;
};
template <int base>
struct power<base,0>{
	static const int value = 1;
};

enum fieldType{SCALAR,VECTOR,TENSOR};

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

template <class T, unsigned int dim, unsigned int dof>
void computeField(fieldType type, unsigned int index, IGAPoint p, const T* U, T* _value=0, T* _grad=0, T* _hess=0);

template<unsigned int DIM, unsigned int DOF>
int initIGA(AppCtx& user, PetscInt p);

#endif
