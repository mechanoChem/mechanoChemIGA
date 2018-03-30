#ifndef coreFunctions_h
#define coreFunctions_h

//extern "C" {
#include "petiga.h"
//}
#include "appCtx.h"

enum fieldType{SCALAR,VECTOR,TENSOR};

/**
 * @defgroup coreFunctions Core functions
 * This is a collection of the core functions of the code. These functions do not require modification
 * by the user.
 */


/**
 * @ingroup coreFunctions
 * Calls quadPtResidual to compute the residual.
 */

template<unsigned int DIM>
PetscErrorCode Residual(IGAPoint p,
			const PetscScalar *U,
			const PetscScalar *Up, 
			const PetscScalar *Upp, 
			PetscScalar *R,
			void *ctx);

/**
 * @ingroup coreFunctions
 * Computes the Jacobian using the residual and Sacado automatic differentiation.
 */

template<unsigned int DIM>
PetscErrorCode Jacobian(IGAPoint p,
			const PetscScalar *U,
			const PetscScalar *Up,
			const PetscScalar *Upp, 
			PetscScalar *K,
			void *ctx);

/**
 * @ingroup coreFunctions
 * Function using PetIGA to retrieve a the values for a field, the gradients, and the hessians.
 */

template <class T, unsigned int dim>
  void ComputeField(fieldType type, unsigned int index, IGAPoint p, const T* U, unsigned int Udof, T* _value=0, T* _grad=0, T* _hess=0);

/**
 * @ingroup coreFunctions
 * Function to setup data structures, etc
 */

template<unsigned int DIM>
int Setup(AppCtx<DIM>& user,Vec *U,Vec *Up,Vec *Upp,SNES &snes);

/**
 * @ingroup coreFunctions
 * Function to initilize PetIGA structures. Called in the "setup" function.
 */

template<unsigned int DIM>
int InitIGA(AppCtx<DIM>& user, IGA &iga, PetscInt p, PetscInt C, PetscInt DOF);

/**
 * @ingroup coreFunctions
 * Form the initial conditions of the IBVP based on user defined functions.
 */

template<unsigned int DIM>
PetscErrorCode FormInitialCondition(IGA iga, Vec U, AppCtx<DIM> *user);

/**
 * @ingroup coreFunctions
 * This computes the residual at the quadrature point level, calling the userFunction residual.
 */

template <class T,unsigned int DIM>
  PetscErrorCode QuadPtResidual(IGAPoint p,
				const T * U,
				const PetscScalar * Up,
				const PetscScalar * Upp,
				T *R,
				void *ctx);

/**
 * @ingroup coreFunctions
 * Compute the residual for the projection fields.
 */

template <int dim>
PetscErrorCode ProjectionResidual(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *ctx);

/**
 * @ingroup coreFunctions
 * Compute the jacobiam for the projection fields using Sacado automatic differentiation.
 */

PetscErrorCode ProjectionJacobian(IGAPoint p, PetscScalar *K, void *ctx);

/**
 * @ingroup coreFunctions
 * Function to project and output desired fields.
 */

template <int dim>
PetscErrorCode ProjectSolution(IGA iga, PetscInt step, Vec U, AppCtx<dim> *user);

/**
 * @ingroup coreFunctions
 * Function to define output.
 */

template <int dim>
PetscErrorCode StepUpdate(PetscInt it_number,PetscReal c_time,Vec U,AppCtx<dim> &user);

/**
 * @ingroup coreFunctions
 * Function defining the workflow.
 */

template <int dim>
PetscErrorCode Run();

/**
 * @ingroup coreFunctions
 * Function that reads the parameter file, if it exists.
 */

template <int dim>
int ReadParameters(AppCtx<dim> &user);

/**
 * @ingroup coreFunctions
 * Function reads the dimension from the parameters file, if it exists.
 */

int ReadDim();

#endif //coreFunctions_h
