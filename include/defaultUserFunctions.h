#ifndef defaultUserFunctions_h
#define defaultUserFunctions_h

//extern "C" {
#include "petiga.h"
//}

#include "tensor.h"
#include "solutionClass.h"
#include "testFunctionClass.h"
#include "appCtx.h"

//include automatic differentiation library
#include <Sacado.hpp>
//#include <Sacado_Fad_BLAS.hpp>
#include <Sacado_Fad_SimpleFad.hpp>


/**
 * @ingroup userFunctions
 * Define the initial conditions of the IBVP for scalar fields. 
 * Return the initial value for scalar solution field with label scalar_i for the node located at the point x.
 *
 * Defaults to a zero initial condition.
 */

template<unsigned int dim>
double scalarInitialConditions(const Tensor<1,dim,double> &x, unsigned int scalar_i, AppCtx<dim> &user){
  return 0;
}

/**
 * @ingroup userFunctions
 * Define the initial conditions of the IBVP for vector fields.
 * Return the initial value (1st order tensor) for vector solution field with label vector_i for the node located at the point x.
 *
 * Defaults to a zero initial condition.
 */

template<unsigned int dim>
Tensor<1,dim,double> vectorInitialConditions(const Tensor<1,dim,double> &x, unsigned int vector_i, AppCtx<dim> &user){
  Tensor<1,dim,double> tmp;
  return tmp;
};

/**
 * @ingroup userFunctions
 * Define the Dirichlet boundary conditions using PetIGA's function. This function defines a uniform value for a given
 * solution field on a face, using the following notation:
 *
 * IGASetBoundaryValue(user.iga,"axis (0,1,2)","boundary (0,1)","dof_index","value")
 *
 * Potentially called multiple times. The input scale can be used to scale the value of the boundary condition
 * if load stepping is being used.
 *
 * Defaults to no Dirichlet boundary conditions.
 */

template<unsigned int dim>
void boundaryConditions(AppCtx<dim>& user, double scale){};

/**
 * @ingroup userFunctions
 * Define any load stepping (i.e. change in Dirichlet b.c.s). Called at each time/load step. This would usually
 * involve calling the boundaryConditions function with an adjusted value for the input scale.
 *
 * Defaults to no load stepping.
 */

template<unsigned int dim>
void loadStep(unsigned int step,AppCtx<dim> &user){};

/**
 * @ingroup userFunctions
 * Specify any adaptive time step schemes (i.e. change the value of dt based on current conditions). Called at each time/load step.
 *
 * Defaults to no adaptive time stepping.
 */

template<unsigned int dim>
void adaptiveTimeStep(unsigned int step,AppCtx<dim> &user){};

/**
 * @ingroup userFunctions
 * Function defining the scalar and vector fields (defined at the quadrature point) to be projected to the nodes 
 * and included in the output file.
 *
 * @param x gives the coordinates of the quadrature point.
 * @param normal gives the unit normal for a surface quadrature point.
 * @param c gives the information (values, gradients, etc.) for the scalar solution fields at the current quadrature point
 * (see documentation for solutionScalars class).
 * @param u gives the information (values, gradients, etc.) for the vector solution fields at the current quadrature point
 * (see documentation for solutionVectors class).
 * @param user is a structure available for parameters related to the initial boundary value problem (e.g. elasticity tensor).
 * @param scalarProjections stores the scalar values calculated at the quadrature points that are to be projected to the nodes.
 * @param vectorProjectionsvector stores the vector values calculated at the quadrature points that are to be projected to the nodes.
 *
 * The stored values and location in the \c scalarProjections and \c vectorProjectionsvector vectors should correspond
 * the name and location defined in the defineParameters function.
 *
 * Defaults to zero valued fields.
 */

template<unsigned int dim>
void projectFields(const Tensor<1,dim,double> &x,
			  const Tensor<1,dim,double> &normal,
			  const solutionScalars<dim,double> &c,
			  const solutionVectors<dim,double> &u,
			  AppCtx<dim> &user,
			  std::vector<double> &scalarProjections,
			  std::vector<Tensor<1,dim,double> > &vectorProjections){};


#endif //defaultUserFunctions_h
