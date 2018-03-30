#ifndef userFunctions_h
#define userFunctions_h

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
 * @defgroup userFunctions User defined functions
 * This is a collection of fuctions that require (viz. defineParameters, residual) or allow user definitions.
 * The optional user defined functions involve defining a function and assigning it to one of the function pointers
 * included in the AppCtx structure.
 */

/**
 * @ingroup userFunctions
 * Function to define any necessary variables in the AppCtx struct. For example, this function is used to set material parameters,
 * define the problem domain and the mesh, set the basis function order, etc. It is also used to set available function pointers to
 * user defined functions. Called once at the beginning of the code.
 *
 * This function must be defined by the user.
 */

template<unsigned int dim>
void defineParameters(AppCtx<dim>& user);

/**
 * @ingroup userFunctions
 * Function defining the residual of the problem, in a format closely following the analytical weak form. The value of
 * the residual is stored in the variable r
 *
 * @param dV is a boolean, "true" if \c residual is being called for the volume integral and 
 * "false" if \c residual is being called for the surface integral.
 * @param dS is a boolean, "false" if \c residual is being called for the volume integral and 
 * "true" if \c residual is being called for the surface integral.
 * @param x gives the coordinates of the quadrature point.
 * @param normal gives the unit normal for a surface quadrature point.
 * @param c gives the information (values, gradients, etc.) for the scalar solution fields at the current quadrature point
 * (see documentation for solutionScalars class).
 * @param u gives the information (values, gradients, etc.) for the vector solution fields at the current quadrature point
 * (see documentation for solutionVectors class).
 * @param w1 gives the information for the scalar test functions.
 * @param w2 gives the information for the vector test functions.
 * @param user is a structure available for parameters related to the initial boundary value problem (e.g. elasticity tensor).
 * @param r stores the scalar value of the residual for the weak form of the PDE which is then used by the core assembly functions.
 *
 * See also the documentation for the Tensor class.
 *
 * As an example, the weak form for finite strain elasticity:
 * \f$\int_\Omega (\nabla{\boldsymbol{w}}:\boldsymbol{P}) d\mathrm{V} - \int_{\partial\Omega} (\boldsymbol{w}\cdot\boldsymbol{h}) d\mathrm{S}\f$
 * would be represented by the following line of code:
 *
 * \code{.cpp} r = double_contract(w2.grad(0),P)*dV - ((w2.val(0)*h)*dS; \endcode
 *
 * where P is a user-defined 2nd order tensor and h is a user-defined 1st order tensor.
 *
 * This function must be defined by the user.
 */

template<unsigned int dim, typename T>
void residual(bool dV,
	      bool dS,
	      const Tensor<1,dim,double> &x,
	      const Tensor<1,dim,double> &normal,
	      const solutionScalars<dim,T> &c,
	      const solutionVectors<dim,T> &u,
	      const testFunctionScalars<dim,T> &w1,
	      const testFunctionVectors<dim,T> &w2,
	      AppCtx<dim> &user,
	      Sacado::Fad::SimpleFad<T> &r);

#endif //userFunctions_h
