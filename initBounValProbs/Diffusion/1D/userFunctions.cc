#include "userFunctions.h"

template<unsigned int dim>
double userScalarInitialConditions(const Tensor<1,dim,double> &x, unsigned int scalar_i, AppCtx<dim> &user)
{

  return x[0];

} //end scalarInitialConditions

template<unsigned int dim>
void defineParameters(AppCtx<dim>& user){
 
  user.N[0] = 100;
  user.L[0] = 1.;

  //Set the domain to be periodic in the x directions
  //user.periodic[0] = PETSC_TRUE;

  user.dtVal = .1;
  user.totalTime = 20;
  user.RESTART_IT = 0;
  user.RESTART_TIME = 0.;
  user.skipOutput = 1;

  user.scalarSolnFields.push_back("c");
  user.polyOrder = 2;
  user.globalContinuity = 1;

  user.scalarInitialConditions = userScalarInitialConditions;

} //end defineParameters

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
	      Sacado::Fad::SimpleFad<T> &r){

  //Chemistry
  double dt = user.dt;
  double jn = 0.;
  double D = 1.; //Mobility
  
  //Fickian diffusion - c
  r = ( w1.val(0)*(c.val(0) - c.valP(0))/dt )*dV;  
  r += D*w1.grad(0)*c.grad(0)*dV;
  r += -w1.val(0)*jn*dS; //Boundary flux
  
} //end residual

#include "userFunctionsInstantiation.h"
