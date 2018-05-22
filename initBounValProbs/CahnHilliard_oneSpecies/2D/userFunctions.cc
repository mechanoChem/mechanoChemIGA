#include "userFunctions.h"

template<unsigned int dim>
double userScalarInitialConditions(const Tensor<1,dim,double> &x, unsigned int scalar_i, const AppCtx<dim> &user)
{

  return 0.5 + 0.3*(0.5 - (double)(rand() % 100 )/100.0); //Random about 0.5

} //end scalarInitialConditions

template<typename T>
T F_cc(T c, double alpha, double ca, double cb){
  
  //Second derivative of the free energy density: f(c) = alpha*(c-ca)^2*(c-cb)^2
  return 2*alpha*(std::pow(2.*c-ca-cb,2) + 2.*(c-ca)*(c-cb));
  
} //end F_cc

template<unsigned int dim>
void defineParameters(AppCtx<dim>& user){
 
  user.N[0] = 100;
  user.N[1] = 100;
  user.L[0] = 1.;
  user.L[1] = 1.;

  //Set the domain to be periodic in the x directions
  user.periodic[0] = PETSC_TRUE;

  //Define some material parameters (can be overwritten by parameters file)
  user.matParam["inFlux"] = 0; //Infux through top
  user.matParam["mobility"] = .1; //Mobility
  user.matParam["kappa"] = .0001; //Gradient energy parameter

  //Define some free energy parameters
  user.matParam["alpha"] = 0.25; //Free energy coefficient
  user.matParam["c_a"] = 0.2; //Composition of phase a
  user.matParam["c_b"] = 0.9; //Composition of phase b

  user.dtVal = .1;
  user.totalTime = 20;
  user.RESTART_IT = 0;
  user.RESTART_TIME = 0.;
  user.skipOutput = 2;

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
  double jn = user.matParam["inFlux"]*(x[1]==user.L[1]); //Influx through the top
  double M = user.matParam["mobility"]; //Mobility
  double kappa = user.matParam["kappa"];
  double tau = 10./user.N[1];
  
  //Get the second derivative of the free energy
  T f_cc;
  f_cc = F_cc(c.val(0),user.matParam["alpha"],user.matParam["c_a"],user.matParam["c_b"]);

  //Cahn-Hilliard with Nitche's method - c
  r = ( w1.val(0)*(c.val(0) - c.valP(0))/dt )*dV;  
  r += M*w1.grad(0)*f_cc*c.grad(0)*dV;
  r += M*kappa*w1.laplacian(0)*c.laplacian(0)*dV;
  
  r += -w1.val(0)*jn*dS; //Boundary flux
  r += -M*kappa*( c.laplacian(0)*(w1.grad(0)*normal) + w1.laplacian(0)*(c.grad(0)*normal) )*dS;
  r += tau*(w1.grad(0)*normal)*(c.grad(0)*normal)*dS;
  
} //end residual

#include "userFunctionsInstantiation.h"
