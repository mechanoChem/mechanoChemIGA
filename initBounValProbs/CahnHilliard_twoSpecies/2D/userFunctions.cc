#include "userFunctions.h"

template<unsigned int dim>
double userScalarInitialConditions(const Tensor<1,dim,double> &x, unsigned int scalar_i, AppCtx<dim> &user)
{

  switch(scalar_i) {
  case 0: return 0.5 + 0.1*(0.5 - (double)(rand() % 100 )/100.0); //Random about 0.5
  case 1: return 0.5 + 0.3*(0.5 - (double)(rand() % 100 )/100.0); //Random about 0.5
  default: return 0.;
  }

} //end scalarInitialConditions

template<typename T>
T F_c1c1(T c1, T c2){
  return 2*c2*(std::pow(2.*c1-1.,2) + 2.*(c1-0.1)*(c1-.9));
}

template<typename T>
T F_c1c2(T c1, T c2){
  return 2.*(c1-0.1)*std::pow(c1-0.9,2) + 2.*std::pow(c1-0.1,2)*(c1-0.9);
}

template<typename T>
T F_c2c2(T c1, T c2){
  return 2.*((c2-0.8)*(2.*c2-1.)+(c2-0.2)*(2.*c2-1.)+2.*(c2-0.2)*(c2-0.8));
} //end free energy derivatives

template<unsigned int dim>
void defineParameters(AppCtx<dim>& user){
 
  user.N[0] = 20;
  user.N[1] = 20;
  
  user.L[0] = 1.;
  user.L[1] = 1.;

  user.dtVal = .1;
  user.totalTime = 20;
  user.RESTART_IT = 0;
  user.RESTART_TIME = 0.;
  user.skipOutput = 2;

  user.scalarSolnFields.push_back("c1");
  user.scalarSolnFields.push_back("c2");

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
  double jn1 = 0, jn2 = 0;
  double M = .1, L = 2.; //Mobility
  double kappa1 = .0005, kappa2 = .0005;
  double tau = 0.1*(user.N[0]/user.L[0]);
  
  //Get chemical potential and derivatives
  T f_c1c1, f_c1c2, f_c2c2;
  f_c1c1 = F_c1c1(c.val(0),c.val(1));
  f_c1c2 = F_c1c2(c.val(0),c.val(1));
  f_c2c2 = F_c2c2(c.val(0),c.val(1));
  
  //Cahn-Hilliard with Nitche's method - c1
  r = ( w1.val(0)*(c.val(0) - c.valP(0))/dt )*dV;  
  r += ( M*w1.grad(0)*(f_c1c1*c.grad(0) + f_c1c2*c.grad(1)))*dV;
  r += M*kappa1*w1.laplacian(0)*c.laplacian(0)*dV;
  
  r += -w1.val(0)*jn1*dS; //Boundary flux
  r += -M*kappa1*( c.laplacian(0)*(w1.grad(0)*normal) + w1.laplacian(0)*(c.grad(0)*normal) )*dS;
  r += tau*(w1.grad(0)*normal)*(c.grad(0)*normal)*dS;
  
  //Cahn-Hilliard with Nitche's method - c2
  r += ( w1.val(1)*(c.val(1) - c.valP(1))/dt )*dV;  
  r += ( M*w1.grad(1)*(f_c1c2*c.grad(0) + f_c2c2*c.grad(1)))*dV;
  r += M*kappa2*w1.laplacian(1)*c.laplacian(1)*dV;
  
  r += -w1.val(1)*jn2*dS; //Boundary flux
  r += -M*kappa2*( c.laplacian(1)*(w1.grad(1)*normal) + w1.laplacian(1)*(c.grad(1)*normal) )*dS;
  r += tau*(w1.grad(1)*normal)*(c.grad(1)*normal)*dS;
  
} //end residual

#include "userFunctionsInstantiation.h"
