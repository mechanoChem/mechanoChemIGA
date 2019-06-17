#include <iostream>
#include <fstream>
#include <mpi.h>
#include <random>
#include "userFunctions.h"
#include "bsplines.h"

template<unsigned int dim>
double userScalarInitialConditions(const Tensor<1,dim,double> &x, unsigned int scalar_i, AppCtx<dim> &user)
{
  std::uniform_real_distribution<> dis(-0.01, 0.01);

  switch(scalar_i) {
  case 0: return 0.707106781186548;
  case 1: return dis(user.gen);
  case 2: return 0.707106781186548;
  case 3: return dis(user.gen);
  default: return 0.;
  }

} //end scalarInitialConditions

template<typename T>
T dgenM(std::vector<T> x){

  //return 4.*x[0]*(std::sqrt(2.) - x[0])*(0.5 - x[1]*x[1]);
  return 16.*x[0]*(std::sqrt(2.) - x[0])*(0.5 - x[1]*x[1])*x[0]*(std::sqrt(2.) - x[0])*(0.5 - x[1]*x[1]);

}

template<typename T>
T M_c(std::vector<T> x){

  //return 4.*(std::sqrt(2.) - 2.*x[0])*(0.5 - x[1]*x[1]);
  return 32.*x[0]*(std::sqrt(2.) - x[0])*(std::sqrt(2.) - 2.*x[0])*(0.5 - x[1]*x[1])*(0.5 - x[1]*x[1]);

}

template<typename T>
T M_eta(std::vector<T> x){

  //return -4.*x[0]*(std::sqrt(2.) - x[0])*(2.*x[1]);
  return -32.*x[0]*(std::sqrt(2.) - x[0])*x[0]*(std::sqrt(2.) - x[0])*(0.5 - x[1]*x[1])*(2.*x[1]);

}


template<unsigned int dim>
void defineParameters(AppCtx<dim>& user){
 
  user.N[0] = 100;
  user.N[1] = 100;
  user.L[0] = 1.;
  user.L[1] = 1.;

  user.dtVal = 0.1;
  user.totalTime = 1e10;
  user.RESTART_IT = 0;
  user.RESTART_TIME = 0.;
  user.skipOutput = 5;

  //Set the domain to be periodic in the x and y directions
  user.periodic[0] = PETSC_TRUE;
  user.periodic[1] = PETSC_TRUE;

  user.scalarSolnFields.push_back("c");
  user.scalarSolnFields.push_back("eta");
  user.scalarSolnFields.push_back("c_b");
  user.scalarSolnFields.push_back("eta_b");

  user.polyOrder = 2;
  user.globalContinuity = 1;

  user.matParam["mobility"] = 0.5;
  user.matParam["kinetic"] = 1.;
  user.matParam["sqrtkappa1"] = 0.005;
  user.matParam["sqrtkappa2"] = 0.005;

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  user.gen.seed(myrank); //For reproducible results

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
  double jn = 0;
  double M = user.matParam["mobility"], L = user.matParam["kinetic"];//10.; //Mobility, kinetic coefficient
  double kappa1 = std::pow(user.matParam["sqrtkappa1"],2), kappa2 = std::pow(user.matParam["sqrtkappa2"],2);
  double tau = 100./user.N[1];
  
  //Get chemical potential and derivatives
  T f_cc_a, f_ceta_a, f_eta_a;
  T f_cc_b, f_ceta_b, f_eta_b;
  T m_a, m_c_a, m_eta_a;
  T m_b, m_c_b, m_eta_b;

  bsplines *freeEn = (bsplines *)user.parameters;
  T y;
  std::vector<T> features(2), dy_dx;
  std::vector<std::vector<T> > ddy_dxdx;

  //Implicit RK4 (2-stage)
  double a11 = 0.25,
    a12 = (3. - 2.*std::sqrt(3.))/12.,
    a21 = (3. + 2.*std::sqrt(3.))/12.,
    a22 = 0.25;

  //Get values for stage a
  features[0] = c.val(0); features[1] = c.val(1)*c.val(1); //We have eta^2 here for symmetry
  freeEn->eval(features,y,dy_dx,ddy_dxdx);
  features[1] = c.val(1); //The mobility just takes in eta, not eta^2
  f_cc_a = ddy_dxdx[0][0];
  f_ceta_a = ddy_dxdx[0][1]*2.*c.val(1);
  f_eta_a = dy_dx[1]*2.*c.val(1);
  m_a = dgenM<T>(features);
  m_c_a = M_c<T>(features);
  m_eta_a = M_eta<T>(features);

  //Get values for stage b
  features[0] = c.val(2); features[1] = c.val(3)*c.val(3); //We have eta^2 here for symmetry
  freeEn->eval(features,y,dy_dx,ddy_dxdx);
  features[1] = c.val(3); //The mobility just takes in eta, not eta^2
  f_cc_b = ddy_dxdx[0][0];
  f_ceta_b = ddy_dxdx[0][1]*2.*c.val(3);
  f_eta_b = dy_dx[1]*2.*c.val(3);
  m_b = dgenM<T>(features);
  m_c_b = M_c<T>(features);
  m_eta_b = M_eta<T>(features);
  /*
  if (c.val(0) < 0. || c.val(0) > std::sqrt(2.) || c.val(1) < -std::sqrt(0.5) || c.val(1) > std::sqrt(0.5)){
    m_a = m_c_a = m_eta_a = 0.;
  }
  if (c.val(2) < 0. || c.val(2) > std::sqrt(2.) || c.val(3) < -std::sqrt(0.5) || c.val(3) > std::sqrt(0.5)){
    m_b = m_c_b = m_eta_b = 0.;
  }
  */

  //Equations 1-2
  //Cahn-Hilliard (periodic BCs)
  r = ( w1.val(0)*(c.val(0) - c.valP(0))/dt )*dV;  
  r += a11*(4.*M*kappa1)*( m_a*f_cc_a*(w1.grad(0)*c.grad(0)) + m_a*f_ceta_a*(w1.grad(0)*c.grad(1)))*dV;
  r += a11*(4.*M*kappa1)*( m_c_a*(w1.grad(0)*c.grad(0)) + m_eta_a*(w1.grad(0)*c.grad(1)))*kappa1*c.laplacian(0)*dV;
  r += a11*(4.*M*kappa1)*m_a*kappa1*w1.laplacian(0)*c.laplacian(0)*dV;
  r += a12*(4.*M*kappa1)*( m_b*f_cc_b*(w1.grad(0)*c.grad(2)) + m_b*f_ceta_b*(w1.grad(0)*c.grad(3)))*dV;
  r += a12*(4.*M*kappa1)*( m_c_b*(w1.grad(0)*c.grad(2)) + m_eta_b*(w1.grad(0)*c.grad(3)))*kappa1*c.laplacian(2)*dV;
  r += a12*(4.*M*kappa1)*m_b*kappa1*w1.laplacian(0)*c.laplacian(2)*dV;
  
  //Allen-Cahn with Nitche's method
  r += ( w1.val(1)*(c.val(1) - c.valP(1))/dt )*dV;
  r += a11*(0.25*L)*( w1.val(1)*m_a*f_eta_a + m_a*kappa2*(w1.grad(1)*c.grad(1)) )*dV;
  r += a11*(0.25*L)*( kappa2*w1.val(1)*(m_c_a*c.grad(0) + m_eta_a*c.grad(1))*c.grad(1) )*dV;
  r += a12*(0.25*L)*( w1.val(1)*m_b*f_eta_b + m_b*kappa2*(w1.grad(1)*c.grad(3)) )*dV;
  r += a12*(0.25*L)*( kappa2*w1.val(1)*(m_c_b*c.grad(2) + m_eta_b*c.grad(3))*c.grad(3) )*dV;
  
  //Equations 3-4
  //Cahn-Hilliard (periodic BCs)
  r += ( w1.val(2)*(c.val(2) - c.valP(0))/dt )*dV;  
  r += a21*(4.*M*kappa1)*( m_a*f_cc_a*(w1.grad(2)*c.grad(0)) + m_a*f_ceta_a*(w1.grad(2)*c.grad(1)))*dV;
  r += a21*(4.*M*kappa1)*( m_c_a*(w1.grad(2)*c.grad(0)) + m_eta_a*(w1.grad(2)*c.grad(1)))*kappa1*c.laplacian(0)*dV;
  r += a21*(4.*M*kappa1)*m_a*kappa1*w1.laplacian(2)*c.laplacian(0)*dV;
  r += a22*(4.*M*kappa1)*( m_b*f_cc_b*(w1.grad(2)*c.grad(2)) + m_b*f_ceta_b*(w1.grad(2)*c.grad(3)))*dV;
  r += a22*(4.*M*kappa1)*( m_c_b*(w1.grad(2)*c.grad(2)) + m_eta_b*(w1.grad(2)*c.grad(3)))*kappa1*c.laplacian(2)*dV;
  r += a22*(4.*M*kappa1)*m_b*kappa1*w1.laplacian(2)*c.laplacian(2)*dV;
  
  //Allen-Cahn with Nitche's method
  r += ( w1.val(3)*(c.val(3) - c.valP(1))/dt )*dV;
  r += a21*(0.25*L)*( w1.val(3)*m_a*f_eta_a + m_a*kappa2*(w1.grad(3)*c.grad(1)) )*dV;
  r += a21*(0.25*L)*( kappa2*w1.val(3)*(m_c_a*c.grad(0) + m_eta_a*c.grad(1))*c.grad(1) )*dV;
  r += a22*(0.25*L)*( w1.val(3)*m_b*f_eta_b + m_b*kappa2*(w1.grad(3)*c.grad(3)) )*dV;
  r += a22*(0.25*L)*( kappa2*w1.val(3)*(m_c_b*c.grad(2) + m_eta_b*c.grad(3))*c.grad(3) )*dV;
  
} //end residual

#include "userFunctionsInstantiation.h"
