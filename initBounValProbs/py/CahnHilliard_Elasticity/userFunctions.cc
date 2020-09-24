#include "userFunctions.h"
#include<iostream>

template<unsigned int dim>
double userScalarInitialConditions(const Tensor<1,dim,double> &x, unsigned int scalar_i, AppCtx<dim> &user)
{

  //  return user.matParam["c_avg"] + 0.1*(0.5 - (double)(rand() % 100 )/100.0); //Random about 0.5
  return user.matParam["c_avg"]*x[0];

} //end scalarInitialConditions

template<typename T>
T fcc(T c, double alpha, double ca, double cb){
  
  //Second derivative of the free energy density: f(c) = alpha*(c-ca)^2*(c-cb)^2
  return 2*alpha*(std::pow(2.*c-ca-cb,2) + 2.*(c-ca)*(c-cb));
  
} //end fcc

template<unsigned int dim,typename T,typename U>
auto custom_contract(const Tensor<2,dim,T> a, const Tensor<3,dim,U> &b) -> Tensor<1,dim,decltype(std::declval<T>()*std::declval<U>())>{

  Tensor<1,dim,decltype(std::declval<T>()*std::declval<U>())> tmp;
  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      for (unsigned int k=0; k<dim; ++k){
	tmp[k] += a[i][j]*b[i][k][j];
      }
    }
  }

  return tmp;

}

template<unsigned int dim>
void userPostParameters(AppCtx<dim> &user)
{

  double E = user.matParam["E"], nu = user.matParam["nu"];
  double lambda = (E*nu)/((1.+nu)*(1.-2.*nu));
  double mu = E/(2.*(1.+nu));
  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      for (unsigned int k=0; k<dim; ++k){
	for (unsigned int l=0; l<dim; ++l){
	  user.C_e[i][j][k][l] = lambda*(i==j)*(k==l) + 2.*0.5*mu*((i==k)*(j==l) + (i==l)*(j==k));
	}
      }
    }
  } //end C_e

  if(user.CahnHilliard || user.Diffusion){
    if(user.scalarSolnFields.size() == 0){
      user.scalarSolnFields.push_back("c");
    }
  }
  if(user.Elasticity){
    if(user.vectorSolnFields.size() == 0){
      user.vectorSolnFields.push_back("displacement");
    }
  }

}

template<unsigned int dim>
void defineParameters(AppCtx<dim>& user){
 
  user.N[0] = 100;
  user.L[0] = 1.;
  if (dim >= 2){
    user.N[1] = 100;
    user.L[1] = 1.;
    if (dim == 3){
      user.N[2] = 10;
      user.L[2] = 0.1;
    }
  }

  //Set the domain to be periodic in all directions
  for (unsigned int i=0; i<dim; ++i){
    user.periodic[i] = PETSC_TRUE;
  }

  //Define some material parameters (can be overwritten by parameters file)
  user.matParam["mobility"] = .1; //Mobility
  user.matParam["kappa"] = .0001; //Gradient energy parameter
  user.matParam["D"] = 20; //diffusivity

  user.matParam["c_avg"] = 0.5;

  user.matParam["Feiga_11"] = 1;
  user.matParam["Feiga_22"] = 1;
  user.matParam["Feiga_33"] = 1;

  user.matParam["Feigb_11"] = 1;//.05;
  user.matParam["Feigb_22"] = 1;//.02;
  user.matParam["Feigb_33"] = 1;

  //Define some free energy parameters
  user.matParam["alpha"] = 0.25; //Free energy coefficient
  user.matParam["c_a"] = 0.2; //Composition of phase a
  user.matParam["c_b"] = 0.9; //Composition of phase b

  //Define elasticity tensor
  user.matParam["E"] = 2.;//0e11;
  user.matParam["nu"] = 0.3;

  user.dtVal = .1;
  user.totalTime = 20;
  user.RESTART_IT = 0;
  user.RESTART_TIME = 0.;
  user.skipOutput = 2;

  user.polyOrder = 2;
  user.globalContinuity = 1;

  user.postParameters = userPostParameters;
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
  double M = user.matParam["mobility"]; //Mobility
  double kappa = user.matParam["kappa"];
  
  //Get the second derivative of the free energy
  T f_cc;
  if(user.CahnHilliard){
    f_cc = fcc(c.val(0),user.matParam["alpha"],user.matParam["c_a"],user.matParam["c_b"]);
  }

  Tensor<2,dim,T> F, Fe, Feig, Feiga, Feigb, invFeig, Ee, E, S, eye;
  Tensor<2,dim,T> E_c, E_cc;
  T psi_cc;
  identity(eye);

  // Deformation gradient, stress, etc.
  if(user.Elasticity){
    F = u.grad(0) + eye;
    if(user.CahnHilliard || user.Diffusion){
      Feiga[0][0] = user.matParam["Feiga_11"];
      Feigb[0][0] = user.matParam["Feigb_11"];
      if (dim >= 2){
	Feiga[1][1] = user.matParam["Feiga_22"];
	Feigb[1][1] = user.matParam["Feigb_22"];
	if (dim == 3){
	  Feiga[2][2] = user.matParam["Feiga_33"];
	  Feigb[2][2] = user.matParam["Feigb_33"];
	}
      }
      //Eigenstrain is linearly interpolated
      Feig = (c.val(0) - user.matParam["c_a"])/(user.matParam["c_b"] - user.matParam["c_a"])*(Feigb - Feiga) + Feiga;
      invFeig = inv(Feig);
      Fe = F*invFeig;
      Ee = 0.5*(trans(Fe)*Fe - eye);
      S = double_contract(user.C_e,Ee); //S = C_e:Ee


      E_c = -1./(user.matParam["c_b"] - user.matParam["c_a"])*(trans(Fe)*Fe*(Feigb - Feiga)*invFeig);
      E_c = 0.5*(E_c + trans(E_c));
      E_cc = 1./pow(user.matParam["c_b"] - user.matParam["c_a"],2)*(trans(invFeig)*(trans(Feigb) - trans(Feiga))*
								    trans(Fe)*Fe*(Feigb - Feiga)*invFeig +
								    2.*trans(Fe)*Fe*(Feigb - Feiga)*invFeig*(Feigb - Feiga)*invFeig);
      E_cc = 0.5*(E_cc + trans(E_cc));

      psi_cc = full_contract(E_c,user.C_e,E_c) + double_contract(S,E_cc);
    }
    else{
      E = 0.5*(trans(F)*F - eye);
      S = double_contract(user.C_e,E);
    }
  }
  double tau = 0.1*(user.N[0]/user.L[0]);

  //Cahn-Hilliard - c
  r = 0;
  if(user.CahnHilliard || user.Diffusion){
    r += ( w1.val(0)*(c.val(0) - c.valP(0))/dt )*dV;
    if (user.CahnHilliard){
      r += M*w1.grad(0)*f_cc*c.grad(0)*dV;
      r += M*kappa*w1.laplacian(0)*c.laplacian(0)*dV;

      //r += -w1.val(0)*jn*dS; //Boundary flux

      //Nitche's method for higher order BC
      r += -M*kappa*( c.laplacian(0)*(w1.grad(0)*normal) + w1.laplacian(0)*(c.grad(0)*normal) )*dS;
      r += tau*(w1.grad(0)*normal)*(c.grad(0)*normal)*dS;
    }
    else{ //Diffusion
      r += user.matParam["D"]*w1.grad(0)*c.grad(0)*dV;
    }
    if(user.Elasticity){
      r += M*w1.grad(0)*psi_cc*c.grad(0)*dV;
      r += M*w1.grad(0)*(custom_contract(Fe*double_contract(user.C_e,E_c)*trans(invFeig),u.hess(0)) + 
			 -1./(user.matParam["c_b"] - user.matParam["c_a"])*
			 custom_contract(Fe*((Feigb - Feiga)*invFeig*S*trans(invFeig) +
					     S*trans(invFeig)*trans(Feigb - Feiga))*trans(invFeig),u.hess(0)))*dV;
      // Elasticity w/ Cahn Hilliard: \int_\Omega (grad{w}:(P*Feig^{-T})) dV
      r += double_contract(w2.grad(0),Fe*S*trans(invFeig))*dV;
    }
  }

  // Elasticity only: \int_\Omega (grad{w}:P) dV
  else if(user.Elasticity){
    r += double_contract(w2.grad(0),F*S)*dV;
  }


} //end residual

#include "userFunctionsInstantiation.h"
