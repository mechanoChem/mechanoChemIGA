#include "userFunctions.h"
#include "json.hpp"
#include <iostream>
#include <map>
#include <string>

template<unsigned int dim>
double userScalarInitialConditions(const Tensor<1,dim,double> &x, unsigned int scalar_i, AppCtx<dim> &user)
{
  nlohmann::json ic = user.param["initial_conditions"];
  double cons = ic["c_avg"].get<double>() - 0.5*ic["c_slope_x"].get<double>()*user.L[0];
  if (dim > 1){
    cons -=  0.5*ic["c_slope_y"].get<double>()*user.L[1];
    if (dim == 3){
      cons -= 0.5*ic["c_slope_z"].get<double>()*user.L[2];
    }
  }
  double val = ic["c_slope_x"].get<double>()*x[0] + cons;
  if (dim > 1){
    val += ic["c_slope_y"].get<double>()*x[1];
    if (dim == 3){
      val += ic["c_slope_z"].get<double>()*x[2];
    }
  }
  val += ic["random_perturb"].get<double>()*(2*(0.5 - (double)(rand() % 100 )/100.0));
  return val;

} //end scalarInitialConditions


template<unsigned int dim>
void userBoundaryConditions(AppCtx<dim>& user, double scale){

  unsigned int nScalars = user.scalarSolnFields.size();

  std::map<std::string,int> bndry2axis = {
    {"x_min",0},
    {"x_max",0},
    {"y_min",1},
    {"y_max",1},
    {"z_min",2},
    {"z_max",2}
  };
  std::map<std::string,int> bndry2side = {
    {"x_min",0},
    {"x_max",1},
    {"y_min",0},
    {"y_max",1},
    {"z_min",0},
    {"z_max",1}
  };
  std::map<std::string,int> field2dof = {
    {"comp",0},
    {"u_x",nScalars},
    {"u_y",nScalars + 1},
    {"u_z",nScalars + 2}
  };

  nlohmann::json bc;
  if (user.param.contains("boundary_conditions")){
    bc = user.param["boundary_conditions"];

    // Loop over boundary conditions
    for (auto field = bc.begin(); field != bc.end(); ++field){
      for (auto bndry = field.value().begin(); bndry != field.value().end(); ++bndry){
	if (bndry.value()["type"].get<std::string>() == "Dirichlet"){
	  //IGASetBoundaryValue(user.iga,"axis","side","dof","val")
	  IGASetBoundaryValue(user.iga,
			      bndry2axis[bndry.key()],
			      bndry2side[bndry.key()],
			      field2dof[field.key()],
			      bndry.value()["value"].get<double>()); 
	}
      }
    }
  }

} //end boundaryConditions

template<typename T>
T fcc(T c, double omega, double ca, double cb){
  
  //Second derivative of the free energy density: f(c) = omega*(c-ca)^2*(c-cb)^2
  return 2*omega*(std::pow(2.*c-ca-cb,2) + 2.*(c-ca)*(c-cb));
  
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
  double E, nu;
  user.param["material_parameters"]["youngs_modulus"].get_to(E);
  user.param["material_parameters"]["poissons_ratio"].get_to(nu);
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

  // Read in which physics will be solved
  if (user.param.contains("physics")){
    user.Elasticity = false;
    user.Diffusion = false;
    user.CahnHilliard = false;
    for (unsigned int i=0; i<user.param["physics"].size(); ++i){
      if(user.param["physics"][i] == "elasticity"){
	user.Elasticity = true;
      }
      else if(user.param["physics"][i] == "diffusion"){
	user.Diffusion = true;
      }
      else if(user.param["physics"][i] == "cahn_hilliard"){
	user.CahnHilliard = true;
      }
    }
  }

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
  //for (unsigned int i=0; i<dim; ++i){
  //  user.periodic[i] = PETSC_TRUE;
  //}

  //Define some material parameters (can be overwritten by parameters file)
  user.param["material_parameters"]["mobility"] = .1; //Mobility
  user.param["material_parameters"]["kappa"] = .0001; //Gradient energy parameter
  user.param["material_parameters"]["diffusivity"] = 20; //diffusivity

  user.param["material_parameters"]["Feiga_11"] = 1;
  user.param["material_parameters"]["Feiga_22"] = 1;
  user.param["material_parameters"]["Feiga_33"] = 1;

  user.param["material_parameters"]["Feigb_11"] = 1;//.05;
  user.param["material_parameters"]["Feigb_22"] = 1;//.02;
  user.param["material_parameters"]["Feigb_33"] = 1;

  //Define some free energy parameters
  user.param["material_parameters"]["omega"] = 0.25; //Free energy coefficient
  user.param["material_parameters"]["c_alpha"] = 0.2; //Composition of phase a
  user.param["material_parameters"]["c_beta"] = 0.9; //Composition of phase b

  user.param["initial_conditions"]["c_avg"] = 0.5;
  user.param["initial_conditions"]["c_slope_x"] = 0.;
  user.param["initial_conditions"]["c_slope_y"] = 0.;
  user.param["initial_conditions"]["c_slope_z"] = 0.;
  user.param["initial_conditions"]["random_perturb"] = 0.;

  //Define elasticity tensor
  user.param["material_parameters"]["youngs_modulus"] = 2.;//0e11;
  user.param["material_parameters"]["poissons_ratio"] = 0.3;

  user.dtVal = .1;
  user.totalTime = 20;
  user.RESTART_IT = 0;
  user.RESTART_TIME = 0.;
  user.skipOutput = 2;

  user.polyOrder = 2;
  user.globalContinuity = 1;

  user.postParameters = userPostParameters;
  user.scalarInitialConditions = userScalarInitialConditions;
  user.boundaryConditions = userBoundaryConditions;

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

  nlohmann::json param = user.param["material_parameters"];

  //Chemistry
  double dt = user.dt;
  double M, kappa, c_alpha, c_beta, D;
  param["mobility"].get_to(M); //Mobility
  param["diffusivity"].get_to(D);
  param["kappa"].get_to(kappa);
  param["c_alpha"].get_to(c_alpha);
  param["c_beta"].get_to(c_beta);

  //Get the second derivative of the free energy
  T f_cc;
  if(user.CahnHilliard){
    f_cc = fcc(c.val(0),param["omega"].get<double>(),c_alpha,c_beta);
  }

  Tensor<2,dim,T> F, Fe, Feig, Feiga, Feigb, invFeig, Ee, E, S, eye;
  Tensor<2,dim,T> E_c, E_cc;
  T psi_cc;
  identity(eye);

  // Deformation gradient, stress, etc.
  if(user.Elasticity){
    F = u.grad(0) + eye;
    if(user.CahnHilliard || user.Diffusion){
      for (unsigned int i=0; i<dim; ++i){
	for(unsigned int j=0; j<dim; ++j){
	  Feiga[i][j] = param["F_eig_alpha"][i][j].get<double>();
	  Feigb[i][j] = param["F_eig_beta"][i][j].get<double>();
	}
      }
      //Eigenstrain is linearly interpolated
      Feig = (c.val(0) - c_alpha)/(c_beta - c_alpha)*(Feigb - Feiga) + Feiga;
      invFeig = inv(Feig);
      Fe = F*invFeig;
      Ee = 0.5*(trans(Fe)*Fe - eye);
      S = double_contract(user.C_e,Ee); //S = C_e:Ee


      E_c = -1./(c_beta - c_alpha)*(trans(Fe)*Fe*(Feigb - Feiga)*invFeig);
      E_c = 0.5*(E_c + trans(E_c));
      E_cc = 1./pow(c_beta - c_alpha,2)*(trans(invFeig)*(trans(Feigb) - trans(Feiga))*
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
      r += D*w1.grad(0)*c.grad(0)*dV;
    }
    if(user.Elasticity){
      r += M*w1.grad(0)*psi_cc*c.grad(0)*dV;
      r += M*w1.grad(0)*(custom_contract(Fe*double_contract(user.C_e,E_c)*trans(invFeig),u.hess(0)) + 
			 -1./(c_beta - c_alpha)*
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
