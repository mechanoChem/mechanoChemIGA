#include "userFunctions.h"

template<unsigned int dim>
void userBoundaryConditions(AppCtx<dim>& user, double scale){

  //IGASetBoundaryValue(user.iga,"axis","side","dof","val")
  IGASetBoundaryValue(user.iga,2,0,0,0.0); 
  IGASetBoundaryValue(user.iga,2,0,1,0.0);
  IGASetBoundaryValue(user.iga,2,0,2,0.0);

} //end boundaryConditions

template<unsigned int dim>
void userProjectFields(const Tensor<1,dim,double> &x,
		   const Tensor<1,dim,double> &normal,
		   const solutionScalars<dim,double> &c,
		   const solutionVectors<dim,double> &u,
		   AppCtx<dim> &user,
		   std::vector<double> &scalarProjections,
		   std::vector<Tensor<1,dim,double> > &vectorProjections){
  
  //Project 1 scalar field and 1 vector fields
  Tensor<2,dim,double> P = get1stPiolaKirchhoff(x,c,u,user,0);
  scalarProjections[0] = std::pow(P[0][0]-P[1][1],2) + std::pow(P[1][1]-P[2][2],2) + std::pow(P[2][2]-P[0][0],2);
  scalarProjections[0] += 6.*(P[0][1]*P[0][1] + P[1][2]*P[1][2] + P[2][0]*P[2][0]);
  scalarProjections[0] = std::sqrt(0.5*scalarProjections[0]);
  
} //end projectFields

template<unsigned int dim,typename T>
Tensor<2,dim,T> get1stPiolaKirchhoff(const Tensor<1,dim,double> &x,
				     const solutionScalars<dim,T> &c,
				     const solutionVectors<dim,T> &u,
				     AppCtx<dim> &user,
				     unsigned int index){
  
  Tensor<2,dim,T> F, E, S, eye;
  identity(eye);
  F = u.grad(index) + eye;
  E = 0.5*(trans(F)*F - eye);
  S = double_contract(user.C_e,E); //S = C_e:E
    
  return F*S;
} //end get1stPiolaKirchhoff

template<unsigned int dim>
void defineParameters(AppCtx<dim>& user){
 
  user.N[0] = 10;
  user.N[1] = 10;
  user.N[2] = 10;
  user.L[0] = 1.;
  user.L[1] = 1.;
  user.L[2] = 1.;

  user.vectorSolnFields.push_back("displacement");
  user.scalarProjectnFields.push_back("vonMises");

  user.polyOrder = 1;
  user.globalContinuity = 0;

  //Define elasticity tensor
  user.matParam["E"] = 2.0e11;
  user.matParam["nu"] = 0.3;
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
  
  // Neumann condition scalar
  user.matParam["h"] = 1.e11;
  user.matParam["h0"] = 0.;

  user.boundaryConditions = userBoundaryConditions;
  user.projectFields = userProjectFields;

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
  
  //Elasticity
  Tensor<2,dim,T> P = get1stPiolaKirchhoff(x,c,u,user,0);
  Tensor<1,dim,T> h;
  h[2] = user.matParam["h"]*x[0]*(x[2] == user.L[2]) + user.matParam["h0"];

  // \int_\Omega (grad{w}:P) dV - \int_{\partial\Omega} (w\cdot h) dS
  r = double_contract(w2.grad(0),P)*dV;
  r += -(w2.val(0)*h)*dS; //Traction

} //end residual

#include "userFunctionsInstantiation.h"
