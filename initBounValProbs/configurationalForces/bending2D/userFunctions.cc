#include "userFunctions.h"

template<unsigned int dim>
void userBoundaryConditions(AppCtx<dim>& user, double scale){

  //IGASetBoundaryValue(user.iga,"axis","side","dof","val")

  double dVal = scale*user.uDirichlet;
  
  //Plane strain bending 
  //Configurational displacements
  IGASetBoundaryValue(user.iga,0,0,0,0.0);
  IGASetBoundaryValue(user.iga,0,0,1,0.0);

  IGASetBoundaryValue(user.iga,0,1,0,0.0);
  IGASetBoundaryValue(user.iga,0,1,1,-0.5*dVal);

  //Total displacements
  IGASetBoundaryValue(user.iga,0,0,2,0.0);
  IGASetBoundaryValue(user.iga,0,0,3,0.0);

  IGASetBoundaryValue(user.iga,0,1,2,0.0);
  IGASetBoundaryValue(user.iga,0,1,3,-dVal);

} //end boundaryConditions

template<unsigned int dim>
void userLoadStep(unsigned int step,AppCtx<dim> &user){

  //Set scale to be the current time
  double scale = user.time;
  user.boundaryConditions(user,scale+user.dt);

} //end loadStep

template<unsigned int dim>
void userProjectFields(const Tensor<1,dim,double> &x,
		   const Tensor<1,dim,double> &normal,
		   const solutionScalars<dim,double> &c,
		   const solutionVectors<dim,double> &u,
		   AppCtx<dim> &user,
		   std::vector<double> &scalarProjections,
		   std::vector<Tensor<1,dim,double> > &vectorProjections){
  
  //Project 1 scalar field
  Tensor<2,dim,double> chi, Theta, eye;
  identity(eye);
  chi = u.grad(0) + eye;
  Theta = 0.5*(trans(chi)*chi - eye);
  scalarProjections[0] = Theta[0][0] - Theta[1][1]; //eta2
  
} //end projectFields

template<unsigned int dim,typename T>
void constitutiveModel(Tensor<3,dim,T> &B,
		       Tensor<2,dim,T> &D,
		       Tensor<2,dim,T> &P,
		       Tensor<2,dim,T> &dpsiS_dchi,
		       Tensor<2,dim,T> &Eshelby,
		       Tensor<2,dim,T> &chi,
		       T &Jchi,
		       const solutionVectors<dim,T> &u){
  
  //Define constitutive model parameters
  double mu = 0.1,
    alphaC = 0.2,
    betaC = 0.1;

  Tensor<3,dim,T> dchi;
  Tensor<2,dim,T> Phi, Theta, eye;
  Tensor<2,dim,T> F, E, S;
  Tensor<1,dim,T> Lambda;
  T psiS;

  //Kinematics
  identity(eye);
  chi = u.grad(0) + eye;
  dchi = u.hess(0);
  Jchi = det(chi);
  Phi = trans(chi)*chi;
  Theta = 0.5*(Phi - eye);
  F = (u.grad(1) + eye)*inv(chi);
  E = 0.5*(trans(F)*F - eye);

  //Stress
  S = betaC*trace(E)*eye + 2*mu*E;
  psiS = 0.5*betaC*trace(E)*trace(E) + mu*double_contract(E,E);
  for (unsigned int i=0; i<dim; ++i){
    Lambda[i] = std::sqrt(Phi[i][i]);
    S[i][i] += (alphaC*Lambda[i] - betaC - 2*mu)*E[i][i];
    psiS += 0.5*(alphaC*Lambda[i] - betaC - 2*mu)*E[i][i]*E[i][i];
    for (unsigned int j=0; j<dim; ++j){
      dpsiS_dchi[j][i] = 0.5*alphaC/Lambda[i]*chi[j][i]*E[i][i]*E[i][i];
    }
  }
  P = F*S;
  Eshelby = psiS*eye - trans(F)*P;

  //The definitions below are only for 2D...
  double Es = 0.1,
    Ed = 1.0,
    El = 0.1;
  T eta1 = Theta[0][0] + Theta[1][1];
  T eta2 = Theta[0][0] - Theta[1][1];
  T eta6 = Theta[0][1];

  T eta2_1=0.0, eta2_2=0.0; 
  for (unsigned int i=0; i<dim; ++i){
    eta2_1+=(chi[i][0]*dchi[i][0][0]-chi[i][1]*dchi[i][1][0]);
    eta2_2+=(chi[i][0]*dchi[i][0][1]-chi[i][1]*dchi[i][1][1]);
  }

  Tensor<2,dim,T> deta1_dchi, deta2_dchi, deta6_dchi;
  Tensor<2,dim,T> deta2_1_dchi, deta2_2_dchi;
  Tensor<3,dim,T> deta2_1_ddchi, deta2_2_ddchi;
  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int J=0; J<dim; ++J){
      deta1_dchi[i][J] = (chi[i][0]*(0==J)+chi[i][1]*(1==J));
      deta2_dchi[i][J] = (chi[i][0]*(0==J)-chi[i][1]*(1==J));
      deta6_dchi[i][J] = (chi[i][1]*(0==J)+chi[i][0]*(1==J))/2.0;
      deta2_1_dchi[i][J] = ((0==J)*dchi[i][0][0]-(1==J)*dchi[i][1][0]);
      deta2_2_dchi[i][J] = ((0==J)*dchi[i][0][1]-(1==J)*dchi[i][1][1]);

      //gradient terms
      for (unsigned int K=0; K<dim; ++K){
	deta2_1_ddchi[i][J][K] = (chi[i][0]*(0==J)-chi[i][1]*(1==J))*(0==K);
	deta2_2_ddchi[i][J][K] = (chi[i][0]*(0==J)-chi[i][1]*(1==J))*(1==K);
      }
    }
  } 
  
  D = 2.*Ed/(Es*Es)*(eta1*deta1_dchi + 2.*(eta2*eta2/(Es*Es) - 1.)*eta2*deta2_dchi + eta6*deta6_dchi);
  D = D + 2*El*El*Ed/(Es*Es)*(eta2_1*deta2_1_dchi + eta2_2*deta2_2_dchi);
  B = 2*El*El*Ed/(Es*Es)*(eta2_1*deta2_1_ddchi + eta2_2*deta2_2_ddchi);
    
} //end constitutiveModel

template<unsigned int dim>
void defineParameters(AppCtx<dim>& user){
 
  user.N[0] = 160;
  user.N[1] = 16;
  user.L[0] = 10.;
  user.L[1] = 1.;

  user.matParam["u_applied"] = 0.196;
  user.uDirichlet = user.matParam["u_applied"];

  user.vectorSolnFields.push_back("configDisplacement");
  user.vectorSolnFields.push_back("totalDisplacement");
  user.scalarProjectnFields.push_back("eta2");

  user.polyOrder = 2;
  user.globalContinuity = 1;

  user.dtVal = 0.005;
  user.totalTime = 1.;
  user.RESTART_IT = 0;
  user.RESTART_TIME = 0;
  
  user.boundaryConditions = userBoundaryConditions;
  user.projectFields = userProjectFields;
  user.loadStep = userLoadStep;

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
  
  Tensor<3,dim,T> B;
  Tensor<2,dim,T> D,P,dpsiS_dchi,Eshelby,chi;
  T Jchi;
  constitutiveModel(B,D,P,dpsiS_dchi,Eshelby,chi,Jchi,u);

  //Weak form
  r = triple_contract(w2.hess(0),B)*dV; //Material
  r += double_contract(w2.grad(0),D+Jchi*(Eshelby*trans(inv(chi)) + dpsiS_dchi))*dV; //Material
  r += double_contract(w2.grad(1),Jchi*P*trans(inv(chi)))*dV; //Standard

} //end residual

#include "userFunctionsInstantiation.h"
