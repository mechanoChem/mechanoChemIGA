//extern "C" {
#include "petiga.h"
//}
#include "constitutive.h"
#include "physicsHeaders.h"
#include "IBVPHeaders.h"
#include "utilsIGAHeaders.h"

//include automatic differentiation library
#include <Sacado.hpp>

//residual function implementation
#undef  __FUNCT__
#define __FUNCT__ "Function"
template <class T,unsigned int DIM, unsigned int DOF>
PetscErrorCode Function(IGAPoint p,PetscReal dt2,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const T * U,
			PetscReal t0,const PetscScalar * U0,
			T *R,void *ctx)
{
  AppCtx *user = (AppCtx *)ctx;

	//Retrieve material parameters
	PetscReal mu = user->matParam["mu"];
	PetscReal betaC = user->matParam["betaC"];
	PetscReal alphaC = user->matParam["alphaC"];
	PetscReal anisoCoeff = user->matParam["anisoCoeff"];

	//Nonconvex free energy parameters
	PetscReal Es = user->matParam["Es"];
	PetscReal Ed = user->matParam["Ed"];
	PetscReal El = user->matParam["El"];

  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  PetscReal *n = p->normal;
  PetscReal c=user->lambda; //used as a load paramter 0<c<1

  //configuration displacement field variable
  T UU[DIM], UUx[DIM][DIM], UUxx[DIM][DIM][DIM];
  computeField<T,DIM,DOF>(VECTOR,0,p,U,&UU[0],&UUx[0][0],&UUxx[0][0][0]);

  //total displacement field variable
  T u[DIM], ux[DIM][DIM], uxx[DIM][DIM][DIM];
  computeField<T,DIM,DOF>(VECTOR,DIM,p,U,&u[0],&ux[0][0],&uxx[0][0][0]);

  //Compute \chi (I+Ux), d\chi (Uxx)
  T chi[DIM][DIM], dchi[DIM][DIM][DIM];
  for (unsigned int i=0; i<DIM; i++) {
    for (unsigned int J=0; J<DIM; J++) {
      chi[i][J]=(i==J)+UUx[i][J];
      for (unsigned int K=0; K<DIM; K++) {
	dchi[i][J][K]=UUxx[i][J][K];
      }
    }
  }

  //Compute J_\chi (the determinant of \chi)
  T J_chi;

  //Compute \chi^{-1}
  T chi_Inv[DIM][DIM];
  if(DIM == 2){
    J_chi = chi[0][0]*chi[1][1] - chi[0][1]*chi[1][0];

    chi_Inv[0][0] = 1./J_chi*chi[1][1];
    chi_Inv[0][1] = -1./J_chi*chi[0][1];
    chi_Inv[1][0] = -1./J_chi*chi[1][0];
    chi_Inv[1][1] = 1./J_chi*chi[0][0];
  }
  else if(DIM == 3){
    J_chi = chi[0][0]*(chi[1][1]*chi[2][2] - chi[1][2]*chi[2][1]) -
      chi[0][1]*(chi[1][0]*chi[2][2] - chi[1][2]*chi[2][0]) +
      chi[0][2]*(chi[1][0]*chi[2][1] - chi[1][1]*chi[2][0]);

    chi_Inv[0][0] = 1./J_chi*(chi[1][1]*chi[2][2] - chi[2][1]*chi[1][2]);
    chi_Inv[0][1] = 1./J_chi*(chi[0][2]*chi[2][1] - chi[0][1]*chi[2][2]);
    chi_Inv[0][2] = 1./J_chi*(chi[0][1]*chi[1][2] - chi[1][1]*chi[0][2]);
    chi_Inv[1][0] = 1./J_chi*(chi[1][2]*chi[2][0] - chi[2][2]*chi[1][0]);
    chi_Inv[1][1] = 1./J_chi*(chi[0][0]*chi[2][2] - chi[2][0]*chi[0][2]);
    chi_Inv[1][2] = 1./J_chi*(chi[0][2]*chi[1][0] - chi[1][2]*chi[0][0]);
    chi_Inv[2][0] = 1./J_chi*(chi[1][0]*chi[2][1] - chi[2][0]*chi[1][1]);
    chi_Inv[2][1] = 1./J_chi*(chi[0][1]*chi[2][0] - chi[2][1]*chi[0][0]);
    chi_Inv[2][2] = 1./J_chi*(chi[0][0]*chi[1][1] - chi[1][0]*chi[0][1]);
  }

  //Check inverse
  T check[DIM][DIM];
  for (unsigned int I=0; I<DIM; I++){
    for (unsigned int J=0; J<DIM; J++){
      check[I][J] = 0.;
      for (unsigned int k=0; k<DIM; k++){
	check[I][J] += chi[I][k]*chi_Inv[k][J];
      }
      if(std::abs(check[I][J] - (I==J)) > 2.0e-15){
	PetscPrintf(PETSC_COMM_WORLD,"inverse is wrong...\n");
      }
    }
  }

  //Compute \Phi=\chi^T*\chi and strain metric, \Xi=0.5*(\Phi-I)
  T Phi[DIM][DIM], Xi[DIM][DIM];
  for (unsigned int I=0; I<DIM; I++){
    for (unsigned int J=0; J<DIM; J++){
      Phi[I][J] = 0.;
      for (unsigned int k=0; k<DIM; k++){
	Phi[I][J] += chi[k][I]*chi[k][J];
      }
      Xi[I][J] = 0.5*(Phi[I][J] - (I==J));
    }
  }

  //Compute normal stretches \Lambda_I = \sqrt{\Phi_{II}}
  T Lambda[DIM];
  for (unsigned int I=0; I<DIM; I++){
    Lambda[I] = sqrt(Phi[I][I]);
  }
  
  //Compute F (I+ux), dF (uxx)
  T F[DIM][DIM], dF[DIM][DIM][DIM];
  for (unsigned int i=0; i<DIM; i++) {
    for (unsigned int J=0; J<DIM; J++) {
      F[i][J] = 0.;
      for (unsigned int k=0; k<DIM; k++){
      	//F[i][J] += ((i==k) + UUx[i][k] + ux[i][k])*chi_Inv[k][J];
	F[i][J] += ((i==k) + ux[i][k])*chi_Inv[k][J]; //let u be the total deformation
      }
    }
  }

  //Compute strain metric, E=0.5*(F^T*F-I)
  T E[DIM][DIM];
  for (unsigned int I=0; I<DIM; I++){
    for (unsigned int J=0; J<DIM; J++){
      E[I][J] = -0.5*(I==J);
      for (unsigned int k=0; k<DIM; k++){
	E[I][J] += 0.5*F[k][I]*F[k][J];
      }
    }
  }
  
  //compute P and Beta
  T P[DIM][DIM];
  T P0[DIM][DIM], Beta0[DIM][DIM][DIM];

  //define alpha and beta tensors
  T alpha[DIM], beta[DIM][DIM];
  for (unsigned int I=0; I<DIM; I++){
    alpha[I] = anisoCoeff*alphaC*(Lambda[I] - (anisoCoeff - 1.)/anisoCoeff);
    for (unsigned int J=0; J<DIM; J++){
      //      beta[I][J] = betaC*Lambda[I]*Lambda[J];
      beta[I][J] = betaC;
    }
  }
  
  //Compute \psi_N (Newtonian strain energy density)
  T psi_N;
  psi_N = 0.;
  for (unsigned int I=0; I<DIM; I++){
    psi_N += 0.5*(alpha[I] - 2*mu)*E[I][I]*E[I][I];
    for (unsigned int J=0; J<DIM; J++){
      if(I != J){
	psi_N += 0.5*beta[I][J]*E[I][I]*E[J][J];
      }
      psi_N += mu*E[I][J]*E[I][J];
    }
  }
  
  //Compute \partial\psi_N/\partial\chi (partial of \psi_N with respect to \chi)
  T dpsi_dchi[DIM][DIM];
  for(unsigned int K=0; K<DIM; K++){
    for(unsigned int L=0; L<DIM; L++){
      dpsi_dchi[K][L] = 0.5*anisoCoeff*(alphaC/Lambda[L])*E[L][L]*E[L][L]*chi[K][L]; //stiffens with elongation
      for (unsigned int I=0; I<DIM; I++){
	if(I != L){
	  //	dpsi_dchi[K][L] += 0.5*betaC*chi[K][L]*E[L][L]*E[I][I]*
	  //										(Lambda[I]/Lambda[L] + Lambda[L]/Lambda[I]);
	}
      }
    }
  }
 
  //2D model
	if(DIM==2){
		//new strain metrics
		T e1=(Xi[0][0]+Xi[1][1]);
		T e2=(Xi[0][0]-Xi[1][1]);
		T e6=(Xi[0][1]);

		//T e_i[3] = {e1, e2, e6};

		T e2_1=0.0, e2_2=0.0; 
		for (unsigned int i=0; i<DIM; ++i){
		  e2_1+=(chi[i][0]*dchi[i][0][0]-chi[i][1]*dchi[i][1][0]);
		  e2_2+=(chi[i][0]*dchi[i][0][1]-chi[i][1]*dchi[i][1][1]);
		}

		for (unsigned int i=0; i<DIM; ++i){
		  for (unsigned int J=0; J<DIM; ++J){
		    T e1_chiiJ=(chi[i][0]*(0==J)+chi[i][1]*(1==J));
		    T e2_chiiJ=(chi[i][0]*(0==J)-chi[i][1]*(1==J));
		    T e6_chiiJ=(chi[i][1]*(0==J)+chi[i][0]*(1==J))/2.0;
		    T e2_1_chiiJ=((0==J)*dchi[i][0][0]-(1==J)*dchi[i][1][0]);
		    T e2_2_chiiJ=((0==J)*dchi[i][0][1]-(1==J)*dchi[i][1][1]);

		    //P
		    P[i][J]=PiJ_2D;
		    P0[i][J]=P0iJ_2D;

		    //gradient terms
		    for (unsigned int K=0; K<DIM; ++K){
		T e2_1_chiiJK=(chi[i][0]*(0==J)-chi[i][1]*(1==J))*(0==K);
		T e2_2_chiiJK=(chi[i][0]*(0==J)-chi[i][1]*(1==J))*(1==K);

		//Beta
		Beta0[i][J][K]=Beta0iJK_2D;
		    }
		  }
		} 
	} 
  //3D model
	else if(DIM==3){
		T e1=(Xi[0][0]+Xi[1][1]+Xi[2][2])/sqrt(3.0);
		T e2=(Xi[0][0]-Xi[1][1])/sqrt(2.0);
		T e3=(Xi[0][0]+Xi[1][1]-2*Xi[2][2])/sqrt(6.0);
		T e4=Xi[1][2], e5=E[2][0], e6=Xi[0][1];

		T e2_1=0.0, e2_2=0.0, e2_3=0.0, e3_1=0.0, e3_2=0.0, e3_3=0.0;
		for (unsigned int i=0; i<DIM; ++i){
		  e2_1+=(chi[i][0]*dchi[i][0][0]-chi[i][1]*dchi[i][1][0])/sqrt(2.0);
		  e2_2+=(chi[i][0]*dchi[i][0][1]-chi[i][1]*dchi[i][1][1])/sqrt(2.0);
		  e2_3+=(chi[i][0]*dchi[i][0][2]-chi[i][1]*dchi[i][1][2])/sqrt(2.0);
		  e3_1+=(chi[i][0]*dchi[i][0][0]+chi[i][1]*dchi[i][1][0]-2*chi[i][2]*dchi[i][2][0])/sqrt(6.0);
		  e3_2+=(chi[i][0]*dchi[i][0][1]+chi[i][1]*dchi[i][1][1]-2*chi[i][2]*dchi[i][2][1])/sqrt(6.0);
		  e3_3+=(chi[i][0]*dchi[i][0][2]+chi[i][1]*dchi[i][1][2]-2*chi[i][2]*dchi[i][2][2])/sqrt(6.0);
		}

		for (unsigned int i=0; i<DIM; ++i){
		  for (unsigned int J=0; J<DIM; ++J){
		    T e1_chiiJ=(chi[i][0]*(0==J)+chi[i][1]*(1==J)+chi[i][2]*(2==J))/sqrt(3.0);
		    T e2_chiiJ=(chi[i][0]*(0==J)-chi[i][1]*(1==J))/sqrt(2.0);
		    T e3_chiiJ=(chi[i][0]*(0==J)+chi[i][1]*(1==J)-2*chi[i][2]*(2==J))/sqrt(6.0);
		    T e4_chiiJ=(chi[i][2]*(1==J)+chi[i][1]*(2==J))/2.0;
		    T e5_chiiJ=(chi[i][0]*(2==J)+chi[i][2]*(0==J))/2.0;
		    T e6_chiiJ=(chi[i][1]*(0==J)+chi[i][0]*(1==J))/2.0;
		    T e2_1_chiiJ=((0==J)*dchi[i][0][0]-(1==J)*dchi[i][1][0])/sqrt(2.0);
		    T e2_2_chiiJ=((0==J)*dchi[i][0][1]-(1==J)*dchi[i][1][1])/sqrt(2.0);
		    T e2_3_chiiJ=((0==J)*dchi[i][0][2]-(1==J)*dchi[i][1][2])/sqrt(2.0);     
		    T e3_1_chiiJ=((0==J)*dchi[i][0][0]+(1==J)*dchi[i][1][0]-2*(2==J)*dchi[i][2][0])/sqrt(6.0);
		    T e3_2_chiiJ=((0==J)*dchi[i][0][1]+(1==J)*dchi[i][1][1]-2*(2==J)*dchi[i][2][1])/sqrt(6.0);
		    T e3_3_chiiJ=((0==J)*dchi[i][0][2]+(1==J)*dchi[i][1][2]-2*(2==J)*dchi[i][2][2])/sqrt(6.0);

		    //P
		    P0[i][J]=P0iJ_3D;  
		    P[i][J]=PiJ_3D; 
		    //P[i][J]=PiJ(mu,alpha,beta,F,E,i,J); 
		    
		    //gradient terms
		    for (unsigned int K=0; K<DIM; ++K){
		T e2_1_chiiJK=(chi[i][0]*(0==J)-chi[i][1]*(1==J))*(0==K)/sqrt(2.0);
		T e2_2_chiiJK=(chi[i][0]*(0==J)-chi[i][1]*(1==J))*(1==K)/sqrt(2.0);
		T e2_3_chiiJK=(chi[i][0]*(0==J)-chi[i][1]*(1==J))*(2==K)/sqrt(2.0);
		T e3_1_chiiJK=(chi[i][0]*(0==J)+chi[i][1]*(1==J)-2*chi[i][2]*(2==J))*(0==K)/sqrt(6.0);
		T e3_2_chiiJK=(chi[i][0]*(0==J)+chi[i][1]*(1==J)-2*chi[i][2]*(2==J))*(1==K)/sqrt(6.0);
		T e3_3_chiiJK=(chi[i][0]*(0==J)+chi[i][1]*(1==J)-2*chi[i][2]*(2==J))*(2==K)/sqrt(6.0);

		//T ei_j_chiiJK[2][3] = {{e2_1_chiiJK, e2_2_chiiJK, e2_3_chiiJK}, {e3_1_chiiJK, e3_2_chiiJK, e3_3_chiiJK}};
		//Beta
		Beta0[i][J][K]=Beta0iJK_3D;
		    }
		  }
		}
	}
	else{
		PetscPrintf(PETSC_COMM_WORLD,"only material models for DIM=2, DIM=3 implemented.... but DIM input is %u\n",DIM); 
		exit(-1);      	
	}

  //Compute the Eshelby stress = \psi_N*\delta_{IJ} - F^T*P
  T Eshelby[DIM][DIM];
  for (unsigned int I=0; I<DIM; I++){
    for (unsigned int J=0; J<DIM; J++){
      Eshelby[I][J] = psi_N*(I==J);
      for (unsigned int k=0; k<DIM; k++){
	Eshelby[I][J] -= F[k][I]*P[k][J];
      }
    }
  }
  

  //get shape function values
  double (*N) = (double (*)) p->shape[0];
  double (*Nx)[DIM] = (double (*)[DIM]) p->shape[1];
  double (*Nxx)[DIM][DIM] = (double (*)[DIM][DIM]) p->shape[2];

  //compute residuals
  bool surfaceFlag=p->atboundary;
  for (unsigned int a=0; a<(unsigned int)nen; a++) {
    double N1[DIM], N2[DIM][DIM];
    for (unsigned int i=0; i<DIM; i++){
      N1[i]=Nx[a][i];
      for (unsigned int j=0; j<DIM; j++){
	N2[i][j]=Nxx[a][i][j];
      }
    }

    if (!surfaceFlag) {
      for (unsigned int i=0; i<DIM; i++){
	T Ru_I=0.0, Ru_i=0.0;
	for (unsigned int J=0; J<DIM; J++){
	  //grad(Na)*P
	  Ru_I += (P0[i][J] + J_chi*dpsi_dchi[i][J])*N1[J];
	  //Ru_I += P0[i][J]*N1[J]; //to replicate
	  for (unsigned int K=0; K<DIM; K++){
	    Ru_I += N2[J][K]*Beta0[i][J][K];
	    Ru_I += J_chi*(P[i][K] + Eshelby[i][K])*chi_Inv[J][K]*N1[J];
	    Ru_i += J_chi*P[i][K]*chi_Inv[J][K]*N1[J];
	  }
	}
	R[a*dof+i] = Ru_I;
	R[a*dof+i+DIM] = Ru_i;
      }
    }
  }

  return 0;
}

template PetscErrorCode Function<PetscReal,2,3>(IGAPoint p,PetscReal dt2,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const PetscReal * U,
			PetscReal t0,const PetscScalar * U0,
			PetscReal *R,void *ctx);
template PetscErrorCode Function<PetscReal,2,4>(IGAPoint p,PetscReal dt2,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const PetscReal * U,
			PetscReal t0,const PetscScalar * U0,
			PetscReal *R,void *ctx);
template PetscErrorCode Function<PetscReal,3,4>(IGAPoint p,PetscReal dt2,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const PetscReal * U,
			PetscReal t0,const PetscScalar * U0,
			PetscReal *R,void *ctx);
template PetscErrorCode Function<PetscReal,3,6>(IGAPoint p,PetscReal dt2,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const PetscReal * U,
			PetscReal t0,const PetscScalar * U0,
			PetscReal *R,void *ctx);

template PetscErrorCode Function<Sacado::Fad::SFad<double,27>,2,3>(IGAPoint p,PetscReal dt2,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const Sacado::Fad::SFad<double,27> * U,
			PetscReal t0,const PetscScalar * U0,
			Sacado::Fad::SFad<double,27> *R,void *ctx);
template PetscErrorCode Function<Sacado::Fad::SFad<double,36>,2,4>(IGAPoint p,PetscReal dt2,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const Sacado::Fad::SFad<double,36> * U,
			PetscReal t0,const PetscScalar * U0,
			Sacado::Fad::SFad<double,36> *R,void *ctx);
template PetscErrorCode Function<Sacado::Fad::SFad<double,108>,3,4>(IGAPoint p,PetscReal dt2,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const Sacado::Fad::SFad<double,108> * U,
			PetscReal t0,const PetscScalar * U0,
			Sacado::Fad::SFad<double,108> *R,void *ctx);
template PetscErrorCode Function<Sacado::Fad::SFad<double,162>,3,6>(IGAPoint p,PetscReal dt2,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const Sacado::Fad::SFad<double,162> * U,
			PetscReal t0,const PetscScalar * U0,
			Sacado::Fad::SFad<double,162> *R,void *ctx);
