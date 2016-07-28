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

  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  PetscReal *n = p->normal;

  //Retrieve parameters
  PetscReal Es = user->matParam["Es"];
  PetscReal Ed = user->matParam["Ed"];
  PetscReal El = user->matParam["El"];
  PetscReal Gl = user->matParam["Gl"];
  PetscReal Cs = user->matParam["Cs"];
  PetscReal Cd = user->matParam["Cd"];
  PetscReal Cl = user->matParam["Cl"];

  PetscReal DVal = user->matParam["DVal"]; //Diffusivity
  PetscReal CVal = user->matParam["CVal"];
  PetscReal gamma = user->matParam["gamma"];

  PetscReal flux[2*DIM];
  flux[0] = user->matParam["flux_xmin"];
  flux[1] = user->matParam["flux_xmax"];
  flux[2] = user->matParam["flux_ymin"];
  flux[3] = user->matParam["flux_ymax"];
  if(DIM==3){
    flux[4] = user->matParam["flux_zmin"];
    flux[5] = user->matParam["flux_zmax"];
  }

  //concentration field variable
  T c, cx[DIM], cxx[DIM][DIM]; PetscReal c0;
  computeField<T,DIM,DOF>(SCALAR,DIM,p,U,&c,&cx[0],&cxx[0][0]);
  computeField<PetscReal,DIM,DOF>(SCALAR,DIM,p,U0,&c0);

  //displacement field variable
  T u[DIM], ux[DIM][DIM], uxx[DIM][DIM][DIM];
  computeField<T,DIM,DOF>(VECTOR,0,p,&U[0],&u[0],&ux[0][0],&uxx[0][0][0]);

  //Compute F (I+Ux), dF (Uxx)
  T F[DIM][DIM], dF[DIM][DIM][DIM];
  for (unsigned int i=0; i<DIM; i++) {
    for (unsigned int J=0; J<DIM; J++) {
      F[i][J]=(i==J)+ux[i][J];
      for (unsigned int K=0; K<DIM; K++) {
	dF[i][J][K]=uxx[i][J][K];
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
  T P[DIM][DIM], Beta[DIM][DIM][DIM];
  T mu_c, mu_eX[DIM-1], eX_x[DIM-1][DIM];
 
  //2D model
  if(DIM==2){  
    //new strain metrics (2D)
    T e1=(E[0][0]+E[1][1]);
    T e2=(E[0][0]-E[1][1]);
    T e6=(E[0][1]);
    T e2_1=0.0, e2_2=0.0; 
    for (unsigned int i=0; i<DIM; ++i){
      e2_1+=(F[i][0]*dF[i][0][0]-F[i][1]*dF[i][1][0]);
      e2_2+=(F[i][0]*dF[i][0][1]-F[i][1]*dF[i][1][1]);
    }
    T e2x[DIM]; 
    eX_x[0][0] = e2_1; eX_x[0][1] = e2_2;
    //
    for (unsigned int i=0; i<DIM; ++i){
      for (unsigned int J=0; J<DIM; ++J){
	T e1_FiJ=(F[i][0]*(0==J)+F[i][1]*(1==J));
	T e2_FiJ=(F[i][0]*(0==J)-F[i][1]*(1==J));
	T e6_FiJ=(F[i][1]*(0==J)+F[i][0]*(1==J))/2.0;
	T e2_1_FiJ=((0==J)*dF[i][0][0]-(1==J)*dF[i][1][0]);
	T e2_2_FiJ=((0==J)*dF[i][0][1]-(1==J)*dF[i][1][1]);
	//P
	P[i][J]=PiJ_2D;

	//gradient terms
	for (unsigned int K=0; K<DIM; ++K){
	  T e2_1_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(0==K);
	  T e2_2_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(1==K);
	  //Beta
	  Beta[i][J][K]=BetaiJK_2D;
	}
      }
    }
    mu_c = mu_c_2D;
    mu_eX[0] = mu_e2_2D;
  }
  
  //3D model
  else if(DIM==3){
    T e1=(E[0][0]+E[1][1]+E[2][2])/sqrt(3.0);
    T e2=(E[0][0]-E[1][1])/sqrt(2.0);
    T e3=(E[0][0]+E[1][1]-2*E[2][2])/sqrt(6.0);
    T e4=E[1][2], e5=E[2][0], e6=E[0][1];
    T e2_1=0.0, e2_2=0.0, e2_3=0.0, e3_1=0.0, e3_2=0.0, e3_3=0.0;
    for (unsigned int i=0; i<DIM; ++i){
      e2_1+=(F[i][0]*dF[i][0][0]-F[i][1]*dF[i][1][0])/sqrt(2.0);
      e2_2+=(F[i][0]*dF[i][0][1]-F[i][1]*dF[i][1][1])/sqrt(2.0);
      e2_3+=(F[i][0]*dF[i][0][2]-F[i][1]*dF[i][1][2])/sqrt(2.0);
      e3_1+=(F[i][0]*dF[i][0][0]+F[i][1]*dF[i][1][0]-2*F[i][2]*dF[i][2][0])/sqrt(6.0);
      e3_2+=(F[i][0]*dF[i][0][1]+F[i][1]*dF[i][1][1]-2*F[i][2]*dF[i][2][1])/sqrt(6.0);
      e3_3+=(F[i][0]*dF[i][0][2]+F[i][1]*dF[i][1][2]-2*F[i][2]*dF[i][2][2])/sqrt(6.0);
    }
    T e2x[DIM], e3x[DIM]; 
    eX_x[0][0] = e2_1; eX_x[0][1] = e2_2; eX_x[0][2] = e2_3;
    eX_x[1][0] = e3_1; eX_x[1][1] = e3_2; eX_x[1][2] = e3_3;
    //
    for (unsigned int i=0; i<DIM; ++i){
      for (unsigned int J=0; J<DIM; ++J){
	T e1_FiJ=(F[i][0]*(0==J)+F[i][1]*(1==J)+F[i][2]*(2==J))/sqrt(3.0);
	T e2_FiJ=(F[i][0]*(0==J)-F[i][1]*(1==J))/sqrt(2.0);
	T e3_FiJ=(F[i][0]*(0==J)+F[i][1]*(1==J)-2*F[i][2]*(2==J))/sqrt(6.0);
	T e4_FiJ=(F[i][2]*(1==J)+F[i][1]*(2==J))/2.0;
	T e5_FiJ=(F[i][0]*(2==J)+F[i][2]*(0==J))/2.0;
	T e6_FiJ=(F[i][1]*(0==J)+F[i][0]*(1==J))/2.0;
	T e2_1_FiJ=((0==J)*dF[i][0][0]-(1==J)*dF[i][1][0])/sqrt(2.0);
	T e2_2_FiJ=((0==J)*dF[i][0][1]-(1==J)*dF[i][1][1])/sqrt(2.0);
	T e2_3_FiJ=((0==J)*dF[i][0][2]-(1==J)*dF[i][1][2])/sqrt(2.0);     
	T e3_1_FiJ=((0==J)*dF[i][0][0]+(1==J)*dF[i][1][0]-2*(2==J)*dF[i][2][0])/sqrt(6.0);
	T e3_2_FiJ=((0==J)*dF[i][0][1]+(1==J)*dF[i][1][1]-2*(2==J)*dF[i][2][1])/sqrt(6.0);
	T e3_3_FiJ=((0==J)*dF[i][0][2]+(1==J)*dF[i][1][2]-2*(2==J)*dF[i][2][2])/sqrt(6.0);
	//P
	P[i][J]=PiJ_3D;   
		    
	//gradient terms
	for (unsigned int K=0; K<DIM; ++K){
	  T e2_1_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(0==K)/sqrt(2.0);
	  T e2_2_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(1==K)/sqrt(2.0);
	  T e2_3_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(2==K)/sqrt(2.0);
	  T e3_1_FiJK=(F[i][0]*(0==J)+F[i][1]*(1==J)-2*F[i][2]*(2==J))*(0==K)/sqrt(6.0);
	  T e3_2_FiJK=(F[i][0]*(0==J)+F[i][1]*(1==J)-2*F[i][2]*(2==J))*(1==K)/sqrt(6.0);
	  T e3_3_FiJK=(F[i][0]*(0==J)+F[i][1]*(1==J)-2*F[i][2]*(2==J))*(2==K)/sqrt(6.0);
	  //Beta
	  Beta[i][J][K]=BetaiJK_3D;
	}
      }
    }
    mu_c = mu_c_3D;
    mu_eX[0] = mu_e2_3D;
    mu_eX[1] = mu_e3_3D;
  }
  else{
    PetscPrintf(PETSC_COMM_WORLD,"only material models for DIM=2, DIM=3 implemented.... but DIM input is %u\n",DIM); 
    exit(-1);    
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

    //mechanics
    if (!surfaceFlag) {
      for (unsigned int i=0; i<DIM; i++){
	T Ru_i=0.0;
	for (unsigned int J=0; J<DIM; J++){
	  //grad(Na)*P
	  Ru_i += N1[J]*P[i][J];
	  for (unsigned int K=0; K<DIM; K++){
	    Ru_i += N2[J][K]*Beta[i][J][K];
	  }
	}
	R[a*dof+i] = Ru_i;
      }
    }
  
    //chemistry
    T laplace_c=0;
    for (unsigned int i=0; i<DIM; i++) laplace_c+=cxx[i][i];
    T Rc=0.0;
    if (!surfaceFlag){
      // Na * c_t
      Rc += N[a] * (c-c0)*(1.0/dt2);
      // grad(Na) . D*(dmuc*grad(C)+dmue2*grad(e2)+dmue3*grad(e3))
      double laplace_N=0.0; T Nxcx=0.0, Nxe2x=0.0;
      for (unsigned int i=0; i<DIM; i++){
	laplace_N+=N2[i][i];
	Nxcx+=N1[i]*cx[i];
	Nxe2x+=N1[i]*eX_x[0][i];
      }
      Rc += DVal*mu_c*Nxcx + DVal*mu_eX[0]*Nxe2x;
      if(DIM==3){ 
	T Nxe3x=0.0;
	for (unsigned int i=0; i<DIM; i++){
	  Nxe3x+=N1[i]*eX_x[1][i];      
	}
	Rc += DVal*mu_eX[1]*Nxe3x;
      }
      // lambda * del2(Na) * D * del2(c)
      Rc += Cl*laplace_N*DVal*laplace_c;
    }
    else{
      // -grad(Na) . (D*del2(c)) n
      T t1 = DVal*laplace_c;
      double laplace_N=0.0;
      for (unsigned int i=0; i<DIM; i++){
	Rc += -N1[i]*t1*n[i];
	laplace_N += N2[i][i];
      }
      // -(gamma*del2(Na)*D)*grad(C).n
      T t2 = gamma*laplace_N*DVal;
      for (unsigned int i=0; i<DIM; i++){
	Rc += -t2*cx[i]*n[i];
      }
      // (C/he)*(grad(Na).n)*D*(grad(C).n)
      double t3=0.0;
      T t4 = (CVal/user->he)*DVal;
      T t5=0.0;
      for (unsigned int i=0; i<DIM; i++){
	t3 += N1[i]*n[i];
	t5 += cx[i]*n[i];
      }
      Rc += t3*t4*t5;
      Rc *= Cl;
      //flux term, Na*J
      for(unsigned int i=0; i<DIM; ++i){
	if(p->boundary_id == i){
	  Rc += -N[a]*flux[i];
	}
      }
    }
    R[a*dof+DIM] = Rc;  
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
