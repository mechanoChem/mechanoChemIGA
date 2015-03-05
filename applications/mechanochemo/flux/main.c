#include <math.h> 
extern "C" {
#include "petiga.h"
}
#include "../../../include/fields.h"
#include "../../../include/petigaksp2.h"
//
#define DIM 2
#define ADSacado
#define numVars 27
typedef struct {
  IGA iga;
  PetscReal Cd, Cg;
  PetscReal Es, Ed, E4, E6, E3, E2, Eii, Eij, Eg, curvature;
  PetscReal dt;
  PetscReal gamma, C, he, cbar, D, flux;
  AppCtxKSP* appCtxKSP;
  PetscReal f0Norm;
} AppCtx;
#include "../../../include/output.h"
#include "../../../include/evaluators.h"
//
#define EXPLICIT
#define finiteStrain
#define EgVAL 1.0
#define dtVAL 1.0e-7
#define NVAL 100
#define FLUX 3
#define bcVAL 2
#define lVAL 0.01
//
//#define PI 3.14159265


#undef  __FUNCT__
#define __FUNCT__ "Function"
template <class T>
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

  //concentration field variable
  T c, cx[DIM], cxx[DIM][DIM]; PetscReal c0;
  computeField<T,DIM,DIM+1>(SCALAR,DIM,p,U,&c,&cx[0],&cxx[0][0]);
  computeField<PetscReal,DIM,DIM+1>(SCALAR,DIM,p,U0,&c0);
  
  //displacement field variable
  T u[DIM], ux[DIM][DIM], uxx[DIM][DIM][DIM];
  computeField<T,DIM,DIM+1>(VECTOR,0,p,&U[0],&u[0],&ux[0][0],&uxx[0][0][0]);

  //problem parameters
  PetscReal Cd=user->Cd;
  PetscReal C4=16*Cd;
  PetscReal C3=-32*Cd;  
  PetscReal C2=16*Cd;
  PetscReal Cg=user->Cg;
  PetscReal Es=user->Es;
  PetscReal Ed=user->Ed;
  PetscReal E6=user->E6;
  PetscReal E4=user->E4;
  PetscReal E3=user->E3;
  PetscReal E2=user->E2;
  PetscReal Eii=user->Eii;
  PetscReal Eij=user->Eij;
  PetscReal Eg=user->Eg;
  PetscReal C=user->C;
  PetscReal he=user->he;
  PetscReal D=user->D;
  PetscReal flux=user->flux;
  PetscReal gamma=user->gamma;
  PetscReal dt=dt2;
  //
  PetscReal Ev=Es*0.1;

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

  //Compute strain metric, E  
  // E=0.5*(F^T*F-I) (finite strain)
  // E=0.5*(Ux+Ux^T)  (small strain)
  T E[DIM][DIM];
  for (unsigned int I=0; I<DIM; I++){
    for (unsigned int J=0; J<DIM; J++){
#ifdef finiteStrain
      E[I][J] = -0.5*(I==J);
      for (unsigned int k=0; k<DIM; k++){
	E[I][J] += 0.5*F[k][I]*F[k][J];
      }
#else
      E[I][J] = 0.5*(F[I][J]+F[J][I]-2*(I==J));
#endif 
   }
  }

  //new strain metrics
#if DIM==3
  T e1=(E[0][0]+E[1][1]+E[2][2])/sqrt(3.0);
  T e2=(E[0][0]-E[1][1])/sqrt(2.0);
  T e3=(E[0][0]+E[1][1]-2*E[2][2])/sqrt(6.0);
  T e4=E[1][2], e5=E[2][0], e6=E[0][1];
  T e2_1=0.0, e2_2=0.0, e2_3=0.0, e3_1=0.0, e3_2=0.0, e3_3=0.0;
#ifdef finiteStrain
  for (unsigned int i=0; i<DIM; ++i){
    e2_1+=(F[i][0]*dF[i][0][0]-F[i][1]*dF[i][1][0])/sqrt(2.0);
    e2_2+=(F[i][0]*dF[i][0][1]-F[i][1]*dF[i][1][1])/sqrt(2.0);
    e2_3+=(F[i][0]*dF[i][0][2]-F[i][1]*dF[i][1][2])/sqrt(2.0);
    e3_1+=(F[i][0]*dF[i][0][0]+F[i][1]*dF[i][1][0]-2*F[i][2]*dF[i][2][0])/sqrt(6.0);
    e3_2+=(F[i][0]*dF[i][0][1]+F[i][1]*dF[i][1][1]-2*F[i][2]*dF[i][2][1])/sqrt(6.0);
    e3_3+=(F[i][0]*dF[i][0][2]+F[i][1]*dF[i][1][2]-2*F[i][2]*dF[i][2][2])/sqrt(6.0);
  }
#else
  e2_1=(dF[0][0][0]-dF[1][1][0])/sqrt(2.0);
  e2_2=(dF[0][0][1]-dF[1][1][1])/sqrt(2.0);
  e2_3=(dF[0][0][2]-dF[1][1][2])/sqrt(2.0);
  e3_1=(dF[0][0][0]+dF[1][1][0]-2*dF[2][2][0])/sqrt(6.0);
  e3_2=(dF[0][0][1]+dF[1][1][1]-2*dF[2][2][1])/sqrt(6.0);
  e3_3=(dF[0][0][2]+dF[1][1][2]-2*dF[2][2][2])/sqrt(6.0);
#endif

  //Pi=E4(e2^2+e3^2)^2+E3*e3*(e3^2-3*e2^2)+E2*(e2^2+e3^2)+Eii*(e1^2)+Eij*(e4^2+e5^2+e6^2)
  //   +Eg(e2_1^2+e2_2^2+e2_3^2+e3_1^2+e3_2^2+e3^2)
  //compute P and Beta
  T P[DIM][DIM], Beta[DIM][DIM][DIM];
#ifdef EXPLICIT
  PetscReal E2c=E2*(5*c0-2.0)/3.0,  E3c=E3*c0, E4c=E4*c0;
#else
  T E2c=E2*(5*c-2.0)/3.0,  E3c=E3*c, E4c=E4*c;
#endif
  //
  for (unsigned int i=0; i<DIM; ++i){
    for (unsigned int J=0; J<DIM; ++J){
#ifdef finiteStrain
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
#else
      //_FiJ here imples _UiJ
      T e1_FiJ=((0==i)*(0==J) + (1==i)*(1==J) + (2==i)*(2==J))/sqrt(3.0);
      T e2_FiJ=((0==i)*(0==J) - (1==i)*(1==J))/sqrt(2.0);
      T e3_FiJ=((0==i)*(0==J) + (1==i)*(1==J) - 2*(2==i)*(2==J))/sqrt(6.0);
      T e4_FiJ=((1==i)*(2==J) + (2==i)*(1==J))/2.0;
      T e5_FiJ=((2==i)*(0==J) + (0==i)*(2==J))/2.0;
      T e6_FiJ=((0==i)*(1==J) + (1==i)*(0==J))/2.0;
      T e2_1_FiJ=0.0; 
      T e2_2_FiJ=0.0; 
      T e2_3_FiJ=0.0; 
      T e3_1_FiJ=0.0; 
      T e3_2_FiJ=0.0; 
      T e3_3_FiJ=0.0; 
#endif
      //P
      P[i][J]=(2*Eii*e1)*e1_FiJ						\
	+(2*E2c*e2 -6*E3c*e2*e3 + 4*E4c*e2*(e2*e2+e3*e3))*e2_FiJ	\
	+(2*E2c*e3 +3*E3c*(e3*e3-e2*e2) + 4*E4c*e3*(e2*e2+e3*e3))*e3_FiJ \
	+(2*Eij*e4)*e4_FiJ						\
	+(2*Eij*e5)*e5_FiJ						\
	+(2*Eij*e6)*e6_FiJ						\
	+ 2*Eg*(e2_1*e2_1_FiJ + e2_2*e2_2_FiJ + e2_3*e2_3_FiJ + e3_1*e3_1_FiJ + e3_2*e3_2_FiJ + e3_3*e3_3_FiJ);
      
      //gradient terms
      for (unsigned int K=0; K<DIM; ++K){
#ifdef finiteStrain
	T e2_1_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(0==K)/sqrt(2.0);
	T e2_2_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(1==K)/sqrt(2.0);
	T e2_3_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(2==K)/sqrt(2.0);
	T e3_1_FiJK=(F[i][0]*(0==J)+F[i][1]*(1==J)-2*F[i][2]*(2==J))*(0==K)/sqrt(6.0);
	T e3_2_FiJK=(F[i][0]*(0==J)+F[i][1]*(1==J)-2*F[i][2]*(2==J))*(1==K)/sqrt(6.0);
	T e3_3_FiJK=(F[i][0]*(0==J)+F[i][1]*(1==J)-2*F[i][2]*(2==J))*(2==K)/sqrt(6.0);
#else
	T e2_1_FiJK=((0==i)*(0==J)-(1==i)*(1==J))*(0==K)/sqrt(2.0);
	T e2_2_FiJK=((0==i)*(0==J)-(1==i)*(1==J))*(1==K)/sqrt(2.0);
	T e2_3_FiJK=((0==i)*(0==J)-(1==i)*(1==J))*(2==K)/sqrt(2.0);
	T e3_1_FiJK=((0==i)*(0==J)+(1==i)*(1==J) -2*(2==i)*(2==J))*(0==K)/sqrt(6.0);
	T e3_2_FiJK=((0==i)*(0==J)+(1==i)*(1==J) -2*(2==i)*(2==J))*(1==K)/sqrt(6.0);
	T e3_3_FiJK=((0==i)*(0==J)+(1==i)*(1==J) -2*(2==i)*(2==J))*(2==K)/sqrt(6.0);
#endif
	//Beta
	Beta[i][J][K]= 	2*Eg*(e2_1*e2_1_FiJK + e2_2*e2_2_FiJK + e2_3*e2_3_FiJK + e3_1*e3_1_FiJK + e3_2*e3_2_FiJK + e3_3*e3_3_FiJK);
      }
    }
  }

#elif DIM==2
  T e1=(E[0][0]+E[1][1]);
  T e2=(E[0][0]-E[1][1]);
  T e3=E[0][1];
  T e2_1=0.0, e2_2=0.0; 
#ifdef finiteStrain
  for (unsigned int i=0; i<DIM; ++i){
    e2_1+=(F[i][0]*dF[i][0][0]-F[i][1]*dF[i][1][0]);
    e2_2+=(F[i][0]*dF[i][0][1]-F[i][1]*dF[i][1][1]);
  }
#else
  e2_1=(dF[0][0][0]-dF[1][1][0]);
  e2_2=(dF[0][0][1]-dF[1][1][1]);
#endif
  //Pi=E4*e2^4+E2*e2^2++Eii*(e1^2)+Eij*(e3^2)
  //   +Eg(e2_1^2+e2_2^2)
  //compute P and Beta
  T P[DIM][DIM], Beta[DIM][DIM][DIM];
#ifdef EXPLICIT
  PetscReal E2c=E2*(2*c0-1.0),  E4c=E4, E6c=E6;  
#else
  T E2c=E2*(2*c-1.0),  E4c=E4, E6c=E6;
#endif
  //
  for (unsigned int i=0; i<DIM; ++i){
    for (unsigned int J=0; J<DIM; ++J){
#ifdef finiteStrain
      T e1_FiJ=(F[i][0]*(0==J)+F[i][1]*(1==J));
      T e2_FiJ=(F[i][0]*(0==J)-F[i][1]*(1==J));
      T e3_FiJ=(F[i][1]*(0==J)+F[i][0]*(1==J))/2.0;
      T e2_1_FiJ=((0==J)*dF[i][0][0]-(1==J)*dF[i][1][0]);
      T e2_2_FiJ=((0==J)*dF[i][0][1]-(1==J)*dF[i][1][1]);
#else
      T e1_FiJ=((0==i)*(0==J) + (1==i)*(1==J));
      T e2_FiJ=((0==i)*(0==J) - (1==i)*(1==J));
      T e3_FiJ=((0==i)*(1==J) + (1==i)*(0==J))/2.0;
      T e2_1_FiJ=0.0; 
      T e2_2_FiJ=0.0; 
#endif
      //P
      P[i][J]=(2*Eii*e1)*e1_FiJ						\
	+(2*E2c*e2 + 4*E4c*e2*e2*e2+ 6*E6c*e2*e2*e2*e2*e2)*e2_FiJ				\
	+(2*Eij*e3)*e3_FiJ						\
	+ 2*Eg*(e2_1*e2_1_FiJ + e2_2*e2_2_FiJ);
      
      //gradient terms
      for (unsigned int K=0; K<DIM; ++K){
#ifdef finiteStrain
	T e2_1_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(0==K);
	T e2_2_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(1==K);
#else
	T e2_1_FiJK=((0==i)*(0==J)-(1==i)*(1==J))*(0==K);
	T e2_2_FiJK=((0==i)*(0==J)-(1==i)*(1==J))*(1==K);
#endif
	//Beta
	Beta[i][J][K]= 	2*Eg*(e2_1*e2_1_FiJK + e2_2*e2_2_FiJK);
      }
    }
  }
#endif 
  
  //phi(c,e)=C4*c^4 + C3*c^3 + C2*c^2 + Cg*(c_1^2+c_2^2) + strain dependent terms 
  //chemical potential
  T mu = 4.0*C4*c*c*c + 3.0*C3*c*c + 2.0*C2*c;
  T dmuc= 12.0*C4*c*c + 6.0*C3*c + 2.0*C2;
  T dmue2=0.0, dmue3=0.0;
#ifndef EXPLICIT
#if DIM==3
  mu+= E4*(e2*e2+e3*e3)*(e2*e2+e3*e3) + E3*e3*(e3*e3-3*e2*e2) + E2*(5.0/3.0)*(e2*e2+e3*e3);
  dmue2+= 4.0*E4*(e2*e2+e3*e3)*e2 - 6.0*E3*e3*e2 + 2.0*E2*(5.0/3.0)*e2;
  dmue3+= 4.0*E4*(e2*e2+e3*e3)*e3 + 3.0*E3*e3*e3 + 2.0*E2*(5.0/3.0)*e3;
#elif DIM==2
  mu+= 2.0*E2*(e2*e2);
  dmue2+= 4.0*E2*e2;
#endif
#endif


  /* //get shape function values */
  double (*N) = (double (*)) p->shape[0];
  double (*Nx)[DIM] = (double (*)[DIM]) p->shape[1];
  double (*Nxx)[DIM][DIM] = (double (*)[DIM][DIM]) p->shape[2];

  //Compute Residual
  bool surfaceFlag=p->atboundary;
  for (unsigned int a=0; a<(unsigned int)nen; a++) {
    double N1[DIM], N2[DIM][DIM];
    for (unsigned int i=0; i<DIM; i++){
      N1[i]=Nx[a][i];
      for (unsigned int j=0; j<DIM; j++){
	N2[i][j]=Nxx[a][i][j];
      }
    }

    //Mechanics
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

    //Chemistry
    T laplace_c=0;
    for (unsigned int i=0; i<DIM; i++) laplace_c+=cxx[i][i];
    T Rc=0.0;
    if (!surfaceFlag){
      // Na * c_t
      Rc += N[a] * (c-c0)*(1.0/dt);
      // grad(Na) . D*(dmuc*grad(C)+dmue2*grad(e2)+dmue3*grad(e3))
      double laplace_N=0.0;
      for (unsigned int i=0; i<DIM; i++){
#if DIM==3
	T e2x=0, e3x=0;
	switch (i) {
	case 0:
	   e2x=e2_1; e3x=e3_1; break;
	case 1:
	   e2x=e2_2; e3x=e3_2; break;
	case 2:
	   e2x=e2_3; e3x=e3_3; break;
	}
	Rc += N1[i]*D*(dmuc*cx[i]+dmue2*e2x+dmue3*e3x);
#elif DIM==2
	T e2x=0;
	switch (i) {
	case 0:
	   e2x=e2_1; break;
	case 1:
	   e2x=e2_2; break;
	}
	Rc += N1[i]*D*(dmuc*cx[i]+dmue2*e2x);
#endif
	laplace_N += N2[i][i];
      }
      // lambda * del2(Na) * D * del2(c)
      Rc += Cg*laplace_N*D*laplace_c;
    }
    else{
      // -grad(Na) . (D*del2(c)) n
      T t1 = D*laplace_c;
      double laplace_N=0.0;
      for (unsigned int i=0; i<DIM; i++){
	Rc += -N1[i]*t1*n[i];
	laplace_N += N2[i][i];
      }
      // -(gamma*del2(Na)*D)*grad(C).n
      T t2 = gamma*laplace_N*D;
      for (unsigned int i=0; i<DIM; i++){
	Rc += -t2*cx[i]*n[i];
      }
      // (C/he)*(grad(Na).n)*D*(grad(C).n)
      double t3=0.0;
      T t4 = (C/he)*D;
      T t5=0.0;
      for (unsigned int i=0; i<DIM; i++){
	t3 += N1[i]*n[i];
	t5 += cx[i]*n[i];
      }
      Rc += t3*t4*t5;
      Rc *= Cg;
      //flux term, Na*J
#if FLUX==3
      Rc +=0.0;
#elif FLUX==0
      Rc += -N[a]*flux;
#elif FLUX==1
      if (n[0]==0.0) Rc += -N[a]*flux;
#elif FLUX==2
      if (n[1]==0.0) Rc += -N[a]*flux;
#else
      PetscPrintf(PETSC_COMM_WORLD,"FLUX key illdefined"); 
      exit(-1);      
#endif
    }
    R[a*dof+DIM] = Rc;
  }
  return 0;
}

typedef struct {
  PetscReal ux, uy;
#if DIM==3  
  PetscReal uz;
#endif
  PetscReal c;
} Field;

PetscErrorCode FormInitialCondition(IGA iga, Vec U, AppCtx *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;
  std::srand(5);
  DM da;
  ierr = IGACreateNodeDM(iga,3,&da);CHKERRQ(ierr);
  Field **u;
  ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
  PetscInt i,j;
  for(i=info.xs;i<info.xs+info.xm;i++){
    for(j=info.ys;j<info.ys+info.ym;j++){
      u[j][i].ux=0.0;
      u[j][i].uy=0.0;
      u[j][i].c= user->cbar + 0.01*(0.5 - (double)(std::rand() % 100 )/100.0);
    }
  }
  ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
  ierr = DMDestroy(&da);;CHKERRQ(ierr); 
  PetscFunctionReturn(0); 
}

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);
  double startTime = MPI_Wtime();
 
  /* Define simulation specific parameters */
  AppCtx user; AppCtxKSP userKSP;

  //problem parameters
  user.Cd=1.0;
  user.Cg=100*pow(lVAL,2.0);
  user.D=10.0;
  user.flux=100.0;
#if FLUX==3 //For Quench
  user.cbar=0.6;
#else
  user.cbar=0.01;
#endif
  user.gamma=1.0;
  user.C=5.0;
  //
  user.Es=0.01; //0.2391; //0.1
  user.Ed=1.0e-3; //1.2847e-04; //1.0
#if DIM==3
  user.E4=1.5*user.Ed/pow(user.Es,2.0);
  user.E3=-user.Ed/user.Es;
  user.E2=-1.5*user.Ed;
#elif DIM==2
  user.curvature=15*user.Ed/pow(user.Es,2.0); //0.0341; //12*user.Ed/pow(user.Es,2.0);
  user.E6=-user.Ed/pow(user.Es,6.0)+user.curvature/(8.0*pow(user.Es,4.0));
  user.E4= 3.0*user.Ed/pow(user.Es,4.0)-user.curvature/(4.0*pow(user.Es,2.0));
  user.E3=0.0;
  user.E2=-3.0*user.Ed/pow(user.Es,2.0) + user.curvature/8.0; 
#endif
  user.Eii=-2*user.E2;
  user.Eij=-2*user.E2;
  user.Eg=EgVAL*(-2*user.E2)*pow(lVAL,2.0);
  //
  user.dt=dtVAL; 
  //
  PetscPrintf(PETSC_COMM_WORLD,"\n\nCg:%8.2e, D:%8.2e, flux:%8.2e, cbar:%8.2e\n",user.Cg,user.D,user.flux,user.cbar);
  PetscPrintf(PETSC_COMM_WORLD,"E6:%8.2e, E4:%8.2e, E2:%8.2e, Eii:%8.2e, Eij:%8.2e, Eg:%8.2e, d:%8.2e, s:%8.2e, curvature:%8.2e\n",user.E6,user.E4,user.E2,user.Eii,user.Eij,user.Eg,user.Ed,user.Es,user.curvature);
  //
  /* Set discretization options */
  PetscInt nsteps = 100;
  PetscInt N=NVAL, p=2, C=PETSC_DECIDE, resStep=0;
  PetscBool output = PETSC_TRUE; 
  PetscBool monitor = PETSC_TRUE; 
  char filePrefix[PETSC_MAX_PATH_LEN] = {0};
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","CahnHilliard2D Options","IGA");CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  if (C == PETSC_DECIDE) C = p-1;
 
  //
  if (p < 2 || C < 0) /* Problem requires a p>=2 C1 basis */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Problem requires minimum of p = 2");
  if (p <= C)         /* Check C < p */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Discretization inconsistent: polynomial order must be greater than degree of continuity");
  
  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDim(iga,DIM);CHKERRQ(ierr);
  ierr = IGASetDof(iga,DIM+1);CHKERRQ(ierr);

  IGAAxis axis0;
  ierr = IGAGetAxis(iga,0,&axis0);CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis0,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis0,N,0.0,1.0,C);CHKERRQ(ierr);
  user.he=1.0/N; 

  IGAAxis axis1;
  ierr = IGAGetAxis(iga,1,&axis1);CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis1,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis1,N,0.0,1.0,C);CHKERRQ(ierr);
#if DIM==3
  IGAAxis axis2;
  ierr = IGAGetAxis(iga,2,&axis2);CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis2,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis2,N,0.0,1.0,C);CHKERRQ(ierr);
#endif
  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  user.iga = iga;
  if (resStep==0){
    char meshfilename[256];
    sprintf(meshfilename, "mesh.dat");
    PetscPrintf(PETSC_COMM_WORLD,"\nWriting mesh file: %s\n", meshfilename);
    ierr = IGAWrite(iga, meshfilename);CHKERRQ(ierr);
  }
  
  Vec U,U0;
  ierr = IGACreateVec(iga,&U);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&U0);CHKERRQ(ierr);
  if (resStep>0){
    MPI_Comm comm;
    PetscViewer viewer;
    char restartfilename[256];
    ierr = PetscObjectGetComm((PetscObject)U0,&comm);CHKERRQ(ierr);
    sprintf(restartfilename,"res%s-%d.dat",filePrefix,resStep);
    ierr = PetscViewerBinaryOpen(comm,restartfilename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecLoad(U0,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);
    PetscPrintf(PETSC_COMM_WORLD,"\nReading solution from restart file: %s\n",restartfilename);
  }
  else{
    ierr = FormInitialCondition(iga, U0, &user); 
  }
  ierr = VecCopy(U0, U);CHKERRQ(ierr);
  
  //
  IGAForm form;
  ierr = IGAGetForm(iga,&form);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,1,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,1,PETSC_TRUE);CHKERRQ(ierr);
#if DIM==3
  ierr = IGAFormSetBoundaryForm (form,2,0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,2,1,PETSC_TRUE);CHKERRQ(ierr);
#endif

  //Dirichlet BC
  double dVal=user.Es*.01;
#if DIM==3
  ierr = IGASetBoundaryValue(iga,0,0,1,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,0,1,1,-dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,0,0,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,1,0,-dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,2,0,2,0.0);CHKERRQ(ierr);  
#elif DIM==2 
#if bcVAL==0
  //Shear BC
  ierr = IGASetBoundaryValue(iga,0,0,1,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,0,1,1,-dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,0,0,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,1,0,-dVal);CHKERRQ(ierr);
#elif bcVAL==1
  //Free BC
  ierr = IGASetBoundaryValue(iga,0,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,0,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,0,1,0,dVal);CHKERRQ(ierr);  
#elif bcVAL==2
  //Fixed BC
  ierr = IGASetBoundaryValue(iga,0,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,0,0,1,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,0,0,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,1,0,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,1,1,0,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,1,1,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,0,1,0,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,0,1,1,0.0);CHKERRQ(ierr);  
#endif
#endif
  //
  ierr = IGASetFormIEFunction(iga,Residual,&user);CHKERRQ(ierr);
  ierr = IGASetFormIEJacobian(iga,Jacobian,&user);CHKERRQ(ierr);
  user.appCtxKSP=&userKSP;
  userKSP.U0=&U;
  userKSP.Es=user.Es;
  //
  TS ts;
  ierr = IGACreateTS(iga,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,100000,1.0);CHKERRQ(ierr);
  ierr = TSSetTime(ts,0.0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,user.dt);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,OutputMonitor<DIM>,&user,NULL);CHKERRQ(ierr);  
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
#if PETSC_VERSION_LE(3,3,0)
  ierr = TSSolve(ts,U,NULL);CHKERRQ(ierr);
#else
  ierr = TSSolve(ts,U);CHKERRQ(ierr);
#endif
  //
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&U0);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = IGADestroy(&iga);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

