#include <math.h> 
extern "C" {
#include "petiga.h"
}
#include "fields.h"
#include "petigaksp2.h"

//set small strain or finite strain
//#define finiteStrain

//include automatic differentiation library
//#define ADSacado
#ifdef ADSacado
#include <Sacado.hpp>
#define numVars 81 //81 for 3D, 18 for 2D
typedef Sacado::Fad::SFad<double,numVars> doubleAD;
//typedef Sacado::Fad::DFad<double> doubleAD;
#else
#include "adept.h"
typedef adept::adouble doubleAD;
#endif

#define DIM 2
#define PI 3.14159265

typedef struct {
  IGA iga;
  PetscReal Es, Ed, E4, E3, E2, Eii, Eij, Eg, El;
  PetscReal c, dt;
  PetscReal C, he;
  AppCtxKSP* appCtxKSP;
  PetscReal f0Norm;
} AppCtx;

#undef  __FUNCT__
#define __FUNCT__ "Function"
template <class T>
PetscErrorCode Function(IGAPoint p,PetscReal dt2,
				PetscReal shift,const PetscScalar *V,
				PetscReal t,const T *U,
				PetscReal t0,const PetscScalar *U0,
				T *R,void *ctx)
{
  AppCtx *user = (AppCtx *)ctx;

  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  PetscReal *n = p->normal;

  //displacement field variables
  T u[DIM], ux[DIM][DIM], uxx[DIM][DIM][DIM];
  computeField<T,DIM,DIM>(VECTOR,0,p,U,&u[0],&ux[0][0],&uxx[0][0][0]);

  //problem parameters

  PetscReal Es=user->Es;
  PetscReal Ed=user->Ed;
  PetscReal E4=user->E4;
  PetscReal E3=user->E3;
  PetscReal E2=user->E2;
  PetscReal Eii=user->Eii;
  PetscReal Eij=user->Eij;
  PetscReal Eg=user->Eg;
  PetscReal El=user->El;
  PetscReal c=user->c+2*t0;
  PetscReal C=user->C;
  PetscReal he=user->he;
  c=1.0;

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
      E[I][J] = 0.5*(ux[I][J]+ux[J][I]);
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
  PetscReal E2c=E2, E3c=E3, E4c=E4;
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
  PetscReal E2c=E2, E4c=E4;
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
	+(2*E2c*e2 + 4*E4c*e2*e2*e2)*e2_FiJ				\
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
  
  /* //get shape function values */
  double (*N) = (double (*)) p->shape[0];
  double (*Nx)[DIM] = (double (*)[DIM]) p->shape[1];
  double (*Nxx)[DIM][DIM] = (double (*)[DIM][DIM]) p->shape[2];

  //Compute Residual
  bool surfaceFlag=p->atboundary;
  T (*Ra)[DIM] = (T (*)[DIM])R;
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
	for (unsigned int j=0; j<DIM; j++){
	  //grad(Na)*P
	  Ru_i += N1[j]*P[i][j];
	  for (unsigned int k=0; k<DIM; k++){
	    Ru_i += N2[j][k]*Beta[i][j][k];
	  }
	}
	Ra[a][i] = Ru_i;
      }
    }
  }
  return 0;
}


#undef  __FUNCT__
#define __FUNCT__ "Residual"
PetscErrorCode Residual(IGAPoint p,PetscReal dt,
                        PetscReal shift,const PetscScalar *V,
                        PetscReal t,const PetscScalar *U,
                        PetscReal t0,const PetscScalar *U0, 
			PetscScalar *R,void *ctx)
{
  Function(p, dt, shift, V, t, U, t0, U0, R, ctx);
  return 0;
}

#undef  __FUNCT__
#define __FUNCT__ "Jacobian"
PetscErrorCode Jacobian(IGAPoint p,PetscReal dt,
				PetscReal shift,const PetscScalar *V,
				PetscReal t,const PetscScalar *U,
				PetscReal t0,const PetscScalar *U0,
				PetscScalar *K,void *ctx)
{
  AppCtx *user = (AppCtx *)ctx;
  const PetscInt nen=p->nen, dof=DIM;
  const PetscReal (*U2)[DIM] = (PetscReal (*)[DIM])U;
#ifdef ADSacado
  if (dof*nen!=numVars) {
    PetscPrintf(PETSC_COMM_WORLD,"\ndof*nen!=numVars.... Set numVars = %u\n",dof*nen); exit(-1);
  }
 std::vector<doubleAD> U_AD(nen*DIM);
  for(int i=0; i<nen*dof; i++){
    U_AD[i]=U[i];
    U_AD[i].diff(i, dof*nen);
  } 
  std::vector<doubleAD> R(nen*dof);
  Function<doubleAD> (p, dt, shift, V, t, &U_AD[0], t0, U0, &R[0], ctx);
  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      for(int n2=0; n2<nen; n2++){
	for(int d2=0; d2<dof; d2++){
      	  K[n1*dof*nen*dof + d1*nen*dof + n2*dof + d2] = R[n1*dof+d1].dx(n2*dof+d2);
	}
      }
    }				
  }
#else
  adept::Stack s;
  std::vector<doubleAD> U_AD(nen*DIM);
  adept::set_values(&U_AD[0],nen*dof,U);
  s.new_recording();
  std::vector<doubleAD> R(nen*dof);
  Function<doubleAD> (p, dt, shift, V, t, &U_AD[0], t0, U0, &R[0], ctx);
  s.independent(&U_AD[0],nen*dof);
  s.dependent(&R[0],nen*dof);
  s.jacobian(K);
#endif  
  return 0;    
}

PetscErrorCode E22System(IGAPoint p, PetscScalar *K, PetscScalar *R, void *ctx)
{	
  AppCtxKSP *user = (AppCtxKSP *)ctx;
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  //displacement field variables
  PetscReal u[DIM], ux[DIM][DIM];
  computeField<PetscReal,DIM,DIM>(VECTOR,0,p,user->localU0,&u[0],&ux[0][0]);
  //Compute F
  PetscReal F[DIM][DIM];
  for (PetscInt i=0; i<DIM; i++) {
    for (PetscInt j=0; j<DIM; j++) {
      F[i][j]=(i==j)+ux[i][j];
    }
  }
  //Compute strain metric, E  (E=0.5*(F^T*F-I))
  PetscReal E[DIM][DIM];
  for (unsigned int I=0; I<DIM; I++){
    for (unsigned int J=0; J<DIM; J++){
#ifdef finiteStrain
      E[I][J] = -0.5*(I==J);
      for (unsigned int k=0; k<DIM; k++){
	E[I][J] += 0.5*F[k][I]*F[k][J];
      }
#else
      E[I][J] = 0.5*(ux[I][J]+ux[J][I]);
#endif 
   }
  }

#if DIM==3  
 //new strain metrics
  PetscReal e1=(E[0][0]+E[1][1]+E[2][2])/sqrt(3.0);
  PetscReal e2=(E[0][0]-E[1][1])/sqrt(2.0);
  PetscReal e3=(E[0][0]+E[1][1]-2*E[2][2])/sqrt(6.0);
  PetscReal e4=E[1][2], e5=E[2][0], e6=E[0][1];
  //PetscPrintf(PETSC_COMM_WORLD,"e1: %8.2e, e2: %8.2e\n",e2,e3);
  
  //compute distance to nearest well
  PetscReal Es=user->Es;
  PetscReal x[3],y[3]; 
  x[0]=0; y[0]=Es; //first well 
  x[1]=-Es*cos(30.0*PI/180.0); y[1]=-Es*sin(30.0*PI/180.0); //second well
  x[2]=+Es*cos(30.0*PI/180.0); y[2]=-Es*sin(30.0*PI/180.0); //third well
  PetscReal dist=sqrt(std::pow(e2-x[0],2.0)+std::pow(e3-y[0],2.0));
  unsigned int wellID=1; 
  for(unsigned int i=1; i<3; i++){
    if(dist>sqrt(pow(e2-x[i],2.0)+pow(e3-y[i],2.0))){
      dist=sqrt(pow(e2-x[i],2.0)+pow(e3-y[i],2.0));
      wellID=i+1;
    }
  }
#elif DIM==2
  //new strain metrics
  PetscReal e1=(E[0][0]+E[1][1]);
  PetscReal e2=(E[0][0]-E[1][1]);
  PetscReal e3=E[0][1];
  
  //compute distance to nearest well
  PetscReal Es=user->Es;
  PetscReal dist=e2-Es;
  unsigned int wellID=1;
#endif
 
  //store L2 projection residual
  const PetscReal (*N) = (PetscReal (*)) p->shape[0];;  
  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      PetscReal val=0.0;
      switch (d1) {
      case 0:
	val=e2; break;
      case 1:
	val=e3; break;
      case 3: //only in 3D with concentration
	val=dist; break;
      case 2: //only in 3D 
	val=(PetscReal) wellID; break;
      } 
      R[n1*dof+d1] = N[n1]*val;
      for(int n2=0; n2<nen; n2++){
	for(int d2=0; d2<dof; d2++){
	  PetscReal val2=0.0;
	  if (d1==d2) {val2 = N[n1] * N[n2];}
	  K[n1*dof*nen*dof + d1*nen*dof + n2*dof + d2] =val2;
	}
      }
    }
  }
  return 0;
}

PetscErrorCode ProjectSolution(IGA iga, PetscInt step, AppCtxKSP *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;
  //Setup linear system for L2 Projection
  Mat A;
  Vec x,b;
  ierr = IGACreateMat(iga,&A);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&x);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&b);CHKERRQ(ierr);
  ierr = IGASetFormSystem(iga,E22System,user);CHKERRQ(ierr);
  ierr = IGAComputeSystem2(iga,A,b);CHKERRQ(ierr);

  //Solver
  KSP ksp;
  ierr = IGACreateKSP(iga,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
  //write solution
  char filename[256];
  sprintf(filename,"./outE%d.dat",step);
  ierr = IGAWriteVec(iga,x,filename);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  PetscFunctionReturn(0); 
}

typedef struct {
  PetscReal ux, uy;
#if DIM==3  
  PetscReal uz;
#endif
} Field;

PetscErrorCode FormInitialCondition(IGA iga, Vec U, AppCtx *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;
  std::srand(5);
  DM da;
  ierr = IGACreateNodeDM(iga,DIM,&da);CHKERRQ(ierr);
#if DIM==3
  Field ***u;
  ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
  PetscInt i,j,k;
  for(i=info.xs;i<info.xs+info.xm;i++){
    for(j=info.ys;j<info.ys+info.ym;j++){
      for(k=info.zs;k<info.zs+info.zm;k++){
	u[k][j][i].ux=0.0;
	u[k][j][i].uy=0.0;
	u[k][j][i].uz=0.0;
      }
    }
  }
#elif DIM==2
  Field **u;
  ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
  PetscInt i,j;
  for(i=info.xs;i<info.xs+info.xm;i++){
    for(j=info.ys;j<info.ys+info.ym;j++){
      u[j][i].ux=0.0;
      u[j][i].uy=0.0;
    }
  }
#endif
  ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
  ierr = DMDestroy(&da);;CHKERRQ(ierr); 
  PetscFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "OutputMonitor"
PetscErrorCode OutputMonitor(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;
  char           filename[256];
  sprintf(filename,"./outU%d.dat",it_number);
  ierr = IGAWriteVec(user->iga,U,filename);CHKERRQ(ierr);
  //PetscPrintf(PETSC_COMM_WORLD,"\nProject Solution\n");
  ProjectSolution(user->iga, it_number, user->appCtxKSP);

  //Set load parameter
  double dVal=user->Es*c_time;
#if DIM==3
  ierr = IGASetBoundaryValue(user->iga,0,0,1,dVal);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user->iga,0,1,1,-dVal);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user->iga,1,0,0,dVal);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user->iga,1,1,0,-dVal);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user->iga,2,0,2,0.0);CHKERRQ(ierr);  
#elif DIM==2 
  ierr = IGASetBoundaryValue(user->iga,0,0,1,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user->iga,0,1,1,-dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user->iga,1,0,0,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user->iga,1,1,0,-dVal);CHKERRQ(ierr);  
#endif
  PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: it_number:%u, c_time:%12.6e, load:%.2e\n", it_number, c_time,dVal);
  PetscFunctionReturn(0);
}

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);
  double startTime = MPI_Wtime();
 
  /* Define simulation specific parameters */
  AppCtx user; AppCtxKSP userKSP;

  //problem parameters
  user.c=1.0;
  user.Es=0.01;
  user.Ed=1.0;
#if DIM==3
  user.E4=1.5*user.Ed/pow(user.Es,4.0);
  user.E3=-user.Ed/pow(user.Es,3.0);
  user.E2=-1.5*user.Ed/pow(user.Es,2.0);
#elif DIM==2
  user.E4= user.Ed/pow(user.Es,4.0);
  user.E3= 0.0;
  user.E2=-2.0*user.Ed/pow(user.Es,2.0);
#endif
  user.Eii=user.Ed/pow(user.Es,2.0);
  user.Eij=user.Ed/pow(user.Es,2.0);
  user.El=0.1;
  user.Eg=pow(user.El,2.0)*user.Ed/pow(user.Es,2.0);
  user.dt=0.01;
  user.C=5.0;

  /* Set discretization options */
  PetscInt nsteps = 100;
  PetscInt N=10, p=2, C=PETSC_DECIDE, resStep=0;
  PetscBool output = PETSC_FALSE; 
  PetscBool monitor = PETSC_FALSE; 
  char filePrefix[PETSC_MAX_PATH_LEN] = {0};
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","CahnHilliard2D Options","IGA");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-N","number of elements (along one DIMension)",__FILE__,N,&N,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-p","polynomial order",__FILE__,p,&p,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-C","global continuity order",__FILE__,C,&C,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-file_prefix","File Prefix",__FILE__,filePrefix,filePrefix,sizeof(filePrefix),PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-ch_output","Enable output files",__FILE__,output,&output,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-ch_monitor","Compute and show statistics of solution",__FILE__,monitor,&monitor,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-nsteps","Number of load steps to take",__FILE__,nsteps,&nsteps,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-dt","time step",__FILE__,user.dt,&user.dt,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-res_step","Restart Step",__FILE__,resStep,&resStep,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  if (C == PETSC_DECIDE) C = p-1;
  //PetscPrintf(PETSC_COMM_WORLD,"\nLambda_u value is: %8.2e\n",user.lambda_u);

  //
  if (p < 2 || C < 0) /* Problem requires a p>=2 C1 basis */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Problem requires minimum of p = 2");
  if (p <= C)         /* Check C < p */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Discretization inconsistent: polynomial order must be greater than degree of continuity");
  
  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDim(iga,DIM);CHKERRQ(ierr);
  ierr = IGASetDof(iga,DIM);CHKERRQ(ierr);

  IGAAxis axis0;
  ierr = IGAGetAxis(iga,0,&axis0);CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis0,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis0,N,0.0,1.0,C);CHKERRQ(ierr);
  user.he=1.0/N; //set he=L/N, by selecting the correct problem length L
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
  double dVal=user.Es*0.01;
#if DIM==3
  ierr = IGASetBoundaryValue(iga,0,0,1,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,0,1,1,-dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,0,0,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,1,0,-dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,2,0,2,0.0);CHKERRQ(ierr);  
#elif DIM==2 
  ierr = IGASetBoundaryValue(iga,0,0,1,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,0,1,1,-dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,0,0,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,1,0,-dVal);CHKERRQ(ierr);  
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
  ierr = TSSetDuration(ts,100,1.0);CHKERRQ(ierr);
  ierr = TSSetTime(ts,0.0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,user.dt);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,OutputMonitor,&user,NULL);CHKERRQ(ierr);
  
  /*
  SNES snes;
  TSGetSNES(ts,&snes);
  SNESLineSearch ls;
  SNESGetLineSearch(snes,&ls);
  SNESLineSearchSetType(ls,SNESLINESEARCHBT);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  SNESLineSearchView(ls,NULL);
  */
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

