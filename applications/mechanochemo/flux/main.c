#include <math.h> 
extern "C" {
#include "petiga.h"
}
#include "../../../include/fields.h"
#include "../../../include/petigaksp2.h"
//general parameters
#define DIM 2
#define ADSacado
#define numVars 27
//physical parameters
#define C4 1.0
#define C3 1.0
#define C2 1.0
#define E4 1.0*(c)
#define E3 1.0*(c)
#define E2 1.0*(c)
#define E4_c 1.0*(c)
#define E3_c 1.0*(c)
#define E2_c 1.0*(c)
#define E4_cc 1.0*(c)
#define E3_cc 1.0*(c)
#define E2_cc 1.0*(c)
#define Eii 1.0
#define Eij 1.0
#define Cl 1.0 //Clambda - constant for gradC.gradC
#define El 1.0 //ELambda - constant for gradE.gradE
#define Gl 1.0 //Glambda - constant for abs(gradC.gradE)
//stress
#define PiJ ((4*E4*e2*e2*e2 + 3*E3*e2*e2 + 2*E2*e2)*e2_FiJ + (e2_1+c_1/2)*e2_1_FiJ + (e2_2+c_2/2)*e2_2_FiJ)
#define BetaiJK ((e2_1+c_1/2)*e2_1_FiJK + (e2_2+c_2/2)*e2_2_FiJK)
//chemistry
#define mu_c (12*C4*c*c+6*C3*c+2*C2+E4_cc*e2*e2*e2*e2+E3_cc*e2*e2*e2+E2_cc*e2*e2)
#define mu_e2 (4*E4_c*e2*e2*e2+3*E3_c*e2*e2+2*E2_c*e2)
//boundary conditions
#define FLUX 0
#define flux 1.0
#define uDirichlet 1.0
//other variables
#define NVal 100
#define Es 0.01
#define DVal 1.0
#define CVal 1.0
#define gamma 1.0
#define cbar 0.05
//time stepping
#define dtVal 1.0e-7

//includes
#include "../../../include/appctx.h"
#include "../../../include/output.h"
#include "../../../include/evaluators.h"
#include "../../../include/initialConditions.h"
#include "../../../include/init.h"

//residual function implementation
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
  T c_1=cx[0], c_2=cx[1]; 
  
  //displacement field variable
  T u[DIM], ux[DIM][DIM], uxx[DIM][DIM][DIM];
  computeField<T,DIM,DIM+1>(VECTOR,0,p,&U[0],&u[0],&ux[0][0],&uxx[0][0][0]);

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
      E[I][J] = -0.5*(I==J);
      for (unsigned int k=0; k<DIM; k++){
	E[I][J] += 0.5*F[k][I]*F[k][J];
      }
    }
  }
  
  //new strain metrics (2D)
  T e1=(E[0][0]+E[1][1]);
  T e2=(E[0][0]-E[1][1]);
  T e3=E[0][1];
  T e2_1=0.0, e2_2=0.0; 
  for (unsigned int i=0; i<DIM; ++i){
    e2_1+=(F[i][0]*dF[i][0][0]-F[i][1]*dF[i][1][0]);
    e2_2+=(F[i][0]*dF[i][0][1]-F[i][1]*dF[i][1][1]);
  }
  //compute P and Beta
  T P[DIM][DIM], Beta[DIM][DIM][DIM];
  //
  for (unsigned int i=0; i<DIM; ++i){
    for (unsigned int J=0; J<DIM; ++J){
      T e1_FiJ=(F[i][0]*(0==J)+F[i][1]*(1==J));
      T e2_FiJ=(F[i][0]*(0==J)-F[i][1]*(1==J));
      T e3_FiJ=(F[i][1]*(0==J)+F[i][0]*(1==J))/2.0;
      T e2_1_FiJ=((0==J)*dF[i][0][0]-(1==J)*dF[i][1][0]);
      T e2_2_FiJ=((0==J)*dF[i][0][1]-(1==J)*dF[i][1][1]);
      //P
      P[i][J]=PiJ;

      //gradient terms
      for (unsigned int K=0; K<DIM; ++K){
	T e2_1_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(0==K);
	T e2_2_FiJK=(F[i][0]*(0==J)-F[i][1]*(1==J))*(1==K);
	//Beta
	Beta[i][J][K]=BetaiJK;
      }
    }
  }
  
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
      Rc += N[a] * (c-c0)*(1.0/dt2);
      // grad(Na) . M*(dmuc*grad(C)+dmue2*grad(e2))
      double laplace_N = N2[0][0] + N2[1][1];
      Rc += DVal*mu_c*(N1[0]*c_1+N1[1]*c_2) + DVal*mu_e2*(N1[0]*e2_1+N1[1]*e2_2);
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


int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);
 
  //application context objects and parameters
  AppCtx user; AppCtxKSP userKSP;
  user.dt=dtVal;
  user.he=1.0/NVal; 
  PetscInt p=2;

  //iga
  IGA iga;
  user.iga = iga;
  Vec U,U0;
  user.appCtxKSP=&userKSP;
  userKSP.U0=&U;

  //initialize
  init(iga, user, NVal, p, 0, U, U0);

  //boundary conditons
  double dVal=uDirichlet;
#if bcVAL==0
  //shear BC
  ierr = IGASetBoundaryValue(iga,0,0,1,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,0,1,1,-dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,0,0,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,1,0,-dVal);CHKERRQ(ierr);
#elif bcVAL==1
  //free BC
  ierr = IGASetBoundaryValue(iga,0,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,0,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,0,1,0,dVal);CHKERRQ(ierr);  
#elif bcVAL==2
  //fixed BC
  ierr = IGASetBoundaryValue(iga,0,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,0,0,1,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,0,0,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,1,0,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,1,1,0,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,1,1,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,0,1,0,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,0,1,1,0.0);CHKERRQ(ierr);  
#endif

  //finalize
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&U0);CHKERRQ(ierr);
  ierr = IGADestroy(&iga);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

