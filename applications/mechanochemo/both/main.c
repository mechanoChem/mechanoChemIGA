#include <math.h> 
extern "C" {
#include "petiga.h"
}
#include "../../../include/fields.h"
//general parameters
#define DIM 2
#define GridScale 0.1
#define ADSacado
#define numVars 27
//physical parameters
#define Cs 1.0
#define Cd -2.0
#define C4 (-16*Cd/std::pow(Cs,4))
#define C3 (32*Cd/std::pow(Cs,3))
#define C2 (-16*Cd/std::pow(Cs,2))
#define Es 0.1
#define Ed -0.1
#define E4 (-Ed/std::pow(Es,4))
#define E3 0.0
#define E2 (2*Ed/std::pow(Es,2))*((2*c-Cs)/Cs)
#define E4_c 0.0
#define E3_c 0.0
#define E2_c (2*Ed/std::pow(Es,2))*(2.0/Cs)
#define E4_cc 0.0
#define E3_cc 0.0
#define E2_cc 0.0
#define Eii (-2*Ed/std::pow(Es,2))
#define Eij (-2*Ed/std::pow(Es,2))
#define Cl  (0.01*GridScale*GridScale) //Clambda - constant for gradC.gradC
//#define El 0.0001 //**ELambda - constant for gradE.gradE
#define Gl 0.0  //Glambda - constant for abs(gradC.gradE)
//stress
#define PiJ (2*Eii*e1*e1_FiJ + 2*Eij*e6*e6_FiJ + (4*E4*e2*e2*e2 + 3*E3*e2*e2 + 2*E2*e2)*e2_FiJ + (El*e2_1+Gl*c_1/2)*e2_1_FiJ + (El*e2_2+Gl*c_2/2)*e2_2_FiJ)
#define BetaiJK ((El*e2_1+Gl*c_1/2)*e2_1_FiJK + (El*e2_2+Gl*c_2/2)*e2_2_FiJK)
//chemistry
#define mu_c (12*C4*c*c+6*C3*c+2*C2+E4_cc*e2*e2*e2*e2+E3_cc*e2*e2*e2+E2_cc*e2*e2)
#define mu_e2 (4*E4_c*e2*e2*e2+3*E3_c*e2*e2+2*E2_c*e2)
//boundary conditions
//#define bcVAL 1 //**
//#define FLUX 3 //**
#define flux 10.0
#define uDirichlet 0.001
//other variables
//#define NVal 60 //**
#define DVal (1.0*GridScale*GridScale)
#define CVal 5.0
#define gamma 1.0
//time stepping
//#define dtVal 1.0e-6 //**
#define skipOutput 10

//other problem headers
#include "../../../include/appctx.h"
#include "../../../include/output.h"
#include "../../../include/evaluators.h"
#include "../../../include/initialConditions.h"
#include "../../../include/init.h"
//physics header
#include "../../../src/mechanochemo/twoD/mechanoChemo2D.h"

//snes convegence test
PetscErrorCode SNESConvergedTest(SNES snes, PetscInt it,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx){
  AppCtx *user  = (AppCtx*) ctx;
  //custom test
  if (NVal>=300){
    if (it>50){
      *reason = SNES_CONVERGED_ITS;
      return(0);
    }
  }
  //default test
  PetscFunctionReturn(SNESConvergedDefault(snes,it,xnorm,snorm,fnorm,reason,ctx));
}

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);
 
  //application context objects and parameters
  AppCtx user;
  user.dt=dtVal;
  user.he=GridScale*1.0/NVal; 
  PetscInt p=2;

  //initialize
  IGA iga;
  Vec U,U0;
  user.iga = iga;
  user.U0=&U0;
  user.U=&U;

  PetscPrintf(PETSC_COMM_WORLD,"initializing...\n");
  init(user, NVal, p);

  //Dirichlet boundary conditons for mechanics
  PetscPrintf(PETSC_COMM_WORLD,"applying bcs...\n");
  
  double dVal=uDirichlet*GridScale;
#if bcVAL==0
  //shear BC
  ierr = IGASetBoundaryValue(user.iga,0,0,1,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,1,1,-dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,0,0,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,1,0,-dVal);CHKERRQ(ierr);
#elif bcVAL==1
  //free BC
  ierr = IGASetBoundaryValue(user.iga,0,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,0,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,0,1,0,dVal);CHKERRQ(ierr);  
#elif bcVAL==2
  //fixed BC
  ierr = IGASetBoundaryValue(user.iga,0,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,0,1,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,0,0,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,1,0,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,1,1,0,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,1,1,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,0,1,0,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,1,1,0.0);CHKERRQ(ierr);  
#endif 

  //time stepping
  TS ts;
  ierr = IGACreateTS(user.iga,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,100000,1.0);CHKERRQ(ierr);
  ierr = TSSetTime(ts,0.0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,user.dt);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,OutputMonitor<DIM>,&user,NULL);CHKERRQ(ierr);  
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  //set snes convergence tests
  SNES snes;
  ierr = TSGetSNES(ts,&snes); CHKERRQ(ierr);
  ierr = SNESSetConvergenceTest(snes,SNESConvergedTest,&user,NULL); CHKERRQ(ierr);

  //run
  PetscPrintf(PETSC_COMM_WORLD,"running...\n");
  user.ts=&ts;
  ierr = TSSolve(ts,*user.U);CHKERRQ(ierr);

  //finalize
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = VecDestroy(user.U);CHKERRQ(ierr);
  ierr = VecDestroy(user.U0);CHKERRQ(ierr);
  ierr = IGADestroy(&user.iga);CHKERRQ(ierr);
  
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

