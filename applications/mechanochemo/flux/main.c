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
#define Gl 0.0 //Glambda - constant for abs(gradC.gradE)
//stress
#define PiJ (2*Eii*e1*e1_FiJ + 2*Eij*e6*e6_FiJ + (4*E4*e2*e2*e2 + 3*E3*e2*e2 + 2*E2*e2)*e2_FiJ + (El*e2_1+Gl*c_1/2)*e2_1_FiJ + (El*e2_2+Gl*c_2/2)*e2_2_FiJ)
#define BetaiJK ((El*e2_1+Gl*c_1/2)*e2_1_FiJK + (El*e2_2+Gl*c_2/2)*e2_2_FiJK)
//chemistry
#define mu_c (12*C4*c*c+6*C3*c+2*C2+E4_cc*e2*e2*e2*e2+E3_cc*e2*e2*e2+E2_cc*e2*e2)
#define mu_e2 (4*E4_c*e2*e2*e2+3*E3_c*e2*e2+2*E2_c*e2)
//boundary conditions
#define bcVAL 1
#define FLUX 0
#define flux 1.0
#define uDirichlet 0.01
//other variables
#define NVal 100
#define Es 0.01
#define DVal 1.0
#define CVal 1.0
#define gamma 1.0
#define cbar 0.05
//time stepping
#define dtVal 1.0e-7

//other problem headers
#include "../../../include/appctx.h"
#include "../../../include/output.h"
#include "../../../include/evaluators.h"
#include "../../../include/initialConditions.h"
#include "../../../include/init.h"
//physics header
#include "../../../src/mechanochemo/twoD/mechanoChemo2D.h"


int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);
 
  //application context objects and parameters
  AppCtx user; AppCtxKSP userKSP;
  user.appCtxKSP=&userKSP;
  user.dt=dtVal;
  user.he=1.0/NVal; 
  PetscInt p=2;

  //initialize
  IGA iga;
  Vec U,U0;
  user.iga = iga;
  userKSP.U0=&U;
  user.U=&U;
  user.U0=&U0;

  PetscPrintf(PETSC_COMM_WORLD,"initializing...\n");
  init(user, userKSP, NVal, p);

  //Dirichlet boundary conditons for mechanics
  PetscPrintf(PETSC_COMM_WORLD,"applying bcs...\n");
  
  double dVal=uDirichlet;
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
  //run
  PetscPrintf(PETSC_COMM_WORLD,"running...\n");
  ierr = TSSolve(ts,*user.U);CHKERRQ(ierr);

  //finalize
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = VecDestroy(user.U);CHKERRQ(ierr);
  ierr = VecDestroy(user.U0);CHKERRQ(ierr);
  ierr = IGADestroy(&user.iga);CHKERRQ(ierr);
  
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

