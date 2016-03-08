#include <math.h> 
//extern "C" {
#include "petiga.h"
//}

//model parameters
#define DIM 3
#define GridScale 1.0
#define ADSacado //enables Sacado for AD instead of Adept
#define numVars 162 //162 in 3D, 36 in 2D for 2*DIM dof

//generic headers
#include "../../../include/fields.h"
#include "../../../include/appctx.h"
#include "../../../include/evaluators.h"
#include "../../../include/init.h"

//2D physical parameters
#if DIM == 2
  #define Es 0.1
  #define Ed 1.0
  #define E4 (Ed/(std::pow(Es,4)))
  #define E2 (Ed/(std::pow(Es,2)))*(2*c-1.0)
  #define Eii (Ed/std::pow(Es,2))
  #define Eij (Ed/std::pow(Es,2))
  #define El .025 //**ELambda - constant for gradE.gradE
#elif DIM == 3
  #define Es 1.0e-1//1.0e-2
  #define Ed -1.0
  #define E4 (-3*Ed/(2*std::pow(Es,4)))
  #define E3 (-Ed/(std::pow(Es,3)))*(c)
  #define E2 (3*Ed/(2*std::pow(Es,2)))*(2*c-1.0)
  #define Eii (-1.5*Ed/std::pow(Es,2))
  #define Eij (-1.5*Ed/std::pow(Es,2))
  #define El 1 //**ELambda - constant for gradE.gradE
#endif
//material model (stress expressions)
//non-gradient St-Venant Kirchoff model with cubic crystal material parameters
#define mu 1e5
#define betaC 1e5
#define alphaC 2e5//(betaC + 2*mu) for isotropic materials
#if DIM == 2
  #define PiJ ((alpha[J]-2*mu-beta[J][J])*F[i][J]*E[J][J] + (beta[J][0]*E[0][0]+beta[J][1]*E[1][1])*F[i][J] + 2*mu*(F[i][0]*E[0][J]+F[i][1]*E[1][J])) //2D
#elif DIM == 3
  #define PiJ ((alpha[J]-2*mu-beta[J][J])*F[i][J]*E[J][J] + (beta[J][0]*E[0][0]+beta[J][1]*E[1][1]+beta[J][2]*E[2][2])*F[i][J] + 2*mu*(F[i][0]*E[0][J]+F[i][1]*E[1][J]+F[i][2]*E[2][J])) //3D
#endif
#define BetaiJK (0.0)
//non-gradient St-Venant Kirchoff model with lambda=mu=1
//#define PiJ ((E[0][0]+E[1][1]+E[2][2])*F[i][J] + 2*(F[i][0]*E[0][J]+F[i][1]*E[1][J]+F[i][2]*E[2][J]))
//#define BetaiJK (0.0)
//gradient model
#if DIM == 2
  #define P0iJ (2*Eii*e1*e1_chiiJ + 2*Eij*e6*e6_chiiJ + (4*E4*e2*e2*e2 - 4*E2*e2)*e2_chiiJ + 2*El*El*Eii*(e2_1*e2_1_chiiJ + e2_2*e2_2_chiiJ))//2D
  #define Beta0iJK  (2*El*El*Eii*(e2_1*e2_1_chiiJK + e2_2*e2_2_chiiJK)) //2D
#elif DIM == 3
  #define P0iJ (2*Eii*e1*e1_chiiJ + 2*Eij*e4*e4_chiiJ + 2*Eij*e5*e5_chiiJ + 2*Eij*e6*e6_chiiJ + (2*E2*e2-6*E3*e2*e3+4*E4*e2*(e2*e2+e3*e3))*e2_chiiJ + (2*E2*e3+3*E3*(e3*e3-e2*e2)+4*E4*e3*(e2*e2+e3*e3))*e3_chiiJ + 2*El*(e2_1*e2_1_chiiJ + e2_2*e2_2_chiiJ + e2_3*e2_3_chiiJ + e3_1*e3_1_chiiJ + e3_2*e3_2_chiiJ + e3_3*e3_3_chiiJ))
  #define Beta0iJK  2*El*(e2_1*e2_1_chiiJK + e2_2*e2_2_chiiJK + e2_3*e2_3_chiiJK + e3_1*e3_1_chiiJK + e3_2*e3_2_chiiJK + e3_3*e3_3_chiiJK)
#endif
//boundary conditions
#define bcVAL 3 //**
#define uDirichlet 0.1
//other variables
#define NVal 10//**
//time stepping
#define dtVal 0.00001 //** // used to set load parameter..so 0<dtVal<1
#define skipOutput 100
//restart file
#define RESTART_IT 400
#define RESTART_TIME 0.004

//physics headers
#include "../../../src/configurationalForces/model.h"
#include "../../../src/configurationalForces/initialConditions.h"
#include "../../../src/configurationalForces/output.h"

//snes convegence test
PetscErrorCode SNESConvergedTest(SNES snes, PetscInt it,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx){
  AppCtx *user  = (AppCtx*) ctx;
  PetscPrintf(PETSC_COMM_WORLD,"xnorm:%12.6e snorm:%12.6e fnorm:%12.6e\n",xnorm,snorm,fnorm);  
  //custom test
  if (NVal>=15){
    if (it>500){
      *reason = SNES_CONVERGED_ITS;
      return(0);
    }
  }
  //default test
  PetscFunctionReturn(SNESConvergedDefault(snes,it,xnorm,snorm,fnorm,reason,ctx));
}


PetscErrorCode SNESConvergedTest_Interactive(SNES snes, PetscInt it,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx)
{
  AppCtx *user  = (AppCtx*) ctx;
  /* Vec R;                                                                                                                
  SNESGetFunction(snes, &R, 0, 0);                                                                                         
  PetscReal normC, normU, normUx, normUy;                                                                                  
  //VecNorm(R, NORM_2, &norm);                                                                                             
  VecSetBlockSize(R,DIM+1);                                                                                                
  VecStrideNorm(R,0,NORM_2,&normUx);                                                                                       
  VecStrideNorm(R,1,NORM_2,&normUy);                                                                                       
  VecStrideNorm(R,2,NORM_2,&normC);                                                                                        
  normU=sqrt(pow(normUx,2)+pow(normUy,2));                                                                                 
  //                                                                                                                       
  if (it==0) user->f0Norm=fnorm;*/
  PetscPrintf(PETSC_COMM_WORLD,"xnorm:%12.6e snorm:%12.6e fnorm:%12.6e\n",xnorm,snorm,fnorm);
  //if ((it>10) && (fnorm/user->f0Norm<1.0e-4)) {                                                                          
  if ((it>5) && (fnorm<1.0e-7)) {
    PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: since it>10 forcefully setting convergence. \n");
    *reason = SNES_CONVERGED_FNORM_ABS;
    return(0);
  }
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
  const unsigned int DOF=2*DIM;

  //initialize
  IGA iga;
  Vec U,U0;
  user.iga = iga;
  user.U0=&U0;
  user.U=&U;
  user.lambda=1.;

  PetscPrintf(PETSC_COMM_WORLD,"initializing...\n");
  init<DOF>(user, NVal, p);

  //initial conditions
  ierr = IGACreateVec(user.iga,user.U);CHKERRQ(ierr);
  ierr = IGACreateVec(user.iga,user.U0);CHKERRQ(ierr);
#if RESTART_IT==0
  ierr = FormInitialCondition(user.iga, *user.U0, &user);CHKERRQ(ierr);
#else
  char           filename[256];
  sprintf(filename,"outU%d.dat",RESTART_IT);  
  ierr = IGAReadVec(user.iga,*user.U0,filename);CHKERRQ(ierr); //Read in vector to restart at step RESTART_IT
#endif 
  ierr = VecCopy(*user.U0, *user.U);CHKERRQ(ierr);

  //Dirichlet boundary conditons for mechanics
  PetscPrintf(PETSC_COMM_WORLD,"applying bcs...\n");
  double dVal=uDirichlet*GridScale;
  //	dVal *= dtVal; //For load stepping
#if bcVAL==0 //unchanged by Greg
  //shear BC
  ierr = IGASetBoundaryValue(user.iga,0,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,1,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,0,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,1,1,1,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,2,0,2,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,2,1,2,0.0);CHKERRQ(ierr); 

  ierr = IGASetBoundaryValue(user.iga,0,0,4,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,1,4,-dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,0,3,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,1,3,-dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,2,0,5,0.0);CHKERRQ(ierr);  
#elif bcVAL==1
  //free BC
  ierr = IGASetBoundaryValue(user.iga,0,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,1,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,0,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,1,1,1,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,2,0,2,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,2,1,2,0.0);CHKERRQ(ierr); 

  ierr = IGASetBoundaryValue(user.iga,0,0,3,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,0,4,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,2,0,5,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,1,3,dVal);CHKERRQ(ierr);  
#elif bcVAL==2 //unchanged by Greg
  //fixed BC
  ierr = IGASetBoundaryValue(user.iga,0,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,1,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,0,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,1,1,1,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,2,0,2,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,2,1,2,0.0);CHKERRQ(ierr); 

  ierr = IGASetBoundaryValue(user.iga,0,0,3,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,0,4,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,0,0,5,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,0,1,3,dVal);CHKERRQ(ierr); 
#elif bcVAL==3
  //bending BC 
  //ierr = IGASetBoundaryValue(user.iga,2,0,2,0.0);CHKERRQ(ierr);  
  //ierr = IGASetBoundaryValue(user.iga,2,1,2,0.0);CHKERRQ(ierr); 

  ierr = IGASetBoundaryValue(user.iga,0,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,0,1,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,0,2,0.0);CHKERRQ(ierr); 

  ierr = IGASetBoundaryValue(user.iga,0,1,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,1,1,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,1,2,0.0);CHKERRQ(ierr);

  /*//Additional
  ierr = IGASetBoundaryValue(user.iga,1,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,0,1,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,0,2,0.0);CHKERRQ(ierr); 

  ierr = IGASetBoundaryValue(user.iga,1,1,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,1,1,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,1,2,0.0);CHKERRQ(ierr);

  ierr = IGASetBoundaryValue(user.iga,2,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,2,0,1,0.0);CHKERRQ(ierr);  
  //ierr = IGASetBoundaryValue(user.iga,2,0,2,0.0);CHKERRQ(ierr); 

  ierr = IGASetBoundaryValue(user.iga,2,1,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,2,1,1,0.0);CHKERRQ(ierr);  
  //ierr = IGASetBoundaryValue(user.iga,2,1,2,0.0);CHKERRQ(ierr);*/

  //plane strain
  //ierr = IGASetBoundaryValue(user.iga,2,0,5,0.0);CHKERRQ(ierr);
  //ierr = IGASetBoundaryValue(user.iga,2,1,5,0.0);CHKERRQ(ierr); 

  ierr = IGASetBoundaryValue(user.iga,0,0,3,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,0,4,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,0,0,5,0.0);CHKERRQ(ierr);  

  ierr = IGASetBoundaryValue(user.iga,0,1,3,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,0,1,4,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,0,1,5,0.0);CHKERRQ(ierr);   
#endif 

  //Set mat type
  ierr = IGASetMatType(user.iga,MATAIJ);CHKERRQ(ierr);
  //ierr = IGASetMatType(user.iga,MATIS);CHKERRQ(ierr);

  //time stepping
  TS ts;
  ierr = IGACreateTS(user.iga,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,100001,1.0);CHKERRQ(ierr);
  ierr = TSSetTime(ts,RESTART_TIME);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,user.dt);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,OutputMonitor<DIM>,&user,NULL);CHKERRQ(ierr);  
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  //set snes convergence tests
  SNES snes;
  ierr = TSGetSNES(ts,&snes); CHKERRQ(ierr);
  //ierr = SNESSetConvergenceTest(snes,SNESConvergedTest,&user,NULL); CHKERRQ(ierr);
  ierr = SNESSetConvergenceTest(snes,SNESConvergedTest_Interactive,(void*)&user,NULL); CHKERRQ(ierr);

  //* //PC
  PC pc;
  KSP ksp;
  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);// */

  /* //PCBDDC
  //ierr = KSPSetType(ksp,KSPCG);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  //ierr = PCSetType(pc,PCBDDC);CHKERRQ(ierr);
  ierr = IGAPreparePCBDDC(iga,pc);CHKERRQ(ierr);// */

  //* //Superlu_dist
  ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
  ierr = PCFactorSetMatSolverPackage(pc,MATSOLVERSUPERLU_DIST);CHKERRQ(ierr);// */

  //run
  PetscPrintf(PETSC_COMM_WORLD,"running...\n");
  user.ts=&ts;
  ierr = TSSolve(ts,*user.U);CHKERRQ(ierr);

  //finalize
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

