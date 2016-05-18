#include <math.h> 
//extern "C" {
#include "petiga.h"
//}

//model parameters
#include "parameters.h"

//generic headers
#include "../../../include/fields.h"
#include "../../../include/appctx.h"
#include "../../../include/evaluators.h"
#include "../../../include/init.h"

//boundary conditions
#include "boundaryConditions.h"

//physics headers
#include "../../../src/configurationalForces/model.h"
#include "../../../src/configurationalForces/initialConditions.h"
#include "../../../src/configurationalForces/output.h"

//problem setup
#include "convergenceTest.h"
#include "setup.h"

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);
 
  //application context objects and parameters
  AppCtx user;

  //initialize
  Vec U,U0;
	ierr = setup<DIM>(user,&U,&U0);

  //Dirichlet boundary conditons for mechanics
  PetscPrintf(PETSC_COMM_WORLD,"applying bcs...\n");
	boundaryConditions<DIM>(user,0.);

  //Set mat type
  ierr = IGASetMatType(user.iga,MATAIJ);CHKERRQ(ierr); //For superlu_dist
  //ierr = IGASetMatType(user.iga,MATIS);CHKERRQ(ierr); //For PCBDDC

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

  //run
  PetscPrintf(PETSC_COMM_WORLD,"running...\n");
  user.ts=&ts;
  ierr = TSSolve(ts,*user.U);CHKERRQ(ierr);

  //finalize
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

