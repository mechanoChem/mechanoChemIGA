#include <math.h> 
//extern "C" {
#include "petiga.h"
//}
//parameters
#include "parameters.h"
//generic functions
#include "../../../include/genericHeaders.h"
//application specific functions
#include "applicationHeaders.h"
//physics functions
#include "../../../src/configurationalForces/physicsHeaders.h"

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);
 
  //Initialize objects and parameters
  AppCtx user;
	//ierr = readParameters(user);
  Vec U,U0;
  TS ts;
  ierr = setup<DIM>(user,&U,&U0,ts);

  //Dirichlet boundary conditons for mechanics
  PetscPrintf(PETSC_COMM_WORLD,"applying bcs...\n");
  ierr = boundaryConditions<DIM>(user,0.);

  //time stepping
  ierr = timeStepSetup(user,ts);

  //set snes convergence tests
  ierr = setConvergenceTest(user,ts);

  //run
  PetscPrintf(PETSC_COMM_WORLD,"running...\n");
  ierr = TSSolve(ts,*user.U);CHKERRQ(ierr);

  //finalize
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

