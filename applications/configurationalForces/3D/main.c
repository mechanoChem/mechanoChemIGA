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
#include "timeStepping.h"

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);
 
  //Initialize objects and parameters
  AppCtx user;
  Vec U,U0;
  TS ts;
	ierr = setup<DIM>(user,&U,&U0,ts);

  //Dirichlet boundary conditons for mechanics
  PetscPrintf(PETSC_COMM_WORLD,"applying bcs...\n");
	ierr = boundaryConditions<DIM>(user,0.);

  //time stepping
	ierr = timeStepping(user,ts);

  //set snes convergence tests
	ierr = setConvergenceTest(user,ts);

  //run
  PetscPrintf(PETSC_COMM_WORLD,"running...\n");
  ierr = TSSolve(ts,*user.U);CHKERRQ(ierr);

  //finalize
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

