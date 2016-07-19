#include <math.h> 
//extern "C" {
#include "petiga.h"
//}
//generic functions
#include "../../../include/genericHeaders.h"
//application specific functions
#include "../applicationHeaders.h"
//physics functions
#include "../../../src/physicsHeaders.h"

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);
 
  //Initialize objects and parameters
  AppCtx user;
	ierr = defineParameters(user);
  Vec U,U0;
  TS ts;

	const unsigned int dim = 3;
	ierr = setup<dim>(user,&U,&U0,ts);

	//Dirichlet boundary conditons for mechanics
	PetscPrintf(PETSC_COMM_WORLD,"applying bcs...\n");
	ierr = boundaryConditions<dim>(user,0.);

	//time stepping
	ierr = timeStepSetup<dim>(user,ts);

  //set snes convergence tests
  ierr = setConvergenceTest(user,ts);

  //run
  PetscPrintf(PETSC_COMM_WORLD,"running...\n");
  ierr = TSSolve(ts,*user.U);CHKERRQ(ierr);

  //finalize
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
