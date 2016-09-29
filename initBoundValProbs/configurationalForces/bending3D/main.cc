#include <math.h> 
//extern "C" {
#include "petiga.h"
//}
//generic functions
#include "utilsIGAHeaders.h"
//application specific functions
#include "IBVPHeaders.h"
//physics functions
#include "physicsHeaders.h"

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);
 
  //Initialize objects and parameters
  AppCtx user;
	ierr = defineParameters(user);
  Vec U,U0;
  TS ts;

	if(user.dim == 2){
		const unsigned int dim = 2;
		const unsigned int dof=2*dim;
		//Setup structures, initial conditions, boundary conditions
		ierr = setup<dim,dof>(user,&U,&U0,ts);
	}
	else if(user.dim == 3){
		const unsigned int dim = 3;
		const unsigned int dof=2*dim;
		//Setup structures, initial conditions, boundary conditions
		ierr = setup<dim,dof>(user,&U,&U0,ts);
	}

  //run
  PetscPrintf(PETSC_COMM_WORLD,"running...\n");
  ierr = TSSolve(ts,*user.U);CHKERRQ(ierr);

  //finalize
  ierr = TSDestroy(user.ts);CHKERRQ(ierr);
  ierr = VecDestroy(user.U);CHKERRQ(ierr);
  ierr = VecDestroy(user.U0);CHKERRQ(ierr);
  ierr = IGADestroy(&user.iga);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
