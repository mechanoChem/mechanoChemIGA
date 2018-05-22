#include <math.h> 
//extern "C" {
#include "petiga.h"
//}
//core functions
#include "coreFunctions.h"
//application specific functions
#include "userFunctions.h"

template <int dim>
PetscErrorCode Run(){
  
  PetscErrorCode  ierr;
  AppCtx<dim> user;
  Vec U,Up,Upp;
  SNES snes;
  unsigned int counter = 0;
  PetscInt n_iter;
  SNESConvergedReason conv_reason;

  //Setup structures, initial conditions, boundary conditions
  ierr = Setup<dim>(user,&U,&Up,&Upp,snes);

  //run
  PetscPrintf(PETSC_COMM_WORLD,"running...\n");
  user.time = user.RESTART_TIME;
  PetscInt step = user.RESTART_IT;
  ierr = StepUpdate<dim>(step,user.time,*(user.U),user);
  while (user.time < user.totalTime){
    PetscPrintf(PETSC_COMM_WORLD,"Step %i, dt %g, time %g\n",step,user.dt,user.time);
    ierr = SNESSolve(snes,NULL,*(user.U));CHKERRQ(ierr);
    ierr = SNESGetIterationNumber(snes,&n_iter);
    ierr = SNESGetConvergedReason(snes,&conv_reason);
    //PetscPrintf(PETSC_COMM_WORLD,"Convergence reason %i\n",conv_reason);
    //Start putting in some infrastructure for adaptive time stepping...
    if(n_iter < 10 && conv_reason > 0){
      //If converged fast enough, move on
      if(n_iter < 3){
	counter++;
	//If converged fast, keep track of it
	if(counter >= 5){
	  //If converged really fast for multiple iterations, scale up time step
	  PetscPrintf(PETSC_COMM_WORLD,"Doubling time step...\n");
	  user.dt *= 2.;
	  counter = 0;
	}
      }
      user.time += user.dt;
      ierr = StepUpdate<dim>(++step,user.time,*(user.U),user);
      ierr = VecCopy(*user.Up, *user.Upp);CHKERRQ(ierr);
      ierr = VecCopy(*user.U, *user.Up);CHKERRQ(ierr);
    }
    else{
      //If it didn't converge fast enough, scale back time step and try again.
      PetscPrintf(PETSC_COMM_WORLD,"Halving time step and trying again...\n");
      user.dt *= 0.5;
      //Reset the current to the previous converged U
      ierr = VecCopy(*user.Up, *user.U);CHKERRQ(ierr);
    }
  }

  //finalize
  ierr = SNESDestroy(user.snes);CHKERRQ(ierr);
  ierr = VecDestroy(user.U);CHKERRQ(ierr);
  ierr = VecDestroy(user.Up);CHKERRQ(ierr);
  ierr = VecDestroy(user.Upp);CHKERRQ(ierr);
  ierr = IGADestroy(&user.iga);CHKERRQ(ierr);
  ierr = IGADestroy(&user.igaProject);CHKERRQ(ierr);  

  return 0;
}

template PetscErrorCode Run<2>();
template PetscErrorCode Run<3>();
