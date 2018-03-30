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
    user.time += user.dt;
    ierr = StepUpdate<dim>(++step,user.time,*(user.U),user);
    ierr = VecCopy(*user.Up, *user.Upp);CHKERRQ(ierr);
    ierr = VecCopy(*user.U, *user.Up);CHKERRQ(ierr);
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
