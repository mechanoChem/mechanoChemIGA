#include <math.h> 
//extern "C" {
#include "petiga.h"
//}
//core functions
#include "coreFunctions.h"
//application specific functions
#include "userFunctions.h"
#include "DNN.h"

#include <chrono>
//#include <time.h>
#include <random>

int main(int argc, char *argv[]) {
 
  //struct timespec start, finish;
  //clock_gettime(CLOCK_MONOTONIC, &start);
  auto start = std::chrono::system_clock::now();
  //srand(0);
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);
 
  //Initialize objects and parameters
  const unsigned int dim = 2;

  //////////////////

  AppCtx<dim> user;
  Vec U,Up,Upp;
  Vec Y1_c, Y1_eta, Y2_c, Y2_eta; //For IRK
  SNES snes;
  unsigned int counter = 0;
  PetscInt n_iter;
  SNESConvergedReason conv_reason;

  DNN freeEnergy;
  freeEnergy.reinit(2);
  user.parameters = &freeEnergy;

  //Setup structures, initial conditions, boundary conditions
  ierr = Setup<dim>(user,&U,&Up,&Upp,snes);

  int size;
  ierr = VecGetLocalSize(*(user.U),&size);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&Y1_c);CHKERRQ(ierr);
  ierr = VecSetSizes(Y1_c,size/4,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(Y1_c);CHKERRQ(ierr);
  ierr = VecDuplicate(Y1_c,&Y1_eta);CHKERRQ(ierr);
  ierr = VecDuplicate(Y1_c,&Y2_c);CHKERRQ(ierr);
  ierr = VecDuplicate(Y1_c,&Y2_eta);CHKERRQ(ierr);

  // Make sure initial conditions on Y2 are the same as Y1
  ierr = VecStrideGather(*(user.U),0,Y1_c,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecStrideGather(*(user.U),1,Y1_eta,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecStrideScatter(Y1_c,2,*user.U,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecStrideScatter(Y1_eta,3,*user.U,INSERT_VALUES);CHKERRQ(ierr);

  ierr = VecCopy(*user.U, *user.Up);CHKERRQ(ierr);
  ierr = VecCopy(*user.Up, *user.Upp);CHKERRQ(ierr);

  //run
  PetscPrintf(PETSC_COMM_WORLD,"running...\n");
  user.time = user.RESTART_TIME;
  PetscInt step = user.RESTART_IT;
  ierr = StepUpdate<dim>(step,user.time,*(user.U),user);

  //clock_gettime(CLOCK_MONOTONIC, &finish);
  //double elapsed = finish.tv_sec + 1.e-9*finish.tv_nsec - start.tv_sec - 1.e-9*start.tv_nsec;
  auto end = std::chrono::system_clock::now();
  auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  double elapsed = 1.e-3*elapsed_ms;
  PetscPrintf(PETSC_COMM_WORLD,"Time elapsed: %f s\n",elapsed);

  while (user.time < user.totalTime){
    PetscPrintf(PETSC_COMM_WORLD,"Step %i, dt %g, time %g\n",step,user.dt,user.time);
    ierr = SNESSolve(snes,NULL,*(user.U));CHKERRQ(ierr);
    ierr = SNESGetIterationNumber(snes,&n_iter);
    ierr = SNESGetConvergedReason(snes,&conv_reason);

    //Start putting in some infrastructure for adaptive time stepping...
    if(n_iter < 10 && conv_reason > 0){
      // Since we're doing IRK4 (Implicit Runge Kutta, 4th order (2 stage),
      // the two fields after solving are Y1 and Y2. We need to correct Y1 to be y_{n+1}
      // using y_{n+1} = y_n + \sqrt{3}*(Y2-Y1)

      // First, extract the 4 subvectors
      ierr = VecStrideGather(*(user.U),0,Y1_c,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecStrideGather(*(user.U),1,Y1_eta,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecStrideGather(*(user.U),2,Y2_c,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecStrideGather(*(user.U),3,Y2_eta,INSERT_VALUES);CHKERRQ(ierr);

      // Now, take \sqrt{3}*(-Y1 + Y2) -> Y1
      ierr = VecAXPBY(Y1_c,std::sqrt(3.),-std::sqrt(3.),Y2_c);CHKERRQ(ierr);
      ierr = VecAXPBY(Y1_eta,std::sqrt(3.),-std::sqrt(3.),Y2_eta);CHKERRQ(ierr);

      // Copy u_n back into u_{n+1}, then add Y1
      ierr = VecCopy(*user.Up, *user.U);CHKERRQ(ierr);
      ierr = VecStrideScatter(Y1_c,0,*user.U,ADD_VALUES);CHKERRQ(ierr);
      ierr = VecStrideScatter(Y1_eta,1,*user.U,ADD_VALUES);CHKERRQ(ierr);
      ierr = VecStrideScatter(Y1_c,2,*user.U,ADD_VALUES);CHKERRQ(ierr);
      ierr = VecStrideScatter(Y1_eta,3,*user.U,ADD_VALUES);CHKERRQ(ierr);

      // Check that c, eta did not go out of the physical space...
      // i.e. 1/sqrt{2} - |c - 1/sqrt{2}| - |eta| > 0
      /*
      ierr = VecStrideGather(*(user.U),0,Y1_c,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecStrideGather(*(user.U),1,Y1_eta,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecAbs(Y1_eta);CHKERRQ(ierr);
      ierr = VecShift(Y1_c,-1./std::sqrt(2.));CHKERRQ(ierr);
      ierr = VecAbs(Y1_c);CHKERRQ(ierr);
      ierr = VecAXPBY(Y1_c,-1.,-1.,Y1_eta);CHKERRQ(ierr);
      ierr = VecShift(Y1_c,1./std::sqrt(2.));CHKERRQ(ierr);
      */
      double min = 1.;
      int i_min;
      //ierr = VecMin(Y1_c,&i_min,&min);CHKERRQ(ierr);
      if (min > 0){

	//If converged fast enough, move on
	if(n_iter < 4){//3){
	  counter++;
	  //If converged fast, keep track of it
	  if(counter >= 5){
          //If converged really fast for multiple iterations, scale up time step
	    PetscPrintf(PETSC_COMM_WORLD,"Increasing time step...\n");
	    user.dt *= std::sqrt(2.);
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
	PetscPrintf(PETSC_COMM_WORLD,"Went outside of physical range. Decreasing time step and trying again...\n");
	user.dt *= 0.5;
	counter = 0.;
	
	//Reset the current to the previous converged U                                                                                                                        
	ierr = VecCopy(*user.Up, *user.U);CHKERRQ(ierr);
      }
    }
    else{
      //If it didn't converge fast enough, scale back time step and try again.                                                                                               
      PetscPrintf(PETSC_COMM_WORLD,"Didn't converge, code %i. Decreasing time step and trying again...\n",conv_reason);
      user.dt *= 0.5;
      counter = 0.;

      //Reset the current to the previous converged U                                                                                                                        
      ierr = VecCopy(*user.Up, *user.U);CHKERRQ(ierr);
    }
    //clock_gettime(CLOCK_MONOTONIC, &finish);
    //elapsed = finish.tv_sec + 1.e-9*finish.tv_nsec - start.tv_sec - 1.e-9*start.tv_nsec;
    auto end = std::chrono::system_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double elapsed = 1.e-3*elapsed_ms;
    PetscPrintf(PETSC_COMM_WORLD,"Time elapsed: %f s\n",elapsed);

  }

  //finalize
  ierr = SNESDestroy(user.snes);CHKERRQ(ierr);
  ierr = VecDestroy(user.U);CHKERRQ(ierr);
  ierr = VecDestroy(user.Up);CHKERRQ(ierr);
  ierr = VecDestroy(user.Upp);CHKERRQ(ierr);
  ierr = VecDestroy(&Y1_c);CHKERRQ(ierr);
  ierr = VecDestroy(&Y1_eta);CHKERRQ(ierr);
  ierr = VecDestroy(&Y2_c);CHKERRQ(ierr);
  ierr = VecDestroy(&Y2_eta);CHKERRQ(ierr);
  ierr = IGADestroy(&user.iga);CHKERRQ(ierr);
  ierr = IGADestroy(&user.igaProject);CHKERRQ(ierr);

  //////////////////
  //ierr = Run<dim>();CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}
