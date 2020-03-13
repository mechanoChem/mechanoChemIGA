#ifndef MECHANOCHEMIGA_H
#define MECHANOCHEMIGA_H

#include <math.h> 
#include <vector>
#include <string>
#include <memory>
//extern "C" {
#include "petiga.h"
//}
//Core functions
#include "coreFunctions.h"
//application specific functions
#include "userFunctions.h"
#include <iostream>

/***************************************************************
 * mechanoChemIGA: a class used in connection with a simple Python
 * interface to allow updating parameters and running from a Python script.
 */
template<unsigned int dim>
class mechanoChemIGA{
 public:

  /**
   * Constructor.
   */
  mechanoChemIGA(int argc, char **argv);

  /**
   * Destructor.
   */
  ~mechanoChemIGA();

  /**
   * Apply initial conditions.
   */
  void apply_initial_conditions();

  /**
   * Apply boundary conditions.
   */
  void apply_boundary_conditions();

  /**
   * Load parameters from parameters file.
   */
  void load_parameters();

  /**
   * Solve using settings in parameters.
   */
  void solve_ibvp();

  /**
   * Output results.
   */
  void output_results();

  /**
   * Get solution.
   */
  std::vector<double> get_solution(){ std::vector<double> a(0,0); return a;};

 private:

  AppCtx<dim> user;
  PetscInt step;
  Vec U,Up,Upp;
  SNES snes;

};

template<unsigned int dim>
mechanoChemIGA<dim>::mechanoChemIGA(int argc, char **argv)
{
  PetscInitialize(&argc,&argv,0,0);

  //Setup structures, initial conditions, boundary conditions
  Setup<dim>(user,&U,&Up,&Upp,snes);
  user.time = user.RESTART_TIME;
  step = user.RESTART_IT;
}

template<unsigned int dim>
mechanoChemIGA<dim>::~mechanoChemIGA()
{
  SNESDestroy(user.snes);
  VecDestroy(user.Upp);
  VecDestroy(user.Up);
  VecDestroy(user.U);
  IGADestroy(&user.iga);
  if (user.scalarProjectnFields.size() + user.vectorProjectnFields.size() > 0){
    IGADestroy(&user.igaProject);
  }
  PetscFinalize();
}

template<unsigned int dim>
void mechanoChemIGA<dim>::apply_initial_conditions()
{
  PetscPrintf(PETSC_COMM_WORLD,"applying ics...\n");

  user.time = user.RESTART_TIME;
  step = user.RESTART_IT;

  if(user.RESTART_IT==0){
    FormInitialCondition(user.iga, *user.Up, &user);
  }
  else{
    char filename[256];
    sprintf(filename,"%s/outU%d.dat",user.outputDir.c_str(),user.RESTART_IT);  
    IGAReadVec(user.iga,*user.Up,filename); //Read in vector to restart at step RESTART_IT
  }
  VecCopy(*user.Up, *user.U);
  VecCopy(*user.Up, *user.Upp); //For now, make U=Up=Upp for initial conditions, but shouldn't need to be this way.
}

template<unsigned int dim>
void mechanoChemIGA<dim>::apply_boundary_conditions()
{
  PetscPrintf(PETSC_COMM_WORLD,"applying bcs...\n");

  for (unsigned int i=0; i<dim; ++i){
    IGAFormClearBoundary(user.iga->form,i,0);
    IGAFormClearBoundary(user.iga->form,i,1);
    IGAFormSetBoundaryForm(user.iga->form,i,0,PETSC_TRUE);
    IGAFormSetBoundaryForm(user.iga->form,i,1,PETSC_TRUE);
  }
  user.boundaryConditions(user,0.);
}

template<unsigned int dim>
void mechanoChemIGA<dim>::load_parameters()
{
  ReadParameters<dim>(user);
  user.dt = user.dtVal;
}

template<unsigned int dim>
void mechanoChemIGA<dim>::solve_ibvp()
{
  unsigned int counter = 0;
  PetscInt n_iter;
  SNESConvergedReason conv_reason;

  PetscPrintf(PETSC_COMM_WORLD,"running...\n");
  StepUpdate<dim>(step,user.time,*(user.U),user);

  int step0 = step;
  while (user.time < user.totalTime && (user.maxTimeSteps == -1 || (step - step0) < user.maxTimeSteps)){
    PetscPrintf(PETSC_COMM_WORLD,"Step %i, dt %g, time %g\n",step,user.dt,user.time);
    SNESSolve(*(user.snes),NULL,*(user.U));
    SNESGetIterationNumber(*(user.snes),&n_iter);
    SNESGetConvergedReason(*(user.snes),&conv_reason);
    //PetscPrintf(PETSC_COMM_WORLD,"Convergence reason %i\n",conv_reason);
    //Start putting in some infrastructure for adaptive time stepping...
    if(n_iter < 10 && conv_reason > 0){
      //If converged fast enough, move on
      if(n_iter < 4){
	counter++;
	//If converged fast, keep track of it
	if(counter >= 5){
	  //If converged really fast for multiple iterations, scale up time step
	  if(user.adapTS){
	    PetscPrintf(PETSC_COMM_WORLD,"Increasing time step...\n");
	    user.dt *= std::pow(2.,1./3.);
	  }
	  counter = 0;
	}
      }
      user.time += user.dt;
      StepUpdate<dim>(++step,user.time,*(user.U),user);
      VecCopy(*user.Up, *user.Upp);
      VecCopy(*user.U, *user.Up);
    }
    else{
      //If it didn't converge fast enough, scale back time step and try again.
      if(user.adapTS){
	PetscPrintf(PETSC_COMM_WORLD,"Halving time step and trying again...\n");
	user.dt *= 0.5;
      }
      counter = 0;
      //Reset the current to the previous converged U
      VecCopy(*user.Up, *user.U);
    }
  }
}

template<unsigned int dim>
void mechanoChemIGA<dim>::output_results()
{
  char filename[256];
  sprintf(filename,"%s/outU%d.dat",user.outputDir.c_str(),step);
  IGAWriteVec(user.iga,*user.U,filename);

  if(user.scalarProjectnFields.size()+user.vectorProjectnFields.size() > 0){
    ProjectSolution<dim>(user.igaProject, step, *user.U, &user); 
  } 

  VecView(U,PETSC_VIEWER_STDOUT_WORLD);

  double *array;
  int size;
  VecGetLocalSize(U,&size);
  VecGetArray(U,&array);
  std::cout << "\nNext way:\n";
  for (int i=0; i<size; ++i){
    std::cout << array[i] << std::endl;
  }

}

#endif //MECHANOCHEMIGA_H
