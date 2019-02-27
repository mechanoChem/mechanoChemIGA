#include <math.h> 
//extern "C" {
#include "petiga.h"
//}
//core functions
#include "coreFunctions.h"
//application specific functions
#include "userFunctions.h"

#include <time.h>

int main(int argc, char *argv[]) {
 
  //struct timespec start, finish;
  //clock_gettime(CLOCK_MONOTONIC, &start);

  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);
 
  //Initialize objects and parameters
  const unsigned int dim = 2;
  ierr = Run<dim>();CHKERRQ(ierr);

  //clock_gettime(CLOCK_MONOTONIC, &finish);
  //double elapsed = finish.tv_sec + 1.e-9*finish.tv_nsec - start.tv_sec - 1.e-9*start.tv_nsec;

  //PetscPrintf(PETSC_COMM_WORLD,"Time elapsed: %f s\n",elapsed);

  ierr = PetscFinalize();CHKERRQ(ierr);

return 0;
}
