#include <math.h> 
//extern "C" {
#include "petiga.h"
//}
//core functions
#include "coreFunctions.h"
//application specific functions
#include "userFunctions.h"
#include "mechanoChemIGA.h"
#include <iostream>

int main(int argc, char *argv[]) {
 
  //Initialize objects and parameters
  const unsigned int dim = 2;
  mechanoChemIGA<dim> prob(argc,argv);

  prob.load_parameters();
  prob.apply_boundary_conditions();
  prob.apply_initial_conditions();
  prob.solve_ibvp();
  prob.output_results();
  /*
  prob.solve_ibvp();

  prob.load_parameters();
  prob.apply_boundary_conditions();
  prob.apply_initial_conditions();
  prob.solve_ibvp();
  */
  return 0;
}
