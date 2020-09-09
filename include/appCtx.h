#ifndef appCtx_h
#define appCtx_h

//extern "C" {
#include "petiga.h"
//}
#include "tensor.h"
#include "solutionClass.h"
//#include "DNN.h"
#include <map>
#include <string>
#include <vector>
#include <random>
#include <stdio.h>

//Forward declaration of structure
template<unsigned int DIM>
struct AppCtx;

/**
 * Function pointer definition for user boundary conditions
 */
template<unsigned int dim>
using BCFunction = void(*)(AppCtx<dim>&, double);

/**
 * Function pointer definition for user initial conditions for scalar fields
 */
template<unsigned int dim>
using scalarICFunction = double(*)(const Tensor<1,dim,double>&, unsigned int, AppCtx<dim>&);

/**
 * Function pointer definition for user initial conditions for vector fields
 */
template<unsigned int dim>
using vectorICFunction = Tensor<1,dim,double>(*)(const Tensor<1,dim,double> &x, unsigned int vector_i, AppCtx<dim> &user);

/**
 * Function pointer definition for user load stepping definition
 */
template<unsigned int dim>
using LSFunction = void(*)(unsigned int step,AppCtx<dim> &user);

/**
 * Function pointer definition for user adaptive time stepping definition
 */
template<unsigned int dim>
using adaptiveTSFunction = void(*)(unsigned int step,AppCtx<dim> &user);

/**
 * Function pointer definition called immmediately after the parameters file is read
 */
template<unsigned int dim>
using postParams = void(*)(AppCtx<dim> &user);

/**
 * Function pointer definition for user definition of fields to project
 */
template<unsigned int dim>
using PFFunction = void(*)(const Tensor<1,dim,double> &x,
			   const Tensor<1,dim,double> &normal,
			   const solutionScalars<dim,double> &c,
			   const solutionVectors<dim,double> &u,
			   AppCtx<dim> &user,
			   std::vector<double> &scalarProjections,
			   std::vector<Tensor<1,dim,double> > &vectorProjections);

/**
 * Structure to hold parameters and variables used by PetIGA, the model, the IBVP, etc. Values used universally are explicitly
 * listed. Otherwise, they should be stored in the std::map matParam.
 */

template<unsigned int DIM>
struct AppCtx{

  /**
   * An IGA object from the PetIGA libary (bitbucket.org/dalcinl/petiga/), used in solving for the solution fields.
   */
  IGA iga;

  /**
   * An IGA object used in solving for the solution fields.
   */
  IGA igaProject;

  /**
   * A pointer to a SNES object from the PETSc library (www.mcs.anl.gov/petsc/), used for nonlinear solves.
   */
  SNES* snes;

  /**
   * A pointer to a PETSc vector that stores the solution vector for the current step.
   */
  Vec *U;

  /**
   * A pointer to a PETSc vector that stores the solution vector for the previous step.
   */
  Vec *Up;

  /**
   * A pointer to a PETSc vector that stores the solution vector for the step before the previous step.
   */
  Vec *Upp;

  /**
   * Value used to scale the domain. Defaults to 1.
   */
  PetscReal GridScale;

  /**
   * Flag specifying the use of the Sacado automatic differentiation package from Trilinos (trilinos.org/packages/sacado/). Defaults to true.
   */
  PetscBool ADSacado;

  /**
   * Value available for use in the definition of Dirchlet boundary conditions.
   */
  PetscReal uDirichlet;

  /**
   * 1st order tensor with components defining the number of elements in each direction.
   */
  Tensor<1,DIM,PetscInt> N;

  /**
   * 1st order tensor with components defining the dimension of the hyperrectangle in each direction.
   */
  Tensor<1,DIM,PetscReal> L;

  /**
   * 1st order tensor with components defining the periodicity of the domain in each direction. Defaults to PETSC_FALSE.
   */
  Tensor<1,DIM,PetscBool> periodic;

  /**
   * Initial time step value. Can be used or modified in the adaptiveTimeStep user function. Defaults to 1.
   */
  PetscReal dtVal;

  /**
   * Current time step value. Can be used or modified in the adaptiveTimeStep user function.
   */
  PetscReal dt;

  /**
   * Value defining the current time in the simulation (should only be read, not defined by the user).
   */
  PetscReal time;

  /**
   * Value defining the total time to be simulated. Defaults to 1.
   */
  PetscReal totalTime;

  /**
   * Integer defining how frequently to output results. Defaults to 1 (i.e. every time step).
   */
  PetscInt skipOutput;

  /**
   * Define the initial iteration number. Defaults to 0.
   */
  PetscInt RESTART_IT;

  /**
   * Define the initial time in the simulation. Defaults to 0.
   */
  PetscReal RESTART_TIME;

  /**
   * The output directory. The code will not create the directory if it does not exist. Defaults to "." (the current directory).
   */
  std::string outputDir;

  /**
   * Vector of strings defining the names of the scalar solution fields. The vector length is used to set the number of scalar solution fields.
   */
  std::vector<std::string> scalarSolnFields;

  /**
   * Vector of strings defining the names of the vector solution fields. The vector length is used to set the number of vector solution fields.
   */
  std::vector<std::string> vectorSolnFields;

  /**
   * Vector of strings defining the names of the scalar projection fields. The vector length is used to set the number of scalar projection fields.
   */
  std::vector<std::string> scalarProjectnFields;

  /**
   * Vector of strings defining the names of the vector projection fields. The vector length is used to set the number of vector projection fields.
   */
  std::vector<std::string> vectorProjectnFields;

  /**
   * Define the polynomial order of the basis functions. Defaults to 2.
   */
  PetscInt polyOrder;

  /**
   * Define the global continuity of the spline basis functions. This value must be less than the polynomial order. Defaults to 1.
   */
  PetscInt globalContinuity;

  /**
   * 4th order elasticity tensor.
   */
  Tensor<4,DIM,double> C_e;

  /**
   * C++ map available for including additional double variables.
   */
  std::map<std::string,PetscReal> matParam;

  /**
   * Void pointer available for additional user objects.
   */
  void *parameters;
  void *parameters2;

  double *par_mat;

  /**
   * Flag specifying whether to solve Cahn-Hilliard equations or not. Defaults to false
   */
  bool CahnHilliard;

  /**
   * Flag specifying whether to solve finite strain elasticity equations or not. Defaults to false
   */
  bool Elasticity;

  /**
   * Object for random number generation
   */
  std::mt19937 gen;

  /**
   * Pointer to function setting Dirichlet boundary conditions.
   */
  BCFunction<DIM> boundaryConditions;

  /**
   * Pointer to function setting initial conditions for scalar fields.
   */
  scalarICFunction<DIM> scalarInitialConditions;

  /**
   * Pointer to function setting initial conditions for vector fields.
   */
  vectorICFunction<DIM> vectorInitialConditions;

  /**
   * Pointer to function defining user load stepping.
   */
  LSFunction<DIM> loadStep;

  /**
   * Pointer to function defining user adaptive time stepping.
   */
  adaptiveTSFunction<DIM> adaptiveTimeStep;

  /**
   * Pointer to function called immediately after reading parameters file.
   */
  postParams<DIM> postParameters;

  /**
   * Pointer to function defining fields to project.
   */
  PFFunction<DIM> projectFields;

  //Define default values
  AppCtx<DIM>() : GridScale(1.), ADSacado(PETSC_TRUE), uDirichlet(0.), dtVal(1.), totalTime(1.), skipOutput(1), RESTART_IT(0), RESTART_TIME(0.), outputDir("."), polyOrder(2), globalContinuity(1), CahnHilliard(false), Elasticity(false) {}
};

#endif
