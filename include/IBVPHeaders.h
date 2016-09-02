#ifndef applicationHeaders_
#define applicationHeaders_

//extern "C" {
#include "petiga.h"
//}

#include "appCtx.h"

/**
 * Define and form the initial conditions of the IBVP.
*/

PetscErrorCode FormInitialCondition(IGA iga, Vec U, AppCtx *user);

/**
 * Setup the data structures, PetIGA structures, PETSc structures, etc. Called once at the beginning of the code.
*/

template<unsigned int DIM,unsigned int DOF>
int setup(AppCtx& user,Vec *U,Vec *U0,TS &ts);

/**
 * Define the Dirichlet boundary conditions. Potentially called multiple times.
*/

int boundaryConditions(AppCtx& user, double scale);

/**
 * Define any load stepping (i.e. change in Dirichlet b.c.s). Called at each time/load step.
 */

PetscErrorCode loadStep(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);

/**
 * Specify any adaptive time step schemes. Called at each time/load step.
*/

PetscErrorCode adaptiveTimeStep(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);

/**
 * Specify the nonliner solve convergence test.
*/

PetscErrorCode SNESConvergedTest(SNES snes, PetscInt it,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx);

/**
 * Function to define function parameters/variables. Called once at the beginning of the code.
*/

int defineParameters(AppCtx& user);

#endif
