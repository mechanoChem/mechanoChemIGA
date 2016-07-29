//extern "C" {
#include "petiga.h"
//}
#include "IBVPHeaders.h"
#include "physicsHeaders.h"

PetscErrorCode loadStep(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;
  //setting load parameter
  user->lambda=1.;
  //user->lambda=c_time;

  PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: load parameter: %6.2e\n",user->lambda);

  double t;
  ierr = TSGetTime(*user->ts,&t);CHKERRQ(ierr);

	//Displacement loading
	double scale = t;
  PetscPrintf(PETSC_COMM_WORLD,"t: %12.6e",t); 
	boundaryConditions(*user,scale);

  PetscFunctionReturn(0);
}

PetscErrorCode adaptiveTimeStep(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;

  //adaptive TS
  double t;
  ierr = TSGetTime(*user->ts,&t);CHKERRQ(ierr);

  double dt=user->dtVal; 

  if (t<0.4) {dt*=10;}

  ierr = TSSetTimeStep(*user->ts,dt);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: initial dt: %12.6e, dt: %12.6e \n",user->dtVal, dt);

  PetscFunctionReturn(0);
}
