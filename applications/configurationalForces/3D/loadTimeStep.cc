#include "applicationHeaders.h"
//extern "C" {
#include "petiga.h"
//}
#include "parameters.h"
#include "../../../src/configurationalForces/physicsHeaders.h"

int timeStepSetup(AppCtx& user, TS& ts){
  PetscErrorCode  ierr;

  ierr = IGACreateTS(user.iga,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,100001,1.0);CHKERRQ(ierr);
  ierr = TSSetTime(ts,RESTART_TIME);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,user.dt);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,OutputMonitor<DIM>,&user,NULL);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,loadStep<DIM>,&user,NULL);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,adaptiveTimeStep<DIM>,&user,NULL);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  return 0;

}

template <int dim>
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
	boundaryConditions<dim>(*user,scale);

  PetscFunctionReturn(0);
}

template int loadStep<2>(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);
template int loadStep<3>(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);

template <int dim>
PetscErrorCode adaptiveTimeStep(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;

  //adaptive TS
  double dt=dtVal; 
  ierr = TSSetTimeStep(*user->ts,dt);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: initial dt: %12.6e, dt: %12.6e \n",dtVal, dt);

  PetscFunctionReturn(0);
}

template int adaptiveTimeStep<2>(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);
template int adaptiveTimeStep<3>(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);
