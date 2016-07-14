#include "../applicationHeaders.h"
//extern "C" {
#include "petiga.h"
//}
//#include "../3D/parameters.h"
#include "../../../src/physicsHeaders.h"

template <int dim>
int timeStepSetup(AppCtx& user, TS& ts){
  PetscErrorCode  ierr;

  ierr = IGACreateTS(user.iga,&ts);CHKERRQ(ierr);
	ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,100001,1.0);CHKERRQ(ierr);
  ierr = TSSetTime(ts,user.RESTART_TIME);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,user.dt);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,OutputMonitor<dim>,&user,NULL);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,loadStep<dim>,&user,NULL);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,adaptiveTimeStep<dim>,&user,NULL);CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  return 0;

}
template int timeStepSetup<2>(AppCtx& user, TS& ts);
template int timeStepSetup<3>(AppCtx& user, TS& ts);

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

template PetscErrorCode loadStep<2>(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);
template PetscErrorCode loadStep<3>(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);

template <int dim>
PetscErrorCode adaptiveTimeStep(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;

  //adaptive TS
  double dt=user->dtVal; 
  ierr = TSSetTimeStep(*user->ts,dt);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: initial dt: %12.6e, dt: %12.6e \n",user->dtVal, dt);

  PetscFunctionReturn(0);
}

template int adaptiveTimeStep<2>(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);
template int adaptiveTimeStep<3>(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);
