#ifndef timeStepping_
#define timeStepping_

int timeStepping(AppCtx& user, TS& ts){
  PetscErrorCode  ierr;

  ierr = IGACreateTS(user.iga,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,100001,1.0);CHKERRQ(ierr);
  ierr = TSSetTime(ts,RESTART_TIME);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,user.dt);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,OutputMonitor<DIM>,&user,NULL);CHKERRQ(ierr);  
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

	return 0;

}

#endif
