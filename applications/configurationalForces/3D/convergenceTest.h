#ifndef convergenceTest_
#define convergenceTest_

//snes convegence test
PetscErrorCode SNESConvergedTest(SNES snes, PetscInt it,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx){
  AppCtx *user  = (AppCtx*) ctx;
  PetscPrintf(PETSC_COMM_WORLD,"xnorm:%12.6e snorm:%12.6e fnorm:%12.6e\n",xnorm,snorm,fnorm);  
  //custom test
  if (NVal>=15){
    if (it>500){
      *reason = SNES_CONVERGED_ITS;
      return(0);
    }
  }
  //default test
  PetscFunctionReturn(SNESConvergedDefault(snes,it,xnorm,snorm,fnorm,reason,ctx));
}

PetscErrorCode SNESConvergedTest_Interactive(SNES snes, PetscInt it,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx)
{
  AppCtx *user  = (AppCtx*) ctx;
  PetscPrintf(PETSC_COMM_WORLD,"xnorm:%12.6e snorm:%12.6e fnorm:%12.6e\n",xnorm,snorm,fnorm);
  if ((it>10) && (fnorm<1.0e-7)) {
    PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: since it>10 forcefully setting convergence. \n");
    *reason = SNES_CONVERGED_FNORM_ABS;
    return(0);
  }
  PetscFunctionReturn(SNESConvergedDefault(snes,it,xnorm,snorm,fnorm,reason,ctx));
}

int setConvergenceTest(AppCtx& user, TS& ts){
  PetscErrorCode  ierr;
  SNES snes;

  ierr = TSGetSNES(ts,&snes); CHKERRQ(ierr);
  ierr = SNESSetConvergenceTest(snes,SNESConvergedTest_Interactive,(void*)&user,NULL); CHKERRQ(ierr);

	return 0;
}

#endif
