#ifndef init_
#define init_

int init(IGA& iga, AppCtx& user, PetscInt N, PetscInt p, PetscInt resStep, Vec&U, Vec&U0){
  PetscErrorCode  ierr;
 
  //set discretization options
  PetscInt nsteps = 100;
  PetscInt C=PETSC_DECIDE;
  PetscBool output = PETSC_TRUE; 
  PetscBool monitor = PETSC_TRUE; 
  char filePrefix[PETSC_MAX_PATH_LEN] = {0};
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","CahnHilliard2D Options","IGA");CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  if (C == PETSC_DECIDE) C = p-1;
 
  //
  if (p < 2 || C < 0) /* Problem requires a p>=2 C1 basis */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Problem requires minimum of p = 2");
  if (p <= C)         /* Check C < p */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Discretization inconsistent: polynomial order must be greater than degree of continuity");
 
  //
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDim(iga,DIM);CHKERRQ(ierr);
  ierr = IGASetDof(iga,DIM+1);CHKERRQ(ierr);

  IGAAxis axis0;
  ierr = IGAGetAxis(iga,0,&axis0);CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis0,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis0,N,0.0,1.0,C);CHKERRQ(ierr);

  IGAAxis axis1;
  ierr = IGAGetAxis(iga,1,&axis1);CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis1,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis1,N,0.0,1.0,C);CHKERRQ(ierr);
  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  //restart
  if (resStep==0){
    char meshfilename[256];
    sprintf(meshfilename, "mesh.dat");
    PetscPrintf(PETSC_COMM_WORLD,"\nWriting mesh file: %s\n", meshfilename);
    ierr = IGAWrite(iga, meshfilename);CHKERRQ(ierr);
  }
  
  //fields 
  ierr = IGACreateVec(iga,&U);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&U0);CHKERRQ(ierr);
  if (resStep>0){
    MPI_Comm comm;
    PetscViewer viewer;
    char restartfilename[256];
    ierr = PetscObjectGetComm((PetscObject)U0,&comm);CHKERRQ(ierr);
    sprintf(restartfilename,"res%s-%d.dat",filePrefix,resStep);
    ierr = PetscViewerBinaryOpen(comm,restartfilename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecLoad(U0,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);
    PetscPrintf(PETSC_COMM_WORLD,"\nReading solution from restart file: %s\n",restartfilename);
  }
  else{
    ierr = FormInitialCondition2D(iga, U0, &user); 
  }
  ierr = VecCopy(U0, U);CHKERRQ(ierr);
  
  //
  IGAForm form;
  ierr = IGAGetForm(iga,&form);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,1,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,1,PETSC_TRUE);CHKERRQ(ierr);

  //assign residual and jacobian functions
  ierr = IGASetFormIEFunction(iga,Residual,&user);CHKERRQ(ierr);
  ierr = IGASetFormIEJacobian(iga,Jacobian,&user);CHKERRQ(ierr);

  //time stepping
  TS ts;
  ierr = IGACreateTS(iga,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,100000,1.0);CHKERRQ(ierr);
  ierr = TSSetTime(ts,0.0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,user.dt);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,OutputMonitor<DIM>,&user,NULL);CHKERRQ(ierr);  
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
#if PETSC_VERSION_LE(3,3,0)
  ierr = TSSolve(ts,U,NULL);CHKERRQ(ierr);
#else
  ierr = TSSolve(ts,U);CHKERRQ(ierr);
#endif
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  //
  return 0;
}

#endif
