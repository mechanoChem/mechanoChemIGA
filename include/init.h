#ifndef init_
#define init_

int init(AppCtx& user, AppCtxKSP& userKSP, PetscInt N, PetscInt p){
  PetscErrorCode  ierr;

  //set discretization options
  PetscInt C=PETSC_DECIDE;
  PetscBool output = PETSC_TRUE; 
  PetscBool monitor = PETSC_TRUE; 
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","CahnHilliard2D Options","IGA");CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  if (C == PETSC_DECIDE) C = p-1;
 
  //
  if (p < 2 || C < 0) /* Problem requires a p>=2 C1 basis */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Problem requires minimum of p = 2");
  if (p <= C)         /* Check C < p */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Discretization inconsistent: polynomial order must be greater than degree of continuity");
 
  //
  ierr = IGACreate(PETSC_COMM_WORLD,&user.iga);CHKERRQ(ierr);
  ierr = IGASetDim(user.iga,DIM);CHKERRQ(ierr);
  ierr = IGASetDof(user.iga,DIM+1);CHKERRQ(ierr);

  IGAAxis axis0;
  ierr = IGAGetAxis(user.iga,0,&axis0);CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis0,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis0,N,0.0,1.0,C);CHKERRQ(ierr);

  IGAAxis axis1;
  ierr = IGAGetAxis(user.iga,1,&axis1);CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis1,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis1,N,0.0,1.0,C);CHKERRQ(ierr);
  ierr = IGASetFromOptions(user.iga);CHKERRQ(ierr);
  ierr = IGASetUp(user.iga);CHKERRQ(ierr);
  //restart
  char meshfilename[256];
  sprintf(meshfilename, "mesh.dat");
  PetscPrintf(PETSC_COMM_WORLD,"\nWriting mesh file: %s\n", meshfilename);
  ierr = IGAWrite(user.iga, meshfilename);CHKERRQ(ierr);
  
  //fields 
  ierr = IGACreateVec(user.iga,user.U);CHKERRQ(ierr);
  ierr = IGACreateVec(user.iga,user.U0);CHKERRQ(ierr);
  ierr = FormInitialCondition2D(user.iga, *user.U0, &user); 
  ierr = VecCopy(*user.U0, *user.U);CHKERRQ(ierr);
  
  //
  IGAForm form;
  ierr = IGAGetForm(user.iga,&form);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,1,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,1,PETSC_TRUE);CHKERRQ(ierr);

  //assign residual and jacobian functions
  ierr = IGASetFormIEFunction(user.iga,Residual,&user);CHKERRQ(ierr);
  ierr = IGASetFormIEJacobian(user.iga,Jacobian,&user);CHKERRQ(ierr);

  return 0;
}

#endif
