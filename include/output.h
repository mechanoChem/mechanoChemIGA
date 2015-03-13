
PetscErrorCode E22System(IGAPoint p, PetscScalar *K, PetscScalar *R, void *ctx)
{	
  AppCtxKSP *user = (AppCtxKSP *)ctx;
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);

  //displacement field variables
  PetscReal u[DIM], ux[DIM][DIM];
  computeField<PetscReal,DIM,DIM+1>(VECTOR,0,p,user->localU0,&u[0],&ux[0][0]);
  //Compute F
  PetscReal F[DIM][DIM];
  for (PetscInt i=0; i<DIM; i++) {
    for (PetscInt j=0; j<DIM; j++) {
      F[i][j]=(i==j)+ux[i][j];
    }
  }

  //Compute strain metric, E  (E=0.5*(F^T*F-I))
  PetscReal E[DIM][DIM];
  for (unsigned int I=0; I<DIM; I++){
    for (unsigned int J=0; J<DIM; J++){
      E[I][J] = -0.5*(I==J);
      for (unsigned int k=0; k<DIM; k++){
	E[I][J] += 0.5*F[k][I]*F[k][J];
      }
    }
  }

  //new strain metrics
  PetscReal e1=(E[0][0]+E[1][1]);
  PetscReal e2=(E[0][0]-E[1][1]);
  PetscReal e6=E[0][1];
  
  //compute distance to nearest well
  PetscReal dist=e2-Es;
  unsigned int wellID=1;
 
  //store L2 projection residual
  const PetscReal (*N) = (PetscReal (*)) p->shape[0];;  
  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      PetscReal val=0.0;
      switch (d1) {
      case 0:
	val=e2; break;
      case 1:
	val=e6; break;
      case 3: //only in 3D with concentration
	val=dist; break;
      case 2: //only in 3D 
	val=(PetscReal) wellID; break;
      } 
      R[n1*dof+d1] = N[n1]*val;
      for(int n2=0; n2<nen; n2++){
	for(int d2=0; d2<dof; d2++){
	  PetscReal val2=0.0;
	  if (d1==d2) {val2 = N[n1] * N[n2];}
	  K[n1*dof*nen*dof + d1*nen*dof + n2*dof + d2] =val2;
	}
      }
    }
  }
  return 0;
}

PetscErrorCode ProjectSolution(IGA iga, PetscInt step, AppCtxKSP *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;
  //Setup linear system for L2 Projection
  Mat A;
  Vec x,b;
  ierr = IGACreateMat(iga,&A);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&x);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&b);CHKERRQ(ierr);
  ierr = IGASetFormSystem(iga,E22System,user);CHKERRQ(ierr);
  ierr = IGAComputeSystem2(iga,A,b);CHKERRQ(ierr);

  //Solver
  KSP ksp;
  ierr = IGACreateKSP(iga,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
  //write solution
  char filename[256];
  sprintf(filename,"./outE%d.dat",step);
  ierr = IGAWriteVec(iga,x,filename);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  PetscFunctionReturn(0); 
}

template <int dim>
PetscErrorCode OutputMonitor(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx)
{
  
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;
  char           filename[256];
  sprintf(filename,"./outU%d.dat",it_number);
  if (it_number%skipOutput==0){
    ierr = IGAWriteVec(user->iga,U,filename);CHKERRQ(ierr);
    ProjectSolution(user->iga, it_number, user->appCtxKSP);
  }
  
  //Check for min(C), and stop if min(C)>1
  PetscReal minC;
  VecSetBlockSize(U,dim+1);  
  VecStrideMin(U,2,NULL,&minC);
  PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: min(C): %12.6e \n",minC);
  if (minC>1.0){
    PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: min(C) > 1.0 so forcefully quitting \n");
    exit(-1);
  }
  
  PetscFunctionReturn(0);
}
