#ifndef output_h
#define output_h
#include <cmath>

PetscErrorCode E22Function(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *ctx)
{	
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);

  //displacement field variables
  PetscReal u[DIM], ux[DIM][DIM];
  computeField<PetscReal,DIM,DIM+1>(VECTOR,0,p,U,&u[0],&ux[0][0]);
  PetscReal c0;
  computeField<PetscReal,DIM,DIM+1>(SCALAR,DIM,p,U,&c0);
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
	val=c0; break;
      case 2: //only in 3D 
	val=c0; break;
      } 
      R[n1*dof+d1] = N[n1]*val;
    }
  }
  return 0;
}

PetscErrorCode E22Jacobian(IGAPoint p, const PetscScalar *U, PetscScalar *K, void *ctx)
{	
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);


  const PetscReal (*N) = (PetscReal (*)) p->shape[0];;  
  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
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

PetscErrorCode ProjectSolution(IGA iga, PetscInt step, Vec U, AppCtx *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;
  //Setup linear system for L2 Projection
  Mat A;
  Vec x,b;
  ierr = IGACreateMat(iga,&A);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&x);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&b);CHKERRQ(ierr);
  ierr = IGASetFormFunction(iga,E22Function,user);CHKERRQ(ierr);
  ierr = IGASetFormJacobian(iga,E22Jacobian,user);CHKERRQ(ierr);
  ierr = IGAComputeFunction(iga,U,b);CHKERRQ(ierr);
  ierr = IGAComputeJacobian(iga,U,A);CHKERRQ(ierr);

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
    ProjectSolution(user->iga, it_number, U, user); 
  }
  
  //Check for min(C), and stop if min(C)>1
  PetscReal minC, maxC;
  VecSetBlockSize(U,dim+1);  
  VecStrideMin(U,2,NULL,&minC);
  VecStrideMax(U,2,NULL,&maxC);
  PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: min(C): %12.6e, max(C): %12.6e \n",minC, maxC);
  if (minC>1.0){
    PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: min(C) > 1.0 so forcefully quitting \n");
    exit(-1);
  }
  //adaptive TS
  double dt=dtVal;
  if (maxC<0.35) dt=dtVal*100;
  else if ((maxC>=0.35)&&(maxC<0.45)) dt=dtVal*10;
  else if (maxC>=0.45) dt=dtVal;  
  ierr = TSSetTimeStep(*user->ts,dt);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: initial dt: %12.6e, dt: %12.6e \n",dtVal, dt); 

  //check for no change in solution between every 100 timesteps and quit
  PetscReal normC;
  VecStrideNorm(U,2,NORM_2,&normC);
  if (it_number%100==0){
    if (it_number>0){
      PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: change in concentration over last 100 time steps: %12.6e, norm: %12.6e, oldnorm: %12.6e\n",std::abs(normC-user->norm), normC, user->norm); 
      if (std::abs(normC-user->norm)<1.0e-6){
	PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: change in concentration below tol (1e-8) over last 100 time steps, solve probable stalled hence quitting\n"); exit(-1);     
      }
    }
    user->norm=normC;
  } 
  //BC
  /*
  double dVal=0.001;
  if ((it_number>100) && (it_number<200)){
    dVal*=(it_number-100);
    PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: dVal: %12.6e \n",dVal); 
    ierr = IGASetBoundaryValue(user->iga,0,0,1,dVal);CHKERRQ(ierr);  
    ierr = IGASetBoundaryValue(user->iga,0,1,1,-dVal);CHKERRQ(ierr);  
    ierr = IGASetBoundaryValue(user->iga,1,0,0,dVal);CHKERRQ(ierr);  
    ierr = IGASetBoundaryValue(user->iga,1,1,0,-dVal);CHKERRQ(ierr);
  }
  */
  // 
  PetscFunctionReturn(0);
}

#endif
