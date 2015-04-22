#ifndef output_h
#define output_h
#include <cmath>

#define PI 3.14159265

#undef  __FUNCT__
#define __FUNCT__ "E22Function"
PetscErrorCode E22Function(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *ctx)
{	
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);

  //displacement field variables
  PetscReal u[DIM], ux[DIM][DIM];
  computeField<PetscReal,DIM,2*DIM>(VECTOR,0,p,U,&u[0],&ux[0][0]);
  //Compute F
  PetscReal chi[DIM][DIM];
  for (PetscInt i=0; i<DIM; i++) {
    for (PetscInt j=0; j<DIM; j++) {
      chi[i][j]=(i==j)+ux[i][j];
    }
  }

  //Compute strain metric, E  (E=0.5*(F^T*F-I))
  PetscReal Xi[DIM][DIM];
  for (unsigned int I=0; I<DIM; I++){
    for (unsigned int J=0; J<DIM; J++){
      Xi[I][J] = -0.5*(I==J);
      for (unsigned int k=0; k<DIM; k++){
	Xi[I][J] += 0.5*chi[k][I]*chi[k][J];
      }
    }
  }

#if DIM==2
  //new strain metrics
  PetscReal e1=(Xi[0][0]+Xi[1][1]);
  PetscReal e2=(Xi[0][0]-Xi[1][1]);
  PetscReal e6=Xi[0][1];
  
  //compute distance to nearest well
  PetscReal dist=e2-Es;
  unsigned int wellID=1;

#elif DIM==3
  //new strain metrics
  PetscReal e1=(Xi[0][0]+Xi[1][1]+Xi[2][2])/sqrt(3.0);
  PetscReal e2=(Xi[0][0]-Xi[1][1])/sqrt(2.0);
  PetscReal e3=(Xi[0][0]+Xi[1][1]-2*Xi[2][2])/sqrt(6.0);
  PetscReal e4=Xi[1][2], e5=Xi[2][0], e6=Xi[0][1];
  //compute distance to nearest well
  PetscReal x[3],y[3]; 
  x[0]=0; y[0]=-Es; //first well 
  x[1]=Es*cos(30.0*PI/180.0); y[1]=Es*sin(30.0*PI/180.0); //second well
  x[2]=-Es*cos(30.0*PI/180.0); y[2]=Es*sin(30.0*PI/180.0); //third well
  PetscReal dist=sqrt(std::pow(e2-x[0],2.0)+std::pow(e3-y[0],2.0));
  unsigned int wellID=1; 
  for(unsigned int i=1; i<3; i++){
    if(dist>sqrt(pow(e2-x[i],2.0)+pow(e3-y[i],2.0))){
      dist=sqrt(pow(e2-x[i],2.0)+pow(e3-y[i],2.0));
      wellID=i+1;
    }
  }
#endif
 
  //store L2 projection residual
  const PetscReal (*N) = (PetscReal (*)) p->shape[0];;  
  for(int n1=0; n1<nen; n1++){
    for(int d1=0; d1<dof; d1++){
      PetscReal val=0.0;
      switch (d1) {
#if DIM==2
      case 0:
	val=e2; break;
      case 1:
	val=e6; break;
#elif DIM==3
      case 0:
	val=e2; break;
      case 1:
	val=e3; break;
      case 2:
	val=wellID; break;
#endif
      } 
      R[n1*dof+d1] = N[n1]*val;
    }
  }
  return 0;
}

#undef  __FUNCT__
#define __FUNCT__ "E22Jacobian"
PetscErrorCode E22Jacobian(IGAPoint p, const PetscScalar *U, PetscScalar *K, void *ctx)
{	
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  //
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

#undef  __FUNCT__
#define __FUNCT__ "ProjectSolution"
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

#undef  __FUNCT__
#define __FUNCT__ "OutputMonitor"
template <int dim>
PetscErrorCode OutputMonitor(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;
  char           filename[256];
  //setting load parameter
  user->lambda=c_time; 
  PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: load parameter: %6.2e\n",c_time);

  //output to file
  sprintf(filename,"./outU%d.dat",it_number);
  if (it_number%skipOutput==0){
    ierr = IGAWriteVec(user->iga,U,filename);CHKERRQ(ierr);
    ProjectSolution(user->iga, it_number, U, user); 
  }
  
  //adaptive TS
  double dt=dtVal;
  double t;
  ierr = TSGetTime(*user->ts,&t);CHKERRQ(ierr);
  if (t<0.4) {dt*=10;}
  ierr = TSSetTimeStep(*user->ts,dt);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: initial dt: %12.6e, dt: %12.6e \n",dtVal, dt); 

  // 
  PetscFunctionReturn(0);
}

#endif
