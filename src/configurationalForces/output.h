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
  PetscReal UU[DIM], UUx[DIM][DIM];
  computeField<PetscReal,DIM,2*DIM>(VECTOR,0,p,U,&UU[0],&UUx[0][0]);
  //Compute \Chi
  PetscReal chi[DIM][DIM];
  for (PetscInt i=0; i<DIM; i++) {
    for (PetscInt j=0; j<DIM; j++) {
      chi[i][j]=(i==j)+UUx[i][j];
    }
  }

  //Compute strain metric, \Xi  (E=0.5*(F^T*F-I))
  PetscReal Phi[DIM][DIM], Xi[DIM][DIM];
  for (unsigned int I=0; I<DIM; I++){
    for (unsigned int J=0; J<DIM; J++){
      Phi[I][J] = 0.;
      for (unsigned int k=0; k<DIM; k++){
	Phi[I][J] += chi[k][I]*chi[k][J];
      }
      Xi[I][J] = 0.5*(Phi[I][J] - (I==J));
    }
  }

  //Compute the normal stretches \Lambda_I = \sqrt{\Phi_{II}} (no sum on I)
  PetscReal Lambda[DIM];
  for (unsigned int I=0; I<DIM; I++){
    Lambda[I] = sqrt(Phi[I][I]);
  }

	//Compute J_\chi (the determinant of \chi)
	PetscReal J_chi;
	J_chi = chi[0][0]*(chi[1][1]*chi[2][2] - chi[1][2]*chi[2][1]) -
					chi[0][1]*(chi[1][0]*chi[2][2] - chi[1][2]*chi[2][0]) +
					chi[0][2]*(chi[1][0]*chi[2][1] - chi[1][1]*chi[2][0]);

	//Compute \chi^{-1}
	PetscReal chi_Inv[DIM][DIM];
	chi_Inv[0][0] = 1./J_chi*(chi[1][1]*chi[2][2] - chi[2][1]*chi[1][2]);
	chi_Inv[0][1] = 1./J_chi*(chi[0][2]*chi[2][1] - chi[0][1]*chi[2][2]);
	chi_Inv[0][2] = 1./J_chi*(chi[0][1]*chi[1][2] - chi[1][1]*chi[0][2]);
	chi_Inv[1][0] = 1./J_chi*(chi[1][2]*chi[2][0] - chi[2][2]*chi[1][0]);
	chi_Inv[1][1] = 1./J_chi*(chi[0][0]*chi[2][2] - chi[2][0]*chi[0][2]);
	chi_Inv[1][2] = 1./J_chi*(chi[0][2]*chi[1][0] - chi[1][2]*chi[0][0]);
	chi_Inv[2][0] = 1./J_chi*(chi[1][0]*chi[2][1] - chi[2][0]*chi[1][1]);
	chi_Inv[2][1] = 1./J_chi*(chi[0][1]*chi[2][0] - chi[2][1]*chi[0][0]);
	chi_Inv[2][2] = 1./J_chi*(chi[0][0]*chi[1][1] - chi[1][0]*chi[0][1]);

  //displacement field variables
  PetscReal u[DIM], ux[DIM][DIM];
  computeField<PetscReal,DIM,2*DIM>(VECTOR,DIM,p,U,&u[0],&ux[0][0]);
  //Compute F (I+ux), dF (uxx)
  PetscReal F[DIM][DIM];
  for (unsigned int i=0; i<DIM; i++) {
    for (unsigned int J=0; J<DIM; J++) {
      F[i][J] = 0.;
      for (unsigned int k=0; k<DIM; k++){
      	F[i][J] += ((i==k) + UUx[i][k] + ux[i][k])*chi_Inv[k][J];
      }
    }
  }

  //Compute strain metric, E=0.5*(F^T*F-I)
  PetscReal E[DIM][DIM];
  for (unsigned int I=0; I<DIM; I++){
    for (unsigned int J=0; J<DIM; J++){
      E[I][J] = -0.5*(I==J);
      for (unsigned int k=0; k<DIM; k++){
	E[I][J] += 0.5*F[k][I]*F[k][J];
      }
    }
  }

  //define alpha and beta tensors
  PetscReal alpha[DIM], beta[DIM][DIM];
  for (unsigned int I=0; I<DIM; I++){
    alpha[I] = alphaC*Lambda[I];
    for (unsigned int J=0; J<DIM; J++){
//      beta[I][J] = betaC*Lambda[I]*Lambda[J];
      beta[I][J] = betaC; //Just use a constant beta for now.
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

  //compute P
  PetscReal P[DIM][DIM];
  for (unsigned int i=0; i<DIM; ++i){
    for (unsigned int J=0; J<DIM; ++J){ 
      P[i][J]=PiJ;
		}
	}


  AppCtx *user = (AppCtx *)ctx;
	if(user->first){
		FILE	*output_file = NULL;
 		PetscFOpen(PETSC_COMM_WORLD,"stress_stretch.txt","a+",&output_file);
		PetscFPrintf(PETSC_COMM_WORLD,output_file,"%g %g %g %g %g",F[0][0],P[0][0],Lambda[0],Lambda[1],Lambda[2]);
		PetscFClose(PETSC_COMM_WORLD,output_file);
		user->first = false;
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
      case 3:
	val=dist; break;
      case 4:
	val=P[0][0]; break;
      case 5:
	val=F[0][0]; break;
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
      if(c_time <1){
	user->lambda=c_time;
  }
  else{
    user->lambda=1.;
  }
  PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: load parameter: %6.2e\n",c_time);
 
  double dt=dtVal;
  double t;
  ierr = TSGetTime(*user->ts,&t);CHKERRQ(ierr);

 //output to file
  user->first = true;
  sprintf(filename,"./outU%d.dat",it_number);
  if (t<=1){
    ierr = IGAWriteVec(user->iga,U,filename);CHKERRQ(ierr);
    ProjectSolution(user->iga, it_number, U, user); 
  }
  if (t>1 && it_number%skipOutput==0){
    ierr = IGAWriteVec(user->iga,U,filename);CHKERRQ(ierr);
    ProjectSolution(user->iga, it_number, U, user); 
  }
//Temporarily reapply initial conditions
//  ierr = VecCopy(*user->U0, *user->U);CHKERRQ(ierr);
  
  //adaptive TS
  if(t<0.1){
    ierr = IGASetBoundaryValue(user->iga,1,1,1,uDirichlet*10.*t);CHKERRQ(ierr);
  }
  else if(t<0.3999) {dt*=10;}
   // else if(t>=1.) {dt*=1;}
   //else if(t>=1.05) {dt*=0.1;}
   //else if(t>=1.) {dt*=1;}
  else if(t>=0.999){dt*=1.;}
  else if(t>=0.6) {dt*=10;}
  // if (t<=.8) {dt*=200;}
  // else if (t<=1.01) {dt*=10;}
  // else if(t<=1.028) {dt *= 1;}
  // else if(t<=1.0296) {dt*=.1;}
  // else if(t<=1.02972){dt *= .01;}
  // else if(t<=1.039) {dt *= .001;}
   //else if(t<=1.079) {dt *= .1;}
   //else if(t<=1.19) {dt *= 1;}
  // else{dt*=10;}
  ierr = TSSetTimeStep(*user->ts,dt);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"USER SIGNAL: initial dt: %12.6e, dt: %12.6e \n",dtVal, dt); 

  double dVal=0.;
  //*
  if(t>1){
 	
    dVal = (t-1)*.1;
  PetscPrintf(PETSC_COMM_WORLD,"t: %12.6e, dVal: %12.6e  \n",t,dVal); 
  /*  ierr = IGASetBoundaryValue(user->iga,0,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user->iga,0,1,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user->iga,1,0,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user->iga,1,1,1,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user->iga,2,0,2,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user->iga,2,1,2,0.0);CHKERRQ(ierr);
*/ 
  // ierr = IGASetBoundaryValue(user->iga,0,0,0,0.0);CHKERRQ(ierr);  
  //ierr = IGASetBoundaryValue(user->iga,1,0,1,0.0);CHKERRQ(ierr);
  //ierr = IGASetBoundaryValue(user->iga,2,0,2,0.0);CHKERRQ(ierr);  
  // ierr = IGASetBoundaryValue(user->iga,1,1,1,1.e-4);CHKERRQ(ierr); 
  //ierr = IGASetBoundaryLoad(user->iga,1,1,1,1.);CHKERRQ(ierr);

  //ierr = IGASetBoundaryValue(user->iga,0,0,3,0.0);CHKERRQ(ierr);  
  //ierr = IGASetBoundaryValue(user->iga,1,0,4,0.0);CHKERRQ(ierr);
  //ierr = IGASetBoundaryValue(user->iga,2,0,5,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user->iga,0,1,3,dVal);CHKERRQ(ierr); 
  }

	FILE	*output_file = NULL;
 	PetscFOpen(PETSC_COMM_WORLD,"stress_stretch.txt","a+",&output_file);
	PetscFPrintf(PETSC_COMM_WORLD,output_file,"\n%g	",dVal);
	PetscFClose(PETSC_COMM_WORLD,output_file);
// */
  PetscFunctionReturn(0);
}

#endif