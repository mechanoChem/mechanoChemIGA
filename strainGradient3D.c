#include <math.h>
#include "petiga.h"
#include "derivatives2.h"

typedef struct{ 
  IGA iga;
  PetscReal dt;
  PetscReal ta, tb, tc;
} AppCtx;

PetscErrorCode computeField(IGAPoint p, const PetscReal* U, PetscReal(*F)[DIM], PetscReal(*dF)[DIM][DIM], PetscReal(*E)[DIM], PetscReal (*dE)[DIM][DIM]){
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  PetscReal (*U2)[DIM] = (PetscReal (*)[DIM])U;
  //F
  PetscReal (*N1)[DIM] = (PetscReal (*)[DIM]) p->shape[1];
  for(unsigned int i=0; i<DIM; i++)
    for(unsigned int J=0; J<DIM; J++){
      F[i][J]=(i==J);
      for(unsigned int n=0; n<nen; n++)
	F[i][J]+= N1[n][J]*U2[n][i];
    }
  //dF
  PetscReal (*N2)[DIM][DIM] = (PetscReal (*)[DIM][DIM]) p->shape[2];
  for(unsigned int i=0; i<DIM; i++)
    for(unsigned int J=0; J<DIM; J++)
      for(unsigned int K=0; K<DIM; K++){
	dF[i][J][K]=0.0;
	for(unsigned int n=0; n<nen; n++)
	  dF[i][J][K]+= N2[n][J][K]*U2[n][i];
      }

  //E
  for (unsigned int I=0; I<DIM; I++)
    for (unsigned int J=0; J<DIM; J++){
      E[I][J] = -0.5*(I==J);
      for (unsigned int k=0; k<DIM; k++)
	E[I][J] += 0.5*F[k][I]*F[k][J];
    }
  
  //dE
  for (unsigned int I=0; I<DIM; I++)
    for (unsigned int J=0; J<DIM; J++)
      for (unsigned int K=0; K<DIM; K++){
	dE[I][J][K]=0.0;
	for (unsigned int i=0; i<DIM; i++)
	  dE[I][J][K] += 0.5*(dF[i][I][K]*F[i][J]+dF[i][J][K]*F[i][I]);
    }
  return 0;
}

#undef  __FUNCT__
#define __FUNCT__ "Residual"
PetscErrorCode Residual(IGAPoint p,PetscReal dt,
                        PetscReal shift,const PetscScalar *V,
                        PetscReal t,const PetscScalar *U,
                        PetscReal t0,const PetscScalar *U0, 
			PetscScalar *R,void *ctx){
  AppCtx *user = (AppCtx *)ctx;
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  PetscReal *n = p->normal;

  //Compute F, dF
  PetscReal F[DIM][DIM], dF[DIM][DIM][DIM];
  PetscReal E[DIM][DIM], dE[DIM][DIM][DIM];
  computeField(p, U, F, dF, E, dE);
 
  //P, Beta
  PetscReal P[DIM][DIM]; 
  PetscReal B[DIM][DIM][DIM];
  residual(F, dF, E, dE, P, B);
  
  /* //get shape function values */
  double (*N) = (double (*)) p->shape[0];
  double (*Nx)[DIM] = (double (*)[DIM]) p->shape[1];
  double (*Nxx)[DIM][DIM] = (double (*)[DIM][DIM]) p->shape[2];
  
  //Compute Residual
  int surfaceFlag=p->atboundary;
  PetscReal (*Ra)[DIM] = (PetscReal (*)[DIM])R;
  for (unsigned int a=0; a<nen; a++) {
    //Mechanics
    if (!surfaceFlag) {
      for (unsigned int i=0; i<DIM; i++){
	PetscReal Ru_i=0.0;
	for (unsigned int j=0; j<DIM; j++){
	  //grad(Na)*P
	  Ru_i += Nx[a][j]*P[i][j];
	  for (unsigned int k=0; k<DIM; k++){
	    Ru_i += Nxx[a][j][k]*B[i][j][k];
	  }
	}
	Ra[a][i] = Ru_i;
      }
    }
  }
  return 0;
}

#undef  __FUNCT__
#define __FUNCT__ "Jacobian"
PetscErrorCode Jacobian(IGAPoint p,PetscReal dt,
			PetscReal shift,const PetscScalar *V,
			PetscReal t,const PetscScalar *U,
			PetscReal t0,const PetscScalar *U0,
			PetscScalar *Kab,void *ctx){
  AppCtx *user = (AppCtx *)ctx;
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  PetscReal *n = p->normal;

  //Compute F, dF
  PetscReal F[DIM][DIM], dF[DIM][DIM][DIM];
  PetscReal E[DIM][DIM], dE[DIM][DIM][DIM];
  computeField(p, U, F, dF, E, dE);
 
  //P, Beta
  PetscReal PbyF[DIM][DIM][DIM][DIM]; 
  PetscReal PbydF[DIM][DIM][DIM][DIM][DIM];
  PetscReal BbyF[DIM][DIM][DIM][DIM][DIM];
  PetscReal BbydF[DIM][DIM][DIM][DIM][DIM]; 
  double ta = MPI_Wtime();
  jacobian(F,dF,E,dE,PbyF,PbydF,BbyF,BbydF);
  user->ta+=MPI_Wtime()-ta;

  /* //get shape function values */
  double (*N) = (double (*)) p->shape[0];
  double (*Nx)[DIM] = (double (*)[DIM]) p->shape[1];
  double (*Nxx)[DIM][DIM] = (double (*)[DIM][DIM]) p->shape[2];
  
  //Compute Residual
  PetscReal PbyFJMil=0.0, PbydFJMNil=0.0, BbyFJKMil=0.0, BbydFJKMNil=0.0;
  double tb = MPI_Wtime();
  for (unsigned int a=0; a<nen; a++) {
    for (unsigned int b=0; b<nen; b++) {
      for (unsigned int i=0; i<DIM; i++){
	for (unsigned int l=0; l<DIM; l++){
	  PetscReal Kabil=0.0;
	  //PbyF
	  for (unsigned int J=0; J<DIM; J++){
	    for (unsigned int M=0; M<DIM; M++){
	      //loop over i,l,a,b
	      PbyFJMil=PbyF[i][l][J][M];
	      Kabil+=Nx[a][J]*PbyFJMil*Nx[b][M]; //PbyF
	      //PbydF
	      for (unsigned int N=0; N<DIM; N++){
		//loop over i,l,a,b
		PbydFJMNil=PbydF[i][l][J][M][N];
		Kabil+=1.0e-15*PbydFJMNil;//Nx[a][J]*PbydFJMNil*Nxx[b][M][N]; //PbydF
	      }
	    }
	    //BbyF
	    for (unsigned int K=0; K<DIM; K++){
	      for (unsigned int M=0; M<DIM; M++){
		//loop over i,l,a,b
		BbyFJKMil=BbyF[i][l][J][K][M];
		Kabil+=1.0e-15*BbyFJKMil;//Nxx[a][J][K]*BbyFJKMil*Nx[b][M]; //BbyF
		//BbydF
		//double tc= MPI_Wtime();
		unsigned int N=K;
		//for (unsigned int N=0; N<DIM; N++){	   
		//loop over i,l,a,b
		BbydFJKMNil=BbydF[i][l][J][K][M];
		Kabil+=1.0e-15*BbydFJKMNil; //Nxx[a][J][K]*BbydFJKMNil*Nxx[b][M][N]; //BbydF
	      }
	    }
	    }
	  //user->tc+=MPI_Wtime()-tc;
	  Kab[a*DIM*nen*DIM + i*nen*DIM + b*DIM + l]=Kabil;
	}
      }
    }
  }
  user->tb+=MPI_Wtime()-tb;
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "OutputMonitor"
PetscErrorCode OutputMonitor(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;
  char           filename[256];
  sprintf(filename,"./outU%d.dat",it_number);
  ierr = IGAWriteVec(user->iga,U,filename);CHKERRQ(ierr);
  //ProjectSolution(user->iga, it_number, user->appCtxKSP);
  PetscFunctionReturn(0);
}

/*
typedef struct {
  PetscReal ux, uy, uz;
} Field;

PetscErrorCode FormInitialCondition(IGA iga, Vec U, AppCtx *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;
  DM da;
  ierr = IGACreateNodeDM(iga,DIM,&da);CHKERRQ(ierr);
  Field ***u;
  ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);

  PetscInt i,j,k;
  for(i=info.xs;i<info.xs+info.xm;i++){
    for(j=info.ys;j<info.ys+info.ym;j++){
      for(k=info.zs;k<info.zs+info.zm;k++){
	u[k][j][i].ux=0.0;
	u[k][j][i].uy=0.0;
	u[k][j][i].uz=0.0;
      }
    }
  }
  ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
  ierr = DMDestroy(&da);;CHKERRQ(ierr); 
  PetscFunctionReturn(0); 
}
*/

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);
  double startTime = MPI_Wtime();
  
  /* Define simulation specific parameters */
  AppCtx user; 
  user.dt       = 1.0e-7; 
  user.ta=0.0;
  user.tb=0.0;
  user.tc=0.0;
  /* Set discretization options */
  PetscInt nsteps = 1;
  PetscInt N=2, p=2, Corder=PETSC_DECIDE, resStep=0;
  PetscBool output = PETSC_FALSE; 
  PetscBool monitor = PETSC_FALSE; 
  char filePrefix[PETSC_MAX_PATH_LEN] = {0};
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","StrainGradient3D Options","IGA");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-N","number of elements (along one dimension)",__FILE__,N,&N,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-p","polynomial order",__FILE__,p,&p,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-C","global continuity order",__FILE__,Corder,&Corder,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-file_prefix","File Prefix",__FILE__,filePrefix,filePrefix,sizeof(filePrefix),PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-ch_output","Enable output files",__FILE__,output,&output,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-ch_monitor","Compute and show statistics of solution",__FILE__,monitor,&monitor,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-nsteps","Number of load steps to take",__FILE__,nsteps,&nsteps,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-dt","time step",__FILE__,user.dt,&user.dt,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-res_step","Restart Step",__FILE__,resStep,&resStep,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  if (Corder == PETSC_DECIDE) Corder = p-1;

  //
  if (p < 2 || Corder < 0) /* Problem requires a p>=2 C1 basis */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Problem requires minimum of p = 2");
  if (p <= Corder)         /* Check C < p */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Discretization inconsistent: polynomial order must be greater than degree of continuity");
  
  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDim(iga,DIM);CHKERRQ(ierr);
  ierr = IGASetDof(iga,DIM);CHKERRQ(ierr);

  IGAAxis axis0;
  ierr = IGAGetAxis(iga,0,&axis0);CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis0,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis0,N,0.0,1.0,Corder);CHKERRQ(ierr);
  IGAAxis axis1;
  ierr = IGAGetAxis(iga,1,&axis1);CHKERRQ(ierr);
  ierr = IGAAxisCopy(axis0,axis1);CHKERRQ(ierr);
#if DIM==3
  IGAAxis axis2;
  ierr = IGAGetAxis(iga,2,&axis2);CHKERRQ(ierr);
  ierr = IGAAxisCopy(axis0,axis2);CHKERRQ(ierr);
#endif
  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  user.iga = iga;
  if (resStep==0){
    char meshfilename[256];
    sprintf(meshfilename, "mesh%s.dat", filePrefix);
    PetscPrintf(PETSC_COMM_WORLD,"\nWriting mesh file: %s\n", meshfilename);
    ierr = IGAWrite(iga, meshfilename);CHKERRQ(ierr);
  }
  
  Vec U,U0;
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
    //ierr = FormInitialCondition(iga, U0, &user); 
  }
  ierr = VecCopy(U0, U);CHKERRQ(ierr);
 
  //
  IGAForm form;
  ierr = IGAGetForm(iga,&form);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,1,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,1,PETSC_TRUE);CHKERRQ(ierr);
#if DIM==3
  ierr = IGAFormSetBoundaryForm (form,2,0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,2,1,PETSC_TRUE);CHKERRQ(ierr);
#endif

  //Dirichlet BC
  double dVal=0.1*1e-15;
  //along z=0 face
  ierr = IGASetBoundaryValue(iga,0,0,0,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(iga,0,0,1,0.0);CHKERRQ(ierr);
#if DIM==3
  ierr = IGASetBoundaryValue(iga,0,0,2,0.0);CHKERRQ(ierr);
#endif
  //along z=1 face
  ierr = IGASetBoundaryValue(iga,0,1,0,-dVal);CHKERRQ(ierr);

  // Setup the nonlinear solver
  ierr = IGASetFormIEFunction(iga,Residual,&user);CHKERRQ(ierr);
  ierr = IGASetFormIEJacobian(iga,Jacobian,&user);CHKERRQ(ierr);
     
   //
  TS ts;
  ierr = IGACreateTS(iga,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  ierr = TSSetDuration(ts,1,1.0);CHKERRQ(ierr);
  ierr = TSSetTime(ts,0.0);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,user.dt);CHKERRQ(ierr);
  if (output) {
    ierr = TSMonitorSet(ts,OutputMonitor,&user,NULL);CHKERRQ(ierr);
  }

  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
#if PETSC_VERSION_LE(3,3,0)
  ierr = TSSolve(ts,U,NULL);CHKERRQ(ierr);
#else
  ierr = TSSolve(ts,U);CHKERRQ(ierr);
#endif
  //
  PetscPrintf(PETSC_COMM_WORLD,"\n\n\nta:%12.4f, tb:%12.4f, tc:%12.4f\n\n\n",user.ta,user.tb,user.tc);
  //
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&U0);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = IGADestroy(&iga);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

