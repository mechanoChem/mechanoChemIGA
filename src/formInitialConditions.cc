//extern "C" {
#include "petiga.h"
//}
#include "coreFunctions.h"

#undef  __FUNCT__
#define __FUNCT__ "FormInitialCondition"
template<unsigned int DIM>
PetscErrorCode FormInitialCondition(IGA iga, Vec U, AppCtx<DIM> *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  srand(myrank);

  DM da;
  unsigned int nSclrs = user->scalarSolnFields.size(), nVctrs = user->vectorSolnFields.size();
  PetscInt dof = nSclrs + DIM*nVctrs;
  ierr = IGACreateNodeDM(iga,dof,&da);CHKERRQ(ierr);

  Tensor<1,DIM,double> x;
  if (DIM == 2){
    PetscScalar ***u;
    ierr = DMDAVecGetArrayDOF(da,U,&u);CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
    PetscInt i,j;
    PetscReal imax = (info.mx-1), Lx = user->L[0]*user->GridScale;
    PetscReal jmax = (info.my-1), Ly = user->L[1]*user->GridScale;
    for(i=info.xs;i<info.xs+info.xm;i++){
      for(j=info.ys;j<info.ys+info.ym;j++){
	x[0] = (i/imax)*Lx;
	x[1] = (j/jmax)*Ly;
	for(unsigned int l=0; l<nSclrs; ++l){
	  u[j][i][l] = user->scalarInitialConditions(x,l,*user);
	}
	for(unsigned int l1=0; l1<nVctrs; ++l1){
	  Tensor<1,DIM,double> vecIC = user->vectorInitialConditions(x,l1,*user);
	  for(unsigned int l2=0; l2<DIM; ++l2){
	    u[j][i][nSclrs+DIM*l1+l2] = vecIC[l2];
	  }
	}
      }
    }  
    ierr = DMDAVecRestoreArrayDOF(da,U,&u);CHKERRQ(ierr);

    ierr = DMDestroy(&da);CHKERRQ(ierr);
  }
  else if (DIM == 3){
    PetscScalar ****u;
    ierr = DMDAVecGetArrayDOF(da,U,&u);CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
    PetscInt i,j,k;
    PetscReal imax = (info.mx-1), Lx = user->L[0]*user->GridScale;
    PetscReal jmax = (info.my-1), Ly = user->L[1]*user->GridScale;
    PetscReal kmax = (info.mz-1), Lz = user->L[2]*user->GridScale;
    for(i=info.xs;i<info.xs+info.xm;i++){
      for(j=info.ys;j<info.ys+info.ym;j++){
	for(k=info.zs;k<info.zs+info.zm;k++){
	  x[0] = (i/imax)*Lx;
	  x[1] = (j/jmax)*Ly;
	  x[2] = (k/kmax)*Lz;
	  for(unsigned int l=0; l<nSclrs; ++l){
	    u[k][j][i][l] = user->scalarInitialConditions(x,l,*user);
	  }
	  for(unsigned int l1=0; l1<nVctrs; ++l1){
	    Tensor<1,DIM,double> vecIC = user->vectorInitialConditions(x,l1,*user);
	    for(unsigned int l2=0; l2<DIM; ++l2){
	      u[k][j][i][nSclrs+DIM*l1+l2] = vecIC[l2];
	    }
	  }
	}
      }
    }  
    ierr = DMDAVecRestoreArrayDOF(da,U,&u);CHKERRQ(ierr);

    ierr = DMDestroy(&da);CHKERRQ(ierr); 
  }

  PetscFunctionReturn(0); 
}

template PetscErrorCode FormInitialCondition<2>(IGA iga, Vec U, AppCtx<2> *user);
template PetscErrorCode FormInitialCondition<3>(IGA iga, Vec U, AppCtx<3> *user);
