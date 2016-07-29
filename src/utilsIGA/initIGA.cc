//extern "C" {
#include "petiga.h"
//}
#include "IBVPHeaders.h"
#include "utilsIGAHeaders.h"

template<unsigned int DIM, unsigned int DOF>
int initIGA(AppCtx& user, PetscInt p){
  PetscErrorCode  ierr;

  //set discretization options
  PetscInt C=PETSC_DECIDE;
  PetscBool output = PETSC_TRUE; 
  PetscBool monitor = PETSC_TRUE; 
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","IGA Options","IGA");CHKERRQ(ierr);
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
  ierr = IGASetDof(user.iga,DOF);CHKERRQ(ierr);

  IGAAxis axis0;
  ierr = IGAGetAxis(user.iga,0,&axis0);CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis0,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis0,user.Nx,0.0,user.Lx*user.GridScale,C);CHKERRQ(ierr); //x direction

  IGAAxis axis1;
  ierr = IGAGetAxis(user.iga,1,&axis1);CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis1,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis1,user.Ny,0.0,user.Ly*user.GridScale,C);CHKERRQ(ierr); //y direction

	if(DIM==3){
		IGAAxis axis2;
		ierr = IGAGetAxis(user.iga,2,&axis2);CHKERRQ(ierr);
		ierr = IGAAxisSetDegree(axis2,p);CHKERRQ(ierr);
		ierr = IGAAxisInitUniform(axis2,user.Nz,0.0,user.Lz*user.GridScale,C);CHKERRQ(ierr); //z direction
	}
  ierr = IGASetFromOptions(user.iga);CHKERRQ(ierr);
  ierr = IGASetUp(user.iga);CHKERRQ(ierr);

  //restart
  char meshfilename[256];
  sprintf(meshfilename, "mesh.dat");
  PetscPrintf(PETSC_COMM_WORLD,"\nWriting mesh file: %s\n", meshfilename);
  ierr = IGAWrite(user.iga, meshfilename);CHKERRQ(ierr);
  
  //
  IGAForm form;
  ierr = IGAGetForm(user.iga,&form);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,1,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,1,1,PETSC_TRUE);CHKERRQ(ierr);
	if(DIM==3){
	  ierr = IGAFormSetBoundaryForm (form,2,0,PETSC_TRUE);CHKERRQ(ierr);
  	ierr = IGAFormSetBoundaryForm (form,2,1,PETSC_TRUE);CHKERRQ(ierr);
	}
  //assign residual and jacobian functions
  ierr = IGASetFormIEFunction(user.iga,Residual<DIM,DOF>,&user);CHKERRQ(ierr);
  ierr = IGASetFormIEJacobian(user.iga,Jacobian<DIM,DOF>,&user);CHKERRQ(ierr);
  //
  return 0;
}

template int initIGA<2,2>(AppCtx& user, PetscInt p);
template int initIGA<2,3>(AppCtx& user, PetscInt p);
template int initIGA<2,4>(AppCtx& user, PetscInt p);
template int initIGA<3,3>(AppCtx& user, PetscInt p);
template int initIGA<3,4>(AppCtx& user, PetscInt p);
template int initIGA<3,6>(AppCtx& user, PetscInt p);
