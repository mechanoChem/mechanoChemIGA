//extern "C" {
#include "petiga.h"
//}
#include "coreFunctions.h"

template<unsigned int DIM>
int InitIGA(AppCtx<DIM>& user, IGA &iga, PetscInt p, PetscInt C, PetscInt DOF){
  PetscErrorCode  ierr;

  //set discretization options
  PetscBool output = PETSC_TRUE; 
  PetscBool monitor = PETSC_TRUE; 
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","IGA Options","IGA");CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  if (C == PETSC_DECIDE) C = p-1;
 
  if (p <= C)         /* Check C < p */
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Discretization inconsistent: polynomial order must be greater than degree of continuity");
 
  //
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDim(iga,DIM);CHKERRQ(ierr);
  ierr = IGASetDof(iga,DOF);CHKERRQ(ierr);

  IGAAxis axis0;
  ierr = IGAGetAxis(iga,0,&axis0);CHKERRQ(ierr);
  ierr = IGAAxisSetPeriodic(axis0,user.periodic[0]);CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis0,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis0,user.N[0],0.0,user.L[0]*user.GridScale,C);CHKERRQ(ierr); //x direction

  if(DIM>=2){
    IGAAxis axis1;
    ierr = IGAGetAxis(iga,1,&axis1);CHKERRQ(ierr);
    ierr = IGAAxisSetPeriodic(axis1,user.periodic[1]);CHKERRQ(ierr);
    ierr = IGAAxisSetDegree(axis1,p);CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axis1,user.N[1],0.0,user.L[1]*user.GridScale,C);CHKERRQ(ierr); //y direction

    if(DIM==3){
      IGAAxis axis2;
      ierr = IGAGetAxis(iga,2,&axis2);CHKERRQ(ierr);
      ierr = IGAAxisSetPeriodic(axis2,user.periodic[2]);CHKERRQ(ierr);
      ierr = IGAAxisSetDegree(axis2,p);CHKERRQ(ierr);
      ierr = IGAAxisInitUniform(axis2,user.N[2],0.0,user.L[2]*user.GridScale,C);CHKERRQ(ierr); //z direction
    }
  }
  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  
  IGAForm form;
  ierr = IGAGetForm(iga,&form);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,0,PETSC_TRUE);CHKERRQ(ierr);
  ierr = IGAFormSetBoundaryForm (form,0,1,PETSC_TRUE);CHKERRQ(ierr);
  if(DIM>=2){
    ierr = IGAFormSetBoundaryForm (form,1,0,PETSC_TRUE);CHKERRQ(ierr);
    ierr = IGAFormSetBoundaryForm (form,1,1,PETSC_TRUE);CHKERRQ(ierr);
    if(DIM==3){
      ierr = IGAFormSetBoundaryForm (form,2,0,PETSC_TRUE);CHKERRQ(ierr);
      ierr = IGAFormSetBoundaryForm (form,2,1,PETSC_TRUE);CHKERRQ(ierr);
    }
  }

  //Set mat type
  ierr = IGASetMatType(iga,MATAIJ);CHKERRQ(ierr); //For superlu_dist (still works for gmres, etc.)

  return 0;
}

template int InitIGA<1>(AppCtx<1>& user, IGA &iga, PetscInt p, PetscInt C, PetscInt DOF);
template int InitIGA<2>(AppCtx<2>& user, IGA &iga, PetscInt p, PetscInt C, PetscInt DOF);
template int InitIGA<3>(AppCtx<3>& user, IGA &iga, PetscInt p, PetscInt C, PetscInt DOF);
