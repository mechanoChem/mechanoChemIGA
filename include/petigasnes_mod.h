#ifndef PETIGASNES_MOD_H
#define PETIGASNES_MOD_H

#include "petiga.h"
#include "appCtx.h"
#include "coreFunctions.h"

PETSC_STATIC_INLINE
PetscBool IGAElementNextFormFunction(IGAElement element,void **ctx)
{
  IGAForm form = element->parent->form;
  if (!IGAElementNextForm(element,form->visit)) return PETSC_FALSE;
  *ctx = form->ops->FunCtx;
  return PETSC_TRUE;
}

PETSC_STATIC_INLINE
PetscBool IGAElementNextFormJacobian(IGAElement element,void **ctx)
{
  IGAForm form = element->parent->form;
  if (!IGAElementNextForm(element,form->visit)) return PETSC_FALSE;
  *ctx = form->ops->JacCtx;
  return PETSC_TRUE;
}

#undef  __FUNCT__
#define __FUNCT__ "IGAComputeFunction_mod"
template<unsigned int DIM>
PetscErrorCode IGAComputeFunction_mod(IGA iga,Vec vecU,Vec vecF)
{
  Vec               localU;
  const PetscScalar *arrayU;
  IGAElement        element;
  IGAPoint          point;
  void              *ctx;
  PetscScalar       *U,*F,*R;
  PetscErrorCode    ierr;

  /* For time stepping */
  Vec vecUp, vecUpp;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(vecU,VEC_CLASSID,2);
  PetscValidHeaderSpecific(vecF,VEC_CLASSID,3);

  /* Get previous solution vectors */
  ctx = iga->form->ops->FunCtx;
  AppCtx<DIM> *user = (AppCtx<DIM> *)ctx;
  vecUp = *(user->Up);
  vecUpp = *(user->Upp);
  Vec               localUp, localUpp;
  const PetscScalar *arrayUp, *arrayUpp;
  PetscScalar       *Up, *Upp;

  /* Clear global vector F*/
  ierr = VecZeroEntries(vecF);CHKERRQ(ierr);

  /* Get local vector U and array */
  ierr = IGAGetLocalVecArray(iga,vecU,&localU,&arrayU);CHKERRQ(ierr);
  ierr = IGAGetLocalVecArray(iga,vecUp,&localUp,&arrayUp);CHKERRQ(ierr);
  ierr = IGAGetLocalVecArray(iga,vecUpp,&localUpp,&arrayUpp);CHKERRQ(ierr);

  /* Element loop */
  ierr = IGABeginElement(iga,&element);CHKERRQ(ierr);
  while (IGANextElement(iga,element)) {
    ierr = IGAElementGetWorkVec(element,&F);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayU,&U);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayUp,&Up);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayUp,&Upp);CHKERRQ(ierr);
    ierr = IGAElementFixValues(element,U);CHKERRQ(ierr); //Apply BCs (I believe...)
    //ierr = IGAElementFixValues(element,Up);CHKERRQ(ierr);
    //ierr = IGAElementFixValues(element,Upp);CHKERRQ(ierr);
    /* FormFunction loop */
    while (IGAElementNextFormFunction(element,&ctx)) {
      /* Quadrature loop */
      ierr = IGAElementBeginPoint(element,&point);CHKERRQ(ierr);
      while (IGAElementNextPoint(element,point)) {
	ierr = IGAPointGetWorkVec(point,&R);CHKERRQ(ierr);
	ierr = Residual<DIM>(point,U,Up,Upp,R,ctx);CHKERRQ(ierr);
	ierr = IGAPointAddVec(point,R,F);CHKERRQ(ierr);
      }
      ierr = IGAElementEndPoint(element,&point);CHKERRQ(ierr);
    }
    ierr = IGAElementFixFunction(element,F);CHKERRQ(ierr);
    ierr = IGAElementAssembleVec(element,F,vecF);CHKERRQ(ierr);
  }
  ierr = IGAEndElement(iga,&element);CHKERRQ(ierr);

  /* Restore local vector U and array */
  ierr = IGARestoreLocalVecArray(iga,vecU,&localU,&arrayU);CHKERRQ(ierr);
  ierr = IGARestoreLocalVecArray(iga,vecUp,&localUp,&arrayUp);CHKERRQ(ierr); //I don't think this is necessary, since we don't change Up or Upp
  ierr = IGARestoreLocalVecArray(iga,vecUpp,&localUpp,&arrayUpp);CHKERRQ(ierr);

  /* Assemble global vector F */
  ierr = VecAssemblyBegin(vecF);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vecF);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGAComputeJacobian_mod"
template<unsigned int DIM>
PetscErrorCode IGAComputeJacobian_mod(IGA iga,Vec vecU,Mat matJ)
{
  Vec               localU;
  const PetscScalar *arrayU;
  IGAElement        element;
  IGAPoint          point;
  void              *ctx;
  PetscScalar       *U,*J,*K;
  PetscErrorCode    ierr;

  /* For time stepping */
  Vec vecUp, vecUpp;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(vecU,VEC_CLASSID,2);
  PetscValidHeaderSpecific(matJ,MAT_CLASSID,3);

  /* Get previous solution vectors */
  ctx = iga->form->ops->JacCtx;
  AppCtx<DIM> *user = (AppCtx<DIM> *)ctx;
  vecUp = *(user->Up);
  vecUpp = *(user->Upp);
  Vec               localUp, localUpp;
  const PetscScalar *arrayUp, *arrayUpp;
  PetscScalar       *Up, *Upp;


  /* Clear global matrix J */
  ierr = MatZeroEntries(matJ);CHKERRQ(ierr);

  /* Get local vector U and array */
  ierr = IGAGetLocalVecArray(iga,vecU,&localU,&arrayU);CHKERRQ(ierr);
  ierr = IGAGetLocalVecArray(iga,vecUp,&localUp,&arrayUp);CHKERRQ(ierr);
  ierr = IGAGetLocalVecArray(iga,vecUpp,&localUpp,&arrayUpp);CHKERRQ(ierr);

  /* Element Loop */
  ierr = IGABeginElement(iga,&element);CHKERRQ(ierr);
  while (IGANextElement(iga,element)) {
    ierr = IGAElementGetWorkMat(element,&J);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayU,&U);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayUp,&Up);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayUpp,&Upp);CHKERRQ(ierr);
    ierr = IGAElementFixValues(element,U);CHKERRQ(ierr);
    //ierr = IGAElementFixValues(element,Up);CHKERRQ(ierr);
    //ierr = IGAElementFixValues(element,Upp);CHKERRQ(ierr);
    /* FormFunction loop */
    while (IGAElementNextFormJacobian(element,&ctx)) {
      /* Quadrature loop */
      ierr = IGAElementBeginPoint(element,&point);CHKERRQ(ierr);
      while (IGAElementNextPoint(element,point)) {
	ierr = IGAPointGetWorkMat(point,&K);CHKERRQ(ierr);
	ierr = Jacobian<DIM>(point,U,Up,Upp,K,ctx);CHKERRQ(ierr);
	ierr = IGAPointAddMat(point,K,J);CHKERRQ(ierr);
      }
      ierr = IGAElementEndPoint(element,&point);CHKERRQ(ierr);
    }
    ierr = IGAElementFixJacobian(element,J);CHKERRQ(ierr);
    ierr = IGAElementAssembleMat(element,J,matJ);CHKERRQ(ierr);
  }
  ierr = IGAEndElement(iga,&element);CHKERRQ(ierr);

  /* Restore local vector U and array */
  ierr = IGARestoreLocalVecArray(iga,vecU,&localU,&arrayU);CHKERRQ(ierr);
  ierr = IGARestoreLocalVecArray(iga,vecUp,&localUp,&arrayUp);CHKERRQ(ierr);
  ierr = IGARestoreLocalVecArray(iga,vecUpp,&localUpp,&arrayUpp);CHKERRQ(ierr);

  /* Assemble global matrix J*/
  ierr = MatAssemblyBegin(matJ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (matJ,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGASNESFormFunction_mod"
template<unsigned int DIM>
PetscErrorCode IGASNESFormFunction_mod(SNES snes,Vec U,Vec F,void *ctx)
{
  IGA            iga = (IGA)ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(snes,SNES_CLASSID,1);
  PetscValidHeaderSpecific(U,VEC_CLASSID,2);
  PetscValidHeaderSpecific(F,VEC_CLASSID,3);
  PetscValidHeaderSpecific(iga,IGA_CLASSID,4);

  ierr = IGAComputeFunction_mod<DIM>(iga,U,F);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGASNESFormJacobian_mod"
template<unsigned int DIM>
PetscErrorCode IGASNESFormJacobian_mod(SNES snes,Vec U,Mat J,Mat P,void *ctx)
{
  IGA            iga = (IGA)ctx;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(snes,SNES_CLASSID,1);
  PetscValidHeaderSpecific(U,VEC_CLASSID,2);
  PetscValidHeaderSpecific(J,MAT_CLASSID,3);
  PetscValidHeaderSpecific(P,MAT_CLASSID,4);
  PetscValidHeaderSpecific(iga,IGA_CLASSID,6);

  ierr = IGAComputeJacobian_mod<DIM>(iga,U,P);CHKERRQ(ierr);
  if (J != P) {
    ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "IGACreateSNES_mod"
template<unsigned int DIM>
PetscErrorCode IGACreateSNES_mod(IGA iga,SNES *snes)
{
  MPI_Comm       comm;
  Vec            F;
  Mat            J;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidPointer(snes,2);

  ierr = IGAGetComm(iga,&comm);CHKERRQ(ierr);
  ierr = SNESCreate(comm,snes);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)*snes,"IGA",(PetscObject)iga);CHKERRQ(ierr);
  ierr = IGASetOptionsHandlerSNES(*snes);CHKERRQ(ierr);

  ierr = IGACreateVec(iga,&F);CHKERRQ(ierr);
  ierr = SNESSetFunction(*snes,F,IGASNESFormFunction_mod<DIM>,iga);CHKERRQ(ierr);
  ierr = VecDestroy(&F);CHKERRQ(ierr);

  ierr = IGACreateMat(iga,&J);CHKERRQ(ierr);
  ierr = SNESSetJacobian(*snes,J,J,IGASNESFormJacobian_mod<DIM>,iga);CHKERRQ(ierr);
  ierr = MatDestroy(&J);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

#endif //PETIGASNES_MOD_H
