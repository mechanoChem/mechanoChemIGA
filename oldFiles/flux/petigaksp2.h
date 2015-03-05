#if !defined(PETIGAKSP2_H)
#define PETIGAKSP2_H
#include "petiga.h"
typedef struct { 
  Vec *U0;
  PetscScalar *localU0;
  PetscInt test;
  PetscReal Es;
} AppCtxKSP;

PETSC_STATIC_INLINE
PetscBool IGAElementNextFormSystem(IGAElement element,IGAFormSystem *sys,void **ctx)
{
  IGAForm form = element->parent->form;
  if (!IGAElementNextForm(element,form->visit)) return PETSC_FALSE;
  *sys = form->ops->System;
  *ctx = form->ops->SysCtx;
  return PETSC_TRUE;
}

#undef  __FUNCT__
#define __FUNCT__ "IGAComputeSystem2"
PetscErrorCode IGAComputeSystem2(IGA iga,Mat matA,Vec vecB)
{
  IGAElement     element;
  IGAPoint       point;
  IGAFormSystem  System;
  AppCtxKSP      *ctx;
  PetscScalar    *A,*B,*U0;
  PetscScalar    *K,*F;
  PetscErrorCode ierr;
  Vec               localU0;
  const PetscScalar *arrayU0;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(matA,MAT_CLASSID,2);
  PetscValidHeaderSpecific(vecB,VEC_CLASSID,3);
  IGACheckSetUp(iga,1);
  IGACheckFormOp(iga,1,System);

  ierr = MatZeroEntries(matA);CHKERRQ(ierr);
  ierr = VecZeroEntries(vecB);CHKERRQ(ierr);

  ierr = PetscLogEventBegin(IGA_FormSystem,iga,matA,vecB,0);CHKERRQ(ierr);
  
  /* Element loop */
  ierr = IGABeginElement(iga,&element);CHKERRQ(ierr);
  IGAForm form = element->parent->form;
  ctx = (AppCtxKSP*) form->ops->SysCtx;
  //PetscPrintf(PETSC_COMM_WORLD,"%d \n", ctx->test);
  ierr = IGAGetLocalVecArray(iga,*(ctx->U0),&localU0,&arrayU0);CHKERRQ(ierr);
  while (IGANextElement(iga,element)) {
    ierr = IGAElementGetWorkMat(element,&A);CHKERRQ(ierr);
    ierr = IGAElementGetWorkVec(element,&B);CHKERRQ(ierr);
    //ierr = IGAElementGetWorkVal(element,&U0);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayU0,&U0);CHKERRQ(ierr);
    ierr = IGAElementFixValues(element,U0);CHKERRQ(ierr);
    /* FormSystem loop */
    while (IGAElementNextFormSystem(element,&System,(void**)&ctx)) {
      ctx->localU0=U0;
      /* Quadrature loop */
      ierr = IGAElementBeginPoint(element,&point);CHKERRQ(ierr);
      while (IGAElementNextPoint(element,point)) {
        ierr = IGAPointGetWorkMat(point,&K);CHKERRQ(ierr);
        ierr = IGAPointGetWorkVec(point,&F);CHKERRQ(ierr);	
        ierr = System(point,K,F,ctx);CHKERRQ(ierr);
        ierr = IGAPointAddMat(point,K,A);CHKERRQ(ierr);
        ierr = IGAPointAddVec(point,F,B);CHKERRQ(ierr);
      }
      ierr = IGAElementEndPoint(element,&point);CHKERRQ(ierr);
    }
    ierr = IGAElementFixSystem(element,A,B);CHKERRQ(ierr);
    ierr = IGAElementAssembleMat(element,A,matA);CHKERRQ(ierr);
    ierr = IGAElementAssembleVec(element,B,vecB);CHKERRQ(ierr);
  }
  ierr = IGAEndElement(iga,&element);CHKERRQ(ierr);

  ierr = IGARestoreLocalVecArray(iga,*(ctx->U0),&localU0,&arrayU0);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(IGA_FormSystem,iga,matA,vecB,0);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(matA,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (matA,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(vecB);CHKERRQ(ierr);
  ierr = VecAssemblyEnd  (vecB);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
#endif/*PETIGAKSP2_H*/
