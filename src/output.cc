/*
Note that the functions IGAElementNextFormFunction and 
IGAComputeProjectionFunction in the file src/output.cc, 
as well as all functions in the file src/petigasnes_mod.h 
were derived from the PetIGA/src/petigasnes.c source code
in the PetIGA library [https://bitbucket.org/dalcinl/petiga/].
Accordingly, we include the license/copyright notice for the 
PetIGA library here in the file LICENSE_PetIGA to apply
to the above functions.
*/

//extern "C" {
#include "petiga.h"
//}
#include <cmath>

#include "coreFunctions.h"
#include <iostream>

#define PI 3.14159265

PETSC_STATIC_INLINE
PetscBool IGAElementNextFormFunction(IGAElement element,IGAFormFunction *fun,void **ctx)
{
  IGAForm form = element->parent->form;
  if (!IGAElementNextForm(element,form->visit)) return PETSC_FALSE;
  *fun = form->ops->Function;
  *ctx = form->ops->FunCtx;
  return PETSC_TRUE;
}

#undef  __FUNCT__
#define __FUNCT__ "IGAComputeProjectionFunction"
PetscErrorCode IGAComputeProjectionFunction(IGA iga,IGA igaProject,Vec vecU,Vec vecF)
{
  Vec               localU;
  const PetscScalar *arrayU;
  IGAElement        element;
  IGAElement        elementP;
  IGAPoint          point;
  IGAFormFunction   Function;
  void              *ctx;
  PetscScalar       *U,*F,*R;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(iga,IGA_CLASSID,1);
  PetscValidHeaderSpecific(igaProject,IGA_CLASSID,1);
  PetscValidHeaderSpecific(vecU,VEC_CLASSID,2);
  PetscValidHeaderSpecific(vecF,VEC_CLASSID,3);
  IGACheckSetUp(iga,1);
  IGACheckSetUp(igaProject,1);
  IGACheckFormOp(igaProject,1,Function);

  /* Clear global vector F*/
  ierr = VecZeroEntries(vecF);CHKERRQ(ierr);

  /* Get local vector U and array */
  ierr = IGAGetLocalVecArray(iga,vecU,&localU,&arrayU);CHKERRQ(ierr);

  /* Element loop */
  ierr = IGABeginElement(iga,&element);CHKERRQ(ierr);
  ierr = IGABeginElement(igaProject,&elementP);CHKERRQ(ierr);
  while (IGANextElement(igaProject,elementP)) {
    IGANextElement(iga,element);

    ierr = IGAElementGetWorkVec(elementP,&F);CHKERRQ(ierr);
    ierr = IGAElementGetValues(element,arrayU,&U);CHKERRQ(ierr);
    /* FormFunction loop */
    while (IGAElementNextFormFunction(elementP,&Function,&ctx)) {
      /* Quadrature loop */
      ierr = IGAElementBeginPoint(elementP,&point);CHKERRQ(ierr);
      while (IGAElementNextPoint(elementP,point)) {
        ierr = IGAPointGetWorkVec(point,&R);CHKERRQ(ierr);
        ierr = Function(point,U,R,ctx);CHKERRQ(ierr);
        ierr = IGAPointAddVec(point,R,F);CHKERRQ(ierr);
      }
      ierr = IGAElementEndPoint(elementP,&point);CHKERRQ(ierr);
    }
    ierr = IGAElementAssembleVec(elementP,F,vecF);CHKERRQ(ierr);
  }
  ierr = IGAEndElement(igaProject,&elementP);CHKERRQ(ierr);
  IGANextElement(iga,element);
  ierr = IGAEndElement(iga,&element);CHKERRQ(ierr);

  /* Restore local vector U and array */
  ierr = IGARestoreLocalVecArray(iga,vecU,&localU,&arrayU);CHKERRQ(ierr);

  /* Assemble global vector F */
  ierr = VecAssemblyBegin(vecF);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(vecF);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "ProjectionResidual"
template <int DIM>
PetscErrorCode ProjectionResidual(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *ctx)
{	
  AppCtx<DIM> *user = (AppCtx<DIM> *)ctx;

  PetscErrorCode ierr;
  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  double (*n) = (double (*)) p->normal;
  double (*X) = (double (*)) p->point;
  Tensor<1,DIM,double> normal, x;
  for(unsigned int i=0; i<DIM; ++i){
    normal[i] = n[i];
    x[i] = X[i];
  } 

  unsigned int nScalars = user->scalarSolnFields.size(), nVectors = user->vectorSolnFields.size();
  unsigned int nSProj = user->scalarProjectnFields.size(), nVProj = user->vectorProjectnFields.size();
  PetscInt UDOF = nScalars + DIM*nVectors;

  solutionScalars<DIM,double> solnS(nScalars);
  solutionVectors<DIM,double> solnV(nVectors);

  //scalar field variables
  double c, cx[DIM], cxx[DIM][DIM];
  for(unsigned int i=0; i<nScalars; ++i){
    ComputeField<PetscReal,DIM>(SCALAR,i,p,U,UDOF,&c,&cx[0],&cxx[0][0]);
    solnS.defineVal(i,c);
    solnS.defineGrad(i,cx);
    solnS.defineHess(i,cxx);
  }

  //vector field variables
  double u[DIM], ux[DIM][DIM], uxx[DIM][DIM][DIM];
  for(unsigned int i=0; i<nVectors; ++i){
    ComputeField<PetscReal,DIM>(VECTOR,nScalars+i*DIM,p,U,UDOF,&u[0],&ux[0][0],&uxx[0][0][0]);
    solnV.defineVal(i,u);
    solnV.defineGrad(i,ux);
    solnV.defineHess(i,uxx);
  }

  //Initialize projection vectors
  std::vector<double> sProjections(nSProj);
  std::vector<Tensor<1,DIM,double> > vProjections(nVProj);

  //store L2 projection residual
  const PetscReal (*N) = (PetscReal (*)) p->shape[0]; 
  for(int n1=0; n1<nen; n1++){
    user->projectFields(x,normal,solnS,solnV,*user,sProjections,vProjections);
    for (unsigned int i=0; i<nSProj; i++){
      R[n1*dof+i] = N[n1]*sProjections[i];
    }
    for (unsigned int i=0; i<nVProj; i++){
      for (unsigned int j=0; j<DIM; j++){
	R[n1*dof+nSProj+i*DIM+j] = N[n1]*vProjections[i][j];
      }
    }
  }

  return 0;
}
template PetscErrorCode ProjectionResidual<1>(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *ctx);
template PetscErrorCode ProjectionResidual<2>(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *ctx);
template PetscErrorCode ProjectionResidual<3>(IGAPoint p, const PetscScalar *U, PetscScalar *R, void *ctx);

#undef  __FUNCT__
#define __FUNCT__ "ProjectionJacobian"
PetscErrorCode ProjectionJacobian(IGAPoint p, PetscScalar *K, void *ctx)
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
template <int dim>
PetscErrorCode ProjectSolution(IGA iga, PetscInt step, Vec U, AppCtx<dim> *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;

  //Setup linear system for L2 Projection
  Mat A;
  Vec x,b;
  ierr = IGACreateMat(iga,&A);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&x);CHKERRQ(ierr);
  ierr = IGACreateVec(iga,&b);CHKERRQ(ierr);
  ierr = IGASetFormFunction(iga,ProjectionResidual<dim>,user);CHKERRQ(ierr);
  ierr = IGASetFormMatrix(iga,ProjectionJacobian,user);CHKERRQ(ierr);
  ierr = IGAComputeProjectionFunction(user->iga,iga,U,b);CHKERRQ(ierr);
  ierr = IGAComputeMatrix(iga,A);CHKERRQ(ierr);

  //*
  //Solver
  KSP ksp;
  ierr = IGACreateKSP(iga,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
  //write solution
  char filename[256];
  sprintf(filename,"%s/outE%d.dat",user->outputDir.c_str(),step);
  ierr = IGAWriteVec(iga,x,filename);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  // */
  PetscFunctionReturn(0); 
}
template PetscErrorCode ProjectSolution<1>(IGA iga, PetscInt step, Vec U, AppCtx<1> *user);
template PetscErrorCode ProjectSolution<2>(IGA iga, PetscInt step, Vec U, AppCtx<2> *user);
template PetscErrorCode ProjectSolution<3>(IGA iga, PetscInt step, Vec U, AppCtx<3> *user);

#undef  __FUNCT__
#define __FUNCT__ "StepUpdate"
template <int dim>
PetscErrorCode StepUpdate(PetscInt it_number,PetscReal c_time,Vec U,AppCtx<dim> &user)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt step = it_number;
 
  //output to file
  char filename[256];
  if (step%user.skipOutput==0){
    sprintf(filename,"%s/outU%d.dat",user.outputDir.c_str(),step);
    ierr = IGAWriteVec(user.iga,U,filename);CHKERRQ(ierr);

    if(user.scalarProjectnFields.size()+user.vectorProjectnFields.size() > 0){
      ProjectSolution<dim>(user.igaProject, step, U, &user); 
    }
  }

  // Call adaptive stepping functions
  user.loadStep(step,user);
  user.adaptiveTimeStep(step,user);

  PetscFunctionReturn(0);
}
template PetscErrorCode StepUpdate<1>(PetscInt it_number,PetscReal c_time,Vec U,AppCtx<1> &user);
template PetscErrorCode StepUpdate<2>(PetscInt it_number,PetscReal c_time,Vec U,AppCtx<2> &user);
template PetscErrorCode StepUpdate<3>(PetscInt it_number,PetscReal c_time,Vec U,AppCtx<3> &user);

