//extern "C" {
#include "petiga.h"
//}
#include "coreFunctions.h"
#include "userFunctions.h"
#include "solutionClass.h"
#include "testFunctionClass.h"

//include automatic differentiation library
#include <Sacado.hpp>
//#include <Sacado_Fad_BLAS.hpp>
#include <Sacado_Fad_SimpleFad.hpp>

//residual function implementation
#undef  __FUNCT__
#define __FUNCT__ "QuadPtResidual"
template <class T,unsigned int DIM>
PetscErrorCode QuadPtResidual(IGAPoint p,
			      const T * U,
			      const PetscScalar * Up,
			      const PetscScalar * Upp,
			      T *R,
			      void *ctx)
{
  typedef Sacado::Fad::SimpleFad<T> AD;

  AppCtx<DIM> *user = (AppCtx<DIM> *)ctx;

  PetscInt nen, dof;
  IGAPointGetSizes(p,0,&nen,&dof);
  unsigned int nScalars = user->scalarSolnFields.size(), nVectors = user->vectorSolnFields.size();

  PetscReal *n = p->normal;
  double (*X) = (double (*)) p->point;
  Tensor<1,DIM,double> normal, x;
  for(unsigned int i=0; i<DIM; ++i){
    normal[i] = n[i];
    x[i] = X[i];
  } 

  solutionScalars<DIM,T> solnS(nScalars);
  solutionVectors<DIM,T> solnV(nVectors);


  //scalar field variables
  T c, cx[DIM], cxx[DIM][DIM];
  double cp, cpx[DIM], cpxx[DIM][DIM];
  double cpp, cppx[DIM], cppxx[DIM][DIM];
  for(unsigned int i=0; i<nScalars; ++i){
    ComputeField<T,DIM>(SCALAR,i,p,U,dof,&c,&cx[0],&cxx[0][0]);
    ComputeField<PetscReal,DIM>(SCALAR,i,p,Up,dof,&cp,&cpx[0],&cpxx[0][0]);
    ComputeField<PetscReal,DIM>(SCALAR,i,p,Upp,dof,&cpp,&cppx[0],&cppxx[0][0]);
    solnS.defineVal(i,c);
    solnS.defineGrad(i,cx);
    solnS.defineHess(i,cxx);
    solnS.defineValP(i,cp);
    solnS.defineGradP(i,cpx);
    solnS.defineHessP(i,cpxx);
    solnS.defineValPP(i,cpp);
    solnS.defineGradPP(i,cppx);
    solnS.defineHessPP(i,cppxx);
  }

  //vector field variables
  T u[DIM], ux[DIM][DIM], uxx[DIM][DIM][DIM];
  double up[DIM], upx[DIM][DIM], upxx[DIM][DIM][DIM];
  double upp[DIM], uppx[DIM][DIM], uppxx[DIM][DIM][DIM];
  for(unsigned int i=0; i<nVectors; ++i){
    ComputeField<T,DIM>(VECTOR,nScalars+i*DIM,p,U,dof,&u[0],&ux[0][0],&uxx[0][0][0]);
    ComputeField<double,DIM>(VECTOR,nScalars+i*DIM,p,Up,dof,&up[0],&upx[0][0],&upxx[0][0][0]);
    ComputeField<double,DIM>(VECTOR,nScalars+i*DIM,p,Upp,dof,&upp[0],&uppx[0][0],&uppxx[0][0][0]);
    solnV.defineVal(i,u);
    solnV.defineGrad(i,ux);
    solnV.defineHess(i,uxx);
    solnV.defineValP(i,up);
    solnV.defineGradP(i,upx);
    solnV.defineHessP(i,upxx);
    solnV.defineValPP(i,upp);
    solnV.defineGradPP(i,uppx);
    solnV.defineHessPP(i,uppxx);
  }

  // Test Functions
  AD W[nen*dof];
  for(int i=0; i<nen*dof; i++){
    W[i]=1.;
    W[i].diff(i, dof*nen);
  } 
  testFunctionScalars<DIM,T> testFnS(nScalars);
  testFunctionVectors<DIM,T> testFnV(nVectors);
  AD residualObject;

  //scalar test functions
  AD w1, w1x[DIM], w1xx[DIM][DIM];
  for(unsigned int i=0; i<nScalars; ++i){
    ComputeField<AD,DIM>(SCALAR,i,p,W,dof,&w1,&w1x[0],&w1xx[0][0]);
    testFnS.defineVal(i,w1);
    testFnS.defineGrad(i,w1x);
    testFnS.defineHess(i,w1xx);
  }

  //vector test functions
  AD w2[DIM], w2x[DIM][DIM], w2xx[DIM][DIM][DIM];
  for(unsigned int i=0; i<nVectors; ++i){
    ComputeField<AD,DIM>(VECTOR,nScalars+i*DIM,p,W,dof,&w2[0],&w2x[0][0],&w2xx[0][0][0]);
    testFnV.defineVal(i,w2);
    testFnV.defineGrad(i,w2x);
    testFnV.defineHess(i,w2xx);
  }

  //compute residuals
  bool surfaceFlag=p->atboundary;
  residual<DIM,T>(!surfaceFlag,surfaceFlag,x,normal,solnS,solnV,testFnS,testFnV,*user,residualObject);

  // Loop over element nodes
  for (unsigned int i=0; i<dof*nen; i++) {
    R[i] = residualObject.fastAccessDx(i);
  }
  return 0;
}

template PetscErrorCode QuadPtResidual<PetscReal,1>(IGAPoint p,
						    const PetscReal * U,
						    const PetscScalar * Up,
						    const PetscScalar * Upp,
						    PetscReal *R,
						    void *ctx);
template PetscErrorCode QuadPtResidual<PetscReal,2>(IGAPoint p,
						    const PetscReal * U,
						    const PetscScalar * Up,
						    const PetscScalar * Upp,
						    PetscReal *R,
						    void *ctx);
template PetscErrorCode QuadPtResidual<PetscReal,3>(IGAPoint p,
						    const PetscReal * U,
						    const PetscScalar * Up,
						    const PetscScalar * Upp,
						    PetscReal *R,
						    void *ctx);
template PetscErrorCode QuadPtResidual<Sacado::Fad::SimpleFad<double>, 1 >(IGAPoint p,
									   const Sacado::Fad::SimpleFad<double> * U,
									   const PetscScalar * Up,
									   const PetscScalar * Upp,
									   Sacado::Fad::SimpleFad<double> *R,
									   void *ctx);
template PetscErrorCode QuadPtResidual<Sacado::Fad::SimpleFad<double>, 2 >(IGAPoint p,
									   const Sacado::Fad::SimpleFad<double> * U,
									   const PetscScalar * Up,
									   const PetscScalar * Upp,
									   Sacado::Fad::SimpleFad<double> *R,
									   void *ctx);
template PetscErrorCode QuadPtResidual<Sacado::Fad::SimpleFad<double>, 3 >(IGAPoint p,
									   const Sacado::Fad::SimpleFad<double> * U,
									   const PetscScalar * Up,
									   const PetscScalar * Upp,
									   Sacado::Fad::SimpleFad<double> *R,
									   void *ctx);
