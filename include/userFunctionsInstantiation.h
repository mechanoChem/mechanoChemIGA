#ifndef userFunctionsInstantiation_h
#define userFunctionsInstantiation_h

#include "tensor.h"
#include "solutionClass.h"
#include "testFunctionClass.h"
#include "appCtx.h"

//include automatic differentiation library
#include <Sacado.hpp>
#include <Sacado_Fad_SimpleFad.hpp>

template void defineParameters<2>(AppCtx<2>& user);
template void defineParameters<3>(AppCtx<3>& user);

template void residual<2,double>(bool dV,
				 bool dS,
				 const Tensor<1,2,double> &x,
				 const Tensor<1,2,double> &normal,
				 const solutionScalars<2,double> &c,
				 const solutionVectors<2,double> &u,
				 const testFunctionScalars<2,double> &w1,
				 const testFunctionVectors<2,double> &w2,
				 AppCtx<2> &user,
				 Sacado::Fad::SimpleFad<double> &r);
template void residual<2,Sacado::Fad::SimpleFad<double> >(bool dV,
				 bool dS,
				 const Tensor<1,2,double> &x,
				 const Tensor<1,2,double> &normal,
				 const solutionScalars<2,Sacado::Fad::SimpleFad<double> > &c,
				 const solutionVectors<2,Sacado::Fad::SimpleFad<double> > &u,
				 const testFunctionScalars<2,Sacado::Fad::SimpleFad<double> > &w1,
				 const testFunctionVectors<2,Sacado::Fad::SimpleFad<double> > &w2,
				 AppCtx<2> &user,
				 Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> > &r);
template void residual<3,double>(bool dV,
				 bool dS,
				 const Tensor<1,3,double> &x,
				 const Tensor<1,3,double> &normal,
				 const solutionScalars<3,double> &c,
				 const solutionVectors<3,double> &u,
				 const testFunctionScalars<3,double> &w1,
				 const testFunctionVectors<3,double> &w2,
				 AppCtx<3> &user,
				 Sacado::Fad::SimpleFad<double> &r);
template void residual<3,Sacado::Fad::SimpleFad<double> >(bool dV,
				 bool dS,
				 const Tensor<1,3,double> &x,
				 const Tensor<1,3,double> &normal,
				 const solutionScalars<3,Sacado::Fad::SimpleFad<double> > &c,
				 const solutionVectors<3,Sacado::Fad::SimpleFad<double> > &u,
				 const testFunctionScalars<3,Sacado::Fad::SimpleFad<double> > &w1,
				 const testFunctionVectors<3,Sacado::Fad::SimpleFad<double> > &w2,
				 AppCtx<3> &user,
				 Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> > &r);

#endif //userFunctionsInstantiation_h
