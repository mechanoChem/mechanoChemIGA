//extern "C" {
#include "petiga.h"
//}
#include <math.h>

#include "coreFunctions.h"
//include automatic differentiation library
#include <Sacado.hpp>
//#include <Sacado_Fad_BLAS.hpp>
#include <Sacado_Fad_SimpleFad.hpp>

template <class T, unsigned int dim>
void ComputeField(fieldType type, unsigned int index, IGAPoint p, const T* U, unsigned int Udof, T* _value, T* _grad, T* _hess)
{
  unsigned int size=pow(dim,(int)type);
  //compute field value
  T (*U2)[Udof] = (T (*)[Udof])U;
  if (_value){
    PetscReal (*N) = (PetscReal (*)) p->shape[0];
    for (unsigned int i=0; i<size; i++) _value[i]=0.0;
    for(unsigned int n=0; n<(unsigned int) p->nen; n++){
      for(unsigned int d=0; d<size; d++) _value[d]+= N[n]*U2[n][index+d];
    }
  }

  //compute field gradient
  if (_grad){
    PetscReal (*N1)[dim] = (PetscReal (*)[dim]) p->shape[1];
    for (unsigned int i=0; i<size*dim; i++) _grad[i]=0.0;
    T (*grad)[dim] = (T (*)[dim])_grad;
    for(unsigned int n=0; n<(unsigned int) p->nen; n++){
      for(unsigned int d1=0; d1<size; d1++){
	for(unsigned int d2=0; d2<dim; d2++) grad[d1][d2]+= N1[n][d2]*U2[n][index+d1];
      }
    }
  }
  
  //compute field hessian
  if (_hess){
    PetscReal (*N2)[dim][dim] = (PetscReal (*)[dim][dim]) p->shape[2];
    for (unsigned int i=0; i<size*dim*dim; i++) _hess[i]=0.0;
    T (*hess)[dim][dim] = (T (*)[dim][dim])_hess;
    for(unsigned int n=0; n<(unsigned int) p->nen; n++){
      for(unsigned int d1=0; d1<size; d1++){
	for(unsigned int d2=0; d2<dim; d2++){
	  for(unsigned int d3=0; d3<dim; d3++)  hess[d1][d2][d3]+= N2[n][d2][d3]*U2[n][index+d1];
	}
      }
    }
  }
}

template void ComputeField<Sacado::Fad::SimpleFad<double>, 1>(fieldType type, unsigned int index, IGAPoint p, 
							 const Sacado::Fad::SimpleFad<double>* U, unsigned int Udof, 
							 Sacado::Fad::SimpleFad<double>* _value=0, 
							 Sacado::Fad::SimpleFad<double>* _grad=0, 
							 Sacado::Fad::SimpleFad<double>* _hess=0);
template void ComputeField<Sacado::Fad::SimpleFad<double>, 2>(fieldType type, unsigned int index, IGAPoint p, 
							 const Sacado::Fad::SimpleFad<double>* U, unsigned int Udof, 
							 Sacado::Fad::SimpleFad<double>* _value=0, 
							 Sacado::Fad::SimpleFad<double>* _grad=0, 
							 Sacado::Fad::SimpleFad<double>* _hess=0);
template void ComputeField<Sacado::Fad::SimpleFad<double>, 3>(fieldType type, unsigned int index, IGAPoint p, 
							 const Sacado::Fad::SimpleFad<double>* U, unsigned int Udof, 
							 Sacado::Fad::SimpleFad<double>* _value=0, 
							 Sacado::Fad::SimpleFad<double>* _grad=0, 
							 Sacado::Fad::SimpleFad<double>* _hess=0);
template void ComputeField<Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> >, 1>(fieldType type, unsigned int index, IGAPoint p, 
									     const Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> >* U, unsigned int Udof, 
									     Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> >* _value=0, 
									     Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> >* _grad=0, 
									     Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> >* _hess=0);
template void ComputeField<Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> >, 2>(fieldType type, unsigned int index, IGAPoint p, 
									     const Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> >* U, unsigned int Udof, 
									     Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> >* _value=0, 
									     Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> >* _grad=0, 
									     Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> >* _hess=0);
template void ComputeField<Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> >, 3>(fieldType type, unsigned int index, IGAPoint p, 
									     const Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> >* U, unsigned int Udof, 
									     Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> >* _value=0, 
									     Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> >* _grad=0, 
									     Sacado::Fad::SimpleFad<Sacado::Fad::SimpleFad<double> >* _hess=0);

template void ComputeField<PetscReal, 1>(fieldType type, unsigned int index, IGAPoint p, const PetscReal* U, unsigned int Udof, PetscReal* _value=0, PetscReal* _grad=0, PetscReal* _hess=0);
template void ComputeField<PetscReal, 2>(fieldType type, unsigned int index, IGAPoint p, const PetscReal* U, unsigned int Udof, PetscReal* _value=0, PetscReal* _grad=0, PetscReal* _hess=0);
template void ComputeField<PetscReal, 3>(fieldType type, unsigned int index, IGAPoint p, const PetscReal* U, unsigned int Udof, PetscReal* _value=0, PetscReal* _grad=0, PetscReal* _hess=0);
