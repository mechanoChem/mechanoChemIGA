#include "genericHeaders.h"
//extern "C" {
#include "petiga.h"
//}
#include <math.h>
//#include "../applications/configurationalForces/3D/parameters.h"

//include automatic differentiation library
#include <Sacado.hpp>

template <class T, unsigned int dim, unsigned int dof>
void computeField(fieldType type, unsigned int index, IGAPoint p, const T* U, T* _value=0, T* _grad=0, T* _hess=0)
{
  unsigned int size=pow(dim,(int)type);
  //compute field value
  T (*U2)[dof] = (T (*)[dof])U;
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

template void computeField<Sacado::Fad::SFad<double,36>, 2, 4>(fieldType type, unsigned int index, IGAPoint p, 
	const Sacado::Fad::SFad<double,36>* U, 
	Sacado::Fad::SFad<double,36>* _value=0, 
	Sacado::Fad::SFad<double,36>* _grad=0, 
	Sacado::Fad::SFad<double,36>* _hess=0);
template void computeField<Sacado::Fad::SFad<double,162>, 3, 6>(fieldType type, unsigned int index, IGAPoint p, 
	const Sacado::Fad::SFad<double,162>* U, 
	Sacado::Fad::SFad<double,162>* _value=0, 
	Sacado::Fad::SFad<double,162>* _grad=0, 
	Sacado::Fad::SFad<double,162>* _hess=0);
template void computeField<Sacado::Fad::SFad<double,27>, 2, 3>(fieldType type, unsigned int index, IGAPoint p, 
	const Sacado::Fad::SFad<double,27>* U, 
	Sacado::Fad::SFad<double,27>* _value=0, 
	Sacado::Fad::SFad<double,27>* _grad=0, 
	Sacado::Fad::SFad<double,27>* _hess=0);
template void computeField<Sacado::Fad::SFad<double,108>, 3, 4>(fieldType type, unsigned int index, IGAPoint p, 
	const Sacado::Fad::SFad<double,108>* U, 
	Sacado::Fad::SFad<double,108>* _value=0, 
	Sacado::Fad::SFad<double,108>* _grad=0, 
	Sacado::Fad::SFad<double,108>* _hess=0);
template void computeField<PetscReal, 2, 4>(fieldType type, unsigned int index, IGAPoint p, const PetscReal* U, PetscReal* _value=0, PetscReal* _grad=0, PetscReal* _hess=0);
template void computeField<PetscReal, 3, 6>(fieldType type, unsigned int index, IGAPoint p, const PetscReal* U, PetscReal* _value=0, PetscReal* _grad=0, PetscReal* _hess=0);
template void computeField<PetscReal, 2, 3>(fieldType type, unsigned int index, IGAPoint p, const PetscReal* U, PetscReal* _value=0, PetscReal* _grad=0, PetscReal* _hess=0);
template void computeField<PetscReal, 3, 4>(fieldType type, unsigned int index, IGAPoint p, const PetscReal* U, PetscReal* _value=0, PetscReal* _grad=0, PetscReal* _hess=0);
