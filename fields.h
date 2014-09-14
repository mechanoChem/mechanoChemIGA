#ifndef FIELDS_H_
#define FIELDS_H_
#include <math.h>
extern "C" {
#include "petiga.h"
}

enum fieldType{SCALAR,VECTOR,TENSOR};

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

#endif /* FIELDS_H_ */

