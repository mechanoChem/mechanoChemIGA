#ifndef TENSOR_H
#define TENSOR_H

#include <vector>
//#include <math.h>
#include <cmath>

//Forward declaration
template<unsigned int order,unsigned int dim,typename T>
  class Tensor;


/**
 * First-order Tensor class
 */
template<unsigned int dim,typename T>
class Tensor <1,dim,T> {
public:
  Tensor();
  Tensor(const Tensor<1,dim,T> &b);

/**
 * Access/set an element of the tensor.
 */
  T& operator[](unsigned int i){ return value[i]; };

/**
 * Access an element of the tensor.
 */
  const T& operator[](unsigned int i) const { return value[i]; };
  Tensor<1,dim,T>& operator=(const Tensor<1,dim,T> &b);
  T norm() const;
  
  std::vector<T> value;
};

/**
 * Default constructor.
 */
template<unsigned int dim,typename T>
Tensor<1,dim,T>::Tensor() : value(dim){}

/**
 * Copy constructor.
 */
template<unsigned int dim,typename T>
Tensor<1,dim,T>::Tensor(const Tensor<1,dim,T> &b){

  *this = b;
  
}

/**
 * Assignment operator.
 */
template<unsigned int dim,typename T>
Tensor<1,dim,T>& Tensor<1,dim,T>::operator=(const Tensor<1,dim,T> &b){
  
  value = b.value;
  
  return *this;
}

/**
 * Compute the L2-norm of a 1st order tensor.
 */
template<unsigned int dim,typename T>
T Tensor<1,dim,T>::norm() const{
  
  T tmp = 0;
  for (unsigned int i=0; i<dim; ++i){
tmp += std::pow(value[i],2);
  }

return std::sqrt(tmp);
  
}

/**
 * Primary Tensor template class
 * @tparam order An unsigned int defining the order of the tensor (1, 2, 3, or 4).
 * @tparam dim An unsigned int defing the spatial dimension of the problem.
 * @tparam T The type that is to be stored in the tensor.
 */
template<unsigned int order,unsigned int dim,typename T>
  class Tensor{
public:
  Tensor();
  Tensor(const Tensor<order,dim,T> &b);

/**
 * Access/set an element of the tensor. Returns a tensor with the order reduced by 1.
 */
  Tensor<order-1,dim,T>& operator[](unsigned int i){ return value[i]; };

/**
 * Access an element of the tensor. Returns a tensor with the order reduced by 1.
 */
  const Tensor<order-1,dim,T>& operator[](unsigned int i) const { return value[i]; };
  Tensor<order,dim,T>& operator=(const Tensor<order,dim,T> &b);
  template<typename U> auto operator*(const Tensor<1,dim,U> &b) -> Tensor<order-1,dim,decltype(std::declval<T>()*std::declval<U>())> const;
  
  std::vector<Tensor<order-1,dim,T> > value;
};

/**
 * Default constructor.
 */
template<unsigned int order,unsigned int dim,typename T>
Tensor<order,dim,T>::Tensor() : value(dim){}

/**
 * Copy constructor.
 */
template<unsigned int order,unsigned int dim,typename T>
Tensor<order,dim,T>::Tensor(const Tensor<order,dim,T> &b){

  *this = b;
  
}

/**
 * Assignment operator.
 */
template<unsigned int order,unsigned int dim,typename T>
Tensor<order,dim,T>& Tensor<order,dim,T>::operator=(const Tensor<order,dim,T> &b){
  
  value = b.value;
  
  return *this;
}

/**
 * Contract the given tensor with a 1st order tensor.
 */
template<unsigned int order,unsigned int dim,typename T>
template<typename U>
  auto Tensor<order,dim,T>::operator*(const Tensor<1,dim,U> &b) -> Tensor<order-1,dim,decltype(std::declval<T>()*std::declval<U>())> const{
  
  Tensor<order-1,dim,decltype(std::declval<T>()*std::declval<U>())> tmp;
  for (unsigned int i=0; i<dim; ++i){
    tmp[i] = this->value[i]*b;
  }
  
  return tmp;
  
}

/**
 * @name Operations related to all Tensor objects:
 */

/**
 * @relates Tensor
 *
 * Scalar multiplication of a Tensor, where the scalar precedes the Tensor object.
 */
template<unsigned int order,unsigned int dim,typename T,typename U> 
  auto operator*(T a, const Tensor<order,dim,U> &b) -> Tensor<order,dim,decltype(std::declval<T>()*std::declval<U>())>{
  
  Tensor<order,dim,decltype(std::declval<T>()*std::declval<U>())> tmp;
  for (unsigned int i=0; i<dim; ++i){
      tmp[i] = a*b[i];
  }
  
  return tmp;
  
}

/**
 * @relates Tensor
 *
 * Scalar multiplication of a Tensor, where the scalar follows the Tensor object.
 */
template<unsigned int order,unsigned int dim,typename T,typename U> 
  auto operator*(const Tensor<order,dim,U> &b,T a) -> Tensor<order,dim,decltype(std::declval<T>()*std::declval<U>())>{
  
  Tensor<order,dim,decltype(std::declval<T>()*std::declval<U>())> tmp;
  for (unsigned int i=0; i<dim; ++i){
      tmp[i] = a*b[i];
  }
  
  return tmp;
  
}

/**
 * @relates Tensor
 *
 * Addition operator for any two Tensors of the same order.
 */
template<unsigned int order,unsigned int dim,typename T,typename U> 
  auto operator+(const Tensor<order,dim,T> &a, const Tensor<order,dim,U> &b) -> Tensor<order,dim,decltype(std::declval<T>()+std::declval<U>())>{
  
  Tensor<order,dim,decltype(std::declval<T>()+std::declval<U>())> tmp;
  for (unsigned int i=0; i<dim; ++i){
      tmp[i] = a[i]+b[i];
  }
  
  return tmp;
  
}

/**
 * @relates Tensor
 *
 * Subtraction operator for any two Tensors of the same order.
 */
template<unsigned int order,unsigned int dim,typename T,typename U> 
  auto operator-(const Tensor<order,dim,T> &a, const Tensor<order,dim,U> &b) -> Tensor<order,dim,decltype(std::declval<T>()-std::declval<U>())>{
  
  Tensor<order,dim,decltype(std::declval<T>()-std::declval<U>())> tmp;
  for (unsigned int i=0; i<dim; ++i){
      tmp[i] = a[i]-b[i];
  }
  
  return tmp;
  
}

/**
 * @name Operations related to 1st order Tensor objects:
 */

/**
 * @relates Tensor
 *
 * Single contraction (i.e. dot product) of two 1st order Tensor objects. Returns a scalar.
 */
template<unsigned int dim,typename T,typename U> 
auto operator*(const Tensor<1,dim,T> &a, const Tensor<1,dim,U> &b) -> decltype(a[0]*b[0]){

  decltype(a[0]*b[0]) tmp(0);
  for (unsigned int i=0; i<dim; ++i){
    tmp += a[i]*b[i];
  }

  return tmp;
  
}

/**
 * @name Operations related to 2nd order Tensor objects:
 */

/**
 * @relates Tensor
 *
 * Double contraction of two 2nd order Tensor objects. Returns a scalar.
 */
template<unsigned int dim,typename T,typename U>
  auto double_contract(const Tensor<2,dim,T> &a, const Tensor<2,dim,U> &b) -> decltype(std::declval<T>()*std::declval<U>()){
  
  decltype(std::declval<T>()*std::declval<U>()) tmp = 0;
  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      tmp += a[i][j]*b[i][j];
    }
  }
  
  return tmp;
  
}

/**
 * @relates Tensor
 *
 * Multiplication operator of two 2nd order Tensor objects (equivalent to matrix multiplication). Returns a 2nd order Tensor object.
 */
template<unsigned int dim,typename T,typename U>
  auto operator*(const Tensor<2,dim,T> a, const Tensor<2,dim,U> &b) -> Tensor<2,dim,decltype(std::declval<T>()*std::declval<U>())>{
  
  Tensor<2,dim,decltype(std::declval<T>()*std::declval<U>())> tmp;
  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      for (unsigned int k=0; k<dim; ++k){
	tmp[i][k] += a[i][j]*b[j][k];
      }
    }
  }
  
  return tmp;
  
}

/**
 * @relates Tensor
 *
 * Converts the 2nd order Tensor object to an isotropic tensor (equivalent to the identity matrix).
 */
template<unsigned int dim,typename T>
  void identity(Tensor<2,dim,T>& a) {
  
  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      a[i][j] = (i == j);
    }
  }
}

/**
 * @relates Tensor
 *
 * Returns the transpose of the 2nd order Tensor object.
 */
template<unsigned int dim,typename T>
Tensor<2,dim,T> trans(const Tensor<2,dim,T>& a) {

  Tensor<2,dim,T> tmp;
  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      tmp[i][j] = a[j][i];
    }
  }

  return tmp;
}

/**
 * @relates Tensor
 *
 * Returns the trace of the 2nd order Tensor object.
 */
template<unsigned int dim,typename T>
T trace(const Tensor<2,dim,T>& a) {

  T tmp(0);
  for (unsigned int i=0; i<dim; ++i){
    tmp += a[i][i];
  }

  return tmp;
}

/**
 * @relates Tensor
 *
 * Returns the determinant of the 2nd order Tensor object (dim==2).
 */
template<typename T>
T det(const Tensor<2,2,T>& a) {
  
  return (a[0][0]*a[1][1] - a[0][1]*a[1][0]);
}

/**
 * @relates Tensor
 *
 * Returns the determinant of the 2nd order Tensor object (dim==3).
 */
template<typename T>
T det(const Tensor<2,3,T>& a) {
  
  return (a[0][0]*(a[1][1]*a[2][2] - a[1][2]*a[2][1])
	  - a[0][1]*(a[1][0]*a[2][2] - a[1][2]*a[2][0])
	  + a[0][2]*(a[1][0]*a[2][1] - a[1][1]*a[2][0]));
}

/**
 * @relates Tensor
 *
 * Returns the inverse of the 2nd order Tensor object (dim==2).
 */
template<typename T>
Tensor<2,2,T> inv(const Tensor<2,2,T>& a) {
    
  T J = det(a);
  Tensor<2,2,T> tmp;
  tmp[0][0] = 1./J*a[1][1];
  tmp[0][1] = -1./J*a[0][1];
  tmp[1][0] = -1./J*a[1][0];
  tmp[1][1] = 1./J*a[0][0];

  return tmp;
}

/**
 * @relates Tensor
 *
 * Returns the inverse of the 2nd order Tensor object (dim==3).
 */
template<typename T>
Tensor<2,3,T> inv(const Tensor<2,3,T>& a) {
    
  T J = det(a);
  Tensor<2,3,T> tmp;
  tmp[0][0] = 1./J*(a[1][1]*a[2][2] - a[2][1]*a[1][2]);
  tmp[0][1] = 1./J*(a[0][2]*a[2][1] - a[0][1]*a[2][2]);
  tmp[0][2] = 1./J*(a[0][1]*a[1][2] - a[1][1]*a[0][2]);
  tmp[1][0] = 1./J*(a[1][2]*a[2][0] - a[2][2]*a[1][0]);
  tmp[1][1] = 1./J*(a[0][0]*a[2][2] - a[2][0]*a[0][2]);
  tmp[1][2] = 1./J*(a[0][2]*a[1][0] - a[1][2]*a[0][0]);
  tmp[2][0] = 1./J*(a[1][0]*a[2][1] - a[2][0]*a[1][1]);
  tmp[2][1] = 1./J*(a[0][1]*a[2][0] - a[2][1]*a[0][0]);
  tmp[2][2] = 1./J*(a[0][0]*a[1][1] - a[1][0]*a[0][1]);

  return tmp;
}

/**
 * @name Operations related to 3rd order Tensor objects:
 */

/**
 * @relates Tensor
 *
 * Triple contraction of two 3rd order Tensor objects. Returns a scalar.
 */
template<unsigned int dim,typename T,typename U>
  auto triple_contract(const Tensor<3,dim,T> &a, const Tensor<3,dim,U> &b) -> decltype(std::declval<T>()*std::declval<U>()){
  
  decltype(std::declval<T>()*std::declval<U>()) tmp = 0;
  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      for (unsigned int k=0; k<dim; ++k){
	tmp += a[i][j][k]*b[i][j][k];
      }
    }
  }
  
  return tmp;
  
}

/**
 * @name Operations related to 4th order Tensor objects:
 */

/**
 * @relates Tensor
 *
 * Performs a full contraction on two 2nd order Tensors, \f$a\f$ and \f$c\f$, and a 4th order Tensor, \f$b\f$, in this manner: \f$a_{ij}b_{ijkl}c_{kl}\f$,
 * with implied summation over all indices. Returns a scalar.
 */
template<unsigned int dim,typename T,typename U,typename V>
auto full_contract(const Tensor<2,dim,T> &a, const Tensor<4,dim,U> &b, const Tensor<2,dim,V> &c) -> decltype(std::declval<T>()*std::declval<U>()*std::declval<V>()){
  
  decltype(std::declval<T>()*std::declval<U>()*std::declval<V>()) tmp = 0;
  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      for (unsigned int k=0; k<dim; ++k){
	for (unsigned int l=0; l<dim; ++l){
	  tmp += a[i][j]*b[i][j][k][l]*c[k][l];
	}
      }
    }
  }
  
  return tmp;
  
}

/**
 * @relates Tensor
 *
 * Performs a double contraction between a 4th order Tensor, \f$a\f$, and a 2nd order Tensor, \f$b\f$, in this manner: \f$a_{ijkl}b_{kl}\f$,
 * with implied summation over \f$k\f$ and \f$l\f$.
 */
template<unsigned int dim,typename T,typename U>
auto double_contract(const Tensor<4,dim,T> &a, const Tensor<2,dim,U> &b) -> Tensor<2,dim,decltype(std::declval<T>()*std::declval<U>())>{
  
  Tensor<2,dim,decltype(std::declval<T>()*std::declval<U>())> tmp;
  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      tmp[i][j] = 0;
      for (unsigned int k=0; k<dim; ++k){
	for (unsigned int l=0; l<dim; ++l){
	  tmp[i][j] += a[i][j][k][l]*b[k][l];
	}
      }
    }
  }
  
  return tmp;
  
}

/**
 * @name Operations converting arrays to Tensor objects:
 */

/**
 * @relates Tensor
 * 
 * Converts 1D array to a 1st order Tensor
 */
template<unsigned int dim,typename T>
  Tensor<1,dim,T> convert(const T b[]){
  
  Tensor<1,dim,T> tmp;
  tmp.value.assign(b,b+dim);
  
  return tmp;
}

/**
 * @relates Tensor
 * 
 * Converts 2D array to a 2nd order Tensor
 */
template<unsigned int dim,typename T>
  Tensor<2,dim,T> convert(const T b[][dim]){
  
  Tensor<2,dim,T> tmp;
  for(unsigned int i=0; i<dim; ++i){
    tmp.value[i] = convert<dim,T>(b[i]);
  }
  
  return tmp;
}

/**
 * @relates Tensor
 * 
 * Converts 3D array to a 3rd order Tensor
 */
template<unsigned int dim,typename T>
  Tensor<3,dim,T> convert(const T b[][dim][dim]){
  
  Tensor<3,dim,T> tmp;
  for(unsigned int i=0; i<dim; ++i){
    tmp.value[i] = convert<dim,T>(b[i]);
  }
  
  return tmp;
}

/**
 * @relates Tensor
 * 
 * Converts 4D array to a 4th order Tensor
 */
template<unsigned int dim,typename T>
  Tensor<4,dim,T> convert(const T b[][dim][dim][dim]){
  
  Tensor<4,dim,T> tmp;
  for(unsigned int i=0; i<dim; ++i){
    tmp.value[i] = convert<dim,T>(b[i]);
  }
  
  return tmp;
}


#endif //TENSOR_H
