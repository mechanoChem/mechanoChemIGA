#ifndef TESTFUNCTIONCLASS_H
#define TESTFUNCTIONCLASS_H

#include "tensor.h"

//include automatic differentiation library
#include <Sacado.hpp>
//#include <Sacado_Fad_BLAS.hpp>
#include <Sacado_Fad_SimpleFad.hpp>

/***************************************************************
 * testFunctionScalars class: contains value, gradient, and hessian of scalar
 * testFunction fields at a given quadrature point.
 */
template<unsigned int dim, typename T>
class testFunctionScalars{
public:

  /**
   * Constructor. Requires the number of scalar fields.
   */
  testFunctionScalars(unsigned int nScalars);

  /**
   * Copy constructor.
   */
  testFunctionScalars(const testFunctionScalars &obj);

  /**
   * Returns the value of the test function associated with scalar field i.
   */
  const Sacado::Fad::SimpleFad<T> val(unsigned int i) const { return value[i]; };

  /**
   * Returns the gradient of the test function associated with scalar field i.
   */
  const Tensor<1,dim,Sacado::Fad::SimpleFad<T> > grad(unsigned int i) const { return gradient[i]; };

  /**
   * Returns the Hessian of the test function associated with scalar field i.
   */
  const Tensor<2,dim,Sacado::Fad::SimpleFad<T> > hess(unsigned int i) const { return hessian[i]; };

  /**
   * Returns the Laplacian of the test function associated with scalar field i.
   */
  const Sacado::Fad::SimpleFad<T> laplacian(unsigned int i) const;
  
  //Input functions

  /**
   * Define the value of the test function associated with scalar field i.
   */
  void defineVal(unsigned int i, Sacado::Fad::SimpleFad<T> val){ value[i] = val; };

  /**
   * Define the gradient of the test function associated with scalar field i.
   */
  void defineGrad(unsigned int i, const Sacado::Fad::SimpleFad<T> grad[]){ gradient[i] = convert<dim,Sacado::Fad::SimpleFad<T> >(grad); };

  /**
   * Define the Hessian of the test function associated with scalar field i.
   */
  void defineHess(unsigned int i, const Sacado::Fad::SimpleFad<T> hess[][dim]){ hessian[i] = convert<dim,Sacado::Fad::SimpleFad<T> >(hess); };
			       
private:
  std::vector<Sacado::Fad::SimpleFad<T> > value;
  std::vector<Tensor<1,dim,Sacado::Fad::SimpleFad<T> > > gradient;
  std::vector<Tensor<2,dim,Sacado::Fad::SimpleFad<T> > > hessian;
  unsigned int size;

};

template<unsigned int dim, typename T>
testFunctionScalars<dim,T>::testFunctionScalars(unsigned int nScalars)
  :
  value(nScalars),
  gradient(nScalars),
  hessian(nScalars)
{

  size = nScalars;

}

template<unsigned int dim, typename T>
testFunctionScalars<dim,T>::testFunctionScalars(const testFunctionScalars &obj)
{

  value = obj.value;
  gradient = obj.gradient;
  hessian = obj.hessian;
  size = obj.size;

}

template<unsigned int dim, typename T>
  const Sacado::Fad::SimpleFad<T> testFunctionScalars<dim,T>::laplacian(unsigned int i) const {

  Sacado::Fad::SimpleFad<T> tmp = 0;
  for (unsigned int j=0; j<dim; ++j){
    tmp += hessian[i][j][j];
  }

  return tmp;

}

/***************************************************************
 * testFunctionVectors class: contains value, gradient, and hessian of vector
 * testFunction fields at a given quadrature point.
 */
template<unsigned int dim, typename T>
class testFunctionVectors{
public:

  /**
   * Constructor. Requires the number of scalar fields.
   */
  testFunctionVectors(unsigned int nVectors);

  /**
   * Copy constructor.
   */
  testFunctionVectors(const testFunctionVectors &obj);

  /**
   * Returns the value of the test function associated with vector field i.
   */
  const Tensor<1,dim,Sacado::Fad::SimpleFad<T> > val(unsigned int i) const { return value[i]; };

  /**
   * Returns the gradient of the test function associated with vector field i.
   */
  const Tensor<2,dim,Sacado::Fad::SimpleFad<T> > grad(unsigned int i) const { return gradient[i]; };

  /**
   * Returns the Hessian of the test function associated with vector field i.
   */
  const Tensor<3,dim,Sacado::Fad::SimpleFad<T> > hess(unsigned int i) const { return hessian[i]; };

  //Input functions

  /**
   * Define the value of the test function associated with vector field i.
   */
  void defineVal(unsigned int i, const Sacado::Fad::SimpleFad<T> val[]){ value[i] = convert<dim,Sacado::Fad::SimpleFad<T> >(val); };

  /**
   * Define the gradient of the test function associated with vector field i.
   */
  void defineGrad(unsigned int i, const Sacado::Fad::SimpleFad<T> grad[][dim]){ gradient[i] = convert<dim,Sacado::Fad::SimpleFad<T> >(grad); };

  /**
   * Define the Hessian of the test function associated with vector field i.
   */
  void defineHess(unsigned int i, const Sacado::Fad::SimpleFad<T> hess[][dim][dim]){ hessian[i] = convert<dim,Sacado::Fad::SimpleFad<T> >(hess); };
			       
private:
  std::vector<Tensor<1,dim,Sacado::Fad::SimpleFad<T> > > value;
  std::vector<Tensor<2,dim,Sacado::Fad::SimpleFad<T> > > gradient;
  std::vector<Tensor<3,dim,Sacado::Fad::SimpleFad<T> > > hessian;
  unsigned int size;

};

template<unsigned int dim, typename T>
testFunctionVectors<dim,T>::testFunctionVectors(unsigned int nVectors)
  :
  value(nVectors),
  gradient(nVectors),
  hessian(nVectors)
{

  size = nVectors;

}

template<unsigned int dim, typename T>
testFunctionVectors<dim,T>::testFunctionVectors(const testFunctionVectors &obj)
{

  value = obj.value;
  gradient = obj.gradient;
  hessian = obj.hessian;
  size = obj.size;

}

#endif //TESTFUNCTIONCLASS_H
