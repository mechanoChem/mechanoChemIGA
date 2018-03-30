#ifndef SOLUTIONCLASS_H
#define SOLUTIONCLASS_H

#include "tensor.h"

/***************************************************************
 * solutionScalars class: contains value, gradient, and hessian of scalar
 * solution fields at a given quadrature point.
 */
template<unsigned int dim, typename T>
class solutionScalars{
public:

  /**
   * Constructor. Requires the number of scalar fields.
   */
  solutionScalars(unsigned int nScalars);

  /**
   * Copy constructor.
   */
  solutionScalars(const solutionScalars &obj);
  
  /**
   * c.val(i) - Returns the value of scalar field i (scalar).
   */
  const T val(unsigned int i) const { return value[i]; };

  /**
   * c.grad(i) - Returns the gradient of scalar field i (1st order tensor).
   */
  const Tensor<1,dim,T> grad(unsigned int i) const { return gradient[i]; };

  /**
   * c.hess(i) - Returns the hessian of scalar field i (2nd order tensor).
   */
  const Tensor<2,dim,T> hess(unsigned int i) const { return hessian[i]; };

  /**
   * c.laplacian(i) - Returns the Laplacian of scalar field i (scalar).
   */
  const T laplacian(unsigned int i) const;
  
  /**
   * c.valP(i) - Returns the value of scalar field i at the previous time step (scalar).
   */
  const double valP(unsigned int i) const { return valueP[i]; };

  /**
   * c.gradP(i) - Returns the gradient of scalar field i at the previous time step (1st order tensor).
   */
  const Tensor<1,dim,double> gradP(unsigned int i) const { return gradientP[i]; };

  /**
   * c.hessP(i) - Returns the hessian of scalar field i at the previous time step (2nd order tensor).
   */
  const Tensor<2,dim,double> hessP(unsigned int i) const { return hessianP[i]; };

  /**
   * c.laplacianP(i) - Returns the Laplacian of scalar field i at the previous time step (scalar).
   */
  const double laplacianP(unsigned int i) const;
  
  /**
   * c.valPP(i) - Returns the value of scalar field i at the step before the previous time step (scalar).
   */
  const double valPP(unsigned int i) const { return valuePP[i]; };

  /**
   * c.gradPP(i) - Returns the gradient of scalar field i at the step before the previous time step (1st order tensor).
   */
  const Tensor<1,dim,double> gradPP(unsigned int i) const { return gradientPP[i]; };

  /**
   * c.hessPP(i) - Returns the hessian of scalar field i at the step before the previous time step (2nd order tensor).
   */
  const Tensor<2,dim,double> hessPP(unsigned int i) const { return hessianPP[i]; };

  /**
   * c.laplacianPP(i) - Returns the Laplacian of scalar field i at the step before the previous time step (scalar).
   */
  const double laplacianPP(unsigned int i) const;
  
  //Input functions

  /**
   * Define the current value of scalar field i.
   */
  void defineVal(unsigned int i, T val){ value[i] = val; };

  /**
   * Define the current gradient of scalar field i with a 1D array.
   */
  void defineGrad(unsigned int i, const T grad[]){ gradient[i] = convert<dim,T>(grad); };

  /**
   * Define the current Hessian of scalar field i with a 2D array.
   */
  void defineHess(unsigned int i, const T hess[][dim]){ hessian[i] = convert<dim,T>(hess); };

  /**
   * Define the previous value of scalar field i.
   */
  void defineValP(unsigned int i, double valP){ valueP[i] = valP; };

  /**
   * Define the previous gradient of scalar field i with a 1D array.
   */
  void defineGradP(unsigned int i, const double gradP[]){ gradientP[i] = convert<dim,double>(gradP); };

  /**
   * Define the previous Hessian of scalar field i with a 2D array.
   */
  void defineHessP(unsigned int i, const double hessP[][dim]){ hessianP[i] = convert<dim,double>(hessP); };

  /**
   * Define the previous previous value of scalar field i.
   */
  void defineValPP(unsigned int i, double valPP){ valuePP[i] = valPP; };

  /**
   * Define the previous previous gradient of scalar field i with a 1D array.
   */
  void defineGradPP(unsigned int i, const double gradPP[]){ gradientPP[i] = convert<dim,double>(gradPP); };

  /**
   * Define the previous previous Hessian of scalar field i with a 2D array.
   */
  void defineHessPP(unsigned int i, const double hessPP[][dim]){ hessianPP[i] = convert<dim,double>(hessPP); };
			       
private:
  std::vector<T> value;
  std::vector<Tensor<1,dim,T> > gradient;
  std::vector<Tensor<2,dim,T> > hessian;
  std::vector<double> valueP;
  std::vector<Tensor<1,dim,double> > gradientP;
  std::vector<Tensor<2,dim,double> > hessianP;
  std::vector<double> valuePP;
  std::vector<Tensor<1,dim,double> > gradientPP;
  std::vector<Tensor<2,dim,double> > hessianPP;
  unsigned int size;

};

template<unsigned int dim, typename T>
solutionScalars<dim,T>::solutionScalars(unsigned int nScalars)
  :
  value(nScalars),
  gradient(nScalars),
  hessian(nScalars),
  valueP(nScalars),
  gradientP(nScalars),
  hessianP(nScalars),
  valuePP(nScalars),
  gradientPP(nScalars),
  hessianPP(nScalars)
{

  size = nScalars;

}

template<unsigned int dim, typename T>
solutionScalars<dim,T>::solutionScalars(const solutionScalars &obj)
{

  value = obj.value;
  gradient = obj.gradient;
  hessian = obj.hessian;
  valueP = obj.valueP;
  gradientP = obj.gradientP;
  hessianP = obj.hessianP;
  valuePP = obj.valuePP;
  gradientPP = obj.gradientPP;
  hessianPP = obj.hessianPP;
  size = obj.size;

}

template<unsigned int dim, typename T>
  const T solutionScalars<dim,T>::laplacian(unsigned int i) const {

  T tmp = 0;
  for (unsigned int j=0; j<dim; ++j){
    tmp += hessian[i][j][j];
  }

  return tmp;

}

template<unsigned int dim, typename T>
  const double solutionScalars<dim,T>::laplacianP(unsigned int i) const {

  double tmp = 0;
  for (unsigned int j=0; j<dim; ++j){
    tmp += hessianP[i][j][j];
  }

  return tmp;

}

template<unsigned int dim, typename T>
  const double solutionScalars<dim,T>::laplacianPP(unsigned int i) const {

  double tmp = 0;
  for (unsigned int j=0; j<dim; ++j){
    tmp += hessianPP[i][j][j];
  }

  return tmp;

}

/***************************************************************
 * solutionVectors class: contains value, gradient, and hessian of vector
 * solution fields at a given quadrature point.
 */
template<unsigned int dim, typename T>
class solutionVectors{
public:

  /**
   * Constructor. Requires the number of vector fields.
   */
  solutionVectors(unsigned int nVectors);

  /**
   * Copy constructor.
   */
  solutionVectors(const solutionVectors &obj);

  /**
   * u.val(i) - Returns the value of vector field i (1st order tensor).
   */
  const Tensor<1,dim,T> val(unsigned int i) const { return value[i]; };

  /**
   * u.grad(i) - Returns the gradient of vector field i (2nd order tensor).
   */
  const Tensor<2,dim,T> grad(unsigned int i) const { return gradient[i]; };

  /**
   * u.hess(i) - Returns the Hessian of vector field i (3rd order tensor).
   */
  const Tensor<3,dim,T> hess(unsigned int i) const { return hessian[i]; };

  /**
   * u.valP(i) - Returns the value of vector field i at the previous time step (1st order tensor).
   */
  const Tensor<1,dim,double> valP(unsigned int i) const { return valueP[i]; };

  /**
   * u.gradP(i) - Returns the gradient of vector field i at the previous time step (2nd order tensor).
   */
  const Tensor<2,dim,double> gradP(unsigned int i) const { return gradientP[i]; };

  /**
   * u.hessP(i) - Returns the Hessian of vector field i at the previous time step (3rd order tensor).
   */
  const Tensor<3,dim,double> hessP(unsigned int i) const { return hessianP[i]; };

  /**
   * u.valPP(i) - Returns the value of vector field i at the step before the previous time step (1st order tensor).
   */
  const Tensor<1,dim,double> valPP(unsigned int i) const { return valuePP[i]; };

  /**
   * u.gradPP(i) - Returns the gradient of vector field i at the step before the previous time step (2nd order tensor).
   */
  const Tensor<2,dim,double> gradPP(unsigned int i) const { return gradientPP[i]; };

  /**
   * u.hessPP(i) - Returns the Hessian of vector field i at the step before the previous time step (3rd order tensor).
   */
  const Tensor<3,dim,double> hessPP(unsigned int i) const { return hessianPP[i]; };

  //Input functions

  /**
   * Define the current value of vector field i with a 1D array.
   */
  void defineVal(unsigned int i, const T val[]){ value[i] = convert<dim,T>(val); };

  /**
   * Define the current gradient of vector field i with a 2D array.
   */
  void defineGrad(unsigned int i, const T grad[][dim]){ gradient[i] = convert<dim,T>(grad); };

  /**
   * Define the current Hessian of vector field i with a 3D array.
   */
  void defineHess(unsigned int i, const T hess[][dim][dim]){ hessian[i] = convert<dim,T>(hess); };

  /**
   * Define the previous value of vector field i with a 1D array.
   */
  void defineValP(unsigned int i, const double valP[]){ valueP[i] = convert<dim,double>(valP); };

  /**
   * Define the previous gradient of vector field i with a 2D array.
   */
  void defineGradP(unsigned int i, const double gradP[][dim]){ gradientP[i] = convert<dim,double>(gradP); };

  /**
   * Define the previous Hessian of vector field i with a 3D array.
   */
  void defineHessP(unsigned int i, const double hessP[][dim][dim]){ hessianP[i] = convert<dim,double>(hessP); };

  /**
   * Define the previous previous value of vector field i with a 1D array.
   */
  void defineValPP(unsigned int i, const double valPP[]){ valuePP[i] = convert<dim,double>(valPP); };

  /**
   * Define the previous previous gradient of vector field i with a 2D array.
   */
  void defineGradPP(unsigned int i, const double gradPP[][dim]){ gradientPP[i] = convert<dim,double>(gradPP); };

  /**
   * Define the previous previous Hessian of vector field i with a 3D array.
   */
  void defineHessPP(unsigned int i, const double hessPP[][dim][dim]){ hessianPP[i] = convert<dim,double>(hessPP); };
    
private:
  std::vector<Tensor<1,dim,T> > value;
  std::vector<Tensor<2,dim,T> > gradient;
  std::vector<Tensor<3,dim,T> > hessian;
  std::vector<Tensor<1,dim,double> > valueP;
  std::vector<Tensor<2,dim,double> > gradientP;
  std::vector<Tensor<3,dim,double> > hessianP;
  std::vector<Tensor<1,dim,double> > valuePP;
  std::vector<Tensor<2,dim,double> > gradientPP;
  std::vector<Tensor<3,dim,double> > hessianPP;
  unsigned int size;

};

template<unsigned int dim, typename T>
solutionVectors<dim,T>::solutionVectors(unsigned int nVectors)
  :
  value(nVectors),
  gradient(nVectors),
  hessian(nVectors),
  valueP(nVectors),
  gradientP(nVectors),
  hessianP(nVectors),
  valuePP(nVectors),
  gradientPP(nVectors),
  hessianPP(nVectors)

{

  size = nVectors;

}

template<unsigned int dim, typename T>
solutionVectors<dim,T>::solutionVectors(const solutionVectors &obj)
{

  value = obj.value;
  gradient = obj.gradient;
  hessian = obj.hessian;
  valueP = obj.valueP;
  gradientP = obj.gradientP;
  hessianP = obj.hessianP;
  valuePP = obj.valuePP;
  gradientPP = obj.gradientPP;
  hessianPP = obj.hessianPP;
  size = obj.size;

}

#endif //SOLUTIONCLASS_H
