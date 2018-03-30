#include "tensor.h"
#include <iostream>

/**
 * @defgroup tests Unit tests
 * This is a collection of unit test functions.
 */

/**
 * @ingroup tests
 * Check if two doubles are equal, within a given tolerance.
 */
bool almost_equal(double a, double b, double eps){
  return (a == b || std::abs(a - b)<eps);
}

/**
 * @ingroup tests
 * Check if two doubles are (almost) equal, then return pass (true) or fail (false).
 */
bool test(double C, double Ctest){

  bool pass = true;
  if (!almost_equal(C,Ctest,1.e-12)){
    std::cout << C << " " << Ctest << " " << C-Ctest << std::endl;
    pass = false;
  }
  
  if(pass){
    std::cout << "Passed.\n";
  }
  else{
    std::cout << "Failed.\n";
  }
  return pass;
}

/**
 * @ingroup tests
 * Check if two 1st order tensors are (almost) equal, then return pass (true) or fail (false).
 */
template<unsigned int dim>
bool test(Tensor<1,dim,double> C, Tensor<1,dim,double> Ctest){

  bool pass = true;
  for (unsigned int i=0; i<dim; ++i){
    if (!almost_equal(C[i],Ctest[i],1.e-12)){
      std::cout << C[i] << " " << Ctest[i] << std::endl;
      pass = false;
    }
  }
  
  if(pass){
    std::cout << "Passed.\n";
  }
  else{
    std::cout << "Failed.\n";
  }
  return pass;
}

/**
 * @ingroup tests
 * Check if two 2nd order tensors are (almost) equal, then return pass (true) or fail (false).
 */
template<unsigned int dim>
bool test(Tensor<2,dim,double> C, Tensor<2,dim,double> Ctest){

  bool pass = true;
  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      if (!almost_equal(C[i][j],Ctest[i][j],1.e-12)){
	std::cout << C[i][j] << " " << Ctest[i][j] << std::endl;
	pass = false;
      }
    }
  }
  
  if(pass){
    std::cout << "Passed.\n";
  }
  else{
    std::cout << "Failed.\n";
  }
  return pass;
}

//Initializations
void initialize_1(Tensor<2,3,double> &A){
  double tmp[3][3] = 
    {
      {1.1,3.2,0.},
      {12.5,9.1,8.3},
      {-0.3,2.4,100.}
    };
  A = convert<3,double>(tmp);
}

void initialize_2(Tensor<2,3,double> &B){
  double tmp[3][3] = 
    {
      {2.6,8.4,0.9},
      {2.3,12.3,9.8},
      {-9.2,-3.6,32.}
    };
  B = convert<3,double>(tmp);
}

void initialize_1(Tensor<2,2,double> &A){
  double tmp[2][2] = 
    {
      {1.1,3.2},
      {12.5,9.1}
    };
  A = convert<2,double>(tmp);
}

void initialize_2(Tensor<2,2,double> &B){
  double tmp[2][2] = 
    {
      {2.6,8.4},
      {2.3,12.3}
    };
  B = convert<2,double>(tmp);
}

void initialize_1(Tensor<1,3,double> &A){
  double tmp[3] = {1.1,3.2,0.};
  A = convert<3,double>(tmp);
}

void initialize_2(Tensor<1,3,double> &B){
  double tmp[3] = {2.6,8.4,0.9};
  B = convert<3,double>(tmp);
}

void initialize_1(Tensor<1,2,double> &A){
  double tmp[2] = {1.1,3.2};
  A = convert<2,double>(tmp);
}

void initialize_2(Tensor<1,2,double> &B){
  double tmp[2] = {2.6,8.4};
  B = convert<2,double>(tmp);
}

void initialize_1(Tensor<4,3,double> &A){
  for(unsigned int i=0; i<3; ++i){
    for(unsigned int j=0; j<3; ++j){
      for(unsigned int k=0; k<3; ++k){
	for(unsigned int l=0; l<3; ++l){
	  A[i][j][k][l] = 11.*(i==j)*(k==l) + 15.*0.5*((i==k)*(j==l)+(i==l)*(j==k));
	}
      }
    }
  }
}

void initialize_1(Tensor<4,2,double> &A){
  for(unsigned int i=0; i<2; ++i){
    for(unsigned int j=0; j<2; ++j){
      for(unsigned int k=0; k<2; ++k){
	for(unsigned int l=0; l<2; ++l){
	  A[i][j][k][l] = 11.*(i==j)*(k==l) + 15.*0.5*((i==k)*(j==l)+(i==l)*(j==k));
	}
      }
    }
  }
}

//Addition tests

/**
 * @ingroup tests
 * Test addition of 2nd order Tensors (dim==3).
 */
bool test_2ndOrder_3D_Add(){
  Tensor<2,3,double> A,B,C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = A+B;

  double tmp[3][3] = 
    {
      {3.7,11.6,0.9},
      {14.8,21.4,18.1},
      {-9.5,-1.2,132.}
    };
  Ctest = convert<3,double>(tmp);

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test addition of 2nd order Tensors (dim==2).
 */
bool test_2ndOrder_2D_Add(){
  Tensor<2,2,double> A,B,C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = A+B;

  double tmp[2][2] = 
    {
      {3.7,11.6},
      {14.8,21.4}
    };
  Ctest = convert<2,double>(tmp);

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test addition of 1st order Tensors (dim==3).
 */
bool test_1stOrder_3D_Add(){
  Tensor<1,3,double> A,B,C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = A+B;

  double tmp[3] = {3.7,11.6,0.9};
  Ctest = convert<3,double>(tmp);

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test addition of 1sd order Tensors (dim==2).
 */
bool test_1stOrder_2D_Add(){
  Tensor<1,2,double> A,B,C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = A+B;

  double tmp[2] = {3.7,11.6};
  Ctest = convert<2,double>(tmp);

  return test(C,Ctest);
}

//Subtraction tests

/**
 * @ingroup tests
 * Test subtraction of 2nd order Tensors (dim==3).
 */
bool test_2ndOrder_3D_Subtract(){
  Tensor<2,3,double> A,B,C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = A-B;

  double tmp[3][3] = 
    {
      {-1.5,-5.2,-0.9},
      {10.2,-3.2,-1.5},
      {8.9,6.,68.}
    };
  Ctest = convert<3,double>(tmp);

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test subtraction of 2nd order Tensors (dim==2).
 */
bool test_2ndOrder_2D_Subtract(){
  Tensor<2,2,double> A,B,C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = A-B;

  double tmp[2][2] = 
    {
      {-1.5,-5.2},
      {10.2,-3.2}
    };
  Ctest = convert<2,double>(tmp);

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test subtraction of 1st order Tensors (dim==3).
 */
bool test_1stOrder_3D_Subtract(){
  Tensor<1,3,double> A,B,C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = A-B;

  double tmp[3] = {-1.5,-5.2,-0.9};
  Ctest = convert<3,double>(tmp);

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test subtraction of 1st order Tensors (dim==2).
 */
bool test_1stOrder_2D_Subtract(){
  Tensor<1,2,double> A,B,C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = A-B;

  double tmp[2] = {-1.5,-5.2};
  Ctest = convert<2,double>(tmp);

  return test(C,Ctest);
}

//Multiplication tests

/**
 * @ingroup tests
 * Test multiplication of 2nd order Tensor with 1st order Tensor (dim==3).
 */
bool test_2ndOrder_3D_multiplication1(){
  Tensor<2,3,double> A;
  Tensor<1,3,double> B,C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = A*B;

  double tmp[3] = {29.74,116.41,109.38};
  Ctest = convert<3,double>(tmp);

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test multiplication of 2nd order Tensor with 1st order Tensor (dim==2).
 */
bool test_2ndOrder_2D_multiplication1(){
  Tensor<2,2,double> A;
  Tensor<1,2,double> B,C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = A*B;

  double tmp[2] = {29.74,108.94};
  Ctest = convert<2,double>(tmp);

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test multiplication of 1st order Tensor with 1st order Tensor (dim==3).
 */
bool test_1stOrder_3D_multiplication1(){
  Tensor<1,3,double> A,B;
  double C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = A*B;

  Ctest = 29.74;

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test multiplication of 2nd order Tensor with 1st order Tensor (dim==2).
 */
bool test_1stOrder_2D_multiplication1(){
  Tensor<1,2,double> A,B;
  double C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = A*B;

  Ctest = 29.74;

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test multiplication of 2nd order Tensors (dim==3).
 */
bool test_2ndOrder_3D_multiplication2(){
  Tensor<2,3,double> A,B,C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = A*B;

  double tmp[3][3] = 
    {
      {10.22,48.6,32.35},
      {-22.93,187.05,366.03},
      {-915.26,-333,3223.25}
    };
  Ctest = convert<3,double>(tmp);

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test multiplication of 2nd order Tensors (dim==2).
 */
bool test_2ndOrder_2D_multiplication2(){
  Tensor<2,2,double> A,B,C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = A*B;

  double tmp[2][2] = 
    {
      {10.22,48.6},
      {53.43,216.93}
    };
  Ctest = convert<2,double>(tmp);

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test double contraction of 2nd order Tensors (dim==3).
 */
bool test_2ndOrder_3D_double_contract(){
  Tensor<2,3,double> A,B;
  double C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = double_contract(A,B);

  Ctest = 3445.88;

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test double contraction of 2nd order Tensors (dim==2).
 */
bool test_2ndOrder_2D_double_contract(){
  Tensor<2,2,double> A,B;
  double C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = double_contract(A,B);

  Ctest = 170.42;

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test double contraction of 4th order Tensor with 2nd order Tensor (dim==3).
 */
bool test_4thOrder_3D_double_contract(){
  Tensor<4,3,double> A;
  Tensor<2,3,double> B,C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = double_contract(A,B);

  double tmp[3][3] = 
    {
      {554.9,80.25,-62.25},
      {80.25,700.4,46.5},
      {-62.25,46.5,995.9}
    };
  Ctest = convert<3,double>(tmp);

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test double contraction of 4th order Tensor with 2nd order Tensor (dim==2).
 */
bool test_4thOrder_2D_double_contract(){
  Tensor<4,2,double> A;
  Tensor<2,2,double> B,C,Ctest;
  initialize_1(A);
  initialize_2(B);
  C = double_contract(A,B);

  double tmp[2][2] = 
    {
      {202.9,80.25},
      {80.25,348.4}
    };
  Ctest = convert<2,double>(tmp);

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test full contraction of 4th order Tensor with 2nd order Tensors (dim==3).
 */
bool test_4thOrder_3D_full_contract(){
  Tensor<4,3,double> D;
  Tensor<2,3,double> A,B;
  double C,Ctest;
  initialize_1(A);
  initialize_1(D);
  initialize_2(B);
  C = full_contract(A,D,B);

  Ctest = 108350.18;

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test full contraction of 4th order Tensor with 2nd order Tensors (dim==2).
 */
bool test_4thOrder_2D_full_contract(){
  Tensor<4,2,double> D;
  Tensor<2,2,double> A,B;
  double C,Ctest;
  initialize_1(A);
  initialize_1(D);
  initialize_2(B);
  C = full_contract(A,D,B);

  Ctest = 4653.555;

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test determinant of 2nd order Tensors (dim==2).
 */
bool test_2ndOrder_2D_determinant(){
  Tensor<2,2,double> A;
  double C,Ctest;
  initialize_1(A);
  C = det(A);
  Ctest = -29.99;

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test determinant of 2nd order Tensors (dim==3).
 */
bool test_2ndOrder_3D_determinant(){
  Tensor<2,3,double> A;
  double C,Ctest;
  initialize_1(A);
  C = det(A);
  Ctest = -3028.88;

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test inverse of 2nd order Tensors (dim==2).
 */
bool test_2ndOrder_2D_inverse(){
  Tensor<2,2,double> A,C,Ctest;
  initialize_1(A);
  C = inv(A);

  double tmp[2][2] = 
    {
      {-0.3034344781593865,0.1067022340780260},
      {0.4168056018672891,-0.0366788929643214}
    };
  Ctest = convert<2,double>(tmp);

  return test(C,Ctest);
}

/**
 * @ingroup tests
 * Test inverse of 2nd order Tensors (dim==3).
 */
bool test_2ndOrder_3D_inverse(){
  Tensor<2,3,double> A,C,Ctest;
  initialize_1(A);
  C = inv(A);

  double tmp[3][3] = 
    {
      {-0.29386439872163966,0.10564961305829218,-0.00876891788383827},
      {0.41351588706056364,-0.03631705448878794,0.00301431552256942},
      {-0.01080597448561845,0.00118855814690579,0.00990134967380682}
    };
  Ctest = convert<3,double>(tmp);

  return test(C,Ctest);
}

