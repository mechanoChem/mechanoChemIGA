#include "test_tensor.h"
#include <iostream>


int main(int argc, char *argv[]) {
 
  unsigned int total_pass = 0;
  bool pass;

  //Test addition (2D & 3D, 1st & 2nd order tensors)
  total_pass += test_2ndOrder_3D_Add();
  total_pass += test_2ndOrder_2D_Add();
  total_pass += test_1stOrder_3D_Add();
  total_pass += test_1stOrder_2D_Add();

  //Test subtraction (2D & 3D, 1st & 2nd order tensors)
  total_pass += test_2ndOrder_3D_Subtract();
  total_pass += test_2ndOrder_2D_Subtract();
  total_pass += test_1stOrder_3D_Subtract();
  total_pass += test_1stOrder_2D_Subtract();

  //Test multiplication
  total_pass += test_2ndOrder_3D_multiplication1();
  total_pass += test_2ndOrder_2D_multiplication1();
  total_pass += test_1stOrder_3D_multiplication1();
  total_pass += test_1stOrder_2D_multiplication1();
  total_pass += test_2ndOrder_3D_multiplication2();
  total_pass += test_2ndOrder_2D_multiplication2();

  //Test contractions
  total_pass += test_2ndOrder_3D_double_contract();
  total_pass += test_2ndOrder_2D_double_contract();
  total_pass += test_4thOrder_3D_double_contract();
  total_pass += test_4thOrder_2D_double_contract();
  total_pass += test_4thOrder_3D_full_contract();
  total_pass += test_4thOrder_2D_full_contract();

  //Test determinant, inverse
  total_pass += test_2ndOrder_2D_determinant();
  total_pass += test_2ndOrder_3D_determinant();
  total_pass += test_2ndOrder_2D_inverse();
  total_pass += test_2ndOrder_3D_inverse();
  
  std::cout << "Total tests passed: " << total_pass << std::endl;

  return 0;
}
