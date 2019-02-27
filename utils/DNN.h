#ifndef DNN_H
#define DNN_H

//#include <Sacado.hpp>
#include <Sacado_Fad_BLAS.hpp>
#include <vector>
#include <fstream>
#include <string>
#include <stdexcept>

extern "C"{
  void dgemv_(char* TRANS, const int* M, const int* N,
	      double* alpha, double* A, const int* LDA, double* X,
	      const int* INCX, double* beta, double* C, const int* INCY);
}

template<typename T>
inline void mult_spec(std::vector<T> &z, std::vector<T> &a, std::vector<double> &weights);

template<>
inline void mult_spec(std::vector<double> &z, std::vector<double> &a, std::vector<double> &weights){

  int mn = weights.size();
  int m = a.size();
  int n = mn/m;
  if(mn%m != 0){
    throw std::invalid_argument( "Incompatible matrix and vector sizes." );
  }
  z.resize(n);

  char tran = 'N';
  double coeff1 = 1., coeff0 = 0.;
  int step = 1;
  dgemv_(&tran,&n,&m,&coeff1,&weights[0],&n,&a[0],&step,&coeff0,&z[0],&step);
}

template<>
inline void mult_spec(std::vector<Sacado::Fad::DFad<double> > &z, std::vector<Sacado::Fad::DFad<double> > &a, std::vector<double> &weights){

  int mn = weights.size();
  int m = a.size();
  int n = mn/m;
  if(mn%m != 0){
    throw std::invalid_argument( "Incompatible matrix and vector sizes." );
  }
  z.resize(n);

  Sacado::Fad::BLAS<int,Sacado::Fad::DFad<double> > blas;
  blas.GEMV(Teuchos::NO_TRANS,n,m,1.,&weights[0],n,&a[0],1,0.,&z[0],1);

}


class DNN {
  
 public:
  DNN(){};
  void reinit(unsigned int nLayers);
  
  template<typename T>
    void eval(std::vector<T> &x,
	      std::vector<T> &y,
	      std::vector<std::vector<T> > &dy_dx,
	      std::vector<std::vector<std::vector<T> > > &ddy_dxdx);
  
 private:
  void readVector(std::string filename,std::vector<double> &vec);
  //void readMatrix(std::string filename,std::vector<std::vector<double> > &mat);
  void readMatrix(std::string filename,std::vector<double> &mat,unsigned int &m);
  
  template<typename T>
    void actFun(std::vector<T> &a, std::vector<T> &z);

  template<typename T>
    void actFunDer(std::vector<T> &a, std::vector<T> &z);

  template<typename T>
    void actFun2ndDer(std::vector<T> &a, std::vector<T> &z);

  template<typename T>
    void mult(std::vector<T> &z, std::vector<T> &a, std::vector<double> &weights);
    //void mult(std::vector<T> &z, const std::vector<T> &a, const std::vector<std::vector<double> > &weights);

  unsigned int n_layers, n_features;
  std::vector<std::vector<double> > bias;
  //std::vector<std::vector<std::vector<double> > > weights;
  std::vector<std::vector<double> > weights;
  std::vector<unsigned int> weights_m;
    
};

inline void DNN::reinit(unsigned int nLayers){
  n_layers = nLayers;
  bias.resize(n_layers);
  weights.resize(n_layers+1);
  weights_m.resize(n_layers+1);
  
  //Read in weights/biases for neural net
  std::string bName = "bias_", wName = "weights_", txt = ".txt", filename;
  for (unsigned int i=0; i<n_layers; ++i){
    filename = bName + std::to_string(i) + txt;
    readVector(filename,bias[i]);
    filename = wName + std::to_string(i) + txt;
    readMatrix(filename,weights[i],weights_m[i]);
  }
  filename = wName + std::to_string(n_layers) + txt;
  readMatrix(filename,weights[n_layers],weights_m[n_layers]);

  n_features = weights_m[0];
}

inline void DNN::readVector(std::string filename,std::vector<double> &vec){
  unsigned int m;
  char c;
  double d;
  std::ifstream file;
  
  file.open(filename);
  file >> c;
  file >> m;
  vec.resize(m);
  for(unsigned int i=0; i<m; ++i){
    file >> d;
    vec[i] = d;
  }
  file.close();
}

//inline void DNN::readMatrix(std::string filename,std::vector<std::vector<double> > &mat){
inline void DNN::readMatrix(std::string filename,std::vector<double> &mat,unsigned int &m){
  unsigned int n;
  char c;
  double d;
  std::ifstream file;
  
  file.open(filename);
  file >> c;
  file >> m; file >> n;
  /*
  mat.resize(m,std::vector<double>(n));
  for(unsigned int i=0; i<m; ++i){
    for(unsigned int j=0; j<n; ++j){
      file >> d;
      mat[i][j] = d;
    }
  }
  */
  mat.resize(m*n);
  for(unsigned int i=0; i<m*n; ++i){
    file >> d;
    mat[i] = d;
  }
  file.close();
}

//Softplus activation
template<typename T>
void DNN::actFun(std::vector<T> &a, std::vector<T> &z){
  a.resize(z.size());
  for(unsigned int i=0; i<a.size(); ++i){
    a[i] = std::log(1. + std::exp(z[i]));
  }
}

//Logisitic (sigmoid) function
template<typename T>
void DNN::actFunDer(std::vector<T> &a, std::vector<T> &z){
  a.resize(z.size());
  for(unsigned int i=0; i<a.size(); ++i){
    a[i] = 1./(1. + std::exp(-z[i]));
  }
}

//Derivative of logisitic function
template<typename T>
void DNN::actFun2ndDer(std::vector<T> &a, std::vector<T> &z){
  a.resize(z.size());
  for(unsigned int i=0; i<a.size(); ++i){
    a[i] = 1./(2. + std::exp(-z[i]) + std::exp(z[i]));
  }
}

//Vector^T*Matrix
template<typename T>
void DNN::mult(std::vector<T> &z, std::vector<T> &a, std::vector<double> &weights){

  mult_spec(z,a,weights);

  /*
void DNN::mult(std::vector<T> &z, const std::vector<T> &a, const std::vector<std::vector<double> > &weights){
  int m = weights.size(), n = weights[0].size();
  if(a.size() != m){
    throw std::invalid_argument( "Incompatible matrix and vector sizes." );
  }
  */
  /*
  int mn = weights.size();
  int m = a.size();
  int n = mn/m;
  if(mn%m != 0){
    throw std::invalid_argument( "Incompatible matrix and vector sizes." );
  }
  z.resize(n);
  */

  /* My implementation
  for(unsigned int i=0; i<n; ++i){
    z[i] = 0.;
    for(unsigned int j=0; j<m; ++j){
      //z[i] += a[j]*weights[j][i];
      z[i] += a[j]*weights[n*j+i];
    }
  }
  // */

  /*
  //Blas
  if(std::is_same<T,double>::value || std::is_same<T,const double>::value){
    char tran = 'N';
    double coeff1 = 1., coeff0 = 0.;
    int step = 1;
    dgemv_(&tran,&n,&m,&coeff1,&weights[0],&n,&a[0],&step,&coeff0,&z[0],&step);
  }
  else{
    Sacado::Fad::BLAS<int,Sacado::Fad::DFad<double> > blas;
    blas.GEMV(Teuchos::NO_TRANS,n,m,1.,&weights[0],n,&a[0],1,0.,&z[0],1);
  }
  // */
}

//Evaluate DNN, its first and second derivatives
template<typename T>
void DNN::eval(std::vector<T> &x,
	       std::vector<T> &y,
	       std::vector<std::vector<T> > &dy_dx,
	       std::vector<std::vector<std::vector<T> > > &ddy_dxdx){

  if(x.size() != n_features){
    throw std::invalid_argument( "Input vector has incompatible number of features." );
  }
  dy_dx.resize(n_features);
  ddy_dxdx.resize(n_features,std::vector<std::vector<T> >(n_features));
  std::vector<T> z, alpha, g1z, g2z, zgamma;
  std::vector<std::vector<T> > beta(n_features), zbeta(n_features);
  std::vector<std::vector<std::vector<T> > > gamma(n_features,std::vector<std::vector<T> >(n_features));

  mult(z,x,weights[0]);
  z += bias[0];
  actFun(alpha,z);
  //For each partial derivative
  for (unsigned int j=0; j<n_features; ++j){
    unsigned int n = weights[0].size()/n_features;
    std::vector<double> weights_0j(&weights[0][j*n],&weights[0][(j+1)*n]);
    actFunDer(beta[j],z);
    //beta[j] *= weights[0][j];
    beta[j] *= weights_0j;
    for (unsigned int k=0; k<n_features; ++k){
      std::vector<double> weights_0k(&weights[0][k*n],&weights[0][(k+1)*n]);
      actFun2ndDer(gamma[j][k],z);
      //gamma[j][k] *= weights[0][j];
      //gamma[j][k] *= weights[0][k];
      gamma[j][k] *= weights_0j;
      gamma[j][k] *= weights_0k;
    }
  }
  for (unsigned int i=1; i<n_layers; ++i){
    mult(z,alpha,weights[i]);
    z += bias[i]; //  z^T = alpha^T*weights + bias^T
    actFun(alpha,z); // alpha = g(z)
    actFunDer(g1z,z); //g'(z)
    actFun2ndDer(g2z,z); //g''(z)

    for (unsigned int j=0; j<n_features; ++j){
      mult(zbeta[j],beta[j],weights[i]); //z_beta_j^T = beta_j^T*weights
      beta[j] = g1z; //actFunDer(beta[j],z);
      beta[j] *= zbeta[j]; //beta_j = g'(z)*z_beta_j
    }

    for (unsigned int j=0; j<n_features; ++j){
      for (unsigned int k=0; k<=j; ++k){
	mult(zgamma,gamma[j][k],weights[i]); //z_gamma_jk^T = gamma_jk^T*weights
	zgamma *= g1z; //g'(z)*z_gamma_jk
	gamma[j][k] = g2z;
	gamma[j][k] *= zbeta[j];
	gamma[j][k] *= zbeta[k];
	gamma[j][k] += zgamma;
	//gamma_jk = g'(z)*z_gamma_jk + g''(z)*z_beta_j*z_beta_k
	gamma[k][j] = gamma[j][k]; //Take advantage of symmetry
      }
    }
  }
  mult(y,alpha,weights[n_layers]); // y = alpha^T*weights
  for (unsigned int j=0; j<n_features; ++j){
    mult(dy_dx[j],beta[j],weights[n_layers]); // dy/dx_j = beta_j^T*weights
    for (unsigned int k=0; k<n_features; ++k){
      mult(ddy_dxdx[j][k],gamma[j][k],weights[n_layers]); // ddy/dx_jdx_k = gamma_jk^T*weights
    }
  }
}

/*
 * Related, nonmember functions
 */

template<typename T, typename U>
  void operator+=(std::vector<T>& a,std::vector<U>& b){
  if(a.size() != b.size()){
    throw std::invalid_argument( "Incompatible vector sizes." );
  }
  for(unsigned int i=0; i<a.size(); ++i){
    a[i] = a[i] + b[i];
  }
}

// This is a component-wise multiplication of the two vectors.
template<typename T, typename U>
  void operator*=(std::vector<T>& a,std::vector<U>& b){
  if(a.size() != b.size()){
    throw std::invalid_argument( "Incompatible vector sizes." );
  }
  for(unsigned int i=0; i<a.size(); ++i){
    a[i] = a[i]*b[i];
  }
}

#endif //DNN_H
