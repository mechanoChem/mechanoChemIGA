#ifndef DNN_H
#define DNN_H

#include <Sacado.hpp>
#include <vector>
#include <fstream>
#include <string>
#include <stdexcept>

class DNN {
  
public:
  DNN(){};
  void reinit(unsigned int nLayers);

  template<typename T>
  void eval(const std::vector<T> &x,
	    std::vector<T> &y,
	    std::vector<std::vector<T> > &dy_dx,
	    std::vector<std::vector<std::vector<T> > > &ddy_dxdx);
  
  //private:
  void readVector(std::string filename,std::vector<double> &vec);
  void readMatrix(std::string filename,std::vector<std::vector<double> > &mat);

  template<typename T>
  void actFun(std::vector<T> &a, const std::vector<T> &z);

  template<typename T>
  void actFunDer(std::vector<T> &a, const std::vector<T> &z);

  template<typename T>
  void actFun2ndDer(std::vector<T> &a, const std::vector<T> &z);

  template<typename T>
  void mult(std::vector<T> &z, const std::vector<T> &a, const std::vector<std::vector<double> > &weights);

  unsigned int n_layers, n_features;
  std::vector<std::vector<double> > bias;
  std::vector<std::vector<std::vector<double> > > weights;
    
};

inline void DNN::reinit(unsigned int nLayers){
  n_layers = nLayers;
  bias.resize(n_layers);
  weights.resize(n_layers+1);
  
  //Read in weights/biases for neural net
  std::string bName = "bias_", wName = "weights_", txt = ".txt", filename;
  for (unsigned int i=0; i<n_layers; ++i){
    filename = bName + std::to_string(i) + txt;
    readVector(filename,bias[i]);
    filename = wName + std::to_string(i) + txt;
    readMatrix(filename,weights[i]);
  }
  filename = wName + std::to_string(n_layers) + txt;
  readMatrix(filename,weights[n_layers]);

  n_features = weights[0].size();
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

inline void DNN::readMatrix(std::string filename,std::vector<std::vector<double> > &mat){
  unsigned int m, n;
  char c;
  double d;
  std::ifstream file;
  
  file.open(filename);
  file >> c;
  file >> m; file >> n;
  mat.resize(m,std::vector<double>(n));
  for(unsigned int i=0; i<m; ++i){
    for(unsigned int j=0; j<n; ++j){
      file >> d;
      mat[i][j] = d;
    }
  }
  file.close();
}

//Softplus activation
template<typename T>
void DNN::actFun(std::vector<T> &a, const std::vector<T> &z){
  a.resize(z.size());
  for(unsigned int i=0; i<a.size(); ++i){
    a[i] = std::log(1. + std::exp(z[i]));
  }
}

//Logisitic (sigmoid) function
template<typename T>
void DNN::actFunDer(std::vector<T> &a, const std::vector<T> &z){
  a.resize(z.size());
  for(unsigned int i=0; i<a.size(); ++i){
    a[i] = 1./(1. + std::exp(-z[i]));
  }
}

//Derivative of logisitic function
template<typename T>
void DNN::actFun2ndDer(std::vector<T> &a, const std::vector<T> &z){
  a.resize(z.size());
  for(unsigned int i=0; i<a.size(); ++i){
    a[i] = 1./(2. + std::exp(-z[i]) + std::exp(z[i]));
  }
}

//Vector^T*Matrix
template<typename T>
void DNN::mult(std::vector<T> &z, const std::vector<T> &a, const std::vector<std::vector<double> > &weights){
  double m = weights.size(), n = weights[0].size();
  if(a.size() != m){
    throw std::invalid_argument( "Incompatible matrix and vector sizes." );
  }
  z.resize(n);
  for(unsigned int i=0; i<n; ++i){
    z[i] = 0.;
    for(unsigned int j=0; j<m; ++j){
      z[i] += a[j]*weights[j][i];
    }
  }
}

//Evaluate DNN, its first and second derivatives
template<typename T>
void DNN::eval(const std::vector<T> &x,
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
    actFunDer(beta[j],z);
    beta[j] *= weights[0][j];
    for (unsigned int k=0; k<n_features; ++k){
      actFun2ndDer(gamma[j][k],z);
      gamma[j][k] *= weights[0][j];
      gamma[j][k] *= weights[0][k];
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
void operator+=(std::vector<T>& a,const std::vector<U>& b){
  if(a.size() != b.size()){
    throw std::invalid_argument( "Incompatible vector sizes." );
  }
  for(unsigned int i=0; i<a.size(); ++i){
    a[i] = a[i] + b[i];
  }
}

// This is a component-wise multiplication of the two vectors.
template<typename T, typename U>
void operator*=(std::vector<T>& a,const std::vector<U>& b){
  if(a.size() != b.size()){
    throw std::invalid_argument( "Incompatible vector sizes." );
  }
  for(unsigned int i=0; i<a.size(); ++i){
    a[i] = a[i]*b[i];
  }
}

#endif //DNN_H
