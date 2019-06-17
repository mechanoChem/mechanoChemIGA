#ifndef bsplines_H
#define bsplines_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
class bsplines{

 public:
  bsplines(){};
  void reinit(unsigned int p_,
	      std::vector<std::vector<double> > X_,
	      std::vector<std::vector<double> > CPts_);

  void reinit(std::string fname);
  /*
  template<typename T,typename U>
  auto divZeroSafe(T a,U b);

  template<typename T>
  T N_ip (unsigned int i, unsigned int j, unsigned int P, T x);

  template<typename T>
  T dkN_dxk (unsigned int k, unsigned int i, unsigned int j, unsigned int P, T x);
  */
  template<typename T>
  void eval(const std::vector<T> &x,
	    T &y,
	    std::vector<T> &dy_dx,
	    std::vector<std::vector<T> > &ddy_dxdx);

 private:
  unsigned int p;
  std::vector<std::vector<double> > X;
  std::vector<std::vector<double> > CPts;
};

inline void bsplines::reinit(unsigned int p_,
			     std::vector<std::vector<double> > X_,
			     std::vector<std::vector<double> > CPts_){

  p = p_;
  X = X_;
  CPts = CPts_;
  
}

inline void bsplines::reinit(std::string fname){

  //Read in knot vectors and control points
  unsigned int m, n;
  char c; //To catch the # sign in the header
  double d;
  std::ifstream file;
  std::string filename = fname+"Cpts.txt";
  file.open(filename);
  file >> c;
  file >> p; file >> m; file >> n;
  CPts.resize(m,std::vector<double>(n));
  for(unsigned int i=0; i<m; ++i){
    for(unsigned int j=0; j<n; ++j){
      file >> d;
      CPts[i][j] = d;
    }
  }
  file.close();

  filename = fname+"knots.txt";
  file.open(filename);
  file >> c;
  file >> p; file >> m; file >> n;
  X.resize(2);
  X[0].resize(m);
  X[1].resize(n);
  for(unsigned int i=0; i<m; ++i){
      file >> d;
      X[0][i] = d;
  }
  for(unsigned int i=0; i<n; ++i){
    file >> d;
    X[1][i] = d;
  }
  file.close();

}
/*
template<typename T,typename U>
auto bsplines::divZeroSafe(T a,U b){
  if(b==0){ return 0.*a; }
  else{ return a/b; }
}

template<typename T>
T bsplines::N_ip (unsigned int i, unsigned int j, unsigned int P, T x){
  T N;
  if(P == 0){
    N = (X[j][i] <= x && x < X[j][i+1]);
  }
  else{
    N = (x - X[j][i])*divZeroSafe(N_ip(i,j,P-1,x),X[j][i+P] - X[j][i]) +
      (X[j][i+P+1] - x)*divZeroSafe(N_ip(i+1,j,P-1,x),(X[j][i+P+1] - X[j][i+1]));
  }
  return N;
}

template<typename T>
T bsplines::dkN_dxk (unsigned int k, unsigned int i, unsigned int j, unsigned int P, T x){
  // k^th derivative at knot i, poly. order P, evaluated at variable x_j = x
  T dkNdxk;
  if (k == 0){
    dkNdxk = N_ip(i,j,P,x);
  }
  else if(P == 0){
    dkNdxk = 0.;
  }
  else{
    dkNdxk = P*divZeroSafe(dkN_dxk(k-1,i,j,P-1,x),X[j][i+P] - X[j][i]) -
      P*divZeroSafe(dkN_dxk(k-1,i+1,j,P-1,x),(X[j][i+P+1] - X[j][i+1]));
  }
  return dkNdxk;
}

template<typename T>
void bsplines::eval(const std::vector<T> &x,
		    T &y,
		    std::vector<T> &dy_dx,
		    std::vector<std::vector<T> > &ddy_dxdx){

  y = 0.;
  dy_dx.resize(2,0.);
  ddy_dxdx.resize(2,std::vector<T>(2,0.));
  unsigned int n1 = X[0].size()-p-1;
  unsigned int n2 = X[1].size()-p-1;
  for(unsigned int j=0; j<n1; ++j){
    for(unsigned int k=0; k<n2; ++k){
      y += CPts[j][k]*N_ip(j,0,p,x[0])*N_ip(k,1,p,x[1]);
      dy_dx[0] += CPts[j][k]*dkN_dxk(1,j,0,p,x[0])*dkN_dxk(0,k,1,p,x[1]);
      dy_dx[1] += CPts[j][k]*dkN_dxk(0,j,0,p,x[0])*dkN_dxk(1,k,1,p,x[1]);
      ddy_dxdx[0][0] += CPts[j][k]*dkN_dxk(2,j,0,p,x[0])*dkN_dxk(0,k,1,p,x[1]);
      ddy_dxdx[1][1] += CPts[j][k]*dkN_dxk(0,j,0,p,x[0])*dkN_dxk(2,k,1,p,x[1]);
      ddy_dxdx[0][1] += CPts[j][k]*dkN_dxk(1,j,0,p,x[0])*dkN_dxk(1,k,1,p,x[1]);
    }
  }
  ddy_dxdx[1][0] = ddy_dxdx[0][1];

}
*/

// Evaluate splines using de Boors' algorithm from de Boor, "On Calculating with B-Splines", 1972.
template<typename T>
void bsplines::eval(const std::vector<T> &x,
		    T &y,
		    std::vector<T> &dy_dx,
		    std::vector<std::vector<T> > &ddy_dxdx){

  // First find out in which knot span the point x is located
  unsigned int i, j, n = X[0].size();
  for (unsigned int k=0; k<n-1; ++k){
    //std::cout << k << " " << X[0][k] << std::endl;
    if ((x[0] >= X[0][k]) && (x[0] < X[0][k+1])){
      i = k;
      break;
    }
  }
  n = X[1].size();
  for (unsigned int k=0; k<n-1; ++k){
    //std::cout << k << " " << X[1][k] << std::endl;
    if ((x[1] >= X[1][k]) && (x[1] < X[1][k+1])){
      j = k;
      break;
    }
  }
  //std::cout << "i: " << i << ", j: " << j << std::endl;

  //Now, compute all necessary basis functions
  std::vector<std::vector <T> > N1(p+1,std::vector<T>(p+1,0.)), N2(p+1,std::vector<T>(p+1,0.));
  std::vector<T> DP1(p,0.), DM1(p,0.), DP2(p,0.), DM2(p,0.);
  N1[0][0] = 1.;
  N2[0][0] = 1.;
  T M1, M2;
  for (unsigned int s=0; s<p; ++s){
    DP1[s] = X[0][i+s+1]-x[0];
    DM1[s] = x[0] - X[0][i-s];
    DP2[s] = X[1][j+s+1]-x[1];
    DM2[s] = x[1] - X[1][j-s];
    N1[0][s+1] = 0;
    N2[0][s+1] = 0; //11*p-1 to here
    for (unsigned int r=0; r<=s; ++r){
      M1 = N1[r][s]/(DP1[r]+DM1[s-r]); //3
      N1[r][s+1] += DP1[r]*M1; //3
      N1[r+1][s+1] = DM1[s-r]*M1; //4
      M2 = N2[r][s]/(DP2[r]+DM2[s-r]); //3
      N2[r][s+1] += DP2[r]*M2; //3
      N2[r+1][s+1] = DM2[s-r]*M2; //4, 21*p*(p+1)/2
    }
  } //Total flop count to this point: //21*p*(p+1)/2 + 11*p - 1

  //Evaluate the function and derivatives
  y = 0.;
  dy_dx.resize(2,0.);
  ddy_dxdx.resize(2,std::vector<T>(2,0.));
  for (unsigned int r1=0; r1<p+1; ++r1){
    for (unsigned int r2=0; r2<p+1; ++r2){
      y += CPts[r1+i-p][r2+j-p]*N1[r1][p]*N2[r2][p]; //7*(p+1)*(p+1) + p*p
    }
  }
  for (unsigned int r1=1; r1<p+1; ++r1){
    for (unsigned int r2=0; r2<p+1; ++r2){
      dy_dx[0] += p*(CPts[r1+i-p][r2+j-p]-CPts[r1+i-p-1][r2+j-p])/(X[0][r1+i]-X[0][r1+i-p])*N1[r1-1][p-1]*N2[r2][p]; 
      dy_dx[1] += p*(CPts[r2+i-p][r1+j-p]-CPts[r2+i-p][r1+j-p-1])/(X[1][r1+j]-X[1][r1+j-p])*N1[r2][p]*N2[r1-1][p-1]; //40*p*(p+1) + p*(p-1)
    }
  }
  for (unsigned int r1=1; r1<p+1; ++r1){
    for (unsigned int r2=1; r2<p+1; ++r2){
      ddy_dxdx[0][1] += p*p*(CPts[r1+i-p][r2+j-p]-CPts[r1+i-p-1][r2+j-p]-CPts[r1+i-p][r2+j-p-1]+CPts[r1+i-p-1][r2+j-p-1])/(X[0][r1+i]-X[0][r1+i-p])/(X[1][r2+j]-X[1][r2+j-p])*N1[r1-1][p-1]*N2[r2-1][p-1]; //42*p*p + (p-1)*(p-1)
    }
  }
  for (unsigned int r1=2; r1<p+1; ++r1){
    for (unsigned int r2=0; r2<p+1; ++r2){
      ddy_dxdx[0][0] += p*(p-1)*(CPts[r1+i-p][r2+j-p]-2*CPts[r1+i-p-1][r2+j-p]+CPts[r1+i-p-2][r2+j-p])/(X[0][r1+i]-X[0][r1+i-p])/(X[0][r1+i-1]-X[0][r1+i-p])*N1[r1-2][p-2]*N2[r2][p]; 
      ddy_dxdx[1][1] += p*(p-1)*(CPts[r2+i-p][r1+j-p]-2*CPts[r2+i-p][r1+j-p-1]+CPts[r2+i-p][r1+j-p-2])/(X[1][r1+j]-X[1][r1+j-p])/(X[1][r1+j-1]-X[1][r1+j-p])*N1[r2][p]*N2[r1-2][p-2]; //72*(p+1)*(p-1) + p*(p-2)
    }
  }
  ddy_dxdx[1][0] = ddy_dxdx[0][1];
} // total: 165*p*p + 70*p ...

#endif //bsplines_H
