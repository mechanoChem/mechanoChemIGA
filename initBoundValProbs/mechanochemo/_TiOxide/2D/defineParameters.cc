#include "../applicationHeaders.h"
//extern "C" {
#include "petiga.h"
//}

//include automatic differentiation library
#include <Sacado.hpp>

int defineParameters(AppCtx& user){

  PetscErrorCode ierr;
 
	user.dim = 2;
	user.GridScale = 1.0;
	user.ADSacado = PETSC_TRUE;
	user.numVars = 36;
	user.uDirichlet = 0.1;
	user.NVal = 5;
	user.dtVal = 0.00001;
	user.skipOutput = 1;
	user.RESTART_IT = 0;
	user.RESTART_TIME = 0.;

	//Anisotropic St. Venant-Kirchhoff
	user.mu = 1e5;
	user.betaC = 1e5;
	user.alphaC = 2e5;

	//Nonconvex free energy parameters
	user.Es = 0.1;
	user.Ed = 1.0;
	user.El = .025;

  return 0;
}

template<class T>
T alphaI(PetscReal alphaC,T LambdaI){
	return alphaC*LambdaI;
}

template Sacado::Fad::SFad<double,27> alphaI<Sacado::Fad::SFad<double,27> >(PetscReal alphaC,Sacado::Fad::SFad<double,27> LambdaI);
template Sacado::Fad::SFad<double,36> alphaI<Sacado::Fad::SFad<double,36> >(PetscReal alphaC,Sacado::Fad::SFad<double,36> LambdaI);
template Sacado::Fad::SFad<double,108> alphaI<Sacado::Fad::SFad<double,108> >(PetscReal alphaC,Sacado::Fad::SFad<double,108> LambdaI);
template Sacado::Fad::SFad<double,162> alphaI<Sacado::Fad::SFad<double,162> >(PetscReal alphaC,Sacado::Fad::SFad<double,162> LambdaI);
template PetscReal alphaI<PetscReal>(PetscReal alphaC,PetscReal LambdaI);

template<class T>
T d_alphaL(PetscReal alphaC,T LambdaL){
	return alphaC/LambdaL;
}

template Sacado::Fad::SFad<double,27> d_alphaL<Sacado::Fad::SFad<double,27> >(PetscReal alphaC,Sacado::Fad::SFad<double,27> LambdaL);
template Sacado::Fad::SFad<double,36> d_alphaL<Sacado::Fad::SFad<double,36> >(PetscReal alphaC,Sacado::Fad::SFad<double,36> LambdaL);
template Sacado::Fad::SFad<double,108> d_alphaL<Sacado::Fad::SFad<double,108> >(PetscReal alphaC,Sacado::Fad::SFad<double,108> LambdaL);
template Sacado::Fad::SFad<double,162> d_alphaL<Sacado::Fad::SFad<double,162> >(PetscReal alphaC,Sacado::Fad::SFad<double,162> LambdaL);
template PetscReal d_alphaL<PetscReal>(PetscReal alphaC,PetscReal LambdaL);

/*template<class T, unsigned int DIM>
T PiJ(double mu, T (&alpha)[DIM], T (&beta)[DIM][DIM], T (&F)[DIM][DIM], T (&E)[DIM][DIM], int i, int J){
	return ((alpha[J]-2*mu-beta[J][J])*F[i][J]*E[J][J] + 
					(beta[J][0]*E[0][0]+beta[J][1]*E[1][1])*F[i][J] + 
					2*mu*(F[i][0]*E[0][J]+F[i][1]*E[1][J]));
}

template Sacado::Fad::SFad<double,36> PiJ<Sacado::Fad::SFad<double,36>,2>(double mu,
				 Sacado::Fad::SFad<double,36> (&alpha)[2],
				 Sacado::Fad::SFad<double,36> (&beta)[2][2],
				 Sacado::Fad::SFad<double,36> (&F)[2][2],
				 Sacado::Fad::SFad<double,36> (&E)[2][2],
				 int i, int J);
template PetscReal PiJ<PetscReal,2>(double mu,
				 PetscReal (&alpha)[2],
				 PetscReal (&beta)[2][2],
				 PetscReal (&F)[2][2],
				 PetscReal (&E)[2][2],
				 int i, int J);*/
