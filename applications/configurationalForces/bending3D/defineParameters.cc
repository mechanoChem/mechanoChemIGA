#include "../applicationHeaders.h"
//extern "C" {
#include "petiga.h"
//}

//include automatic differentiation library
#include <Sacado.hpp>

int defineParameters(AppCtx& user){

  PetscErrorCode ierr;
 
	user.dim = 3;
	user.GridScale = 1.0;
	user.ADSacado = PETSC_TRUE;
	user.numVars = 162;
	user.uDirichlet = 0.1;
	user.NVal = 10;
	user.dtVal = 0.0001;
	user.skipOutput = 1;
	user.RESTART_IT = 0;
	user.RESTART_TIME = 0.;

	user.beamRatio = 10; //Ratio of length to width of beam. Keep it an integer to make setting the number of elements easy.

	//Anisotropic St. Venant-Kirchhoff
	user.mu = 1e-1;
	user.betaC = 1e-1;
	user.alphaC = 2e-1;

	//Nonconvex free energy parameters
	user.Es = 0.1;
	user.Ed = 1.0;
	user.El = .25;

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
