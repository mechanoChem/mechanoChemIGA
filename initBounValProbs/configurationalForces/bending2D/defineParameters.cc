//extern "C" {
#include "petiga.h"
//}

#include "applicationHeaders.h"
//include automatic differentiation library
#include <Sacado.hpp>

int defineParameters(AppCtx& user){

  PetscErrorCode ierr;
 
	user.dim = 2;
	user.GridScale = 1.0;
	user.ADSacado = PETSC_TRUE;
	user.numVars = 36;
	user.uDirichlet = 0.1;
	user.NVal = 80;
	user.dtVal = 0.001;
	user.skipOutput = 5;
	user.RESTART_IT = 0;
	user.RESTART_TIME = 0.;
	user.beamRatio = 10; //Ratio of length to width of beam. Keep it an integer to make setting the number of elements easy.

	//Anisotropic St. Venant-Kirchhoff
	user.matParam["mu"] = 1e-1;
	user.matParam["betaC"] = 1e-1;
	user.matParam["alphaC"] = 2e-1;
	user.matParam["anisoCoeff"] = 1.;

	//Nonconvex free energy parameters
	user.matParam["Es"] = 0.1;
	user.matParam["Ed"] = 1.0;
	user.matParam["El"] = .1;

  return 0;
}
