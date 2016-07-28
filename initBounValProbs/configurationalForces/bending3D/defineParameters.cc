//extern "C" {
#include "petiga.h"
//}

#include "IBVPHeaders.h"
//include automatic differentiation library
#include <Sacado.hpp>

int defineParameters(AppCtx& user){

  PetscErrorCode ierr;
 
	user.dim = 3;
	user.GridScale = 1.0;
	user.ADSacado = PETSC_TRUE;
	user.numVars = 162;
	user.uDirichlet = 0.1;
	user.Nx = 100;
	user.Ny = 10;
	user.Nz = 10;
	user.Lx = 10.;
	user.Ly = 1.;
	user.Lz = 1.;
	user.dtVal = 0.001;
	user.skipOutput = 5;
	user.RESTART_IT = 0;
	user.RESTART_TIME = 0.;

	//Anisotropic St. Venant-Kirchhoff
	user.matParam["mu"] = 1e-1;
	user.matParam["betaC"] = 1e-1;
	user.matParam["alphaC"] = 2e-1;
	user.matParam["anisoCoeff"] = 1.;

	//Nonconvex free energy parameters
	user.matParam["Es"] = 0.1;
	user.matParam["Ed"] = 1.0;
	user.matParam["El"] = .25;

  return 0;
}
