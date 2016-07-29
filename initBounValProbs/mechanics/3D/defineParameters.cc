//extern "C" {
#include "petiga.h"
//}

#include "IBVPHeaders.h"
//include automatic differentiation library
#include <Sacado.hpp>

int defineParameters(AppCtx& user){

  PetscErrorCode ierr;
 
	user.dim = 3;
	user.GridScale = 1.;
	user.ADSacado = PETSC_TRUE;
	user.numVars = 81;
	user.uDirichlet = 0.001;
	user.Nx = 50;
	user.Ny = 50;
	user.Nz = 50;
	user.Lx = 1.;
	user.Ly = 1.;
	user.Lz = 1.;
	user.dtVal = 1.0e-2;
	user.skipOutput = 1;
	user.RESTART_IT = 0;
	user.RESTART_TIME = 0.;

	//Elastic free energy parameters
	user.matParam["Es"] = 0.01;
	user.matParam["Ed"] = 1.;
	user.matParam["El"] = 0.1;

  return 0;
}
