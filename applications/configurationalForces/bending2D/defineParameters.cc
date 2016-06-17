#include "../applicationHeaders.h"
//extern "C" {
#include "petiga.h"
//}

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

  return 0;
}
