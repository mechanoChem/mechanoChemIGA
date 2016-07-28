#ifndef appCtx_h
#define appCtx_h

//extern "C" {
#include "petiga.h"
//}
#include <map>
#include <string>

struct AppCtx{
  IGA iga;
  TS* ts;
  Vec *U, *U0;
  PetscReal dt;
  PetscReal he;
  PetscReal norm;
  PetscReal lambda; //load parameter
	PetscInt dim;
	PetscReal GridScale;
	PetscBool ADSacado;
	PetscInt numVars;
	PetscReal uDirichlet;
	PetscInt Nx, Ny, Nz; //Number of elements in each direction (note that Nz is not used if dim = 2)
	PetscReal Lx, Ly, Lz; //Dimensions of body
	PetscReal dtVal;
	PetscInt skipOutput;
	PetscInt RESTART_IT;
	PetscReal RESTART_TIME;
	PetscInt beamRatio; //Ratio of length to width of beam. Keep it an integer to make setting the number of elements easy.

	std::map<std::string,PetscReal> matParam;
};

#endif
