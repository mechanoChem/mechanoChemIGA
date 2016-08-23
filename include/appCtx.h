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
  PetscReal he;
  PetscReal norm;
  PetscReal lambda; //load parameter
	PetscInt dim;
	PetscReal GridScale;
	PetscBool ADSacado;
	PetscReal uDirichlet;
	PetscInt Nx, Ny, Nz; //Number of elements in each direction (note that Nz is not used if dim = 2)
	PetscReal Lx, Ly, Lz; //Dimensions of body
	PetscReal dtVal;
	PetscInt skipOutput;
	PetscInt RESTART_IT;
	PetscReal RESTART_TIME;

	std::map<std::string,PetscReal> matParam;

	//Define default values
	AppCtx() : dim(3), GridScale(1.), ADSacado(PETSC_TRUE), uDirichlet(0.), Nx(2), Ny(2), Nz(2), Lx(1.), Ly(1.), Lz(1.), dtVal(1.e-5), skipOutput(1), RESTART_IT(0), RESTART_TIME(0.) {}
};

#endif
