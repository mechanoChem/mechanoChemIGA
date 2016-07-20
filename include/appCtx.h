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
	PetscInt NVal;
	PetscReal dtVal;
	PetscInt skipOutput;
	PetscInt RESTART_IT;
	PetscReal RESTART_TIME;
	PetscInt beamRatio; //Ratio of length to width of beam. Keep it an integer to make setting the number of elements easy.

	std::map<std::string,PetscReal> matParam;
/*
	//For measuring anisotropy
	PetscReal F00; //F_{11} - Deformation gradient
	PetscReal P00; //P_{11} - 1st Piola-Kirchhoff stress
	PetscReal Lambda1; //Configurational stretch in e_1 direction

	//Anisotropic St. Venant-Kirchhoff
	PetscReal mu;
	PetscReal alphaC;
	PetscReal betaC;
	PetscReal anisoCoeff;

	//Nonconvex free energy parameters
	PetscReal Es;
	PetscReal Ed;
	PetscReal El;
*/
};

#endif
