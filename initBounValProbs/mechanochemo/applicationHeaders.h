#ifndef applicationHeaders_
#define applicationHeaders_

//extern "C" {
#include "petiga.h"
//}

typedef struct {
  IGA iga;
  TS* ts;
  Vec *U, *U0;
  PetscReal dt;
  PetscReal he;
  PetscReal norm;
  PetscReal lambda; //load parameter
	Vec *totalEnergy;
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

	//Anisotropic St. Venant-Kirchhoff
	PetscReal mu;
	PetscReal alphaC;
	PetscReal betaC;

	//Nonconvex free energy parameters
	PetscReal Es;
	PetscReal Ed;
	PetscReal El;

} AppCtx;

template<unsigned int DIM>
int setup(AppCtx& user,Vec *U,Vec *U0,Vec *totalEnergy,TS &ts);

template<unsigned int dim>
int boundaryConditions(AppCtx& user, double scale);

template<int dim>
int timeStepSetup(AppCtx& user, TS& ts);

template <int dim>
PetscErrorCode loadStep(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);

template <int dim>
PetscErrorCode adaptiveTimeStep(TS ts,PetscInt it_number,PetscReal c_time,Vec U,void *mctx);

PetscErrorCode SNESConvergedTest(SNES snes, PetscInt it,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx);

PetscErrorCode SNESConvergedTest_Interactive(SNES snes, PetscInt it,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx);

int setConvergenceTest(AppCtx& user, TS& ts);

int defineParameters(AppCtx& user);

template<class T>
T alphaI(PetscReal alphaC,T LambdaI);

template<class T>
T d_alphaL(PetscReal alphaC,T LambdaL);

/*template<class T, unsigned int DIM>
T PiJ(double mu, T (&alpha)[DIM], T (&beta)[DIM][DIM], T (&F)[DIM][DIM], T (&E)[DIM][DIM], int i, int J);*/

#endif
