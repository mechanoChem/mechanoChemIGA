#include "../applicationHeaders.h"
//extern "C" {
#include "petiga.h"
//}
#include "../../../include/genericHeaders.h"
#include "../../../src/physicsHeaders.h"

template<unsigned int DIM>
int setup(AppCtx& user,Vec *U,Vec *U0,TS &ts){

  PetscErrorCode ierr;
 
  //application context objects and parameters
  user.dt=user.dtVal;
  user.he=user.GridScale*1.0/user.NVal;
  user.lambda=1.;
  PetscInt p=2;
  const unsigned int DOF=2*DIM;

  user.ts=&ts;
  user.U0=U0;
  user.U=U;
  IGA iga;
  user.iga = iga;

  PetscPrintf(PETSC_COMM_WORLD,"initializing...\n");
  init<DIM,DOF>(user, user.NVal, p);

  //initial conditions
  ierr = IGACreateVec(user.iga,user.U);CHKERRQ(ierr);
  ierr = IGACreateVec(user.iga,user.U0);CHKERRQ(ierr);

	if(user.RESTART_IT==0){
	  ierr = FormInitialCondition(user.iga, *user.U0, &user);CHKERRQ(ierr);
	}
	else{
		char filename[256];
		sprintf(filename,"outU%d.dat",user.RESTART_IT);  
		ierr = IGAReadVec(user.iga,*user.U0,filename);CHKERRQ(ierr); //Read in vector to restart at step RESTART_IT
	}
  ierr = VecCopy(*user.U0, *user.U);CHKERRQ(ierr);

  //Set mat type
  ierr = IGASetMatType(user.iga,MATAIJ);CHKERRQ(ierr); //For superlu_dist (still works for gmres, etc.)

  return 0;
}

template int setup<2>(AppCtx& user,Vec *U,Vec *U0,TS &ts);
template int setup<3>(AppCtx& user,Vec *U,Vec *U0,TS &ts);
