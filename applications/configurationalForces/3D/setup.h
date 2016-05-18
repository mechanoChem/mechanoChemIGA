#ifndef setup_
#define setup_

template<unsigned int dim>
inline int setup(AppCtx& user,Vec *U,Vec *U0){

  PetscErrorCode ierr;
 
  //application context objects and parameters
  user.dt=dtVal;
  user.he=GridScale*1.0/NVal;
  user.lambda=1.;
  PetscInt p=2;
  const unsigned int DOF=2*dim;

  user.U0=U0;
  user.U=U;
  IGA iga;
  user.iga = iga;

  PetscPrintf(PETSC_COMM_WORLD,"initializing...\n");
  init<DOF>(user, NVal, p);

  //initial conditions
  ierr = IGACreateVec(user.iga,user.U);CHKERRQ(ierr);
  ierr = IGACreateVec(user.iga,user.U0);CHKERRQ(ierr);
#if RESTART_IT==0
  ierr = FormInitialCondition(user.iga, *user.U0, &user);CHKERRQ(ierr);
#else
  char           filename[256];
  sprintf(filename,"outU%d.dat",RESTART_IT);  
  ierr = IGAReadVec(user.iga,*user.U0,filename);CHKERRQ(ierr); //Read in vector to restart at step RESTART_IT
#endif 
  ierr = VecCopy(*user.U0, *user.U);CHKERRQ(ierr);

	return 0;
}

#endif
