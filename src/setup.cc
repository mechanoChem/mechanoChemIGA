
//extern "C" {
#include "petiga.h"
//}
#include "coreFunctions.h"
#include "userFunctions.h"
#include "defaultUserFunctions.h"
#include "petigasnes_mod.h"

PetscErrorCode SNESConvergedTest(SNES snes, PetscInt it,PetscReal xnorm, PetscReal snorm, PetscReal fnorm, SNESConvergedReason *reason, void *ctx)
{
  PetscPrintf(PETSC_COMM_WORLD,"xnorm:%12.6e snorm:%12.6e fnorm:%12.6e\n",xnorm,snorm,fnorm);
  //default test
  PetscFunctionReturn(SNESConvergedDefault(snes,it,xnorm,snorm,fnorm,reason,ctx));
}

template<unsigned int DIM>
int Setup(AppCtx<DIM>& user,Vec *U,Vec *Up,Vec *Upp,SNES &snes){

  PetscErrorCode ierr;

  //Set default user functions (these may be overridden by the defineParameters function
  user.boundaryConditions = boundaryConditions;
  user.scalarInitialConditions = scalarInitialConditions;
  user.vectorInitialConditions = vectorInitialConditions;
  user.loadStep = loadStep;
  user.adaptiveTimeStep = adaptiveTimeStep;
  user.projectFields = projectFields;

  //Set periodicity to be false by default (can be overridden by user)
  for (unsigned int i=0; i<DIM; ++i){
    user.periodic[i] = PETSC_FALSE;
  }

  //Define user parameters and functions
  defineParameters<DIM>(user);
  ReadParameters<DIM>(user);
  user.dt = user.dtVal;
  //Check that the output directory does not end with a "/"
  if (user.outputDir.back()=='/'){
    user.outputDir.pop_back();
  }

  //define number of fields base on vector of names
  PetscInt nScalars = user.scalarSolnFields.size();
  PetscInt nVectors = user.vectorSolnFields.size();
  PetscInt nScalarProjections = user.scalarProjectnFields.size();
  PetscInt nVectorProjections = user.vectorProjectnFields.size();

  //Write out file to be used in postprocessing
  PetscViewer viewer;
  PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
  PetscViewerSetType(viewer, PETSCVIEWERASCII);
  PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
  char filename[256];
  sprintf(filename,"%s/fieldInfo.txt",user.outputDir.c_str(),user.RESTART_IT);
  PetscViewerFileSetName(viewer, filename);
  PetscViewerASCIIPrintf(viewer, "%d",DIM);
  PetscViewerASCIIPrintf(viewer, "\n%d",nScalars);
  for(unsigned int i=0; i<nScalars; ++i){
    PetscViewerASCIIPrintf(viewer, " %s",user.scalarSolnFields[i].c_str());
  }
  PetscViewerASCIIPrintf(viewer, "\n%d",nVectors);
  for(unsigned int i=0; i<nVectors; ++i){
    PetscViewerASCIIPrintf(viewer, " %s",user.vectorSolnFields[i].c_str());
  }
  PetscViewerASCIIPrintf(viewer, "\n%d",nScalarProjections);
  for(unsigned int i=0; i<nScalarProjections; ++i){
    PetscViewerASCIIPrintf(viewer, " %s",user.scalarProjectnFields[i].c_str());
  }
  PetscViewerASCIIPrintf(viewer, "\n%d",nVectorProjections);
  for(unsigned int i=0; i<nVectorProjections; ++i){
    PetscViewerASCIIPrintf(viewer, " %s",user.vectorProjectnFields[i].c_str());
  }
  PetscViewerDestroy(&viewer);

  //application context objects and parameters
  user.snes=&snes;
  user.Upp=Upp;
  user.Up=Up;
  user.U=U;
  IGA iga, igaProject;
  user.iga = iga;
  user.igaProject = igaProject;
  PetscInt DOF = user.scalarSolnFields.size()+ DIM*user.vectorSolnFields.size();
  PetscInt DOFproject = user.scalarProjectnFields.size() + DIM*user.vectorProjectnFields.size();

  PetscPrintf(PETSC_COMM_WORLD,"initializing...\n");
  InitIGA<DIM>(user, user.iga, user.polyOrder, user.globalContinuity, DOF);
  if(nScalarProjections+nVectorProjections > 0){
    InitIGA<DIM>(user, user.igaProject, user.polyOrder, user.globalContinuity, DOFproject);
  }

  //initial conditions
  ierr = IGACreateVec(user.iga,user.U);CHKERRQ(ierr);
  ierr = IGACreateVec(user.iga,user.Up);CHKERRQ(ierr);
  ierr = IGACreateVec(user.iga,user.Upp);CHKERRQ(ierr);

  if(user.RESTART_IT==0){
    ierr = FormInitialCondition(user.iga, *user.Up, &user);CHKERRQ(ierr);

    char meshfilename[256];
    sprintf(meshfilename, "%s/mesh.dat",user.outputDir.c_str());
    PetscPrintf(PETSC_COMM_WORLD,"\nWriting mesh file: %s\n", meshfilename);
    ierr = IGAWrite(user.iga, meshfilename);CHKERRQ(ierr);
  }
  else{
    char filename[256];
    sprintf(filename,"%s/outU%d.dat",user.outputDir.c_str(),user.RESTART_IT);  
    ierr = IGAReadVec(user.iga,*user.Up,filename);CHKERRQ(ierr); //Read in vector to restart at step RESTART_IT
  }

  ierr = VecCopy(*user.Up, *user.U);CHKERRQ(ierr);
  ierr = VecCopy(*user.Up, *user.Upp);CHKERRQ(ierr); //For now, make U=Up=Upp for initial conditions, but shouldn't need to be this way.

  //Dirichlet boundary conditons for mechanics
  PetscPrintf(PETSC_COMM_WORLD,"applying bcs...\n");
  user.boundaryConditions(user,0.);

  //Assign residual and jacobian functions
  user.iga->form->ops->FunCtx= &user;
  user.iga->form->ops->JacCtx= &user;

  //Nonlinear solver
  ierr = IGACreateSNES_mod<DIM>(user.iga,&snes);CHKERRQ(ierr);
  //PetscReal atol, rtol, stol;
  //PetscInt maxit, maxf;
  //ierr = SNESGetTolerances(snes,&atol,&rtol,&stol,&maxit,&maxf);
  //std::cout << atol << " " << rtol << " " << stol << " " << maxit << " " << maxf << std::endl;
  ierr = SNESSetTolerances(snes,1.e-50,1.e-8,1.e-8,10,1000);
  ierr = SNESSetConvergenceTest(snes,SNESConvergedTest,(void*)&user,NULL); CHKERRQ(ierr);
  ierr = SNESSetType(snes,"newtonls");

  KSP ksp;
  ierr = SNESGetKSP(snes,&ksp);
#if defined(PETSC_HAVE_SUPERLU_DIST)
  PC pc;
  ierr = KSPGetPC(ksp,&pc);
  ierr = KSPSetType(ksp,"preonly");
  ierr = PCSetType(pc,"lu");
  ierr = PCFactorSetMatSolverPackage(pc,"superlu_dist");
  //ierr = PCFactorSetMatSolverType(pc,"superlu_dist"); //Changed name with petsc 3.9
  //PetscPrintf(PETSC_COMM_WORLD,"\nUsing superlu_dist\n");
#else
  ierr = KSPSetType(ksp,"gmres");
  //PetscPrintf(PETSC_COMM_WORLD,"\nUsing gmres\n");
#endif
  //MatSolverType stype;
  //PCFactorGetMatSolverType(pc,&stype);
  //PetscPrintf(PETSC_COMM_WORLD,"\nUsing %s\n",stype);

  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  return 0;
}

template int Setup<2>(AppCtx<2>& user,Vec *U,Vec *Up,Vec *Upp,SNES &snes);
template int Setup<3>(AppCtx<3>& user,Vec *U,Vec *Up,Vec *Upp,SNES &snes);
