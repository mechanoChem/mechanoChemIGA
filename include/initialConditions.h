#ifndef initialconditions_
#define initialconditions_

typedef struct {
  PetscReal ux, uy;
  PetscReal c;
} Field;


PetscErrorCode FormInitialCondition2D(IGA iga, Vec U, AppCtx *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;
  std::srand(5);
  DM da;
  ierr = IGACreateNodeDM(iga,3,&da);CHKERRQ(ierr);
  Field **u;
  ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
#if FLUX==3 //Quench
  double c_average=0.45;
#else       //Flux 
  double c_average=0.05;
#endif
  PetscInt i,j;
  for(i=info.xs;i<info.xs+info.xm;i++){
    for(j=info.ys;j<info.ys+info.ym;j++){
      u[j][i].ux=0.0;
      u[j][i].uy=0.0;
      u[j][i].c= c_average + 0.01*(0.5 - (double)(std::rand() % 100 )/100.0);
    }
  }
  ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
  ierr = DMDestroy(&da);;CHKERRQ(ierr); 
  PetscFunctionReturn(0); 
}

#endif
