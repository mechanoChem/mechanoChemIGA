#ifndef initialconditions_
#define initialconditions_

typedef struct {
  PetscReal ux, uy;
  PetscReal c;
#if DIM==3
  PetscReal uz;
#endif
} Field;


PetscErrorCode FormInitialCondition(IGA iga, Vec U, AppCtx *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;
  std::srand(5);
  DM da;
  ierr = IGACreateNodeDM(iga,DIM+1,&da);CHKERRQ(ierr);
#if DIM==2
  Field **u;
#elif DIM==3
  Field ***u;
#endif
  ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
#if FLUX==3 //Quench
  double c_average=0.46;
#else       //Flux 
  double c_average=0.05;
#endif
#if DIM==2
  PetscInt i,j;
  for(i=info.xs;i<info.xs+info.xm;i++){
    for(j=info.ys;j<info.ys+info.ym;j++){
      u[j][i].ux=0.0;
      u[j][i].uy=0.0;
      u[j][i].c= c_average + 0.01*(0.5 - (double)(std::rand() % 100 )/100.0);
    }
  }
#elif DIM==3
  PetscInt i,j,k;
  for(i=info.xs;i<info.xs+info.xm;i++){
    for(j=info.ys;j<info.ys+info.ym;j++){
      for(k=info.zs;k<info.zs+info.zm;k++){
	u[k][j][i].ux=0.0;
	u[k][j][i].uy=0.0;
	u[k][j][i].uz=0.0;
	u[k][j][i].c= c_average + 0.01*(0.5 - (double)(std::rand() % 100 )/100.0);
      }
    }
  }    
#endif
  ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
  ierr = DMDestroy(&da);;CHKERRQ(ierr); 
  PetscFunctionReturn(0); 
}

#endif
