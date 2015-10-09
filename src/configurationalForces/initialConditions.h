#ifndef initialconditions_
#define initialconditions_

typedef struct {
  PetscReal Ux, Uy;
#if DIM==3
  PetscReal Uz;
#endif
  PetscReal ux, uy;
#if DIM==3
  PetscReal uz;
#endif
} Field;


#undef  __FUNCT__
#define __FUNCT__ "FormInitialCondition"
PetscErrorCode FormInitialCondition(IGA iga, Vec U, AppCtx *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;
  DM da;
  ierr = IGACreateNodeDM(iga,2*DIM,&da);CHKERRQ(ierr);
  //ierr = IGACreateNodeDM(iga,DIM,&da);CHKERRQ(ierr);
#if DIM==2
  Field **u;
#elif DIM==3
  Field ***u;
#endif
  ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
  DMDALocalInfo info;
  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
  //
#if DIM==2
  PetscInt i,j;
  for(i=info.xs;i<info.xs+info.xm;i++){
    for(j=info.ys;j<info.ys+info.ym;j++){
      u[j][i].Ux=0.0;
      u[j][i].Uy=0.0;
      u[j][i].ux=0.0;
      u[j][i].uy=0.0;
    }
  }
#elif DIM==3
  PetscInt i,j,k;
  for(i=info.xs;i<info.xs+info.xm;i++){
    for(j=info.ys;j<info.ys+info.ym;j++){
      for(k=info.zs;k<info.zs+info.zm;k++){
	u[k][j][i].Ux=0.0;
	u[k][j][i].Uy=0.0;
	u[k][j][i].Uz=0.0;
	u[k][j][i].ux=0.0;
	u[k][j][i].uy=0.0;
	u[k][j][i].uz=0.0;
      }
    }
  }    
#endif
  ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
  ierr = DMDestroy(&da);;CHKERRQ(ierr); 
  PetscFunctionReturn(0); 
}

#endif
