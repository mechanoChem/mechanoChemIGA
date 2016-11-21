//extern "C" {
#include "petiga.h"
//}
#include "IBVPHeaders.h"
#include "physicsHeaders.h"

template<unsigned int dim>
struct Field;

template<>
struct Field<2>{
  PetscReal ux, uy, c;
};

template<>
struct Field<3>{
  PetscReal ux, uy, uz, c;
};

#undef  __FUNCT__
#define __FUNCT__ "FormInitialCondition"
PetscErrorCode FormInitialCondition(IGA iga, Vec U, AppCtx *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;
  srand(5);
  DM da;
  ierr = IGACreateNodeDM(iga,user->dim+1,&da);CHKERRQ(ierr);
  double c_average=0.05;

  if(user->dim==2){
    Field<2> **u;
    ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
    PetscInt i,j;
    PetscReal jmax = (info.my-1), Ly = user->GridScale;
    for(i=info.xs;i<info.xs+info.xm;i++){
      for(j=info.ys;j<info.ys+info.ym;j++){
	PetscReal y = (j/jmax)*Ly;
      	u[j][i].ux=0.0;
      	u[j][i].uy=0.0;
      	u[j][i].c= 0.498*y/Ly + 0.001; //Uniform gradient initially
      }
    }
    ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr);
  }
  else if(user->dim==3){
    Field<3> ***u;
    ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
    PetscInt i,j,k;
    PetscReal jmax = (info.my-1), Ly = user->GridScale;
    for(i=info.xs;i<info.xs+info.xm;i++){
      for(j=info.ys;j<info.ys+info.ym;j++){
	for(k=info.zs;k<info.zs+info.zm;k++){
	  PetscReal y = (j/jmax)*Ly;
	  u[k][j][i].ux=0.0;
	  u[k][j][i].uy=0.0;
	  u[k][j][i].uz=0.0;
	  u[k][j][i].c= 0.498*y/Ly + 0.001; //Uniform gradient initially
	}
      }
    }    
    ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
  }

  ierr = DMDestroy(&da);CHKERRQ(ierr); 

  PetscFunctionReturn(0); 
}
