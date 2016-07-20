//extern "C" {
#include "petiga.h"
//}
#include "applicationHeaders.h"
#include "physicsHeaders.h"

template<unsigned int dim>
struct Field;

template<>
struct Field<2>{
  PetscReal Ux, Uy;
  PetscReal ux, uy;
};

template<>
struct Field<3>{
  PetscReal Ux, Uy, Uz;
  PetscReal ux, uy, uz;
};

#undef  __FUNCT__
#define __FUNCT__ "FormInitialCondition"
PetscErrorCode FormInitialCondition(IGA iga, Vec U, AppCtx *user)
{	
  PetscErrorCode ierr;
  PetscFunctionBegin;
  DM da;
  ierr = IGACreateNodeDM(iga,2*user->dim,&da);CHKERRQ(ierr);
  //ierr = IGACreateNodeDM(iga,DIM,&da);CHKERRQ(ierr);

	if(user->dim==2){
		Field<2> **u;
		ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
		DMDALocalInfo info;
		ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
		PetscInt i,j;
		for(i=info.xs;i<info.xs+info.xm;i++){
		  for(j=info.ys;j<info.ys+info.ym;j++){
		    u[j][i].Ux=0.0;
		    u[j][i].Uy=0.0;
		    u[j][i].ux=0.0;
		    u[j][i].uy=0.0;
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
  	ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
	}

	ierr = DMDestroy(&da);;CHKERRQ(ierr); 

  PetscFunctionReturn(0); 
}
