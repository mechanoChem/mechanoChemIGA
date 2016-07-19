#include "../applicationHeaders.h"
//extern "C" {
#include "petiga.h"
//}

template<unsigned int dim>
int boundaryConditions(AppCtx& user, double scale){
  PetscErrorCode  ierr;

	double t = scale, scale1, scale2;
	if(t<1){
		scale1 = scale;
		scale2 = 0.;
	}
	else{
		scale1 = 1;
		scale2 = t-1+user.dt;
	}
  double dVal=scale2*user.uDirichlet*user.GridScale; //In this case, this is just the applied standard displacmenet

  PetscPrintf(PETSC_COMM_WORLD,"  dVal: %12.6e  \n",dVal);
	FILE	*output_file = NULL;
	PetscFOpen(PETSC_COMM_WORLD,"stress_stretch.txt","a+",&output_file);
	PetscFPrintf(PETSC_COMM_WORLD,output_file,"%g ",dVal);
	PetscFClose(PETSC_COMM_WORLD,output_file);

  //Configurational displacements (forced into e1 tetragonal)
  ierr = IGASetBoundaryValue(user.iga,0,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,0,1,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,2,0,2,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,1,1,1,(-1. + sqrt(1. + sqrt(8./3.)*user.Es))*scale1);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,0,1,0,(-1. + sqrt(1. - sqrt(2./3.)*user.Es))*scale1);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,2,1,2,(-1. + sqrt(1. - sqrt(2./3.)*user.Es))*scale1);CHKERRQ(ierr);

  //Total displacements (configurational displacement + uniaxial tension)
  ierr = IGASetBoundaryValue(user.iga,0,0,3,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,0,4,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,2,0,5,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,0,1,3,dVal + (-1. + sqrt(1. - sqrt(2./3.)*user.Es))*scale1);CHKERRQ(ierr);  

  return 0;
}

template int boundaryConditions<2>(AppCtx& user, double scale);
template int boundaryConditions<3>(AppCtx& user, double scale);
