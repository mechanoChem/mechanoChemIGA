//extern "C" {
#include "petiga.h"
//}
#include "IBVPHeaders.h"

int boundaryConditions(AppCtx& user, double scale){
  PetscErrorCode  ierr;

  double dVal=scale*user.uDirichlet*user.GridScale;
  PetscPrintf(PETSC_COMM_WORLD,"  dVal: %12.6e  \n",dVal);
	/*
  //shear BC
  ierr = IGASetBoundaryValue(iga,0,0,1,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,0,1,1,-dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,0,0,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,1,1,0,-dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(iga,2,0,2,0.0);CHKERRQ(ierr);  
	
  //free BC
  ierr = IGASetBoundaryValue(user.iga,0,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,1,0,1,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,2,0,2,0.0);CHKERRQ(ierr);
  ierr = IGASetBoundaryValue(user.iga,0,1,0,dVal);CHKERRQ(ierr);  
	*/
  //fixed BC
  ierr = IGASetBoundaryValue(user.iga,0,0,0,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,0,1,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,0,2,0.0);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,1,1,dVal);CHKERRQ(ierr);  
  ierr = IGASetBoundaryValue(user.iga,0,1,2,dVal);CHKERRQ(ierr);   

  return 0;
}
