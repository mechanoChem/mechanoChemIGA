#include "applicationHeaders.h"
//extern "C" {
#include "petiga.h"
//}
#include "parameters.h"

template<unsigned int dim>
int boundaryConditions(AppCtx& user, double scale){
  PetscErrorCode  ierr;

  double dVal=scale*uDirichlet*GridScale;
  PetscPrintf(PETSC_COMM_WORLD,"  dVal: %12.6e  \n",dVal);

  if(dim == 2){
    //Plane strain bending 
    //Configurational displacements
    ierr = IGASetBoundaryValue(user.iga,0,0,0,0.0);CHKERRQ(ierr);  
    ierr = IGASetBoundaryValue(user.iga,0,0,1,0.0);CHKERRQ(ierr);  

    ierr = IGASetBoundaryValue(user.iga,0,1,0,0.0);CHKERRQ(ierr);  
    ierr = IGASetBoundaryValue(user.iga,0,1,1,-0.5*dVal);CHKERRQ(ierr);  

    //Total displacements
    ierr = IGASetBoundaryValue(user.iga,0,0,2,0.0);CHKERRQ(ierr);  
    ierr = IGASetBoundaryValue(user.iga,0,0,3,0.0);CHKERRQ(ierr);

    ierr = IGASetBoundaryValue(user.iga,0,1,2,0.0);CHKERRQ(ierr);
    ierr = IGASetBoundaryValue(user.iga,0,1,3,-dVal);CHKERRQ(ierr); 
  }

  if(dim == 3){
    //Plane strain bending 
    //Configurational displacements
    ierr = IGASetBoundaryValue(user.iga,2,0,2,0.0);CHKERRQ(ierr);  
    ierr = IGASetBoundaryValue(user.iga,2,1,2,0.0);CHKERRQ(ierr); 

    ierr = IGASetBoundaryValue(user.iga,0,0,0,0.0);CHKERRQ(ierr);  
    ierr = IGASetBoundaryValue(user.iga,0,0,1,0.0);CHKERRQ(ierr);  
    ierr = IGASetBoundaryValue(user.iga,0,0,2,0.0);CHKERRQ(ierr); 

    ierr = IGASetBoundaryValue(user.iga,0,1,0,0.0);CHKERRQ(ierr);  
    ierr = IGASetBoundaryValue(user.iga,0,1,1,-0.5*dVal);CHKERRQ(ierr);  
    ierr = IGASetBoundaryValue(user.iga,0,1,2,0.0);CHKERRQ(ierr);

    //Total displacements
    ierr = IGASetBoundaryValue(user.iga,2,0,5,0.0);CHKERRQ(ierr);
    ierr = IGASetBoundaryValue(user.iga,2,1,5,0.0);CHKERRQ(ierr); 

    ierr = IGASetBoundaryValue(user.iga,0,0,3,0.0);CHKERRQ(ierr);  
    ierr = IGASetBoundaryValue(user.iga,0,0,4,0.0);CHKERRQ(ierr);
    ierr = IGASetBoundaryValue(user.iga,0,0,5,0.0);CHKERRQ(ierr);  

    ierr = IGASetBoundaryValue(user.iga,0,1,3,0.0);CHKERRQ(ierr);
    ierr = IGASetBoundaryValue(user.iga,0,1,4,-dVal);CHKERRQ(ierr);
    ierr = IGASetBoundaryValue(user.iga,0,1,5,0.0);CHKERRQ(ierr);   
  }

  return 0;
}

template int boundaryConditions<2>(AppCtx& user, double scale);
template int boundaryConditions<3>(AppCtx& user, double scale);
