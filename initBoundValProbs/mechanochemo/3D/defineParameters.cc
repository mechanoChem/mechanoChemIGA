//extern "C" {
#include "petiga.h"
//}

#include "IBVPHeaders.h"
//include automatic differentiation library
#include <Sacado.hpp>

int defineParameters(AppCtx& user){

  PetscErrorCode ierr;
 
	user.dim = 3;
	user.GridScale = 1.;
	user.ADSacado = PETSC_TRUE;
	user.uDirichlet = 0.01;
	user.Nx = 80;
	user.Ny = 80;
	user.Nz = 80;
	user.Lx = 1.;
	user.Ly = 1.;
	user.Lz = 1.;
	user.dtVal = 1.0e-5;
	user.skipOutput = 1;
	user.RESTART_IT = 0;
	user.RESTART_TIME = 0.;

	//Elastic free energy parameters
	user.matParam["Es"] = 0.01;
	user.matParam["Ed"] = -0.001;
	user.matParam["El"] = 0.01;

	user.matParam["Gl"] = 0.0;

	//Chemical free energy parameters
	user.matParam["Cs"] = 1.0;
	user.matParam["Cd"] = -1.0;
	user.matParam["Cl"] = 0.01*user.GridScale*user.GridScale;

	//Influx
	user.matParam["flux_xmin"] = 0.0; //flux in side x_min
	user.matParam["flux_xmax"] = 0.0; //flux in side x_max
	user.matParam["flux_ymin"] = 0.0; //in y_min
	user.matParam["flux_ymax"] = 0.0; //in y_max
	user.matParam["flux_zmin"] = 0.0; //in z_min
	user.matParam["flux_zmax"] = 0.0; //in z_max

	user.matParam["DVal"] = 10.0*user.GridScale*user.GridScale; //Diffusivity
	user.matParam["CVal"] = 5.0;
	user.matParam["gamma"] = 1.0;

  return 0;
}
