//extern "C" {
#include "petiga.h"
//}

#include "IBVPHeaders.h"
//include automatic differentiation library
#include <Sacado.hpp>
#include <fstream>

int defineParameters(AppCtx& user){

  PetscErrorCode ierr;
 
  user.dim = 3;
  user.GridScale = 1;
  user.ADSacado = PETSC_TRUE;
  user.uDirichlet = 0;
  user.Nx = 100;
  user.Ny = 100;
  user.Ny = 100;
  user.Lx = 2.;
  user.Ly = 2.;
  user.Lz = 2.;
  user.dtVal = 10.;
  user.skipOutput = 5;
  user.RESTART_IT = 0;
  user.RESTART_TIME = 0.;

  //Elastic free energy parameters
  user.matParam["Es"] = 0.1;
  user.matParam["Ed"] = -0.1;
  user.matParam["El"] = 0.0001;

  user.matParam["Gl"] = 0.0;

  //Chemical free energy parameters
  user.matParam["Cs"] = 1.0; //Not used in this example
  user.matParam["Cd"] = -2.0; //Not used in this example
  user.matParam["Cl"] = 2.5e-4;

  //Influx
  user.matParam["flux_xmin"] = 0.0; //flux in side x_min
  user.matParam["flux_xmax"] = 0.0; //flux in side x_max
  user.matParam["flux_ymin"] = 0.0; //in y_min
  user.matParam["flux_ymax"] = 0.0; //in y_max
  //user.matParam["flux_zmin"] = 0.0; //in z_min
  //user.matParam["flux_zmax"] = 0.0; //in z_max

  user.matParam["DVal"] = 9.e-6; //Diffusivity, um^2/s
  user.matParam["CVal"] = 5.0;
  user.matParam["gamma"] = 1.0;

  user.matParam["splineOrRK"] = 1.; //1 for spline, 2 for RK poly (place coefficient & breaks files in the same directory as defineParameters.cc
  //Note that currently this assumes log terms scaled to be over 0 and 1/2

  //Read in curve fit data
  if(user.matParam["splineOrRK"] == 1 || user.matParam["splineOrRK"] == 2){
    double data;
    std::ifstream file;

    //Load breaks vector
    if(user.matParam["splineOrRK"] == 1){
      file.open("splineBreaks.txt");
      user.breaks.clear();
      while(file >> data){
	user.breaks.push_back(data);
      }
      file.close();
    }

    //Load coefficients matrix
    if(user.matParam["splineOrRK"] == 1){
      file.open("splineCoeffs.txt");
    }
    else{
      file.open("polyCoeffs.txt");
    }
    user.coeffs.clear();
    std::string line;
    while(std::getline(file, line)){
      std::vector<double> row;
      std::istringstream Row(line);
      while (Row >> data)
	{
	  row.push_back(data);
	}
      user.coeffs.push_back(row);
    }
    file.close();
  }

  return 0;
}
