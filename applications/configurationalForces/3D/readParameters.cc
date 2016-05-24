#include "applicationHeaders.h"
//extern "C" {
#include "petiga.h"
//}
#include<fstream>
#include<string>
#include<map>

int readParameters(AppCtx& user){
	std::fstream file;
	std::string tempVar;
	double tempVal;	
	std::map<std::string,double> parameters;

	file.open("parameters.txt");
	while(file.good()){
		if(file.peek() != '#' && file.peek() != '\n'){
			file >> tempVar >> tempVal;
			parameters[tempVar] = tempVal;
		}
		file.ignore(256,'\n');
		if(file.fail() && !file.eof()){
			PetscPrintf(PETSC_COMM_WORLD,"Error: Input file format violated. Should be a variable name followed by number. Comments designated by '#'.\n");
			return(1);
		}
	}
  file.close();

	user.DIM = (PetscInt) parameters["DIM"];
	user.GridScale = (PetscReal) parameters["GridScale"];
	user.ADSacado = (PetscBool) parameters["ADSacado"];
  user.numVars = (PetscInt) parameters["numVars"];
	user.uDirichlet = (PetscReal) parameters["uDirichlet"];
	user.NVal = (PetscInt) parameters["NVal"];
	user.dtVal = (PetscReal) parameters["dtVal"];
	user.skipOutput = (PetscInt) parameters["skipOutput"];
	user.RESTART_IT = (PetscInt) parameters["RESTART_IT"];
	user.RESTART_TIME = (PetscReal) parameters["RESTART_TIME"];

}
