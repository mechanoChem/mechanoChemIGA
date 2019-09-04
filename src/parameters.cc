#include <map>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include "coreFunctions.h"
#include "tensor.h"
#include "appCtx.h"

int ReadDim(){
  int dim = 0;

  std::fstream file("parameters.prm");
  if (!file){
    file.close();
    file.open("../parameters.prm");
  }
  if(file){
    std::string line, temp1;

    while( std::getline(file,line) ){
      if (line.find("=") != -1){
	//Remove equals signs
	line.replace(line.find("="),1," ");
	std::istringstream iss(line);
	iss >> temp1;
	if(temp1=="set"){
	  iss >> temp1;
	  if(temp1 == "dim" || temp1 == "DIM"){
	    iss >> dim;
	  }
	}
      }
    }
    file.close();
  }

  return dim;
}

template<int dim>
int ReadParameters(AppCtx<dim> &user){
  
  std::fstream file("parameters.prm");
  if (!file){
    file.close();
    file.open("../parameters.prm");
  }
  if(file){
    std::string line, temp1;
  
    std::map<std::string,double*> doubles;
    std::map<std::string,int*> ints;
    std::map<std::string,std::string*> strings;
    //std::map<std::string,bool*> bools;
    std::map<std::string,Tensor<1,dim,int>* > tensorsInt;
    std::map<std::string,Tensor<1,dim,double>* > tensorsDouble;
    doubles["dtVal"] = &user.dtVal;
    doubles["totalTime"] = &user.totalTime;
    doubles["RESTART_TIME"] = &user.RESTART_TIME;
    ints["RESTART_IT"] = &user.RESTART_IT;
    ints["skipOutput"] = &user.skipOutput;
    ints["polyOrder"] = &user.polyOrder;
    ints["globalContinuity"] = &user.globalContinuity;
    tensorsInt["N"] = &user.N;
    tensorsDouble["L"] = &user.L;
    strings["outputDir"] = &user.outputDir;

    while( std::getline(file,line) ){
      if (line.find("=") != -1){
	//Remove equals signs
	line.replace(line.find("="),1," ");
	//Remove commas
	while (line.find(",") != -1){
	  line.replace(line.find(","),1," ");
	}
	//Remove brackets
	while (line.find("[") != -1){
	  line.replace(line.find("["),1," ");
	}
	while (line.find("]") != -1){
	  line.replace(line.find("]"),1," ");
	}
	std::istringstream iss(line);
	iss >> temp1;
	if(temp1=="set"){
	  iss >> temp1;
	  if(doubles.count(temp1) == 1){
	    iss >> *doubles[temp1];
	  }
	  else if(ints.count(temp1) == 1){
	    iss >> *ints[temp1];
	  }
	  //else if(bools.count(temp1) == 1){
	  //  iss >> std::boolalpha >> *bools[temp1];
	  //}
	  else if(strings.count(temp1) == 1){
	    iss >> *strings[temp1];
	  }
	  else if(tensorsInt.count(temp1) == 1){
	    //If a tensor, go back and remove "x" from the line
	    while (line.find("x") != -1){
	      line.replace(line.find("x"),1," ");
	    }
	    std::istringstream iss2(line);
	    iss2 >> temp1;
	    iss2 >> temp1;
	    for(unsigned int i=0; i<dim; ++i){
	      iss2 >> (*tensorsInt[temp1])[i];
	    }
	  }
	  else if(tensorsDouble.count(temp1) == 1){
	    //If a tensor, go back and remove "x" from the line
	    while (line.find("x") != -1){
	      line.replace(line.find("x"),1," ");
	    }
	    std::istringstream iss2(line);
	    iss2 >> temp1;
	    iss2 >> temp1;
	    for(unsigned int i=0; i<dim; ++i){
	      iss2 >> (*tensorsDouble[temp1])[i];
	    }
	  }
	  //Otherwise, put it in the matParam map
	  else if(temp1 != "DIM" && temp1 != "dim"){
	    iss >> user.matParam[temp1];
	    //PetscPrintf(PETSC_COMM_WORLD,"%s is not a valid variable.\n",temp1.c_str());
	  }
	}
      }
    }
    file.close();
  }

  return 0;
}

template int ReadParameters<1>(AppCtx<1> &user);
template int ReadParameters<2>(AppCtx<2> &user);
template int ReadParameters<3>(AppCtx<3> &user);
