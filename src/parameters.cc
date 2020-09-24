#include <map>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include "coreFunctions.h"
#include "tensor.h"
#include "appCtx.h"
#include "json.hpp"

int ReadDim(){
  int dim = 0;

  nlohmann::json j_tmp;

  std::fstream file("parameters.json");
  if (!file){
    file.close();
    file.open("../parameters.json");
  }
  if(file){
    file >> j_tmp;
    file.close();

    if (j_tmp.contains("dim")){
      dim = j_tmp["dim"].get<int>();
    }
  }
  return dim;
}

template<typename T>
void copy_if_contains(nlohmann::json j, std::string name, T &val){
  if (j.contains(name)){
    j[name].get_to(val);
  }
}

template<int dim>
int ReadParameters(AppCtx<dim> &user){

  std::fstream file("parameters.json");
  if (!file){
    file.close();
    file.open("../parameters.json");
  }
  if(file){
    nlohmann::json j_tmp;
    file >> j_tmp;
    user.param.update(j_tmp); //Update values from parameters file
    file.close();

    copy_if_contains(user.param,"Dtval",user.dtVal);
    copy_if_contains(user.param,"totalTime",user.totalTime);
    copy_if_contains(user.param,"RESTART_TIME",user.RESTART_TIME);
    copy_if_contains(user.param,"maxTimeSteps",user.maxTimeSteps);
    copy_if_contains(user.param,"RESTART_IT",user.RESTART_IT);
    copy_if_contains(user.param,"skipOutput",user.skipOutput);
    copy_if_contains(user.param,"polyOrder",user.polyOrder);
    copy_if_contains(user.param,"globalContinuity",user.globalContinuity);
    copy_if_contains(user.param,"adapTS",user.adapTS);
    copy_if_contains(user.param,"outputDir",user.outputDir);
    copy_if_contains(user.param,"Elasticity",user.Elasticity);
    copy_if_contains(user.param,"CahnHilliard",user.CahnHilliard);
    copy_if_contains(user.param,"Diffusion",user.Diffusion);
    copy_if_contains(user.param,"adapTS",user.adapTS);
    if (user.param.contains("N")){
      for (unsigned int i = 0; i<dim; ++i){
	user.param["N"][i].get_to(user.N[i]);
      }
    }
    if (user.param.contains("L")){
      for (unsigned int i = 0; i<dim; ++i){
	user.param["L"][i].get_to(user.L[i]);
      }
    }
  }

  return 0;
}

template int ReadParameters<1>(AppCtx<1> &user);
template int ReadParameters<2>(AppCtx<2> &user);
template int ReadParameters<3>(AppCtx<3> &user);
