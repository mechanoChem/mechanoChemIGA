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
    nlohmann::json j_tmp, param_flat;
    file >> j_tmp;
    param_flat = user.param.flatten();
    param_flat.update(j_tmp.flatten()); //Update values from parameters file
    user.param = param_flat.unflatten();

    file.close();

    copy_if_contains(user.param,"dt",user.dtVal);
    copy_if_contains(user.param,"max_time",user.totalTime);
    copy_if_contains(user.param,"restart_time",user.RESTART_TIME);
    copy_if_contains(user.param,"max_timesteps",user.maxTimeSteps);
    copy_if_contains(user.param,"restart_it",user.RESTART_IT);
    copy_if_contains(user.param,"skip_output",user.skipOutput);
    copy_if_contains(user.param,"poly_order",user.polyOrder);
    copy_if_contains(user.param,"global_continuity",user.globalContinuity);
    copy_if_contains(user.param,"adaptive_ts",user.adapTS);
    copy_if_contains(user.param,"output_dir",user.outputDir);
    if (user.param.contains("mesh")){
      for (unsigned int i = 0; i<dim; ++i){
	user.param["mesh"][i].get_to(user.N[i]);
      }
    }
    if (user.param.contains("box_geom")){
      for (unsigned int i = 0; i<dim; ++i){
	user.param["box_geom"][i].get_to(user.L[i]);
      }
    }
  }

  return 0;
}

template int ReadParameters<1>(AppCtx<1> &user);
template int ReadParameters<2>(AppCtx<2> &user);
template int ReadParameters<3>(AppCtx<3> &user);
