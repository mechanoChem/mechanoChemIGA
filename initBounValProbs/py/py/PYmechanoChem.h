#include "setup_simulator.h"
#define DIMS 2

class PYmechanoChem
{
 public:
  PYmechanoChem();
	
  void setup_mechanoChemIGA();
  void simulate();
  void setParameter();
  //
  std::shared_ptr<mechanoChemIGA<DIMS> > IGAsimulator;
  int argc=0;
  char **argv;
};

PYmechanoChem::PYmechanoChem()
{
  std::cout<<"PYmechanoChem initialated"<<std::endl;
}
	

void PYmechanoChem::setup_mechanoChemIGA()
{
  std::cout<<"setup_mechanoChemIGA"<<std::endl;
  setup_simulator(IGAsimulator,argc,argv);
  IGAsimulator->load_parameters();
  IGAsimulator->apply_boundary_conditions();
  IGAsimulator->apply_initial_conditions();
}

void PYmechanoChem::simulate()
{
  IGAsimulator->solve_ibvp();
IGAsimulator->output_results();
}

void PYmechanoChem::setParameter()
{
  IGAsimulator->load_parameters();
  IGAsimulator->apply_boundary_conditions();
  IGAsimulator->apply_initial_conditions();
}
