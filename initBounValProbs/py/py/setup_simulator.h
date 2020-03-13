#include "mechanoChemIGA.h"
#define DIMS 2
void setup_simulator(std::shared_ptr<mechanoChemIGA<DIMS> >& IGAsimulator,int argc,char **argv)
{
  std::shared_ptr<mechanoChemIGA<DIMS> > simulator_tem(new mechanoChemIGA<DIMS>(argc,argv));	
  IGAsimulator = simulator_tem;	
}
