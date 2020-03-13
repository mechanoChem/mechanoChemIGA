#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


#include "PYmechanoChem.h"

namespace py = pybind11;

PYBIND11_MODULE(PYmechanoChem,m)
{
  // Declares class Simulation, specifying the constructor, the greet function
  // and declaring msg as an attribute visible from python, which can be queried
	// and set through the C++ get and set functions.
  py::class_< PYmechanoChem >(m,"PYmechanoChem")
          .def(py::init())
	  .def("setup_mechanoChemIGA", &PYmechanoChem::setup_mechanoChemIGA)
	  .def("simulate", &PYmechanoChem::simulate)
	  .def("setParameter", &PYmechanoChem::setParameter);
}
