/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_FLOW_TYPEDEFS_HH_
#define AMANZI_FLOW_TYPEDEFS_HH_

#include <boost/array.hpp>

namespace Amanzi {
namespace AmanziFlow {

typedef boost::array<double, 2> bc_tuple;
typedef std::pair<double, double> dt_tuple;

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

