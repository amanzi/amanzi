/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATORS_TYPEDEFS_HH_
#define AMANZI_OPERATORS_TYPEDEFS_HH_

#include <boost/array.hpp>

namespace Amanzi {
namespace Operators {

typedef boost::array<double, 2> bc_tuple;
typedef std::pair<double, double> dt_tuple;

} // namespace Operators
} // namespace Amanzi

#endif
