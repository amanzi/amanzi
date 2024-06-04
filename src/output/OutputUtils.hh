/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! OutputUtils: Utility functions for I/O.
#ifndef AMANZI_OUTPUT_UTILS_HH
#define AMANZI_OUTPUT_UTILS_HH

#include "AmanziTypes.hh"

namespace Amanzi {
namespace OutputUtils {

// helper for writing maps
Vector_type_<GO>
asVector(const Map_ptr_type& map);

std::vector<std::string>
getNames(const Teuchos::ParameterList& attrs, std::size_t count);

} // namespace OutputUtils
} // namespace Amanzi


#endif
