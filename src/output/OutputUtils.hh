//! OutputUtils: Utility functions for I/O.
/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

#ifndef AMANZI_OUTPUT_UTILS_HH
#define AMANZI_OUTPUT_UTILS_HH

#include "AmanziTypes.hh"

namespace Amanzi {

Map_ptr_type GetNaturalMap(const Map_ptr_type& ghosted_map, const Map_ptr_type& owned_map);


} // namespace Amanzi


#endif
