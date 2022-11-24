/*
   State

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Some basic typedefs for State and company.
*/

#ifndef AMANZI_STATE_DEFS_HH_
#define AMANZI_STATE_DEFS_HH_

#include <string>
#include <vector>
#include <set>

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "Tag.hh"

namespace Amanzi {

typedef bool NullFactory; // placeholder object for no factory required

namespace Tags {
static const Tag DEFAULT("");
static const Tag CURRENT("current");
static const Tag INTER("inter");
static const Tag NEXT(""); // an alias used by ATS
static const Tag COPY("copy");
} // namespace Tags

enum class EvaluatorType { PRIMARY, SECONDARY, INDEPENDENT, OTHER };

} // namespace Amanzi

#endif
