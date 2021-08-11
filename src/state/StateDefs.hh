/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Some basic typedefs for State and company.
   ------------------------------------------------------------------------- */

#ifndef AMANZI_STATE_DEFS_HH_
#define AMANZI_STATE_DEFS_HH_

#include <string>
#include <vector>
#include <set>

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"

namespace Amanzi {

// Fields
typedef enum { NULL_FIELD_TYPE, COMPOSITE_VECTOR_FIELD, CONSTANT_VECTOR, CONSTANT_SCALAR } FieldType;

typedef enum { UNKNOWN, PRIMARY, SECONDARY, INDEPENDENT } EvaluatorType;

} // namespace

#endif
