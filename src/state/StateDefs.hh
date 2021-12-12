/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
   State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
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

namespace Amanzi {

// Fields
typedef enum { NULL_FIELD_TYPE, COMPOSITE_VECTOR_FIELD, CONSTANT_VECTOR, CONSTANT_SCALAR } FieldType;

enum class DataType { NULL_TYPE,
                      COMPOSITE_VECTOR,
                      CONSTANT_DOUBLE_VECTOR,
                      CONSTANT_DOUBLE };

typedef enum { UNKNOWN, PRIMARY, SECONDARY, INDEPENDENT } EvaluatorType;

typedef bool NullFactory; // placeholder object for no factory required

} // namespace

#endif
