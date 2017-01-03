/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
#include <vector>

#include "OperatorDefs.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

void Schema::Init(int i) { 
  base_ = i & (OPERATOR_SCHEMA_BASE_NODE + OPERATOR_SCHEMA_BASE_EDGE 
             + OPERATOR_SCHEMA_BASE_FACE + OPERATOR_SCHEMA_BASE_CELL);

  schema_.clear();

  SchemaItem item;
  if (i & OPERATOR_SCHEMA_DOFS_NODE) {
    item.set(OPERATOR_SCHEMA_DOFS_NODE, SCHEMA_DOFS_SCALAR, 1);
    schema_.push_back(item);
  }
  if (i & OPERATOR_SCHEMA_DOFS_EDGE) {
    item.set(OPERATOR_SCHEMA_DOFS_EDGE, SCHEMA_DOFS_SCALAR, 1);
    schema_.push_back(item);
  }
  if (i & OPERATOR_SCHEMA_DOFS_FACE) {
    item.set(OPERATOR_SCHEMA_DOFS_FACE, SCHEMA_DOFS_SCALAR, 1);
    schema_.push_back(item);
  }
  if (i & OPERATOR_SCHEMA_DOFS_CELL) {
    item.set(OPERATOR_SCHEMA_DOFS_CELL, SCHEMA_DOFS_SCALAR, 1);
    schema_.push_back(item);
  }
}


int Schema::OldSchema() const {
  int i(base_);
  for (auto it = schema_.begin(); it != schema_.end(); ++it) i += it->location; 
  return i;
}

}  // namespace Operators
}  // namespace Amanzi


