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

/* ******************************************************************
* Constructor takes the old schema as input.
****************************************************************** */
void Schema::Init(int i) { 
  base_ = i & (OPERATOR_SCHEMA_BASE_NODE + OPERATOR_SCHEMA_BASE_EDGE 
             + OPERATOR_SCHEMA_BASE_FACE + OPERATOR_SCHEMA_BASE_CELL);

  items_.clear();

  SchemaItem item;
  if (i & OPERATOR_SCHEMA_DOFS_NODE) {
    item.set(OPERATOR_SCHEMA_DOFS_NODE, SCHEMA_DOFS_SCALAR, 1);
    items_.push_back(item);
  }
  if (i & OPERATOR_SCHEMA_DOFS_EDGE) {
    item.set(OPERATOR_SCHEMA_DOFS_EDGE, SCHEMA_DOFS_SCALAR, 1);
    items_.push_back(item);
  }
  if (i & OPERATOR_SCHEMA_DOFS_FACE) {
    item.set(OPERATOR_SCHEMA_DOFS_FACE, SCHEMA_DOFS_SCALAR, 1);
    items_.push_back(item);
  }
  if (i & OPERATOR_SCHEMA_DOFS_CELL) {
    item.set(OPERATOR_SCHEMA_DOFS_CELL, SCHEMA_DOFS_SCALAR, 1);
    items_.push_back(item);
  }
}


/* ******************************************************************
* Compatibility: returns old schema
****************************************************************** */
int Schema::OldSchema() const
{
  int i(base_);
  for (auto it = items_.begin(); it != items_.end(); ++it) i += it->location; 
  return i;
}


/* ******************************************************************
* Returns standard name for geometric location of DOF.
****************************************************************** */
std::string Schema::LocationName(int loc) const 
{
  if (loc == OPERATOR_SCHEMA_DOFS_NODE) {
    return "node";
  } else if (loc == OPERATOR_SCHEMA_DOFS_EDGE) {
    return "edge";
  } else if (loc == OPERATOR_SCHEMA_DOFS_FACE) {
    return "face";
  } else if (loc == OPERATOR_SCHEMA_DOFS_CELL) {
    return "cell";
  }
}


/* ******************************************************************
* Returns standard mesh id for geometric location of DOF.
****************************************************************** */
AmanziMesh::Entity_kind Schema::LocationMeshID(int loc) const 
{
  if (loc == OPERATOR_SCHEMA_DOFS_NODE) {
    return AmanziMesh::NODE;
  } else if (loc == OPERATOR_SCHEMA_DOFS_EDGE) {
    return AmanziMesh::EDGE;
  } else if (loc == OPERATOR_SCHEMA_DOFS_FACE) {
    return AmanziMesh::FACE;
  } else if (loc == OPERATOR_SCHEMA_DOFS_CELL) {
    return AmanziMesh::CELL;
  }
}


/* ******************************************************************
* Auxiliary routine creates new name.
****************************************************************** */
std::string Schema::CreateUniqueName() const
{
  std::string name;
  for (auto it = items_.begin(); it != items_.end(); ++it) {
    name.append(LocationName(it->location)); 
    name.append("+");
  }
  return name;
}

}  // namespace Operators
}  // namespace Amanzi


