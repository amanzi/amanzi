/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <algorithm>
#include <iostream>
#include <vector>

#include "OperatorDefs.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
 * Constructor takes the old schema as input.
 ****************************************************************** */
void
Schema::Init(Teuchos::ParameterList& plist,
             Teuchos::RCP<const AmanziMesh::Mesh> mesh)
{
  Errors::Message msg;

  if (plist.isParameter("base")) {
    base_ = StringToKind(plist.get<std::string>("base"));
  } else {
    msg << "Parameter schema->base is missing.";
    Exceptions::amanzi_throw(msg);
  }

  std::vector<std::string> name;
  if (plist.isParameter("location")) {
    name = plist.get<Teuchos::Array<std::string>>("location").toVector();
  } else {
    msg << "Parameter schema->location is missing.";
    Exceptions::amanzi_throw(msg);
  }

  std::vector<std::string> type;
  if (plist.isParameter("type")) {
    type = plist.get<Teuchos::Array<std::string>>("type").toVector();
  } else {
    msg << "Parameter schema->type is missing.";
    Exceptions::amanzi_throw(msg);
  }

  std::vector<int> ndofs;
  if (plist.isParameter("number")) {
    ndofs = plist.get<Teuchos::Array<int>>("number").toVector();
  } else {
    msg << "Parameter schema->number is missing.";
    Exceptions::amanzi_throw(msg);
  }

  // Populate schema and save it.
  items_.clear();
  for (int i = 0; i < name.size(); i++) {
    AddItem(StringToKind(name[i]), StringToType(type[i]), ndofs[i]);
  }

  Finalize(mesh);
}


/* ******************************************************************
 * Backward compatibility: takes the old schema as input.
 ****************************************************************** */
void
Schema::Init(int i)
{
  base_ = AmanziMesh::CELL; // default

  if (i & OPERATOR_SCHEMA_BASE_NODE) {
    base_ = AmanziMesh::NODE;
  } else if (i & OPERATOR_SCHEMA_BASE_EDGE) {
    base_ = AmanziMesh::EDGE;
  } else if (i & OPERATOR_SCHEMA_BASE_FACE) {
    base_ = AmanziMesh::FACE;
  } else if (i & OPERATOR_SCHEMA_BASE_CELL) {
    base_ = AmanziMesh::CELL;
  }

  items_.clear();

  SchemaItem item;
  if (i & OPERATOR_SCHEMA_DOFS_NODE) {
    item.set(AmanziMesh::NODE, DOF_Type::SCALAR, 1);
    items_.push_back(item);
  }
  if (i & OPERATOR_SCHEMA_DOFS_EDGE) {
    item.set(AmanziMesh::EDGE, DOF_Type::SCALAR, 1);
    items_.push_back(item);
  }
  if (i & OPERATOR_SCHEMA_DOFS_FACE) {
    item.set(AmanziMesh::FACE, DOF_Type::SCALAR, 1);
    items_.push_back(item);
  }
  if (i & OPERATOR_SCHEMA_DOFS_CELL) {
    item.set(AmanziMesh::CELL, DOF_Type::SCALAR, 1);
    items_.push_back(item);
  }
  if (i & OPERATOR_SCHEMA_DOFS_BNDFACE) {
    item.set(AmanziMesh::BOUNDARY_FACE, DOF_Type::SCALAR, 1);
    items_.push_back(item);
  }
}


/* ******************************************************************
 * Backward compatibility: takes kind as input.
 ****************************************************************** */
void
Schema::Init(AmanziMesh::Entity_kind kind, int nvec)
{
  base_ = kind;

  SchemaItem item;
  item.set(kind, DOF_Type::SCALAR, nvec);

  items_.clear();
  items_.push_back(item);
}


/* ******************************************************************
 * Compute offsets (starting position of DOF ids).
 ****************************************************************** */
void
Schema::Finalize(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
{
  offset_.clear();

  int m(0);
  for (auto it = items_.begin(); it != items_.end(); ++it) {
    offset_.push_back(m);
    int nent = mesh->num_entities(it->kind, AmanziMesh::Parallel_type::OWNED);
    m += nent * it->num;
  }
}


/* ******************************************************************
 * Compute local (cell-based) offsets
 ****************************************************************** */
void
Schema::ComputeOffset(int c, Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                      std::vector<int>& offset) const
{
  AmanziMesh::Entity_ID_List nodes, edges, faces;

  offset.clear();
  offset.push_back(0);

  int ndofs;
  for (auto it = items_.begin(); it != items_.end(); ++it) {
    if (it->kind == AmanziMesh::NODE) {
      mesh->cell_get_nodes(c, &nodes);
      ndofs = nodes.size();
    } else if (it->kind == AmanziMesh::EDGE) {
      mesh->cell_get_nodes(c, &edges);
      ndofs = edges.size();
    } else if (it->kind == AmanziMesh::FACE) {
      ndofs = mesh->cell_get_num_faces(c);
    } else if (it->kind == AmanziMesh::CELL) {
      ndofs = 1;
    }

    offset.push_back(ndofs * it->num);
  }
}


/* ******************************************************************
 * Compatibility: returns old schema
 ****************************************************************** */
int
Schema::OldSchema() const
{
  int i(0);

  // convert base
  if (base_ == AmanziMesh::NODE) {
    i = OPERATOR_SCHEMA_BASE_NODE;
  } else if (base_ == AmanziMesh::EDGE) {
    i = OPERATOR_SCHEMA_BASE_EDGE;
  } else if (base_ == AmanziMesh::FACE) {
    i = OPERATOR_SCHEMA_BASE_FACE;
  } else if (base_ == AmanziMesh::CELL) {
    i = OPERATOR_SCHEMA_BASE_CELL;
  }

  for (auto it = items_.begin(); it != items_.end(); ++it) {
    if (it->kind == AmanziMesh::NODE) {
      i += OPERATOR_SCHEMA_DOFS_NODE;
    } else if (it->kind == AmanziMesh::EDGE) {
      i += OPERATOR_SCHEMA_DOFS_EDGE;
    } else if (it->kind == AmanziMesh::FACE) {
      i += OPERATOR_SCHEMA_DOFS_FACE;
    } else if (it->kind == AmanziMesh::CELL) {
      i += OPERATOR_SCHEMA_DOFS_CELL;
    } else if (it->kind == AmanziMesh::BOUNDARY_FACE) {
      i += OPERATOR_SCHEMA_DOFS_BNDFACE;
    }
  }
  return i;
}


/* ******************************************************************
 * Returns standard name for geometric location of DOF.
 ****************************************************************** */
std::string
Schema::KindToString(AmanziMesh::Entity_kind kind) const
{
  if (kind == AmanziMesh::NODE) {
    return "node";
  } else if (kind == AmanziMesh::EDGE) {
    return "edge";
  } else if (kind == AmanziMesh::FACE) {
    return "face";
  } else if (kind == AmanziMesh::CELL) {
    return "cell";
  } else if (kind == AmanziMesh::BOUNDARY_FACE) {
    return "boundary_face";
  }
  return "null";
}


/* ******************************************************************
 * Returns standard mesh id for geometric location of DOF.
 ****************************************************************** */
AmanziMesh::Entity_kind
Schema::StringToKind(std::string& name) const
{
  if (name == "node") {
    return AmanziMesh::NODE;
  } else if (name == "edge") {
    return AmanziMesh::EDGE;
  } else if (name == "face") {
    return AmanziMesh::FACE;
  } else if (name == "cell") {
    return AmanziMesh::CELL;
  }
}


/* ******************************************************************
 * Returns standard mesh id for geometric location of DOF.
 ****************************************************************** */
DOF_Type
Schema::StringToType(std::string& name) const
{
  if (name == "scalar") {
    return DOF_Type::SCALAR;
  } else if (name == "vector") {
    return DOF_Type::VECTOR;
  } else if (name == "point") {
    return DOF_Type::POINT;
  } else if (name == "normal component") {
    return DOF_Type::NORMAL_COMPONENT;
  } else if (name == "moment") {
    return DOF_Type::MOMENT;
  }
  return DOF_Type::SCALAR;
}


/* ******************************************************************
 * Auxiliary routine creates new name.
 ****************************************************************** */
std::string
Schema::CreateUniqueName() const
{
  std::string name(KindToString(base_)), c("_");
  for (auto it = items_.begin(); it != items_.end(); ++it) {
    name.append(c);
    name.append(KindToString(it->kind));
    name.append(std::to_string(it->num));
    c = "+";
  }

  std::transform(name.begin(), name.end(), name.begin(), ::toupper);
  return name;
}

} // namespace Operators
} // namespace Amanzi
