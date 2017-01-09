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
    item.set(AmanziMesh::NODE, SCHEMA_DOFS_SCALAR, 1);
    items_.push_back(item);
  }
  if (i & OPERATOR_SCHEMA_DOFS_EDGE) {
    item.set(AmanziMesh::EDGE, SCHEMA_DOFS_SCALAR, 1);
    items_.push_back(item);
  }
  if (i & OPERATOR_SCHEMA_DOFS_FACE) {
    item.set(AmanziMesh::FACE, SCHEMA_DOFS_SCALAR, 1);
    items_.push_back(item);
  }
  if (i & OPERATOR_SCHEMA_DOFS_CELL) {
    item.set(AmanziMesh::CELL, SCHEMA_DOFS_SCALAR, 1);
    items_.push_back(item);
  }
}


/* ******************************************************************
* Compute offsets (starting position of DOF ids).
****************************************************************** */
void Schema::Finalize(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
{
  offset_.clear();

  int m(0);
  for (auto it = items_.begin(); it != items_.end(); ++it) {
    offset_.push_back(m);
    int nent = mesh->num_entities(it->kind, AmanziMesh::OWNED);
    m += nent * it->num;
  }
}


/* ******************************************************************
* Compute local (cell-based) offsets
****************************************************************** */
void Schema::ComputeOffset(int c, Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                           std::vector<int>& offset)
{
  AmanziMesh::Entity_ID_List nodes, edges, faces;

  offset.clear();
  offset.push_back(0);

  int ndofs;
  for (auto it = items_.begin(); it != items_.end(); ++it) {
    if (it->kind == AmanziMesh::NODE) {
      mesh->cell_get_nodes(c, &nodes);
      ndofs = nodes.size();
    }
    else if (it->kind == AmanziMesh::EDGE) {
      mesh->cell_get_nodes(c, &edges);
      ndofs = edges.size();
    }
    else if (it->kind == AmanziMesh::FACE) {
      ndofs = mesh->cell_get_num_faces(c);
    }
    else if (it->kind == AmanziMesh::CELL) {
      ndofs = 1;
    }

    offset.push_back(ndofs * it->num);
  }
}


/* ******************************************************************
* Compatibility: returns old schema
****************************************************************** */
int Schema::OldSchema() const
{
  int i(base_);
  for (auto it = items_.begin(); it != items_.end(); ++it) {
    if (it->kind == AmanziMesh::NODE) {
      i += OPERATOR_SCHEMA_DOFS_NODE; 
    } else if (it->kind == AmanziMesh::EDGE) {
      i += OPERATOR_SCHEMA_DOFS_EDGE; 
    } else if (it->kind == AmanziMesh::FACE) {
      i += OPERATOR_SCHEMA_DOFS_FACE; 
    } else if (it->kind == AmanziMesh::CELL) {
      i += OPERATOR_SCHEMA_DOFS_CELL; 
    }
  }
  return i;
}


/* ******************************************************************
* Returns standard name for geometric location of DOF.
****************************************************************** */
std::string Schema::KindToString(AmanziMesh::Entity_kind kind) const 
{
  if (kind == AmanziMesh::NODE) {
    return "node";
  } else if (kind == AmanziMesh::EDGE) {
    return "edge";
  } else if (kind == AmanziMesh::FACE) {
    return "face";
  } else if (kind == AmanziMesh::CELL) {
    return "cell";
  }
}


/* ******************************************************************
* Returns standard mesh id for geometric location of DOF.
****************************************************************** */
AmanziMesh::Entity_kind Schema::StringToKind(std::string& name) const 
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
* Auxiliary routine creates new name.
****************************************************************** */
std::string Schema::CreateUniqueName() const
{
  std::string name;
  for (auto it = items_.begin(); it != items_.end(); ++it) {
    name.append(KindToString(it->kind)); 
    name.append("+");
  }
  return name;
}

}  // namespace Operators
}  // namespace Amanzi


