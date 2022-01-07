/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <iostream>
#include <vector>

#include "OperatorDefs.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Constructor from a bilinear form
****************************************************************** */
void Schema::Init(Teuchos::RCP<const WhetStone::BilinearForm> form, 
                  Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                  AmanziMesh::Entity_kind base)
{ 
  base_ = base;
  items_ = form->schema();
  Finalize(mesh);
}


/* ******************************************************************
* Backward compatibility: takes the old schema as input.
****************************************************************** */
void Schema::Init(int i)
{ 
  base_ = AmanziMesh::Entity_kind::CELL;  // default

  if (i & OPERATOR_SCHEMA_BASE_NODE) {
    base_ = AmanziMesh::Entity_kind::NODE;
  } else if (i & OPERATOR_SCHEMA_BASE_EDGE) {
    base_ = AmanziMesh::Entity_kind::EDGE;
  } else if (i & OPERATOR_SCHEMA_BASE_FACE) {
    base_ = AmanziMesh::Entity_kind::FACE;
  } else if (i & OPERATOR_SCHEMA_BASE_CELL) {
    base_ = AmanziMesh::Entity_kind::CELL;
  }

  items_.clear();

  if (i & OPERATOR_SCHEMA_DOFS_NODE) {
    items_.push_back(std::make_tuple(AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::SCALAR, 1));
  }
  if (i & OPERATOR_SCHEMA_DOFS_EDGE) {
    items_.push_back(std::make_tuple(AmanziMesh::Entity_kind::EDGE, WhetStone::DOF_Type::SCALAR, 1));
  }
  if (i & OPERATOR_SCHEMA_DOFS_FACE) {
    items_.push_back(std::make_tuple(AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::SCALAR, 1));
  }
  if (i & OPERATOR_SCHEMA_DOFS_CELL) {
    items_.push_back(std::make_tuple(AmanziMesh::Entity_kind::CELL, WhetStone::DOF_Type::SCALAR, 1));
  }
  if (i & OPERATOR_SCHEMA_DOFS_BNDFACE) {
    items_.push_back(std::make_tuple(AmanziMesh::Entity_kind::BOUNDARY_FACE, WhetStone::DOF_Type::SCALAR, 1));
  }
}


/* ******************************************************************
* Backward compatibility: takes kind as input.
****************************************************************** */
void Schema::Init(AmanziMesh::Entity_kind kind, int nvec) 
{
  base_ = kind;

  items_.clear();
  items_.push_back(std::make_tuple(kind, WhetStone::DOF_Type::SCALAR, nvec));
}


/* ******************************************************************
* Compute offsets (starting position of DOF ids).
****************************************************************** */
void Schema::Finalize(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
{
  offset_.clear();

  int m(0);
  for (auto it = items_.begin(); it != items_.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    offset_.push_back(m);
    int nent = mesh->getNumEntities(kind, AmanziMesh::Parallel_type::OWNED);
    m += nent * num;
  }
}


/* ******************************************************************
* Compute local (cell-based) offsets
****************************************************************** */
void Schema::ComputeOffset(int c, Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                           std::vector<int>& offset) const
{
  offset.clear();
  offset.push_back(0);

  int ndofs;
  for (auto it = items_.begin(); it != items_.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    if (kind == AmanziMesh::Entity_kind::NODE) {
      AmanziMesh::Entity_ID_List nodes;
      mesh->getCellNodes(c, nodes);
      ndofs = nodes.size();
    }
    else if (kind == AmanziMesh::Entity_kind::EDGE) {
      const auto& edges = mesh->getCellEdges(c);
      ndofs = edges.size();
    }
    else if (kind == AmanziMesh::Entity_kind::FACE) {
      ndofs = mesh->getCellNumFaces(c);
    }
    else if (kind == AmanziMesh::Entity_kind::CELL) {
      ndofs = 1;
    }

    offset.push_back(ndofs * num);
  }
}


/* ******************************************************************
* Compatibility: returns old schema
****************************************************************** */
int Schema::OldSchema() const
{
  int i(0);

  // convert base
  if (base_ == AmanziMesh::Entity_kind::NODE) {
    i = OPERATOR_SCHEMA_BASE_NODE;
  } else if (base_ == AmanziMesh::Entity_kind::EDGE) {
    i = OPERATOR_SCHEMA_BASE_EDGE;
  } else if (base_ == AmanziMesh::Entity_kind::FACE) {
    i = OPERATOR_SCHEMA_BASE_FACE;
  } else if (base_ == AmanziMesh::Entity_kind::CELL) {
    i = OPERATOR_SCHEMA_BASE_CELL;
  }

  for (auto it = items_.begin(); it != items_.end(); ++it) {
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, std::ignore) = *it;

    if (kind == AmanziMesh::Entity_kind::NODE) {
      i += OPERATOR_SCHEMA_DOFS_NODE; 
    } else if (kind == AmanziMesh::Entity_kind::EDGE) {
      i += OPERATOR_SCHEMA_DOFS_EDGE; 
    } else if (kind == AmanziMesh::Entity_kind::FACE) {
      i += OPERATOR_SCHEMA_DOFS_FACE; 
    } else if (kind == AmanziMesh::Entity_kind::CELL) {
      i += OPERATOR_SCHEMA_DOFS_CELL;     
    } else if (kind == AmanziMesh::Entity_kind::BOUNDARY_FACE) {
      i += OPERATOR_SCHEMA_DOFS_BNDFACE;     
    }
  }
  return i;
}


/* ******************************************************************
* Returns standard name for geometric location of DOF.
****************************************************************** */
std::string Schema::KindToString(AmanziMesh::Entity_kind kind) const 
{
  if (kind == AmanziMesh::Entity_kind::NODE) {
    return "node";
  } else if (kind == AmanziMesh::Entity_kind::EDGE) {
    return "edge";
  } else if (kind == AmanziMesh::Entity_kind::FACE) {
    return "face";
  } else if (kind == AmanziMesh::Entity_kind::CELL) {
    return "cell";
  } else if (kind == AmanziMesh::Entity_kind::BOUNDARY_FACE) {
    return "boundary_face";
  }
  return "null";
}


/* ******************************************************************
* Returns standard mesh id for geometric location of DOF.
****************************************************************** */
AmanziMesh::Entity_kind Schema::StringToKind(const std::string& name) const 
{
  return AmanziMesh::createEntityKind(name);
}


/* ******************************************************************
* Returns standard mesh id for geometric location of DOF.
****************************************************************** */
WhetStone::DOF_Type Schema::StringToType(const std::string& name) const 
{
  if (name == "scalar") {
    return WhetStone::DOF_Type::SCALAR;
  } else if (name == "vector") {
    return WhetStone::DOF_Type::VECTOR;
  } else if (name == "point") {
    return WhetStone::DOF_Type::POINT;
  } else if (name == "normal component") {
    return WhetStone::DOF_Type::NORMAL_COMPONENT;
  } else if (name == "moment") {
    return WhetStone::DOF_Type::MOMENT;
  }
  return WhetStone::DOF_Type::SCALAR;
}


/* ******************************************************************
* Auxiliary routine creates new name.
****************************************************************** */
std::string Schema::CreateUniqueName() const
{
  std::string name(KindToString(base_)), c("_");
  for (auto it = items_.begin(); it != items_.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    name.append(c);
    name.append(KindToString(kind)); 
    name.append(std::to_string(num)); 
  }

  std::transform(name.begin(), name.end(), name.begin(), ::toupper);
  return name;
}

}  // namespace Operators
}  // namespace Amanzi


