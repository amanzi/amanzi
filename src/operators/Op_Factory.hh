//! Helper factory for storing Ops in State

/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OP_FACTORY_HH_
#define AMANZI_OP_FACTORY_HH_

#include "Op.hh"
#include "Op_Cell_Cell.hh"
#include "Op_Face_Cell.hh"
#include "Op_Cell_FaceCell.hh"

namespace Amanzi {
namespace Operators {

class Op_Factory {
 public:
  Op_Factory() {};

  void set_mesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) { mesh_ = mesh; }
  Teuchos::RCP<const AmanziMesh::Mesh> mesh() const { return mesh_; }
  
  void set_name(const std::string& name) { name_ = name; }
  void set_schema(const Schema& schema_row, const Schema& schema_col) {
    schema_row_ = schema_row;
    schema_col_ = schema_col;
  }
  void set_schema(const Schema& schema) {
    schema_row_ = schema;
    schema_col_ = schema;
    set_name(schema.CreateUniqueName());
  }
  const Schema& schema_col() const { return schema_col_; }
  const Schema& schema_row() const { return schema_row_; }
  

  Teuchos::RCP<Op> Create() const {
    AMANZI_ASSERT(mesh_.get());
    AMANZI_ASSERT(schema_row_.size() > 0);
    AMANZI_ASSERT(schema_row_.OldSchema() == schema_col_.OldSchema());

    if (schema_row_.OldSchema() == (OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL)) {
      return Teuchos::rcp(new Op_Cell_Cell(name_, mesh_));
    } else if (schema_row_.OldSchema() == 
               (OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_FACE | OPERATOR_SCHEMA_DOFS_CELL)) {
      return Teuchos::rcp(new Op_Cell_FaceCell(name_, mesh_));
    } else if (schema_row_.OldSchema() ==
               (OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL)) {
      return Teuchos::rcp(new Op_Face_Cell(name_, mesh_));
    } else {
      Errors::Message message("Unimeplemented Op Schema requested in Op_Factory");
      throw(message);
    }
    return Teuchos::null;
  }
  
 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::string name_;
  Schema schema_row_, schema_col_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

