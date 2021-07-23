//! A factory for creating Operator objects (the global version)

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

/*!

*/

#include "Operator_Factory.hh"
#include "Operator_Cell.hh"
#include "Operator_FaceCell.hh"

namespace Amanzi {
namespace Operators {

Teuchos::RCP<Operator>
Operator_Factory::Create() {
  if (!plist_.get())
    plist_ = Teuchos::rcp(new Teuchos::ParameterList("operator"));

  // deduce the type
  // -- first choice: type is provided in the plist
  std::string operator_type;
  if (plist_->isParameter("operator type")) {
    operator_type = plist_->get<std::string>("operator type");

    if (operator_type == "Operator_Cell") {
      // build the CVS from the global schema
      CompositeVectorSpace cvs;
      cvs.SetMesh(mesh_)->SetGhosted(true);
      cvs.AddComponent("cell", AmanziMesh::CELL, 1);
      return Teuchos::rcp(new Operator_Cell(cvs.CreateSpace(), *plist_, OPERATOR_SCHEMA_DOFS_CELL));

    } else {
      Errors::Message msg;
      msg << "Operator_Factory: unsupported operator type.";
      Exceptions::amanzi_throw(msg);
    }

  } else if (cvs_row_->size() != 0) {
    if (cvs_row_->HasComponent("cell")) {
      if (cvs_row_->HasComponent("face")) {
        return Teuchos::rcp(new Operator_FaceCell(cvs_row_, *plist_));
      } else {
        return Teuchos::rcp(new Operator_Cell(cvs_row_, *plist_, OPERATOR_SCHEMA_DOFS_CELL));
      }
    } else {
      Errors::Message msg;
      msg << "Operator_Factory: unsupported CompositeVector's component.";
      Exceptions::amanzi_throw(msg);
    }
  }
  return Teuchos::null;
}

}  // namespace Operators
}  // namespace Amanzi
