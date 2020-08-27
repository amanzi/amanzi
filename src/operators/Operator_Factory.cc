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

#include <tuple>

#include "UniqueHelpers.hh"

#include "Operator_Cell.hh"
#include "Operator_FaceCell.hh"
#include "Operator_FaceCellSff.hh"
#include "Operator_Factory.hh"
#include "Operator_Edge.hh"
#include "Operator_Node.hh"
#include "Operator_Schema.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Original version
****************************************************************** */
Teuchos::RCP<Operator>
Operator_Factory::Create()
{
  if (!plist_.get())
    plist_ = Teuchos::rcp(new Teuchos::ParameterList("operator"));

  // deduce the type
  // -- first choice: type is provided in the plist
  std::string pde_type;
  if (plist_->isParameter("operator type")) {
    pde_type = plist_->get<std::string>("operator type");

    if (pde_type == "Operator_Cell") {
      // build the CVS from the global schema
      Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
      cvs->SetMesh(mesh_)->SetGhosted(true);
      cvs->AddComponent("cell", AmanziMesh::CELL, 1);

      return Teuchos::rcp(new Operator_Cell(cvs, *plist_, OPERATOR_SCHEMA_DOFS_CELL));

    } else {
      Errors::Message msg;
      msg << "Operator_Factory: unsupported operator type.";
      Exceptions::amanzi_throw(msg);
    }

  } else if (cvs_row_.size() != 0) {
    if (cvs_row_.HasComponent("cell")) {
      if (cvs_row_.HasComponent("face")) {
        auto cvs_row = Teuchos::rcp(new CompositeVectorSpace(cvs_row_));
        return Teuchos::rcp(new Operator_FaceCell(cvs_row, *plist_));
      } else {
        auto cvs_row = Teuchos::rcp(new CompositeVectorSpace(cvs_row_));
        return Teuchos::rcp(new Operator_Cell(cvs_row, *plist_, OPERATOR_SCHEMA_DOFS_CELL));
      }
    } else {
      Errors::Message msg;
      msg << "Operator_Factory: unsupported CompositeVector's component.";
      Exceptions::amanzi_throw(msg);
    }
  }
  return Teuchos::null;
}


/* ******************************************************************
* Schema-based version
****************************************************************** */
Teuchos::RCP<Operator>
Operator_Factory::CreateFromSchema()
{
  AMANZI_ASSERT(schema_row_.size() > 0 && schema_col_.size() > 0);
  AMANZI_ASSERT(schema_row_.base() == schema_col_.base());
  AMANZI_ASSERT(schema_row_.base() == AmanziMesh::CELL);

  if (!plist_.get())
    plist_ = Teuchos::rcp(new Teuchos::ParameterList("operator"));

  auto base = schema_row_.base();
  int size = schema_row_.size();
  int num1 = std::get<2>(schema_row_[0]);  // number of dofs

  auto cvs1 = Teuchos::rcp(new CompositeVectorSpace(cvsFromSchema(schema_row_, mesh_, true)));

  // named operator is the best choise for a square operator
  if (schema_row_ == schema_col_ && size == 1 && num1 == 1) {
    int old_schema = schema_row_.OldSchema();
    auto entity = std::get<0>(schema_row_[0]);  

    if (entity == AmanziMesh::CELL) {
      return Teuchos::rcp(new Operator_Cell(cvs1, *plist_, old_schema));
    } else if (entity == AmanziMesh::FACE) {
      cvs1->AddComponent("cell", AmanziMesh::CELL, 1);
      return Teuchos::rcp(new Operator_FaceCellSff(cvs1, *plist_));
    } else if (entity == AmanziMesh::EDGE) {
      return Teuchos::rcp(new Operator_Edge(cvs1, *plist_));
    } else if (entity == AmanziMesh::NODE) {
      return Teuchos::rcp(new Operator_Node(cvs1, *plist_));
    }
  } 

  // named operator is the best choise for a square operator
  if (schema_row_ == schema_col_ && size == 2) {
    int num2 = std::get<2>(schema_row_[1]);
  
    auto ent1 = std::get<0>(schema_row_[0]);  
    auto ent2 = std::get<0>(schema_row_[1]);
    if (num1 == num2 && num1 == 1) {
      if ((ent1 == AmanziMesh::CELL && ent2 == AmanziMesh::FACE) ||
          (ent2 == AmanziMesh::CELL && ent1 == AmanziMesh::FACE))
      return Teuchos::rcp(new Operator_FaceCell(cvs1, *plist_));
    }
  }

  // abstract operator is the best choice for multiple dofs, more
  // then two items in schema 
  num1 = 0;
  for (int i = 0; i < schema_row_.size(); ++i)
    num1 = std::max(num1, std::get<2>(schema_row_[i]));
  for (int i = 0; i < schema_col_.size(); ++i)
    num1 = std::max(num1, std::get<2>(schema_col_[i]));

  if (num1 > 1) {
    auto cvs2 = Teuchos::rcp(new CompositeVectorSpace(cvsFromSchema(schema_col_, mesh_, true)));
    return Teuchos::rcp(new Operator_Schema(cvs1, cvs2, *plist_, schema_row_, schema_col_));
  }
  
  Errors::Message msg;
  msg << "Operator factory failed: unsupported combination of schema parameters";
  Exceptions::amanzi_throw(msg);
  return Teuchos::null;
}

}  // namespace Operators
}  // namespace Amanzi
