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

#ifndef AMANZI_OPERATOR_FACTORY_HH_
#define AMANZI_OPERATOR_FACTORY_HH_

#include "Schema.hh"
#include "Operator.hh"

namespace Amanzi {
namespace Operators {

class Operator_Factory {
 public:
  Operator_Factory() {};

  void set_mesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) { mesh_ = mesh; }
  void set_cvs(const CompositeVectorSpace& cvs) {
    cvs_row_ = cvs;
    cvs_col_ = cvs;
  }
  void set_cvs(const CompositeVectorSpace& cvs_row, const CompositeVectorSpace& cvs_col) {
    cvs_row_ = cvs_row;
    cvs_col_ = cvs_col;
  }
  void set_schema(const Schema& schema) {
    schema_row_ = schema;
    schema_col_ = schema;
  }
  void set_schema(const Schema& schema_row, const Schema& schema_col) {
    schema_row_ = schema_row;
    schema_col_ = schema_col;
  }
  void set_plist(const Teuchos::RCP<Teuchos::ParameterList>& plist) {
    plist_ = plist;
  }
  
  Teuchos::RCP<Operator> Create();

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  CompositeVectorSpace cvs_row_, cvs_col_;
  Schema schema_row_, schema_col_;

  Teuchos::RCP<Teuchos::ParameterList> plist_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif
