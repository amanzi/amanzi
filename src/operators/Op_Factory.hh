//! Helper factory for storing Ops in State

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OP_FACTORY_HH_
#define AMANZI_OP_FACTORY_HH_

namespace Amanzi {
namespace Operators {

template<class Op_t>
class Op_Factory {
 public:
  Op_Factory() {}

  void set_mesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) { mesh_ = mesh; }
  void set_name(const std::string& name) { name_ = name; }
  Teuchos::RCP<Op_t> Create() const {
    return Teuchos::rcp(new Op_t(name_, mesh_));
  }

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::string name_;

};

}  // namespace Operators
}  // namespace Amanzi


#endif

