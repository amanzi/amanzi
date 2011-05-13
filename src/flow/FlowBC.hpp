#ifndef __FLOWBC_H__
#define __FLOWBC_H__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"

class FlowBC {

public:

  enum bc_types {
    PRESSURE_CONSTANT = 1,
    NO_FLOW,
    DARCY_CONSTANT,
    STATIC_HEAD
  };

  struct bc_spec {
    bc_types Type;
    unsigned int SetID;
    std::vector<unsigned int> Faces;
    std::vector<double> Aux;
    double Value;
  };

public:
  FlowBC(Teuchos::ParameterList &params, const Teuchos::RCP<Amanzi::AmanziMesh::Mesh> &mesh);
  ~FlowBC() {}

  const int NumBC () const { return bc_.size(); }

  bc_spec& operator [] (int index) { return bc_[index]; }
  // can't do this because we need to write/read from the Aux member
  //const bc_spec& operator [] (int index) const { return bc_[index]; }

private:

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_;

  std::vector<bc_spec> bc_;

};

#endif
