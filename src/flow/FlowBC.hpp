#ifndef __FLOWBC_H__
#define __FLOWBC_H__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh_maps_base.hh"

class FlowBC {

public:

  enum bc_types {
    PRESSURE_CONSTANT = 1,
    NO_FLOW,
    DARCY_CONSTANT,
    STATIC_HEAD,
    TIME_DEPENDENT_PRESSURE_CONSTANT
  };

  struct bc_spec {
    bc_types Type;
    unsigned int SetID;
    std::vector<unsigned int> Faces;
    std::vector<double> Aux;
    double Value;
    double InitialValue;
    double FinalTime;
    double InitialTime;
  };

public:
  FlowBC(Teuchos::ParameterList &params, const Teuchos::RCP<Mesh_maps_base> &mesh);
  ~FlowBC() {}

  const int NumBC () const { return bc_.size(); }

  bc_spec& operator [] (int index) { return bc_[index]; }
  // can't do this because we need to write/read from the Aux member
  //const bc_spec& operator [] (int index) const { return bc_[index]; }

private:

  Teuchos::RCP<Mesh_maps_base> mesh_;

  std::vector<bc_spec> bc_;

};

#endif
