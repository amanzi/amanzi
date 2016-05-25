/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatsky
*/

#ifndef AMANZI_TRANSPORT_BOUNDARY_FUNCTION_COUPLER_HH_
#define AMANZI_TRANSPORT_BOUNDARY_FUNCTION_COUPLER_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "MultiFunction.hh"
#include "TransportBoundaryFunction.hh"

namespace Amanzi {
namespace Transport {

class TransportBoundaryFunction_Coupler : public TransportBoundaryFunction {
 public:
  TransportBoundaryFunction_Coupler(const Teuchos::RCP<const State>& S, 
                                    const Teuchos::RCP<const AmanziMesh::Mesh> &mesh) :
    finalize_(false),
    TransportBoundaryFunction(S, mesh) {};
  ~TransportBoundaryFunction_Coupler() {};
  
  void Compute(double time);

  void Define(const std::vector<std::string> &regions, Teuchos::ParameterList& spec);


  void Define(std::string region, Teuchos::ParameterList& spec);


 private:
  void Finalize_();
  bool finalize_;
  Teuchos::ParameterList spec_list_;
};

}  // namespace Transport
}  // namespace Amanzi


#endif
