/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov
*/

#ifndef AMANZI_TRANSPORT_BOUNDARY_FUNCTION_TRACER_HH_
#define AMANZI_TRANSPORT_BOUNDARY_FUNCTION_TRACER_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "MultiFunction.hh"
#include "TransportBoundaryFunction.hh"

namespace Amanzi {
namespace Transport {

class TransportBoundaryFunction_Tracer : public TransportBoundaryFunction {
 public:
  TransportBoundaryFunction_Tracer(const Teuchos::RCP<const State>& S,
                                   const Teuchos::RCP<const AmanziMesh::Mesh> &mesh) :
    TransportBoundaryFunction(S,mesh) {};
  ~TransportBoundaryFunction_Tracer() {};
  
  void Compute(double time);

  void Define(const std::vector<std::string> &regions,
              const Teuchos::RCP<const MultiFunction> &f);

  void Define(std::string region,
              const Teuchos::RCP<const MultiFunction> &f);

 private:
  void Finalize_();
};

}  // namespace Transport
}  // namespace Amanzi


#endif
