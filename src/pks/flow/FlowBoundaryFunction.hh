/*
  Flow PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1)
           Ethan Coon (version 2)
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)

  Function applied to a mesh component with at most one function 
  application per entity.
*/

#ifndef AMANZI_FLOW_BOUNDARY_FUNCTION_HH_
#define AMANZI_FLOW_BOUNDARY_FUNCTION_HH_

#include <map>
#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "PK_BoundaryFunction.hh"

namespace Amanzi {
namespace Flow {

class FlowBoundaryFunction : public PK_BoundaryFunction {
 public:
  FlowBoundaryFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
                       PK_BoundaryFunction(mesh) {};
  
  void ComputeShift(double T, double* shift);

  // access / set
  double reference_pressure() { return reference_pressure_; }
  void set_reference_pressure(double p0) { reference_pressure_ = p0; }

 private:
  double reference_pressure_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif
