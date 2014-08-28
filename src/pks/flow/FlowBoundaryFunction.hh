/*
  This is the flow component of the Amanzi code.
 
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

#include "Mesh.hh"
#include "MultiFunction.hh"
#include "unique_mesh_function.hh"

namespace Amanzi {
namespace Functions {

const int BOUNDARY_FUNCTION_ACTION_NONE = 0;
const int BOUNDARY_FUNCTION_ACTION_HEAD_RELATIVE = 1;

typedef std::pair<std::string, int> Action;

class FlowBoundaryFunction : public UniqueMeshFunction {
 public:
  FlowBoundaryFunction(const Teuchos::RCP<const AmanziMesh::Mesh> &mesh) :
      UniqueMeshFunction(mesh),
      finalized_(false) {}
  
  void Define(const std::vector<std::string> &regions,
              const Teuchos::RCP<const MultiFunction> &f, 
              int method);
  void Define(std::string region,
              const Teuchos::RCP<const MultiFunction> &f,
              int method);

  void Compute(double time);
  void ComputeShift(double T, double* shift);

  void Finalize();

  // access / set
  const std::vector<Action>& actions() { return actions_; } 
  double reference_pressure() { return reference_pressure_; }
  void set_reference_pressure(double p0) { reference_pressure_ = p0; }

  // iterator methods
  typedef std::map<int,double>::const_iterator Iterator;
  Iterator begin() const { return value_.begin(); }
  Iterator end() const  { return value_.end(); }
  Iterator find(const int j) const { return value_.find(j); }
  std::map<int,double>::size_type size() { return value_.size(); }

 protected:
  std::map<int,double> value_;
  bool finalized_;

 private:
  std::vector<Action> actions_;
  double reference_pressure_;
};

}  // namespace Functions
}  // namespace Amanzi

#endif
