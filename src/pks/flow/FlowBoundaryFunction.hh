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

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "MultiFunction.hh"
#include "unique_mesh_function.hh"

namespace Amanzi {
namespace Flow {

class FlowBoundaryFunction : public Functions::UniqueMeshFunction {
 public:
  FlowBoundaryFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      UniqueMeshFunction(mesh),
      finalized_(false), 
      global_size_(0) {};
  
  void Define(const std::vector<std::string>& regions,
              const Teuchos::RCP<const MultiFunction>& f, 
              int method);
  void Define(std::string& region,
              const Teuchos::RCP<const MultiFunction>& f,
              int method);

  void Compute(double time);
  void ComputeShift(double T, double* shift);

  void Finalize();

  // access / set
  const std::vector<CommonDefs::Action>& actions() { return actions_; } 
  double reference_pressure() { return reference_pressure_; }
  int global_size() { return global_size_; }
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
  std::vector<CommonDefs::Action> actions_;
  double reference_pressure_;
  int global_size_;  // number of data stored on all processors
};

}  // namespace Flow
}  // namespace Amanzi

#endif
