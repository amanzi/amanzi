/*
  Process Kernels
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Function applied to a mesh component with at most one function 
  application per entity.
*/

#ifndef AMANZI_PK_BOUNDARY_FUNCTION_HH_
#define AMANZI_PK_BOUNDARY_FUNCTION_HH_

#include <map>
#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "MultiFunction.hh"
#include "UniqueMeshFunction.hh"

namespace Amanzi {

typedef std::pair<std::string, int> Action;

class PK_BoundaryFunction : public Functions::UniqueMeshFunction {
 public:
  PK_BoundaryFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      UniqueMeshFunction(mesh),
      finalized_(false), 
      global_size_(0) {};
  ~PK_BoundaryFunction() {};
  
  virtual void Define(const std::vector<std::string>& regions,
                      const Teuchos::RCP<const MultiFunction>& f, 
                      int method);
  virtual void Define(std::string& region,
                      const Teuchos::RCP<const MultiFunction>& f,
                      int method);

  virtual void Compute(double time);
  virtual void Finalize();

  // access / set
  const std::vector<Action>& actions() { return actions_; } 
  int global_size() { return global_size_; }

  // iterator methods
  typedef std::map<int,double>::const_iterator Iterator;
  Iterator begin() const { return value_.begin(); }
  Iterator end() const  { return value_.end(); }
  Iterator find(const int j) const { return value_.find(j); }
  std::map<int, double>::size_type size() { return value_.size(); }

 protected:
  std::map<int,double> value_;
  bool finalized_;

  std::vector<Action> actions_;
  int global_size_;  // number of data stored on all processors
};

}  // namespace Amanzi

#endif
