/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef __AMANZI_DOMAIN_FUNCTION_HH__
#define __AMANZI_DOMAIN_FUNCTION_HH__

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "function.hh"
#include "unique_mesh_function.hh"

namespace Amanzi {
namespace Functions {

class DomainFunction : public UniqueMeshFunction {
 public:
  explicit DomainFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : 
      UniqueMeshFunction(mesh),
      finalized_(false) {}

  void Define(const std::vector<std::string>& regions,
              const Teuchos::RCP<const VectorFunction>& f);
  void Define(std::string region,
              const Teuchos::RCP<const VectorFunction> &f);

  void Compute(double time);

  void Finalize() {};
 
  // iterator methods
  typedef std::map<int,double>::const_iterator Iterator;
  Iterator begin() const { return value_.begin(); }
  Iterator end() const  { return value_.end(); }
  Iterator find(const int j) const { return value_.find(j); }
  std::map<int,double>::size_type size() { return value_.size(); }
 

 private:
  std::vector<int> actions_;


 protected:
  std::map<int,double> value_;
  bool finalized_;

};


}  // namespace Functions
}  // namespace Amanzi

#endif  // AMANZI_DOMAIN_FUNCTION_HH_
