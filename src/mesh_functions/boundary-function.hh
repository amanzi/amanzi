#ifndef AMANZI_BOUNDARY_FUNCTION_HH_
#define AMANZI_BOUNDARY_FUNCTION_HH_

#include <vector>
#include <set>
#include <map>
#include <utility>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "function.hh"

namespace Amanzi {

class BoundaryFunction {
  
  typedef std::set<AmanziMesh::Entity_ID> Domain;
  typedef std::pair<Domain,Teuchos::RCP<const Function> > Spec;
  typedef std::vector<Spec> SpecList;
  
 public:
  BoundaryFunction(const Teuchos::RCP<const AmanziMesh::Mesh> &mesh) : mesh_(mesh) {}
  
  void Define(const std::vector<int> &regions, const Teuchos::RCP<const Function> &f);
  void Define(int region, const Teuchos::RCP<const Function> &f);
  
  typedef const std::map<int,double>& Value;
  typedef std::map<int,double>::const_iterator ValueIterator;
  Value operator() (double t) { Compute(t); return value_; }
  
  std::map<int,double>::size_type size() { return value_.size(); }
  
  // another possible scheme
  void Compute(double);
  typedef std::map<int,double>::const_iterator Iterator;
  Iterator begin() { return value_.begin(); }
  Iterator end()   { return value_.end(); }
  
 private:
  const Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  SpecList spec_list_;
  std::map<int,double> value_;
};

} // namespace Amanzi

#endif // AMANZI_BOUNDARY_FUNCTION_HH_
