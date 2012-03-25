#ifndef AMANZI_MESH_FUNCTION_HH_
#define AMANZI_MESH_FUNCTION_HH_

#include <vector>
#include <set>
#include <map>
#include <utility>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "function.hh"

namespace Amanzi {

typedef std::set<AmanziMesh::Entity_ID> Domain;
typedef std::pair<Domain,Teuchos::RCP<const Function> > Spec;
typedef std::vector<Spec> SpecList;
typedef std::map<int,double>::const_iterator Iterator;

class MeshFunction {  
 public:
  MeshFunction() {};
  virtual ~MeshFunction() {};

  virtual void Define(const std::vector<std::string>& regions, const Teuchos::RCP<const Function>& f) = 0;
  void DefineFromString(const std::string region, const Teuchos::RCP<const Function>& f);
  
  void Compute(double);
  Iterator begin() const { return value_.begin(); }
  Iterator end() const  { return value_.end(); }
  Iterator find(const int j) const { return value_.find(j); }
  
  std::map<int,double>::size_type size() { return value_.size(); }
  
 public:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  SpecList spec_list_;
  std::map<int,double> value_;
};

} // namespace Amanzi

#endif  // AMANZI_MESH_FUNCTION_HH_
