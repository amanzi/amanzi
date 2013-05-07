/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $AMANZI_DIR/COPYRIGHT
Author (v1): Neil Carlson
       (v2): Ethan Coon

Function applied to a mesh component with at most one function application per
entity.

------------------------------------------------------------------------- */


#ifndef AMANZI_TRANSPORT_BOUNDARY_FUNCTION_HH_
#define AMANZI_TRANSPORT_BOUNDARY_FUNCTION_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "MultiFunction.hh"
#include "unique_mesh_function.hh"

namespace Amanzi {
namespace Functions {

typedef std::pair<std::string, int> Action;

class TransportBoundaryFunction : public UniqueMeshFunction {
  
public:
  TransportBoundaryFunction(const Teuchos::RCP<const AmanziMesh::Mesh> &mesh) :
      UniqueMeshFunction(mesh),
      finalized_(false) {}
  
  void Define(const std::vector<std::string> &regions,
              const Teuchos::RCP<const MultiFunction> &f);

  void Define(std::string region,
              const Teuchos::RCP<const MultiFunction> &f);

  void Compute(double time);

  void ComputeShift(double T, double* shift);

  void Finalize();

  // iterator methods
  typedef std::map<int,double>::const_iterator Iterator;
  Iterator begin() const { return value_.begin(); }
  Iterator end() const  { return value_.end(); }
  Iterator find(const int j) const { return value_.find(j); }
  std::map<int,double>::size_type size() { return value_.size(); }

protected:
  std::map<int,double> value_;
  bool finalized_;
};

} // namespace
} // namespace


#endif
