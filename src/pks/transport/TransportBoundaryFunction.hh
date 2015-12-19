/*
  Transport PK 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson
           Ethan Coon
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_TRANSPORT_BOUNDARY_FUNCTION_HH_
#define AMANZI_TRANSPORT_BOUNDARY_FUNCTION_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "UniqueMeshFunction.hh"

namespace Amanzi {
namespace Transport {

class TransportBoundaryFunction : public Functions::UniqueMeshFunction {
 public:
  TransportBoundaryFunction(const Teuchos::RCP<const AmanziMesh::Mesh> &mesh) :
      UniqueMeshFunction(mesh) {};
  virtual ~TransportBoundaryFunction() {};
  
  virtual void Compute(double time) = 0;

  // access
  std::vector<std::string>& tcc_names() { return tcc_names_; }
  std::vector<int>& tcc_index() { return tcc_index_; }

  std::vector<int>& faces() { return faces_; }
  std::vector<std::vector<double> >& values() { return values_; }

 protected:
  std::vector<std::string> tcc_names_;  // list of component names
  std::vector<int> tcc_index_;  // index of component in the global list

  std::vector<int> faces_;  // list of boundary faces 
  std::vector<std::vector<double> > values_;  // component values on boundary faces
};

}  // namespace Transport
}  // namespace Amanzi


#endif
