/*
This is the mimetic discretization component of the Amanzi code. 
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage:  
*/

#ifndef __WHETSTONE_hpp__
#define __WHETSTONE_hpp__

/*
This is the discretization package, release alpha.
*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "tensor.hpp"
#include "Point.hh"
#include "Mesh.hh"


namespace Amanzi {
namespace WhetStone {

class NLFV { 
 public:
  NLFV(Teuchos::RCP<AmanziMesh::Mesh> mesh) { mesh_ = mesh; };
  ~NLFV();
 private:
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

