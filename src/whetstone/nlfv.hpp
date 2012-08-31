/*
This is the mimetic discretization component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Release name: aka-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
*/

#ifndef __NLFV_HPP__
#define __NLFV_HPP__

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

