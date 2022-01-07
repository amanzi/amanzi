/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Mimetic discretization of elliptic operator using edge-based
  degrees of freedom shows flexibility of the discretization framework.
 
  Usage:
 
  1. Include base class for the mimetic methods and the factory of
     discretiation methods:

  #include "MFD3D.hh"
  #include "BilinearFormFactory.hh"

  2. Add variable for the static registry. In this example MyMethod 
  is Diffusion_Edge:

  static RegisteredFactory<MFD3D_My_Method> factory_;

  3. Explicitly instantiate the static registry (in .cc file) with
     the unique name "my unique method". The factory takes this name
     from the parameter list.

  RegisteredFactory<MFD3D_MyMethod> MFD3D_MyMethod::factory_("my unique method");
*/

#ifndef AMANZI_WHETSTONE_MFD3D_DIFFUSION_EDGE_HH_
#define AMANZI_WHETSTONE_MFD3D_DIFFUSION_EDGE_HH_

#include <tuple>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_Diffusion_Edge : public MFD3D {
 public:
  MFD3D_Diffusion_Edge(const Teuchos::ParameterList& plist,
                       const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : MFD3D(mesh) {};

  // main methods 
  // -- symmetric schema
  virtual std::vector<SchemaItem> schema() const override {
    return std::vector<SchemaItem>(1, std::make_tuple(AmanziMesh::Entity_kind::EDGE, DOF_Type::SCALAR, 1));
  }

  // -- stiffness matrix
  int H1consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac);
  virtual int StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A) override;

 private:
  static RegisteredFactory<MFD3D_Diffusion_Edge> factory_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

