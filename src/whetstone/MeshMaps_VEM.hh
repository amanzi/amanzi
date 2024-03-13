/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Maps between mesh objects located on different meshes, e.g. two
  states of a deformable mesh: virtual element implementation.
*/

#ifndef AMANZI_WHETSTONE_MESH_MAPS_VEM_HH_
#define AMANZI_WHETSTONE_MESH_MAPS_VEM_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "MeshMapsBase.hh"
#include "Polynomial.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class MeshMaps_VEM : public MeshMapsBase {
 public:
  MeshMaps_VEM(Teuchos::RCP<const AmanziMesh::Mesh> mesh, const Teuchos::ParameterList& plist)
    : MeshMapsBase(mesh)
  {
    ParseInputParameters_(plist);
  }
  MeshMaps_VEM(Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
               Teuchos::RCP<const AmanziMesh::Mesh> mesh1,
               const Teuchos::ParameterList& plist)
    : MeshMapsBase(mesh0, mesh1)
  {
    ParseInputParameters_(plist);
  }
  ~MeshMaps_VEM(){};

  // remap pseudo velocity
  virtual void VelocityFace(int f, VectorPolynomial& vf) const override;
  virtual void VelocityCell(int c,
                            const std::vector<VectorPolynomial>& ve,
                            const std::vector<VectorPolynomial>& vf,
                            VectorPolynomial& vc) const override;

 private:
  void LeastSquareProjector_Cell_(int order,
                                  int c,
                                  const std::vector<VectorPolynomial>& vf,
                                  VectorPolynomial& vc) const;

  // io
  void ParseInputParameters_(const Teuchos::ParameterList& plist);

 private:
  std::string method_, projector_;
  int order_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
