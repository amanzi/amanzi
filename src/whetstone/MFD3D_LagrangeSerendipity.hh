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

  Serendipity Lagrange-type element: degrees of freedom are nodal values
  and moments on edges, faces and inside cell. The number of later is
  reduced significantly for polygonal cells.
*/

#ifndef AMANZI_MFD3D_LAGRANGE_SERENDIPITY_HH_
#define AMANZI_MFD3D_LAGRANGE_SERENDIPITY_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MFD3D_LagrangeAnyOrder.hh"
#include "Polynomial.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_LagrangeSerendipity : public MFD3D_LagrangeAnyOrder {
 public:
  MFD3D_LagrangeSerendipity(const Teuchos::ParameterList& plist,
                            const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- stiffness matrix
  int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- l2 projectors
  virtual void L2Cell(int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments,
                      Polynomial& uc) override
  {
    ProjectorCell_(mesh_, c, ve, vf, ProjectorType::L2, moments, uc);
  }

  virtual void L2Face(int f,
                      const std::vector<Polynomial>& ve,
                      const Polynomial* moments,
                      Polynomial& uf) override
  {
    ProjectorFace_(f, ve, ProjectorType::L2, moments, uf);
  }

  // -- h1 projectors
  virtual void H1Cell(int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments,
                      Polynomial& uc) override
  {
    ProjectorCell_(mesh_, c, ve, vf, ProjectorType::H1, moments, uc);
  }

  virtual void H1Face(int f,
                      const std::vector<Polynomial>& ve,
                      const Polynomial* moments,
                      Polynomial& uf) override
  {
    ProjectorFace_(f, ve, ProjectorType::H1, moments, uf);
  }

  // other methods
  void L2Cell_LeastSquare(int c,
                          const std::vector<Polynomial>& vf,
                          const Polynomial* moments,
                          Polynomial& uc)
  {
    ProjectorCell_(mesh_, c, vf, vf, ProjectorType::LS, moments, uc);
  }

 private:
  void ProjectorCell_(const Teuchos::RCP<const AmanziMesh::Mesh>& mymesh,
                      int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const ProjectorType type,
                      const Polynomial* moments,
                      Polynomial& uc);

  void ProjectorFace_(int f,
                      const std::vector<Polynomial>& ve,
                      const ProjectorType type,
                      const Polynomial* moments,
                      Polynomial& uf);

  void CalculateDOFsOnBoundary_(const Teuchos::RCP<const AmanziMesh::Mesh>& mymesh,
                                int c,
                                const std::vector<Polynomial>& ve,
                                const std::vector<Polynomial>& vf,
                                DenseVector& vdof);

 private:
  static RegisteredFactory<MFD3D_LagrangeSerendipity> reg_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
