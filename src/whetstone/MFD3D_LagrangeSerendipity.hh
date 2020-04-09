/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

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
#include "MFD3D_Lagrange.hh"
#include "Polynomial.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_LagrangeSerendipity : public MFD3D_Lagrange {
 public:
  MFD3D_LagrangeSerendipity(const Teuchos::ParameterList& plist,
                            const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~MFD3D_LagrangeSerendipity(){};

  // required methods
  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N,
                            DenseMatrix& Ac) override;
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- cell projectors
  virtual void L2Cell(int c, const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override
  {
    ProjectorCell_(c, vf, ProjectorType::L2, moments, uc);
  }

  virtual void H1Cell(int c, const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override
  {
    ProjectorCell_(c, vf, ProjectorType::H1, moments, uc);
  }

  // -- face projectors
  virtual void H1Face(int f, const std::vector<Polynomial>& ve,
                      const Polynomial* moments, Polynomial& uf) override
  {
    ProjectorFace_(f, ve, ProjectorType::H1, moments, uf);
  }

  // other methods
  void L2Cell_LeastSquare(int c, const std::vector<Polynomial>& vf,
                          const Polynomial* moments, Polynomial& uc)
  {
    ProjectorCell_(c, vf, ProjectorType::LS, moments, uc);
  }

 private:
  void ProjectorCell_(int c, const std::vector<Polynomial>& vf,
                      const ProjectorType type, const Polynomial* moments,
                      Polynomial& uc);
  void ProjectorFace_(int f, const std::vector<Polynomial>& ve,
                      const ProjectorType type, const Polynomial* moments,
                      Polynomial& uf)
  {
    AMANZI_ASSERT(false);
  }

  void CalculateDOFsOnBoundary_(int c, const std::vector<Polynomial>& vf,
                                DenseVector& vdof);

 private:
  static RegisteredFactory<MFD3D_LagrangeSerendipity> factory_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
