/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  This is the energy component of the Amanzi code.

*/

#ifndef AMANZI_ENERGY_ANALYTIC_BASE_HH_
#define AMANZI_ENERGY_ANALYTIC_BASE_HH_

#include "Mesh.hh"
#include "Tensor.hh"

class AnalyticBase {
 public:
  AnalyticBase(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : mesh_(mesh){};
  ~AnalyticBase(){};

  // problem coefficients: conductivity and fluid velocity
  virtual Amanzi::WhetStone::Tensor
  Conductivity(int c, const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  virtual Amanzi::AmanziGeometry::Point
  FluidVelocity(int c, const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // primary variable is temperature: provide its derivatives
  virtual double temperature_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  virtual Amanzi::AmanziGeometry::Point
  flux_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // error calculation: L2 and maximum norms
  void ComputeCellError(const Epetra_MultiVector& temp,
                        double t,
                        double& l2_norm,
                        double& l2_err,
                        double& inf_err)
  {
    l2_norm = 0.0;
    l2_err = 0.0;
    inf_err = 0.0;

    int ncells =
      mesh_->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
    for (int c = 0; c < ncells; c++) {
      const Amanzi::AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
      double tmp = temperature_exact(xc, t);
      double volume = mesh_->getCellVolume(c);

      double err = std::fabs(tmp - temp[0][c]);
      l2_err += std::pow(err, 2.0) * volume;
      inf_err = std::max(inf_err, err);
      l2_norm += std::pow(tmp, 2.0) * volume;
      // std::cout << c << " " << tmp << " " << temp[0][c] << " err=" << err << std::endl;
    }
#ifdef HAVE_MPI
    double tmp = l2_norm;
    mesh_->getComm()->SumAll(&tmp, &l2_norm, 1);
    tmp = l2_err;
    mesh_->getComm()->SumAll(&tmp, &l2_err, 1);
    tmp = inf_err;
    mesh_->getComm()->MaxAll(&tmp, &inf_err, 1);
#endif
    l2_norm = sqrt(l2_norm);
    l2_err = sqrt(l2_err);
  }

 protected:
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_;
};

#endif
