/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATOR_ANALYTIC_ELASTICITY_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_ELASTICITY_BASE_HH_

#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"
#include "Mesh.hh"

class AnalyticElasticityBase {
 public:
  AnalyticElasticityBase(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
    : mesh_(mesh)
  {
    nnodes_owned = mesh_->num_entities(
      Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::OWNED);
    ncells_owned = mesh_->num_entities(
      Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

    nnodes_wghost = mesh_->num_entities(Amanzi::AmanziMesh::NODE,
                                        Amanzi::AmanziMesh::Parallel_type::ALL);
    ncells_wghost = mesh_->num_entities(Amanzi::AmanziMesh::NODE,
                                        Amanzi::AmanziMesh::Parallel_type::ALL);
  };
  ~AnalyticElasticityBase(){};

  // analytic solution for elasticity-type problem
  // -- stiffness/elasticity tensor T
  virtual Amanzi::WhetStone::Tensor
  Tensor(const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  // -- analytic solution
  virtual Amanzi::AmanziGeometry::Point
  velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  virtual double
  pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  // -- source term
  virtual Amanzi::AmanziGeometry::Point
  source_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // error calculation
  // -- velocity
  void ComputeNodeError(Amanzi::CompositeVector& u, double t, double& unorm,
                        double& l2_err, double& inf_err)
  {
    unorm = 0.0;
    l2_err = 0.0;
    inf_err = 0.0;

    // calculate nodal volumes
    Amanzi::AmanziMesh::Entity_ID_List nodes;

    Teuchos::RCP<Amanzi::CompositeVectorSpace> cvs =
      Teuchos::rcp(new Amanzi::CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true)->AddComponent(
      "node", Amanzi::AmanziMesh::NODE, 1);

    Amanzi::CompositeVector vol(*cvs);
    Epetra_MultiVector& vol_node = *vol.ViewComponent("node", true);
    vol.putScalar(0.0);

    for (int c = 0; c != ncells_owned; ++c) {
      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      for (int i = 0; i < nnodes; i++) {
        vol_node[0][nodes[i]] += mesh_->cell_volume(c, false) / nnodes;
      }
    }
    vol.GatherGhostedToMaster("node");

    // calculate errors
    int d = mesh_->space_dimension();
    Amanzi::AmanziGeometry::Point p(d), ucalc(d);
    Epetra_MultiVector& u_node = *u.ViewComponent("node");

    for (int v = 0; v < nnodes_owned; ++v) {
      mesh_->node_get_coordinates(v, &p);

      const Amanzi::AmanziGeometry::Point& uexact = velocity_exact(p, t);
      for (int i = 0; i < d; ++i) ucalc[i] = u_node[i][v];

      double tmp = L22(ucalc - uexact);
      l2_err += tmp * vol_node[0][v];
      inf_err = std::max(inf_err, sqrt(tmp));
      unorm += L22(uexact) * vol_node[0][v];
      // std::cout << v << " uh=" << ucalc << " ex=" << uexact << std::endl;
    }
#ifdef HAVE_MPI
    double tmp = unorm;
    mesh_->get_comm()->SumAll(&tmp, &unorm, 1);
    tmp = l2_err;
    mesh_->get_comm()->SumAll(&tmp, &l2_err, 1);
    tmp = inf_err;
    mesh_->get_comm()->MaxAll(&tmp, &inf_err, 1);
#endif
    unorm = sqrt(unorm);
    l2_err = sqrt(l2_err);
  }

  // -- pressure
  void ComputeCellError(Amanzi::CompositeVector& p, double t, double& pnorm,
                        double& l2_err, double& inf_err)
  {
    pnorm = 0.0;
    l2_err = 0.0;
    inf_err = 0.0;

    Epetra_MultiVector& p_cell = *p.ViewComponent("cell");

    for (int c = 0; c < ncells_owned; ++c) {
      const Amanzi::AmanziGeometry::Point& xm = mesh_->cell_centroid(c);
      double area = mesh_->cell_volume(c, false);

      double pexact = pressure_exact(xm, t);
      double tmp = fabs(p_cell[0][c] - pexact);

      l2_err += tmp * tmp * area;
      inf_err = std::max(inf_err, tmp);
      pnorm += pexact * pexact * area;
      // std::cout << c << " ph=" << p_cell[0][c] << " ex=" << pexact <<
      // std::endl;
    }
#ifdef HAVE_MPI
    double tmp = pnorm;
    mesh_->get_comm()->SumAll(&tmp, &pnorm, 1);
    tmp = l2_err;
    mesh_->get_comm()->SumAll(&tmp, &l2_err, 1);
    tmp = inf_err;
    mesh_->get_comm()->MaxAll(&tmp, &inf_err, 1);
#endif
    pnorm = sqrt(pnorm);
    l2_err = sqrt(l2_err);
  }

 protected:
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_;

  int nnodes_owned, ncells_owned, nnodes_wghost, ncells_wghost;
};

#endif
