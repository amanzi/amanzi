/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for testing elasticity-type problems.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_ELASTICITY_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_ELASTICITY_BASE_HH_

#include "CompositeVector.hh"
#include "Mesh.hh"

class AnalyticElasticityBase {
 public:
  AnalyticElasticityBase(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : mesh_(mesh) {};
  ~AnalyticElasticityBase() {};

  // analytic solution for elasticity-type problem
  // -- stiffness/elasticity tensor T
  virtual Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  // -- analytic solution 
  virtual Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  virtual double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  // -- source term
  virtual Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // error calculation
  void ComputeNodeError(Amanzi::CompositeVector& u, double t, double& unorm, double& l2_err, double& inf_err) {
    unorm = 0.0;
    l2_err = 0.0;
    inf_err = 0.0;

    int d = mesh_->space_dimension(); 
    Amanzi::AmanziGeometry::Point p(d), ucalc(d);
    Epetra_MultiVector& u_node = *u.ViewComponent("node");

    int nnodes = mesh_->num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::OWNED);
    for (int v = 0; v < nnodes; ++v) {
      mesh_->node_get_coordinates(v, &p);

      const Amanzi::AmanziGeometry::Point& uexact = velocity_exact(p, t);
      for (int i = 0; i < d; ++i) ucalc[i] = u_node[i][v];

      l2_err += norm(ucalc - uexact);
      std::cout << v << " uh=" << ucalc << " ex=" << uexact << std::endl;
      unorm += norm(uexact);
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

 protected:
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_;
};

#endif

