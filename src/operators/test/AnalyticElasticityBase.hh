/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for testing elasticity and Stokes-type problems.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_ELASTICITY_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_ELASTICITY_BASE_HH_

#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"
#include "Mesh.hh"

class AnalyticElasticityBase {
 public:
  AnalyticElasticityBase(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : mesh_(mesh) {
    nnodes_owned = mesh_->num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::OWNED);
    ncells_owned = mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

    nnodes_wghost = mesh_->num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::ALL);
    ncells_wghost = mesh_->num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::ALL);
  };
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
  // -- velocity or displacement
  void VectorNodeError(Amanzi::CompositeVector& u, double t, double& unorm, double& l2_err, double& inf_err);
  void VectorCellError(Amanzi::CompositeVector& u, double t, double& unorm, double& l2_err, double& inf_err);
  // -- pressure
  void ScalarCellError(Amanzi::CompositeVector& p, double t, double& pnorm, double& l2_err, double& inf_err);

  // communications
  void GlobalOp(std::string op, double* val, int n);

 protected:
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_;

  int nnodes_owned, ncells_owned, nnodes_wghost, ncells_wghost;
};


/* ******************************************************************
* Error in displacement.
****************************************************************** */
inline
void AnalyticElasticityBase::VectorNodeError(
    Amanzi::CompositeVector& u, double t, double& unorm,
    double& l2_err, double& inf_err)
{
  unorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  // calculate nodal volumes
  Amanzi::AmanziMesh::Entity_ID_List nodes;

  Teuchos::RCP<Amanzi::CompositeVectorSpace> cvs = Teuchos::rcp(new Amanzi::CompositeVectorSpace());
  cvs->SetMesh(mesh_)->SetGhosted(true)->AddComponent("node", Amanzi::AmanziMesh::NODE, 1);

  Amanzi::CompositeVector vol(*cvs);
  Epetra_MultiVector& vol_node = *vol.ViewComponent("node", true);
  vol.PutScalar(0.0);

  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    for (int i = 0; i < nnodes; i++) {
      vol_node[0][nodes[i]] += mesh_->cell_volume(c) / nnodes; 
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
    // std::cout << v << " uh=" << ucalc << " ex=" << uexact << " xv=" << p << std::endl;
  }
#ifdef HAVE_MPI
  GlobalOp("sum", &unorm, 1);
  GlobalOp("sum", &l2_err, 1);
  GlobalOp("max", &inf_err, 1);
#endif
  unorm = sqrt(unorm);
  l2_err = sqrt(l2_err);
}


/* ******************************************************************
* Error in velocity or displacement
****************************************************************** */
inline
void AnalyticElasticityBase::VectorCellError(
    Amanzi::CompositeVector& u, double t, double& unorm,
    double& l2_err, double& inf_err)
{
  unorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  int d = mesh_->space_dimension(); 
  Amanzi::AmanziGeometry::Point ucalc(d);
  Epetra_MultiVector& u_cell = *u.ViewComponent("cell");

  for (int c = 0; c < ncells_owned; ++c) {
    const Amanzi::AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    double vol = mesh_->cell_volume(c);

    const Amanzi::AmanziGeometry::Point& uexact = velocity_exact(xc, t);
    for (int i = 0; i < d; ++i) ucalc[i] = u_cell[i][c];

    double tmp = L22(ucalc - uexact);
    l2_err += tmp * vol;
    inf_err = std::max(inf_err, sqrt(tmp));
    unorm += L22(uexact) * vol;
    // std::cout << c << " uh=" << ucalc << " ex=" << uexact << " l2=" << l2_err << std::endl;
  }
#ifdef HAVE_MPI
  GlobalOp("sum", &unorm, 1);
  GlobalOp("sum", &l2_err, 1);
  GlobalOp("max", &inf_err, 1);
#endif
  unorm = sqrt(unorm);
  l2_err = sqrt(l2_err);
}


/* ******************************************************************
* Error in pressure
****************************************************************** */
inline
void AnalyticElasticityBase::ScalarCellError(
    Amanzi::CompositeVector& p, double t, double& pnorm,
    double& l2_err, double& inf_err)
{
  pnorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  Epetra_MultiVector& p_cell = *p.ViewComponent("cell");

  for (int c = 0; c < ncells_owned; ++c) {
    const Amanzi::AmanziGeometry::Point& xm = mesh_->cell_centroid(c);
    double area = mesh_->cell_volume(c);

    double pexact = pressure_exact(xm, t);
    double tmp = fabs(p_cell[0][c] - pexact);

    l2_err += tmp * tmp * area;
    inf_err = std::max(inf_err, tmp);
    pnorm += pexact * pexact * area;
    // std::cout << c << " ph=" << p_cell[0][c] << " ex=" << pexact << std::endl;
  }
#ifdef HAVE_MPI
  GlobalOp("sum", &pnorm, 1);
  GlobalOp("sum", &l2_err, 1);
  GlobalOp("max", &inf_err, 1);
#endif
  pnorm = sqrt(pnorm);
  l2_err = sqrt(l2_err);
}


/* ******************************************************************
* Collective communications.
****************************************************************** */
inline
void AnalyticElasticityBase::GlobalOp(std::string op, double* val, int n)
{
  double* val_tmp = new double[n];
  for (int i = 0; i < n; ++i) val_tmp[i] = val[i];

  if (op == "sum") 
    mesh_->get_comm()->SumAll(val_tmp, val, n);
  else if (op == "max") 
    mesh_->get_comm()->MaxAll(val_tmp, val, n);

  delete[] val_tmp;
}

#endif

