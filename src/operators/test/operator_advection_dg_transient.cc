/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  DG methods for linear advection equations.
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "UnitTest++.h"

// Amanzi
#include "Explicit_TI_RK.hh"
#include "CompositeVector.hh"
#include "DG_Modal.hh"
#include "MeshFactory.hh"
#include "MeshMapsFactory.hh"
#include "NumericalIntegration.hh"
#include "OperatorUtils.hh"
#include "OutputXDMF.hh"
#include "Tensor.hh"
#include "VectorObjects.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "LimiterCellDG.hh"
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "PDE_AdvectionRiemann.hh"
#include "PDE_Reaction.hh"
#include "ReconstructionCellLinear.hh"

#include "AnalyticDG02b.hh"
#include "AnalyticDG06.hh"
#include "AnalyticDG06b.hh"
#include "AnalyticDG06c.hh"
#include "AnalyticDG07.hh"
#include "AnalyticDG07b.hh"
#include "AnalyticDG08.hh"
#include "MeshDeformation.hh"

// global variables
bool exact_solution_expected = false;

namespace Amanzi {

/* *****************************************************************
* Base class: analytic velocities
***************************************************************** */
template <class Analytic>
class AdvectionFn : public Explicit_TI::fnBase<CompositeVector> {
 public:
  AdvectionFn(Teuchos::ParameterList& plist, int nx, double dt0,
              const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
              Teuchos::RCP<WhetStone::DG_Modal> dg,  
              bool conservative_form, std::string weak_form);

  // functional in dy/dt = F(y)
  void FunctionalTimeDerivative(double t, const CompositeVector& u, CompositeVector& f) override;

  // modify time step
  void set_dt(double dt) { dt_ = dt; }

  // calculate cell-center and face-centered velocities using analytic formulas
  virtual void ComputeVelocities(
      double t, double dt, const CompositeVector& u,
      const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> >& velc,
      const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& velf,
      const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& divc);

  // limit solution
  void ApplyLimiter(std::string& limiter, CompositeVector& u);

  // name
  virtual std::string name() { return "analytic velocity"; }

 public:
  double l2norm, dt_stable_min;
  double limiter_min, limiter_mean;

  Teuchos::RCP<Operators::PDE_AdvectionRiemann> op_flux;
  Teuchos::RCP<Operators::PDE_Abstract> op_adv;
  Teuchos::RCP<Operators::PDE_Abstract> op_mass;
  Teuchos::RCP<Operators::PDE_Abstract> op_reac;

 protected:
  Analytic ana_;
  Teuchos::RCP<Operators::Operator> global_op_;
  const Teuchos::ParameterList plist_;

  int nx_;
  double dt_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int ncells_wghost_, nfaces_wghost_;

  int order_;
  Teuchos::RCP<WhetStone::DG_Modal> dg_;  

  double weak_sign_;
  bool conservative_form_, divergence_term_, setup_;
  std::string pk_name_;
};


/* *****************************************************************
* Advection: velocities are computed using L2 or H1 projection
***************************************************************** */
template <class Analytic>
class AdvectionFn_Projection : public AdvectionFn<Analytic> {
 public:
  AdvectionFn_Projection(Teuchos::ParameterList& plist, int nx, double dt0,
                         const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                         Teuchos::RCP<WhetStone::DG_Modal> dg,  
                         bool conservative_form, std::string weak_form)
    : AdvectionFn<Analytic>(plist, nx, dt0, mesh, dg,  conservative_form, weak_form) {};

  // calculate cell-center and face-centered velocities using L2 or H1 projection
  virtual void ComputeVelocities(
      double t, double dt, const CompositeVector& u,
      const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> >& velc,
      const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& velf,
      const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& divc) override;

  // name
  virtual std::string name() override { return "high order"; }

 public:
  using AdvectionFn<Analytic>::plist_;
  using AdvectionFn<Analytic>::order_;

  using AdvectionFn<Analytic>::mesh_;
  using AdvectionFn<Analytic>::ncells_wghost_;
  using AdvectionFn<Analytic>::nfaces_wghost_;

  using AdvectionFn<Analytic>::divergence_term_;

  using AdvectionFn<Analytic>::ana_;
  using AdvectionFn<Analytic>::dt_stable_min;
};


/* *****************************************************************
* Advection: velocities are computed using L2 or H1 projection
***************************************************************** */
template <class Analytic>
class AdvectionFn_LevelSet : public AdvectionFn<Analytic> {
 public:
  AdvectionFn_LevelSet(Teuchos::ParameterList& plist, int nx, double dt0,
                       const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                       Teuchos::RCP<WhetStone::DG_Modal> dg,  
                       bool conservative_form, std::string weak_form)
    : AdvectionFn<Analytic>(plist, nx, dt0, mesh, dg,  conservative_form, weak_form) {};

  // calculate cell-center and face-centered velocities using L2 or H1 projection
  virtual void ComputeVelocities(
      double t, double dt, const CompositeVector& u,
      const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> >& velc,
      const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& velf,
      const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& divc) override;

  // name
  virtual std::string name() override { return "level set"; }

 public:
  using AdvectionFn<Analytic>::plist_;

  using AdvectionFn<Analytic>::dg_;
  using AdvectionFn<Analytic>::order_;

  using AdvectionFn<Analytic>::mesh_;
  using AdvectionFn<Analytic>::ncells_wghost_;
  using AdvectionFn<Analytic>::nfaces_wghost_;

  using AdvectionFn<Analytic>::divergence_term_;

  using AdvectionFn<Analytic>::ana_;
  using AdvectionFn<Analytic>::dt_stable_min;
};


/* *****************************************************************
* Constructor
***************************************************************** */
template <class Analytic>
AdvectionFn<Analytic>::AdvectionFn(
    Teuchos::ParameterList& plist, int nx, double dt, 
    const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
    Teuchos::RCP<WhetStone::DG_Modal> dg,  
    bool conservative_form, std::string weak_form)
    : dt_stable_min(1e+99),
      limiter_min(-1.0),
      limiter_mean(-1.0),
      ana_(mesh, dg->get_order(), true),
      plist_(plist),
      nx_(nx),
      dt_(dt), 
      mesh_(mesh),
      dg_(dg), 
      conservative_form_(conservative_form)
{
  divergence_term_ = !conservative_form_;
  if (weak_form == "dual") {
    weak_sign_ = 1.0;
    pk_name_ = "PK operator";
  } else if (weak_form == "primal") {
    weak_sign_ = -1.0;
    pk_name_ = "PK operator: primal";
    divergence_term_ = false;
  } else if (weak_form == "gauss points") {
    weak_sign_ = -1.0;
    pk_name_ = "PK operator: gauss points";
    divergence_term_ = false;
  }

  // create global operator 
  // -- upwind flux term
  Teuchos::ParameterList op_list = plist.sublist(pk_name_).sublist("flux operator");
  op_flux = Teuchos::rcp(new Operators::PDE_AdvectionRiemann(op_list, mesh));
  global_op_ = op_flux->global_operator();

  // -- volumetric advection term
  op_list = plist.sublist(pk_name_).sublist("advection operator");
  op_adv = Teuchos::rcp(new Operators::PDE_Abstract(op_list, global_op_));

  // -- accumulation term
  op_list = plist.sublist(pk_name_).sublist("inverse mass operator");
  op_mass = Teuchos::rcp(new Operators::PDE_Abstract(op_list, mesh_));

  // -- reaction term
  if (divergence_term_) {
    op_list = plist.sublist(pk_name_).sublist("reaction operator");
    op_reac = Teuchos::rcp(new Operators::PDE_Abstract(op_list, global_op_));
  }

  order_ = dg_->get_order();

  // mesh dimensions
  nfaces_wghost_ = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  ncells_wghost_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  // cotrol variables
  setup_ = true;
}


/* *****************************************************************
* Functional F in dy/dt - F(y) = 0.
***************************************************************** */
template <class Analytic>
void AdvectionFn<Analytic>::FunctionalTimeDerivative(
    double t, const CompositeVector& u, CompositeVector& func)
{
  int d = mesh_->space_dimension();
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // update velocity coefficient
  WhetStone::VectorPolynomial v;
  auto velc = Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_wghost));
  auto velf = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(nfaces_wghost));
  auto divc = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(ncells_wghost));

  ComputeVelocities(t, dt_, u, velc, velf, divc);

  // update problem coefficients
  // -- accumulation
  auto K = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(ncells_wghost));
  WhetStone::Polynomial Kc(d, 0);
  Kc(0, 0) = 1.0;
  for (int c = 0; c < ncells_wghost; c++) {
    (*K)[c] = Kc;
  }

  // -- source term
  int nk = WhetStone::PolynomialSpaceDimension(d, order_);

  WhetStone::Polynomial sol, src, pc(d, order_);
  WhetStone::DenseVector data(pc.size());
  WhetStone::NumericalIntegration numi(mesh_);

  CompositeVector& rhs = *global_op_->rhs();
  Epetra_MultiVector& rhs_c = *rhs.ViewComponent("cell");
  rhs_c.PutScalar(0.0);

  for (int c = 0; c < ncells; ++c) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    ana_.SourceTaylor(xc, t, src);

    for (auto it = pc.begin(); it < pc.end(); ++it) {
      int n = it.PolynomialPosition();

      WhetStone::Polynomial cmono(d, it.multi_index(), 1.0);
      cmono.set_origin(xc);      

      WhetStone::Polynomial tmp = src * cmono;      

      data(n) = numi.IntegratePolynomialCell(c, tmp);
    }

    // -- convert moment to my basis
    dg_->cell_basis(c).LinearFormNaturalToMy(data);
    for (int n = 0; n < pc.size(); ++n) {
      rhs_c[n][c] = data(n);
    }
  }

  // -- boundary data
  auto bc = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::VECTOR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<std::vector<double> >& bc_value = bc->bc_value_vector(nk);

  bool flag;
  WhetStone::Polynomial coefs;

  for (int f = 0; f < nfaces_wghost; f++) {
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 ||
        fabs(xf[1]) < 1e-6 || fabs(xf[1] - 1.0) < 1e-6 ||
        fabs(xf[d - 1]) < 1e-6 || fabs(xf[d - 1] - 1.0) < 1e-6) {
      AmanziGeometry::Point vp = ana_.VelocityExact(xf, t); 

      // const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      const AmanziGeometry::Point& normal = ana_.face_normal_exterior(f, &flag);

      if (vp * normal < -1e-12) {
        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        if (weak_sign_ < 0.0) bc_model[f] = Operators::OPERATOR_BC_DIRICHLET_TYPE2;

        ana_.SolutionTaylor(xf, t, coefs);
        data = coefs.coefs();

        for (int i = 0; i < nk; ++i) {
          bc_value[f][i] = data(i);
        }
      } else if (weak_sign_ < 0.0) {
        bc_model[f] = Operators::OPERATOR_BC_REMOVE;
      }
    }
  }

  // populate the global operator
  op_flux->SetBCs(bc, bc);
  op_flux->Setup(velc, velf);
  op_flux->UpdateMatrices(*velf);
  op_flux->local_op()->Rescale(weak_sign_);
  op_flux->ApplyBCs(true, true, true);

  op_adv->Setup(velc, false);
  op_adv->UpdateMatrices();
  op_adv->local_op()->Rescale(-weak_sign_);

  if (divergence_term_) {
    op_reac->Setup(divc, false);
    op_reac->UpdateMatrices();
    op_reac->local_op()->Rescale(-weak_sign_);
  }

  if (setup_) {
    op_mass->Setup(K, false);
    op_mass->UpdateMatrices(Teuchos::null, Teuchos::null);
    setup_ = false;
  }

  // calculate functional
  CompositeVector u1(u);
  global_op_->ComputeResidual(u, u1);

  // invert vector
  op_mass->global_operator()->Apply(u1, func);

  // statistics: l2 norm of solution
  u.Dot(u, &l2norm);
}


/* *****************************************************************
* Definition of velocities: analytic velocities
***************************************************************** */
template <class Analytic>
void AdvectionFn<Analytic>::ComputeVelocities(
    double t, double dt, const CompositeVector& u,
    const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> >& velc,
    const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& velf,
    const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& divc)
{
  double dt_stable(1e+99), alpha(1.0), vmag;
  WhetStone::VectorPolynomial v;

  for (int c = 0; c < ncells_wghost_; ++c) {
    ana_.VelocityTaylor(mesh_->cell_centroid(c), t, v); 
    (*velc)[c] = v;
  }

  for (int f = 0; f < nfaces_wghost_; ++f) {
    double area = mesh_->face_area(f);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

    ana_.VelocityTaylor(xf, t, v); 
    (*velf)[f] = v * mesh_->face_normal(f);

    v.Value(xf).Norm2(&vmag);
    dt_stable = std::min(dt_stable, area / vmag);
  }
  dt_stable *= alpha / (2 * order_ + 1);
  dt_stable_min = std::min(dt_stable_min, dt_stable);

  if (divergence_term_) {
    for (int c = 0; c < ncells_wghost_; ++c) {
      (*divc)[c] = Divergence((*velc)[c]);
    }
  }
}


/* *****************************************************************
* Change original definitions of velocities: projection algorithm
***************************************************************** */
template <class Analytic>
void AdvectionFn_Projection<Analytic>::ComputeVelocities(
    double t, double dt, const CompositeVector& u,
    const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> >& velc,
    const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& velf,
    const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& divc)
{
  int d = mesh_->space_dimension();
  double dt_stable(1e+99), alpha(1.0), vmag;

  // create a mesh map at time t
  const Teuchos::ParameterList& map_list = plist_.sublist("maps");
  WhetStone::MeshMapsFactory maps_factory;
  auto maps = maps_factory.Create(map_list, mesh_, mesh_);

  // calculate approximate velocities
  AmanziMesh::Entity_ID_List edges;
  int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  double dtfac = 1.0 / dt;
  for (int c = 0; c < ncells_wghost; ++c) {
    const auto& faces = mesh_->cell_get_faces(c);
    int nfaces = faces.size();

    WhetStone::VectorPolynomial v;
    std::vector<WhetStone::VectorPolynomial> vvf, vve;

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      double area = mesh_->face_area(f);
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

      ana_.VelocityTaylor(xf, t, v);
      vvf.push_back(v * dt);
      (*velf)[f] = v * mesh_->face_normal(f);

      v.Value(xf).Norm2(&vmag);
      dt_stable = std::min(dt_stable, area / vmag);
    }

    if (d == 3) {
      mesh_->cell_get_edges(c, &edges);
      int nedges = edges.size();

      for (int n = 0; n < nedges; ++n) {
        int e = edges[n];
        ana_.VelocityTaylor(mesh_->edge_centroid(e), t, v);
        vve.push_back(v * dt);
      }
    }

    maps->VelocityCell(c, vve, vvf, (*velc)[c]);
    (*velc)[c] *= dtfac;

    if (divergence_term_) (*divc)[c] = Divergence((*velc)[c]);
  }

  dt_stable *= alpha / (2 * order_ + 1);
  dt_stable_min = std::min(dt_stable_min, dt_stable);
}


/* *****************************************************************
* Change original definitions of velocities: level set algorithm
***************************************************************** */
template <class Analytic>
void AdvectionFn_LevelSet<Analytic>::ComputeVelocities(
    double t, double dt, const CompositeVector& u,
    const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial> >& velc,
    const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& velf,
    const Teuchos::RCP<std::vector<WhetStone::Polynomial> >& divc)
{
  u.ScatterMasterToGhosted();
  const Epetra_MultiVector& u_c = *u.ViewComponent("cell", true);

  int dim = mesh_->space_dimension();
  int nk = u_c.NumVectors();
  WhetStone::DenseVector data(nk);

  AMANZI_ASSERT(dim == 2);

  int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int nfaces_owned  = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  // cell-based velocity is constant for dGP1
  // we approximate it with a linear function for dGP2
  AmanziMesh::Entity_ID_List faces, cells;
  AmanziGeometry::Point zero(dim);

  // -- normalized cell-centered velocity
  for (int c = 0; c < ncells_wghost; ++c) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    for (int i = 0; i < nk; ++i) data(i) = u_c[i][c];
    dg_->cell_basis(c).ChangeBasisMyToNatural(data);
    WhetStone::Polynomial poly(dim, order_, data); 
    poly.set_origin(xc);

    (*velc)[c] = GradientOnUnitSphere(poly, order_ - 1);

    if (divergence_term_) (*divc)[c] = Divergence((*velc)[c]);
  }

  // -- normalized face-based velocities
  WhetStone::VectorPolynomial vvf(dim, dim, order_ - 1);

  int mk = WhetStone::PolynomialSpaceDimension(dim, order_ - 1);
  auto cvs = Operators::CreateCompositeVectorSpace(mesh_, "face", AmanziMesh::FACE, dim * mk, true);

  CompositeVector vecf(*cvs);
  Epetra_MultiVector vecf_f = *vecf.ViewComponent("face", true);
  vecf.PutScalar(0.0);

  for (int f = 0; f < nfaces_owned; ++f) {
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    WhetStone::Polynomial poly(dim, order_);
    poly.set_origin(xf);

    for (int n = 0; n < ncells; ++n) {
      int c = cells[n];
      for (int i = 0; i < nk; ++i) data(i) = u_c[i][c];
      dg_->cell_basis(c).ChangeBasisMyToNatural(data);
      WhetStone::Polynomial tmp(dim, order_, data); 
      tmp.set_origin(mesh_->cell_centroid(c));

      tmp.ChangeOrigin(xf);
      poly += tmp;
    }

    vvf = GradientOnUnitSphere(poly, order_ - 1);

    for (int i = 0; i < 2; ++i) {
      for (int m = 0; m < mk; ++m) {
        vecf_f[i * mk + m][f] = vvf[i](m);
      }
    }
  } 
    
  vecf.ScatterMasterToGhosted();

  // face-based fluxes scaled by area
  for (int f = 0; f < nfaces_wghost; ++f) {
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    for (int i = 0; i < 2; ++i) {
      for (int m = 0; m < mk; ++m) {
        vvf[i](m) = vecf_f[i * mk + m][f];
      }
    }

    (*velf)[f] = vvf * normal;
    (*velf)[f].set_origin(xf);
  }

  // update CFL condition
  double dt_stable(1.0e+99);
  for (int f = 0; f < nfaces_wghost; ++f) {
    dt_stable = std::min(dt_stable, mesh_->face_area(f));
  }
  dt_stable_min = std::min(dt_stable_min, dt_stable / (2 * order_ + 1));
}


/* *****************************************************************
* Limit gradient
***************************************************************** */
template <class Analytic>
void AdvectionFn<Analytic>::ApplyLimiter(std::string& name, CompositeVector& u)
{
  if (name == "none") return;

  const Epetra_MultiVector& u_c = *u.ViewComponent("cell", true);
  int dim = mesh_->space_dimension();

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost, 0.0);

  // initialize limiter
  Teuchos::ParameterList plist;
  plist.set<std::string>("limiter", name);
  plist.set<int>("polynomial_order", 1);
  plist.set<bool>("limiter extension for transport", false);
  plist.set<std::string>("limiter stencil", "cell to all cells");
  plist.set<int>("limiter points", 3);

  // limit gradient and save it to solution
  Operators::LimiterCellDG limiter(mesh_);
  limiter.Init(plist);

  u.ScatterMasterToGhosted();

  if (name == "Barth-Jespersen dg") {
    limiter.ApplyLimiterDG(u.ViewComponent("cell", true), *dg_, bc_model, bc_value);
  } else {
    int nk = u_c.NumVectors();
    WhetStone::DenseVector data(nk);

    CompositeVectorSpace cvs;
    cvs.SetMesh(mesh_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, dim);
    auto grad = Teuchos::rcp(new CompositeVector(cvs));
    Epetra_MultiVector& grad_c = *grad->ViewComponent("cell", true);

    for (int c = 0; c < ncells_owned; ++c) {
      // mean value is preserved automatically for the partially orthogonalized basis
      // otherwise, a more complicated routine is needed
      AMANZI_ASSERT(dg_->cell_basis(c).id() == WhetStone::TAYLOR_BASIS_NORMALIZED_ORTHO);

      for (int i = 0; i < nk; ++i) data(i) = u_c[i][c];
      for (int i = dim + 1; i < nk; ++i) data(i) = 0.0;

      dg_->cell_basis(c).ChangeBasisMyToNatural(data);
      for (int i = 0; i < dim; ++i) grad_c[i][c] = data(i + 1);
      u_c[0][c] = data(0);
    }

    grad->ScatterMasterToGhosted("cell");

    auto lifting = Teuchos::rcp(new Operators::ReconstructionCellLinear(mesh_, grad));
    limiter.ApplyLimiter(u.ViewComponent("cell", true), 0, lifting, bc_model, bc_value);

    for (int c = 0; c < ncells_owned; ++c) {
      data(0) = u_c[0][c];
      for (int i = 0; i < dim; ++i) data(i + 1) = grad_c[i][c];
      for (int i = dim + 1; i < nk; ++i) data(i) = 0.0;

      dg_->cell_basis(c).ChangeBasisNaturalToMy(data);
      for (int i = 0; i < nk; ++i) u_c[i][c] = data(i);
    }
  }

  // statistics
  auto& lim = *limiter.limiter();
  double tmp1(1.0), tmp2(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    tmp1 = std::min(tmp1, lim[c]);
    tmp2 += lim[c];
  }
  mesh_->get_comm()->MinAll(&tmp1, &limiter_min, 1);
  mesh_->get_comm()->SumAll(&tmp2, &limiter_mean, 1);
  limiter_mean /= lim.Map().MaxAllGID() + 1;
}

}  // namespace Amanzi


/* *****************************************************************
* This tests the transient advection scheme.
***************************************************************** */
// support functions for the level set algorithm
bool inside1(const Amanzi::AmanziGeometry::Point& p) {
  Amanzi::AmanziGeometry::Point c(0.5, 0.5);
  return (norm(p - c) < 0.06); 
}
bool inside2(const Amanzi::AmanziGeometry::Point& p) {
  Amanzi::AmanziGeometry::Point c(1.0, 0.0);
  return (norm(p) < 0.06 || norm(p - c) < 0.06); 
}

// support function for visualization: extrapolation to mesh nodes
Teuchos::RCP<Epetra_MultiVector> InterpolateCellToNode(
    Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh,
    const Amanzi::WhetStone::DG_Modal& dg,
    const Epetra_MultiVector& uc)
{
  int order = dg.get_order();
  int nk = uc.NumVectors();
  Amanzi::AmanziGeometry::Point xv(mesh->space_dimension());
  Amanzi::AmanziMesh::Entity_ID_List cells;

  int nnodes_owned = mesh->num_entities(Amanzi::AmanziMesh::NODE, Amanzi::AmanziMesh::Parallel_type::OWNED);
  auto un = Teuchos::rcp(new Epetra_MultiVector(mesh->node_map(false), 1));

  for (int v = 0; v < nnodes_owned; ++v) {
    mesh->node_get_coordinates(v, &xv);
    mesh->node_get_cells(v, Amanzi::AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    double value(0.0);
    Amanzi::WhetStone::DenseVector data(nk);

    for (int n = 0; n < ncells; ++n) {
      int c = cells[n];
      for (int i = 0; i < nk; ++i) data(i) = uc[i][c];
      auto poly = dg.cell_basis(c).CalculatePolynomial(mesh, c, order, data);
      value += poly.Value(xv);
    }
    (*un)[0][v] = value / ncells;
  } 

  return un;
}

template <class Analytic, class Advection>
void Transient(std::string filename, int nx, int ny, int nz,
               double dt0, double tend,
               bool conservative_form = true, 
               std::string weak_form = "dual",
               std::string limiter = "none")
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  // read parameter list
  std::string xmlFileName = "test/operator_advection_dg_transient.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  std::string pk_name = "PK operator";
  if (weak_form == "primal") pk_name = "PK operator: primal";
  if (weak_form == "gauss points") pk_name = "PK operator: gauss points";

  int order = plist.sublist(pk_name)
                   .sublist("flux operator")
                   .sublist("schema").get<int>("method order");
  Teuchos::ParameterList dg_list = plist.sublist(pk_name)
                                        .sublist("flux operator").sublist("schema");
  std::string basis = dg_list.get<std::string>("dg basis");

  { 
    std::string problem = (conservative_form) ? "conservative" : "non-conservative";
    if (MyPID == 0) {
      std::cout << "\nTest: dG transient advection: " << filename 
                << ", order=" << order << ", dt=" << dt0
                << "\n      PDE=" << problem << ", basis=" << basis 
                << "\n      weak formulation=\"" << weak_form << "\""
                << ", limiter=\"" << limiter << "\"" << std::endl;
    }
  }

  // create a mesh framework
  MeshFactory meshfactory(comm, Teuchos::null);
  // meshfactory.set_partitioner(AmanziMesh::Partitioner_type::ZOLTAN_RCB);
  meshfactory.set_preference(Preference({Framework::MSTK, Framework::STK}));
  RCP<Mesh> mesh;
  if (nx == 0 || ny == 0)
    mesh = meshfactory.create(filename, true, true);
  else if (nz == 0) 
    mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, nx, ny);
  else 
    mesh = meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, nx, ny, nz, true, true);

  // DeformMesh(mesh, deform, 0.0);

  // initialize I/O
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh, true, false);

  // create main advection class
  auto dg = Teuchos::rcp(new WhetStone::DG_Modal(dg_list, mesh));
  Advection fn(plist, nx, dt0, mesh, dg, conservative_form, weak_form);

  if (fn.name() == "high order" && MyPID == 0) {
    const auto& map_list = plist.sublist("maps");
    int vel_order = map_list.get<int>("method order");
    std::string vel_method = map_list.get<std::string>("method");
    std::string vel_projector = map_list.get<std::string>("projector");
      
    std::cout << "      face velocity: order=" << vel_order 
              << ", projector=" << vel_projector 
              << ", method=\"" << vel_method << "\"" << std::endl;
  }

  // create initial guess
  Analytic ana(mesh, order, true);

  CompositeVector& rhs = *fn.op_flux->global_operator()->rhs();
  CompositeVector sol(rhs), sol_next(rhs);
  sol.PutScalar(0.0);

  Epetra_MultiVector& sol_c = *sol.ViewComponent("cell");
  ana.InitialGuess(*dg, sol_c, 0.0);

  int nstep(0);
  double dt(dt0), t(0.0), tio(dt0);
  Explicit_TI::RK<CompositeVector> rk(fn, Explicit_TI::tvd_3rd_order, sol);

  while(std::fabs(t < tend) < dt/4 || dt > 1e-12) {
    fn.set_dt(dt);
    fn.ApplyLimiter(limiter, sol);
    rk.TimeStep(t, dt, sol, sol_next);

    sol = sol_next;
    t += dt;
    nstep++;

    // visualization
    if (std::fabs(t - tio) < dt/4) {
      tio = std::min(tio + 0.1, tend); 
      ana.GlobalOp("min", &fn.dt_stable_min, 1);
      if (MyPID == 0)
        printf("t=%9.6f |p|=%12.8g  CFL=%8.6f  limiter min/mean: %8.4g %8.4g\n",
            t, fn.l2norm, fn.dt_stable_min, fn.limiter_min, fn.limiter_mean);

      const Epetra_MultiVector& pc = *sol.ViewComponent("cell");
      auto pn = InterpolateCellToNode(mesh, *dg, pc);
  
      io.InitializeCycle(t, nstep, "");
      io.WriteVector(*pc(0), "solution", AmanziMesh::CELL);
      io.WriteVector(*(*pn)(0), "interpolation", AmanziMesh::NODE);
      io.FinalizeCycle();
    }

    dt = std::min(dt0, tend - t);

    // overwrite solution at the origin
    if (fn.name() == "level set") {
      ana.InitialGuess(*dg, sol_c, t, inside1);
    }
  }

  // compute solution error
  sol.ScatterMasterToGhosted();
  Epetra_MultiVector& p = *sol.ViewComponent("cell", false);

  double pnorm, pl2_err, pinf_err, pl2_mean, pinf_mean, pl2_int;
  ana.ComputeCellError(*dg, p, tend, pnorm, pl2_err, pinf_err, pl2_mean, pinf_mean, pl2_int);

  double pface_inf, grad_pface_inf;
  ana.ComputeFaceError(*dg, p, tend, pface_inf, grad_pface_inf);

  if (MyPID == 0) {
    printf("nx=%3d CELL: (mean) L2(p)=%9.6g  Inf(p)=%9.6g\n", nx, pl2_mean, pinf_mean);
    printf("            (total) L2(p)=%9.6g  Inf(p)=%9.6g\n", pl2_err, pinf_err);
    printf("         (integral) L2(p)=%9.6g\n", pl2_int);
    printf("       FACE:        Inf(p)=%9.6g  Inf(grad p)=%9.6g\n", pface_inf, grad_pface_inf);
    if (exact_solution_expected) 
      CHECK(pl2_mean < 1e-10);
    else if (fn.name() == "level set") 
      CHECK(pl2_mean < 0.3 / nx);
    else if (limiter == "none") 
      CHECK(pl2_mean < 0.12 / nx);
    else
      CHECK(pl2_mean < 0.3 / nx);
  }
}


TEST(OPERATOR_ADVECTION_TRANSIENT_DG) {
  double dT(0.1), T1(1.0);
  exact_solution_expected = true;
  Transient<AnalyticDG02b, Amanzi::AdvectionFn_Projection<AnalyticDG02b> >("square", 4,4,0, dT,T1, false);
  Transient<AnalyticDG02b, Amanzi::AdvectionFn_Projection<AnalyticDG02b> >("square", 4,4,0, dT,T1, false, "gauss points");

  exact_solution_expected = false;
  Transient<AnalyticDG06b, Amanzi::AdvectionFn_Projection<AnalyticDG06b> >("square", 4,4,0, dT,T1);
  Transient<AnalyticDG06,  Amanzi::AdvectionFn_Projection<AnalyticDG06> > ("square", 4,4,0, dT,T1, false, "primal");
  Transient<AnalyticDG06,  Amanzi::AdvectionFn_Projection<AnalyticDG06> > ("square", 4,4,0, dT,T1, false, "dual");
  Transient<AnalyticDG06,  Amanzi::AdvectionFn_Projection<AnalyticDG06> > ("square", 4,4,0, dT,T1, false, "gauss points", "Barth-Jespersen dg");
  Transient<AnalyticDG06c, Amanzi::AdvectionFn_Projection<AnalyticDG06c> >("cube",   2,2,1, dT,T1);

  dT = 0.01;
  Transient<AnalyticDG07, Amanzi::AdvectionFn_LevelSet<AnalyticDG07> >("square", 6, 6, 0, dT,  T1, false, "primal");

  /*
  double dT(0.05), T1(1.0);
  Transient<AnalyticDG06c, Amanzi::AdvectionFn_Projection<AnalyticDG06c> >("cube", 8, 8, 8, dT / 2, T1);
  Transient<AnalyticDG06c, Amanzi::AdvectionFn_Projection<AnalyticDG06c> >("cube",16,16,16, dT / 4, T1);
  */

  /*
  double dT(0.01), T1(1.0);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("square",  16, 16, 0, dT,   T1);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("square",  32, 32, 0, dT/2, T1);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("square",  64, 64, 0, dT/4, T1);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("square", 128,128, 0, dT/8, T1);
  */

  /*
  double dT(0.01), T1(1.0);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("test/triangular8.exo",    8,0,0, dT,   T1);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("test/triangular16.exo",  16,0,0, dT/2, T1);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("test/triangular32.exo",  32,0,0, dT/4, T1);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("test/triangular64.exo",  64,0,0, dT/8, T1);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("test/triangular128.exo",128,0,0, dT/16,T1);
  */

  /*
  double dT(0.002), T1(1.0);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("test/median15x16.exo",   16,0,0, dT,  T1);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("test/median32x33.exo",   32,0,0, dT/2,T1);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("test/median63x64.exo",   64,0,0, dT/4,T1);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("test/median127x128.exo",128,0,0, dT/8,T1);
  */

  /*
  double dT(0.01), T1(1.0);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("test/mesh_poly20x20.exo",   20,0,0, dT,  T1);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("test/mesh_poly40x40.exo",   40,0,0, dT/2,T1);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("test/mesh_poly80x80.exo",   80,0,0, dT/4,T1);
  Transient<AnalyticDG06, Amanzi::AdvectionFn_Projection<AnalyticDG06> >("test/mesh_poly160x160.exo",160,0,0, dT/8,T1);
  */

  /*
  double dT(0.001), T1(1.0);
  Transient<AnalyticDG07, Amanzi::AdvectionFn_LevelSet<AnalyticDG07> >("square", 16, 16, 0, dT,  T1, false, "primal");
  Transient<AnalyticDG07, Amanzi::AdvectionFn_LevelSet<AnalyticDG07> >("square", 32, 32, 0, dT/2,T1, false, "primal");
  Transient<AnalyticDG07, Amanzi::AdvectionFn_LevelSet<AnalyticDG07> >("square", 64, 64, 0, dT/4,T1, false, "primal");
  Transient<AnalyticDG07, Amanzi::AdvectionFn_LevelSet<AnalyticDG07> >("square",128,128, 0, dT/8,T1, false, "primal");
  */

  /*
  double dT(0.001), T1(1.0);
  Transient<AnalyticDG07, Amanzi::AdvectionFn_LevelSet<AnalyticDG07> >("test/median15x16.exo",   16,0,0, dT,  T1, false, "primal");
  Transient<AnalyticDG07, Amanzi::AdvectionFn_LevelSet<AnalyticDG07> >("test/median32x33.exo",   32,0,0, dT/2,T1, false, "primal");
  Transient<AnalyticDG07, Amanzi::AdvectionFn_LevelSet<AnalyticDG07> >("test/median63x64.exo",   64,0,0, dT/4,T1, false, "primal");
  Transient<AnalyticDG07, Amanzi::AdvectionFn_LevelSet<AnalyticDG07> >("test/median127x128.exo",128,0,0, dT/8,T1, false, "primal");
  */

  /*
  double dT(0.001), T1(0.8);
  Transient<AnalyticDG07b, Amanzi::AdvectionFn_LevelSet<AnalyticDG07b> >("square", 128, 128, 0, dT/6,T1, false, "primal");
  Transient<AnalyticDG07b, Amanzi::AdvectionFn_LevelSet<AnalyticDG07b> >("test/median127x128.exo", 128,0,0, dT/6,T1, false, "primal");
  */

  /*
  double dT(0.01), T1(6.2832);
  Transient<AnalyticDG08, Amanzi::AdvectionFn_Projection<AnalyticDG08> >("square", 128,128,0, dT/8,T1, true, "dual", "none");
  Transient<AnalyticDG08, Amanzi::AdvectionFn_Projection<AnalyticDG08> >("test/median127x128.exo", 128,0,0, dT/8,T1, true, "dual", "none");
  */
}


