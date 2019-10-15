/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

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
#include "OutputXDMF.hh"
#include "CompositeVector.hh"
#include "DG_Modal.hh"
#include "LinearOperatorGMRES.hh"
#include "MeshFactory.hh"
#include "MeshMapsFactory.hh"
#include "NumericalIntegration.hh"
#include "OperatorUtils.hh"
#include "Tensor.hh"
#include "VectorPolynomial.hh"

// Operators
#include "AnalyticDG02b.hh"
#include "AnalyticDG06.hh"
#include "AnalyticDG06b.hh"
#include "AnalyticDG07.hh"
#include "AnalyticDG07b.hh"
#include "AnalyticDG08.hh"

#include "LimiterCell.hh"
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "PDE_AdvectionRiemann.hh"
#include "PDE_Reaction.hh"

// global variables
bool exact_solution_expected = false;

namespace Amanzi {

template <class AnalyticDG>
class AdvectionFn : public Explicit_TI::fnBase<CompositeVector> {
 public:
  AdvectionFn(Teuchos::ParameterList& plist, int nx, double dt0,
              const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
              Teuchos::RCP<WhetStone::DG_Modal> dg, bool conservative_form,
              std::string weak_form);

  // functional in dy/dt = F(y)
  void FunctionalTimeDerivative(double t, const CompositeVector& u,
                                CompositeVector& f) override;

  // modify time step
  void set_dt(double dt) { dt_ = dt; }

  // modify cell and face velocities using a mesh map
  void ApproximateVelocity_Projection(
    double t, double dt,
    const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial>>& velc,
    const Teuchos::RCP<std::vector<WhetStone::Polynomial>>& velf,
    const Teuchos::RCP<std::vector<WhetStone::Polynomial>>& divc);

  void ApproximateVelocity_LevelSet(
    double t, const CompositeVector& u,
    const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial>>& velc,
    const Teuchos::RCP<std::vector<WhetStone::Polynomial>>& velf,
    const Teuchos::RCP<std::vector<WhetStone::Polynomial>>& divc);

  // limit solution
  void ApplyLimiter(std::string& limiter, CompositeVector& u);

 public:
  double l2norm;
  double limiter_min, limiter_mean;

  Teuchos::RCP<Operators::PDE_AdvectionRiemann> op_flux;
  Teuchos::RCP<Operators::PDE_Abstract> op_adv;
  Teuchos::RCP<Operators::PDE_Abstract> op_mass;
  Teuchos::RCP<Operators::PDE_Abstract> op_reac;

 private:
  AnalyticDG ana_;
  Teuchos::RCP<Operators::Operator> global_op_;
  const Teuchos::ParameterList plist_;

  int nx_;
  double dt_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<AmanziMesh::Mesh> mesh_new_;

  int order_;
  Teuchos::RCP<WhetStone::DG_Modal> dg_;

  double weak_sign_;
  bool conservative_form_, divergence_term_;
  bool setup_, high_order_velf_, level_set_velf_;
  std::string pk_name_;
};


/* *****************************************************************
 * Constructor
 ***************************************************************** */
template <class AnalyticDG>
AdvectionFn<AnalyticDG>::AdvectionFn(
  Teuchos::ParameterList& plist, int nx, double dt,
  const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
  Teuchos::RCP<WhetStone::DG_Modal> dg, bool conservative_form,
  std::string weak_form)
  : plist_(plist),
    nx_(nx),
    dt_(dt),
    mesh_(mesh),
    dg_(dg),
    ana_(mesh, dg->order(), true),
    conservative_form_(conservative_form),
    limiter_min(-1.0),
    limiter_mean(-1.0)
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
  Teuchos::ParameterList op_list =
    plist.sublist(pk_name_).sublist("flux operator");
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

  order_ = dg_->order();

  // create auxiliary mesh
  auto comm = Amanzi::getDefaultComm();
  AmanziMesh::MeshFactory meshfactory(comm, mesh_->geometric_model());
  // factory.set_partitioner(AmanziMesh::Partitioner_type::ZOLTAN_RCB);
  meshfactory.set_preference(
    AmanziMesh::Preference({ AmanziMesh::Framework::MSTK }));

  std::string name = plist.get<std::string>("file name");
  if (name == "square")
    mesh_new_ = meshfactory.create(0.0, 0.0, 1.0, 1.0, nx_, nx_);
  else
    mesh_new_ = meshfactory.create(name, true, true);

  // cotrol variables
  name = plist.get<std::string>("face velocity method");
  high_order_velf_ = (name == "high order");
  level_set_velf_ = (name == "level set");

  setup_ = true;
}


/* *****************************************************************
 * Functional F in dy/dt = F(y)
 ***************************************************************** */
template <class AnalyticDG>
void
AdvectionFn<AnalyticDG>::FunctionalTimeDerivative(double t,
                                                  const CompositeVector& u,
                                                  CompositeVector& f)
{
  int d = mesh_->space_dimension();
  int ncells =
    mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost =
    mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost =
    mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // update velocity coefficient
  WhetStone::VectorPolynomial v;
  auto velc =
    Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_wghost));
  auto velf =
    Teuchos::rcp(new std::vector<WhetStone::Polynomial>(nfaces_wghost));
  auto divc =
    Teuchos::rcp(new std::vector<WhetStone::Polynomial>(ncells_wghost));

  for (int c = 0; c < ncells_wghost; ++c) {
    ana_.VelocityTaylor(mesh_->cell_centroid(c), t, v);
    (*velc)[c] = v;
    (*velc)[c] *= -weak_sign_;
  }

  for (int f = 0; f < nfaces_wghost; ++f) {
    ana_.VelocityTaylor(mesh_->face_centroid(f), t, v);
    (*velf)[f] = v * (mesh_->face_normal(f) * weak_sign_);
  }

  if (divergence_term_) {
    for (int c = 0; c < ncells_wghost; ++c) {
      (*divc)[c] = Divergence((*velc)[c]);
    }
  }

  // modify analytic Taylor expansions
  if (level_set_velf_) {
    ApproximateVelocity_LevelSet(t, u, velc, velf, divc);
  } else {
    ApproximateVelocity_Projection(t, dt_, velc, velf, divc);
  }

  // update problem coefficients
  // -- accumulation
  auto K = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(ncells_wghost));
  WhetStone::Polynomial Kc(d, 0);
  Kc(0, 0) = 1.0;
  for (int c = 0; c < ncells_wghost; c++) { (*K)[c] = Kc; }

  // -- source term
  int nk = (order_ + 1) * (order_ + 2) / 2;

  WhetStone::Polynomial sol, src, pc(2, order_);
  WhetStone::DenseVector data(pc.size());
  WhetStone::NumericalIntegration numi(mesh_);

  CompositeVector& rhs = *global_op_->rhs();
  Epetra_MultiVector& rhs_c = *rhs.ViewComponent("cell");
  rhs_c.putScalar(0.0);

  for (int c = 0; c < ncells; ++c) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    double volume = mesh_->cell_volume(c, false);

    ana_.SourceTaylor(xc, t, src);

    for (auto it = pc.begin(); it < pc.end(); ++it) {
      int n = it.PolynomialPosition();
      int k = it.MonomialSetOrder();

      WhetStone::Polynomial cmono(2, it.multi_index(), 1.0);
      cmono.set_origin(xc);

      WhetStone::Polynomial tmp = src * cmono;

      data(n) = numi.IntegratePolynomialCell(c, tmp);
    }

    // -- convert moment to my basis
    dg_->cell_basis(c).LinearFormNaturalToMy(data);
    for (int n = 0; n < pc.size(); ++n) { rhs_c[n][c] = data(n); }
  }

  // -- boundary data
  auto bc = Teuchos::rcp(
    new Operators::BCs(mesh_, AmanziMesh::FACE, Operators::DOF_Type::VECTOR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<std::vector<double>>& bc_value = bc->bc_value_vector(nk);

  WhetStone::Polynomial coefs;

  for (int f = 0; f < nfaces_wghost; f++) {
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    if (fabs(xf[0]) < 1e-6 || fabs(xf[0] - 1.0) < 1e-6 || fabs(xf[1]) < 1e-6 ||
        fabs(xf[1] - 1.0) < 1e-6) {
      AmanziGeometry::Point vp = ana_.VelocityExact(xf, t);

      if (vp * normal < -1e-12) {
        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        if (weak_sign_ < 0.0)
          bc_model[f] = Operators::OPERATOR_BC_DIRICHLET_TYPE2;

        ana_.SolutionTaylor(xf, t, coefs);
        data = coefs.coefs();

        for (int i = 0; i < nk; ++i) { bc_value[f][i] = data(i); }
      } else if (weak_sign_ < 0.0) {
        bc_model[f] = Operators::OPERATOR_BC_REMOVE;
      }
    }
  }

  // populate the global operator
  op_flux->SetBCs(bc, bc);
  op_flux->Setup(velc, velf);
  op_flux->UpdateMatrices(velf.ptr());
  op_flux->ApplyBCs(true, true, true);

  op_adv->SetupPolyVector(velc);
  op_adv->UpdateMatrices();

  if (divergence_term_) {
    op_reac->SetupPoly(divc);
    op_reac->UpdateMatrices();
  }

  if (setup_) {
    op_mass->SetupPoly(K);
    op_mass->UpdateMatrices(Teuchos::null);
    setup_ = false;
  }

  // calculate functional
  CompositeVector u1(u);
  global_op_->ComputeResidual(u, u1);

  // invert vector
  op_mass->global_operator()->apply(u1, f);

  // statistics: l2 norm of solution
  l2norm = u.dot(u);
}


/* *****************************************************************
 * Change original definitions of velocities: projection algorithm
 ***************************************************************** */
template <class AnalyticDG>
void
AdvectionFn<AnalyticDG>::ApproximateVelocity_Projection(
  double t, double dt,
  const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial>>& velc,
  const Teuchos::RCP<std::vector<WhetStone::Polynomial>>& velf,
  const Teuchos::RCP<std::vector<WhetStone::Polynomial>>& divc)
{
  // deform auxiliary mesh for time interval (t, t + dt)
  int dim = mesh_->space_dimension();
  AmanziGeometry::Point xv(dim), yv(dim), sv(dim);
  AmanziMesh::Entity_ID_List faces, nodeids;
  AmanziGeometry::Point_List new_positions, final_positions;

  int ncells_wghost =
    mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nnodes_wghost =
    mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  for (int v = 0; v < nnodes_wghost; ++v) {
    mesh_->node_get_coordinates(v, &xv);
    yv = ana_.VelocityExact(xv, t);
    xv += yv * dt;

    nodeids.push_back(v);
    new_positions.push_back(xv);
  }
  mesh_new_->deform(nodeids, new_positions, false, &final_positions);

  // create a mesh map at time t
  const Teuchos::ParameterList& map_list = plist_.sublist("maps");
  WhetStone::MeshMapsFactory maps_factory;
  auto maps = maps_factory.Create(map_list, mesh_, mesh_new_);

  // calculate approximate velocities
  double dtfac = 1.0 / dt;
  for (int c = 0; c < ncells_wghost; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    WhetStone::VectorPolynomial v;
    std::vector<WhetStone::VectorPolynomial> vvf;

    if (high_order_velf_) { // use high-order face velocity
      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        ana_.VelocityTaylor(mesh_->face_centroid(f), t, v);
        vvf.push_back(v * dt);
        (*velf)[f] = v * (mesh_->face_normal(f) * weak_sign_);
      }
    } else {
      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        maps->VelocityFace(f, v);
        vvf.push_back(v);
        (*velf)[f] = (v * mesh_->face_normal(f)) * dtfac * weak_sign_;
      }
    }

    maps->VelocityCell(c, vvf, (*velc)[c]);
    (*velc)[c] *= -dtfac * weak_sign_;

    if (divergence_term_) (*divc)[c] = Divergence((*velc)[c]);
  }
}


/* *****************************************************************
 * Change original definitions of velocities: level set algorithm
 ***************************************************************** */
template <class AnalyticDG>
void
AdvectionFn<AnalyticDG>::ApproximateVelocity_LevelSet(
  double t, const CompositeVector& u,
  const Teuchos::RCP<std::vector<WhetStone::VectorPolynomial>>& velc,
  const Teuchos::RCP<std::vector<WhetStone::Polynomial>>& velf,
  const Teuchos::RCP<std::vector<WhetStone::Polynomial>>& divc)
{
  u.ScatterMasterToGhosted();
  const Epetra_MultiVector& u_c = *u.ViewComponent("cell", true);

  int dim = mesh_->space_dimension();
  int nk = u_c.getNumVectors();
  WhetStone::DenseVector data(nk);

  int ncells_wghost =
    mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost =
    mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int nfaces_owned =
    mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  // cell-based velocity is constant for dGP1
  // we approximate it with a linear function for dGP2
  AmanziMesh::Entity_ID_List faces, cells;
  AmanziGeometry::Point zero(dim);

  // -- normalized cell-centered velocity
  for (int c = 0; c < ncells_wghost; ++c) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);

    // ana_.VelocityTaylor(xc, t, (*velc)[c]);
    // (*velc)[c] *= -weak_sign_;  // velocity is negative of normalized
    // gradient

    for (int i = 0; i < nk; ++i) data(i) = u_c[i][c];
    dg_->cell_basis(c).ChangeBasisMyToNatural(data);
    WhetStone::Polynomial poly(dim, order_, data);
    poly.set_origin(xc);

    poly *= weak_sign_;
    (*velc)[c] = GradientOnUnitSphere(poly, order_ - 1);

    if (divergence_term_) (*divc)[c] = Divergence((*velc)[c]);
  }

  // -- normalized face-based velocities
  WhetStone::VectorPolynomial vvf(dim, dim, order_ - 1);

  int mk = WhetStone::PolynomialSpaceDimension(dim, order_ - 1);
  auto cvs = Operators::CreateCompositeVectorSpace(
    mesh_, "face", AmanziMesh::FACE, dim * mk, true);

  CompositeVector vecf(*cvs);
  Epetra_MultiVector vecf_f = *vecf.ViewComponent("face", true);
  vecf.putScalar(0.0);

  for (int f = 0; f < nfaces_owned; ++f) {
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

    // ana_.VelocityTaylor(xf, t, vvf);
    // vvf *= weak_sign_;  // velocity is negative of normalized gradient

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
      poly -= tmp;
    }

    poly *= weak_sign_;
    vvf = GradientOnUnitSphere(poly, order_ - 1);

    for (int i = 0; i < 2; ++i) {
      for (int m = 0; m < mk; ++m) { vecf_f[i * mk + m][f] = vvf[i](m); }
    }
  }

  vecf.ScatterMasterToGhosted();

  // face-based fluxes scaled by area
  for (int f = 0; f < nfaces_wghost; ++f) {
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    double area = mesh_->face_area(f);

    for (int i = 0; i < 2; ++i) {
      for (int m = 0; m < mk; ++m) { vvf[i](m) = vecf_f[i * mk + m][f]; }
    }

    (*velf)[f] = vvf * normal;
    (*velf)[f].set_origin(xf);
  }
}


/* *****************************************************************
 * Limit gradient
 ***************************************************************** */
template <class AnalyticDG>
void
AdvectionFn<AnalyticDG>::ApplyLimiter(std::string& name, CompositeVector& u)
{
  if (name == "none") return;

  const Epetra_MultiVector& u_c = *u.ViewComponent("cell", true);
  int dim = mesh_->space_dimension();

  int ncells_owned =
    mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost =
    mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

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
  Operators::LimiterCell limiter(mesh_);
  limiter.Init(plist);

  u.ScatterMasterToGhosted();

  if (name == "Barth-Jespersen dg") {
    limiter.ApplyLimiter(
      u.ViewComponent("cell", true), *dg_, bc_model, bc_value);
  } else {
    int nk = u_c.getNumVectors();
    WhetStone::DenseVector data(nk);

    CompositeVectorSpace cvs;
    cvs.SetMesh(mesh_)->SetGhosted(true)->AddComponent(
      "cell", AmanziMesh::CELL, dim);
    auto grad = Teuchos::rcp(new CompositeVector(cvs));
    Epetra_MultiVector& grad_c = *grad->ViewComponent("cell", true);

    for (int c = 0; c < ncells_owned; ++c) {
      // mean value is preserved automatiacally for the partially orthogonalized
      // basis otherwise, a more complicated routine is needed
      AMANZI_ASSERT(dg_->cell_basis(c).id() ==
                    WhetStone::TAYLOR_BASIS_NORMALIZED_ORTHO);

      for (int i = 0; i < nk; ++i) data(i) = u_c[i][c];
      for (int i = dim + 1; i < nk; ++i) data(i) = 0.0;

      dg_->cell_basis(c).ChangeBasisMyToNatural(data);
      for (int i = 0; i < dim; ++i) grad_c[i][c] = data(i + 1);
      u_c[0][c] = data(0);
    }

    grad->ScatterMasterToGhosted("cell");

    limiter.ApplyLimiter(
      u.ViewComponent("cell", true), 0, grad, bc_model, bc_value);

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
  Teuchos::reduceAll(
    *mesh_->get_comm(), Teuchos::REDUCE_MIN, 1, &tmp1, &limiter_min);
  Teuchos::reduceAll(
    *mesh_->get_comm(), Teuchos::REDUCE_SUM, 1, &tmp2, &limiter_mean);
  limiter_mean /= lim.getMap().getMaxGlobalIndex() + 1;
}

} // namespace Amanzi


/* *****************************************************************
 * This tests exactness of the transient advection scheme for
 * dp/dt + div(v p) = f.
 ***************************************************************** */
bool
inside1(const Amanzi::AmanziGeometry::Point& p)
{
  Amanzi::AmanziGeometry::Point c(0.5, 0.5);
  return (norm(p - c) < 0.06);
}
bool
inside2(const Amanzi::AmanziGeometry::Point& p)
{
  Amanzi::AmanziGeometry::Point c(1.0, 0.0);
  return (norm(p) < 0.06 || norm(p - c) < 0.06);
}

template <class AnalyticDG>
void
AdvectionTransient(std::string filename, int nx, int ny, double dt0,
                   double tend, const Amanzi::Explicit_TI::method_t& rk_method,
                   bool conservative_form = true,
                   std::string weak_form = "dual",
                   std::string face_velocity_method = "high order",
                   std::string limiter = "none")
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int getRank = comm->getRank();

  // read parameter list
  std::string xmlFileName = "test/operator_advection_dg_transient.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  ParameterList plist = xmlreader.getParameters();

  std::string pk_name = "PK operator";
  if (weak_form == "primal") pk_name = "PK operator: primal";
  if (weak_form == "gauss points") pk_name = "PK operator: gauss points";

  int order =
    plist.sublist(pk_name).sublist("flux operator").get<int>("method order");
  Teuchos::ParameterList dg_list =
    plist.sublist(pk_name).sublist("flux operator");
  std::string basis = dg_list.get<std::string>("dg basis");

  {
    std::string problem = (conservative_form) ? ", conservative PDE" : "";
    if (getRank == 0) {
      std::cout << "\nTest: 2D dG transient advection: " << filename
                << ", order=" << order << problem << ", basis=" << basis
                << "\n      weak formulation=\"" << weak_form << "\""
                << ", face_velocity=\"" << face_velocity_method << "\""
                << ", limiter=\"" << limiter << "\"" << std::endl;
      if (face_velocity_method == "high order") {
        const auto& map_list = plist.sublist("maps");
        int vel_order = map_list.get<int>("method order");
        std::string vel_method = map_list.get<std::string>("method");
        std::string vel_projector = map_list.get<std::string>("projector");

        std::cout << "      face velocity: order=" << vel_order
                  << ", projector=" << vel_projector << ", method=\""
                  << vel_method << "\"" << std::endl;
      }
    }
  }

  // create a mesh framework
  ParameterList region_list = plist.sublist("regions");
  Teuchos::RCP<GeometricModel> gm =
    Teuchos::rcp(new GeometricModel(2, region_list, *comm));

  MeshFactory meshfactory(comm, gm);
  // meshfactory.set_partitioner(AmanziMesh::Partitioner_type::ZOLTAN_RCB);
  meshfactory.set_preference(Preference({ Framework::MSTK, Framework::STK }));
  RCP<const Mesh> mesh;
  if (nx == 0 || ny == 0)
    mesh = meshfactory.create(filename, true, true);
  else
    mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, nx, ny);

  int ncells =
    mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces =
    mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost =
    mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost =
    mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // initialize I/O
  Teuchos::ParameterList iolist;
  iolist.get<std::string>("file name base", "plot");
  OutputXDMF io(iolist, mesh, true, false);

  // create main advection class
  plist.set<std::string>("file name", filename);
  plist.set<std::string>("face velocity method", face_velocity_method);

  auto dg = Teuchos::rcp(new WhetStone::DG_Modal(dg_list, mesh));
  AdvectionFn<AnalyticDG> fn(
    plist, nx, dt0, mesh, dg, conservative_form, weak_form);

  // create initial guess
  AnalyticDG ana(mesh, order, true);

  CompositeVector& rhs = *fn.op_flux->global_operator()->rhs();
  CompositeVector sol(rhs), sol_next(rhs);
  sol.putScalar(0.0);

  Epetra_MultiVector& sol_c = *sol.ViewComponent("cell");
  ana.InitialGuess(*dg, sol_c, 0.0);

  int nstep(0);
  double dt(dt0), t(0.0), tio(dt0);
  Explicit_TI::RK<CompositeVector> rk(fn, rk_method, sol);

  while (std::fabs(t < tend) < dt / 4 || dt > 1e-12) {
    fn.set_dt(dt);
    fn.ApplyLimiter(limiter, sol);
    rk.TimeStep(t, dt, sol, sol_next);

    sol = sol_next;
    t += dt;
    nstep++;

    // visualization
    if (std::fabs(t - tio) < dt / 4) {
      tio = std::min(tio + 0.1, tend);
      if (getRank == 0)
        printf("t=%9.6f |p|=%12.8g limiter min/mean: %8.4g %8.4g\n",
               t,
               fn.l2norm,
               fn.limiter_min,
               fn.limiter_mean);

      const Epetra_MultiVector& p = *sol.ViewComponent("cell");

      io.InitializeCycle(t, nstep);
      io.WriteVector(*p(0), "solution");
      io.FinalizeCycle();
    }

    dt = std::min(dt0, tend - t);

    // overwrite solution at the origin
    if (face_velocity_method == "level set") {
      ana.InitialGuess(*dg, sol_c, t, inside2);
    }
  }

  // compute solution error
  sol.ScatterMasterToGhosted();
  Epetra_MultiVector& p = *sol.ViewComponent("cell", false);

  double pnorm, pl2_err, pinf_err, pl2_mean, pinf_mean, pl2_int;
  ana.ComputeCellError(
    *dg, p, tend, pnorm, pl2_err, pinf_err, pl2_mean, pinf_mean, pl2_int);

  if (getRank == 0) {
    printf(
      "nx=%3d (mean) L2(p)=%9.6g  Inf(p)=%9.6g\n", nx, pl2_mean, pinf_mean);
    printf("      (total) L2(p)=%9.6g  Inf(p)=%9.6g\n", pl2_err, pinf_err);
    printf("   (integral) L2(p)=%9.6g\n", pl2_int);
    if (exact_solution_expected)
      CHECK(pl2_mean < 1e-10);
    else if (limiter == "none")
      CHECK(pl2_mean < 0.1 / nx);
    else
      CHECK(pl2_mean < 0.3 / nx);
  }
}


TEST(OPERATOR_ADVECTION_TRANSIENT_DG)
{
  double dT(0.1), T1(1.0);
  auto rk_order = Amanzi::Explicit_TI::tvd_3rd_order;
  exact_solution_expected = true;
  AdvectionTransient<AnalyticDG02b>("square", 4, 4, dT, T1, rk_order, false);
  AdvectionTransient<AnalyticDG02b>(
    "square", 4, 4, dT, T1, rk_order, false, "gauss points");

  exact_solution_expected = false;
  AdvectionTransient<AnalyticDG06b>("square", 4, 4, dT, T1, rk_order);
  AdvectionTransient<AnalyticDG06>(
    "square", 4, 4, dT, T1, rk_order, false, "primal");
  AdvectionTransient<AnalyticDG06>(
    "square", 4, 4, dT, T1, rk_order, false, "dual", "high order");
  AdvectionTransient<AnalyticDG06>("square",
                                   4,
                                   4,
                                   dT,
                                   T1,
                                   rk_order,
                                   false,
                                   "gauss points",
                                   "high order",
                                   "Barth-Jespersen dg");

  /*
  double dT(0.01), T1(1.0);
  auto rk_order = Amanzi::Explicit_TI::tvd_3rd_order;
  AdvectionTransient<AnalyticDG06>("square",  16, 16, dT,  T1, rk_order);
  AdvectionTransient<AnalyticDG06>("square",  32, 32, dT/2,T1, rk_order);
  AdvectionTransient<AnalyticDG06>("square",  64, 64, dT/4,T1, rk_order);
  AdvectionTransient<AnalyticDG06>("square", 128,128, dT/8,T1, rk_order);
  */

  /*
  double dT(0.01), T1(1.0);
  auto rk_order = Amanzi::Explicit_TI::tvd_3rd_order;
  AdvectionTransient<AnalyticDG06>("test/triangular8.exo",    8,0, dT,   T1,
  rk_order); AdvectionTransient<AnalyticDG06>("test/triangular16.exo",  16,0,
  dT/2, T1, rk_order); AdvectionTransient<AnalyticDG06>("test/triangular32.exo",
  32,0, dT/4, T1, rk_order);
  AdvectionTransient<AnalyticDG06>("test/triangular64.exo",  64,0, dT/8, T1,
  rk_order); AdvectionTransient<AnalyticDG06>("test/triangular128.exo",128,0,
  dT/16,T1  rk_order);
  */

  /*
  double dT(0.001), T1(1.0);
  auto rk_order = Amanzi::Explicit_TI::tvd_3rd_order;
  AdvectionTransient<AnalyticDG06>("test/median15x16.exo",   16,0, dT,  T1,
  rk_order); AdvectionTransient<AnalyticDG06>("test/median32x33.exo",   32,0,
  dT/2,T1, rk_order); AdvectionTransient<AnalyticDG06>("test/median63x64.exo",
  64,0, dT/4,T1, rk_order);
  AdvectionTransient<AnalyticDG06>("test/median127x128.exo",128,0, dT/8,T1,
  rk_order);
  */

  /*
  double dT(0.01), T1(1.0);
  auto rk_order = Amanzi::Explicit_TI::tvd_3rd_order;
  AdvectionTransient<AnalyticDG06>("test/mesh_poly20x20.exo",   20,0, dT,  T1,
  rk_order); AdvectionTransient<AnalyticDG06>("test/mesh_poly40x40.exo",   40,0,
  dT/2,T1, rk_order);
  AdvectionTransient<AnalyticDG06>("test/mesh_poly80x80.exo",   80,0, dT/4,T1,
  rk_order); AdvectionTransient<AnalyticDG06>("test/mesh_poly160x160.exo",160,0,
  dT/8,T1, rk_order);
  */

  /*
  double dT(0.001), T1(1.0);
  auto rk_order = Amanzi::Explicit_TI::tvd_3rd_order;
  AdvectionTransient<AnalyticDG07>("square", 20, 20, dT,  T1, rk_order, false,
  "primal", "level set"); AdvectionTransient<AnalyticDG07>("square", 40, 40,
  dT/2,T1, rk_order, false, "primal", "level set");
  AdvectionTransient<AnalyticDG07>("square", 80, 80, dT/4,T1, rk_order, false,
  "primal", "level set"); AdvectionTransient<AnalyticDG07>("square",160,160,
  dT/8,T1, rk_order, false, "primal", "level set");
  */

  /*
  double dT(0.001), T1(1.0);
  auto rk_order = Amanzi::Explicit_TI::tvd_3rd_order;
  AdvectionTransient<AnalyticDG07>("test/median15x16.exo",   16,0, dT,  T1,
  rk_order, false, "primal", "level set");
  AdvectionTransient<AnalyticDG07>("test/median32x33.exo",   32,0, dT/2,T1,
  rk_order, false, "primal", "level set");
  AdvectionTransient<AnalyticDG07>("test/median63x64.exo",   64,0, dT/4,T1,
  rk_order, false, "primal", "level set");
  AdvectionTransient<AnalyticDG07>("test/median127x128.exo",128,0, dT/8,T1,
  rk_order, false, "primal", "level set");
  */

  /*
  double dT(0.001), T1(0.8);
  auto rk_order = Amanzi::Explicit_TI::tvd_3rd_order;
  AdvectionTransient<AnalyticDG07b>("square", 128, 128, dT/6,T1, rk_order,
  false, "primal", "level set");
  AdvectionTransient<AnalyticDG07b>("test/median127x128.exo", 128,0, dT/6,T1,
  rk_order, false, "primal", "level set");
  */

  /*
  double dT(0.01), T1(6.2832);
  auto rk_order = Amanzi::Explicit_TI::tvd_3rd_order;
  AdvectionTransient<AnalyticDG08>("square", 160,160, dT/8,T1, rk_order, true,
  "dual", "high order", "none");
  AdvectionTransient<AnalyticDG08>("test/median127x128.exo", 128,0, dT/8,T1,
  rk_order, true, "dual", "high order", "none");
  */
}
