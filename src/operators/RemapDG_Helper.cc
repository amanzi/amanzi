/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  The helper advection-based base class for various remap methods. It
  provides support of time integration and calculation of various static
  and dynamic geometric quantities. The actual time-step loop could be
  implemented differently by an application.
*/

#include "Epetra_Vector.h"

#include "ReconstructionCellLinear.hh"
#include "RemapDG.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace Operators {

/* *****************************************************************
* Initialization of remap: operarot and face velocity.
***************************************************************** */
RemapDG_Helper::RemapDG_Helper(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
                               const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
                               Teuchos::ParameterList& plist)
  : mesh0_(mesh0), mesh1_(mesh1), dim_(mesh0->getSpaceDimension()), plist_(plist)
{
  // mesh data
  ncells_owned_ = mesh0_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  ncells_wghost_ = mesh0_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  nfaces_owned_ = mesh0_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  nfaces_wghost_ = mesh0_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  if (mesh0_->hasEdges()) {
    nedges_owned_ = mesh0_->getNumEntities(AmanziMesh::Entity_kind::EDGE, AmanziMesh::Parallel_kind::OWNED);
    nedges_wghost_ = mesh0_->getNumEntities(AmanziMesh::Entity_kind::EDGE, AmanziMesh::Parallel_kind::ALL);
  }

  auto& pklist = plist_.sublist("PK operator");
  order_ = pklist.sublist("flux operator").sublist("schema").template get<int>("method order");

  // other control variable
  bc_type_ = OPERATOR_BC_NONE;
  std::string name = pklist.template get<std::string>("boundary conditions");
  if (name == "remove") bc_type_ = OPERATOR_BC_REMOVE;

  // initialize limiter
  auto limlist = plist_.sublist("limiter");
  is_limiter_ = (limlist.template get<std::string>("limiter") != "none");

  if (is_limiter_) {
    smoothness_ = limlist.template get<std::string>("smoothness indicator", "none");
    limiter_ = Teuchos::rcp(new LimiterCellDG(mesh0_));
    limiter_->Init(limlist);
  }

  // miscallateous
  nfun_ = 0;
  sharp_ = 0.0;
}


/* *****************************************************************
* Initialization of operators
***************************************************************** */
void
RemapDG_Helper::InitializeOperators(const Teuchos::RCP<WhetStone::DG_Modal> dg)
{
  dg_ = dg;

  // create right-hand side operator
  // -- flux
  auto oplist = plist_.sublist("PK operator").sublist("flux operator");
  op_flux_ = Teuchos::rcp(new PDE_AdvectionRiemann(oplist, mesh0_));
  auto global_op = op_flux_->global_operator();

  // -- advection
  oplist = plist_.sublist("PK operator").sublist("advection operator");
  op_adv_ = Teuchos::rcp(new PDE_Abstract(oplist, global_op));

  // create left-hand side operator
  oplist = plist_.sublist("PK operator").sublist("reaction operator");
  op_reac_ = Teuchos::rcp(new PDE_Reaction(oplist, mesh0_));

  // boundary data
  int nk = WhetStone::PolynomialSpaceDimension(dim_, order_);
  auto bc = Teuchos::rcp(new BCs(mesh0_, AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::VECTOR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<std::vector<double>>& bc_value = bc->bc_value_vector(nk);

  const auto& fmap = mesh0_->getMap(AmanziMesh::Entity_kind::FACE,true);
  const auto& bmap = mesh0_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,true);
  for (int bf = 0; bf < bmap.NumMyElements(); ++bf) {
    int f = fmap.LID(bmap.GID(bf));
    for (int i = 0; i < nk; ++i) bc_value[f][i] = 0.0;
    bc_model[f] = bc_type_;
  }
  op_flux_->SetBCs(bc, bc);

  // memory allocation for velocities
  velf_ = Teuchos::rcp(new std::vector<WhetStone::SpaceTimePolynomial>(nfaces_wghost_));
  velc_ = Teuchos::rcp(new std::vector<WhetStone::VectorSpaceTimePolynomial>(ncells_owned_));
  det_ = Teuchos::rcp(new std::vector<WhetStone::SpaceTimePolynomial>(ncells_owned_));

  // memory allocation for non-conservative field
  field_ = Teuchos::rcp(new CompositeVector(*op_reac_->global_operator()->rhs()));

  // memory allocation for new features
  jac_ = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(ncells_owned_));
}


/* *****************************************************************
* Initialization of static edge and face velocities
***************************************************************** */
void
RemapDG_Helper::StaticEdgeFaceVelocities()
{
  auto map_list = plist_.sublist("maps");
  WhetStone::MeshMapsFactory maps_factory;
  maps_ = maps_factory.Create(map_list, mesh0_, mesh1_);

  velf_vec_.resize(nfaces_wghost_);
  for (int f = 0; f < nfaces_wghost_; ++f) { maps_->VelocityFace(f, velf_vec_[f]); }

  if (mesh0_->hasEdges()) {
    vele_vec_.resize(nedges_wghost_);
    for (int e = 0; e < nedges_wghost_; ++e) { maps_->VelocityEdge(e, vele_vec_[e]); }
  }
}


/* *****************************************************************
* Initialization of the constant cell velocity
***************************************************************** */
void
RemapDG_Helper::StaticCellVelocity()
{
  uc_.resize(ncells_owned_);

  for (int c = 0; c < ncells_owned_; ++c) {
    // faces are always included
    const auto& faces = mesh0_->getCellFaces(c);

    std::vector<WhetStone::VectorPolynomial> vve, vvf;
    for (int n = 0; n < faces.size(); ++n) { vvf.push_back(velf_vec_[faces[n]]); }

    // edges are included in 3D only
    if (dim_ == 3) {
      auto edges = mesh0_->getCellEdges(c);

      for (int n = 0; n < edges.size(); ++n) { vve.push_back(vele_vec_[edges[n]]); }
    }

    maps_->VelocityCell(c, vve, vvf, uc_[c]);
  }
}


/* *****************************************************************
* Initialization of space-time co-velocity v = u * (j J^{-t} N)
***************************************************************** */
void
RemapDG_Helper::StaticFaceCoVelocity()
{
  WhetStone::VectorSpaceTimePolynomial cn;
  for (int f = 0; f < nfaces_wghost_; ++f) {
    WhetStone::VectorSpaceTimePolynomial map(dim_, dim_, 1), tmp(dim_, dim_, 0);
    const auto& origin = velf_vec_[f][0].get_origin();

    for (int i = 0; i < dim_; ++i) {
      map[i][0].Reshape(dim_, std::max(1, order_), true);
      map[i][0](1, i) = 1.0; // map = x
      map[i][0].set_origin(origin);
      map[i][1] = velf_vec_[f][i]; // map = x + t * u

      tmp[i][0] = velf_vec_[f][i];
    }

    maps_->NansonFormula(f, map, cn);
    (*velf_)[f] = tmp * cn;
  }
}


/* *****************************************************************
* Initialization of the constant cell velocity
***************************************************************** */
void
RemapDG_Helper::StaticCellCoVelocity()
{
  for (int c = 0; c < ncells_owned_; ++c) {
    WhetStone::MatrixPolynomial Jc;
    maps_->Jacobian(uc_[c], Jc);

    // space-time cell velocity: v = -j J^{-1} u = -C^t u
    WhetStone::MatrixSpaceTimePolynomial Jt(dim_, dim_, dim_, 1), Ct;
    WhetStone::VectorSpaceTimePolynomial tmp(dim_, dim_, 0);
    const auto& origin = uc_[c][0].get_origin();

    for (int i = 0; i < dim_; ++i) {
      for (int j = 0; j < dim_; ++j) {
        Jt(i, j)[0].Reshape(dim_, 0, true);
        Jt(i, j)[0].set_origin(origin);
        Jt(i, j)[1] = Jc(i, j); // Jt = 1 + t * J
      }
      Jt(i, i)[0](0) = 1.0;
      tmp[i][0] = uc_[c][i];
    }

    maps_->Cofactors(Jt, Ct);

    tmp *= -1.0;
    Ct.Multiply(tmp, (*velc_)[c], true);

    maps_->Determinant(Jt, (*det_)[c]);
  }
}


/* *****************************************************************
* Limit non-conservative field x
***************************************************************** */
void
RemapDG_Helper::ApplyLimiter(double t, CompositeVector& x)
{
  auto& x_c = *x.ViewComponent("cell", true);
  int nk = x_c.NumVectors();

  // create list of cells where to apply limiter
  double L(-1.0);
  double threshold = -4.0 * std::log10((double)order_) - L;
  AmanziMesh::Entity_ID_View ids("ids", ncells_owned_);
  int ids_ct = 0; 

  for (int c = 0; c < ncells_owned_; ++c) {
    if (smoothness_ == "high order term" && order_ > 1) {
      double honorm(0.0);
      for (int i = dim_ + 1; i < nk; ++i) honorm += x_c[i][c] * x_c[i][c];

      double xnorm = honorm;
      for (int i = 0; i <= dim_; ++i) xnorm += x_c[i][c] * x_c[i][c];

      if (xnorm > 0.0 && std::log10(honorm / xnorm) > threshold) ids[ids_ct++] = c;
    } else {
      ids[ids_ct++] = c;
    }
  }
  Kokkos::resize(ids, ids_ct); 
  int nids, itmp = ids.size();
  mesh0_->getComm()->SumAll(&itmp, &nids, 1);
  sharp_ = std::max(sharp_, 100.0 * nids / x.ViewComponent("cell")->GlobalLength());

  // apply limiter
  std::vector<int> bc_model(nfaces_wghost_, OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost_, 0.0);

  x.ScatterMasterToGhosted("cell");

  if (limiter_->get_type() == OPERATOR_LIMITER_BARTH_JESPERSEN_DG ||
      limiter_->get_type() == OPERATOR_LIMITER_MICHALAK_GOOCH_DG ||
      limiter_->get_type() == OPERATOR_LIMITER_BARTH_JESPERSEN_DG_HIERARCHICAL) {
    limiter_->ApplyLimiterDG(ids, x.ViewComponent("cell", true), *dg_, bc_model, bc_value);
  } else {
    // -- create gradient in the natural basis
    WhetStone::DenseVector data(nk);

    CompositeVectorSpace cvs;
    cvs.SetMesh(mesh0_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::Entity_kind::CELL, dim_);
    auto grad = Teuchos::rcp(new CompositeVector(cvs));
    Epetra_MultiVector& grad_c = *grad->ViewComponent("cell", true);

    // -- mean value is preserved automatically for the partially orthogonalized basis
    //    otherwise, a more complicated algorithm is needed
    AMANZI_ASSERT(nk > dim_ || dg_->cell_basis(0).id() == WhetStone::TAYLOR_BASIS_NORMALIZED_ORTHO);

    for (int c = 0; c < ncells_wghost_; ++c) {
      for (int i = 0; i < nk; ++i) data(i) = x_c[i][c];
      for (int i = dim_ + 1; i < nk; ++i) data(i) = 0.0;

      dg_->cell_basis(c).ChangeBasisMyToNatural(data);
      for (int i = 0; i < dim_; ++i) grad_c[i][c] = data(i + 1);
      x_c[0][c] = data(0);
    }

    // -- limit gradient and save it to solution
    //    Reconstruction object does nothing but keeping poiter to gradient
    auto lifting = Teuchos::rcp(new ReconstructionCellLinear(mesh0_, grad));
    limiter_->ApplyLimiter(ids, x.ViewComponent("cell", true), 0, lifting, bc_model, bc_value);

    for (int n = 0; n < ids.size(); ++n) {
      int c = ids[n];
      data(0) = x_c[0][c];
      for (int i = 0; i < dim_; ++i) data(i + 1) = grad_c[i][c];
      for (int i = dim_ + 1; i < nk; ++i) data(i) = 0.0;

      dg_->cell_basis(c).ChangeBasisNaturalToMy(data);
      for (int i = 0; i < nk; ++i) x_c[i][c] = data(i);
    }
  }
}

} // namespace Operators
} // namespace Amanzi
