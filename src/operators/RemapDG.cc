/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "Epetra_Vector.h"

#include "RemapDG.hh"

namespace Amanzi {
namespace Operators {

/* *****************************************************************
 * Initialization of remap: operarot and face velocity.
 ***************************************************************** */
RemapDG::RemapDG(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
                 const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
                 Teuchos::ParameterList& plist)
  : mesh0_(mesh0), mesh1_(mesh1), plist_(plist), dim_(mesh0->space_dimension())
{
  // mesh data
  ncells_owned_ =
    mesh0_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost_ =
    mesh0_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  nfaces_owned_ =
    mesh0_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost_ =
    mesh0_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  auto& pklist = plist_.sublist("PK operator");
  order_ = pklist.sublist("flux operator").template get<int>("method order");

  // other control variable
  bc_type_ = OPERATOR_BC_NONE;
  std::string name = pklist.template get<std::string>("boundary conditions");
  if (name == "remove") bc_type_ = OPERATOR_BC_REMOVE;

  name = pklist.template get<std::string>("jacobian determinant method");
  if (name == "VEM")
    det_method_ = OPERATOR_DETERMINANT_VEM;
  else if (name == "exact time integration")
    det_method_ = OPERATOR_DETERMINANT_EXACT_TI;
  else if (name == "monotone")
    det_method_ = OPERATOR_DETERMINANT_MONOTONE;

  // initialize limiter
  auto limlist = plist_.sublist("limiter");
  is_limiter_ = (limlist.template get<std::string>("limiter") != "none");

  if (is_limiter_) {
    smoothness_ =
      limlist.template get<std::string>("smoothness indicator", "none");
    limiter_ = Teuchos::rcp(new LimiterCell(mesh0_));
    limiter_->Init(limlist);
  }

  // miscallateous
  nfun_ = 0;
  sharp_ = 0.0;
}


/* *****************************************************************
 * Initialization of opertors
 ***************************************************************** */
void
RemapDG::InitializeOperators(const Teuchos::RCP<WhetStone::DG_Modal> dg)
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
  auto bc = Teuchos::rcp(
    new BCs(mesh0_, AmanziMesh::FACE, Operators::DOF_Type::VECTOR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<std::vector<double>>& bc_value = bc->bc_value_vector(nk);

  const auto& fmap = mesh0_->face_map(true);
  const auto& bmap = mesh0_->exterior_face_map(true);
  for (int bf = 0; bf < bmap.getNodeNumElements(); ++bf) {
    int f = fmap.getLocalElement(bmap.getGlobalElement(bf));
    for (int i = 0; i < nk; ++i) bc_value[f][i] = 0.0;
    bc_model[f] = bc_type_;
  }
  op_flux_->SetBCs(bc, bc);

  // memory allocation for velocities
  velf_ = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(nfaces_wghost_));
  velc_ =
    Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_owned_));
  det_ =
    Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_owned_));

  // memory allocation for non-conservative field
  field_ =
    Teuchos::rcp(new CompositeVector(*op_reac_->global_operator()->rhs()));
}


/* *****************************************************************
 * Initialization of velocities and deformation tensor
 ***************************************************************** */
void
RemapDG::InitializeFaceVelocity()
{
  auto map_list = plist_.sublist("maps");
  WhetStone::MeshMapsFactory maps_factory;
  maps_ = maps_factory.Create(map_list, mesh0_, mesh1_);

  velf_vec_.resize(nfaces_wghost_);
  for (int f = 0; f < nfaces_wghost_; ++f) {
    maps_->VelocityFace(f, velf_vec_[f]);
  }
}


/* *****************************************************************
 * Initialization of the deformation tensor
 ***************************************************************** */
void
RemapDG::InitializeJacobianMatrix()
{
  WhetStone::Entity_ID_List faces;
  J_.resize(ncells_owned_);
  uc_.resize(ncells_owned_);

  for (int c = 0; c < ncells_owned_; ++c) {
    mesh0_->cell_get_faces(c, &faces);

    std::vector<WhetStone::VectorPolynomial> vvf;
    for (int n = 0; n < faces.size(); ++n) {
      vvf.push_back(velf_vec_[faces[n]]);
    }

    maps_->VelocityCell(c, vvf, uc_[c]);
    maps_->Jacobian(uc_[c], J_[c]);
  }
}


/* *****************************************************************
 * Main routine: evaluation of functional at time t
 ***************************************************************** */
void
RemapDG::FunctionalTimeDerivative(double t, const CompositeVector& u,
                                  CompositeVector& f)
{
  // -- populate operators
  //    geometric data were updated during solution modification
  op_adv_->SetupPolyVector(velc_);
  op_adv_->UpdateMatrices();

  op_flux_->Setup(velc_, velf_);
  op_flux_->UpdateMatrices(velf_.ptr());
  op_flux_->ApplyBCs(true, true, true);

  // -- calculate right-hand_side
  op_flux_->global_operator()->apply(*field_, f);

  nfun_++;
}


/* *****************************************************************
 * Limiting the non-conservative field at time t
 ***************************************************************** */
void
RemapDG::ModifySolution(double t, CompositeVector& u)
{
  DynamicFaceVelocity(t);
  DynamicCellVelocity(t);

  // populate operators
  op_reac_->Setup(det_);
  op_reac_->UpdateMatrices(Teuchos::null);

  // solve the problem with the mass matrix
  auto& matrices = op_reac_->local_matrices()->matrices;
  for (int n = 0; n < matrices.size(); ++n) matrices[n].InverseSPD();
  op_reac_->global_operator()->apply(u, *field_);

  // limit non-conservative field and update the conservative field
  if (is_limiter_) {
    // -- save original field and limit it
    auto& field_c = *field_->ViewComponent("cell");
    auto orig_c = field_c;

    ApplyLimiter(t, *field_);

    // -- recover original mass matrices FIXME (lipnikov@lanl.gov)
    for (int n = 0; n < matrices.size(); ++n) matrices[n].InverseSPD();

    // -- shift mean values
    auto& climiter = *limiter_->limiter();
    auto& u_c = *u.ViewComponent("cell");
    int nk = u_c.getNumVectors();

    for (int c = 0; c < ncells_owned_; ++c) {
      double a = climiter[c];
      if (a < 1.0) {
        double mass(0.0);
        for (int i = 0; i < nk; ++i) {
          mass += matrices[c](i, 0) * orig_c[i][c];
        }

        field_c[0][c] = a * orig_c[0][c] + (1.0 - a) * mass / matrices[c](0, 0);
      }
    }

    // -- update conservative field
    op_reac_->global_operator()->apply(*field_, u);
  }
}


/* *****************************************************************
 * Calculates various geometric quantaties on intermediate meshes.
 ***************************************************************** */
void
RemapDG::DynamicJacobianMatrix(int c, double t,
                               const WhetStone::MatrixPolynomial& J,
                               WhetStone::MatrixPolynomial& Jt)
{
  int nJ = J.NumRows();
  Jt = J * t;

  for (int i = 0; i < nJ; ++i) { Jt(i, i % dim_)(0) += 1.0; }
}


/* *****************************************************************
 * Calculate face co-velocity in reference coordinates
 ***************************************************************** */
void
RemapDG::DynamicFaceVelocity(double t)
{
  WhetStone::VectorPolynomial tmp, fmap, cn;

  for (int f = 0; f < nfaces_wghost_; ++f) {
    // cn = j J^{-t} N dA
    fmap = t * velf_vec_[f];
    maps_->NansonFormula(f, fmap, cn);
    (*velf_)[f] = velf_vec_[f] * cn;
  }
}


/* *****************************************************************
 * Cell co-velocity in reference coordinates and Jacobian determinant
 ***************************************************************** */
void
RemapDG::DynamicCellVelocity(double t)
{
  WhetStone::MatrixPolynomial Jt, C;
  for (int c = 0; c < ncells_owned_; ++c) {
    DynamicJacobianMatrix(c, t, J_[c], Jt);
    maps_->Cofactors(Jt, C);
    maps_->Determinant(Jt, (*det_)[c]);

    // negative co-velocity, v = -C^t u
    int nC = C.NumRows();
    (*velc_)[c].resize(nC);

    int kC = nC / dim_;
    for (int n = 0; n < kC; ++n) {
      int m = n * dim_;
      for (int i = 0; i < dim_; ++i) {
        (*velc_)[c][m + i].Reshape(dim_, 0, true);
        (*velc_)[c][m + i].set_origin(uc_[c][0].origin());

        for (int k = 0; k < dim_; ++k) {
          (*velc_)[c][m + i] -= C(m + k, i) * uc_[c][m + k];
        }
      }
    }
  }
}


/* *****************************************************************
 * Limit non-conservative field x
 ***************************************************************** */
void
RemapDG::ApplyLimiter(double t, CompositeVector& x)
{
  const Epetra_MultiVector& x_c = *x.ViewComponent("cell", true);
  int nk = x_c.getNumVectors();

  // create list of cells where to apply limiter
  double L(-1.0);
  double threshold = -4.0 * std::log10((double)order_) - L;
  AmanziMesh::Entity_ID_List ids;

  for (int c = 0; c < ncells_owned_; ++c) {
    if (smoothness_ == "high order term" && order_ > 1) {
      double honorm(0.0);
      for (int i = dim_ + 1; i < nk; ++i) honorm += x_c[i][c] * x_c[i][c];

      double xnorm = honorm;
      for (int i = 0; i <= dim_; ++i) xnorm += x_c[i][c] * x_c[i][c];

      if (xnorm > 0.0 && std::log10(honorm / xnorm) > threshold)
        ids.push_back(c);
    } else {
      ids.push_back(c);
    }
  }

  int nids, itmp = ids.size();
  mesh0Teuchos::reduceAll(*_->get_comm(), Teuchos::REDUCE_SUM, 1, &itmp, &nids);
  sharp_ =
    std::max(sharp_, 100.0 * nids / x.ViewComponent("cell")->getGlobalLength());

  // apply limiter
  std::vector<int> bc_model(nfaces_wghost_, OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost_, 0.0);

  x.ScatterMasterToGhosted("cell");

  if (limiter_->type() == OPERATOR_LIMITER_BARTH_JESPERSEN_DG) {
    limiter_->ApplyLimiter(
      ids, x.ViewComponent("cell", true), *dg_, bc_model, bc_value);
  } else {
    // -- create gradient in the natural basis
    WhetStone::DenseVector data(nk);

    CompositeVectorSpace cvs;
    cvs.SetMesh(mesh0_)->SetGhosted(true)->AddComponent(
      "cell", AmanziMesh::CELL, dim_);
    auto grad = Teuchos::rcp(new CompositeVector(cvs));
    Epetra_MultiVector& grad_c = *grad->ViewComponent("cell", true);

    // -- mean value is preserved automatiacally for the partially
    // orthogonalized basis
    //    otherwise, a more complicated algorithm is needed
    AMANZI_ASSERT(dg_->cell_basis(0).id() ==
                  WhetStone::TAYLOR_BASIS_NORMALIZED_ORTHO);

    for (int c = 0; c < ncells_wghost_; ++c) {
      for (int i = 0; i < nk; ++i) data(i) = x_c[i][c];
      for (int i = dim_ + 1; i < nk; ++i) data(i) = 0.0;

      dg_->cell_basis(c).ChangeBasisMyToNatural(data);
      for (int i = 0; i < dim_; ++i) grad_c[i][c] = data(i + 1);
      x_c[0][c] = data(0);
    }

    // -- limit gradient and save it to solution
    limiter_->ApplyLimiter(
      ids, x.ViewComponent("cell", true), 0, grad, bc_model, bc_value);

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


/* *****************************************************************
 * Change between conservative and non-conservative variables.
 ***************************************************************** */
void
RemapDG::NonConservativeToConservative(double t, const CompositeVector& u,
                                       CompositeVector& v)
{
  DynamicFaceVelocity(t);
  DynamicCellVelocity(t);

  op_reac_->Setup(det_);
  op_reac_->UpdateMatrices(Teuchos::null);
  op_reac_->global_operator()->apply(u, v);
}

void
RemapDG::ConservativeToNonConservative(double t, const CompositeVector& u,
                                       CompositeVector& v)
{
  DynamicFaceVelocity(t);
  DynamicCellVelocity(t);

  op_reac_->Setup(det_);
  op_reac_->UpdateMatrices(Teuchos::null);

  auto& matrices = op_reac_->local_matrices()->matrices;
  for (int n = 0; n < matrices.size(); ++n) matrices[n].InverseSPD();

  op_reac_->global_operator()->apply(u, v);
}

} // namespace Operators
} // namespace Amanzi
