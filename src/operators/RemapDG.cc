/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The helper advection-based base class for various remap methods. It
  provides support of time integration and calculation of various static
  and dynamic geometric quantities. The actual time-step loop could be
  implemented differently by an application.
*/

#include "RemapDG.hh"

namespace Amanzi {

/* *****************************************************************
* Initialization of remap: operarot and face velocity.
***************************************************************** */
RemapDG::RemapDG(
    const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
    const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
    Teuchos::ParameterList& plist) 
  : mesh0_(mesh0),
    mesh1_(mesh1),
    plist_(plist),
    dim_(mesh0->space_dimension())
{
  // mesh data
  ncells_owned_ = mesh0_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost_ = mesh0_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  nfaces_owned_ = mesh0_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost_ = mesh0_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  order_ = plist_.sublist("PK operator")
                 .sublist("flux operator").template get<int>("method order");

  // other control variable
  bc_type_ = Operators::OPERATOR_BC_NONE;
  std::string name = plist_.sublist("PK operator").template get<std::string>("boundary conditions");
  if (name == "remove")
    bc_type_ = Operators::OPERATOR_BC_REMOVE;

  consistent_jac_ = plist_.sublist("PK operator").template get<bool>("consistent jacobian");
  smoothness_ = plist_.sublist("limiter").template get<std::string>("smoothness indicator", "none");
  is_limiter_ = (plist_.sublist("limiter").template get<std::string>("limiter") != "none");

  // miscallateous
  nfun_ = 0;
  sharp_ = 0.0;

  t_adv_ = -1.0;
  t_flux_ = -1.0;
  t_reac_inv_ = -1.0;
}


/* *****************************************************************
* Initialization of opertors
***************************************************************** */
void RemapDG::InitializeOperators(const Teuchos::RCP<WhetStone::DG_Modal> dg)
{
  dg_ = dg;

  // create right-hand side operator
  // -- flux
  auto oplist = plist_.sublist("PK operator").sublist("flux operator");
  op_flux_ = Teuchos::rcp(new Operators::PDE_AdvectionRiemann(oplist, mesh0_));
  auto global_op = op_flux_->global_operator();

  // -- advection
  oplist = plist_.sublist("PK operator").sublist("advection operator");
  op_adv_ = Teuchos::rcp(new Operators::PDE_Abstract(oplist, global_op));

  // create left-hand side operator 
  oplist = plist_.sublist("PK operator").sublist("reaction operator");
  op_reac_ = Teuchos::rcp(new Operators::PDE_Reaction(oplist, mesh0_));

  // boundary data
  int nk = WhetStone::PolynomialSpaceDimension(dim_, order_);
  auto bc = Teuchos::rcp(new Operators::BCs(mesh0_, AmanziMesh::FACE, Operators::DOF_Type::VECTOR));
  std::vector<int>& bc_model = bc->bc_model();
  std::vector<std::vector<double> >& bc_value = bc->bc_value_vector(nk);

  const auto& fmap = mesh0_->face_map(true);
  const auto& bmap = mesh0_->exterior_face_map(true);
  for (int bf = 0; bf < bmap.NumMyElements(); ++bf) {
    int f = fmap.LID(bmap.GID(bf));
    for (int i = 0; i < nk; ++i) bc_value[f][i] = 0.0;
    bc_model[f] = bc_type_;
  }
  op_flux_->SetBCs(bc, bc);

  // memory allocation for velocities
  velf_ = Teuchos::rcp(new std::vector<WhetStone::Polynomial>(nfaces_wghost_));
  velc_ = Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_owned_));
  jac_ = Teuchos::rcp(new std::vector<WhetStone::VectorPolynomial>(ncells_owned_));
}


/* *****************************************************************
* Initialization of velocities and deformation tensor
***************************************************************** */
void RemapDG::InitializeFaceVelocity()
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
void RemapDG::InitializeJacobianMatrix()
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
* Initialization of the consistent jacobian
***************************************************************** */
void RemapDG::InitializeConsistentJacobianDeterminant()
{
  // constant part of determinant
  DynamicFaceVelocity(0.0);
  DynamicCellVelocity(0.0, false);

  op_adv_->SetupPolyVector(velc_);
  op_adv_->UpdateMatrices();

  op_reac_->Setup(jac_);
  op_reac_->UpdateMatrices(Teuchos::null);

  auto& matrices = op_reac_->local_matrices()->matrices;
  for (int n = 0; n < matrices.size(); ++n) matrices[n].Inverse();

  op_flux_->Setup(velc_, velf_);
  op_flux_->UpdateMatrices(velf_.ptr());
  op_flux_->ApplyBCs(true, true, true);

  CompositeVector& tmp = *op_reac_->global_operator()->rhs();
  CompositeVector one(tmp), u0(tmp), u1(tmp);
  Epetra_MultiVector& one_c = *one.ViewComponent("cell", true);

  one.PutScalarMasterAndGhosted(0.0);
  for (int c = 0; c < ncells_wghost_; ++c) one_c[0][c] = 1.0;

  op_flux_->global_operator()->Apply(one, tmp);
  op_reac_->global_operator()->Apply(tmp, u0);

  // linear part of determinant
  double dt(0.01);
  DynamicFaceVelocity(dt);
  DynamicCellVelocity(dt, false);
 
  op_adv_->SetupPolyVector(velc_);
  op_adv_->UpdateMatrices();

  op_flux_->Setup(velc_, velf_);
  op_flux_->UpdateMatrices(velf_.ptr());
  op_flux_->ApplyBCs(true, true, true);

  op_flux_->global_operator()->Apply(one, tmp);
  op_reac_->global_operator()->Apply(tmp, u1);
  u1.Update(-1.0/dt, u0, 1.0/dt);

  // save as polynomials
  int nk = one_c.NumVectors();
  Amanzi::WhetStone::DenseVector data(nk);
  Epetra_MultiVector& u0c = *u0.ViewComponent("cell", true);
  Epetra_MultiVector& u1c = *u1.ViewComponent("cell", true);

  jac0_.resize(ncells_owned_);
  jac1_.resize(ncells_owned_);

  for (int c = 0; c < ncells_owned_; ++c) {
    const auto& basis = dg_->cell_basis(c);

    for (int i = 0; i < nk; ++i) data(i) = u0c[i][c];
    jac0_[c] = basis.CalculatePolynomial(mesh0_, c, order_, data);

    for (int i = 0; i < nk; ++i) data(i) = u1c[i][c];
    jac1_[c] = basis.CalculatePolynomial(mesh0_, c, order_, data);
  }

  // t_adv_ = dt;
  // t_flux_ = dt;
  // t_reac_inv_ = 0.0;
}


/* *****************************************************************
* Main routine: evaluation of functional
***************************************************************** */
void RemapDG::FunctionalTimeDerivative(
    double t, const CompositeVector& u, CompositeVector& f)
{
  CompositeVector uu(u);
  if (is_limiter_) ApplyLimiter(t, uu);

  DynamicFaceVelocity(t);
  DynamicCellVelocity(t, consistent_jac_);

  // -- populate operators
  if (t_adv_ != t) { 
    op_adv_->SetupPolyVector(velc_);
    op_adv_->UpdateMatrices();
    // t_adv_ = t;
  }

  if (t_flux_ != t) {
    op_flux_->Setup(velc_, velf_);
    op_flux_->UpdateMatrices(velf_.ptr());
    op_flux_->ApplyBCs(true, true, true);
    // t_flux_ = t;
  }

  if (t_reac_inv_ != t) {
    op_reac_->Setup(jac_);
    op_reac_->UpdateMatrices(Teuchos::null);

    auto& matrices = op_reac_->local_matrices()->matrices;
    for (int n = 0; n < matrices.size(); ++n) matrices[n].Inverse();
    // t_reac_inv_ = t;
  }

  // -- solve the problem with mass matrix
  auto& x = *op_reac_->global_operator()->rhs(); // re-use memory
  op_reac_->global_operator()->Apply(uu, x);

  // -- calculate right-hand_side
  op_flux_->global_operator()->Apply(x, f);

  nfun_++;
}


/* *****************************************************************
* Calculates various geometric quantaties on intermediate meshes.
***************************************************************** */
void RemapDG::DynamicJacobianMatrix(
    int c, double t, const WhetStone::MatrixPolynomial& J, WhetStone::MatrixPolynomial& Jt)
{
  int nJ = J.NumRows();
  Jt = J * t;

  for (int i = 0; i < nJ; ++i) {
    Jt(i, i % dim_)(0) += 1.0;
  }
}


/* *****************************************************************
* Calculate face co-velocity in reference coordinates
***************************************************************** */
void RemapDG::DynamicFaceVelocity(double t)
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
void RemapDG::DynamicCellVelocity(double t, bool consistent_det)
{
  WhetStone::MatrixPolynomial Jt, C;
  for (int c = 0; c < ncells_owned_; ++c) {
    DynamicJacobianMatrix(c, t, J_[c], Jt);
    maps_->Cofactors(Jt, C);
    if (consistent_det) {
      double tmp = t * t / 2;
      (*jac_)[c][0] = t * jac0_[c] + tmp * jac1_[c];
      (*jac_)[c][0](0) += 1.0;
    } else {
      maps_->Determinant(Jt, (*jac_)[c]);
    }
    
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
* Limit conservative field u
***************************************************************** */
void RemapDG::ApplyLimiter(double t, CompositeVector& u)
{
  // calculate non-conservative field x
  auto& x = *op_reac_->global_operator()->rhs(); // re-use memory
  ConservativeToNonConservative(t, u, x);

  const Epetra_MultiVector& x_c = *x.ViewComponent("cell", true);
  int nk = x_c.NumVectors();

  // create list of cells where to apply limiter
  double threshold = -4.0;
  AmanziMesh::Entity_ID_List ids;

  for (int c = 0; c < ncells_owned_; ++c) {
    if (smoothness_ == "high order term" && order_ > 1) {
      double honorm(0.0);
      for (int i = dim_ + 1; i < nk; ++i)
        honorm += x_c[i][c] * x_c[i][c];

      double xnorm = honorm;
      for (int i = 0; i <= dim_; ++i)
        xnorm += x_c[i][c] * x_c[i][c];

      if (xnorm > 0.0 && std::log10(honorm / xnorm) > threshold)
        ids.push_back(c);
    } else {
      ids.push_back(c);
    }
  }

  int nids, itmp = ids.size();
  mesh0_->get_comm()->SumAll(&itmp, &nids, 1);
  sharp_ = std::max(sharp_, 100.0 * nids / x.ViewComponent("cell")->GlobalLength());

  // initialize limiter
  auto limlist = plist_.sublist("limiter");
  Operators::LimiterCell limiter(mesh0_);
  limiter.Init(limlist);

  std::string name = limlist.template get<std::string>("limiter");
  
  // apply limiter
  std::vector<int> bc_model(nfaces_wghost_, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_wghost_, 0.0);

  x.ScatterMasterToGhosted("cell");

  if (name == "Barth-Jespersen dg") {
    limiter.ApplyLimiter(ids, x.ViewComponent("cell", true), *dg_, bc_model, bc_value);
  } else {
    // -- create gradient in the natural basis
    WhetStone::DenseVector data(nk);

    CompositeVectorSpace cvs;
    cvs.SetMesh(mesh0_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, dim_);
    auto grad = Teuchos::rcp(new CompositeVector(cvs));
    Epetra_MultiVector& grad_c = *grad->ViewComponent("cell", true);

    // -- mean value is preserved automatiacally for the partially orthogonalized basis
    //    otherwise, a more complicated algorithm is needed
    AMANZI_ASSERT(dg_->cell_basis(0).id() == WhetStone::TAYLOR_BASIS_NORMALIZED_ORTHO);

    for (int c = 0; c < ncells_wghost_; ++c) {
      for (int i = 0; i < nk; ++i) data(i) = x_c[i][c];
      for (int i = dim_ + 1; i < nk; ++i) data(i) = 0.0;

      dg_->cell_basis(c).ChangeBasisMyToNatural(data);
      for (int i = 0; i < dim_; ++i) grad_c[i][c] = data(i + 1);
      x_c[0][c] = data(0);
    }

    // -- limit gradient and save it to solution
    limiter.ApplyLimiter(ids, x.ViewComponent("cell", true), 0, grad, bc_model, bc_value);

    for (int n = 0; n < ids.size(); ++n) {
      int c = ids[n];
      data(0) = x_c[0][c];
      for (int i = 0; i < dim_; ++i) data(i + 1) = grad_c[i][c];
      for (int i = dim_ + 1; i < nk; ++i) data(i) = 0.0;

      dg_->cell_basis(c).ChangeBasisNaturalToMy(data);
      for (int i = 0; i < nk; ++i) x_c[i][c] = data(i);
    }
  }

  NonConservativeToConservative(t, x, u);
}


/* *****************************************************************
* Change between conservative and non-conservative variables.
***************************************************************** */
void RemapDG::NonConservativeToConservative(
    double t, const CompositeVector& u, CompositeVector& v)
{
  DynamicFaceVelocity(t);
  DynamicCellVelocity(t, consistent_jac_);

  op_reac_->Setup(jac_);
  op_reac_->UpdateMatrices(Teuchos::null);
  op_reac_->global_operator()->Apply(u, v);
  t_reac_inv_ = -1.0;
}

void RemapDG::ConservativeToNonConservative(
    double t, const CompositeVector& u, CompositeVector& v)
{
  DynamicFaceVelocity(t);
  DynamicCellVelocity(t, consistent_jac_);

  if (t_reac_inv_ != t) {
    op_reac_->Setup(jac_);
    op_reac_->UpdateMatrices(Teuchos::null);

    auto& matrices = op_reac_->local_matrices()->matrices;
    for (int n = 0; n < matrices.size(); ++n) matrices[n].Inverse();
    // t_reac_inv_ = t;
  }

  op_reac_->global_operator()->Apply(u, v);
}

} // namespace Amanzi

