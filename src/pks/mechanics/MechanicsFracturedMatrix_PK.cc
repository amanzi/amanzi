/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Mechanics PK

*/

#include <vector>

// Amanzi
#include "UniqueLocalIndex.hh"

// Amanzi::Mechanics
#include "MechanicsFracturedMatrix_PK.hh"

namespace Amanzi {
namespace Mechanics {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* New constructor: extracts lists and requires fields.
****************************************************************** */
MechanicsFracturedMatrix_PK::MechanicsFracturedMatrix_PK(
  Teuchos::ParameterList& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& glist,
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<TreeVector>& soln)
  : MechanicsSmallStrain_PK(pk_tree, glist, S, soln)
{
  S_ = S;

  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  ec_list_ = Teuchos::sublist(pk_list, pk_name, true);

  // We also need iscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_solver_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(ec_list_, "time integrator", true);

  // domain and primary evaluators
  domain_ = ec_list_->get<std::string>("domain name", "domain");
  displacement_key_ = Keys::getKey(domain_, "displacement");
  AddDefaultPrimaryEvaluator(S_, displacement_key_);

  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = ec_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("MechanicsMatrixFracture", vlist));
}


/* ******************************************************************
* Contribution of this PK for DAG.
****************************************************************** */
void
MechanicsFracturedMatrix_PK::Setup()
{
  mesh_ = S_->GetMesh("domain");
  mesh_fracture_ = S_->GetMesh("fracture");
  int d = mesh_->getSpaceDimension();

  aperture_key_ = Keys::getKey("fracture", "aperture");
  ref_aperture_key_ = Keys::getKey("fracture", "ref_aperture");
  aperture_stiffness_key_ = Keys::getKey("fracture", "aperture_stiffness");
  pressure_key_ = Keys::getKey("fracture", "pressure");

  // displacement field
  std::vector<WhetStone::SchemaItem> items;
  items.push_back(std::make_tuple(AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::POINT, d));
  items.push_back(
    std::make_tuple(AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::NORMAL_COMPONENT, 1));

  std::vector<std::string> names = ec_list_->sublist("operators")
                                     .sublist("elasticity operator")
                                     .sublist("schema")
                                     .get<Teuchos::Array<std::string>>("fracture")
                                     .toVector();
  AMANZI_ASSERT(names.size() == 1);
  auto cvs = Operators::CreateFracturedMatrixCVS(mesh_, mesh_fracture_, names[0], items);

  if (!S_->HasRecord(displacement_key_)) {
    *S_->Require<CV_t, CVS_t>(displacement_key_, Tags::DEFAULT).SetMesh(mesh_)->SetGhosted(true) =
      *cvs;

    eval_ = Teuchos::rcp_static_cast<EvaluatorPrimary<CV_t, CVS_t>>(
      S_->GetEvaluatorPtr(displacement_key_, Tags::DEFAULT));
  }

  // aperture fields
  if (!S_->HasRecord(aperture_key_)) {
    S_->Require<CV_t, CVS_t>(aperture_key_, Tags::DEFAULT)
      .SetMesh(mesh_fracture_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  if (!S_->HasRecord(ref_aperture_key_)) {
    S_->Require<CV_t, CVS_t>(ref_aperture_key_, Tags::DEFAULT)
      .SetMesh(mesh_fracture_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  if (!S_->HasRecord(aperture_stiffness_key_)) {
    S_->Require<CV_t, CVS_t>(aperture_stiffness_key_, Tags::DEFAULT)
      .SetMesh(mesh_fracture_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  if (!S_->HasRecord(pressure_key_)) {
    S_->Require<CV_t, CVS_t>(pressure_key_, Tags::DEFAULT)
      .SetMesh(mesh_fracture_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  MechanicsSmallStrain_PK::Setup();
}


/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
****************************************************************** */
void
MechanicsFracturedMatrix_PK::FunctionalResidual(double t_old,
                                                double t_new,
                                                Teuchos::RCP<TreeVector> u_old,
                                                Teuchos::RCP<TreeVector> u_new,
                                                Teuchos::RCP<TreeVector> f)
{
  UpdateSourceBoundaryData(t_old, t_new);

  S_->GetEvaluator(shear_modulus_key_).Update(*S_, passwd_);

  op_matrix_elas_->global_operator()->Init();

  // add external forces
  auto rhs = op_matrix_->rhs();
  if (use_gravity_) AddGravityTerm(*rhs);
  if (poroelasticity_) AddPressureGradient(*rhs);
  if (thermoelasticity_) AddTemperatureGradient(*rhs);

  // assemble local matrices for fractures
  op_matrix_elas_->UpdateMatrices();
  AddFractureMatrices_(*rhs);

  op_matrix_elas_->ApplyBCs(true, true, true);

  // compute negative residual, A u - f
  op_matrix_->ComputeNegativeResidual(*u_new->Data(), *f->Data());
}


/* *******************************************************************
* Updates after successful time steps
******************************************************************* */
void
MechanicsFracturedMatrix_PK::CommitStep(double t_old, double t_new, const Tag& tag)
{
  MechanicsSmallStrain_PK::Setup();

  // compute aperture as diference of twin face DoFs
  auto& a_c = *S_->GetW<CV_t>(aperture_key_, Tags::DEFAULT, "").ViewComponent("cell");
  const auto& u_f = *S_->Get<CV_t>(displacement_key_).ViewComponent("face", true);
  const auto& fmap = *S_->Get<CV_t>(displacement_key_).ComponentMap("face", true);

  auto kind = AmanziMesh::Entity_kind::CELL;
  int ncells = mesh_fracture_->getNumEntities(kind, AmanziMesh::Parallel_kind::OWNED);

  for (int c = 0; c < ncells; ++c) {
    int f = mesh_fracture_->getEntityParent(kind, c);
    const auto& cells = mesh_->getFaceCells(f);

    int c1 = cells[0];
    if (c1 > ncells_owned_) c1 = cells[1];
    const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c1);

    int pos = std::distance(faces.begin(), std::find(faces.begin(), faces.end(), f));
    int g = fmap.FirstPointInElement(f);
    int shift = Operators::UniqueIndexFaceToCells(*mesh_, f, c1);

    a_c[0][c] = (u_f[0][g + shift] - u_f[0][g + 1 - shift]) * dirs[pos];
  }
}


/* ******************************************************************
* Calculate f(u, du/dt) = a d(s(u))/dt + A*u - rhs.
****************************************************************** */
void
MechanicsFracturedMatrix_PK::AddFractureMatrices_(CompositeVector& rhs)
{
  const auto& E_c = *S_->Get<CV_t>(aperture_stiffness_key_).ViewComponent("cell");
  const auto& a0_c = *S_->Get<CV_t>(ref_aperture_key_).ViewComponent("cell");
  const auto& p_c = *S_->Get<CV_t>(pressure_key_).ViewComponent("cell");
  const auto& rhs_f = *rhs.ViewComponent("face");

  const auto& fmap = *rhs.ComponentMap("face", true);
  rhs.PutScalarGhosted(0.0);

  auto kind = AmanziMesh::Entity_kind::CELL;
  int ncells = mesh_fracture_->getNumEntities(kind, AmanziMesh::Parallel_kind::OWNED);

  int pos, np, i1, i2;
  for (int c = 0; c < ncells; ++c) {
    int f = mesh_fracture_->getEntityParent(kind, c);
    const auto& cells = mesh_->getFaceCells(f);

    int c1 = cells[0];
    if (c1 >= ncells_owned_) c1 = cells[1];
    const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c1);
    int nfaces = faces.size();

    auto& Acell = op_matrix_elas_->local_op()->matrices[c1];
    int nrows = Acell.NumRows();

    np = 0;
    pos = 0;
    for (int i = 0; i < nfaces; ++i) {
      if (f == faces[i]) {
        pos = i;
        i1 = np;
      }
      np += fmap.ElementSize(faces[i]);
    }
    i1 = nrows - np + i1;
    i2 = i1 + 1;

    double area = mesh_fracture_->getCellVolume(c);
    double force = area * E_c[0][c] / a0_c[0][c];
    Acell(i1, i1) += force;
    Acell(i2, i2) += force;

    Acell(i1, i2) -= force;
    Acell(i2, i1) -= force;

    // update RHS using fluid pressure in fracture
    int g = fmap.FirstPointInElement(f);
    int shift = Operators::UniqueIndexFaceToCells(*mesh_, f, c1);

    rhs_f[0][g + shift] += area * p_c[0][c] * dirs[pos];
    rhs_f[0][g + 1 - shift] -= area * p_c[0][c] * dirs[pos];
  }

  rhs.GatherGhostedToMaster("face", Add);
}

} // namespace Mechanics
} // namespace Amanzi
