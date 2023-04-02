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

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"

#include "EvaluatorMultiplicativeReciprocal.hh"
#include "EvaluatorPrimary.hh"
#include "Mesh.hh"
#include "Mesh_Algorithms.hh"
#include "PK_DomainFunctionFactory.hh"
#include "State.hh"
#include "WhetStoneDefs.hh"

#include "Energy_PK.hh"
#include "EnthalpyEvaluator.hh"

namespace Amanzi {
namespace Energy {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* Default constructor for Energy PK.
****************************************************************** */
Energy_PK::Energy_PK(Teuchos::ParameterList& pk_tree,
                     const Teuchos::RCP<Teuchos::ParameterList>& glist,
                     const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<TreeVector>& soln)
  : PK_PhysicalBDF(pk_tree, glist, S, soln), glist_(glist), passwd_(""), flow_on_manifold_(false)
{
  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  ep_list_ = Teuchos::sublist(pk_list, pk_name, true);

  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  ti_list_ = Teuchos::sublist(ep_list_, "time integrator");

  // domain name
  domain_ = ep_list_->get<std::string>("domain name", "domain");

  // create verbosity object
  S_ = S;
  mesh_ = S->GetMesh(domain_);
  dim = mesh_->space_dimension();

  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // workflow can be affected by the list of models
  auto physical_models = Teuchos::sublist(ep_list_, "physical models and assumptions");
  flow_on_manifold_ = physical_models->get<bool>("flow and transport in fractures", false);
}


/* ******************************************************************
* Construction of PK global variables.
****************************************************************** */
void
Energy_PK::Setup()
{
  temperature_key_ = Keys::getKey(domain_, "temperature");

  energy_key_ = Keys::getKey(domain_, "energy");
  prev_energy_key_ = Keys::getKey(domain_, "prev_energy");
  enthalpy_key_ = Keys::getKey(domain_, "enthalpy");
  conductivity_key_ = Keys::getKey(domain_, "thermal_conductivity");
  particle_density_key_ = Keys::getKey(domain_, "particle_density");

  aperture_key_ = Keys::getKey(domain_, "aperture");
  prev_aperture_key_ = Keys::getKey(domain_, "prev_aperture");
  conductivity_eff_key_ = Keys::getKey(domain_, "thermal_conductivity_effective");
  conductivity_gen_key_ = (!flow_on_manifold_) ? conductivity_key_ : conductivity_eff_key_;

  ie_liquid_key_ = Keys::getKey(domain_, "internal_energy_liquid");
  ie_gas_key_ = Keys::getKey(domain_, "internal_energy_gas");
  ie_rock_key_ = Keys::getKey(domain_, "internal_energy_rock");

  mol_density_liquid_key_ = Keys::getKey(domain_, "molar_density_liquid");
  mass_density_liquid_key_ = Keys::getKey(domain_, "mass_density_liquid");

  mol_density_gas_key_ = Keys::getKey(domain_, "molar_density_gas");
  x_gas_key_ = Keys::getKey(domain_, "molar_fraction_gas");

  vol_flowrate_key_ = Keys::getKey(domain_, "volumetric_flow_rate");
  mol_flowrate_key_ = Keys::getKey(domain_, "molar_flow_rate");
  sat_liquid_key_ = Keys::getKey(domain_, "saturation_liquid");
  Key pressure_key = Keys::getKey(domain_, "pressure");

  // require first-requested state variables
  if (!S_->HasRecord("atmospheric_pressure")) {
    S_->Require<double>("atmospheric_pressure", Tags::DEFAULT, "state");
  }

  // require primary state variables
  std::vector<std::string> names({ "cell", "face" });
  std::vector<int> ndofs(2, 1);
  std::vector<AmanziMesh::Entity_kind> locations({ AmanziMesh::CELL, AmanziMesh::FACE });

  S_->Require<CV_t, CVS_t>(temperature_key_, Tags::DEFAULT)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponents(names, locations, ndofs);

  if (!S_->HasEvaluator(temperature_key_, Tags::DEFAULT)) {
    Teuchos::ParameterList elist(temperature_key_);
    elist.set<std::string>("evaluator name", temperature_key_);
    temperature_eval_ = Teuchos::rcp(new EvaluatorPrimary<CV_t, CVS_t>(elist));
    S_->SetEvaluator(temperature_key_, Tags::DEFAULT, temperature_eval_);
  } else {
    temperature_eval_ = Teuchos::rcp_static_cast<EvaluatorPrimary<CV_t, CVS_t>>(
      S_->GetEvaluatorPtr(temperature_key_, Tags::DEFAULT));
  }

  // conserved quantity from the last time step.
  if (!S_->HasRecord(prev_energy_key_)) {
    S_->Require<CV_t, CVS_t>(prev_energy_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->GetRecordW(prev_energy_key_, passwd_).set_io_vis(false);
  }

  // other fields
  // -- energies (liquid, rock)
  S_->Require<CV_t, CVS_t>(ie_liquid_key_, Tags::DEFAULT, ie_liquid_key_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S_->RequireEvaluator(ie_liquid_key_, Tags::DEFAULT);

  S_->RequireDerivative<CV_t, CVS_t>(
      ie_liquid_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT, ie_liquid_key_)
    .SetGhosted();

  S_->Require<CV_t, CVS_t>(ie_rock_key_, Tags::DEFAULT, ie_rock_key_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S_->RequireEvaluator(ie_rock_key_, Tags::DEFAULT);

  S_->RequireDerivative<CV_t, CVS_t>(
      ie_rock_key_, Tags::DEFAULT, temperature_key_, Tags::DEFAULT, ie_rock_key_)
    .SetGhosted();

  // -- densities
  S_->Require<CV_t, CVS_t>(mol_density_liquid_key_, Tags::DEFAULT, mol_density_liquid_key_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S_->RequireEvaluator(mol_density_liquid_key_, Tags::DEFAULT);

  if (S_->GetEvaluator(mol_density_liquid_key_)
        .IsDifferentiableWRT(*S_, temperature_key_, Tags::DEFAULT)) {
    S_->RequireDerivative<CV_t, CVS_t>(mol_density_liquid_key_,
                                       Tags::DEFAULT,
                                       temperature_key_,
                                       Tags::DEFAULT,
                                       mol_density_liquid_key_)
      .SetGhosted();
  }

  S_->Require<CV_t, CVS_t>(mass_density_liquid_key_, Tags::DEFAULT, mass_density_liquid_key_)
    .SetMesh(mesh_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S_->RequireEvaluator(mass_density_liquid_key_, Tags::DEFAULT);

  if (S_->GetEvaluator(mass_density_liquid_key_)
        .IsDifferentiableWRT(*S_, temperature_key_, Tags::DEFAULT)) {
    S_->RequireDerivative<CV_t, CVS_t>(mass_density_liquid_key_,
                                       Tags::DEFAULT,
                                       temperature_key_,
                                       Tags::DEFAULT,
                                       mass_density_liquid_key_)
      .SetGhosted();
  }

  // -- volumetric and molar flow rates
  if (!S_->HasRecord(vol_flowrate_key_)) {
    S_->Require<CV_t, CVS_t>(vol_flowrate_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }
  if (!S_->HasRecord(mol_flowrate_key_)) {
    auto cvs = S_->Require<CV_t, CVS_t>(vol_flowrate_key_, Tags::DEFAULT, passwd_);
    *S_->Require<CV_t, CVS_t>(mol_flowrate_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true) = cvs;
  }

  // -- effective fracture conductivity
  if (flow_on_manifold_) {
    S_->Require<CV_t, CVS_t>(conductivity_eff_key_, Tags::DEFAULT, conductivity_eff_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    S_->Require<CV_t, CVS_t>(aperture_key_, Tags::DEFAULT, aperture_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator(aperture_key_, Tags::DEFAULT);

    Teuchos::ParameterList elist(conductivity_eff_key_);
    std::vector<std::string> listm(
      { Keys::getVarName(aperture_key_), Keys::getVarName(conductivity_key_) });
    elist.set<std::string>("my key", conductivity_eff_key_)
      .set<Teuchos::Array<std::string>>("multiplicative dependencies", listm)
      .set<std::string>("tag", "");
    auto eval = Teuchos::rcp(new EvaluatorMultiplicativeReciprocal(elist));
    S_->SetEvaluator(conductivity_eff_key_, Tags::DEFAULT, eval);
  }

  // if flow is missing, we need more fields
  // -- saturation
  if (!S_->HasRecord(sat_liquid_key_)) {
    S_->Require<CV_t, CVS_t>(sat_liquid_key_, Tags::DEFAULT, sat_liquid_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    AddDefaultIndependentEvaluator(S_, sat_liquid_key_, Tags::DEFAULT, 1.0);
  }

  // -- pressure
  if (!S_->HasRecord(pressure_key)) {
    S_->Require<CV_t, CVS_t>(pressure_key, Tags::DEFAULT, pressure_key)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    AddDefaultIndependentEvaluator(S_, pressure_key, Tags::DEFAULT, 101325.0);
  }

  // -- fracture aperture
  if (flow_on_manifold_) {
    S_->Require<CV_t, CVS_t>(aperture_key_, Tags::DEFAULT, aperture_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator(aperture_key_, Tags::DEFAULT);

    S_->Require<CV_t, CVS_t>(prev_aperture_key_, Tags::DEFAULT)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
}


/* ******************************************************************
* Basic initialization of energy classes.
****************************************************************** */
void
Energy_PK::Initialize()
{
  // Create BCs objects
  // -- memory
  op_bc_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op_bc_enth_ =
    Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  auto bc_list =
    Teuchos::rcp(new Teuchos::ParameterList(ep_list_->sublist("boundary conditions", false)));

  // -- temperature
  if (bc_list->isSublist("temperature")) {
    PK_DomainFunctionFactory<PK_DomainFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("temperature");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);
        bc_temperature_.push_back(
          bc_factory.Create(spec, "boundary temperature", AmanziMesh::FACE, Teuchos::null));
      }
    }
  }

  // -- energy flux
  if (bc_list->isSublist("energy flux")) {
    PK_DomainFunctionFactory<PK_DomainFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("energy flux");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);
        bc_flux_.push_back(
          bc_factory.Create(spec, "outward energy flux", AmanziMesh::FACE, Teuchos::null));
      }
    }
  }


  if (ep_list_->isSublist("source terms")) {
    PK_DomainFunctionFactory<PK_DomainFunction> factory(mesh_, S_);
    auto src_list = ep_list_->sublist("source terms");
    for (auto it = src_list.begin(); it != src_list.end(); ++it) {
      std::string name = it->first;
      if (src_list.isSublist(name)) {
        Teuchos::ParameterList& spec = src_list.sublist(name);
        srcs_.push_back(factory.Create(spec, "source", AmanziMesh::CELL, Teuchos::null));
      }
    }
  }

  // initialized fields
  InitializeFields_();

  // other parameters
  prec_include_enthalpy_ =
    ep_list_->sublist("operators").get<bool>("include enthalpy in preconditioner", true);
}


/* ****************************************************************
* This completes initialization of missed fields in the state.
* This is useful for unit tests.
**************************************************************** */
void
Energy_PK::InitializeFields_()
{
  InitializeCVField(S_, *vo_, temperature_key_, Tags::DEFAULT, passwd_, 298.0);
  InitializeCVField(S_, *vo_, vol_flowrate_key_, Tags::DEFAULT, passwd_, 0.0);
  InitializeCVField(S_, *vo_, mol_flowrate_key_, Tags::DEFAULT, passwd_, 0.0);
}


/* ******************************************************************
* Converts scalar conductivity to a tensorial field: not used yet.
****************************************************************** */
bool
Energy_PK::UpdateConductivityData(const Teuchos::Ptr<State>& S)
{
  bool update = S->GetEvaluator(conductivity_gen_key_).Update(*S, passwd_);
  if (update) {
    const auto& conductivity = *S->Get<CV_t>(conductivity_gen_key_).ViewComponent("cell");
    WhetStone::Tensor Ktmp(dim, 1);

    K.clear();
    for (int c = 0; c < ncells_owned; c++) {
      Ktmp(0, 0) = conductivity[0][c];
      K.push_back(Ktmp);
    }
  }
  return update;
}


/* ******************************************************************
* A wrapper for updating boundary conditions.
****************************************************************** */
void
Energy_PK::UpdateSourceBoundaryData(double t_old, double t_new, const CompositeVector& u)
{
  for (int i = 0; i < bc_temperature_.size(); ++i) { bc_temperature_[i]->Compute(t_old, t_new); }

  for (int i = 0; i < bc_flux_.size(); ++i) { bc_flux_[i]->Compute(t_old, t_new); }

  for (int i = 0; i < srcs_.size(); ++i) { srcs_[i]->Compute(t_old, t_new); }

  ComputeBCs(u);
}


/* ******************************************************************
* Add source and sink terms.
****************************************************************** */
void
Energy_PK::AddSourceTerms(CompositeVector& rhs)
{
  Epetra_MultiVector& rhs_cell = *rhs.ViewComponent("cell");

  for (int i = 0; i < srcs_.size(); ++i) {
    for (auto it = srcs_[i]->begin(); it != srcs_[i]->end(); ++it) {
      int c = it->first;
      rhs_cell[0][c] += mesh_->cell_volume(c) * it->second[0];
    }
  }
}


/* ******************************************************************
* Add a boundary marker to used faces.
* WARNING: we can skip update of ghost boundary faces, b/c they
* should be always owned.
****************************************************************** */
void
Energy_PK::ComputeBCs(const CompositeVector& u)
{
  std::vector<int>& bc_model = op_bc_->bc_model();
  std::vector<double>& bc_value = op_bc_->bc_value();
  std::vector<double>& bc_mixed = op_bc_->bc_mixed();

  for (int n = 0; n < bc_model.size(); n++) {
    bc_model[n] = Operators::OPERATOR_BC_NONE;
    bc_value[n] = 0.0;
    bc_mixed[n] = 0.0;
  }

  for (int i = 0; i < bc_temperature_.size(); ++i) {
    for (auto it = bc_temperature_[i]->begin(); it != bc_temperature_[i]->end(); ++it) {
      int f = it->first;
      bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value[f] = it->second[0];
    }
  }

  for (int i = 0; i < bc_flux_.size(); ++i) {
    for (auto it = bc_flux_[i]->begin(); it != bc_flux_[i]->end(); ++it) {
      int f = it->first;
      bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
      bc_value[f] = it->second[0];
    }
  }

  dirichlet_bc_faces_ = 0;
  for (int f = 0; f < nfaces_owned; ++f) {
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) dirichlet_bc_faces_++;
  }
#ifdef HAVE_MPI
  int tmp = dirichlet_bc_faces_;
  mesh_->get_comm()->SumAll(&tmp, &dirichlet_bc_faces_, 1);
#endif

  // additional boundary conditions
  // -- copy essential conditions to primary variables
  // BoundaryDataToFaces(op_bc_, *S_->GetFieldData(temperature_key_, passwd_));

  // -- populate BCs
  S_->GetEvaluator(enthalpy_key_).Update(*S_, passwd_);
  const auto& enth = *S_->Get<CV_t>(enthalpy_key_).ViewComponent("boundary_face", true);

  std::vector<int>& bc_model_enth_ = op_bc_enth_->bc_model();
  std::vector<double>& bc_value_enth_ = op_bc_enth_->bc_value();

  for (int n = 0; n < bc_model.size(); ++n) {
    bc_model_enth_[n] = Operators::OPERATOR_BC_NONE;
    bc_value_enth_[n] = 0.0;
  }

  int nbfaces = enth.MyLength();
  for (int bf = 0; bf < nbfaces; ++bf) {
    int f = getBoundaryFaceFace(*mesh_, bf);
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
      bc_model_enth_[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value_enth_[f] = enth[0][bf];
    }
  }
}


/* ******************************************************************
* Return a pointer to a local operator
****************************************************************** */
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
Energy_PK::ModifyCorrection(double dt,
                            Teuchos::RCP<const TreeVector> res,
                            Teuchos::RCP<const TreeVector> u,
                            Teuchos::RCP<TreeVector> du)
{
  int ntemp_clipped(0);

  // clipping temperature
  for (auto comp = u->Data()->begin(); comp != u->Data()->end(); ++comp) {
    const auto& uc = *u->Data()->ViewComponent(*comp);
    auto& duc = *du->Data()->ViewComponent(*comp);

    int ncomp = u->Data()->size(*comp, false);
    for (int i = 0; i < ncomp; ++i) {
      double tmp0 = uc[0][i] / 20.0;
      double tmp1 = std::min(tmp0, uc[0][i] - 273.65);
      if (duc[0][i] < -tmp0) {
        ntemp_clipped++;
        duc[0][i] = -tmp0;
      } else if (duc[0][i] > tmp1) {
        ntemp_clipped++;
        duc[0][i] = tmp1;
      }
    }
  }

  int tmp(ntemp_clipped);
  u->Data()->Comm()->SumAll(&tmp, &ntemp_clipped, 1); // find the global clipping

  return (ntemp_clipped) > 0 ? AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED :
                               AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}


/* ******************************************************************
* Return a pointer to a local operator
****************************************************************** */
Teuchos::RCP<Operators::Operator>
Energy_PK::my_operator(const Operators::OperatorType& type)
{
  if (type == Operators::OPERATOR_MATRIX)
    return op_matrix_;
  else if (type == Operators::OPERATOR_PRECONDITIONER_RAW)
    return op_preconditioner_;
  return Teuchos::null;
}

} // namespace Energy
} // namespace Amanzi
