/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"

#include "GMVMesh.hh"
#include "Mesh.hh"
#include "PK_DomainFunctionFactory.hh"
#include "primary_variable_field_evaluator.hh"
#include "State.hh"
#include "WhetStoneDefs.hh"

#include "Energy_PK.hh"
#include "EnthalpyEvaluator.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Default constructor for Energy PK.
****************************************************************** */
Energy_PK::Energy_PK(Teuchos::ParameterList& pk_tree,
                     const Teuchos::RCP<Teuchos::ParameterList>& glist,
                     const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<TreeVector>& soln) :
    PK_PhysicalBDF(pk_tree, glist, S, soln),
    glist_(glist),
    passwd_("thermal")
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

  // keys
  temperature_key_ = Keys::getKey(domain_, "temperature"); 

  energy_key_ = Keys::getKey(domain_, "energy"); 
  prev_energy_key_ = Keys::getKey(domain_, "prev_energy");
  enthalpy_key_ = Keys::getKey(domain_, "enthalpy");
  conductivity_key_ = Keys::getKey(domain_, "thermal_conductivity");
  particle_density_key_ = Keys::getKey(domain_, "particle_density");

  ie_liquid_key_ = Keys::getKey(domain_, "internal_energy_liquid");
  ie_rock_key_ = Keys::getKey(domain_, "internal_energy_rock");

  mol_density_liquid_key_ = Keys::getKey(domain_, "molar_density_liquid");
  mass_density_liquid_key_ = Keys::getKey(domain_, "mass_density_liquid");

  darcy_flux_key_ = Keys::getKey(domain_, "darcy_flux"); 
}


/* ******************************************************************
* Construction of PK global variables.
****************************************************************** */
void Energy_PK::Setup(const Teuchos::Ptr<State>& S)
{
  // require first-requested state variables
  if (!S->HasField("atmospheric_pressure")) {
    S->RequireScalar("atmospheric_pressure", passwd_);
  }

  // require primary state variables
  if (!S->HasField(temperature_key_)) {
    std::vector<std::string> names({"cell", "face"});
    std::vector<int> ndofs(2, 1);
    std::vector<AmanziMesh::Entity_kind> locations({AmanziMesh::CELL, AmanziMesh::FACE});
 
    S->RequireField(temperature_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", temperature_key_);
    temperature_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator(temperature_key_, temperature_eval_);
  } else {
    temperature_eval_ = 
        Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(S->GetFieldEvaluator(temperature_key_));
  }

  // conserved quantity from the last time step.
  if (!S->HasField(prev_energy_key_)) {
    S->RequireField(prev_energy_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->GetField(prev_energy_key_, passwd_)->set_io_vis(false);
  }

  // other fields
  // -- energies
  S->RequireField(ie_liquid_key_)->SetMesh(mesh_)->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S->RequireFieldEvaluator(ie_liquid_key_);

  // -- densities
  S->RequireField(mol_density_liquid_key_)->SetMesh(mesh_)->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S->RequireFieldEvaluator(mol_density_liquid_key_);

  S->RequireField(mass_density_liquid_key_)->SetMesh(mesh_)->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S->RequireFieldEvaluator(mass_density_liquid_key_);

  // -- darcy flux
  if (!S->HasField(darcy_flux_key_)) {
    S->RequireField(darcy_flux_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }
}


/* ******************************************************************
* Basic initialization of energy classes.
****************************************************************** */
void Energy_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  // Create BCs objects
  // -- memory
  op_bc_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op_bc_enth_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  auto bc_list = Teuchos::rcp(new Teuchos::ParameterList(ep_list_->sublist("boundary conditions", false)));

  // -- temperature
  if (bc_list->isSublist("temperature")) {
    PK_DomainFunctionFactory<PK_DomainFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("temperature");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);
        bc_temperature_.push_back(bc_factory.Create(
            spec, "boundary temperature", AmanziMesh::FACE, Teuchos::null));
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
        bc_flux_.push_back(bc_factory.Create(
            spec, "outward energy flux", AmanziMesh::FACE, Teuchos::null));
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
  prec_include_enthalpy_ = ep_list_->sublist("operators")
                                   .get<bool>("include enthalpy in preconditioner", true);
}


/* ****************************************************************
* This completes initialization of missed fields in the state.
* This is useful for unit tests.
**************************************************************** */
void Energy_PK::InitializeFields_()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  if (!S_->GetField(temperature_key_, passwd_)->initialized()) {
    S_->GetFieldData(temperature_key_, passwd_)->PutScalar(298.0);
    S_->GetField(temperature_key_, passwd_)->set_initialized();

    if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
        *vo_->os() << "initialized temperature to default value 298 K." << std::endl;  
  }

  if (S_->GetField(darcy_flux_key_)->owner() == passwd_) {
    if (!S_->GetField(darcy_flux_key_, passwd_)->initialized()) {
      S_->GetFieldData(darcy_flux_key_, passwd_)->PutScalar(0.0);
      S_->GetField(darcy_flux_key_, passwd_)->set_initialized();

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
          *vo_->os() << "initialized darcy_flux to default value 0.0" << std::endl;  
    }
  }
}


/* ******************************************************************
* Converts scalar conductivity to a tensorial field: not used yet.
****************************************************************** */
bool Energy_PK::UpdateConductivityData(const Teuchos::Ptr<State>& S)
{
  bool update = S->GetFieldEvaluator(conductivity_key_)->HasFieldChanged(S, passwd_);
  if (update) {
    const auto& conductivity = *S->GetFieldData(conductivity_key_)->ViewComponent("cell");
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
void Energy_PK::UpdateSourceBoundaryData(double t_old, double t_new, const CompositeVector& u)
{
  for (int i = 0; i < bc_temperature_.size(); ++i) {
    bc_temperature_[i]->Compute(t_old, t_new);
  }

  for (int i = 0; i < bc_flux_.size(); ++i) {
    bc_flux_[i]->Compute(t_old, t_new);
  }

  for (int i = 0; i < srcs_.size(); ++i) {
    srcs_[i]->Compute(t_old, t_new);
  }

  ComputeBCs(u);
}


/* ******************************************************************
* Add source and sink terms.                                   
****************************************************************** */
void Energy_PK::AddSourceTerms(CompositeVector& rhs)
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
void Energy_PK::ComputeBCs(const CompositeVector& u)
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
  const auto& temp = *S_->GetFieldData(temperature_key_)->ViewComponent("face", true);
  const auto& n_l = *S_->GetFieldData(mol_density_liquid_key_)->ViewComponent("boundary_face", true);

  const Epetra_Map& ext_face_map = mesh_->exterior_face_map(true);
  const Epetra_Map& face_map = mesh_->face_map(true);

  int nbfaces = n_l.MyLength();
  for (int bf = 0; bf < nbfaces; ++bf) {
    int f = face_map.LID(ext_face_map.GID(bf));
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) temp[0][f] = bc_value[f];
  }

  // -- populate BCs
  S_->GetFieldEvaluator(enthalpy_key_)->HasFieldChanged(S_.ptr(), passwd_);
  const auto& enth = *S_->GetFieldData(enthalpy_key_)->ViewComponent("boundary_face", true);

  std::vector<int>& bc_model_enth_ = op_bc_enth_->bc_model();
  std::vector<double>& bc_value_enth_ = op_bc_enth_->bc_value();

  for (int n = 0; n < bc_model.size(); ++n) {
    bc_model_enth_[n] = Operators::OPERATOR_BC_NONE;
    bc_value_enth_[n] = 0.0;
  }

  for (int bf = 0; bf < nbfaces; ++bf) {
    int f = face_map.LID(ext_face_map.GID(bf));
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
      bc_model_enth_[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value_enth_[f] = enth[0][bf] * n_l[0][bf];
    }
  }
}


/* ******************************************************************
* Return a pointer to a local operator
****************************************************************** */
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
    Energy_PK::ModifyCorrection(
        double dt, Teuchos::RCP<const TreeVector> res,
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
  u->Data()->Comm()->SumAll(&tmp, &ntemp_clipped, 1);  // find the global clipping

  return (ntemp_clipped) > 0 ? AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED :
      AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}


/* ******************************************************************
* Return a pointer to a local operator
****************************************************************** */
Teuchos::RCP<Operators::Operator> Energy_PK::my_operator(
    const Operators::OperatorType& type)
{
  if (type == Operators::OPERATOR_MATRIX) return op_matrix_;
  else if (type == Operators::OPERATOR_PRECONDITIONER_RAW) return op_preconditioner_;
  return Teuchos::null;
}

}  // namespace Energy
}  // namespace Amanzi

