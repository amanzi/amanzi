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

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Default constructor for Energy PK.
****************************************************************** */
Energy_PK::Energy_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                     Teuchos::RCP<State> S) :
    glist_(glist),
    passwd_("thermal")
{
  S_ = S;
  mesh_ = S->GetMesh();
  dim = mesh_->space_dimension();

  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  energy_key_ = "energy";
  prev_energy_key_ = "prev_energy";
  enthalpy_key_ = "enthalpy";
  conductivity_key_ = "thermal_conductivity";
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
  std::vector<std::string> names(2);
  names[0] = "cell";
  names[1] = "face";
 
  std::vector<AmanziMesh::Entity_kind> locations(2);
  locations[0] = AmanziMesh::CELL;
  locations[1] = AmanziMesh::FACE;
 
  std::vector<int> ndofs(2, 1);
  
  temperature_key_ = Keys::getKey(domain_, "temperature"); 
  if (!S->HasField(temperature_key_)) {
    S->RequireField(temperature_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", temperature_key_);
    temperature_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator(temperature_key_, temperature_eval_);
  }

  // conserved quantity from the last time step.
  if (!S->HasField("prev_energy")) {
    S->RequireField("prev_energy", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->GetField("prev_energy", passwd_)->set_io_vis(false);
  }

  // Fields for energy as independent PK
  if (!S->HasField("darcy_flux")) {
    S->RequireField("darcy_flux", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }
}


/* ******************************************************************
* Basic initialization of energy classes.
****************************************************************** */
void Energy_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  // Energy list has only one sublist
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist_, "PKs", true);
  Teuchos::RCP<Teuchos::ParameterList> ep_list = Teuchos::sublist(pk_list, "energy", true);

  // Create BCs objects
  // -- memory
  op_bc_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op_bc_enth_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  Teuchos::RCP<Teuchos::ParameterList>
      bc_list = Teuchos::rcp(new Teuchos::ParameterList(ep_list->sublist("boundary conditions", true)));

  // -- temperature
  if (bc_list->isSublist("temperature")) {
    PK_DomainFunctionFactory<PK_DomainFunction> bc_factory(mesh_);

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
    PK_DomainFunctionFactory<PK_DomainFunction> bc_factory(mesh_);

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

  // initialized fields
  InitializeFields_();

  // other parameters
  prec_include_enthalpy_ = ep_list->sublist("operators")
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

  if (S_->GetField("darcy_flux")->owner() == passwd_) {
    if (!S_->GetField("darcy_flux", passwd_)->initialized()) {
      S_->GetFieldData("darcy_flux", passwd_)->PutScalar(0.0);
      S_->GetField("darcy_flux", passwd_)->set_initialized();

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
    const Epetra_MultiVector& conductivity = *S->GetFieldData(conductivity_key_)->ViewComponent("cell");
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

  ComputeBCs(u);
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
  int flag_essential_bc = (dirichlet_bc_faces_ > 0) ? 1 : 0;

  // verify that the algebraic problem is consistent
#ifdef HAVE_MPI
  int flag = flag_essential_bc;
  mesh_->get_comm()->MaxAll(&flag, &flag_essential_bc, 1);  // find the global maximum
#endif
  if (! flag_essential_bc && vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "WARNING: no essential boundary conditions, solver may fail" << std::endl;
  }

  // additional boundary conditions
  AmanziMesh::Entity_ID_List cells;
  const auto& enth = *S_->GetFieldData(enthalpy_key_)->ViewComponent("cell");

  std::vector<int>& bc_model_enth_ = op_bc_enth_->bc_model();
  std::vector<double>& bc_value_enth_ = op_bc_enth_->bc_value();

  for (int n = 0; n < bc_model.size(); ++n) {
    bc_model_enth_[n] = Operators::OPERATOR_BC_NONE;
    bc_value_enth_[n] = 0.0;
  }

  for (int f = 0; f < bc_model.size(); ++f) {
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

      bc_model_enth_[f] = Operators::OPERATOR_BC_DIRICHLET;
      bc_value_enth_[f] = enth[0][cells[0]];
    }
  }
}

}  // namespace Energy
}  // namespace Amanzi

