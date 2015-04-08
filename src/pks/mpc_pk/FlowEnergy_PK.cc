/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov

  Process kernel for coupling Flow PK with Energy PK.
*/

#include "FlowEnergy_PK.hh"
#include "MPCStrong.hh"

namespace Amanzi {

/* ******************************************************************* 
* Constructor
******************************************************************* */
FlowEnergy_PK::FlowEnergy_PK(Teuchos::ParameterList& pk_tree,
                             const Teuchos::RCP<Teuchos::ParameterList>& glist,
                             const Teuchos::RCP<State>& S,
                             const Teuchos::RCP<TreeVector>& soln) :
    glist_(glist),
    Amanzi::MPCStrong<FnTimeIntegratorPK>(pk_tree, glist, S, soln)
{
  Teuchos::ParameterList vlist;
  vo_ = new VerboseObject("FlowEnergy_PK", vlist); 
}


/* ******************************************************************* 
* Physics-based setup of PK.
******************************************************************* */
void FlowEnergy_PK::Setup()
{
  mesh_ = S_->GetMesh();
  int dim = mesh_->space_dimension();

  Teuchos::ParameterList& elist = S_->FEList();

  // solid
  if (!S_->HasField("density_rock")) {
    S_->RequireField("density_rock", "density_rock")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireFieldEvaluator("density_rock");
  }

  if (!S_->HasField("internal_energy_rock")) {
    elist.sublist("internal_energy_rock")
         .set<std::string>("field evaluator type", "iem")
         .set<std::string>("internal energy key", "internal_energy_rock");
    elist.sublist("internal_energy_rock").sublist("IEM parameters")
         .set<std::string>("IEM type", "linear")
         .set<double>("heat capacity [J/mol-K]", 620.0);
  }

  // gas
  if (!S_->HasField("internal_energy_gas")) {
    elist.sublist("internal_energy_gas")
         .set<std::string>("field evaluator type", "iem water vapor")
         .set<std::string>("internal energy key", "internal_energy_gas");
  }

  if (!S_->HasField("molar_density_gas")) {
    elist.sublist("molar_density_gas")
         .set<std::string>("field evaluator type", "eos")
         .set<std::string>("EOS basis", "molar")
         .set<std::string>("molar density key", "molar_density_gas");
    elist.sublist("molar_density_gas").sublist("EOS parameters")
         .set<std::string>("EOS type", "vapor in gas");
    elist.sublist("molar_density_gas").sublist("EOS parameters")
         .sublist("gas EOS parameters")
         .set<std::string>("EOS type", "ideal gas");
  }

  if (!S_->HasField("molar_fraction_gas")) {
    elist.sublist("molar_fraction_gas")
         .set<std::string>("field evaluator type", "molar fraction gas")
         .set<std::string>("molar fraction key", "molar_fraction_gas");
    elist.sublist("molar_fraction_gas")
         .sublist("vapor pressure model parameters")
         .set<std::string>("vapor pressure model type", "water vapor over water/ice");
  }

  // liquid
  if (!S_->HasField("internal_energy_liquid")) {
    elist.sublist("internal_energy_liquid")
         .set<std::string>("field evaluator type", "iem")
         .set<std::string>("internal energy key", "internal_energy_liquid");
    elist.sublist("internal_energy_liquid")
         .sublist("IEM parameters")
         .set<std::string>("IEM type", "linear")
         .set<double>("heat capacity [J/mol-K]", 76.0);
  }

  if (!S_->HasField("molar_density_liquid")) {
    elist.sublist("molar_density_liquid")
         .set<std::string>("field evaluator type", "eos")
         .set<std::string>("EOS basis", "both")
         .set<std::string>("molar density key", "molar_density_liquid")
         .set<std::string>("mass density key", "mass_density_liquid");
    elist.sublist("molar_density_liquid").sublist("EOS parameters")
         .set<std::string>("EOS type", "liquid water");
  }

  if (!S_->HasField("mass_density_liquid")) {
    S_->RequireField("mass_density_liquid", "mass_density_liquid")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireFieldEvaluator("mass_density_liquid");
  }

  // other
  if (!S_->HasField("effective_pressure")) {
    elist.sublist("effective_pressure")
         .set<std::string>("field evaluator type", "effective_pressure");
  }

  if (!S_->HasField("mass_density_liquid")) {
    S_->RequireFieldEvaluator("mass_density_liquid");
  }

  // inform other PKs about strong coupling
  // -- flow
  Teuchos::ParameterList& flow = glist_->sublist("PKs").sublist("Flow")
                                        .sublist("Richards problem")
                                        .sublist("physical models and assumptions");
  flow.set("vapor diffusion", true);
  flow.set<std::string>("water content model", "generic");

  // -- energy
  Teuchos::ParameterList& energy = glist_->sublist("PKs").sublist("Energy")
                                          .sublist("physical models and assumptions");
  energy.set("vapor diffusion", true);

  // process other PKs.
  MPCStrong<FnTimeIntegratorPK>::Setup();
}


/* ******************************************************************* 
* Performs one time step.
******************************************************************* */
bool FlowEnergy_PK::AdvanceStep(double t_old, double t_new)
{
  // flow
  // -- swap saturations
  S_->GetFieldEvaluator("saturation_liquid")->HasFieldChanged(S_.ptr(), "flow");
  const CompositeVector& sat = *S_->GetFieldData("saturation_liquid");
  CompositeVector& sat_prev = *S_->GetFieldData("prev_saturation_liquid", "flow");

  CompositeVector sat_prev_copy(sat_prev);
  sat_prev = sat;

  // -- swap water_contents
  S_->GetFieldEvaluator("water_content")->HasFieldChanged(S_.ptr(), "flow");
  CompositeVector& wc = *S_->GetFieldData("water_content", "water_content");
  CompositeVector& wc_prev = *S_->GetFieldData("prev_water_content", "flow");

  CompositeVector wc_prev_copy(wc_prev);
  wc_prev = wc;

  // energy
  // -- swap conserved energys
  S_->GetFieldEvaluator("energy")->HasFieldChanged(S_.ptr(), "thermal");
  const CompositeVector& e = *S_->GetFieldData("energy");
  CompositeVector& e_prev = *S_->GetFieldData("prev_energy", "thermal");

  CompositeVector e_prev_copy(e_prev);
  e_prev = e;
 
  // try a step
  bool fail = MPCStrong<FnTimeIntegratorPK>::AdvanceStep(t_old, t_new);

  if (fail) {
    // revover the original conserved quantaties
    *S_->GetFieldData("prev_saturation_liquid", "flow") = sat_prev_copy;
    *S_->GetFieldData("prev_water_content", "flow") = wc_prev_copy;
    *S_->GetFieldData("prev_energy", "thermal") = e_prev_copy;

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Step failed. Restored prev_saturation_liquid, prev_water_content, prev_energy" << std::endl;
  }

  return fail;
}

}  // namespace Amanzi

