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
#include "PK_MPCStrong.hh"

namespace Amanzi {

/* ******************************************************************* 
* Constructor
******************************************************************* */
FlowEnergy_PK::FlowEnergy_PK(Teuchos::ParameterList& pk_tree,
                             const Teuchos::RCP<Teuchos::ParameterList>& glist,
                             const Teuchos::RCP<State>& S,
                             const Teuchos::RCP<TreeVector>& soln) :
    glist_(glist),
    Amanzi::PK_MPC<PK_BDF>(pk_tree, glist, S, soln),
    Amanzi::PK_MPCStrong<PK_BDF>(pk_tree, glist, S, soln)
{
  Teuchos::ParameterList vlist;
  vo_ =  Teuchos::rcp(new VerboseObject("FlowEnergy_PK", vlist)); 
}


/* ******************************************************************* 
* Physics-based setup of PK.
******************************************************************* */
void FlowEnergy_PK::Setup(const Teuchos::Ptr<State>& S)
{
  mesh_ = S->GetMesh();
  int dim = mesh_->space_dimension();

  Teuchos::ParameterList& elist = S->FEList();

  // Fields for solids
  // -- rock
  if (!S->HasField("particle_density")) {
    S->RequireField("particle_density", "particle_density")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator("particle_density");
  }

  if (!S->HasField("internal_energy_rock")) {
    elist.sublist("internal_energy_rock")
         .set<std::string>("field evaluator type", "iem")
         .set<std::string>("internal energy key", "internal_energy_rock");
    elist.sublist("internal_energy_rock").sublist("IEM parameters")
         .set<std::string>("iem type", "linear")
         .set<double>("heat capacity [J/mol-K]", 620.0);
  }

  // Fields for gas
  // -- internal energy
  if (!S->HasField("internal_energy_gas")) {
    elist.sublist("internal_energy_gas")
         .set<std::string>("field evaluator type", "iem water vapor")
         .set<std::string>("internal energy key", "internal_energy_gas");
  }

  // -- molar density
  if (!S->HasField("molar_density_gas")) {
    elist.sublist("molar_density_gas")
         .set<std::string>("field evaluator type", "eos")
         .set<std::string>("eos basis", "molar")
         .set<std::string>("molar density key", "molar_density_gas");
    elist.sublist("molar_density_gas").sublist("EOS parameters")
         .set<std::string>("eos type", "vapor in gas");
    elist.sublist("molar_density_gas").sublist("EOS parameters")
         .sublist("gas EOS parameters")
         .set<std::string>("eos type", "ideal gas")
         .set<double>("molar mass of gas", 28.9647e-03);  // dry air
  }

  // -- molar fraction
  if (!S->HasField("molar_fraction_gas")) {
    elist.sublist("molar_fraction_gas")
         .set<std::string>("field evaluator type", "molar fraction gas")
         .set<std::string>("molar fraction key", "molar_fraction_gas");
    elist.sublist("molar_fraction_gas")
         .sublist("vapor pressure model parameters")
         .set<std::string>("vapor pressure model type", "water vapor over water/ice");
  }

  // Fields for liquid
  // -- internal energy
  if (!S->HasField("internal_energy_liquid")) {
    elist.sublist("internal_energy_liquid")
         .set<std::string>("field evaluator type", "iem")
         .set<std::string>("internal energy key", "internal_energy_liquid");
    elist.sublist("internal_energy_liquid")
         .sublist("IEM parameters")
         .set<std::string>("iem type", "linear")
         .set<double>("heat capacity [J/mol-K]", 76.0);
  }

  // -- molar and mass density
  if (!S->HasField("molar_density_liquid")) {
    elist.sublist("molar_density_liquid")
         .set<std::string>("field evaluator type", "eos")
         .set<std::string>("eos basis", "both")
         .set<std::string>("molar density key", "molar_density_liquid")
         .set<std::string>("mass density key", "mass_density_liquid");
    elist.sublist("molar_density_liquid").sublist("EOS parameters")
         .set<std::string>("eos type", "liquid water");
    elist.sublist("molar_density_liquid")
         .sublist("verbose object").set<std::string>("verbosity level", "medium");

    S->RequireField("molar_density_liquid", "molar_density_liquid")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator("molar_density_liquid");
    S->RequireFieldEvaluator("mass_density_liquid");
  }

  // -- viscosity model
  if (!S->HasField("viscosity_liquid")) {
    elist.sublist("viscosity_liquid")
         .set<std::string>("field evaluator type", "viscosity")
         .set<std::string>("viscosity key", "viscosity_liquid")
         .sublist("viscosity model parameters")
         .set<std::string>("viscosity relation type", "liquid water");
    elist.sublist("viscosity_liquid")
         .sublist("verbose object").set<std::string>("verbosity level", "high");

    S->RequireField("viscosity_liquid", "viscosity_liquid")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator("viscosity_liquid");
  }

  // Other fields
  if (!S->HasField("effective_pressure")) {
    elist.sublist("effective_pressure")
         .set<std::string>("field evaluator type", "effective_pressure");

    S->RequireField("effective_pressure", "effective_pressure")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator("effective_pressure");
    S->GetField("effective_pressure", "effective_pressure")->set_io_vis(false);
  }

  // inform other PKs about strong coupling
  // -- flow
  Teuchos::ParameterList& flow = glist_->sublist("PKs").sublist("flow")
                                        .sublist("Richards problem")
                                        .sublist("physical models and assumptions");
  flow.set("vapor diffusion", true);
  flow.set<std::string>("water content model", "generic");

  // -- energy
  Teuchos::ParameterList& energy = glist_->sublist("PKs").sublist("energy")
                                          .sublist("physical models and assumptions");
  energy.set("vapor diffusion", true);

  // process other PKs.
  PK_MPCStrong<PK_BDF>::Setup(S);
}


/* ******************************************************************* 
* Performs one time step.
******************************************************************* */
bool FlowEnergy_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  // flow
  // -- swap saturations (current and previous)
  S_->GetFieldEvaluator("saturation_liquid")->HasFieldChanged(S_.ptr(), "flow");
  const CompositeVector& sat = *S_->GetFieldData("saturation_liquid");
  CompositeVector& sat_prev = *S_->GetFieldData("prev_saturation_liquid", "flow");

  CompositeVector sat_prev_copy(sat_prev);
  sat_prev = sat;

  // -- swap water_contents (current and previous)
  S_->GetFieldEvaluator("water_content")->HasFieldChanged(S_.ptr(), "flow");
  CompositeVector& wc = *S_->GetFieldData("water_content", "water_content");
  CompositeVector& wc_prev = *S_->GetFieldData("prev_water_content", "flow");

  CompositeVector wc_prev_copy(wc_prev);
  wc_prev = wc;

  // energy
  // -- swap conserved energies (current and previous)
  S_->GetFieldEvaluator("energy")->HasFieldChanged(S_.ptr(), "thermal");
  const CompositeVector& e = *S_->GetFieldData("energy");
  CompositeVector& e_prev = *S_->GetFieldData("prev_energy", "thermal");

  CompositeVector e_prev_copy(e_prev);
  e_prev = e;
 
  // try a step
  bool fail = PK_MPCStrong<PK_BDF>::AdvanceStep(t_old, t_new, reinit);

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

