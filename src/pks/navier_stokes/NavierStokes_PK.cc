/*
  Navier Stokes PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// Amanzi::NavierStokes
#include "NavierStokes_PK.hh"

namespace Amanzi {
namespace NavierStokes {

/* ******************************************************************
* New constructor: extracts lists and requires fields.
****************************************************************** */
NavierStokes_PK::NavierStokes_PK(Teuchos::ParameterList& pk_tree,
                                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                 const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<TreeVector>& soln) :
  soln_(soln)
{
  S_ = S;

  std::string pk_name = pk_tree.name();
  const char* result = pk_name.data();

  boost::iterator_range<std::string::iterator> res = boost::algorithm::find_last(pk_name,"->"); 
  if (res.end() - pk_name.end() != 0) boost::algorithm::erase_head(pk_name,  res.end() - pk_name.begin());

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  ns_list_ = Teuchos::sublist(pk_list, pk_name, true);

  // We also need iscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_solver_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(ns_list_, "time integrator", true);
}


/* ******************************************************************
* Old constructor for unit tests.
****************************************************************** */
NavierStokes_PK::NavierStokes_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                 const std::string& pk_list_name,
                                 Teuchos::RCP<State> S,
                                 const Teuchos::RCP<TreeVector>& soln) :
    glist_(glist),
    soln_(soln)
{
  S_ = S;

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  ns_list_ = Teuchos::sublist(pk_list, "Navier Stokes", true);
 
  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_solver_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(ns_list_, "time integrator");

  vo_ = Teuchos::null;
}


/* ******************************************************************
* Define structure of this PK. We request physical fields and their
* evaluators. Selection of a few models is available and driven by
* model factories, evaluator factories, and parameters of the list
* "physical models and assumptions".
****************************************************************** */
void NavierStokes_PK::Setup(const Teuchos::Ptr<State>& S)
{
  dt_ = 0.0;
  mesh_ = S->GetMesh();
  dim = mesh_->space_dimension();

  // primary fields
  // -- pressure
  if (!S->HasField("pressure")) {
    S->RequireField("pressure", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", "pressure");
    pressure_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator("pressure", pressure_eval_);
  }

  // -- velocity
  std::vector<std::string> names = {"node", "cell"};
  std::vector<AmanziMesh::Entity_kind> locations = {AmanziMesh::NODE, AmanziMesh::FACE};
  std::vector<int> ndofs = {dim, 1};

  if (!S->HasField("fluid_velocity")) {
    S->RequireField("fluid_velocity", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", "fluid_velocity");
    fluid_velocity_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator("fluid_velocity", fluid_velocity_eval_);
  }
}


/* ******************************************************************
* This is a long but simple routine. It goes through flow parameter
* list and initializes various objects including those created during 
* the setup step.
****************************************************************** */
void NavierStokes_PK::Initialize(const Teuchos::Ptr<State>& S)
{
}

}  // namespace NavierStokes
}  // namespace Amanzi
