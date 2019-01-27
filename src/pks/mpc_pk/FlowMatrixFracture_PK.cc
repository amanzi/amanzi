/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy

  Process kernel that couples Flow flow in matrix and fractures.
*/

#include "PDE_DiffusionFracturedMatrix.hh"
#include "primary_variable_field_evaluator.hh"

#include "FlowMatrixFracture_PK.hh"
#include "PK_MPCStrong.hh"

namespace Amanzi {

/* ******************************************************************* 
* Constructor
******************************************************************* */
FlowMatrixFracture_PK::FlowMatrixFracture_PK(Teuchos::ParameterList& pk_tree,
                                               const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                               const Teuchos::RCP<State>& S,
                                               const Teuchos::RCP<TreeVector>& soln) :
    glist_(glist),
    Amanzi::PK_MPC<PK_BDF>(pk_tree, glist, S, soln),
    Amanzi::PK_MPCStrong<PK_BDF>(pk_tree, glist, S, soln)
{
  Teuchos::ParameterList vlist;
  vo_ =  Teuchos::rcp(new VerboseObject("MatrixFracture_PK", vlist)); 
}


/* ******************************************************************* 
* Physics-based setup of PK.
******************************************************************* */
void FlowMatrixFracture_PK::Setup(const Teuchos::Ptr<State>& S)
{
  mesh_domain_ = S->GetMesh();
  mesh_fracture_ = S->GetMesh("fracture");
  int dim = mesh_domain_->space_dimension();

  Teuchos::ParameterList& elist = S->FEList();

  // primary and secondary fields for matrix affected by non-uniform
  // distribution of DOFs
  // -- pressure
  auto cvs = Operators::CreateFracturedMatrixCVS(mesh_domain_, mesh_fracture_);
  if (!S->HasField("pressure")) {
    *S->RequireField("pressure", "flow")->SetMesh(mesh_domain_)->SetGhosted(true) = *cvs;
  }

  // -- darcy flux
  if (!S->HasField("darcy_flux")) {
    std::string name("face");
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S->RequireField("darcy_flux", "flow")->SetMesh(mesh_domain_)->SetGhosted(true) 
      ->SetComponent(name, mmap, gmap, 1);
  }

  // Require additional fields and evaluators

  // inform dependent PKs about coupling
  // -- flow (matrix)
  Teuchos::ParameterList& mflow = glist_->sublist("PKs").sublist("flow matrix")
                                         .sublist("Darcy problem")
                                         .sublist("physical models and assumptions");
  mflow.set<std::string>("coupled matrix fracture flow", "matrix");

  // -- flow (fracture)
  Teuchos::ParameterList& fflow = glist_->sublist("PKs").sublist("flow fracture")
                                         .sublist("Darcy problem")
                                         .sublist("physical models and assumptions");
  fflow.set<std::string>("coupled matrix fracture flow", "fracture");

  // process other PKs.
  PK_MPCStrong<PK_BDF>::Setup(S);
}


/* ******************************************************************* 
* Performs one time step.
******************************************************************* */
bool FlowMatrixFracture_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = PK_MPCStrong<PK_BDF>::AdvanceStep(t_old, t_new, reinit);

  if (fail) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Step failed." << std::endl;
  }

  return fail;
}

}  // namespace Amanzi

