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
  mesh_ = S->GetMesh();
  int dim = mesh_->space_dimension();

  Teuchos::ParameterList& elist = S->FEList();

  // Fields for coupling terms
  // -- bc for matrix
  /*
  if (!S->HasField("matrix_fracture_bc")) {
    elist.sublist("matrix_fracture_bc")
         .set<std::string>("field evaluator type", "matrix_fracture_bc");

    S->RequireField("matrix_fracture_bc", "matrix_fracture_bc")->SetMesh(mesh_)
      ->SetGhosted(true)->SetComponent("boundary face", AmanziMesh::FACE, 2);
    // S->RequireFieldEvaluator("matrix_fracture_bc");
    S->GetField("matrix_fracture_bc", "matrix_fracture_bc")->set_io_vis(false);
  }
  */

  // -- source for fracture

  // inform other PKs about strong coupling
  // -- flow (matrix)
  Teuchos::ParameterList& mflow = glist_->sublist("PKs").sublist("flow matrix")
                                         .sublist("Flow problem")
                                         .sublist("physical models and assumptions");
  mflow.set("coupled matrix fracture flow", true);

  // -- flow (fracture)
  Teuchos::ParameterList& fflow = glist_->sublist("PKs").sublist("flow fracture")
                                         .sublist("Flow problem")
                                         .sublist("physical models and assumptions");
  fflow.set("coupled matrix fracture flow", true);

  // process other PKs.
  PK_MPCStrong<PK_BDF>::Setup(S);
}


/* ******************************************************************* 
* Performs one time step.
******************************************************************* */
bool FlowMatrixFracture_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  // try a step
  bool fail = PK_MPCStrong<PK_BDF>::AdvanceStep(t_old, t_new, reinit);

  if (fail) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Step failed." << std::endl;
  }

  return fail;
}

}  // namespace Amanzi

