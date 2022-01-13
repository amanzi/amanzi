/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy

  Process kernel that couples Transport in matrix and fracture.
*/

#include "PDE_DiffusionFracturedMatrix.hh"
#include "EvaluatorPrimary.hh"

#include "TransportMatrixFracture_PK.hh"
#include "PK_MPCStrong.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************* 
* Constructor
******************************************************************* */
TransportMatrixFracture_PK::TransportMatrixFracture_PK(Teuchos::ParameterList& pk_tree,
                                                       const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                                       const Teuchos::RCP<State>& S,
                                                       const Teuchos::RCP<TreeVector>& soln) :
    glist_(glist),
    Amanzi::PK(pk_tree, glist, S, soln),
    Amanzi::PK_MPCWeak(pk_tree, glist, S, soln)
{
  Teuchos::ParameterList vlist;
  vo_ = Teuchos::rcp(new VerboseObject("TransportMatrixFracture_PK", vlist)); 
}


/* ******************************************************************* 
* Physics-based setup of PK.
******************************************************************* */
void TransportMatrixFracture_PK::Setup()
{
  mesh_domain_ = S_->GetMesh();
  mesh_fracture_ = S_->GetMesh("fracture");

  // darcy fluxes use non-uniform distribution of DOFs
  // -- darcy flux for matrix
  if (!S_->HasData("darcy_flux")) {
    auto cvs = Operators::CreateFracturedMatrixCVS(mesh_domain_, mesh_fracture_);
    auto mmap = cvs->Map("face", false);
    auto gmap = cvs->Map("face", true);
    S_->Require<CV_t, CVS_t>("darcy_flux", Tags::DEFAULT, "transport")
      .SetMesh(mesh_domain_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, mmap, gmap, 1);
  }

  // -- darcy flux for fracture
  if (!S_->HasData("fracture-darcy_flux")) {
    auto cvs = Operators::CreateNonManifoldCVS(mesh_fracture_);
    *S_->Require<CV_t, CVS_t>("fracture-darcy_flux", Tags::DEFAULT, "transport")
      .SetMesh(mesh_fracture_)->SetGhosted(true) = *cvs;
  }

  // add boundary condition to transport in matrix list
  auto pks = glist_->sublist("PKs").sublist(name_).get<Teuchos::Array<std::string> >("PKs order").toVector();
  Teuchos::ParameterList& bclist = glist_->sublist("PKs")
      .sublist(pks[0]).sublist("boundary conditions")
      .sublist("concentration").sublist("coupling").sublist("BC coupling");

   Teuchos::Array<std::string> regs;
   regs.push_back("FRACTURE_NETWORK_INTERNAL");
   bclist.set<std::string>("spatial distribution method", "domain coupling")
         .set<std::string>("submodel", "field")
         .set<Teuchos::Array<std::string> >("regions", regs);

   bclist.sublist("boundary concentration")
         .set<std::string>("external field key", "fracture-total_component_concentration");

  // add source term to transport in fracture list
  Teuchos::ParameterList& srclist = glist_->sublist("PKs")
      .sublist(pks[1]).sublist("source terms")
      .sublist("concentration").sublist("coupling").sublist("fracture");

   regs.clear();
   regs.push_back("All");
   srclist.set<std::string>("spatial distribution method", "domain coupling")
          .set<std::string>("submodel", "rate")
          .set<Teuchos::Array<std::string> >("regions", regs);

   srclist.sublist("sink")
         .set<std::string>("external field key", "total_component_concentration")
         .set<std::string>("flux key", "darcy_flux");

  // setup the sub-PKs
  PK_MPCWeak::Setup();
}


/* ******************************************************************* 
* Reduce stable dt to avoid the 2-cycle behavior of the transport PK
******************************************************************* */
double TransportMatrixFracture_PK::get_dt() {
  return 0.999 * PK_MPCWeak::get_dt();
}

  
/* ******************************************************************* 
* Performs one time step.
******************************************************************* */
bool TransportMatrixFracture_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = PK_MPCWeak::AdvanceStep(t_old, t_new, reinit);

  if (fail) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Step failed." << std::endl;
  }

  return fail;
}

}  // namespace Amanzi

