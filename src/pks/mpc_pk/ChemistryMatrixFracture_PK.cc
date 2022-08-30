/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy

  Process kernel that couples chemistry PKs in matrix and fracture.
*/

#include "ChemistryMatrixFracture_PK.hh"
#include "PK_MPCWeak.hh"

namespace Amanzi {

/* ******************************************************************* 
* Constructor
******************************************************************* */
ChemistryMatrixFracture_PK::ChemistryMatrixFracture_PK(Teuchos::ParameterList& pk_tree,
                                                       const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                                       const Teuchos::RCP<State>& S,
                                                       const Teuchos::RCP<TreeVector>& soln) :
    Amanzi::PK(pk_tree, glist, S, soln),
    Amanzi::PK_MPCWeak(pk_tree, glist, S, soln)
{
  Teuchos::ParameterList vlist;
  vo_ = Teuchos::rcp(new VerboseObject("ChemistryMatrixFracture_PK", vlist)); 
}


/* ******************************************************************* 
* Physics-based setup of PK.
******************************************************************* */
void ChemistryMatrixFracture_PK::Setup()
{
  mesh_domain_ = S_->GetMesh();
  mesh_fracture_ = S_->GetMesh("fracture");

  // setup the sub-PKs
  PK_MPCWeak::Setup();
}

}  // namespace Amanzi

