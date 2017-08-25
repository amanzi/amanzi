/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Explicit.hh"
#include "PK_MixinLeaf.hh"

using namespace Amanzi;

class PK_ODE_Explicit : public PK_MixinLeaf<PK_Explicit<PK> > {

};

Teuchos::RCP<PK>
create() {
  Teuchos::RCP<Teuchos::ParameterList> pk_tree = Teuchos::rcp(new Teuchos::ParameterList("my pk tree"));
  Teuchos::RCP<Teuchos::ParameterList> global_list = Teuchos::rcp(new Teuchos::ParameterList("main"));
  Teuchos::RCP<State> S = Teuchos::rcp(new State());
  Teuchos::RCP<TreeVector> soln = Teuchos::rcp(new TreeVector());
  return Teuchos::rcp(new PK_ODE_Explicit(pk_tree, global_list, S, soln));
}


SUITE(PK_ODE_EXPLICIT) {

  TEST_FIXTURE(Construct) {
    create();
  }

}

