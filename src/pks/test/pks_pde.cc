/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

/* Test basic implicit and explicit PDEs

At this point PKs manage memory and interface time integrators with the DAG.
These tests that functionality with a series of ODEs.
   
*/

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Default.hh"
#include "PK_MixinLeaf.hh"
#include "PK_MixinExplicit.hh"
#include "PK_MixinExplicitSubcycled.hh"
#include "PK_MixinImplicit.hh"
#include "PK_MixinImplicitSubcycled.hh"
#include "PK_MixinPredictorCorrector.hh"
#include "PK_Adaptors.hh"


#include "test_pks.hh"
#include "pks_test_harness.hh"

using namespace Amanzi;



SUITE(PKS_PDE) {

  // Forward Euler tests with each of 3 PKs
  TEST(A_FORWARD_EULER) {
    auto run = createExplicit("A", "forward euler");
    auto nsteps = run_test(run->S, run->pk);
    CHECK_CLOSE(2.0, (*run->S->Get<CompositeVector>("primaryA")
                      .ViewComponent("cell",false))[0][0], 1.e-10);
    CHECK_EQUAL(10, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }
  
}

