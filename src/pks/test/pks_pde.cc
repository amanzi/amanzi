/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

/* Test basic implicit and explicit PDEs

At this point PKs manage memory and interface time integrators with the DAG.
These tests that functionality with a series of ODEs.

Solves the pseudo-1D PDE:

du/dt - div k grad u = q

on (-1,1)

with initial data u = 0.75,
boundary data u = 0
source data q = 0.2 * (pi/2)^2 * cos(pi/2 x)
coefficient k = 0.2

At steady state, the solution is:

u = cos(pi/2 x)

*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "CompositeVector.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Adaptors.hh"
#include "PK_Default.hh"
#include "PK_MixinExplicit.hh"
#include "PK_MixinExplicitSubcycled.hh"
#include "PK_MixinImplicit.hh"
#include "PK_MixinImplicitSubcycled.hh"
#include "PK_MixinLeaf.hh"
#include "PK_MixinPredictorCorrector.hh"

#include "pks_test_harness.hh"
#include "test_pk_pde.hh"

using namespace Amanzi;

// static const double PI_2 = 1.5707963267948966;

SUITE(PKS_PDE)
{
  // // Forward Euler tests with each of 3 PKs
  // TEST(DIFFUSION_FE_EXPLICIT)
  // {
  //   // run the simulation
  //   using PK_t = PK_Explicit_Adaptor<PK_PDE_Explicit<
  //     PK_MixinExplicit<PK_MixinLeafCompositeVector<PK_Default>>>>;
  //   auto run = createRunPDE<PK_t>("diffusion FE explicit", "test/pks_pde.xml");
  //   auto nsteps = run_test(run->S, run->pk, 20.0);

  //   // print final solution
  //   std::cout << "Final solution" << std::endl;
  //   run->S->Get<CompositeVector>("u").Print(std::cout);

  //   // check error
  //   auto m = run->S->GetMesh();
  //   int ncells =
  //     m->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  //   {
  //     auto u = run->S->Get<CompositeVector>("u").ViewComponent<AmanziDefaultHost>("cell", false);
  //     for (int c = 0; c != ncells; ++c) {
  //       auto p = m->cell_centroid(c);
  //       double val = std::cos(PI_2 * p[0]);
  //       CHECK_CLOSE(val, u(c,0), 1.e-3);
  //     }
  //   }

  //   // check timesteps
  //   CHECK_EQUAL(20. / 0.001, nsteps.first);
  //   CHECK_EQUAL(0, nsteps.second);
  // }

  // Forward Euler tests with each of 3 PKs
  TEST(DIFFUSION_FE_IMPLICIT)
  {
    // run the simulation
    using PK_t = PK_Implicit_Adaptor<PK_PDE_Implicit<
      PK_MixinImplicit<PK_MixinLeafCompositeVector<PK_Default>>>>;
    auto run = createRunPDE<PK_t>("diffusion FE implicit", "test/pks_pde.xml");
    auto nsteps = run_test(run->S, run->pk, 20.0);

    // print final solution
    std::cout << "Final solution" << std::endl;
    run->S->Get<CompositeVector>("u").Print(std::cout);

    // check error
    auto m = run->S->GetMesh();
    int ncells =
      m->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    {
      auto u = run->S->Get<CompositeVector>("u").ViewComponent<AmanziDefaultHost>("cell", false);
      for (int c = 0; c != ncells; ++c) {
        auto p = m->cell_centroid(c);
        double val = std::cos(PI_2 * p[0]);
        CHECK_CLOSE(val, u(c,0), 1.e-3);
      }
    }
    
    // check timesteps
    CHECK_EQUAL(20. / 1, nsteps.first);
    CHECK_EQUAL(0, nsteps.second);
  }
}
