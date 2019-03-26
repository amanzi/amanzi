/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   Amanzi

   License: see $AMANZI_DIR/COPYRIGHT
   Author: Ethan Coon

   Unit tests for the communication patterns used in CompositeVector
------------------------------------------------------------------------- */
#include <vector>

#include "UnitTest++.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"

#include "AmanziComm.hh"
#include "AmanziMap.hh"
#include "AmanziVector.hh"

using namespace Amanzi;


TEST(COMMUNICATION_PATTERN_DISTINCT_VECTORS) {
  // having some trouble with communication -- lets make sure we understand
  // how to use tpetra
  auto comm = getDefaultComm();

  int size = comm->getSize(); CHECK_EQUAL(2,size);
  int rank = comm->getRank();

  Map_ptr_type owned_map, ghost_map;
  if (rank == 0) {
    std::vector<int> gids_owned{0,1};
    owned_map = Teuchos::rcp(new Map_type(4, gids_owned.data(), 2, 0, comm));

    std::vector<int> gids_used{0,1,2};
    ghost_map = Teuchos::rcp(new Map_type(6, gids_used.data(), 3, 0, comm));
  } else {
    std::vector<int> gids_owned{2,3};
    owned_map = Teuchos::rcp(new Map_type(4, gids_owned.data(), 2, 0, comm));

    std::vector<int> gids_used{2,3,1};
    ghost_map = Teuchos::rcp(new Map_type(6, gids_used.data(), 3, 0, comm));
  }

  Import_type importer(owned_map, ghost_map);

  auto ghost = Teuchos::rcp(new MultiVector_type(ghost_map, 1));
  //  auto owned = ghost->offsetViewNonConst(owned_map, 0);
  auto owned = Teuchos::rcp(new MultiVector_type(owned_map, 1));

  owned->putScalar((double)(rank+1));
  ghost->putScalar(0.);
  ghost->doImport(*owned, importer, Tpetra::INSERT);

  // manually get a host view and check
  ghost->sync_host();
  {
    // NOTE INDICES HERE -- they are backwards relative to Epetra!
    auto ghost_v = ghost->getLocalViewHost();
    if (rank == 0) {
      CHECK_CLOSE(1.0, ghost_v(0,0), 1.e-6);
      CHECK_CLOSE(2.0, ghost_v(2,0), 1.e-6);
    } else {
      CHECK_CLOSE(2.0, ghost_v(0,0), 1.e-6);
      CHECK_CLOSE(1.0, ghost_v(2,0), 1.e-6);
    }
  }
}


