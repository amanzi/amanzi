/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "UnitTest++.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "BilinearFormFactory.hh"

/* ****************************************************************
* Test of face centroids
**************************************************************** */
TEST(FACTORY_DISCRETIZATION_METHODS) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: factory of discretization methods" << std::endl;
 
  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
  MeshFactory factory(comm);
  factory.preference(pref);
  // Teuchos::RCP<const Mesh> mesh = factory(0.0, 0.0, 1.0, 1.0, 1, 1); 
  Teuchos::RCP<const Mesh> mesh = factory("test/one_pentagon.exo");

  std::string names[8] = {"diffusion", "diffusion generalized", "elasticity",
                          "CrouzeixRaviart", "BernardiRaugel", "Lagrange serendipity",
                          "Lagrange", "dg modal"};

  for (int i = 0; i < 8; ++i) {
    Teuchos::ParameterList plist;
    plist.set<std::string>("method", names[i])
         .set<int>("method order", 1)
         .set<std::string>("dg basis", "normalized");
    Teuchos::RCP<BilinearForm> form = BilinearFormFactory::Create(plist, mesh); 
    CHECK(form != Teuchos::null);
  }

  delete comm;
}


