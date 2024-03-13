/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, version 2.0
  Release name: naka-to.

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
TEST(FACTORY_DISCRETIZATION_METHODS)
{
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::WhetStone;

  std::cout << "Test: meshfactory of discretization methods" << std::endl;

  Preference pref;
  pref.clear();
  pref.push_back(Framework::MSTK);

  auto comm = Amanzi::getDefaultComm();
  MeshFactory meshfactory(comm);
  meshfactory.set_preference(pref);
  // Teuchos::RCP<const Mesh> mesh = meshfactory.create(0.0, 0.0, 1.0, 1.0, 1, 1);
  Teuchos::RCP<const Mesh> mesh = meshfactory.create("test/one_pentagon.exo");

  std::string names[8] = { "diffusion",      "diffusion generalized",
                           "elasticity",     "CrouzeixRaviart",
                           "BernardiRaugel", "Lagrange serendipity",
                           "Lagrange",       "dg modal" };

  for (int i = 0; i < 8; ++i) {
    Teuchos::ParameterList plist;
    plist.set<std::string>("method", names[i])
      .set<int>("method order", 1)
      .set<std::string>("base", "cell")
      .set<std::string>("dg basis", "normalized")
      .set<int>("quadrature order", 1);
    Teuchos::RCP<BilinearForm> form = BilinearFormFactory::Create(plist, mesh);
    CHECK(form != Teuchos::null);
  }
}
