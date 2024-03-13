/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy(dasvyat@lanl.gov)
*/

/*
  Operators

*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"
#include "GMVMesh.hh"

#include "OperatorDefs.hh"
#include "BoundaryFlux.hh"

namespace Amanzi {

class Model {
 public:
  Model(){};
  ~Model(){};

  // main members
  double Value(double pc) const { return analytic(pc); }

  double analytic(double pc) const { return 1 + pc; }

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};

typedef double (Model::*NonLinFn)(double pc) const;

} // namespace Amanzi


/* *****************************************************************
* This test replaces diffusion tensor and boundary conditions by
* continuous functions.
* **************************************************************** */
TEST(BOUNDARYFLUX)
{
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;
  using namespace Amanzi::Operators;

  auto comm = Amanzi::getDefaultComm();
  int MyPID = comm->MyPID();

  if (MyPID == 0) std::cout << "\nTest: BoundaryFluxSolver" << std::endl;

  // create model of nonlinearity
  Teuchos::RCP<const Model> model = Teuchos::rcp(new Model());
  NonLinFn func = &Model::Value;

  double trans_f = 1;
  double g_f = 2.;
  double lmd = 1;
  double cell_val = 1;
  double bnd_flux = 1;
  double patm = 0.;
  int dir = 1;
  double max_val = 100;
  double min_val = 0;
  double eps = 1e-6;

  BoundaryFaceSolver<Model> BndFaceSolver(
    trans_f, g_f, cell_val, lmd, bnd_flux, dir, patm, max_val, min_val, eps, model, func);

  double face_value;
  face_value = BndFaceSolver.FaceValue();

  if (MyPID == 0) std::cout << "Face value " << face_value << "\n";

  CHECK(fabs(face_value - 0.585786) < eps);
}
