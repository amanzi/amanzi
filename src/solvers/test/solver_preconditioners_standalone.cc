/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include <string>

#include "Teuchos_GlobalMPISession.hpp"
#include "VerboseObject_objs.hh"


#include "Teuchos_RCP.hpp"

#include "exceptions.hh"
#include "IterativeMethodPCG.hh"


#include "AmanziComm.hh"

#include "Preconditioner.hh"
#include "InverseFactory.hh"


const int N = 125;
using namespace Amanzi;
using namespace Amanzi::AmanziSolvers;

inline Teuchos::RCP<Matrix_type>
matrix(const Teuchos::RCP<Map_type>& map)
{
  auto A_ = Teuchos::rcp(new Matrix_type(map, map, 3));

  double v0[2] = { 1.0, -1.0 };
  int inds0[2] = { 0, 1 };
  A_->insertLocalValues(0, 2, v0, inds0);

  for (int i = 1; i < map->getLocalNumElements() - 1; i++) {
    int indices[3];
    double values[3] = { double(-i), double(2 * i + 1), double(-i - 1) };
    for (int k = 0; k < 3; k++) indices[k] = i + k - 1;
    A_->insertLocalValues(i, 3, values, indices);
  }

  int i = map->getLocalNumElements() - 1;
  double vN[2] = { double(-i), double(2 * i + 1) };
  int indsN[2] = { i - 1, i };
  A_->insertLocalValues(i, 2, vN, indsN);

  A_->fillComplete(map, map);

  return A_;
}


inline Teuchos::RCP<Amanzi::AmanziSolvers::Preconditioner>
preconditioner(const std::string& name, const Teuchos::RCP<Matrix_type>& mat)
{
  Teuchos::ParameterList plist;
  plist.set<std::string>("preconditioning method", name);
  Teuchos::ParameterList& tmp = plist.sublist(name + " parameters");

  tmp.sublist("verbose object").set<std::string>("verbosity level", "high");
  // if (name == "ml") {
  //   tmp.set("coarse: max size", 5);
  //   tmp.set("cycle applications", 2);
  //   tmp.set("ML output", 0);
  // } else if (name == "boomer amg") {
  //   tmp.set("max coarse size", 5);
  //   tmp.set("cycle applications", 1);
  // }

  // static_assert(Amanzi::AmanziSolvers::Impl::is_assembled<Matrix_type>::value,
  // "Matrix_type is assembled?");
  auto pc = createAssembledMethod<>(name, plist);
  pc->set_matrix(mat);
  return pc;
};

inline Teuchos::RCP<
  IterativeMethodPCG<Matrix_type, Amanzi::AmanziSolvers::Preconditioner, Vector_type, Map_type>>
get_solver(const std::string& name, const Teuchos::RCP<Matrix_type>& m)
{
  auto pc = preconditioner(name, m);

  Teuchos::ParameterList plist;
  plist.set("error tolerance", 1.e-12);
  plist.set("maximum number of iterations", 200);
  auto inv = Teuchos::rcp(new IterativeMethodPCG<Matrix_type,
                                                 Amanzi::AmanziSolvers::Preconditioner,
                                                 Vector_type,
                                                 Map_type>());
  inv->set_inverse_parameters(plist);
  inv->set_matrices(m, pc);
  inv->initializeInverse();
  inv->computeInverse();
  return inv;
}

int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  std::cout << "\nComparison of preconditioners for N=125" << std::endl;

  auto comm = getDefaultComm();
  auto map = Teuchos::rcp(new Map_type(N, 0, comm));
  auto m = matrix(map);


  Vector_type u(m->getRangeMap());
  {
    auto uv = u.getLocalViewHost();
    for (int i = 0; i < N; i++) uv(i, 0) = 1.0 / (i + 2.0);
  }
  u.sync_device();

  Vector_type v(m->getDomainMap());

  for (const auto& prec_name : {
         "identity", "diagonal", "ifpack2: ILUT", "ifpack2: FAST_ILU"
         //"ifpack2: KSPILUK"
         //"boomer amg",
         //"euclid",
         //"ml"
       }) {
    auto solver = get_solver(prec_name, m);
    v.putScalar(0.0);

    std::cout << "Preconditioner: " << prec_name << std::endl
              << "-------------------------------------------" << std::endl;
    solver->applyInverse(u, v);

    v.sync_host();
    auto vv = v.getLocalViewHost();
    assert(fabs(11.03249773994628 - vv(0, 0)) < 1e-6);
    assert(fabs(10.53249773994628 - vv(1, 0)) < 1e-6);
  }
}
