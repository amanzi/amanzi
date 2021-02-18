#include <iostream>
#include <string>

#ifdef _OPENMP
#include "omp.h"
#endif

#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "exceptions.hh"
#include "IterativeMethodPCG.hh"


#include "AmanziComm.hh"

#include "Preconditioner.hh"
#include "InverseFactory.hh"

SUITE(SOLVERS) {
const int N = 125;
using namespace Amanzi;
using namespace Amanzi::AmanziSolvers;

inline Teuchos::RCP<Matrix_type> matrix(const Teuchos::RCP<Map_type>& map) {
    auto A_ = Teuchos::rcp(new Matrix_type(map, map, 3));

    double v0[2] = {1.0, -1.0};
    int inds0[2] = {0,1};
    A_->insertLocalValues(0, 2, v0, inds0);

    for (int i = 1; i < map->getNodeNumElements()-1; i++) {
      int indices[3];
      double values[3] = { double(-i), double(2 * i + 1), double(-i - 1) };
      for (int k = 0; k < 3; k++) indices[k] = i + k - 1;
      A_->insertLocalValues(i, 3, values, indices);
    }

    int i = map->getNodeNumElements()-1;
    double vN[2] = {double(-i), double(2*i+1)};
    int indsN[2] = {i-1,i};
    A_->insertLocalValues(i, 2, vN, indsN);

    A_->fillComplete(map, map);

  return A_;
}


inline Teuchos::RCP<Amanzi::AmanziSolvers::Preconditioner>
preconditioner(const std::string& name,
               const Teuchos::RCP<Matrix_type>& mat)
{
  Teuchos::ParameterList plist;
  plist.set<std::string>("preconditioning method", name);
  Teuchos::ParameterList& tmp = plist.sublist(name+" parameters");

  tmp.sublist("verbose object").set<std::string>("verbosity level", "high");

  if (name == "ifpack2: RILUK") {
    tmp.set("fact: type", "KSPILUK"); // get the kokkos kernels variant
  }

  //static_assert(Amanzi::AmanziSolvers::Impl::is_assembled<Matrix_type>::value, "Matrix_type is assembled?");
  auto pc = createAssembledMethod<>(name, plist);
  pc->set_matrix(mat);
  return pc;
};

inline Teuchos::RCP<IterativeMethodPCG<Matrix_type,Amanzi::AmanziSolvers::Preconditioner,Vector_type,Map_type>>
    get_solver(const std::string& name, const Teuchos::RCP<Matrix_type>& m)
{
  auto pc = preconditioner(name, m);

  Teuchos::ParameterList plist;
  plist.set("error tolerance", 1.e-12);
  plist.set("maximum number of iterations", 200);
  auto inv = Teuchos::rcp(new IterativeMethodPCG<Matrix_type,Amanzi::AmanziSolvers::Preconditioner,Vector_type,Map_type>());
  inv->set_inverse_parameters(plist);
  inv->set_matrices(m, pc);
  inv->initializeInverse();
  inv->computeInverse();
  return inv;
}


TEST(PRECONDITIONERS) {
  std::cout << "\nComparison of preconditioners for N=125" << std::endl;

  auto comm = getDefaultComm();
  auto map = Teuchos::rcp(new Map_type(N, 0, comm));
  auto m = matrix(map);


  Vector_type u(m->getRangeMap());
  {
    auto uv = u.getLocalViewHost();
    for (int i = 0; i < N; i++) uv(i,0) = 1.0 / (i + 2.0);
  }
  u.sync_device();

  Vector_type v(m->getDomainMap());

  for (const auto& prec_name : {
    "identity",
    "diagonal",
    "ifpack2: ILUT",
    "ifpack2: RILUK",
    "hypre: boomer amg",
    "hypre: euclid",
    "muelu"
  }) {
    auto solver = get_solver(prec_name, m);
    v.putScalar(0.0);

    std::cout << "Preconditioner: " << prec_name << std::endl
              << "-------------------------------------------" << std::endl;
    solver->applyInverse(u, v);

    v.sync_host();
    auto vv = v.getLocalViewHost();
    CHECK_CLOSE(11.03249773994628, vv(0,0), 1e-6);
    CHECK_CLOSE(10.53249773994628, vv(1,0), 1e-6);
    std::cout << std::endl << std::endl;
  }
};

#if 0 
#ifdef _OPENMP
TEST(PRECONDITIONERS_OMP) {
  std::cout << "\nComparison of preconditioners for N=125" << std::endl;

  auto comm = getDefaultComm();
  auto map = Teuchos::rcp(new Map_type(N, 0, comm));

  double cpu0 = omp_get_wtime();

#pragma omp parallel for shared(map) num_threads(2)
  for (const auto& prec_name : {"identity", "diagonal"}) {
    auto m = matrix(map);

    Vector_type u(*map), v(*map);
    for (int i = 0; i < N; i++) u[i] = 1.0 / (i + 2.0);

    auto solver = get_solver(prec_name, m);
    v.putScalar(0.0);

    solver->applyInverse(u, v);
    CHECK_CLOSE(11.03249773994628, v[0], 1e-6);
    CHECK_CLOSE(10.53249773994628, v[1], 1e-6);
  }


  int nthreads = omp_get_max_threads();
  double cpu1 = omp_get_wtime();
  std::cout << "CPU (parallel): " << cpu1 - cpu0 << " [sec]  threads=" << nthreads << std::endl;

  // serial run
  cpu0 = omp_get_wtime();
  for (const auto& prec_name : {"identity", "diagonal"}) {
    auto m = matrix(map);

    Vector_type u(*map), v(*map);
    for (int i = 0; i < N; i++) u[i] = 1.0 / (i + 2.0);

    auto solver = get_solver(prec_name, m);
    v.putScalar(0.0);

    solver->applyInverse(u, v);
    CHECK_CLOSE(11.03249773994628, v[0], 1e-6);
    CHECK_CLOSE(10.53249773994628, v[1], 1e-6);
  }


  int nthreads = omp_get_max_threads();
  double cpu1 = omp_get_wtime();
  std::cout << "CPU (serial): " << cpu1 - cpu0 << " [sec]  threads=" << nthreads << std::endl;

};
#endif
#endif 
}




