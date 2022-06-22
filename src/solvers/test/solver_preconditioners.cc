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

    for (int i = 0; i < map->getLocalNumElements()-1; i++) {
      int indices[3];
      double values[3] = { double(-i), double(2 * i + 1), double(-i - 1) };
      for (int k = 0; k < 3; k++) indices[k] = i + k - 1;
      A_->insertLocalValues(i, 3, values, indices);
    }

    int i = map->getLocalNumElements()-1;
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

  if(name == "ifpack2: ILUT"){
    tmp.set("fact: ilut level-of-fill", 1.0)
    .set("fact: drop tolerance",0.0); 
  }
  if (name == "ifpack2: RILUK") {
    tmp.set("fact: type", "KSPILUK") // get the kokkos kernels variant
    .set<int>("fact: iluk level-of-fill", 1)
    .set<double>("fact: drop tolerance", 0.0);
  }
#ifdef KOKKOS_ENABLE_CUDA
  if(name == "hypre: boomer amg"){
    tmp.set<int>("verbosity", 1)
    .set<double>("strong threshold", 0.5)
    .set<int>("cycle applications", 2)
    .set<int>("smoother sweeps", 3)
    .set<int>("coarsening type", 8) /* 8: PMIS */
    .set<int>("interpolation type", 6) 
    .set<int>("relaxation order", 0) /* must be false */
    .set<int>("relaxation type", 6); 
  }
#else
  if(name == "hypre: boomer amg"){
    tmp.set<int>("verbosity", 1);
  }
#endif  

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
  std::cout << "\nComparison of preconditioners for N="<<N << std::endl;

  auto comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_SELF));
  auto map = Teuchos::rcp(new Map_type(N, 0, comm));
  auto m = matrix(map);

  Vector_type u(map);
  {
    auto uv = u.getLocalViewHost(Tpetra::Access::ReadWrite);
    for (int i = 0; i < N; i++) uv(i,0) = 1.0 / (i + 2.0);
  }

  Vector_type v(map);

#ifdef KOKKOS_ENABLE_CUDA
  static std::vector<std::string> prec_name = { 
    "identity",
    "diagonal",
    "hypre: boomer amg",
    //"ifpack2: SCHWARZ",
    //"ifpack2: ILUT",
    //"ifpack2: RILUK",
    //"ifpack2: FAST_ILU"
  };
#else 
  static std::vector<std::string> prec_name = { 
    "identity",
    "diagonal",
    "ifpack2: ILUT",
    "hypre: boomer amg",
    "hypre: euclid",
    "muelu"
  };
#endif 

  for (const auto& prec_name : prec_name) {
    auto solver = get_solver(prec_name, m);
    v.putScalar(0.0);

    std::cout << "Preconditioner: " << prec_name << std::endl
              << "-------------------------------------------" << std::endl;
    solver->applyInverse(u, v);

    auto vv = v.getLocalViewHost(Tpetra::Access::ReadOnly);
    CHECK_CLOSE(11.03249773994628, vv(0,0), 1e-6);
    CHECK_CLOSE(10.53249773994628, vv(1,0), 1e-6);
    std::cout << std::endl << std::endl;
  }
};

#if 0 
#ifdef _OPENMP
TEST(PRECONDITIONERS_OMP) {
  std::cout << "\nComparison of preconditioners for N="<<N << std::endl;

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




