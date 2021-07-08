#include <iostream>
#include <string>

#ifdef _OPENMP
#include "omp.h"
#endif

#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "UnitTest++.h"

#include "exceptions.hh"
#include "IterativeMethodPCG.hh"

#include "Preconditioner.hh"
#include "InverseFactory.hh"

SUITE(SOLVERS) {
const int N = 125;
using namespace Amanzi;
using namespace Amanzi::AmanziSolvers;

inline Teuchos::RCP<Epetra_CrsMatrix> matrix(const Epetra_Map& map) {
  auto A = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map, map, 3));
  for (int i = 0; i < N; i++) {
    int indices[3];
    double values[3] = {double(-i), double(2 * i + 1), double(-i - 1)};
    for (int k = 0; k < 3; k++) indices[k] = i + k - 1; 
    A->InsertMyValues(i, 3, values, indices);
  }
  A->FillComplete(map, map);
  return A;
}
  

inline Teuchos::RCP<Amanzi::AmanziSolvers::Preconditioner>
preconditioner(const std::string& name,
               const Teuchos::RCP<Epetra_CrsMatrix>& mat)
{
  Teuchos::ParameterList plist;
  plist.set<std::string>("preconditioning method", name);
  Teuchos::ParameterList& tmp = plist.sublist(name+" parameters");

  tmp.sublist("verbose object").set<std::string>("verbosity level", "high");
  if (name == "ml") {
    tmp.set("coarse: max size", 5);
    tmp.set("cycle applications", 2);
    tmp.set("ML output", 0);
  } else if (name == "muelu") {
    tmp.set<int>("coarse: max size", 5)
       .set<std::string>("coarse: type", "SuperLU_dist")
       .set<std::string>("verbosity", "low")
       .set<std::string>("multigrid algorithm", "sa")
       .set<std::string>("smoother: type", "RELAXATION").sublist("smoother: params")
         .set<std::string>("relaxation: type", "symmetric Gauss-Seidel")
         .set<int>("relaxation: sweeps", 1)
         .set<double>("relaxation: damping factor", 0.9);
  } else if (name == "boomer amg") {
    tmp.set("max coarse size", 5);
    tmp.set("cycle applications", 1);
  }

  static_assert(Amanzi::AmanziSolvers::Impl::is_assembled<Epetra_CrsMatrix>::value, "Epetra_CrsMatrix is assembled?");
  auto pc = createAssembledMethod<>(name, plist);
  pc->set_matrix(mat);
  return pc;
};    

inline Teuchos::RCP<IterativeMethodPCG<Epetra_CrsMatrix,Amanzi::AmanziSolvers::Preconditioner,Epetra_Vector,Epetra_Map>>
    get_solver(const std::string& name, const Teuchos::RCP<Epetra_CrsMatrix>& m)
{
  auto pc = preconditioner(name, m);
  
  Teuchos::ParameterList plist;
  plist.set("error tolerance", 1.e-12);
  plist.set("maximum number of iterations", 200);

  auto inv = Teuchos::rcp(new IterativeMethodPCG<Epetra_CrsMatrix,Amanzi::AmanziSolvers::Preconditioner,Epetra_Vector,Epetra_Map>());
  inv->set_inverse_parameters(plist);
  inv->set_matrices(m, pc);
  inv->InitializeInverse();
  inv->ComputeInverse();

  return inv;
}


TEST(PRECONDITIONERS) {
  std::cout << "\nComparison of preconditioners for N=125" << std::endl;

  Epetra_MpiComm comm(MPI_COMM_SELF);
  Epetra_Map map(N, 0, comm);
  auto m = matrix(map);

  Epetra_Vector u(map), v(map);
  for (int i = 0; i < N; i++) u[i] = 1.0 / (i + 2.0);
  
  for (const auto& prec_name : {
    "identity",
    "diagonal",
    "block ilu",
    "boomer amg",
    "euclid",
    "ml"
#if defined(HAVE_MUELU_EPETRA)
    , "muelu"
#endif
  }) {
    auto solver = get_solver(prec_name, m);
    v.PutScalar(0.0);

    std::cout << "Preconditioner: " << prec_name << std::endl
              << "-------------------------------------------" << std::endl;
    solver->ApplyInverse(u, v);
    CHECK_CLOSE(11.03249773994628, v[0], 1e-6);
    CHECK_CLOSE(10.53249773994628, v[1], 1e-6);
  }
};

#if 0 
#ifdef _OPENMP
TEST(PRECONDITIONERS_OMP) {
  std::cout << "\nComparison of preconditioners for N=125" << std::endl;

  Epetra_MpiComm comm(MPI_COMM_SELF);
  auto map = Teuchos::rcp(new Epetra_Map(N, 0, comm));

  double cpu0 = omp_get_wtime();

#pragma omp parallel for shared(map) num_threads(2)
  for (const auto& prec_name : {"identity", "diagonal"}) {
    auto m = matrix(*map);

    Epetra_Vector u(*map), v(*map);
    for (int i = 0; i < N; i++) u[i] = 1.0 / (i + 2.0);
  
    auto solver = get_solver(prec_name, m);
    v.PutScalar(0.0);

    solver->ApplyInverse(u, v);
    CHECK_CLOSE(11.03249773994628, v[0], 1e-6);
    CHECK_CLOSE(10.53249773994628, v[1], 1e-6);
  }


  int nthreads = omp_get_max_threads();
  double cpu1 = omp_get_wtime();
  std::cout << "CPU (parallel): " << cpu1 - cpu0 << " [sec]  threads=" << nthreads << std::endl;

  // serial run
  cpu0 = omp_get_wtime();
  for (const auto& prec_name : {"identity", "diagonal"}) {
    auto m = matrix(*map);

    Epetra_Vector u(*map), v(*map);
    for (int i = 0; i < N; i++) u[i] = 1.0 / (i + 2.0);
  
    auto solver = get_solver(prec_name, m);
    v.PutScalar(0.0);

    solver->ApplyInverse(u, v);
    CHECK_CLOSE(11.03249773994628, v[0], 1e-6);
    CHECK_CLOSE(10.53249773994628, v[1], 1e-6);
  }


  nthreads = omp_get_max_threads();
  cpu1 = omp_get_wtime();
  std::cout << "CPU (serial): " << cpu1 - cpu0 << " [sec]  threads=" << nthreads << std::endl;
  
};
#endif
#endif 

}




