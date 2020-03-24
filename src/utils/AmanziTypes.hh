/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! Typedefs to make forward declarations and interfaces a bit easier.

/*!

  Forward declarations of types for use in more generic code.

*/


#ifndef AMANZI_TYPES_HH_
#define AMANZI_TYPES_HH_

#include "Teuchos_RCPDecl.hpp"

#define TRILINOS_TPETRA_STACK 1

#ifdef TRILINOS_TPETRA_STACK

#  include "Teuchos_Comm.hpp"
#  include "Teuchos_DefaultMpiComm.hpp"
#  include "Tpetra_Map_fwd.hpp"
#  include "Tpetra_Import_fwd.hpp"
#  include "Tpetra_Vector_fwd.hpp"
#  include "Tpetra_MultiVector_fwd.hpp"
#  include "Tpetra_RowMatrix_fwd.hpp"

#  include "Kokkos_Core.hpp"

// // host execution space
// #ifdef KOKKOS_ENANLE_OPENMP
// using AmanziHostExecutionSpace = Kokkos::OpenMP;
// #elif KOKKOS_ENABLE_THREADS
// using AmanziHostExecutionSpace = Kokkos::Threads;
// #else
// using AmanziHostExecutionSpace = Kokkos::Serial;
// #endif // KOKKOS_HOST

// #ifdef KOKKOS_ENABLE_CUDA
// // Default CUDA using UVM
// using AmanziDefaultDevice = Kokkos::Device<Kokkos::Cuda,
// Kokkos::CudaUVMSpace>; using AmanziDefaultHost =
// Kokkos::Device<AmanziHostExecutionSpace, Kokkos::CudaUVMSpace>;
// //Kokkos::Serial; // ???? #else using AmanziDefaultHost =
// AmanziHostExecutionSpace; #endif // KOKKOS_ENABLE_CUDA

using AmanziDefaultDevice =
  Kokkos::Device<Kokkos::DefaultExecutionSpace,
                 Kokkos::DefaultExecutionSpace::memory_space>;
using AmanziDefaultHost =
  Kokkos::Device<Kokkos::DefaultExecutionSpace,
                 Kokkos::DefaultExecutionSpace::memory_space>;


#else  // TRILINOS_TPETRA_STACK
class Epetra_Comm;
class Epetra_MpiComm;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_Vector;
class Epetra_IntVector;
class Epetra_MultiVector;
// class Epetra_MultiIntVector; // defined in trilinos > 12.??
#endif // TRILINOS_TPETRA_STACK


namespace Amanzi {

#ifdef TRILINOS_TPETRA_STACK

// -- default specializations -- NOTE: this is the ony place here and in data
//    structures where these types are defined!
using double_type = double;
using int_type = int;
using LO = int;
using GO = int;

// Tpetra uses Teuchos Comm
typedef Teuchos::Comm<int> Comm_type;
#  ifdef HAVE_MPI
typedef Teuchos::MpiComm<int> MpiComm_type;
#  endif

// Tpetra maps and importers
typedef Tpetra::Map<LO,GO> Map_type;
typedef Tpetra::Map<LO,GO> BlockMap_type; // is there a Tpetra block map?
typedef Tpetra::Import<LO,GO> Import_type;

// Tpetra vectors
// -- alias
template <typename Scalar>
using Vector_type_ = Tpetra::Vector<Scalar,LO,GO>;

template <typename Scalar>
using MultiVector_type_ = Tpetra::MultiVector<Scalar,LO,GO>;

template <typename Scalar>
using Vector_ptr_type_ = Teuchos::RCP<Vector_type_<Scalar>>;
template <typename Scalar>
using cVector_ptr_type_ = Teuchos::RCP<const Vector_type_<Scalar>>;
template <typename Scalar>
using MultiVector_ptr_type_ = Teuchos::RCP<MultiVector_type_<Scalar>>;
template <typename Scalar>
using cMultiVector_ptr_type_ = Teuchos::RCP<const MultiVector_type_<Scalar>>;

using Vector_type = Vector_type_<double_type>;
using MultiVector_type = MultiVector_type_<double_type>;
using IntVector_type = Vector_type_<int_type>;
using IntMultiVector_type = MultiVector_type_<int_type>;

using Matrix_type = Tpetra::RowMatrix<double_type,LO,GO>;
using Matrix_ptr_type = Teuchos::RCP<Matrix_type>;
using cMatrix_ptr_type = Teuchos::RCP<const Matrix_type>;

// Kokkos Views into vectors
template <class Device_type, typename Scalar>
using VectorView_type_ =
  Kokkos::View<Scalar*, Kokkos::LayoutLeft,
               Device_type>; // MH: layout depends on Tpetra fix later

template <class Device_type, typename Scalar>
using cVectorView_type_ =
  Kokkos::View<const Scalar*, Kokkos::LayoutLeft, Device_type>;

template <class Device_type, typename Scalar>
using MultiVectorView_type_ =
  Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device_type>;

template <class Device_type, typename Scalar>
using cMultiVectorView_type_ =
  Kokkos::View<const Scalar**, Kokkos::LayoutLeft, Device_type>;

// and some defaults
template <class Device_type>
using VectorView_type = VectorView_type_<Device_type, double_type>;
template <class Device_type>
using cVectorView_type = cVectorView_type_<Device_type, double_type>;
template <class Device_type>
using MultiVectorView_type = MultiVectorView_type_<Device_type, double_type>;
template <class Device_type>
using cMultiVectorView_type = cMultiVectorView_type_<Device_type, double_type>;

template <class Device_type>
using IntVectorView_type = VectorView_type_<Device_type, int_type>;
template <class Device_type>
using cIntVectorView_type = cVectorView_type_<Device_type, int_type>;
template <class Device_type>
using IntMultiVectorView_type = MultiVectorView_type_<Device_type, int_type>;
template <class Device_type>
using cIntMultiVectorView_type = cMultiVectorView_type_<Device_type, int_type>;


#else // Trilinos Epetra stack

// Epetra Comm
typedef Epetra_Comm Comm_type;
#  ifdef HAVE_MPI
typedef Epetra_MpiComm MpiComm_type;
#  endif

// Epetra maps
typedef Epetra_Map Map_type;
typedef Epetra_BlockMap BlockMap_type;
typedef Epetra_Import Import_type;

// Epetra vectors
typedef Epetra_Vector Vector_type;
typedef Epetra_IntVector IntVector_type;
typedef Epetra_MultiVector MultiVector_type;
// typedef Epetra_MultiIntVector IntMultiVector_type; // defined in trilinos
// > 12.??


#endif

// derived pointer types
typedef Teuchos::RCP<const Comm_type> Comm_ptr_type;
typedef Teuchos::RCP<const Map_type> Map_ptr_type;
typedef Teuchos::RCP<const BlockMap_type> BlockMap_ptr_type;
typedef Teuchos::RCP<const Import_type> Import_ptr_type;

// non-const
typedef Teuchos::RCP<Vector_type> Vector_ptr_type;
typedef Teuchos::RCP<MultiVector_type> MultiVector_ptr_type;

// const
typedef Teuchos::RCP<const Vector_type> cVector_ptr_type;
typedef Teuchos::RCP<const MultiVector_type> cMultiVector_ptr_type;


} // namespace Amanzi

#endif
