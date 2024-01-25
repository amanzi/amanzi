/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Typedefs to make forward declarations and interfaces a bit easier.
/*!

  Forward declarations of types for use in more generic code.

*/

#pragma once

#include <memory>
#include "Teuchos_RCPDecl.hpp"
#include "Kokkos_Core.hpp"

#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Tpetra_Map_fwd.hpp"
#include "Tpetra_Import_fwd.hpp"
#include "Tpetra_Export_fwd.hpp"
#include "Tpetra_Vector_fwd.hpp"
#include "Tpetra_MultiVector_fwd.hpp"
#include "Tpetra_CrsGraph_fwd.hpp"
#include "Tpetra_CrsMatrix_fwd.hpp"

#include "Kokkos_Core.hpp"

// #define OUTPUT_CUDA
// #undef OUTPUT_CUDA

namespace Amanzi {

using DefaultExecutionSpace = Kokkos::DefaultExecutionSpace;
using DefaultHostExecutionSpace = Kokkos::DefaultHostExecutionSpace;

using DefaultMemorySpace = DefaultExecutionSpace::memory_space;
using DefaultHostMemorySpace = DefaultHostExecutionSpace::memory_space;

using HostSpace = Kokkos::HostSpace;
using DefaultDevice = Kokkos::Device<DefaultExecutionSpace, DefaultMemorySpace>;
using DefaultHost = Kokkos::Device<DefaultHostExecutionSpace, DefaultHostMemorySpace>;


// KOKKOS UVM
#ifdef KOKKOS_ENABLE_CUDA_UVM
#  define NOTSERIAL
// If we are using UVM, we sometimes prefer to turn it off to avoid unnecessary
// syncing when we want to stay fully on the device and never come back to the
// host.
using DeviceOnlyMemorySpace = Kokkos::CudaSpace;
using HostSpaceSpecial = Kokkos::CudaUVMSpace;
using DeviceSpaceSpecial = Kokkos::CudaUVMSpace;

using DeviceSpecial =
  Kokkos::Device<Kokkos::HostSpace::execution_space, Kokkos::Cuda::memory_space>;

// If we are using UVM, dual views are not synced to DefaultHostMemorySpace,
// but to CudaUVMSpace even on the host_mirror
using MirrorHost = Kokkos::Device<DefaultExecutionSpace, Kokkos::CudaUVMSpace>;
//using NT = Kokkos::Compat::KokkosCudaWrapperNode;
using Layout = Kokkos::LayoutLeft;
#endif

// KOKKOS CUDA
#ifdef KOKKOS_ENABLE_CUDA
#  ifndef KOKKOS_ENABLE_CUDA_UVM
#    define NOTSERIAL

using DeviceOnlyMemorySpace = Kokkos::CudaSpace;
// using MirrorHost =
//   Kokkos::Device<Kokkos::DefaultHostExecutionSpace, DefaultMemorySpace>;
using MirrorHost = Kokkos::Device<DefaultExecutionSpace, DefaultHostMemorySpace>;
// using MirrorHost =
//   Kokkos::Device<Kokkos::DefaultHostExecutionSpace, DefaultHostMemorySpace>;
// using MirrorHost =
//   Kokkos::Device<Kokkos::DefaultExecutionSpace, DefaultMemorySpace>;
//using NT = Kokkos::Compat::KokkosCudaWrapperNode;
using Layout = Kokkos::LayoutLeft;

using HostSpaceSpecial = Kokkos::HostSpace;
using DeviceSpaceSpecial = Amanzi::DeviceOnlyMemorySpace;


using DeviceSpecial =
  Kokkos::Device<Kokkos::HostSpace::execution_space, Kokkos::HostSpace::memory_space>;


#  endif

// KOKKOS OPENMP
#ifdef KOKKOS_ENABLE_OPENMP
#  define NOTSERIAL
using DeviceOnlyMemorySpace = Kokkos::OMPSpace;
// using MirrorHost =
//   Kokkos::Device<Kokkos::DefaultHostExecutionSpace, DefaultMemorySpace>;
using MirrorHost = Kokkos::Device<DefaultExecutionSpace, DefaultHostMemorySpace>;
// using MirrorHost =
//   Kokkos::Device<Kokkos::DefaultHostExecutionSpace, DefaultHostMemorySpace>;
// using MirrorHost =
//   Kokkos::Device<Kokkos::DefaultExecutionSpace, DefaultMemorySpace>;
//using NT = Kokkos::Compat::KokkosCudaWrapperNode;
using Layout = Kokkos::LayoutRight

  using HostSpaceSpecial = Kokkos::HostSpace;
using DeviceSpaceSpecial = Amanzi::DeviceOnlyMemorySpace;
using DeviceSpecial =
  Kokkos::Device<Kokkos::HostSpace::execution_space, Kokkos::Cuda::memory_space>;

#  endif
#endif

// KOKKOS SERIAL
#ifndef NOTSERIAL

using HostSpaceSpecial = Kokkos::HostSpace;
using DeviceSpaceSpecial = Kokkos::HostSpace;
//using NT = Kokkos::Compat::KokkosSerialWrapperNode;
using DeviceOnlyMemorySpace = DefaultExecutionSpace::memory_space;
using MirrorHost = DefaultHost;
using Layout = Kokkos::LayoutRight;
using DeviceSpecial =
  Kokkos::Device<Kokkos::HostSpace::execution_space, Kokkos::HostSpace::memory_space>;

#endif


enum class MemSpace_kind { HOST, DEVICE };

template <MemSpace_kind MEM>
using MemoryLocation = std::
  conditional_t<MEM == MemSpace_kind::DEVICE, Kokkos::DefaultExecutionSpace, Kokkos::HostSpace>;


} // namespace Amanzi

namespace Teuchos {
template <typename T>
inline Teuchos::RCP<T>
rcp(std::unique_ptr<T>&& in)
{
  return Teuchos::rcp(in.release());
}
} // namespace Teuchos

namespace Amanzi {


// -- default specializations -- NOTE: this is the ony place here and in data
//    structures where these types are defined!
using double_type = double;
using int_type = int;
using LO = int;
using GO = int;

// Tpetra uses Teuchos Comm
typedef Teuchos::Comm<int> Comm_type;
typedef Teuchos::MpiComm<int> MpiComm_type;

using ParameterList_ptr_type = Teuchos::RCP<Teuchos::ParameterList>;

// Tpetra maps and importers
using Map_type = Tpetra::Map<LO, GO>;
using BlockMap_type = Tpetra::Map<LO, GO>;
using Import_type = Tpetra::Import<LO, GO>;
using Export_type = Tpetra::Export<LO, GO>;

// Tpetra vectors
// -- alias
template <typename Scalar>
using Vector_type_ = Tpetra::Vector<Scalar, LO, GO>;

template <typename Scalar>
using MultiVector_type_ = Tpetra::MultiVector<Scalar, LO, GO>;

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

// Kokkos Views into vectors
template <class Device_type, typename Scalar>
using VectorView_type_ = Kokkos::View<Scalar*,
                                      Kokkos::LayoutLeft,
                                      Device_type>; // MH: layout depends on Tpetra fix later

template <class Device_type, typename Scalar>
using cVectorView_type_ = Kokkos::View<const Scalar*, Kokkos::LayoutLeft, Device_type>;

template <class Device_type, typename Scalar>
using MultiVectorView_type_ = Kokkos::View<Scalar**, Kokkos::LayoutLeft, Device_type>;

template <class Device_type, typename Scalar>
using cMultiVectorView_type_ = Kokkos::View<const Scalar**, Kokkos::LayoutLeft, Device_type>;

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

// Graphs and Matrices
using Graph_type = Tpetra::CrsGraph<LO, GO>;
using Matrix_type = Tpetra::CrsMatrix<double_type, LO, GO>;

// derived pointer types
typedef Teuchos::RCP<const Comm_type> Comm_ptr_type;
typedef Teuchos::RCP<const Map_type> Map_ptr_type;
typedef Teuchos::RCP<const BlockMap_type> BlockMap_ptr_type;
using Import_ptr_type = Teuchos::RCP<const Import_type>;
using Export_ptr_type = Teuchos::RCP<const Export_type>;

// non-const
typedef Teuchos::RCP<Vector_type> Vector_ptr_type;
typedef Teuchos::RCP<MultiVector_type> MultiVector_ptr_type;

// const
typedef Teuchos::RCP<const Vector_type> cVector_ptr_type;
typedef Teuchos::RCP<const MultiVector_type> cMultiVector_ptr_type;

using Graph_ptr_type = Teuchos::RCP<Graph_type>;
using Matrix_ptr_type = Teuchos::RCP<Matrix_type>;
using cMatrix_ptr_type = Teuchos::RCP<const Matrix_type>;

} // namespace Amanzi


