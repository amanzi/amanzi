/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Typedefs to make forward declarations and interfaces a bit easier.
/*!

  Forward declarations of types for use in more generic code.

*/


#ifndef AMANZI_TYPES_HH_
#define AMANZI_TYPES_HH_

#include <memory>
#include "Teuchos_RCPDecl.hpp"
#include "Kokkos_Core.hpp"

#ifdef TRILINOS_TPETRA_STACK
#  include "Teuchos_Comm.hpp"
#  include "Teuchos_MpiComm.hpp"
#  include "Teuchos_SerialComm.hpp"
#else

#  include "Epetra_MpiComm.h"
#  include "Epetra_SerialComm.h"
#  include "Epetra_Map.h"
#  include "Epetra_Import.h"
#endif

namespace Teuchos {
template <typename T>
inline Teuchos::RCP<T>
rcp(std::unique_ptr<T>&& in)
{
  return Teuchos::rcp(in.release());
}
} // namespace Teuchos

namespace Amanzi {

using DefaultExecutionSpace = Kokkos::DefaultExecutionSpace;
using DefaultHostExecutionSpace = Kokkos::DefaultHostExecutionSpace;

using DefaultMemorySpace = DefaultExecutionSpace::memory_space;
using DefaultHostMemorySpace = DefaultHostExecutionSpace::memory_space;

using HostSpace = Kokkos::HostSpace;
using DefaultDevice = Kokkos::Device<DefaultExecutionSpace, DefaultMemorySpace>;
using DefaultHost = Kokkos::Device<DefaultHostExecutionSpace, DefaultHostMemorySpace>;


#ifdef TRILINOS_TPETRA_STACK
typedef Teuchos::Comm<int> Comm_type;
#  ifdef HAVE_MPI
typedef Teuchos::MpiComm<int> MpiComm_type;
typedef Teuchos::SerialComm<int> SerialComm_type;
#  endif
#else
typedef Epetra_Comm Comm_type;
#  ifdef HAVE_MPI
typedef Epetra_MpiComm MpiComm_type;
typedef Epetra_SerialComm SerialComm_type;
#  endif
#endif

typedef Teuchos::RCP<const Comm_type> Comm_ptr_type;
using Map_type = Epetra_Map;
using Map_ptr_type = Teuchos::RCP<Map_type>;
using Import_type = Epetra_Import;
using Import_ptr_type = Teuchos::RCP<Import_type>;


enum class MemSpace_kind { HOST, DEVICE };

template <MemSpace_kind MEM>
using MemoryLocation = std::
  conditional_t<MEM == MemSpace_kind::DEVICE, Kokkos::DefaultExecutionSpace, Kokkos::HostSpace>;


} // namespace Amanzi

#endif
