/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! Typedefs to make forward declarations and interfaces a bit easier.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!

  Forward declarations of types for use in more generic code.

*/


#ifndef AMANZI_TYPES_HH_
#define AMANZI_TYPES_HH_

#include "Teuchos_RCPDecl.hpp"

#define TRILINOS_TPETRA_STACK 1

#ifdef TRILINOS_TPETRA_STACK

#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Tpetra_Map_fwd.hpp"
#include "Tpetra_Import_fwd.hpp"
#include "Tpetra_Vector_fwd.hpp"
#include "Tpetra_MultiVector_fwd.hpp"

#include "Kokkos_Core.hpp"
#ifdef HAVE_CUDA
using AmanziDefaultDevice = Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>;
#else
using AmanziDefaultDevice = Kokkos::Serial;
#endif

#else

class Epetra_Comm;
class Epetra_MpiComm;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_Vector;
class Epetra_IntVector;
class Epetra_MultiVector;
//class Epetra_MultiIntVector; // defined in trilinos > 12.??

#endif



namespace Amanzi {

#ifdef TRILINOS_TPETRA_STACK

// Tpetra uses Teuchos Comm
typedef Teuchos::Comm<int> Comm_type;
#ifdef HAVE_MPI
typedef Teuchos::MpiComm<int> MpiComm_type;
#endif

// Tpetra maps and importers
typedef Tpetra::Map<> Map_type;
typedef Tpetra::Map<> BlockMap_type; // is there a Tpetra block map?
typedef Tpetra::Import<> Import_type;

// Tpetra vectors
typedef Tpetra::Vector<> Vector_type;
typedef Tpetra::Vector<int> IntVector_type;
typedef Tpetra::MultiVector<> MultiVector_type;
typedef Tpetra::MultiVector<int> IntMultiVector_type;


#else // Trilinos Epetra stack

// Epetra Comm
typedef Epetra_Comm Comm_type;
#ifdef HAVE_MPI
typedef Epetra_MpiComm MpiComm_type;
#endif

// Epetra maps
typedef Epetra_Map Map_type;
typedef Epetra_BlockMap BlockMap_type;
typedef Epetra_Import Import_type;

// Epetra vectors
typedef Epetra_Vector Vector_type;
typedef Epetra_IntVector IntVector_type;
typedef Epetra_MultiVector MultiVector_type;
//typedef Epetra_MultiIntVector IntMultiVector_type; // defined in trilinos > 12.??


#endif

// derived pointer types
typedef Teuchos::RCP<const Comm_type> Comm_ptr_type;
typedef Teuchos::RCP<const Map_type> Map_ptr_type;
typedef Teuchos::RCP<const BlockMap_type> BlockMap_ptr_type;
typedef Teuchos::RCP<const Import_type> Import_ptr_type;

// non-const
typedef Teuchos::RCP<Vector_type> Vector_ptr_type;
typedef Teuchos::RCP<MultiVector_type> MultiVector_ptr_type;
typedef Teuchos::RCP<IntVector_type> IntVector_ptr_type;
//typedef Teuchos::RCP<MultiIntVector_type> MultiIntVector_ptr_type;

// const
typedef Teuchos::RCP<const Vector_type> cVector_ptr_type;
typedef Teuchos::RCP<const MultiVector_type> cMultiVector_ptr_type;
typedef Teuchos::RCP<const IntVector_type> cIntVector_ptr_type;
//typedef Teuchos::RCP<const MultiIntVector_type> cMultiIntVector_ptr_type;


} // namespace Amanzi

#endif
