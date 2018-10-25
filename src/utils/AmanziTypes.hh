/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Type-defs for use throughout Amanzi

#ifndef AMANZI_TYPES_HH_
#define AMANZI_TYPES_HH_

#include "Teuchos_RCPDecl.hpp"
#include "Tpetra_Map_fwd.hpp"
#include "Tpetra_Import_fwd.hpp"
#include "Tpetra_Vector_fwd.hpp"
#include "Tpetra_MultiVector_fwd.hpp"
#include "Tpetra_CrsGraph_fwd.hpp"
#include "Tpetra_CrsMatrix_fwd.hpp"

namespace Teuchos {
template<typename Ordinal>
class Comm;
}


namespace Amanzi {

// data structure typedefs
typedef Teuchos::Comm<int> Comm_type;
typedef Teuchos::RCP<const Comm_type> Comm_ptr_type;

typedef Tpetra::Map<> Map_type;
typedef Teuchos::RCP<const Map_type> Map_ptr_type;

typedef Tpetra::Import<> Import_type;
typedef Teuchos::RCP<const Import_type> Import_ptr_type;

typedef Tpetra::Vector<> Vector_type;
typedef Teuchos::RCP<const Vector_type> Vector_ptr_type;

typedef Tpetra::Vector<int> IntVector_type;
typedef Teuchos::RCP<const IntVector_type> IntVector_ptr_type;

typedef Tpetra::MultiVector<> MultiVector_type;
typedef Teuchos::RCP<MultiVector_type> MultiVector_ptr_type;

typedef Tpetra::MultiVector<int> IntMultiVector_type;
typedef Teuchos::RCP<IntMultiVector_type> IntMultiVector_ptr_type;

typedef Tpetra::CrsGraph<> CrsGraph_type;
typedef Teuchos::RCP<CrsGraph_type> CrsGraph_ptr_type;

typedef Tpetra::CrsMatrix<> CrsMatrix_type;
typedef Teuchos::RCP<CrsMatrix_type> CrsMatrix_ptr_type;

} // namespace Amanzi

#endif
