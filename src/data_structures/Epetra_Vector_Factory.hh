//! Helper factory for storing Epetra_Vectors  in State

/*

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatsky (dasvyat@lanl.gov)

  A simple wrapper for creating Epetra_Vectors  in State.
*/

#ifndef AMANZI_EPETRA_VECTOR_FACTORY_HH_
#define AMANZI_EPETRA_VECTOR_FACTORY_HH_

#include "dbc.hh"
#include "errors.hh"
#include "Epetra_SerialComm.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"

namespace Amanzi {
//
// A factory for putting these into state.
// -----------------------------------------------------------------------------
class Epetra_Vector_Factory {
 public:
  Epetra_Vector_Factory() :
      size_(0)
  {};

  int get_size() const { return size_; }
  void set_size(int size) { size_ = size; }
  
  Teuchos::RCP<Epetra_Vector> Create() const {
    AMANZI_ASSERT(size_ > 0);
    Epetra_SerialComm comm;
    Epetra_LocalMap map(size_, 0, comm);
    return Teuchos::rcp(new Epetra_Vector(map));
  }

 private:
  int size_;
};

}  // namespace Amanzi

#endif
