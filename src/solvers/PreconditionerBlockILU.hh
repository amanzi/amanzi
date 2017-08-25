/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! PreconditionerBlockILU:   Incomplete LU preconditioner.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

The internal parameters for block ILU are as follows:

Example:

.. code-block:: xml

  <ParameterList name="block ilu parameters">
    <Parameter name="fact: relax value" type="double" value="1.0"/>
    <Parameter name="fact: absolute threshold" type="double" value="0.0"/>
    <Parameter name="fact: relative threshold" type="double" value="1.0"/>
    <Parameter name="fact: level-of-fill" type="int" value="0"/>
    <Parameter name="overlap" type="int" value="0"/>
    <Parameter name="schwarz: combine mode" type="string" value="Add"/>
    </ParameterList>
  </ParameterList>

*/


#ifndef AMANZI_PRECONDITIONER_BLOCK_ILU_HH_
#define AMANZI_PRECONDITIONER_BLOCK_ILU_HH_

#include "Ifpack.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

class PreconditionerBlockILU : public Preconditioner {
 public:
  PreconditionerBlockILU() {};
  ~PreconditionerBlockILU() {};

  void Init(const std::string& name, const Teuchos::ParameterList& list);
  void Update(const Teuchos::RCP<Epetra_RowMatrix>& A);
  void Destroy();

  int ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv);

  int returned_code() { return returned_code_; }

 private:
  Teuchos::ParameterList list_;
  Teuchos::RCP<Ifpack_Preconditioner> IfpILU_;

  bool initialized_;
  int returned_code_;
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi



#endif
