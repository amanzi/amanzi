/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! PreconditionerBlockILU:   Incomplete LU preconditioner.

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

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

class PreconditionerBlockILU
  : public Preconditioner<Epetra_RowMatrix, Epetra_MultiVector> {
 public:
  PreconditionerBlockILU(){};
  ~PreconditionerBlockILU(){};

  void
  Init(const std::string& name, const Teuchos::ParameterList& list) override;
  void Update(const Teuchos::RCP<const Epetra_RowMatrix>& A) override;
  void Destroy() override;

  int ApplyInverse(const Epetra_MultiVector& v,
                   Epetra_MultiVector& hv) const override;

  int returned_code() override { return returned_code_; }

 private:
  Teuchos::ParameterList list_;
  Teuchos::RCP<Ifpack_Preconditioner> IfpILU_;

  bool initialized_;
  mutable int returned_code_;
};

} // namespace AmanziPreconditioners
} // namespace Amanzi


#endif
