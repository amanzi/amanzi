/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/
//! Incomplete LU preconditioner.

/*!

Incomplete LU is an approximate scheme based on partial factorization.  The
implementation here is that provided in the Ifpack package of Trilinos.  This
approach is a block solve that performs the ILU on each MPI process and uses
Additive Schwarz to combine the blocks.

This is provided when using the `"preconditioner type`"=`"block ilu`" in the
`Preconditioner`_ spec.

.. _preconditioner-typed-block-ilu-spec:
.. admonition:: preconditioner-typed-block-ilu-spec:

    * `"fact: relax value`" ``[double]`` **1.0**
    * `"fact: absolute threshold`" ``[double]`` **0.0**
    * `"fact: relative threshold`" ``[double]`` **1.0**
    * `"fact: level-of-fill`" ``[int]`` **0**
    * `"overlap`" ``[int]`` **0** Overlap of the combination.
    * `"schwarz: combine mode`" ``[string]`` **Add** 

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
