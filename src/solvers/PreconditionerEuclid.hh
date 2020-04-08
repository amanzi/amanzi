/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/
//! HYPRE's parallel ILU as a preconditioner.

/*!

Euclid is a Parallel Incomplete LU, provided as part of the HYPRE project
through the Ifpack interface.

This is provided when using the `"preconditioner type`"=`"euclid`" in the
`Preconditioner`_ spec.

.. _preconditioner-typed-euclid-spec:
.. admonition:: preconditioner-typed-euclid-spec:

    * `"ilu(k) fill level`" ``[int]`` **1** The factorization level.
    * `"ilut drop tolerance`" ``[double]`` **0** Defines a drop tolerance relative to the largest absolute value of any entry in the row being factored.
    * `"rescale row`" ``[bool]`` **false** If true, values are scaled prior to factorization so that largest value in any row is +1 or -1. Note that this can destroy matrix symmetry.
    * `"verbosity`" ``[int]`` **0** Prints a summary of runtime settings and timing information to stdout.

*/

#ifndef AMANZI_PRECONDITIONER_EUCLID_HH_
#define AMANZI_PRECONDITIONER_EUCLID_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Ifpack.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

class PreconditionerEuclid : public Preconditioner {
 public:
  PreconditionerEuclid() {};
  ~PreconditionerEuclid() {};

  void Init(const std::string& name, const Teuchos::ParameterList& list);
  void Update(const Teuchos::RCP<Epetra_RowMatrix>& A);
  void Destroy() {};

  int ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv);

  int returned_code() { return returned_code_; }

 private:
  Teuchos::ParameterList plist_;
  std::vector<Teuchos::RCP<FunctionParameter> > funcs_;

  int returned_code_;
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_;

};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi



#endif
