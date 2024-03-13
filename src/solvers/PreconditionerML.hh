/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! Trilinos ML smoothed aggregation multigrid.
/*!

This is provided when using the `"preconditioning method`"=`"ml`" in the
`Preconditioner`_ spec.

.. warning:: no input spec defined

See also: https://trilinos.github.io/pdfs/mlguide5.pdf

Example:

.. code-block:: xml

   <ParameterList name="ml parameters">
     <Parameter name="ML output" type="int" value="0"/>
     <Parameter name="aggregation: damping factor" type="double" value="1.33"/>
     <Parameter name="aggregation: nodes per aggregate" type="int" value="3"/>
     <Parameter name="aggregation: threshold" type="double" value="0.0"/>
     <Parameter name="aggregation: type" type="string" value="Uncoupled"/>
     <Parameter name="coarse: type" type="string" value="Amesos-KLU"/>
     <Parameter name="coarse: max size" type="int" value="128"/>
     <Parameter name="coarse: damping factor" type="double" value="1.0"/>
     <Parameter name="cycle applications" type="int" value="2"/>
     <Parameter name="eigen-analysis: iterations" type="int" value="10"/>
     <Parameter name="eigen-analysis: type" type="string" value="cg"/>
     <Parameter name="max levels" type="int" value="40"/>
     <Parameter name="prec type" type="string" value="MGW"/>
     <Parameter name="smoother: damping factor" type="double" value="1.0"/>
     <Parameter name="smoother: pre or post" type="string" value="both"/>
     <Parameter name="smoother: sweeps" type="int" value="2"/>
     <Parameter name="smoother: type" type="string" value="Gauss-Seidel"/>
   </ParameterList>

 */

#ifndef AMANZI_PRECONDITIONER_ML_HH_
#define AMANZI_PRECONDITIONER_ML_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "ml_MultiLevelPreconditioner.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziSolvers {

class PreconditionerML : public Preconditioner {
 public:
  PreconditionerML() : Preconditioner(), initialized_(false){};

  virtual ~PreconditionerML()
  {
    // unclear whether this is needed or not...  It seems that it
    // ought to be, but destructors crash occassionally if it is used.
    //    if (ML_.get()) ML_->DestroyPreconditioner();
  }

  virtual void set_matrices(const Teuchos::RCP<Epetra_CrsMatrix>& m,
                            const Teuchos::RCP<Epetra_CrsMatrix>& h) override final;

  virtual void set_inverse_parameters(Teuchos::ParameterList& list) override final;
  virtual void InitializeInverse() override final;
  virtual void ComputeInverse() override final;
  virtual int ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const override final;

  virtual int returned_code() const override final { return returned_code_; }
  virtual std::string returned_code_string() const override final
  {
    if (returned_code_ == 0) return "success";
    return "PreconditionerML: unknown error";
  }

 private:
  Teuchos::ParameterList list_;
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> ML_;

  mutable int returned_code_;
  bool initialized_;
};

} // namespace AmanziSolvers
} // namespace Amanzi


#endif
