/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (amklinv@sandia.gov)
           Ethan Coon (coonet@ornl.gov)
*/

//! Direct solvers via Trilinos.
/*!

Amesos library of Trilinos package provides interfaces to a few direct solvers.
List `"amesos parameters`" contains parameters that understood by this library.
These parameters may violate the camel-case convention employed by this spec.
Additional parameters are:

* `"solver name`" [string] declares name of one of the supported direct solvers.
  Available options are `"klu`", `"superludist`", `"basker`", etc, see Amesos and
  Amesos2 manuals for details. The default value is serial solver `"klu`".

* `"amesos version`" [int] specifies version of Amesos. Available options are 1 and 2.
  The default value is 1.

.. code-block:: xml

  <ParameterList name="_AMESOS KLU">  <!-- parent list -->
  <Parameter name="direct method" type="string" value="amesos"/>
  <ParameterList name="amesos parameters">
    <Parameter name="solver name" type="string" value="klu"/>
    <Parameter name="amesos version" type="int" value="1"/>
  </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_AMESOS_OPERATOR_HH_
#define AMANZI_AMESOS_OPERATOR_HH_

#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "Inverse.hh"


class Epetra_LinearProblem;
class Amesos_BaseSolver;

namespace Amanzi {

class VerboseObject;

namespace AmanziSolvers {

/* ******************************************************************
* Auxiliary base class.
****************************************************************** */
class DirectMethodAmesos
  : public Inverse<Epetra_CrsMatrix, Epetra_CrsMatrix, Epetra_Vector, Epetra_Map> {
 private:
  using Inv = Inverse<Epetra_CrsMatrix, Epetra_CrsMatrix, Epetra_Vector, Epetra_Map>;

 public:
  DirectMethodAmesos() : Inv(), inited_(false), updated_(false), computed_(false){};

  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) override;
  virtual void InitializeInverse() override;
  virtual void ComputeInverse() override;
  virtual int ApplyInverse(const Epetra_Vector&, Epetra_Vector&) const override;

  virtual int returned_code() const override { return returned_code_; }
  virtual std::string returned_code_string() const override;

 protected:
  using Inv::m_;
  using Inv::h_;

  mutable Teuchos::ParameterList plist_;
  Teuchos::RCP<VerboseObject> vo_;

  mutable Teuchos::RCP<Epetra_LinearProblem> problem_;
  mutable Teuchos::RCP<Amesos_BaseSolver> solver_;

  std::string solver_name_;
  mutable int returned_code_;

  bool inited_;
  bool updated_;
  bool computed_;
};


} // namespace AmanziSolvers
} // namespace Amanzi

#endif
