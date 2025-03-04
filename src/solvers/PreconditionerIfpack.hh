/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! Ifpack suite of preconditioners, including block ILU.
/*!

The Ifpack (Incomplete Factorization Package) from Trilinos provides
additive-Schwarz-based incomplete factorization methods, including ILU and many
others.

The method is provided in the parent list, e.g. `"preconditioning method`" =
`"ifpack: METHOD`" when calling the Preconditioner factory.

Valid methods include:

- `"point relaxation`" : Additive Schwarz + relaxation
- `"block relaxation`" : Additive Schwarz + block relaxation
- `"Amesos`" : Additive Schwarz with Amesos on block
- `"IC`" : Additive Schwarz + Incomplete Cholesky (symmetry required)
- `"ICT`" : Additive Schwarz + tolerance based Incomplete Cholesky (symmetry required)
- `"ILU`" : Additive Schwarz + Incomplete LU
- `"ILUT`" : Additive Schwarz + tolerance based Incomplete LU
- `"block ilu`" : is the same as `"ILU`", as this is the legacy Amanzi name.

Note that all of these can be used without Additive Schwarz by appending
"stand-alone".

The full list of relevant parameters is somewhat method-dependent, and is
documented extensively in the `Trilinos Documentation
<https://docs.trilinos.org/dev/packages/ifpack/doc/html/index.html>`_.

Here we document a subset of the most frequently used parameters -- advanced
users should read the Ifpack User Guide above to see all options.

.. _preconditioner-ifpack-ilu-spec:
.. admonition:: preconditioner-ifpack-ilu-spec:

    * `"schwarz: combine mode`" ``[string]`` **Add** Note that `"Zero`" may
      perform better for nonsymmetric cases.
    * `"overlap`" ``[int]`` **0** overlap of the Additive Schwarz
    * `"fact: relax value`" ``[double]`` **0.0** If nonzero, dropped values are added to the diagonal (times this factor).
    * `"fact: absolute threshold`" ``[double]`` **0.0** Defines the value to
      add to each diagonal element (times the sign of the actual diagonal
      element).
    * `"fact: relative threshold`" ``[double]`` **1.0** Multiplies the
      diagonal by this value before checking the threshold.
    * `"fact: level-of-fill`" ``[int]`` **0**

.. _preconditioner-ifpack-relaxation-spec:
.. admonition:: preconditioner-ifpack-relaxation-spec:

    * `"schwarz: combine mode`" ``[string]`` **Add** Note that `"Zero`" may
      perform better for nonsymmetric cases.
    * `"overlap`" ``[int]`` **0** overlap of the Additive Schwarz
    * `"relaxation: type`" ``[string]`` **Jacobi**
    * `"relaxation: sweeps`" ``[int]`` **1**
    * `"relaxation: damping factor`" ``[double]`` **1.0**

.. _preconditioner-ifpack-amesos-spec:
.. admonition:: preconditioner-amesos-relaxation-spec:

    * `"schwarz: combine mode`" ``[string]`` **Add** Note that `"Zero`" may
      perform better for nonsymmetric cases.
    * `"overlap`" ``[int]`` **0** overlap of the Additive Schwarz
    * `"amesos: solver type`" ``[string]`` **Amesos_Klu**

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


#ifndef AMANZI_PRECONDITIONER_IFPACK_HH_
#define AMANZI_PRECONDITIONER_IFPACK_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Preconditioner.hh"

class Ifpack_Preconditioner;

namespace Amanzi {

class VerboseObject;

namespace AmanziSolvers {

class PreconditionerIfpack : public Preconditioner {
 public:
  PreconditionerIfpack() : Preconditioner(), initialized_(false){};

  virtual void set_inverse_parameters(Teuchos::ParameterList& list) override final;
  virtual void InitializeInverse() override final;
  virtual void ComputeInverse() override final;

  virtual int ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const override final;

  virtual int returned_code() const override final { return returned_code_; }
  virtual std::string returned_code_string() const override final;

 private:
  Teuchos::ParameterList plist_;
  Teuchos::RCP<Ifpack_Preconditioner> IfpILU_;
  Teuchos::RCP<VerboseObject> vo_;

  bool initialized_;
  mutable int returned_code_;
};

} // namespace AmanziSolvers
} // namespace Amanzi


#endif
