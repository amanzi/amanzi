/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
      Ethan Coon (coonet@ornl.gov)
*/

//! PreconditionerBoomerAMG: HYPRE's multigrid preconditioner.

/*!
Internal parameters for Boomer AMG include

* `"tolerance`" ``[double]`` if is not zero, the preconditioner is dynamic
  and approximate the inverse matrix with the prescribed tolerance (in
  the energy norm ???).

* `"smoother sweeps`" ``[int]`` **3** defines the number of smoothing loops.
Default is 3.

* `"cycle applications`" ``[int]`` **5** defines the number of V-cycles.

* `"strong threshold`" ``[double]`` **0.5** defines the number of V-cycles.
Default is 5.

* `"relaxation type`" ``[int]`` **6** defines the smoother to be used. Default
is 6 which specifies a symmetric hybrid Gauss-Seidel / Jacobi hybrid method.
TODO: add others!

* `"coarsen type`" ``[int]`` **0** defines the coarsening strategy to be used.
Default is 0 which specifies a Falgout method. TODO: add others!

* `"max multigrid levels`" ``[int]`` optionally defined the maximum number of
multigrid levels.

* `"use block indices`" ``[bool]`` **false** If true, uses the `"systems of
    PDEs`" code with blocks given by the SuperMap, or one per DoF per entity
    type.

* `"number of functions`" ``[int]`` **1** Any value > 1 tells Boomer AMG to
  use the `"systems of PDEs`" code with strided block type.  Note that, to use
  this approach, unknowns must be ordered with DoF fastest varying (i.e. not
  the native Epetra_MultiVector order).  By default, it uses the `"unknown`"
  approach in which each equation is coarsened and interpolated independently.

* `"nodal strength of connection norm`" ``[int]`` tells AMG to coarsen such
    that each variable has the same coarse grid - sometimes this is more
    "physical" for a particular problem. The value chosen here for nodal
    determines how strength of connection is determined between the
    coupled system.  I suggest setting nodal = 1, which uses a Frobenius
    norm.  This does NOT tell AMG to use nodal relaxation.
    Default is 0.

* `"verbosity`" ``[int]`` **0** prints a summary of run time settings and
  timing information to stdout.  `"1`" prints coarsening info, `"2`" prints
  smoothing info, and `"3`'" prints both.

Example:

.. code-block:: xml

  <ParameterList name="boomer amg parameters">
    <Parameter name="tolerance" type="double" value="0.0"/>
    <Parameter name="smoother sweeps" type="int" value="3"/>
    <Parameter name="cycle applications" type="int" value="5"/>
    <Parameter name="strong threshold" type="double" value="0.5"/>
    <Parameter name="coarsen type" type="int" value="0"/>
    <Parameter name="relaxation type" type="int" value="3"/>
    <Parameter name="verbosity" type="int" value="0"/>
    <Parameter name="number of functions" type="int" value="1"/>
  </ParameterList>

*/

#ifndef AMANZI_PRECONDITIONER_BOOMERAMG_HH_
#define AMANZI_PRECONDITIONER_BOOMERAMG_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Ifpack.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

class PreconditionerBoomerAMG
  : public Preconditioner<Epetra_RowMatrix, Epetra_MultiVector> {
 public:
  PreconditionerBoomerAMG()
    : num_blocks_(0), block_indices_(Teuchos::null), IfpHypre_(Teuchos::null)
  {}
  ~PreconditionerBoomerAMG(){};

  void
  Init(const std::string& name, const Teuchos::ParameterList& list) override;
  void Update(const Teuchos::RCP<const Epetra_RowMatrix>& A) override;
  void Destroy() override{};

  int ApplyInverse(const Epetra_MultiVector& v,
                   Epetra_MultiVector& hv) const override;

  int returned_code() override { return returned_code_; }

 private:
  Teuchos::ParameterList plist_;
  std::vector<Teuchos::RCP<FunctionParameter>> funcs_;
  Teuchos::RCP<std::vector<int>> block_indices_;
  int num_blocks_;
  int block_index_function_index_;

  mutable int returned_code_;
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_;
};

} // namespace AmanziPreconditioners
} // namespace Amanzi


#endif
