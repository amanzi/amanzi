/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (coonet@ornl.gov)
*/

//! Hypre based preconditioners include Algebraic MultiGrid and global ILU
/*!

Boomer AMG is a HYPRE product consisting of a variety of Algebraic Multigrid
methods.  It is accessed through Ifpack.

This is provided when using the `"preconditioning method`"=`"boomer amg`" or
`"preconditioning method`" = `"hypre: boomer amg`" in the `Preconditioner`_
spec.

.. _preconditioner-typed-boomer-amg-spec:
.. admonition:: preconditioner-typed-boomer-amg-spec:

    * `"tolerance`" ``[double]`` **0.** If is not zero, the preconditioner is dynamic 
      and approximate the inverse matrix with the prescribed tolerance (in
      the energy norm ???).

    * `"smoother sweeps`" ``[int]`` **3** defines the number of smoothing loops. Default is 3.

    * `"cycle applications`" ``[int]`` **5** defines the number of V-cycles.

    * `"strong threshold`" ``[double]`` **0.5** defines the number of V-cycles. Default is 5.

    * `"relaxation type`" ``[int]`` **6** defines the smoother to be used. Default is 6 
      which specifies a symmetric hybrid Gauss-Seidel / Jacobi hybrid method. TODO: add others!

    * `"coarsen type`" ``[int]`` **0** defines the coarsening strategy to be used. Default is 0 
      which specifies a Falgout method. TODO: add others!

    * `"max multigrid levels`" ``[int]`` optionally defined the maximum number of multigrid levels.

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
      determines how strength of connection is determined between the coupled
      system.  I suggest setting nodal = 1, which uses a Frobenius norm.  This
      does NOT tell AMG to use nodal relaxation.  Default is 0.

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


Euclid is a Parallel Incomplete LU, provided as part of the HYPRE project
through the Ifpack interface.

This is provided when using the `"preconditioning method`"=`"euclid`" or
=`"hypre: euclid`" in the `Preconditioner`_ spec.

.. _preconditioner-typed-euclid-spec:
.. admonition:: preconditioner-typed-euclid-spec:

    * `"ilu(k) fill level`" ``[int]`` **1** The factorization level.
    * `"ilut drop tolerance`" ``[double]`` **0** Defines a drop tolerance relative to the largest absolute value of any entry in the row being factored.
    * `"rescale row`" ``[bool]`` **false** If true, values are scaled prior to factorization so that largest value in any row is +1 or -1. Note that this can destroy matrix symmetry.
    * `"verbosity`" ``[int]`` **0** Prints a summary of runtime settings and timing information to stdout.

  
*/

#ifndef AMANZI_PRECONDITIONER_BOOMERAMG_HH_
#define AMANZI_PRECONDITIONER_BOOMERAMG_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Ifpack_Hypre.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {

class VerboseObject;

namespace AmanziSolvers {

class PreconditionerHypre : public Amanzi::AmanziSolvers::Preconditioner {
 public:
  PreconditionerHypre() :
      Amanzi::AmanziSolvers::Preconditioner(),
      num_blocks_(0),
      block_indices_(Teuchos::null),
      IfpHypre_(Teuchos::null),
      returned_code_(0)
  {}

  virtual void InitializeInverse(Teuchos::ParameterList& list) override final;
  virtual void UpdateInverse() override final;
  virtual void ComputeInverse() override final;
  virtual int ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const override final;

  virtual int returned_code() const override final { return returned_code_; }
  virtual std::string returned_code_string() const override final {
    if (returned_code_ == 0) return "not yet applied.";
    return "success";
  }
  

 private:

  void Init_();
  void InitBoomer_();
  void InitEuclid_();
  
  Teuchos::ParameterList plist_;
  Teuchos::RCP<VerboseObject> vo_;
  
  Hypre_Solver method_;
  std::vector<Teuchos::RCP<FunctionParameter> > funcs_;
  Teuchos::RCP<std::vector<int> > block_indices_;
  int num_blocks_;
  int block_index_function_index_;
  
  mutable int returned_code_;
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_;
  Teuchos::RCP<Epetra_RowMatrix> A_;
};

}  // namespace AmanziSolvers
}  // namespace Amanzi



#endif
