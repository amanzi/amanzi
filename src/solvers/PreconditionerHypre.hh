/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (coonet@ornl.gov)
*/

//! Hypre based preconditioners include Algebraic MultiGrid and global ILU
/*!

Boomer AMG is a HYPRE product consisting of a variety of Algebraic Multigrid
methods.  It is accessed through Ifpack.

This is provided when using the `"preconditioning method`"=`"boomer amg`" or
`"preconditioning method`" = `"hypre: boomer amg`" in the `Preconditioner`_
spec.

.. _preconditioner-boomer-amg-spec:
.. admonition:: preconditioner-boomer-amg-spec:

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


ILU is a Parallel Incomplete LU, provided as part of the HYPRE project
through the Ifpack interface.

This is provided when using the `"preconditioning method`"=`"ILU`" or
=`"hypre: ILU`" in the `Preconditioner`_ spec.

.. _preconditioner-ILU-spec:
.. admonition:: preconditioner-ILU-spec:

    * `"ilu(k) fill level`" ``[int]`` **1** The factorization level.
    * `"ilut drop tolerance`" ``[double]`` **0** Defines a drop tolerance relative to the largest absolute value of any entry in the row being factored.
    * `"verbosity`" ``[int]`` **0** Prints a summary of runtime settings and timing information to stdout.


*/

#ifndef AMANZI_PRECONDITIONER_BOOMERAMG_HH_
#define AMANZI_PRECONDITIONER_BOOMERAMG_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Ifpack_Hypre.h"

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "AmanziTypes.hh"

#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "krylov.h"
#include "_hypre_parcsr_mv.h"
#include "_hypre_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {

class VerboseObject;

namespace AmanziSolvers {

class PreconditionerHypre : public AmanziSolvers::Preconditioner {
  enum HyprePreconditioners { Boomer, ILU, MGR, AMS };

 public:
  PreconditionerHypre()
    : AmanziSolvers::Preconditioner(),
      block_indices_(Teuchos::null),
      num_blocks_(0),
      returned_code_(0)
  {}

  ~PreconditionerHypre()
  {
    HYPRE_IJVectorDestroy(XHypre_);
    HYPRE_IJVectorDestroy(YHypre_);
    HYPRE_IJMatrixDestroy(HypreA_);
    if (method_type_ == Boomer) {
      HYPRE_BoomerAMGDestroy(method_);
    } else if (method_type_ == ILU) {
      HYPRE_ILUDestroy(method_);
    } else if (method_type_ == AMS) {
      HYPRE_AMSDestroy(method_);
      if(plist_.isParameter("graph coordinates") &&
         plist_.isType<Teuchos::RCP<Epetra_MultiVector>>("graph coordinates")) {
        HYPRE_IJVectorDestroy(xHypre_);
        HYPRE_IJVectorDestroy(yHypre_);
        HYPRE_IJVectorDestroy(zHypre_);
      }
      if (plist_.isParameter("discrete gradient operator") &&
          plist_.isType<Teuchos::RCP<Epetra_CrsMatrix>>("discrete gradient operator")) {
        HYPRE_IJMatrixDestroy(HypreG_);
      }
    } else if (method_type_ == MGR) {
      HYPRE_MGRDestroy(method_);
    }
  }

  virtual void set_inverse_parameters(Teuchos::ParameterList& list) override final;
  virtual void InitializeInverse() override final;
  virtual void ComputeInverse() override final;
  virtual int ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const override final;

  virtual int returned_code() const override final { return returned_code_; }
  virtual std::string returned_code_string() const override final { return "success"; }

 private:
  void copy_matrix_();


  void InitBoomer_();
  void InitILU_();
  void InitAMS_();
  void Init_(){};
  void InitMGR_();

  int SetCoordinates_(Teuchos::RCP<Epetra_MultiVector>);
  int SetDiscreteGradient_(Teuchos::RCP<const Epetra_CrsMatrix>);
  Teuchos::RCP<const Epetra_Map>
  MakeContiguousColumnMap_(Teuchos::RCP<const Epetra_RowMatrix>&) const;

  Teuchos::ParameterList plist_;
  Teuchos::RCP<VerboseObject> vo_;

  HYPRE_Solver method_;
  Teuchos::RCP<std::vector<int>> block_indices_;
  int num_blocks_;

  mutable int returned_code_;
  Teuchos::RCP<Epetra_RowMatrix> A_;

  Teuchos::RCP<const Map_type> GloballyContiguousRowMap_;
  Teuchos::RCP<const Map_type> GloballyContiguousColMap_;

  Teuchos::RCP<const Map_type> GloballyContiguousNodeRowMap_;
  Teuchos::RCP<const Map_type> GloballyContiguousNodeColMap_;

  HYPRE_ParCSRMatrix ParMatrix_;
  HYPRE_IJMatrix HypreA_;
  HYPRE_ParVector ParX_;
  HYPRE_ParVector ParY_;
  HYPRE_IJVector XHypre_;
  HYPRE_IJVector YHypre_;

  HYPRE_ParCSRMatrix ParMatrixG_;
  HYPRE_IJMatrix HypreG_;

  HYPRE_IJVector xHypre_;
  HYPRE_IJVector yHypre_;
  HYPRE_IJVector zHypre_;
  HYPRE_ParVector xPar_;
  HYPRE_ParVector yPar_;
  HYPRE_ParVector zPar_;

  Teuchos::RCP<hypre_ParVector> XVec_;
  Teuchos::RCP<hypre_ParVector> YVec_;

  static bool inited;
  HyprePreconditioners method_type_ = Boomer;

  Teuchos::RCP<const Epetra_CrsMatrix> G_;
};

} // namespace AmanziSolvers
} // namespace Amanzi


#endif
