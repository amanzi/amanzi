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


Euclid is a Parallel Incomplete LU, provided as part of the HYPRE project
through the Ifpack interface.

This is provided when using the `"preconditioning method`"=`"euclid`" or
=`"hypre: euclid`" in the `Preconditioner`_ spec.

.. _preconditioner-euclid-spec:
.. admonition:: preconditioner-euclid-spec:

    * `"ilu(k) fill level`" ``[int]`` **1** The factorization level.
    * `"ilut drop tolerance`" ``[double]`` **0** Defines a drop tolerance relative to the largest absolute value of any entry in the row being factored.
    * `"rescale row`" ``[bool]`` **false** If true, values are scaled prior to factorization so that largest value in any row is +1 or -1. Note that this can destroy matrix symmetry.
    * `"verbosity`" ``[int]`` **0** Prints a summary of runtime settings and timing information to stdout.


*/

#ifndef AMANZI_PRECONDITIONER_BOOMERAMG_HH_
#define AMANZI_PRECONDITIONER_BOOMERAMG_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Tpetra_RowMatrix_decl.hpp"

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ArrayRCP.hpp"

#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "krylov.h"
#include "_hypre_parcsr_mv.h"
#include "_hypre_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE.h"

#include "cuda_decl.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

#define HAVE_IFPACK2_HYPRE

namespace Amanzi {

class VerboseObject;

namespace AmanziSolvers {

class PreconditionerHypre : public Preconditioner {

  enum HyprePreconditioners {Boomer, Euclid};
  
  using RowMatrix_type = Tpetra::RowMatrix<double_type,LO,GO>;

 public:

  ~PreconditionerHypre(){
    HYPRE_IJVectorDestroy(XHypre_);
    HYPRE_IJVectorDestroy(YHypre_);
    if(PrecondType == Boomer){
      HYPRE_BoomerAMGDestroy(HyprePrecond_);
    } else if(PrecondType == Euclid){
      HYPRE_EuclidDestroy(HyprePrecond_); 
    }
  }

  static void init(){
    nvtxRangePush("HP: init");
    if(!inited){
      HYPRE_Init(); 
      HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE);
      HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);
      HYPRE_SetSpGemmUseCusparse(true);
      HYPRE_SetUseGpuRand(true);
      //if (useHypreGpuMemPool)
      //{
        /* use hypre's GPU memory pool */
        //HYPRE_SetGPUMemoryPoolSize(bin_growth, min_bin, max_bin, max_bytes);
      //}
      //else if (useUmpireGpuMemPool)
      //{
        /* or use Umpire GPU memory pool */
      //  HYPRE_SetUmpireUMPoolName("HYPRE_UM_POOL_TEST");
      //  HYPRE_SetUmpireDevicePoolName("HYPRE_DEVICE_POOL_TEST");
      //}
      inited = true; 
    }
    nvtxRangePop(); 
  }

  PreconditionerHypre() :
      Preconditioner(),
      num_blocks_(0),
      block_indices_(Teuchos::null),
      HyprePrecond_(),
      returned_code_(0)
  {
    init(); 
  }

  virtual void set_inverse_parameters(Teuchos::ParameterList& list) override final;
  virtual void initializeInverse() override final;
  virtual void computeInverse() override final;
  virtual int applyInverse(const Vector_type& v, Vector_type& hv) const override final;

  virtual int returned_code() const override final { return returned_code_; }
  virtual std::string returned_code_string() const override final {
    return "success";
  }  
  
  // Need to be public for Kokkos::parallel_for
  void copy_matrix_(); 


 private:
  void Init_(){};
  void InitBoomer_();
  void InitEuclid_();
  
  Teuchos::RCP<const Map_type> make_contiguous_(Teuchos::RCP<const RowMatrix_type> &Matrix); 

  Teuchos::ParameterList plist_;
  Teuchos::RCP<VerboseObject> vo_;
  
  Teuchos::RCP<std::vector<int> > block_indices_;
  int num_blocks_;
  mutable int returned_code_;

  Teuchos::RCP<const Map_type> GloballyContiguousRowMap_;
  Teuchos::RCP<const Map_type> GloballyContiguousColMap_;

  HYPRE_Solver HyprePrecond_;
  HYPRE_ParCSRMatrix ParMatrix_;
  HYPRE_IJMatrix HypreA_;
  HYPRE_ParVector ParX_;
  HYPRE_ParVector ParY_;  
  HYPRE_IJVector XHypre_;
  HYPRE_IJVector YHypre_;
  Teuchos::RCP<hypre_ParVector> XVec_;
  Teuchos::RCP<hypre_ParVector> YVec_;

  Teuchos::RCP<RowMatrix_type> h_row;

  static bool inited; 
  HyprePreconditioners PrecondType = Boomer; 

};

}  // namespace AmanziSolvers
}  // namespace Amanzi



#endif
