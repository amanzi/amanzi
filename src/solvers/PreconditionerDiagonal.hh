/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/
//! Diagonal preconditioner.

/*!

Simply applys the pointwise inverse of the diagonal of the matrix as an
extremely cheap matrix.

This is provided when using the `"preconditioning method`"=`"diagonal`" in the
`Preconditioner`_ spec.

No parameters are required.

*/


#ifndef AMANZI_PRECONDITIONER_DIAGONAL_HH_
#define AMANZI_PRECONDITIONER_DIAGONAL_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
//#include "Epetra_MultiVector.h"
//#include "Epetra_RowMatrix.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

#include "nvToolsExt.h"


namespace Amanzi {
namespace AmanziSolvers {

class PreconditionerDiagonal : public Preconditioner {
 public:
  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) override final {};
  virtual void initializeInverse() override final {}
  virtual void computeInverse() override final {
    nvtxRangePushA("PreconditionerDiagonal::computeInverse");    
    #ifdef OUTPUT_CUDA 
    
    Teuchos::RCP<Teuchos::FancyOStream> out = getFancyOStream (Teuchos::rcpFromRef (std::cout));
    h_->describe(*out,Teuchos::VERB_EXTREME);
    #endif 
    diagonal_ = Teuchos::rcp(new Vector_type(h_->getRowMap()));
    diagonal_->putScalar(0.);
    h_->getLocalDiagCopy(*diagonal_);

    #ifdef OUTPUT_CUDA
    {
      std::cout<<"getLocalDiagCopy: "; 
    auto diag = diagonal_->getLocalViewDevice ();
    Kokkos::parallel_for(
      "",
      diag.extent(0), 
      KOKKOS_LAMBDA(const int i){
        printf("%.4f - ",diag(i,0)); 
      }
    );  
    Kokkos::fence(); 
    std::cout<<std::endl;
    }
    #endif 
    diagonal_->reciprocal(*diagonal_);

    nvtxRangePop();
  };

  virtual int applyInverse(const Vector_type& v, Vector_type& hv) const override final {
    nvtxRangePushA("PreconditionerDiagonal::applyInverse");

    AMANZI_ASSERT(diagonal_.get()); // Compute called
    auto diag = diagonal_->getLocalViewDevice ();
    #ifdef OUTPUT_CUDA
    std::cout<<"Diagonal : applyInverse"<<std::endl;
    Kokkos::parallel_for(
      "",
      diag.extent(0), 
      KOKKOS_LAMBDA(const int i){
        printf("%.4f - ",diag(i,0)); 
      }
    );  
    Kokkos::fence(); 
    std::cout<<std::endl;
    #endif 
    hv.elementWiseMultiply(1., v, *diagonal_, 0.);
    nvtxRangePop();

    return 0;
  }

  virtual int returned_code() const override final  { return returned_code_; }
  virtual std::string returned_code_string() const override final {
    if (returned_code_ == 0) return "success";
    return "failed ReciprocalMultiply()";
  }

 private:
  Teuchos::RCP<Vector_type> diagonal_;
  mutable int returned_code_;
};

}  // namespace AmanziSolvers
}  // namespace Amanzi



#endif
