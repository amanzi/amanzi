/*
  License: BSD
  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov) (ATS version)
  MatrixMFD provides a mimetic discretization for the elliptic operator div K grad u.

*/

#ifndef BLOCK_MATRIX_HH_
#define BLOCK_MATRIX_HH_

#include <strings.h>

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"


#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include "CompositeVectorSpace.hh"
#include "CompositeVector.hh"
#include "CompositeMatrix.hh"
#include "EpetraMatrix.hh"

namespace Amanzi {
namespace Operators {

class BlockMatrix : public CompositeMatrix{
 public:
  BlockMatrix(Teuchos::RCP<CompositeVectorSpace> space){ space_ = space;};

  BlockMatrix(const BlockMatrix& other){
    space_ = other.space_;
    blocks.resize( other.blocks.size() );
    for (int i=0; i<blocks.size(); i++){
       blocks[i]    = other.blocks[i];
    }
    for (int i=0; i<blocks_pc.size(); i++){
      blocks_pc[i] = other.blocks_pc[i];
    }
  };

  ~BlockMatrix(){};

 // CompositeMatrix methods
  virtual const CompositeVectorSpace& DomainMap() const {
    return *space_; }

  virtual const CompositeVectorSpace& RangeMap() const {
    return *space_; }

  virtual Teuchos::RCP<CompositeMatrix> Clone() const {
    return Teuchos::rcp(new BlockMatrix(*this)); }

  virtual int Apply(const CompositeVector& X,
                     CompositeVector& Y) const
  {
    int i=0;
    int err = 0;
    for (CompositeVector::name_iterator comp=X.begin(); comp!=X.end(); ++comp){
      const Epetra_MultiVector& temp_x = *(X.ViewComponent(*comp,false));
      Epetra_MultiVector& temp_y = *(Y.ViewComponent(*comp,false));

      err = blocks[i]->Apply(temp_x, temp_y);      
      i++;
    }
    return err;

  };

  virtual int ApplyInverse(const CompositeVector& X,
                            CompositeVector& Y) const
  {
    int i=0;
    int err = 0;


    i=0;
    for (CompositeVector::name_iterator comp=X.begin(); comp!=X.end(); ++comp, ++i){
    //for (CompositeVector::name_iterator comp=X.end(); comp!=X.begin(); --comp){
      //--comp;

      const Epetra_MultiVector& temp_x = *(X.ViewComponent(*comp,false));
      Epetra_MultiVector& temp_y = *(Y.ViewComponent(*comp,false));
      //std::cout<<i<<" "<<*comp<<"\n";
      err = blocks_pc[i]->ApplyInverse(temp_x, temp_y);      
    }

    return err;

  };


  void SetNumberofBlocks(int n){  
    blocks.resize(n);
    blocks_pc.resize(n);
  }
  int GetNumberofBlocks(){ return blocks.size();}

  void SetBlock(int i, Teuchos::RCP<Epetra_FECrsMatrix> A){
    ASSERT(i < blocks.size());
    blocks[i] = A;
  }

  void SetPrec(int i,  Teuchos::RCP<AmanziPreconditioners::Preconditioner> A_pc){
    ASSERT(i < blocks_pc.size());
    blocks_pc[i] = A_pc;
  }

  void InitPreconditioner(const Teuchos::ParameterList& plist_){
    if (plist_.isSublist("preconditioner")) {
      Teuchos::ParameterList pc_list = plist_.sublist("preconditioner");
      AmanziPreconditioners::PreconditionerFactory pc_fac;
      for (int i=0; i < blocks_pc.size(); i++){
	blocks_pc[i] = pc_fac.Create(pc_list);
      }
    }
  }

  void UpdatePreconditioner() {
    for (int i=0; i < blocks_pc.size(); i++){
      if (blocks_pc[i] == Teuchos::null) {
	Errors::Message msg("BlockMFD::UpdatePreconditioner block is not initialized");
	Exceptions::amanzi_throw(msg);
      }
      blocks_pc[i]->Destroy();
      blocks_pc[i]->Update(blocks[i]);
    }
  }



   // VectorSpace describing both domain and range
  Teuchos::RCP<CompositeVectorSpace> space_;

  Teuchos::ParameterList prec_list;

  protected:
  mutable std::vector<Teuchos::RCP<Epetra_FECrsMatrix> > blocks;
  mutable std::vector<Teuchos::RCP<AmanziPreconditioners::Preconditioner> > blocks_pc;

};

}
}

#endif
