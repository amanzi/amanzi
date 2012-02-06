/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Implementation of TreeVector, a nested, hierarchical data structure for PK
   hierarchies.  This nested vector allows each physical PK to push back
   Epetra_MultiVectors to store their solution, and allows MPCs to push back
   TreeVectors in a nested format.  It also provides an implementation of the
   Vector interface for use with time integrators/nonlinear solvers.
   ------------------------------------------------------------------------- */

#include "dbc.hh"
#include "TreeVector.hh"

namespace Amanzi {

  // Basic constructor of an empty TreeVector
  TreeVector::TreeVector(std::string name) : name_(name) {
  };

  // TreeVector::TreeVector(std::string name, const Teuchos::RCP<Epetra_MultiVector> &data) :
  //     name_(name) {
  //   Teuchos::RCP<Epetra_MultiVector> new_EMV
  //       = Teuchos::rcp(new Epetra_MultiVector( data->Map(), data->NumVectors() ));
  //   data_.push_back(new_EMV);
  // };

  // TreeVector::TreeVector(std::string name, const Teuchos::RCP<Epetra_Vector> &data) :
  //     name_(name) {
  //   Teuchos::RCP<Epetra_MultiVector> new_EMV
  //       = Teuchos::rcp(new Epetra_MultiVector( data->Map(), data->NumVectors() ));
  //   data_.push_back(new_EMV);
  // };

  // Copy constructors
  TreeVector::TreeVector(std::string name, const Teuchos::RCP<TreeVector> &tv) : name_(name) {
    for (std::vector< Teuchos::RCP<Epetra_MultiVector> >::const_iterator datum = tv->data_.begin();
         datum != tv->data_.end(); ++datum) {

      Teuchos::RCP<Epetra_MultiVector> new_EMV
        = Teuchos::rcp(new Epetra_MultiVector( (*datum)->Map(), (*datum)->NumVectors()  ) );
      data_.push_back(new_EMV);
    }
    for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator sv = tv->subvecs_.begin();
         sv != tv->subvecs_.end(); ++sv) {

      Teuchos::RCP<Amanzi::TreeVector> new_TV
        = Teuchos::rcp(new Amanzi::TreeVector((*sv)->name_, *sv) );
      subvecs_.push_back(new_TV);
    }
  };

  TreeVector::TreeVector(std::string name, const TreeVector &tv) : name_(name) {
    for (std::vector< Teuchos::RCP<Epetra_MultiVector> >::const_iterator datum = tv.data_.begin();
         datum != tv.data_.end(); ++datum) {

      Teuchos::RCP<Epetra_MultiVector> new_EMV
        = Teuchos::rcp(new Epetra_MultiVector( (*datum)->Map(), (*datum)->NumVectors()  ) );
      data_.push_back(new_EMV);
    }
    for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator sv = tv.subvecs_.begin();
         sv != tv.subvecs_.end(); ++sv) {

      Teuchos::RCP<Amanzi::TreeVector> new_TV
        = Teuchos::rcp(new Amanzi::TreeVector((*sv)->name_, *sv) );
      subvecs_.push_back(new_TV);
    }
  };

  TreeVector::TreeVector(const TreeVector &tv) : name_(std::string("")) {
    for (std::vector< Teuchos::RCP<Epetra_MultiVector> >::const_iterator datum = tv.data_.begin();
         datum != tv.data_.end(); ++datum) {

      Teuchos::RCP<Epetra_MultiVector> new_EMV
        = Teuchos::rcp(new Epetra_MultiVector( (*datum)->Map(), (*datum)->NumVectors()  ) );
      data_.push_back(new_EMV);
    }
    for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator sv = tv.subvecs_.begin();
         sv != tv.subvecs_.end(); ++sv) {

      Teuchos::RCP<Amanzi::TreeVector> new_TV
        = Teuchos::rcp(new Amanzi::TreeVector((*sv)->name_, *sv) );
      subvecs_.push_back(new_TV);
    }
  };

  TreeVector::TreeVector(const Teuchos::RCP<TreeVector>& tv) : name_(std::string("")) {
    for (std::vector< Teuchos::RCP<Epetra_MultiVector> >::const_iterator datum = tv->data_.begin();
         datum != tv->data_.end(); ++datum) {

      Teuchos::RCP<Epetra_MultiVector> new_EMV
        = Teuchos::rcp(new Epetra_MultiVector( (*datum)->Map(), (*datum)->NumVectors()  ) );
      data_.push_back(new_EMV);
    }
    for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator sv = tv->subvecs_.begin();
         sv != tv->subvecs_.end(); ++sv) {

      Teuchos::RCP<Amanzi::TreeVector> new_TV
        = Teuchos::rcp(new Amanzi::TreeVector((*sv)->name_, *sv) );
      subvecs_.push_back(new_TV);
    }
  };

  TreeVector& TreeVector::operator=(double value) {
    // Set data of the vector.
    for (std::vector< Teuchos::RCP<Epetra_MultiVector> >::iterator datum = data_.begin();
         datum != data_.end(); ++datum) {
      (*datum)->PutScalar(value);
    }
    for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
         subvec != subvecs_.end(); ++subvec) {
      *(*subvec) = value;
    }
    return *this;
  };

  TreeVector& TreeVector::operator=(const Epetra_Vector &value) {
    // Set data at this node.
    ASSERT(subvecs_.size() == 0);
    ASSERT(data_.size() == 1);

    Epetra_MultiVector datum = *data_[0];
    ASSERT(datum.NumVectors() == 1);
    *(datum(0)) = value;
    return *this;
  };

  TreeVector& TreeVector::operator=(const Epetra_MultiVector &value) {
    ASSERT(subvecs_.size() == 0);
    ASSERT(data_.size() == 1);

    *data_[0] = value;
    return *this;
  };

  TreeVector& TreeVector::operator=(const TreeVector &value) {
    // Copy values at this node and all child nodes.
    ASSERT(subvecs_.size() == value.subvecs_.size());
    ASSERT(data_.size() == value.data_.size());

    std::vector< Teuchos::RCP<Epetra_MultiVector> >::const_iterator value_datum = value.data_.begin();
    for (std::vector< Teuchos::RCP<Epetra_MultiVector> >::iterator datum = data_.begin();
         datum != data_.end(); ++datum, ++value_datum) {
      *(*datum) = *(*value_datum);
    }

    std::vector< Teuchos::RCP<TreeVector> >::const_iterator value_subvec = value.subvecs_.begin();
    for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
         subvec != subvecs_.end(); ++subvec, value_subvec) {
      *(*subvec) = *(*value_subvec);
    }
    return *this;
  }

  int TreeVector::PutScalar(double scalar) {
    // Set all data of this node and all child nodes to scalar.
    for (std::vector< Teuchos::RCP<Epetra_MultiVector> >::iterator datum = data_.begin();
         datum != data_.end(); ++datum) {
      (*datum)->PutScalar(scalar);
    }
    for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
         subvec != subvecs_.end(); ++subvec) {
      (*subvec)->PutScalar(scalar);
    }
    return 0;
  }

  int TreeVector::NormInf(double* ninf) const {
    // Take the L_Inf norm of this.
    if (ninf == NULL) return 1;

    bool got_one_already(false);
    double ninf_loc;
    for (std::vector< Teuchos::RCP<Epetra_MultiVector> >::const_iterator datum = data_.begin();
         datum != data_.end(); ++datum)
      for (int i=0; i<(*datum)->NumVectors(); i++) {
        (*(*datum))(i)->NormInf(&ninf_loc);
        if (got_one_already) {
          if (ninf_loc > *ninf) {
            *ninf = ninf_loc;
          }
        } else {
          *ninf = ninf_loc;
          got_one_already = true;
        }
      }

    for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator subvec = subvecs_.begin();
         subvec != subvecs_.end(); ++subvec) {
      (*subvec)->NormInf(&ninf_loc);
      if (got_one_already) {
        if (ninf_loc > *ninf) {
          *ninf = ninf_loc;
        }
      } else {
        *ninf = ninf_loc;
        got_one_already = true;
      }
    }

    if (! got_one_already) {
      return 1;
    } else {
      return 0;
    }
  }

  void TreeVector::Print(ostream &os) const {
    // Print data to ostream for this node and all children.
    os << name_ << std::endl;

    for (std::vector< Teuchos::RCP<Epetra_MultiVector> >::const_iterator datum = data_.begin();
         datum != data_.end(); ++datum) {
      (*datum)->Print(os);
    }
    for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator subvec = subvecs_.begin();
         subvec != subvecs_.end(); ++subvec) {
      (*subvec)->Print(os);
    }
  }

  void TreeVector::Scale(double value) {
    // this <- value*this
    for (std::vector< Teuchos::RCP<Epetra_MultiVector> >::iterator datum = data_.begin();
         datum != data_.end(); ++datum) {
      (*datum)->Scale(value);
    }
    for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
         subvec != subvecs_.end(); ++subvec) {
      (*subvec)->Scale(value);
    }
  };

  void TreeVector::Shift(double scalarA) {
    // this <- this + scalarA
    for (unsigned int i = 0; i != data_.size(); ++i) {
      Epetra_MultiVector work(*data_[i]);
      work.PutScalar(1.0);
      data_[i]->Update(scalarA, work, 1.0);
    }
    for (unsigned int i = 0; i != subvecs_.size(); ++i) {
      subvecs_[i]->Shift(scalarA);
    }
  };

  int TreeVector::Dot(const Vector& other, double* result) const {
    // compute the dot product of all components of the tree vector
    // viewed as one flat vector

    const TreeVector *other_ptr = static_cast<const TreeVector *> (&other);
    if (!other_ptr) return 1;
    if (other_ptr->subvecs_.size() != subvecs_.size()) return 2;
    if (other_ptr->data_.size() != data_.size()) return 3;

    *result = 0.0;
    for (unsigned int i = 0; i != data_.size(); ++i) {
      double intermediate_result[(other_ptr->data_[i])->NumVectors()];
      int ierr = data_[i]->Dot( *(other_ptr->data_[i]),intermediate_result);
      if (ierr) return ierr;
      for (unsigned int j=0; j<(other_ptr->data_[i])->NumVectors(); ++j) {
        *result += intermediate_result[j];
      }
    }
    for (unsigned int i = 0; i != subvecs_.size(); ++i) {
      double intermediate_result;
      int ierr = subvecs_[i]->Dot(*(other_ptr->subvecs_[i]), &intermediate_result);
      if (ierr) return ierr;
      *result += intermediate_result;
    }

    return 0;
  };

  // this <- scalarA*A + scalarThis*this
  Vector& TreeVector::Update(double scalarA, const Vector& A, double scalarThis) {
    const TreeVector *A_ptr = static_cast<const TreeVector *> (&A);
    ASSERT(A_ptr);
    ASSERT(A_ptr->subvecs_.size() == subvecs_.size());
    ASSERT(A_ptr->data_.size() == data_.size());

    for (unsigned int i = 0; i != data_.size(); ++i) {
      data_[i]->Update(scalarA, *(A_ptr->data_[i]), scalarThis);
    }
    for (unsigned int i = 0; i != subvecs_.size(); ++i) {
      subvecs_[i]->Update(scalarA, *A_ptr->subvecs_[i], scalarThis);
    }

    return *this;
  };

  Vector& TreeVector::Update(double scalarA, const Vector& A,
                             double scalarB, const Vector& B, double scalarThis) {
    const TreeVector *A_ptr = static_cast<const TreeVector *> (&A);
    ASSERT(A_ptr);
    ASSERT(A_ptr->subvecs_.size() == subvecs_.size());
    ASSERT(A_ptr->data_.size() == data_.size());

    const TreeVector *B_ptr = dynamic_cast<const TreeVector *> (&B);
    ASSERT(B_ptr);
    ASSERT(B_ptr->subvecs_.size() == subvecs_.size());
    ASSERT(B_ptr->data_.size() == data_.size());

    for (unsigned int i = 0; i != data_.size(); ++i) {
      data_[i]->Update(scalarA, *(A_ptr->data_[i]), scalarB, *(B_ptr->data_[i]), scalarThis);
    }
    for (unsigned int i = 0; i != subvecs_.size(); ++i) {
      subvecs_[i]->Update(scalarA, *(A_ptr->subvecs_[i]), scalarB, *(B_ptr->subvecs_[i]),
                          scalarThis);
    }
    return *this;
  };

  int TreeVector::Multiply(double scalarAB, const Vector& A, const Vector& B, double scalarThis) {
    // this <- scalarAB * A@B + scalarThis*this  (@ is the elementwise product
    const TreeVector *A_ptr = static_cast<const TreeVector *> (&A);
    const TreeVector *B_ptr = static_cast<const TreeVector *> (&B);

    ASSERT(A_ptr);
    ASSERT(A_ptr->subvecs_.size() == subvecs_.size());
    ASSERT(A_ptr->data_.size() == data_.size());

    ASSERT(B_ptr);
    ASSERT(B_ptr->subvecs_.size() == subvecs_.size());
    ASSERT(B_ptr->data_.size() == data_.size());

    for (unsigned int i = 0; i != data_.size(); ++i) {
      data_[i]->Multiply(scalarAB, *(A_ptr->data_[i]), *(B_ptr->data_[i]), scalarThis);
    }
    for (unsigned int i = 0; i != subvecs_.size(); ++i) {
      subvecs_[i]->Multiply(scalarAB, *A_ptr->subvecs_[i], *B_ptr->subvecs_[i], scalarThis);
    }

    return 0;
  };


  int TreeVector::SubVector(std::string subname, Teuchos::RCP<TreeVector> &vec) {
    // Get a pointer to the sub-vector "subname".
    for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
         subvec != subvecs_.end(); ++subvec) {
      if (subname == (*subvec)->Name()) {
        vec = *subvec;
        return 0;
      }
    }
    for (std::vector< Teuchos::RCP<TreeVector> >::iterator subvec = subvecs_.begin();
         subvec != subvecs_.end(); ++subvec) {
      int ierr = (*subvec)->SubVector(subname, vec);
      if (!ierr) return ierr;
    }
    return 1;
  };

  int TreeVector::SubVector(std::string subname,
                            Teuchos::RCP<const TreeVector> &vec) const {
    // Get a pointer to the sub-vector "subname".
    for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator subvec = subvecs_.begin();
         subvec != subvecs_.end(); ++subvec) {
      if (subname == (*subvec)->Name()) {
        vec = *subvec;
        return 0;
      }
    }
    for (std::vector< Teuchos::RCP<TreeVector> >::const_iterator subvec = subvecs_.begin();
         subvec != subvecs_.end(); ++subvec) {
      int ierr = (*subvec)->SubVector(subname, vec);
      if (!ierr) return ierr;
    }
    return 1;
  };

  void TreeVector::PushBack(Teuchos::RCP<TreeVector>& subvec) {
    subvecs_.push_back(subvec);
  };

  void TreeVector::PushBack(Teuchos::RCP<Epetra_MultiVector>& data) {
    data_.push_back(data);
  };

} // namespace
