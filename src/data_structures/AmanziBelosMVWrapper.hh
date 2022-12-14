/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Alicia Klinvex (amklinv@sandia.gov)
*/

#ifndef AMANZI_BELOS_MV_WRAPPER_HH_
#define AMANZI_BELOS_MV_WRAPPER_HH_

#include <BelosMultiVec.hpp>

namespace Amanzi {
// This is not an optimal implementation, but we can improve it later
template <class Vector>
class CompositeMultiVector : public Belos::MultiVec<double> {
 private:
  typedef double ScalarType;
  typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

 public:
  // Constructors
  CompositeMultiVector(Teuchos::RCP<Vector> cmv)
  {
    cmv_.resize(1);
    cmv_[0] = cmv;
  }
  CompositeMultiVector(std::vector<Teuchos::RCP<Vector>>& cmv) { cmv_ = cmv; }

  // Creation methods for new multivectors
  CompositeMultiVector* Clone(const int numvecs) const;
  CompositeMultiVector* CloneCopy() const;
  CompositeMultiVector* CloneCopy(const std::vector<int>& index) const;
  CompositeMultiVector* CloneViewNonConst(const std::vector<int>& index);
  const CompositeMultiVector* CloneView(const std::vector<int>& index) const;

  // Dimension information methods
  ptrdiff_t GetGlobalLength() const { return cmv_[0]->GlobalLength(); } // TODO: Fix this
  int GetNumberVecs() const { return cmv_.size(); }

  // Update methods
  void MvTimesMatAddMv(double alpha,
                       const Belos::MultiVec<double>& A,
                       const Teuchos::SerialDenseMatrix<int, double>& B,
                       double beta);
  void MvAddMv(ScalarType alpha,
               const Belos::MultiVec<double>& A,
               ScalarType beta,
               const Belos::MultiVec<double>& B);
  void MvScale(ScalarType alpha);
  void MvScale(const std::vector<ScalarType>& alpha);
  void MvTransMv(ScalarType alpha,
                 const Belos::MultiVec<double>& A,
                 Teuchos::SerialDenseMatrix<int, ScalarType>& B) const;
  void MvDot(const Belos::MultiVec<double>& A, std::vector<ScalarType>& b) const;

  // Norm method
  void MvNorm(std::vector<MagnitudeType>& normvec, Belos::NormType type = Belos::TwoNorm) const;

  // Initialization methods
  void SetBlock(const Belos::MultiVec<double>& A, const std::vector<int>& index);
  void MvRandom();
  void MvInit(ScalarType alpha);

  // Print method
  void MvPrint(std::ostream& os) const;

  Teuchos::RCP<Vector> getVector(int i) const { return cmv_[i]; }

 private:
  std::vector<Teuchos::RCP<Vector>> cmv_;
}; // end class CompositeMultiVector

template <class Vector>
CompositeMultiVector<Vector>*
CompositeMultiVector<Vector>::Clone(const int numvecs) const
{
  std::vector<Teuchos::RCP<Vector>> newcmv(numvecs);
  for (int i = 0; i < numvecs; i++) newcmv[i] = Teuchos::rcp(new Vector(*cmv_[0]));
  CompositeMultiVector* temp = new CompositeMultiVector(newcmv);
  temp->MvInit(0);

  return temp;
}

template <class Vector>
CompositeMultiVector<Vector>*
CompositeMultiVector<Vector>::CloneCopy() const
{
  std::vector<Teuchos::RCP<Vector>> newcmv(cmv_.size());
  for (int i = 0; i < cmv_.size(); i++) newcmv[i] = Teuchos::rcp(new Vector(*cmv_[i]));
  CompositeMultiVector* temp = new CompositeMultiVector(newcmv);

  return temp;
}

template <class Vector>
CompositeMultiVector<Vector>*
CompositeMultiVector<Vector>::CloneCopy(const std::vector<int>& index) const
{
  std::vector<Teuchos::RCP<Vector>> newcmv(index.size());
  for (int i = 0; i < index.size(); i++) newcmv[i] = Teuchos::rcp(new Vector(*cmv_[index[i]]));
  CompositeMultiVector* temp = new CompositeMultiVector(newcmv);

  return temp;
}

template <class Vector>
CompositeMultiVector<Vector>*
CompositeMultiVector<Vector>::CloneViewNonConst(const std::vector<int>& index)
{
  std::vector<Teuchos::RCP<Vector>> newcmv(index.size());
  for (int i = 0; i < index.size(); i++) newcmv[i] = cmv_[index[i]];
  CompositeMultiVector* temp = new CompositeMultiVector(newcmv);

  return temp;
}

template <class Vector>
const CompositeMultiVector<Vector>*
CompositeMultiVector<Vector>::CloneView(const std::vector<int>& index) const
{
  std::vector<Teuchos::RCP<Vector>> newcmv(index.size());
  for (int i = 0; i < index.size(); i++) newcmv[i] = cmv_[index[i]];
  CompositeMultiVector* temp = new CompositeMultiVector(newcmv);

  return temp;
}

template <class Vector>
void
CompositeMultiVector<Vector>::MvTimesMatAddMv(double alpha,
                                              const Belos::MultiVec<double>& A,
                                              const Teuchos::SerialDenseMatrix<int, double>& B,
                                              double beta)
{
  const CompositeMultiVector* cmvA = dynamic_cast<const CompositeMultiVector*>(&A);
  TEUCHOS_TEST_FOR_EXCEPTION(
    cmvA == 0,
    std::invalid_argument,
    "CompositeMultiVector::MvTimesMatAddMv: A must be a CompositeMultiVector");

  MvScale(beta);
  for (int r = 0; r < B.numRows(); r++) {
    for (int c = 0; c < B.numCols(); c++) {
      // mv[c] += alpha b(r,c) A[r]
      cmv_[c]->Update(alpha * B(r, c), *cmvA->cmv_[r], 1);
    }
  }
}

template <class Vector>
void
CompositeMultiVector<Vector>::MvAddMv(double alpha,
                                      const Belos::MultiVec<double>& A,
                                      double beta,
                                      const Belos::MultiVec<double>& B)
{
  const CompositeMultiVector* cmvA = dynamic_cast<const CompositeMultiVector*>(&A);
  TEUCHOS_TEST_FOR_EXCEPTION(cmvA == 0,
                             std::invalid_argument,
                             "CompositeMultiVector::MvAddMv: A must be a CompositeMultiVector");
  const CompositeMultiVector* cmvB = dynamic_cast<const CompositeMultiVector*>(&B);
  TEUCHOS_TEST_FOR_EXCEPTION(cmvB == 0,
                             std::invalid_argument,
                             "CompositeMultiVector::MvAddMv: B must be a CompositeMultiVector");

  for (int i = 0; i < cmv_.size(); i++)
    cmv_[i]->Update(alpha, *cmvA->cmv_[i], beta, *cmvB->cmv_[i], 0);
}

template <class Vector>
void
CompositeMultiVector<Vector>::MvScale(ScalarType alpha)
{
  for (int i = 0; i < cmv_.size(); i++) cmv_[i]->Scale(alpha);
}

template <class Vector>
void
CompositeMultiVector<Vector>::MvScale(const std::vector<ScalarType>& alpha)
{
  for (int i = 0; i < cmv_.size(); i++) cmv_[i]->Scale(alpha[i]);
}

template <class Vector>
void
CompositeMultiVector<Vector>::MvTransMv(ScalarType alpha,
                                        const Belos::MultiVec<double>& A,
                                        Teuchos::SerialDenseMatrix<int, ScalarType>& B) const
{
  const CompositeMultiVector* cmvA = dynamic_cast<const CompositeMultiVector*>(&A);
  TEUCHOS_TEST_FOR_EXCEPTION(cmvA == 0,
                             std::invalid_argument,
                             "CompositeMultiVector::MvTransMv: A must be a CompositeMultiVector");

  for (int r = 0; r < cmvA->GetNumberVecs(); r++) {
    for (int c = 0; c < GetNumberVecs(); c++) cmv_[c]->Dot(*cmvA->cmv_[r], &B(r, c));
  }
}

template <class Vector>
void
CompositeMultiVector<Vector>::MvDot(const Belos::MultiVec<double>& A,
                                    std::vector<ScalarType>& b) const
{
  const CompositeMultiVector* cmvA = dynamic_cast<const CompositeMultiVector*>(&A);
  TEUCHOS_TEST_FOR_EXCEPTION(cmvA == 0,
                             std::invalid_argument,
                             "CompositeMultiVector::MvAddMv: A must be a CompositeMultiVector");

  for (int i = 0; i < cmv_.size(); i++) { cmv_[i]->Dot(*cmvA->cmv_[i], &b[i]); }
}

template <class Vector>
void
CompositeMultiVector<Vector>::MvNorm(std::vector<MagnitudeType>& normvec,
                                     Belos::NormType type) const
{
  if (type == Belos::TwoNorm) {
    for (int i = 0; i < cmv_.size(); i++) cmv_[i]->Norm2(&normvec[i]);
  } else if (type == Belos::OneNorm) {
    for (int i = 0; i < cmv_.size(); i++) cmv_[i]->Norm1(&normvec[i]);
  } else // if (type == Belos::InfNorm)
  {
    for (int i = 0; i < cmv_.size(); i++) cmv_[i]->NormInf(&normvec[i]);
  }
}

template <class Vector>
void
CompositeMultiVector<Vector>::SetBlock(const Belos::MultiVec<double>& A,
                                       const std::vector<int>& index)
{
  const CompositeMultiVector* cmvA = dynamic_cast<const CompositeMultiVector*>(&A);
  TEUCHOS_TEST_FOR_EXCEPTION(cmvA == 0,
                             std::invalid_argument,
                             "CompositeMultiVector::SetBlock: A must be a CompositeMultiVector");

  for (int i = 0; i < index.size(); i++) *cmv_[index[i]] = *cmvA->cmv_[i];
}

template <class Vector>
void
CompositeMultiVector<Vector>::MvRandom()
{
  for (int i = 0; i < cmv_.size(); i++) cmv_[i]->Random();
}

template <class Vector>
void
CompositeMultiVector<Vector>::MvInit(ScalarType alpha)
{
  for (int i = 0; i < cmv_.size(); i++) cmv_[i]->PutScalar(alpha);
}

template <class Vector>
void
CompositeMultiVector<Vector>::MvPrint(std::ostream& os) const
{
  for (int i = 0; i < cmv_.size(); i++) cmv_[i]->Print(os);
}

} // namespace Amanzi

#endif
