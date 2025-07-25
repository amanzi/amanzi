/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry

  Simple class template
*/

#ifndef AMANZI_CHEMISTRY_LU_SOLVER_HH_
#define AMANZI_CHEMISTRY_LU_SOLVER_HH_

#include <ostream>
#include <vector>

#include "errors.hh"

#include "MatrixBlock.hh"

namespace Amanzi {
namespace AmanziChemistry {

#define PREFIX
#define F77_LAPACK_MANGLE(lcase, UCASE) lcase##_
#define DGESV_F77 F77_LAPACK_MANGLE(dgesv, DGESV)

#ifdef __cplusplus
extern "C"
{
#endif

  void PREFIX
  DGESV_F77(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);

#ifdef __cplusplus
}
#endif


class LUSolver {
 public:
  LUSolver()
    : size_(0)
  {
    ipiv_.clear();
  }
  ~LUSolver() {};

  void Initialize(int size)
  {
    size_ = size;
    ipiv_.resize(size);
  }

  void Solve(MatrixBlock* A, std::vector<double>* b)
  {
    // scale A x = b
    int ierr, nrhs(1), n(A->size());
    for (int i = 0; i < n; ++i) {
      double big(0.0);
      for (int j = 0; j < n; ++j) {
        big = std::max(big, fabs((*A)(i, j)));
      }
      for (int j = 0; j < n; ++j) (*A)(i, j) /= big;
      (*b)[i] /= big;
    }

    // run solver
    DGESV_F77(&n, &nrhs, A->GetValues(), &n, &(ipiv_[0]), &(*b)[0], &n, &ierr);

    if (ierr != 0) {
      std::string msg("LUSolver::Decomposition() : Singular matrix.");
      Exceptions::amanzi_throw(Errors::Message(msg));
    }
  }

 private:
  int size_;
  std::vector<int> ipiv_;
};

} // namespace AmanziChemistry
} // namespace Amanzi

#endif
