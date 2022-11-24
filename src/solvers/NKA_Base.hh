/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Neil N. Carlson <neil.n.carlson@gmail.com>

  NONLINEAR_KRYLOV_ACCELERATOR

  This code implements the nonlinear Krylov accelerator introduced in [1]
  for inexact Newton's (IN) method, where the correction equation of
  Newton's method is only approximately solved because the Jacobian matrix
  is approximated and/or the linear system is not solved exactly.  Placed
  in the iteration loop, this black-box accelerator listens to the sequence
  of inexact corrections and replaces them with accelerated corrections;
  the resulting method is a type of accelerated inexact Newton (AIN) method.
  Note that an IN iteration is merely a standard fixed point iteration for
  a preconditioned system, and so this accelerator is more generally
  applicable to fixed point iterations.

  This code is a straightforward translation of the original Fortran 95
  implementation into C.

  [1] N.N.Carlson and K.Miller, "Design and application of a gradient-
      weighted moving finite element code I: in one dimension", SIAM J.
      Sci. Comput;, 19 (1998), pp. 728-765.  See section 9.

  Copyright (c) 2009  Neil N. Carlson

  Permission is hereby granted, free of charge, to any person obtaining a
  copy of this software and associated documentation files (the "Software"),
  to deal in the Software without restriction, including without limitation
  the rights to use, copy, modify, merge, publish, distribute, sublicense,
  and/or sell copies of the Software, and to permit persons to whom the
  Software is furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included
  in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
  DEALINGS IN THE SOFTWARE.

  NOTE:
   The following code was adapted from Neil Carlson's  NLKAIN code that is on
   Sourceforge, see http://sourceforge.net/projects/nlkain/
*/

#ifndef AMANZI_NKA_BASE_HH_
#define AMANZI_NKA_BASE_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseMatrix.h"

#include "VerboseObject.hh"
#include "dbc.hh"

namespace Amanzi {
namespace AmanziSolvers {

#define NKA_EOL -1

template <class Vector, class VectorSpace>
class NKA_Base {
 public:
  NKA_Base(int mvec, double vtol, const VectorSpace& map);

  void Init(Teuchos::ParameterList& plist)
  {
    vo_ = Teuchos::rcp(new VerboseObject("NKA_Base", plist));
  }


  void Relax();
  void Restart();
  void
  Correction(const Vector& f, Vector& v, const Teuchos::Ptr<const Vector>& oldv = Teuchos::null);

 private:
  bool subspace_; // boolean: a nonempty subspace
  bool pending_;  // contains pending vectors -- boolean
  int mvec_;      // maximum number of subspace vectors
  double vtol_;   // vector drop tolerance

  std::vector<Teuchos::RCP<Vector>> v_; // subspace storage
  std::vector<Teuchos::RCP<Vector>> w_; // function difference vectors
  Epetra_SerialDenseMatrix h_;          // matrix of w vector inner products

  // Linked-list organization of the vector storage.
  int first_v_;             // index of first_v subspace vector
  int last_v_;              // index of last_v subspace vector
  int free_v_;              // index of the initial vector in free storage linked list
  std::vector<int> next_v_; // next_v index link field
  std::vector<int> prev_v_; // previous index link field in doubly-linked subspace v

  Teuchos::RCP<VerboseObject> vo_;
};


/* ******************************************************************
 * Allocate memory
 ***************************************************************** */
template <class Vector, class VectorSpace>
NKA_Base<Vector, VectorSpace>::NKA_Base(int mvec, double vtol, const VectorSpace& map)
  : subspace_(false), pending_(false), mvec_(std::max(mvec, 1)), vtol_(vtol)
{
  v_.resize(mvec_ + 1);
  w_.resize(mvec_ + 1);

  for (int i = 0; i < mvec_ + 1; i++) {
    v_[i] = Teuchos::rcp(new Vector(map));
    v_[i]->PutScalar(0.0);
    w_[i] = Teuchos::rcp(new Vector(map));
    w_[i]->PutScalar(0.0);
  }

  h_.Shape(mvec_ + 1, mvec_ + 1);

  next_v_.resize(mvec_ + 1, NKA_EOL);
  prev_v_.resize(mvec_ + 1, NKA_EOL);
};


/* ******************************************************************
 * TBW
 ***************************************************************** */
template <class Vector, class VectorSpace>
void
NKA_Base<Vector, VectorSpace>::Relax()
{
  if (pending_) {
    // Drop the initial slot where the pending_ vectors are stored.
    AMANZI_ASSERT(first_v_ >= 0);

    int new_v = first_v_;
    first_v_ = next_v_[first_v_];
    if (first_v_ == NKA_EOL) {
      last_v_ = NKA_EOL;
    } else {
      prev_v_[first_v_] = NKA_EOL;
    }

    // Update the free storage list.
    next_v_[new_v] = free_v_;
    free_v_ = new_v;
    pending_ = false;
  }
}


/* ******************************************************************
 * TBW
 ***************************************************************** */
template <class Vector, class VectorSpace>
void
NKA_Base<Vector, VectorSpace>::Restart()
{
  // No vectors are stored.
  first_v_ = NKA_EOL;
  last_v_ = NKA_EOL;
  subspace_ = false;
  pending_ = false;

  // Initialize the free storage linked list.
  free_v_ = 0;
  for (int k = 0; k < mvec_; k++) { next_v_[k] = k + 1; }
  next_v_[mvec_] = NKA_EOL;
}


/* ******************************************************************
 * TBW
 ***************************************************************** */
template <class Vector, class VectorSpace>
void
NKA_Base<Vector, VectorSpace>::Correction(const Vector& f,
                                          Vector& dir,
                                          const Teuchos::Ptr<const Vector>& old_dir)
{
  // UPDATE THE ACCELERATION SUBSPACE
  if (pending_) {
    AMANZI_ASSERT(first_v_ != NKA_EOL);

    auto wp = w_[first_v_];
    auto vp = v_[first_v_];

    // next_v_ function difference w_1
    wp->Update(-1., f, 1.);

    // If the function difference is 0, we can't update the subspace with
    // this data; so we toss it out and continue.  In this case it is likely
    // that the outer iterative solution procedure has gone badly awry
    // (unless the function value is itself 0), and we merely want to do
    // something reasonable here and hope that situation is detected on the
    // outside.
    double s = 0.;
    double sv = 0.;
    wp->Norm2(&s);
    vp->Norm2(&sv);
    double sf = 0.;
    f.Norm2(&sf);

    if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "L2 norms (unscaled): du_prev, du_NKA, du = " << s << "," << sv << "," << sf
                 << std::endl;
    }

    bool too_close = (s == 0);

    if (!too_close) {
      // update the NKA space if the previously returned direction had been
      // modified before being applied
      if (old_dir != Teuchos::null) *vp = *old_dir;

      // Normalize w_1 and apply same factor to v_1.
      wp->Scale(1.0 / s);
      vp->Scale(1.0 / s);

      s = 0.;
      wp->Norm2(&s);
      vp->Norm2(&sv);
    }

    if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "Dot products = " << s << "," << sv << std::endl;
    }
    if (too_close) {
      // if our rescaled difference is too small,
      if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "Dot product confused = " << s << "," << sv << " tossing iterate"
                   << std::endl;
      }
      // nka_relax sets pending to false
      Relax();
    }
  }

  if (pending_) {
    AMANZI_ASSERT(first_v_ != NKA_EOL);
    auto wp = w_[first_v_];
    double s = 0.;
    wp->Norm2(&s);

    // Update H.
    for (int k = next_v_[first_v_]; k != NKA_EOL; k = next_v_[k]) {
      double hk = 0.;
      w_[k]->Norm2(&hk);
      hk = 0.;
      int ierr = wp->Dot(*w_[k], &hk);
      AMANZI_ASSERT(!ierr);
      h_[first_v_][k] = hk;
      hk = 0.;
      w_[k]->Dot(*wp, &hk);
    }

    // CHOLESKI FACTORIZATION OF H = W^t W
    // original matrix kept in the upper triangle (implicit unit diagonal)
    // lower triangle holds the factorization

    // Trivial initial factorization stage.
    int nvec = 1;
    h_[first_v_][first_v_] = 1.0;

    // The accelerated correction
    for (int k = next_v_[first_v_]; k != NKA_EOL; k = next_v_[k]) {
      // Maintain at most MVEC_ vectors.
      if (++nvec > mvec_) {
        // Drop the last_v_ vector and update the free storage list.
        AMANZI_ASSERT(last_v_ == k);
        next_v_[last_v_] = free_v_;
        free_v_ = k;
        last_v_ = prev_v_[k];
        next_v_[last_v_] = NKA_EOL;
        break;
      }

      // Single stage of Choleski factorization.
      double hkk = 1.0;
      for (int j = first_v_; j != k; j = next_v_[j]) {
        double hkj = h_[j][k];
        for (int i = first_v_; i != j; i = next_v_[i]) { hkj -= h_[k][i] * h_[j][i]; }
        hkj /= h_[j][j];
        h_[k][j] = hkj;
        hkk -= hkj * hkj;
      }

      if (hkk > std::pow(vtol_, 2)) {
        h_[k][k] = std::sqrt(hkk);
      } else {
        if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
          Teuchos::OSTab tab = vo_->getOSTab();
          *vo_->os() << "Vectors are linearly dependent, hkk=" << hkk << ", tossing iterate"
                     << std::endl;
        }

        // The current w nearly lies in the span of the previous vectors k vectors.
        // Drop the kth vector and keep the new one.
        AMANZI_ASSERT(prev_v_[k] !=
                      NKA_EOL); // when k is the only vector, its diagonal is 1.  Internal error!
        next_v_[prev_v_[k]] = next_v_[k];
        if (next_v_[k] == NKA_EOL) {
          last_v_ = prev_v_[k];
        } else {
          prev_v_[next_v_[k]] = prev_v_[k];
        }
        // update the free storage list,
        next_v_[k] = free_v_;
        free_v_ = k;

        // back-up and move on to the next_v_ vector.
        k = prev_v_[k];
        nvec--;
      }
    }

    AMANZI_ASSERT(first_v_ != NKA_EOL);
    subspace_ = true; // the acceleration subspace_ isn't empty
  }

  // ACCELERATED CORRECTION
  dir = f;
  Vector dir_update(dir);

  // Locate storage for the new vectors.
  AMANZI_ASSERT(free_v_ != NKA_EOL);
  int new_v = free_v_;
  free_v_ = next_v_[free_v_];

  // Save the original f for the next_v_ call.
  *w_[new_v] = f;

  if (subspace_) {
    std::vector<double> c(mvec_ + 1, 0.0);

    // Project f onto the span of the w vectors:
    // forward substitution
    for (int j = first_v_; j != NKA_EOL; j = next_v_[j]) {
      double cj = 0.;
      int ierr = dir.Dot(*w_[j], &cj);
      AMANZI_ASSERT(!ierr);

      for (int i = first_v_; i != j; i = next_v_[i]) { cj -= h_[j][i] * c[i]; }
      c[j] = cj / h_[j][j];
    }
    // backward substitution
    for (int j = last_v_; j != NKA_EOL; j = prev_v_[j]) {
      double cj = c[j];
      for (int i = last_v_; i != j; i = prev_v_[i]) { cj -= h_[i][j] * c[i]; }
      c[j] = cj / h_[j][j];
    }

    // The accelerated correction
    for (int k = first_v_; k != NKA_EOL; k = next_v_[k]) {
      dir_update.Update(1., *v_[k], -1., *w_[k], 0.);
      dir.Update(c[k], dir_update, 1.0);
    }
  }

  // Save the accelerated correction for the next_v_ call.
  *v_[new_v] = dir;


  // Prepend the new vectors to the list.
  prev_v_[new_v] = NKA_EOL;
  next_v_[new_v] = first_v_;
  if (first_v_ == NKA_EOL) {
    last_v_ = new_v;
  } else {
    prev_v_[first_v_] = new_v;
  }
  first_v_ = new_v;

  // The original f and accelerated correction are cached for the next_v_ call.
  pending_ = true;
};

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
