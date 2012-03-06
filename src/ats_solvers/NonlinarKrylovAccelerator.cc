//   NONLINEAR_KRYLOV_ACCELERATOR

//   Neil N. Carlson <neil.n.carlson@gmail.com>

//   This code implements the nonlinear Krylov accelerator introduced in [1]
//   for inexact Newton's (IN) method, where the correction equation of
//   Newton's method is only approximately solved because the Jacobian matrix
//   is approximated and/or the linear system is not solved exactly.  Placed
//   in the iteration loop, this black-box accelerator listens to the sequence
//   of inexact corrections and replaces them with accelerated corrections;
//   the resulting method is a type of accelerated inexact Newton (AIN) method.
//   Note that an IN iteration is merely a standard fixed point iteration for
//   a preconditioned system, and so this accelerator is more generally
//   applicable to fixed point iterations.

//   This code is a straightforward translation of the original Fortran 95
//   implementation into C.

//   [1] N.N.Carlson and K.Miller, "Design and application of a gradient-
//       weighted moving finite element code I: in one dimension", SIAM J.
//       Sci. Comput;, 19 (1998), pp. 728-765.  See section 9.

//   ************************************************************************

//   Copyright (c) 2009  Neil N. Carlson

//   Permission is hereby granted, free of charge, to any person obtaining a
//   copy of this software and associated documentation files (the "Software"),
//   to deal in the Software without restriction, including without limitation
//   the rights to use, copy, modify, merge, publish, distribute, sublicense,
//   and/or sell copies of the Software, and to permit persons to whom the
//   Software is furnished to do so, subject to the following conditions:

//   The above copyright notice and this permission notice shall be included
//   in all copies or substantial portions of the Software.

//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//   DEALINGS IN THE SOFTWARE.

#include "NonlinarKrylovAccelerator.hh"

///////////////////////////////////////////////////////////////////////

namespace Amanzi {

NonlinearKrylovAccelerator::NonlinearKrylovAccelerator (int mvec, double vtol, const Amanzi::TreeVector& initvec)
{
  mvec_ = std::max(mvec,1); // we cannot have mvec_ < 1
  vtol_ = vtol;

  v_ = new Teuchos::RCP<Amanzi::TreeVector> [mvec_+1];
  w_ = new Teuchos::RCP<Amanzi::TreeVector> [mvec_+1];

  for (int i=0; i<mvec_+1; i++)
  {
    v_[i] = Teuchos::rcp( new  Amanzi::TreeVector( std::string("v_"), initvec ) );
    w_[i] = Teuchos::rcp( new  Amanzi::TreeVector( std::string("w_"), initvec ) );
  }

  h_ =  new double* [mvec_+1];
  for (int j = 0; j < mvec_+1; j++)
  {
    h_[j] = new double[mvec_+1];
  }

  next_v_ = new int [mvec_ + 1];
  prev_v_ = new int [mvec_ + 1];

};

///////////////////////////////////////////////////////////////////////

NonlinearKrylovAccelerator::~NonlinearKrylovAccelerator()
{
  
  delete [] v_;
  delete [] w_;
  for (int j=0;j<mvec_+1; ++j) {
    delete [] h_[j];
  }
  delete [] h_;
  delete [] next_v_;
  delete [] prev_v_;
};

///////////////////////////////////////////////////////////////////////

void NonlinearKrylovAccelerator::nka_relax()
{
  int new_v;

  if (pending_)
  {
    // Drop the initial slot where the pending_ vectors are stored.
    assert(first_v_ >= 0);
    new_v = first_v_;
    first_v_ = next_v_[first_v_];
    if (first_v_ == NKAEOL) {
      last_v_ = NKAEOL;
    } else {
      prev_v_[first_v_] = NKAEOL;
    }
    // Update the free storage list.
    next_v_[new_v] = free_v_;
    free_v_ = new_v;
    pending_ = NKAFALSE;
  }
}

///////////////////////////////////////////////////////////////////////

void NonlinearKrylovAccelerator::nka_restart ()
{
  int k;

  // No vectors are stored.
  first_v_    = NKAEOL;
  last_v_     = NKAEOL;
  subspace_ = NKAFALSE;
  pending_  = NKAFALSE;

  // Initialize the free storage linked list.
  free_v_ = 0;
  for (k = 0; k < mvec_; k++) {
    next_v_[k] = k + 1;
  }
  next_v_[mvec_] = NKAEOL;
}

///////////////////////////////////////////////////////////////////////

void NonlinearKrylovAccelerator::nka_correction (const Amanzi::TreeVector& f, Amanzi::TreeVector &dir)
{
  int i, j, k, nvec, new_v;
  double s, hkk, hkj, cj;
  double  *hk, *hj, *c;



  Teuchos::RCP<Amanzi::TreeVector> vp, wp;
  Teuchos::RCP<Amanzi::TreeVector> ff = Teuchos::rcp(new Amanzi::TreeVector(std::string("ff"), f));
  *ff = f;

  // UPDATE THE ACCELERATION SUBSPACE_

  if (pending_) {
    // next_v_ function difference w_1
    wp = w_[first_v_];
    
    wp->Update(-1.0, *ff, 1.0);

    if ( wp->Dot(*wp, &s) != 0) throw "Dot problem(1)!";
    s = sqrt(s);
    
    // If the function difference is 0, we can't update the subspace_ with
    // this data; so we toss it out and continue.  In this case it is likely
    // that the outer iterative solution procedure has gone badly awry
    // (unless the function value is itself 0), and we merely want to do
    // something reasonable here and hope that situation is detected on the
    // outside.
    if (s == 0.0) {
      // nka_relax sets pending_ to NKAFALSE
      nka_relax();
    }
  }

  if (pending_)  {
    // Normalize w_1 and apply same factor to v_1.
    vp = v_[first_v_];

    vp->Scale(1.0/s);
    wp->Scale(1.0/s);


    // Update H.
    for (k = next_v_[first_v_]; k != NKAEOL; k = next_v_[k]) {
      // h[first_v_][k] = wp->innerProduct(*w[k]);
      if (wp->Dot(*w_[k], &h_[first_v_][k]) != 0) throw "Dot problem(2)!";
    }


    //  CHOLESKI FACTORIZATION OF H = W^t W
    //  original matrix kept in the upper triangle (implicit unit diagonal)
    //  lower triangle holds the factorization

    // Trivial initial factorization stage.
    nvec = 1;
    h_[first_v_][first_v_] = 1.0;

    for (k = next_v_[first_v_]; k != NKAEOL; k = next_v_[k]) {

      // Maintain at most MVEC_ vectors.
      if (++nvec > mvec_) {
        // Drop the last_v_ vector and update the free storage list.
        assert(last_v_ == k);
        next_v_[last_v_] = free_v_;
        free_v_ = k;
        last_v_ = prev_v_[k];
        next_v_[last_v_] = NKAEOL;
        break;
      }
      // Single stage of Choleski factorization.
      hk = h_[k];   // row k of H
      hkk = 1.0;
      for (j = first_v_; j != k; j = next_v_[j]) {
        hj = h_[j];   // row j of H
        hkj = hj[k];
        for (i = first_v_; i != j; i = next_v_[i]) {
          hkj -= hk[i] * hj[i];
        }
        hkj /= hj[j];
        hk[j] = hkj;
        hkk -= hkj*hkj;
      }
      if (hkk > vtol_*vtol_) {
        hk[k] = sqrt(hkk);
      } else {
        
        // The current w nearly lies in the span of the previous vectors:
        // Drop this vector,
        assert(prev_v_[k] != NKAEOL);
        next_v_[prev_v_[k]] = next_v_[k];
        if (next_v_[k] == NKAEOL) {
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

    assert(first_v_ != NKAEOL);
    subspace_ = NKATRUE; // the acceleration subspace_ isn't empty
  }


  

  //  ACCELERATED CORRECTION

  // Locate storage for the new vectors.
  assert(free_v_ != NKAEOL);
  new_v = free_v_;
  free_v_ = next_v_[free_v_];

  // Save the original f for the next_v_ call.
  *w_[new_v] = *ff;


  if (subspace_) {
    c = new double[mvec_ + 1];

    assert(c != NULL);
    // Project f onto the span of the w vectors:
    // forward substitution
    for (j = first_v_; j != NKAEOL; j = next_v_[j]) {
      // cj = (*ff).innerProduct(*w[j]);
      if (ff->Dot(*w_[j], &cj) != 0) throw "Dot problem(3)!";

      for (i = first_v_; i != j; i = next_v_[i]) {
        cj -= h_[j][i] * c[i];
      }
      c[j] = cj / h_[j][j];
    }
    // backward substitution
    for (j = last_v_; j != NKAEOL; j = prev_v_[j]) {
      cj = c[j];
      for (i = last_v_; i != j; i = prev_v_[i]) {
        cj -= h_[i][j] * c[i];
      }
      c[j] = cj / h_[j][j];
    }

    // The accelerated correction
    for (k = first_v_; k != NKAEOL; k = next_v_[k]) {
      wp = w_[k];
      vp = v_[k];

      ff->Update(c[k], *vp, -c[k], *wp, 1.0);
    }
    delete c;
  }

  // Save the accelerated correction for the next_v_ call.
  v_[new_v] = ff;


  // Prepend the new vectors to the list.
  prev_v_[new_v] = NKAEOL;
  next_v_[new_v] = first_v_;
  if (first_v_ == NKAEOL) {
    last_v_ = new_v;
  } else {
    prev_v_[first_v_] = new_v;
  }
  first_v_ = new_v;

  // The original f and accelerated correction are cached for the next_v_ call.
  pending_ = NKATRUE;

  // pass back the accelerated correction vector
  dir = *ff;

};

}
