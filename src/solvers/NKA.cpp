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


//   The following code was adapted from Neil Carlson's  NLKAIN code that is on 
//   Sourceforge, see http://sourceforge.net/projects/nlkain/   
//   Markus Berndt, CCS-2




#include "NKA.H"

///////////////////////////////////////////////////////////////////////

nka::nka (int mvec_, double vtol_, const NOX::Abstract::Vector &initvec)
{
  mvec = max(mvec_,1); // we cannot have mvex < 1
  vtol = vtol_; 
  
  v = new Teuchos::RCP<NOX::Abstract::Vector> [mvec+1];
  w = new Teuchos::RCP<NOX::Abstract::Vector> [mvec+1];
  
  for (int i=0; i<mvec+1; i++) 
    {
      v[i] = initvec.clone(NOX::ShapeCopy);
      w[i] = initvec.clone(NOX::ShapeCopy);
    }
  
  h =  new double* [mvec+1];
  for (int j = 0; j < mvec+1; j++) 
    {
      h[j] = new double[mvec+1];
    }
  
  next_v = new int [mvec + 1];
  prev_v = new int [mvec + 1];

};

///////////////////////////////////////////////////////////////////////

nka::~nka()
{
  delete [] v;
  delete [] w;
  for (int j=0;j<mvec+1; ++j) {
    delete [] h[j];
  }
  delete [] h;
  delete [] next_v;
  delete [] prev_v;
};

///////////////////////////////////////////////////////////////////////

void nka::nka_relax()
{
  int new_v;
  
  if (pending) 
    {
      // Drop the initial slot where the pending vectors are stored.
      assert(first_v >= 0);
      new_v = first_v;
      first_v = next_v[first_v];
      if (first_v == NKAEOL) 
	{
	  last_v = NKAEOL;
	} 
      else 
	{
	  prev_v[first_v] = NKAEOL;
	}
      // Update the free storage list.
      next_v[new_v] = free_v;
      free_v = new_v;
      pending = NKAFALSE;
    }    
}

///////////////////////////////////////////////////////////////////////

void nka::nka_restart ()
{
  int k;
  
  // No vectors are stored. 
  first_v    = NKAEOL;
  last_v     = NKAEOL;
  subspace = NKAFALSE;
  pending  = NKAFALSE;
  
  // Initialize the free storage linked list. 
  free_v = 0;
  for (k = 0; k < mvec; k++) 
    {
      next_v[k] = k + 1;
    }
  next_v[mvec] = NKAEOL;
}

///////////////////////////////////////////////////////////////////////

void nka::nka_correction (NOX::Abstract::Vector &dir, 
			  const Teuchos::RCP<NOX::Abstract::Vector> f,
			  double damp)
{
  int i, j, k, nvec, new_v;
  double s, hkk, hkj, cj;
  double  *hk, *hj, *c;
  
  Teuchos::RCP<NOX::Abstract::Vector> vp, wp;
  Teuchos::RCP<NOX::Abstract::Vector> ff = (*f.get()).clone(NOX::DeepCopy);

  // UPDATE THE ACCELERATION SUBSPACE

  if (pending) 
    {
      // next_v function difference w_1 
      wp = w[first_v];
      


      wp->update(-1.0, *ff.get(), 1.0);
      
      s = wp->innerProduct(*wp);
      s = sqrt(s);
      
      // If the function difference is 0, we can't update the subspace with
      // this data; so we toss it out and continue.  In this case it is likely
      // that the outer iterative solution procedure has gone badly awry
      // (unless the function value is itself 0), and we merely want to do
      // something reasonable here and hope that situation is detected on the
      // outside. 
      if (s == 0.0) 
	{
	  // nka_relax sets pending to NKAFALSE
	  nka_relax();
	}
    }
  
  if (pending)
    {
      // Normalize w_1 and apply same factor to v_1. 
      vp = v[first_v];
     
      // first damp if requested, the user will supply a damping factor != 1.0, if 
      // the previously suggested NKA update was damped
      if (damp != 1.0) vp->scale(damp);
 
      vp->scale(1.0/s);
      wp->scale(1.0/s);
      
      
      // Update H. 
      for (k = next_v[first_v]; k != NKAEOL; k = next_v[k])
	{
	  h[first_v][k] = wp->innerProduct(*w[k]);
	}
      
      
      //  CHOLESKI FACTORIZATION OF H = W^t W
      //  original matrix kept in the upper triangle (implicit unit diagonal)
      //  lower triangle holds the factorization
      
      // Trivial initial factorization stage. 
      nvec = 1;
      h[first_v][first_v] = 1.0;
      
      for (k = next_v[first_v]; k != NKAEOL; k = next_v[k]) 
	{
	  
	  // Maintain at most MVEC vectors. 
	  if (++nvec > mvec) 
	    {
	      // Drop the last_v vector and update the free storage list. 
	      assert(last_v == k);
	      next_v[last_v] = free_v;
	      free_v = k;
	      last_v = prev_v[k];
	      next_v[last_v] = NKAEOL;
	      break;
	    }
	  
	  // Single stage of Choleski factorization. 
	  hk = h[k];   // row k of H 
	  hkk = 1.0;
	  for (j = first_v; j != k; j = next_v[j]) 
	    {
	      hj = h[j];   // row j of H 
	      hkj = hj[k];
	      for (i = first_v; i != j; i = next_v[i])
		{
		  hkj -= hk[i] * hj[i];
		}
	      hkj /= hj[j];
	      hk[j] = hkj;
	      hkk -= hkj*hkj;
	    }
	  
	  if (hkk > pow(vtol,2)) 
	    {
	      hk[k] = sqrt(hkk);
	    } 
	  else  
	    {
	      // The current w nearly lies in the span of the previous vectors: 
	      // Drop this vector, 
	      assert(prev_v[k] != NKAEOL);
	      next_v[prev_v[k]] = next_v[k];
	      if (next_v[k] == NKAEOL)
		{
		  last_v = prev_v[k];
		}
	      else
		{
		  prev_v[next_v[k]] = prev_v[k];
		}
	      // update the free storage list, 
	      next_v[k] = free_v;
	      free_v = k;
	      // back-up and move on to the next_v vector. 
	      k = prev_v[k];
	      nvec--;
	    }
	}
      
      assert(first_v != NKAEOL);
      subspace = NKATRUE; // the acceleration subspace isn't empty 
      
    }
  

  //  ACCELERATED CORRECTION
  
  // Locate storage for the new vectors. 
  assert(free_v != NKAEOL);
  new_v = free_v;
  free_v = next_v[free_v];
  
  // Save the original f for the next_v call. 
  *w[new_v] = *ff;
  
  if (subspace) 
    {
      c = new double[mvec + 1];
      
      assert(c != NULL);
      // Project f onto the span of the w vectors: 
      // forward substitution 
      for (j = first_v; j != NKAEOL; j = next_v[j]) 
	{
	  cj = (*ff).innerProduct(*w[j]);
	  
	  for (i = first_v; i != j; i = next_v[i]) 
	    {
	      cj -= h[j][i] * c[i];
	    }
	  c[j] = cj / h[j][j];
	}
      // backward substitution 
      for (j = last_v; j != NKAEOL; j = prev_v[j]) 
	{
	  cj = c[j];
	  for (i = last_v; i != j; i = prev_v[i])
	    {
	      cj -= h[i][j] * c[i];
	    }
	  c[j] = cj / h[j][j];
	}
      // The accelerated correction 
      for (k = first_v; k != NKAEOL; k = next_v[k]) 
	{
	  wp = w[k];
	  vp = v[k];

	  (*ff).update(c[k], *vp, -c[k], *wp, 1.0); 
	}
      delete c;
    }
  
  // Save the accelerated correction for the next_v call. 
  v[new_v] = ff;
  
  
  // Prepend the new vectors to the list.
  prev_v[new_v] = NKAEOL;
  next_v[new_v] = first_v;
  if (first_v == NKAEOL) 
    {
      last_v = new_v;
    } 
  else 
    {
      prev_v[first_v] = new_v;
    }
  first_v = new_v;
  
  // The original f and accelerated correction are cached for the next_v call. 
  pending = NKATRUE;
  
  // pass back the accelerated correction vector
  dir = *ff;

};



