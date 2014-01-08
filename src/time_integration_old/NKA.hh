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



#ifndef __NKA_H__
#define __NKA_H__

#include "NOX.H"
#include "NOX_Abstract_Vector.H"

#define NKATRUE 1
#define NKAFALSE 0
#define NKAEOL -1   // end-of-list marker 

class nka {
 
public:

  nka (int, double, const NOX::Abstract::Vector&);
  ~nka();
  void nka_relax();
  void nka_restart ();
  void nka_correction (NOX::Abstract::Vector&, const Teuchos::RCP<NOX::Abstract::Vector>, 
		       Teuchos::RCP<NOX::Abstract::Vector> = Teuchos::null);


private:

  int subspace;       // boolean: a nonempty subspace 
  int pending;        // contains pending vectors -- boolean 
  int mvec;           // maximum number of subspace vectors 
  double vtol;        // vector drop tolerance 

  Teuchos::RCP<NOX::Abstract::Vector> *v;  // subspace storage
  
  Teuchos::RCP<NOX::Abstract::Vector> *w;  // function difference vectors 
  
  double **h;   // matrix of w vector inner products */

  // Linked-list organization of the vector storage. 
  int first_v;  // index of first_v subspace vector 
  int last_v;   // index of last_v subspace vector 
  int free_v;   // index of the initial vector in free storage linked list 
  int *next_v;  // next_v index link field 
  int *prev_v;  // previous index link field in doubly-linked subspace v 

};


#endif
