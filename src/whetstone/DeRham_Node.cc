/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Derham complex: mimetic inner products for nodal DOFs.
*/

#include "Mesh.hh"

#include "DeRham_Node.hh"
#include "WhetStone_typedefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Efficient implementation is possible in 2D. Hence, we fork the code.
* Non-symmetric tensor is not yet used.
****************************************************************** */
int DeRham_Node::L2consistency(int c, const Tensor& T,
                               DenseMatrix& N, DenseMatrix& Mc, bool symmetry)
{
  int ok;
  return ok;
}

}  // namespace WhetStone
}  // namespace Amanzi

