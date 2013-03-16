/*
This is the mimetic discretization component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Release name: ara-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
*/

#include <cmath>
#include <vector>

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Mesh.hh"
#include "Point.hh"
#include "errors.hh"

#include "mfd3d.hh"
#include "tensor.hh"


namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Consistency condition for stifness matrix in geomechanics. 
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int MFD3D::H1consistencyElasticity(int cell, const Tensor& T,
                                   Teuchos::SerialDenseMatrix<int, double>& N,
                                   Teuchos::SerialDenseMatrix<int, double>& Ac)
{
  return WHETSTONE_ELEMENTAL_MATRIX_OK;  // (lipnikov@lanl.gov)
}

}  // namespace WhetStone
}  // namespace Amanzi



