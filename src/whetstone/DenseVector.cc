/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Replacement of linear vectors. It may go away after code upgrade.
*/

#include "DenseVector.hh"

namespace Amanzi {
namespace WhetStone {

#if 0 
/* ******************************************************************
 * Assignment operator
 ****************************************************************** */
template<class MEMSPACE> 
DenseVector<MEMSPACE>&
DenseVector<MEMSPACE>::operator=(const DenseVector<MEMSPACE>& B)
{
  if (this != &B) {
    if (mem_ < B.m_) {
      if (data_ != NULL) { delete[] data_; }
      data_ = new double[B.m_];
      mem_ = B.m_;
    }
    m_ = B.m_;
    const double* b = B.Values();
    for (int i = 0; i < m_; ++i) data_[i] = b[i];
  }
  return (*this);
}
#endif 
} // namespace WhetStone
} // namespace Amanzi
