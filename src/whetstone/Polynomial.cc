/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Operations with polynomials.
*/

#include "Polynomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructor
****************************************************************** */
Polynomial::Polynomial(int d, int order) : d_(d), order_(order)
{
  size_ = 0;
  coefs_.resize(order_ + 1);
  for (int i = 0; i <= order_; ++i) {
    coefs_[i] = Monomial(d_, i);
    size_ += coefs_[i].size();
  }
}


/* ******************************************************************
* Iterators: 2D algorithm
****************************************************************** */
void Polynomial::IteratorReset()
{
  it_k_ = 0;
  it_m_ = 0;
  it_index_[0] = 0;
  it_index_[1] = 0;
  it_index_[2] = 0;
}


void Polynomial::IteratorNext() 
{
  if (it_index_[0] == 0) {
    it_k_++;
    it_m_ = 0;
    it_index_[0] = it_k_;
    it_index_[1] = 0;
    it_index_[2] = 0;
  } else {
    it_m_++;
    it_index_[0]--;
    it_index_[1] = it_k_ - it_index_[0];
  }
}


/* ******************************************************************
* Global position of monomial defined by multi index. 2D algorithm.
****************************************************************** */
int Polynomial::MonomialPosition(const int* multi_index) const
{
  return multi_index[1];
}

}  // namespace WhetStone
}  // namespace Amanzi


