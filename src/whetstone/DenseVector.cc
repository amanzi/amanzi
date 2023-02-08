/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Replacement of linear vectors. It may go away after code upgrade.
*/

#include "DenseVector.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructors
****************************************************************** */
DenseVector::DenseVector(int mrow) : m_(mrow), mem_(mrow)
{
  data_ = new double[mem_];
}


DenseVector::DenseVector(int mrow, double* data)
{
  m_ = mrow;
  mem_ = mrow;
  data_ = new double[mem_];
  for (int i = 0; i < m_; i++) data_[i] = data[i];
}


DenseVector::DenseVector(const DenseVector& B)
{
  m_ = B.NumRows();
  mem_ = m_;
  data_ = NULL;
  if (m_ > 0) {
    data_ = new double[m_];
    const double* dataB = B.Values();
    for (int i = 0; i < m_; i++) data_[i] = dataB[i];
  }
}


DenseVector::DenseVector(const AmanziMesh::Double_List& B)
{
  m_ = B.size();
  mem_ = m_;
  data_ = NULL;
  if (m_ > 0) {
    data_ = new double[m_];
    for (int i = 0; i < m_; i++) data_[i] = B[i];
  }
}


/* ******************************************************************
* Smart memory management preserving data
****************************************************************** */
void
DenseVector::Reshape(int mrow, double val)
{
  int m0 = m_;
  m_ = mrow;

  if (mem_ < m_) {
    double* data_tmp = new double[m_];
    if (data_ != NULL) {
      for (int i = 0; i < mem_; i++) data_tmp[i] = data_[i];
      delete[] data_;
    }
    mem_ = m_;
    data_ = data_tmp;
  }

  if (val != 0.0) {
    for (int i = m0; i < m_; ++i) data_[i] = val; 
  }
}


/* ******************************************************************
* Assignment operator
****************************************************************** */
DenseVector&
DenseVector::operator=(const DenseVector& B)
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


/* ******************************************************************
* Vector based initialization. The size of the vector is not changed!
****************************************************************** */
void
DenseVector::PutVector(const DenseVector& v, double val)
{
  int mmin = std::min(m_, v.NumRows());
  const double* vdata = v.Values();
  for (int i = 0; i < mmin; ++i) data_[i] = vdata[i];
  for (int i = mmin; i < m_; ++i) data_[i] = val;
}

} // namespace WhetStone
} // namespace Amanzi
