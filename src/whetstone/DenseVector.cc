/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Ring algebra of linear vectors and more.
*/

#include "DenseVector.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Constructors
****************************************************************** */
DenseVector::DenseVector(int mrow) : m_(mrow), mem_(mrow) {
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


DenseVector::DenseVector(const std::vector<double>& B)
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
void DenseVector::Reshape(int mrow)
{
  m_ = mrow;

  if (mem_ < m_) {
    int mem_next_ = 2 * m_;
    double* data_tmp = new double[mem_next_];
    if (data_ != NULL) {
      // it is dangerous to not initialize reshaped data to zero
      for (int i = 0; i < mem_; ++i) data_tmp[i] = data_[i];
      for (int i = mem_; i < m_; ++i) data_tmp[i] = 0.0;
      delete [] data_;
    }
    mem_ = mem_next_;
    data_ = data_tmp;
  }
}


/* ******************************************************************
* Consequitive groups of size stride1 are copied to groups of size
* stride2. Gaps are filled with 0. Incomplete last group is lost.
****************************************************************** */
void DenseVector::Regroup(int stride1, int stride2) 
{
  if (data_ == NULL) return;

  int ngroups = m_ / stride1;
  m_ = stride2 * ngroups;

  double* data_tmp;
  if (mem_ < m_)  {
    data_tmp = new double[m_];
  } else {
    data_tmp = data_;
  }

  int stride = std::min(stride1, stride2);
  for (int n = 0; n < ngroups; ++n) {
    int i1 = n * stride1;
    int i2 = n * stride2;

    for (int i = 0; i < stride; ++i) data_tmp[i2 + i] = data_[i1 + i]; 
    for (int i = stride; i < stride2; ++i) data_tmp[i2 + i] = 0.0; 
  }

  if (mem_ < m_) {
    delete [] data_;
    data_ = data_tmp;
    mem_ = m_;
  }
}


/* ******************************************************************
* Assignment operators
****************************************************************** */
DenseVector& DenseVector::operator=(const DenseVector& other)
{
  if (this != &other) {
    if (mem_ < other.m_) {
      if (data_ != NULL) {
        delete [] data_;
      }
      data_ = new double[other.m_];
      mem_ = other.m_;
    }
    m_ = other.m_;
    const double *b = other.Values();
    for (int i = 0; i < m_; ++i) data_[i] = b[i];
  }
  return *this;
}


DenseVector& DenseVector::operator=(DenseVector&& other) noexcept
{
  if (this != &other) {
    m_ = other.m_;
    mem_ = other.mem_;
    data_ = other.data_;
    other.data_ = NULL;
  }
  return *this;
}


/* ******************************************************************
* Ring algebra
****************************************************************** */
DenseVector operator+(const DenseVector& v1, const DenseVector& v2)
{
  DenseVector tmp(v1);
  double* data = tmp.Values();
  const double* data1 = v1.Values();
  const double* data2 = v2.Values();

  for (int i = 0; i < tmp.NumRows(); ++i) {
    *data = *data1 + *data2;
    data++;
    data1++;
    data2++;
  }

  return tmp;
}


DenseVector operator-(const DenseVector& v1, const DenseVector& v2)
{
  DenseVector tmp(v1);
  double* data = tmp.Values();
  const double* data1 = v1.Values();
  const double* data2 = v2.Values();

  for (int i = 0; i < tmp.NumRows(); ++i) {
    *data = *data1 - *data2;
    data++;
    data1++;
    data2++;
  }

  return tmp;
}


/* ******************************************************************
* Vector based initialization. The size of the vector is not changed!
****************************************************************** */
void DenseVector::PutVector(const DenseVector& v, double val)
{
  int mmin = std::min(m_, v.NumRows());
  const double* vdata = v.Values();
  for (int i = 0; i < mmin; ++i) data_[i] = vdata[i];   
  for (int i = mmin; i < m_; ++i) data_[i] = val;
}

}  // namespace WhetStone
}  // namespace Amanzi

