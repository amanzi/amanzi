/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Simple class template
*/

#include <cmath>
#include <iostream>

#include "exceptions.hh"

#include "chemistry_exception.hh"
#include "lu_solver.hh"
#include "matrix_block.hh"

namespace Amanzi {
namespace AmanziChemistry {

const double LUSolver::kSmallNumber = 1.0e-20;

LUSolver::LUSolver(void)
    : system_size_(0),
      row_interchange_(0.0),
      factored_(false) {
  pivoting_indices_.clear();
  row_scaling_.clear();
}


void LUSolver::Initialize(const int size) {
  set_system_size(size);
  pivoting_indices_.resize(system_size());
  row_scaling_.resize(system_size());
}


void LUSolver::Solve(MatrixBlock* A, std::vector<double>* b) {
  // TODO(bandre): check that the sizes of A, b and 'system_size' are
  // equal.
  Decomposition(A);
  BackSolve(A, b);
}


void LUSolver::Decomposition(MatrixBlock* A) {
  double** a = A->GetValues();
  row_interchange_ = 1.0;
  for (int i = 0; i < system_size(); ++i) {
    double big = 0.0;
    for (int j = 0; j < system_size(); ++j) {
      double temp = std::fabs(a[i][j]);
      if (temp > big) {
        big = temp;
      }
    }
    if (big == 0.0) {
      std::stringstream error_stream;
      error_stream << "LUSolver::Decomposition() : Singular matrix." << std::endl;
      Exceptions::amanzi_throw(ChemistryUnrecoverableError(error_stream.str()));
    }
    row_scaling_.at(i) = 1.0 / big;
  }

  int imax = 0;
  for (int j = 0; j < system_size(); ++j) {
    double temp;
    for (int i = 0; i < j; ++i) {
      double sum = a[i][j];
      for (int k = 0; k < i; k++) {
        sum -= a[i][k] * a[k][j];
      }
      a[i][j] = sum;
    }
    double big = 0.0;
    for (int i = j; i < system_size(); ++i) {
      double sum =  a[i][j];
      for (int k = 0; k < j; ++k) {
        sum -= a[i][k] * a[k][j];
      }
      a[i][j] = sum;
      temp = row_scaling_.at(i) * std::fabs(sum);
      if (temp >= big) {
        big = temp;
        imax = i;
      }
    }
    if (j != imax) {
      for (int k = 0; k < system_size(); ++k) {
        temp = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = temp;
      }
      row_interchange_ = -row_interchange_;
      row_scaling_.at(imax) = row_scaling_.at(j);
    }
    pivoting_indices_.at(j) = imax;
    if (a[j][j] == 0.0) {
      a[j][j] = kSmallNumber;
    }
    if (j != system_size() - 1) {
      temp = 1.0 / (a[j][j]);
      for (int i = j + 1; i < system_size(); ++i) {
        a[i][j] *= temp;
      }
    }
  }

  a = NULL;
  factored_ = true;
}


void LUSolver::BackSolve(MatrixBlock* A, std::vector<double>* b) {
  double** a = A->GetValues();
  int ii = 0;
  for (int i = 0; i < system_size(); i++) {
    int ip = pivoting_indices_.at(i);
    double sum = b->at(ip);
    b->at(ip) = b->at(i);
    if (ii != 0) {
      for (int j = ii - 1; j < i; ++j) {
        sum -= a[i][j] * b->at(j);
      }
    } else if (sum != 0.0) {
      ii = i + 1;
    }
    b->at(i) = sum;
  }

  for (int i = system_size() - 1; i >= 0; --i) {
    double sum = b->at(i);
    for (int j = i + 1; j < system_size(); ++j) {
      sum -= a[i][j] * b->at(j);
    }
    b->at(i) = sum / a[i][i];
  }
  a = NULL;
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
