/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  Tabulated equations of state: rows/columns are p/T
*/

#include <fstream>
#include <string>

#include "errors.hh"
#include "LookupTable.hh"

namespace Amanzi {
namespace AmanziEOS {

/* ******************************************************************
* Function evaluation
****************************************************************** */
double
LookupTable::Function(double T, double p, int* ierr)
{
  int ip, jp;
  *ierr = FindBox_(T, p, &ip, &jp);
  if (*ierr > 0) return 0.0;

  // bilinear interpolation
  double a = (axisP_[ip] - p) / (axisP_[ip] - axisP_[ip - 1]);
  double b = (axisT_[jp] - T) / (axisT_[jp] - axisT_[jp - 1]);

  double val = (1.0 - a) * ((1.0 - b) * F_[ip][jp] + b * F_[ip][jp - 1]) +
               a * ((1.0 - b) * F_[ip - 1][jp] + b * F_[ip - 1][jp - 1]);
  return val;
}


/* ******************************************************************
* Derivative in p evaluation
****************************************************************** */
double
LookupTable::DFunctionDp(double T, double p, int* ierr)
{
  int ip, jp;
  *ierr = FindBox_(T, p, &ip, &jp);
  if (*ierr > 0) return 0.0;

  // bilinear interpolation
  double a = (axisP_[ip] - p) / (axisP_[ip] - axisP_[ip - 1]);
  double b = (axisT_[jp] - T) / (axisT_[jp] - axisT_[jp - 1]);

  double v00 = DerivativeP_(ip - 1, jp - 1);
  double v01 = DerivativeP_(ip - 1, jp);
  double v10 = DerivativeP_(ip, jp - 1);
  double v11 = DerivativeP_(ip, jp);

  double val = (1.0 - a) * ((1.0 - b) * v11 + b * v10) + a * ((1.0 - b) * v01 + b * v00);
  return val;
}


/* ******************************************************************
* Derivative in T evaluation
****************************************************************** */
double
LookupTable::DFunctionDT(double T, double p, int* ierr)
{
  int ip, jp;
  *ierr = FindBox_(T, p, &ip, &jp);
  if (*ierr > 0) return 0.0;

  // bilinear interpolation
  double a = (axisP_[ip] - p) / (axisP_[ip] - axisP_[ip - 1]);
  double b = (axisT_[jp] - T) / (axisT_[jp] - axisT_[jp - 1]);

  double v00 = DerivativeT_(ip - 1, jp - 1);
  double v01 = DerivativeT_(ip - 1, jp);
  double v10 = DerivativeT_(ip, jp - 1);
  double v11 = DerivativeT_(ip, jp);

  double val = (1.0 - a) * ((1.0 - b) * v11 + b * v10) + a * ((1.0 - b) * v01 + b * v00);
  return val;
}


/* ******************************************************************
* Derivative in pressure
****************************************************************** */
double
LookupTable::DerivativeP_(int i, int j)
{
  int n(0);
  int nP = axisP_.size();
  double val(0.0);

  if (i > 0) {
    val += (F_[i][j] - F_[i - 1][j]) / (axisP_[i] - axisP_[i - 1]);
    n++;
  }
  if (i < nP - 1) {
    val += (F_[i + 1][j] - F_[i][j]) / (axisP_[i + 1] - axisP_[i]);
    n++;
  }

  return val / n;
}


/* ******************************************************************
* Derivative in temperature
****************************************************************** */
double
LookupTable::DerivativeT_(int i, int j)
{
  int n(0);
  int nT = axisT_.size();
  double val(0.0);

  if (j > 0) {
    val += (F_[i][j] - F_[i][j - 1]) / (axisT_[j] - axisT_[j - 1]);
    n++;
  }
  if (j < nT - 1) {
    val += (F_[i][j + 1] - F_[i][j]) / (axisT_[j + 1] - axisT_[j]);
    n++;
  }

  return val / n;
}


/* ******************************************************************
* Returns right-top corner of the table box
****************************************************************** */
int
LookupTable::FindBox_(double T, double p, int* ip, int* jp)
{
  int nT = axisT_.size();
  int nP = axisP_.size();

  if (T < axisT_[0] || T > axisT_[nT - 1] || p < axisP_[0] || p > axisP_[nP - 1]) return 1;

  *ip = 0;
  *jp = 0;
  for (int i = 1; i < nP; ++i) {
    if (axisP_[i] >= p) {
      *ip = i;
      break;
    }
  }
  for (int i = 1; i < nT; ++i) {
    if (axisT_[i] >= T) {
      *jp = i;
      break;
    }
  }

  return 0;
}


/* ******************************************************************
* Create error message
****************************************************************** */
std::string
LookupTable::ErrorMessage(double T, double p)
{
  std::stringstream ss;
  ss << "out of bounds values: T=" << T << " or p=" << p;
  return ss.str();
}

} // namespace AmanziEOS
} // namespace Amanzi
