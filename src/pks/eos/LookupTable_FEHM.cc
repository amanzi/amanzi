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

#include "dbc.hh"
#include "errors.hh"
#include "LookupTable_FEHM.hh"

namespace Amanzi {
namespace AmanziEOS {

/* ******************************************************************
* Populate internal structures
****************************************************************** */
LookupTable_FEHM::LookupTable_FEHM(Teuchos::ParameterList& plist)
  : LookupTable(plist)
{
  std::string filename = plist.get<std::string>("table name");
  field_ = plist.get<std::string>("field name");

  std::ifstream ifs;
  ifs.open(filename, std::ifstream::in);

  Errors::Message msg;
  msg << "\nFailed to open/read data from file: " << filename;
  if (ifs.fail() ) Exceptions::amanzi_throw(msg);

  char line[200];
  ifs.getline(line, 100);
  if (ifs.fail() ) Exceptions::amanzi_throw(msg);

  int nT, nP, nB, nK;
  ifs >> nT >> nP >> nB;
  ifs.getline(line, 100);
  if (ifs.fail() ) Exceptions::amanzi_throw(msg);

  // primary variables are T and p
  ifs >> scaleT_ >> shiftT_;
  ifs.getline(line, 100);
  if (ifs.fail() ) Exceptions::amanzi_throw(msg);

  ifs >> scaleP_ >> shiftP_;
  ifs.getline(line, 100);
  if (ifs.fail() ) Exceptions::amanzi_throw(msg);

  ifs >> nK;
  ifs.getline(line, 100);
  if (ifs.fail() ) Exceptions::amanzi_throw(msg);

  // T and p axis
  double value;
  ifs.getline(line, 100);
  axisT_.resize(nT);
  for (int i = 0; i < nT; ++i) {
    ifs >> value;
    axisT_[i] = scaleT_ * value + shiftT_;
  }
  ifs.getline(line, 100);

  ifs.getline(line, 100);
  axisP_.resize(nP);
  for (int i = 0; i < nP; ++i) {
    ifs >> value;
    axisP_[i] = scaleP_ * value + shiftP_;
  }
  ifs.getline(line, 100);

  AMANZI_ASSERT(axisT_[0] > 273.0 || axisP_[nP - 1] > 1.0e+5); // scaling is missing

  // skip lines
  for (int i = 0; i < 10; ++i) ifs.getline(line, 100);

  // allocate memory for table
  F_.resize(nP);
  for (int i = 0; i < nP; ++i) F_[i].resize(nT);
  dFdP_.resize(nP);
  for (int i = 0; i < nP; ++i) dFdP_[i].resize(nT);
  dFdT_.resize(nP);
  for (int i = 0; i < nP; ++i) dFdT_[i].resize(nT);

  // We assume the following order: density, enthalpy, viscosity
  std::string label;
  if (field_ == "density") {
    label = ReadBlock_(ifs, nP, nT, true);
    ReadBlock_(ifs, nP, nT, false);
    ReadBlock_(ifs, nP, nT, false);
  } else if (field_ == "viscosity") {
    ReadBlock_(ifs, nP, nT, false);
    ReadBlock_(ifs, nP, nT, false);
    ReadBlock_(ifs, nP, nT, true);
  } else if (field_ == "internal_energy") {
    ReadBlock_(ifs, nP, nT, true);
    auto density = F_;
    auto DdensityDT = dFdT_;
    auto DdensityDP = dFdP_;

    ReadBlock_(ifs, nP, nT, true);
    for (int i = 0; i < nP; ++i) {
      for (int j = 0; j < nT; ++j) {
        F_[i][j] = F_[i][j] * 1000.0 - axisP_[i] / density[i][j];

        double tmp = axisP_[i] / std::pow(density[i][j], 2);
        dFdT_[i][j] = dFdT_[i][j] * 1000.0 + tmp * DdensityDT[i][j];
        dFdP_[i][j] = dFdP_[i][j] * 1000.0 + tmp * DdensityDP[i][j] - 1.0 / density[i][j];
      }
    }

    ReadBlock_(ifs, nP, nT, false);
  }

  // Read saturation line point
  ifs.getline(line, 100);
  while (strncmp(line, " Saturation line closeness", 26) != 0) ifs.getline(line, 200);

  int n, nS;
  for (int i = 0; i < nP * nT + 1; ++i) ifs >> nS;
  satP_.resize(nS);
  satT_.resize(nS);

  ifs.getline(line, 100);
  ifs.getline(line, 100);

  for (int i = 0; i < nS; ++i) {
    ifs >> satP_[i];
    ifs >> satT_[i];
    ifs >> value;
    ifs >> n;
    ifs >> n;
  }

  // Read physical values on liquid size
  ReadBlockSat_(ifs, nS, satFl_);
  ReadBlockSat_(ifs, nS, satFg_);

  // Optimization: downward map
  ComputeDownwardMap_(nP, nT, nS);
}


/* ******************************************************************
* Function evaluation
****************************************************************** */
double
LookupTable_FEHM::Function(double T, double p, int* ierr)
{
  return Function_(T, p, ierr, F_);
}


/* ******************************************************************
* Derivative evaluation
****************************************************************** */
double
LookupTable_FEHM::DFunctionDT(double T, double p, int* ierr)
{
  return Function_(T, p, ierr, dFdT_);
}


/* ******************************************************************
* Derivative evaluation
****************************************************************** */
double
LookupTable_FEHM::DFunctionDp(double T, double p, int* ierr)
{
  return Function_(T, p, ierr, dFdP_);
}


/* ******************************************************************
* Generic function evaluation
****************************************************************** */
double
LookupTable_FEHM::Function_(double T,
                            double p,
                            int* ierr,
                            const std::vector<std::vector<double>>& F)
{
  int ip, jt, n, m(0);
  *ierr = FindBox_(T, p, &ip, &jt);
  if (*ierr > 0) return 0.0;

  double f00, f01, f10, f11, g00, g11;
  f00 = g00 = F[ip - 1][jt - 1];
  f01 = F[ip - 1][jt];
  f10 = F[ip][jt - 1];
  f11 = g11 = F[ip][jt];

  n = map_[ip - 1][jt - 1];
  if (n >= 0) {
    f00 = satFl_[n][0];
    g00 = satFg_[n][0];
    m++;
  }

  n = map_[ip - 1][jt];
  if (n >= 0) {
    f01 = satFg_[n][0];
    m++;
  }

  n = map_[ip][jt - 1];
  if (n >= 0) {
    f10 = satFl_[n][0];
    m++;
  }

  n = map_[ip][jt];
  if (n >= 0) {
    f11 = satFl_[n][0];
    g11 = satFg_[n][0];
    m++;
  }

  // bilinear and linear interpolations
  double a, b, dp, dT, val;
  if (m < 2) {
    a = (axisP_[ip] - p) / (axisP_[ip] - axisP_[ip - 1]);
    b = (axisT_[jt] - T) / (axisT_[jt] - axisT_[jt - 1]);

    val = (1.0 - a) * ((1.0 - b) * f11 + b * f10) + a * ((1.0 - b) * f01 + b * f00);
  } else if (m == 2) {
    dp = axisP_[ip] - axisP_[ip - 1];
    dT = axisT_[jt] - axisT_[jt - 1];

    a = (p - axisP_[ip - 1]) / dp;
    b = (T - axisT_[jt - 1]) / dT;

    if (a + b <= 1.0) {
      val = a * f11 + b * f00 + (1.0 - a - b) * f10;
    } else {
      val = (1.0 - a) * g00 + (1.0 - b) * g11 + (a + b - 1.0) * f01;
    }
  } else {
    AMANZI_ASSERT(false);
    val = 0.0;
  }
  return val;
}


/* ******************************************************************
* Location
****************************************************************** */
int
LookupTable_FEHM::Location(double T, double p, int* ierr)
{
  int ip, jt, m(0);
  *ierr = FindBox_(T, p, &ip, &jt);
  if (*ierr > 0) return EOS_TABLE_UNKNOWN;

  if (map_[ip - 1][jt - 1] >= 0) m++;
  if (map_[ip][jt] >= 0) m++;
  if (m < 2) return map_[ip - 1][jt - 1];

  double dp, dT, a, b, tol;
  dp = axisP_[ip] - axisP_[ip - 1];
  dT = axisT_[jt] - axisT_[jt - 1];

  a = (p - axisP_[ip - 1]) / dp;
  b = (T - axisT_[jt - 1]) / dT;
  tol = b - a;

  if (tol > EOS_TABLE_TOL) return EOS_TABLE_GAS;
  else if (tol < -EOS_TABLE_TOL) return EOS_TABLE_LIQUID;

  return EOS_TABLE_SATURATION_CURVE;
}


/* ******************************************************************
* Reading field block
****************************************************************** */
std::string
LookupTable_FEHM::ReadBlock_(std::ifstream& ifs, int nP, int nT, bool flag)
{
  double value;
  char line[100];
  std::string label;

  ifs.getline(line, 100);
  label = line;
  if (flag) {
    for (int i = 0; i < nP; ++i)
      for (int j = 0; j < nT; ++j) ifs >> F_[i][j];
  } else {
    for (int i = 0; i < nT * nP; ++i) ifs >> value;
  }
  ifs.getline(line, 100);

  // -- derivative wrt temperature
  ifs.getline(line, 100);
  if (flag) {
    for (int i = 0; i < nP; ++i)
      for (int j = 0; j < nT; ++j) ifs >> dFdT_[i][j];
  } else {
    for (int i = 0; i < nT * nP; ++i) ifs >> value;
  }
  ifs.getline(line, 100);

  // -- derivative wrt pressure
  ifs.getline(line, 100);
  if (flag) {
    for (int i = 0; i < nP; ++i) {
      for (int j = 0; j < nT; ++j) {
        ifs >> dFdP_[i][j];
        dFdP_[i][j] /= scaleP_;
      }
    }
  } else {
    for (int i = 0; i < nT * nP; ++i) ifs >> value;
  }
  ifs.getline(line, 100);

  return label;
}


/* ******************************************************************
* Reading saturation block
****************************************************************** */
void
LookupTable_FEHM::ReadBlockSat_(std::ifstream& ifs, int nS, std::vector<std::vector<double>>& satF)
{
  char line[200];
  double data[9];

  satF.resize(nS);

  while (strncmp(line, " Density", 8) != 0) ifs.getline(line, 200);
  for (int i = 0; i < nS; ++i) {
    for (int k = 0; k < 9; ++k) ifs >> data[k];

    satF[i].resize(3);
    if (field_ == "density") {
      for (int k = 0; k < 3; ++k) satF[i][k] = data[k];
    } else if (field_ == "viscosity") {
      for (int k = 0; k < 3; ++k) satF[i][k] = data[6 + k];
    } else if (field_ == "internal_energy") {
      for (int k = 0; k < 3; ++k) satF[i][k] = data[3 + k];
    }
  }
}


/* ******************************************************************
* Downward map: table node -> saturation point or -1.
****************************************************************** */
void
LookupTable_FEHM::ComputeDownwardMap_(int nP, int nT, int nS)
{
  map_.resize(nP);
  for (int i = 0; i < nP; ++i) {
    map_[i].resize(nT);
    for (int j = 0; j < nT; ++j) map_[i][j] = EOS_TABLE_UNKNOWN;
  }

  int i0, j0;
  double p, T;
  for (int n = 0; n < nS; ++n) {
    p = scaleP_ * satP_[n] + shiftP_;
    T = scaleT_ * satT_[n] + shiftT_;

    i0 = j0 = -1;
    for (int i = 0; i < nT; ++i) {
      if (std::fabs(T - axisT_[i]) / scaleT_ < EOS_TABLE_TOL) {
        i0 = i;
        break;
      }
    }
    for (int i = 0; i < nP; ++i) {
      if (std::fabs(p - axisP_[i]) / scaleP_ < EOS_TABLE_TOL) {
        j0 = i;
        break;
      }
    }
    if (i0 < 0 || j0 < 0) AMANZI_ASSERT(false);

    map_[j0][i0] = n;
  }

  // mark gases
  int phase = EOS_TABLE_LIQUID;
  for (int i = 0; i < nP; ++i) {
    for (int j = 0; j < nT; ++j) {
      if (map_[i][j] == EOS_TABLE_UNKNOWN) {
        map_[i][j] = phase;
      } else if (map_[i][j] >= 0) {
        phase = EOS_TABLE_GAS;
      }
    }
  }
}

} // namespace AmanziEOS
} // namespace Amanzi
