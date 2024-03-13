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
#include "LookupTable_FEHM.hh"

namespace Amanzi {
namespace AmanziEOS {

/* ******************************************************************
* Populate internal structures
****************************************************************** */
LookupTable_FEHM::LookupTable_FEHM(Teuchos::ParameterList& plist) : LookupTable(plist)
{
  std::string filename = plist.get<std::string>("table name");
  std::string field = plist.get<std::string>("field name");
  M_ = plist.get<double>("molar weight");

  std::ifstream ifs;
  ifs.open(filename, std::ifstream::in);

  Errors::Message msg;
  msg << "\nFailed to open/read data from file: " << filename;
  if (ifs.fail()) Exceptions::amanzi_throw(msg);

  char line[100];
  ifs.getline(line, 100);
  if (ifs.fail()) Exceptions::amanzi_throw(msg);

  int nT, nP, nB, nK;
  ifs >> nT >> nP >> nB;
  ifs.getline(line, 100);
  if (ifs.fail()) Exceptions::amanzi_throw(msg);

  // primary variables are T and p
  ifs >> scaleT_ >> shiftT_;
  ifs.getline(line, 100);
  if (ifs.fail()) Exceptions::amanzi_throw(msg);

  ifs >> scaleP_ >> shiftP_;
  ifs.getline(line, 100);
  if (ifs.fail()) Exceptions::amanzi_throw(msg);

  ifs >> nK;
  ifs.getline(line, 100);
  if (ifs.fail()) Exceptions::amanzi_throw(msg);

  // T and p axis
  double value;
  ifs.getline(line, 100);
  axisT_.resize(nT);
  for (int i = 0; i < nT; ++i) {
    ifs >> value;
    axisT_[i] = scaleT_ * value + shiftT_ + 273.15;
  }
  ifs.getline(line, 100);

  ifs.getline(line, 100);
  axisP_.resize(nP);
  for (int i = 0; i < nP; ++i) {
    ifs >> value;
    axisP_[i] = (scaleP_ * value + shiftP_) * 1.0e+6;
  }
  ifs.getline(line, 100);

  // skip lines
  for (int i = 0; i < 10; ++i) ifs.getline(line, 100);

  // allocate memory for table
  F_.resize(nP);
  for (int i = 0; i < nP; ++i) F_[i].resize(nT);

  // We assume the following order: density, enthalpy, viscosity
  std::string label;
  if (field == "density") {
    label = ReadBlock_(ifs, nP, nT);
  } else if (field == "viscosity") {
    for (int i = 0; i < 3; ++i) ReadBlock_(ifs, nP, nT);
  } else if (field == "internal_energy") {
    ReadBlock_(ifs, nP, nT);
    auto density = F_;

    ReadBlock_(ifs, nP, nT);
    for (int i = 0; i < nP; ++i)
      for (int j = 0; j < nT; ++j) F_[i][j] = F_[i][j] * 1000.0 - axisP_[i] / density[i][j];
  }
}


/* ******************************************************************
* Function evaluation
****************************************************************** */
std::string
LookupTable_FEHM::ReadBlock_(std::ifstream& ifs, int nP, int nT)
{
  char line[100];
  std::string label;

  ifs.getline(line, 100);
  label = line;
  for (int i = 0; i < nP; ++i)
    for (int j = 0; j < nT; ++j) ifs >> F_[i][j];
  ifs.getline(line, 100);

  // -- derivative wrt temperature
  double value;
  ifs.getline(line, 100);
  for (int i = 0; i < nT * nP; ++i) ifs >> value;
  ifs.getline(line, 100);

  // -- derivative wrt pressure
  ifs.getline(line, 100);
  for (int i = 0; i < nT * nP; ++i) ifs >> value;
  ifs.getline(line, 100);

  return label;
}

} // namespace AmanziEOS
} // namespace Amanzi
