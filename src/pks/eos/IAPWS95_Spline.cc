/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Maitri Dalal (mdalal@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  Splined equation of state based on the IAPWS95 model.
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "SplineCubicNotAKnot2D.hh"

#include "IAPWS95_Spline.hh"

namespace Amanzi {
namespace AmanziEOS {

/* ******************************************************************
* Populate internal structures
****************************************************************** */
IAPWS95_Spline::IAPWS95_Spline(Teuchos::ParameterList& plist)
  : IAPWS95(plist)
{
  std::string filename = plist.get<std::string>("table name");
  ReadTable_(filename);

  int ok = spline_.build(delta_, tau_, values_);
}
  

/* ******************************************************************
* Populate internal structures
****************************************************************** */
std::array<double, 6>
IAPWS95_Spline::ResidualPart(double rho, double T)
{
  double tau = TC / T;
  double delta = rho / RHOC;
  return spline_.evaluate(delta, tau); 
}


/* ******************************************************************
* Populate internal structures
****************************************************************** */
void
IAPWS95_Spline::ReadTable_(const std::string& filename)
{
  Errors::Message msg;

  std::ifstream file(filename);
  if (!file.is_open()) {
    msg << "\nFailed to open/read data from file: " << filename;
    Exceptions::amanzi_throw(msg);
  }

  std::string line;
  if (!std::getline(file, line)) {
    msg << "\nEmpty CSV file";
    Exceptions::amanzi_throw(msg);
  }

  auto headers = SplitLine_(line);
  int tau_col(0), delta_col(1), alpha_col(2);

  std::set<double> tau_set, delta_set;
  std::map<std::pair<double, double>, double> data;

  while (std::getline(file, line)) {
    if (line.empty()) continue;

    auto fields = SplitLine_(line);
    if (fields.size() <= 2) {
      msg << "\nMalformed CSV row";
      Exceptions::amanzi_throw(msg);
    }

    double tau = std::stod(fields[tau_col]);
    double delta = std::stod(fields[delta_col]);
    double alpha_r = std::stod(fields[alpha_col]);

    tau_set.insert(tau);
    delta_set.insert(delta);

    data[{tau, delta}] = alpha_r;
  }

  tau_.assign(tau_set.begin(), tau_set.end());
  delta_.assign(delta_set.begin(), delta_set.end());

  values_.resize(tau_.size(), std::vector<double>(delta_.size()));

  for (int i = 0; i < tau_.size(); ++i) {
    for (int j = 0; j < delta_.size(); ++j) {
      auto key = std::make_pair(tau_[i], delta_[j]);
      auto it = data.find(key);

      if (it == data.end()) {
        msg << "\nMissing EOS data";
        Exceptions::amanzi_throw(msg);
      }

      values_[i][j] = it->second;
    }
  }

  // verify table compatibility with cubic splines
  nx_ = tau_.size();
  ny_ = delta_.size();

  if (nx_ < 4 || ny_ < 4) {
    msg << "\nCubic interpolation requires at least 4 points in each direction.";
    Exceptions::amanzi_throw(msg);
  }

  bool ok1 = CheckStrictMonotonicity_(tau_);
  bool ok2 = CheckStrictMonotonicity_(delta_);

  if (!ok1 || !ok2) {
    msg << "\nGrid coordinate must be strictly increasing.";
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
* Helper function
****************************************************************** */
std::vector<std::string>
IAPWS95_Spline::SplitLine_(const std::string& line)
{
  std::vector<std::string> fields;
  std::string field;
  std::stringstream ss(line);

  while (std::getline(ss, field, ',')) {
    fields.push_back(field);
  }

  return fields;
}


/* ******************************************************************
* Check monotonicity of coordinates
****************************************************************** */
bool
IAPWS95_Spline::CheckStrictMonotonicity_(const std::vector<double>& x)
{
  for (int i = 1; i < x.size(); ++i) {
    if (x[i] <= x[i - 1]) return false;
  }
  return true;
}

} // namespace AmanziEOS
} // namespace Amanzi

