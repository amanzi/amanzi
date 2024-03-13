/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef AMANZI_MPC_UTILS_HH
#define AMANZI_MPC_UTILS_HH

#include <iostream>
#include "stdlib.h"
#include "math.h"

std::vector<std::pair<double, double>>
getObservations(const std::string& filename, int skip)
{
  std::vector<std::pair<double, double>> data;

  std::ifstream ifs(filename.c_str(), std::ios::in);

  // skip two lines
  char line[80];
  ifs.getline(line, 80);
  ifs.getline(line, 80);

  do {
    std::string word;
    ifs >> word;
    if (ifs.eof()) break;

    for (int i = 0; i < skip; ++i) ifs >> word;
    double val1 = std::atof(word.c_str());
    ifs >> word;
    double val2 = std::atof(word.c_str());
    data.push_back(std::make_pair(val1, val2));
    if (ifs.eof()) break;
  } while (true);

  return data;
}

#endif
