/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/


#include "UnitTest++.h"

#include "Key.hh"
#include "Formatter.hh"

using namespace Amanzi;
using namespace Amanzi::Utils;

SUITE(DEBUG_STRING_FORMATTING) {

TEST(HEADER) {
  int header_width = 20;
  int cellnum_width = 5;
  Formatter f(15, 7, header_width, cellnum_width);

  std::vector<std::string> names = { "", "s", "longer", "very_long_header_name" };
  std::vector<int> cells = { 0, 1, 10, 100, 99999};

  for (auto& name : names) {
    for (auto& c : cells) {
      auto fh = f.formatHeader(name, c);
      std::cout << name << "+" << c << " \t\t = \"" << fh << "\"" << std::endl;

      //  header_width
      CHECK_EQUAL(header_width, fh.size());
      CHECK(Keys::ends_with(fh, "): "));
    }
  }

  // GIDs > 5 sig digits are long
  for (auto& name : names) {
    int c = 100000;
    auto fh = f.formatHeader(name, c);
    std::cout << name << "+" << c << " \t\t = \"" << fh << "\"" << std::endl;

    //  header_width
    CHECK_EQUAL(header_width+1, fh.size());
    CHECK(Keys::ends_with(fh, "): "));
  }

}


TEST(NUMBERS) {
  int width = 15;
  int precision = 7;
  Formatter f(width, precision, 20, 6);

  std::vector<double> vals = {0.0, -0.0,
    1.01, -1.01, 1.01234567890123456789, -1.01234567890123456789,
    1.01e1, -1.01e1, 1.01234567890123456789e1, -1.01234567890123456789e1,
    1.01e2, -1.01e2, 1.01234567890123456789e2, -1.01234567890123456789e2,
    101325.01, -101325.01, 101325.01234567890123456789, -101325.01234567890123456789,
    1.01e-1, -1.01e-1, 1.01234567890123456789e-1, -1.01234567890123456789e-1,
    1.01e-2, -1.01e-2, 1.01234567890123456789e-2, -1.01234567890123456789e-2,
    1.01e-4, -1.01e-4, 1.01234567890123456789e-4, -1.01234567890123456789e-4,
    1.1e-7, -1.1e-7, 1.01234567890123456789e-7, -1.01234567890123456789e-7,
    1.1e7, -1.1e7, 1.01234567890123456789e7, -1.01234567890123456789e7};


  int i = 0;
  for (auto& val : vals) {
    auto vs = f.format(val);
    int n_tabs = ((i < 2) || ((i+2) % 4) > 1) ? 1 : 2;
    std::cout << std::setprecision(8) << val;
    for (int j=0; j!=n_tabs; ++j) std::cout << "\t";
    std::cout << " = \"" << vs << "\"" << std::endl;
    CHECK_EQUAL(width, vs.size());

    if (i < 26) CHECK_EQUAL('.', vs[width - precision - 1]);
    else CHECK_EQUAL('.', vs[width - precision - 4 - 1]);
    i++;
  }
}


TEST(NUMBERS2) {
  int width = 13;
  int precision = 6;
  Formatter f(width, precision, 20, 6);

  std::vector<double> vals = {0.0, -0.0,
    1.01, -1.01, 1.01234567890123456789, -1.01234567890123456789,
    1.01e1, -1.01e1, 1.01234567890123456789e1, -1.01234567890123456789e1,
    1.01e2, -1.01e2, 1.01234567890123456789e2, -1.01234567890123456789e2,
    101325.01, -101325.01, 101325.01234567890123456789, -101325.01234567890123456789,
    1.01e-1, -1.01e-1, 1.01234567890123456789e-1, -1.01234567890123456789e-1,
    1.01e-2, -1.01e-2, 1.01234567890123456789e-2, -1.01234567890123456789e-2,
    1.01e-4, -1.01e-4, 1.01234567890123456789e-4, -1.01234567890123456789e-4,
    1.1e-7, -1.1e-7, 1.01234567890123456789e-7, -1.01234567890123456789e-7,
    1.1e7, -1.1e7, 1.01234567890123456789e7, -1.01234567890123456789e7};


  int i = 0;
  for (auto& val : vals) {
    auto vs = f.format(val);
    int n_tabs = ((i < 2) || ((i+2) % 4) > 1) ? 1 : 2;
    std::cout << std::setprecision(8) << val;
    for (int j=0; j!=n_tabs; ++j) std::cout << "\t";
    std::cout << " = \"" << vs << "\"" << std::endl;
    CHECK_EQUAL(width, vs.size());

    // if (i < 26) CHECK_EQUAL('.', vs[width - precision - 1]);
    // else CHECK_EQUAL('.', vs[width - precision - 4 - 1]);
    i++;
  }
}

}
