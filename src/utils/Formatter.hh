/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/
#pragma once

#include <string>
#include <cstdarg>
#include <vector>
#include <cmath>

namespace Amanzi {
namespace Utils {

//
// Code from StackExchange Network answer
// https://stackoverflow.com/a/49812018/1322752
// Licensed under CC BY-SA 3.0
//
inline const std::string vformat(const char * const zcFormat, ...)
{
  // initialize use of the variable argument array
  std::va_list vaArgs;
  va_start(vaArgs, zcFormat);

  // reliably acquire the size
  // from a copy of the variable argument array
  // and a functionally reliable call to mock the formatting
  std::va_list vaArgsCopy;
  va_copy(vaArgsCopy, vaArgs);
  const int iLen = std::vsnprintf(NULL, 0, zcFormat, vaArgsCopy);
  va_end(vaArgsCopy);

  // return a formatted string without risking memory mismanagement
  // and without assuming any compiler or platform specific behavior
  std::vector<char> zc(iLen + 1);
  std::vsnprintf(zc.data(), zc.size(), zcFormat, vaArgs);
  va_end(vaArgs);
  return std::string(zc.data(), iLen);
}


class Formatter {
 public:
  Formatter(int width,
            int precision,
            int header_width,
            int cellnum_width);

  int getWidth() const { return width_; }
  // NOTE, this may also change precision
  void setWidth(int width);

  int getPrecision() const { return precision_; }
  void setPrecision(int prec) { precision_ = prec; }

  std::string format(double dat) const;
  std::string formatHeader(std::string header, int c) const;

private:
  int width_;
  int precision_;
  int header_width_;
  int cellnum_width_;
};

} // Utils
} // Amanzi
