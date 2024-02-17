/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include <sstream>
#include <iomanip>
#include <iostream>

#include "dbc.hh"
#include "Formatter.hh"

namespace Amanzi {
namespace Utils {

Formatter::Formatter(int width, int precision, int header_width, int cellnum_width)
  : width_(width), precision_(precision), header_width_(header_width), cellnum_width_(cellnum_width)
{
  AMANZI_ASSERT(width_ >= (precision_ + 7));
  AMANZI_ASSERT(header_width > (cellnum_width + 6));
}


void
Formatter::setWidth(int width)
{
  width_ = width;
  if (width < (precision_ + 7)) setPrecision(width - 7);
}

std::string
Formatter::format(double dat) const
{
  std::stringstream datastream;
  if (dat == 0.) {
    // how many spaces do we have to the left of the 0?
    // one for dot, one for zero
    int left_spaces = width_ - precision_ - 2;

    // AMANZI_ASSERT(width >= precision_ + 2);
    datastream << std::string(left_spaces, ' ') << "0." << std::string(precision_, ' ');
  } else {
    // can we fit in non-scientific notation?  How many digits to the left of
    // the dot?
    int mag = std::floor(std::log10(std::abs(dat)));

    // and how much room to the left of the dot do we have for those digits?
    // subtract one for dot and one for sign
    int left_of_dot = width_ - precision_ - 2;

    // use fixed format if mag > 0 and it fits
    // mag == 1 is 10, requires mag+1 slots
    //
    // also use fixed format if 0.01 or bigger
    if ((mag + 1) <= left_of_dot && mag >= -2) { // fixed format
      // how many spaces?
      int space_width = mag < 0 ? left_of_dot - 1 : left_of_dot - (mag + 1);

      std::stringstream formatstream2;
      formatstream2 << Utils::vformat("%% #.%df", precision_);
      datastream << std::string(space_width, ' ')
                 << Utils::vformat(formatstream2.str().c_str(), dat);

    } else { // sci format
      // 4 digits of exponent + dot + precision
      int pre_dot_width = width_ - precision_ - 5;

      // one digit plus sign
      int space_width = pre_dot_width - 2;

      std::stringstream formatstream2;
      formatstream2 << Utils::vformat("%% #.%de", precision_);
      datastream << std::string(space_width, ' ')
                 << Utils::vformat(formatstream2.str().c_str(), dat);
    }
  }
  return datastream.str();
}


std::string
Formatter::formatHeader(std::string header, int c) const
{
  int header_width = header_width_ - cellnum_width_ - 4;
  if (header.size() > header_width) {
    header.erase(header_width);
  } else if (header.size() < header_width) {
    header.append(header_width - header.size(), ' ');
  }

  std::stringstream headerstream;
  headerstream << header << "(" << std::setw(cellnum_width_) << std::right << c << "): ";
  return headerstream.str();
}

} // namespace Utils
} // namespace Amanzi
