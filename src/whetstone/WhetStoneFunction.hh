/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Virtual class for continuous functions. It has limited interface
  compared to class Amanzi::Function
*/

#ifndef AMANZI_WHETSTONE_FUNCTION_HH_
#define AMANZI_WHETSTONE_FUNCTION_HH_

#include "Point.hh"

namespace Amanzi {
namespace WhetStone {

class WhetStoneFunction {
 public:
  WhetStoneFunction() {};
  virtual ~WhetStoneFunction() {};

  virtual double Value(const AmanziGeometry::Point& xp) const = 0;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
