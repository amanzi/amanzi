/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "UniqueHelpers.hh"
#include "UnitTest++.h"
#include <boost/any.hpp>
#include <vector>

TEST(ANY)
{
  std::vector<boost::any> stuff;

  boost::any thing1 = std::make_unique<double>(1.1);
  stuff.push_back(thing1);

  boost::any thing2 = std::make_unique<bool>(true);
  stuff.push_back(thing2);
}
