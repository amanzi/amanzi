/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

//  Tests the DataFactory object, which stores a factory and generates data as
//  part of the container functionality of State.


#include "UnitTest++.h"
#include "data/DataFactory.hh"

using namespace Amanzi;

#include "Vec.hh"

TEST(NULL_FACTORY)
{
  DataFactory fac = dataFactory<double, NullFactory>();
  auto s = fac.Create();
  s.Set(1.1);
  CHECK_EQUAL(1.1, s.Get<double>());
}

TEST(NULL_FACTORY_MISDIRECTED)
{
  DataFactory fac = dataFactory<double, NullFactory>();
  auto s = data<double>();
  fac.Create(s);
  s.Set(1.1);
  CHECK_EQUAL(1.1, s.Get<double>());
}

TEST(VEC_FACTORY)
{
  g_constructor_calls_default = 0;
  g_constructor_calls_main = 0;
  g_constructor_calls_copy = 0;
  DataFactory fac = dataFactory<Vec, VecFactory>();
  fac.GetW<Vec, VecFactory>().set_size(2);

  auto s = fac.Create();
  s.GetW<Vec>().v[0] = 1.1;
  CHECK_EQUAL(1.1, s.Get<Vec>().v[0]);
  CHECK_EQUAL(1, g_constructor_calls_main);
  CHECK_EQUAL(0, g_constructor_calls_copy);
  CHECK_EQUAL(0, g_constructor_calls_default);
}

TEST(VEC_FACTORY_MISDIRECTED)
{
  g_constructor_calls_default = 0;
  g_constructor_calls_main = 0;
  g_constructor_calls_copy = 0;
  DataFactory fac = dataFactory<Vec, VecFactory>();
  fac.GetW<Vec, VecFactory>().set_size(2);

  auto s = data<Vec>();
  fac.Create(s);
  s.GetW<Vec>().v[1] = 1.1;
  CHECK_EQUAL(1.1, s.Get<Vec>().v[1]);
  CHECK_EQUAL(1, g_constructor_calls_main);
  CHECK_EQUAL(0, g_constructor_calls_copy);
  CHECK_EQUAL(0, g_constructor_calls_default);
}
