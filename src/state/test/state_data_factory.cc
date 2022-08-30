#include "UnitTest++.h"
#include "data/DataFactory.hh"

using namespace Amanzi;

#include "Vec.hh"

TEST(NULL_FACTORY) {
  DataFactory fac = dataFactory<double, NullFactory>();
  CHECK(fac.ValidType<double>());
  CHECK(!fac.ValidType<int>());
  CHECK(!fac.ValidType<Vec>());

  auto s = fac.Create();
  CHECK(fac.ValidType<double>());
  CHECK(!fac.ValidType<int>());
  CHECK(!fac.ValidType<Vec>());

  s.Assign(1.1);
  CHECK_EQUAL(1.1, s.Get<double>());
}

TEST(NULL_FACTORY_MISDIRECTED) {
  DataFactory fac = dataFactory<double, NullFactory>();
  auto s = data<double>();
  fac.Create(s);
  s.Assign(1.1);
  CHECK_EQUAL(1.1, s.Get<double>());
}

TEST(VEC_FACTORY) {
  g_constructor_calls_default = 0;
  g_constructor_calls_main = 0;
  g_constructor_calls_copy = 0;
  DataFactory fac = dataFactory<Vec, VecFactory>();

  bool valid = fac.ValidType<Vec>(); CHECK(valid);
  valid = fac.ValidType<Vec,VecFactory>(); CHECK(valid);
  valid = !fac.ValidType<double>(); CHECK(valid);
  valid = !fac.ValidType<Vec,double>(); CHECK(valid);
  valid = !fac.ValidType<double,VecFactory>(); CHECK(valid);

  fac.GetW<Vec, VecFactory>().set_size(2);

  auto s = fac.Create();
  valid = fac.ValidType<Vec>(); CHECK(valid);
  valid = fac.ValidType<Vec,VecFactory>(); CHECK(valid);
  valid = !fac.ValidType<double>(); CHECK(valid);
  valid = !fac.ValidType<Vec,double>(); CHECK(valid);
  valid = !fac.ValidType<double,VecFactory>(); CHECK(valid);

  s.GetW<Vec>().v[0] = 1.1;
  CHECK_EQUAL(1.1, s.Get<Vec>().v[0]);
  CHECK_EQUAL(1, g_constructor_calls_main);
  CHECK_EQUAL(0, g_constructor_calls_copy);
  CHECK_EQUAL(0, g_constructor_calls_default);
}

TEST(VEC_FACTORY_MISDIRECTED) {
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
