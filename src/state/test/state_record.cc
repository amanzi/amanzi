/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

// Tests the Record and RecordSet objects, which store meta-data about data.

#include "UniqueHelpers.hh"
#include "UnitTest++.h"
#include "data/RecordSet.hh"
#include "errors.hh"

#include "Vec.hh"

using namespace Amanzi;

std::unique_ptr<RecordSet>
getRecordDouble()
{
  auto f = std::make_unique<RecordSet>("my_field");
  f->RequireRecord("", "my_owner");
  f->SetType<double>();
  f->CreateData();
  return f;
}

std::unique_ptr<RecordSet>
getRecordVec()
{
  auto f = std::make_unique<RecordSet>("my_field");
  f->RequireRecord("", "my_owner");
  auto& fac = f->SetType<Vec, VecFactory>();
  fac.set_size(2);
  f->CreateData();
  return f;
}

TEST(FIELD_MAKE_EMPTY)
{
  RecordSet f{ "my_field" };
  CHECK(!f.HasType());
}

TEST(FIELD_MAKE_NULLFACTORY)
{
  auto f = getRecordDouble();
}

TEST(FIELD_SET_NULLFACTORY)
{
  auto f = getRecordDouble();
  f->Set("my_owner", 1.1);
  CHECK_EQUAL(1.1, f->Get<double>());
}

TEST(FIELD_MAKE_FACTORY)
{
  auto f = getRecordVec();
}

TEST(FIELD_SET_FACTORY)
{
  auto f = getRecordVec();
  f->GetW<Vec>("my_owner").v[0] = 1.1;
  CHECK_EQUAL(1.1, f->Get<Vec>().v[0]);
}

TEST(FIELD_SET_TWICE_OK)
{
  auto f1 = getRecordVec();
  f1->SetType<Vec, VecFactory>();
  // f1->SetType<double,VecFactory>(); // NOTE: This should not compile!
  // f1->SetType<Vec,double>(); // NOTE: This should not compile!
  CHECK_THROW(f1->SetType<Vec>(), Errors::Message);

  auto f2 = getRecordDouble();
  f2->SetType<double>();
  // f2->SetType<double,VecFactory>();  // NOTE: This should not compile!
  // f2->SetType<Vec,double>();  // NOTE: This should not compile!
  CHECK_THROW(f2->SetType<Vec>(), Errors::Message);
}

TEST(FIELD_THROWS_ON_BAD_OWNER)
{
  auto f = getRecordDouble();
  CHECK_THROW(f->Set("not_my_owner", 1.1), Errors::Message);
}

TEST(FIELD_THROWS_ON_BAD_TYPE)
{
  auto f = getRecordDouble();
  f->Set("my_owner", 1.1);
  CHECK_THROW(f->Get<bool>(), Errors::Message);
}

TEST(FIELD_COPY)
{
  auto f = getRecordDouble();
  f->RequireRecord("prev", "my_other_owner");
  f->CreateData();

  f->Set("my_owner", 1.1);
  CHECK_EQUAL(1.1, f->Get<double>());
  f->Set("prev", "my_other_owner", 2.1);
  CHECK_EQUAL(1.1, f->Get<double>());
  CHECK_EQUAL(2.1, f->Get<double>("prev"));
}

// TEST(FIELD_SET_FROM_OTHER) {
//   auto f1 = getRecordDouble()->GetRecord("");
//   auto f2 = getRecordDouble()->GetRecord("");

//   f2.Set("my_owner", 0.0);
//   f1.Set("my_owner", 1.1);

//   f2.SetFromOther("my_owner", f1);
//   f1.Set("my_owner", 2.1);

//   CHECK_EQUAL(1.1, f2.Get<double>());
//   CHECK_EQUAL(2.1, f1.Get<double>());
// }
