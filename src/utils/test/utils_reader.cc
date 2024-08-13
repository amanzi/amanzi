/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "ReaderHDF5.hh"
#include "ReaderNetCDF.hh"
#include "UnitTest++.h"

using namespace Amanzi;

using Vec = Teuchos::Array<double>;
using IVec = Teuchos::Array<int>;
using Mat = Teuchos::SerialDenseMatrix<int,double>;

void test_vec1(const Vec& v) {
  CHECK_EQUAL(3, v.length());
  CHECK_CLOSE(0.0, v[0], 1.e-10);
  CHECK_CLOSE(1.0, v[1], 1.e-10);
  CHECK_CLOSE(2.0, v[2], 1.e-10);
}

void test_vec2(const Vec& v) {
  CHECK_EQUAL(3, v.length());
  CHECK_CLOSE(3.0, v[0], 1.e-10);
  CHECK_CLOSE(4.0, v[1], 1.e-10);
  CHECK_CLOSE(5.0, v[2], 1.e-10);
}


void test_ivec1(const IVec& v) {
  CHECK_EQUAL(3, v.length());
  CHECK_EQUAL(0, v[0]);
  CHECK_EQUAL(1, v[1]);
  CHECK_EQUAL(2, v[2]);
}

void test_mat1(const Mat& v) {
  CHECK_EQUAL(2, v.numRows());
  CHECK_EQUAL(2, v.numCols());
  CHECK_CLOSE(0.0, v(0,0), 1.e-10);
  CHECK_CLOSE(1.0, v(1,0), 1.e-10);
  CHECK_CLOSE(2.0, v(0,1), 1.e-10);
  CHECK_CLOSE(3.0, v(1,1), 1.e-10);
}

void test_mat2(const Mat& v) {
  CHECK_EQUAL(2, v.numRows());
  CHECK_EQUAL(2, v.numCols());
  CHECK_CLOSE(2.0, v(0,0), 1.e-10);
  CHECK_CLOSE(3.0, v(1,0), 1.e-10);
  CHECK_CLOSE(4.0, v(0,1), 1.e-10);
  CHECK_CLOSE(5.0, v(1,1), 1.e-10);
}

void test_mat3(const Mat& v) {
  CHECK_EQUAL(3, v.numRows());
  CHECK_EQUAL(2, v.numCols());
  CHECK_CLOSE(0.0, v(0,0), 1.e-10);
  CHECK_CLOSE(1.0, v(1,0), 1.e-10);
  CHECK_CLOSE(2.0, v(2,0), 1.e-10);
  CHECK_CLOSE(3.0, v(0,1), 1.e-10);
  CHECK_CLOSE(4.0, v(1,1), 1.e-10);
  CHECK_CLOSE(5.0, v(2,1), 1.e-10);
}


void testReader(const Amanzi::Reader& r) {
  // test valid things
  // -------------------
  CHECK(r.hasVariableOrGroup("vec1"));
  CHECK(r.hasVariableOrGroup("/vec1"));
  CHECK(r.hasVariableOrGroup("/group1/vec3"));

  // test vec1
  Vec v1;
  r.read("vec1", v1);
  test_vec1(v1);

  // test int vec1
  IVec iv1;
  r.read("int_vec1", iv1);
  test_ivec1(iv1);

  // test vec2
  Vec v2;
  r.read("vec2", v2, 0);
  test_vec1(v2);
  r.read("vec2", v2, 1);
  test_vec2(v2);

  // test mat1
  Mat m1;
  r.read("mat1", m1);
  test_mat1(m1);

  // test mat2
  Mat m2;
  r.read("mat2", m2, 0);
  test_mat1(m2);
  r.read("mat2", m2, 1);
  test_mat2(m2);

  // test vec3
  Vec v3;
  r.read("/group1/vec3", v3);
  test_vec1(v3);

  // test mat3
  Mat m3;
  r.read("/mat3", m3);
  test_mat3(m3);

  // test invalid things
  // -----------------------
  CHECK(!r.hasVariableOrGroup("vec17"));
  CHECK(!r.hasVariableOrGroup("/group1/thing3"));

  // can't read thing that doesn't exist
  CHECK_THROW(r.read("vec17", v1), Errors::Message);

  // can't read past end of UNLIMITED dimension
  CHECK_THROW(r.read("mat2", m2, 3), Errors::Message);

  // can't read a vector that is actually a matrix
  // CHECK_THROW(r.read("mat2", v2, 0), Errors::Message);

  // can't read a matrix that is actually a vector
  CHECK_THROW(r.read("vec1", m1), Errors::Message);
}


TEST(UTILS_HDF5_READER) {
  auto r = std::make_shared<ReaderHDF5>("test/test.h5");
  testReader(*r);
}

TEST(UTILS_NETCDF_READER) {
  auto r = std::make_shared<ReaderNetCDF>("test/test.nc");
  testReader(*r);
}





