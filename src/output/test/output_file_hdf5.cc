/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

// Tests FileHDF5 writer

#include "UnitTest++.h"

#include "AmanziComm.hh"
#include "AmanziVector.hh"
#include "FileHDF5.hh"

#include "Mesh_MSTK.hh"

#include <fstream>
#include <iterator>
#include <string>
#include <algorithm>



SUITE(FILE_HDF5) {

  using namespace Amanzi;

TEST(ATTR_WRITEREAD) {
  // write
  double sixone(6.1);
  int seven(7);
  bool truthy(true);
  std::string name("hello");

  {
    Amanzi::FileHDF5 out(Amanzi::getDefaultComm(), "test.h5", FILE_CREATE);
    out.WriteAttribute(sixone, "six_one");
    out.WriteAttribute(seven, "seven");
    out.WriteAttribute(truthy, "truthy");
    out.WriteAttribute(name, "name");
  }


  {
    Amanzi::FileHDF5 in(Amanzi::getDefaultComm(), "test.h5", FILE_READONLY);
    CHECK_CLOSE(sixone, in.template ReadAttribute<double>("six_one"), 1.e-10);
    CHECK_EQUAL(seven, in.template ReadAttribute<int>("seven"));
    CHECK_EQUAL(truthy, in.template ReadAttribute<bool>("truthy"));
    CHECK_EQUAL(name, in.template ReadAttribute<std::string>("name"));
  }
}  


TEST(ATTR_WRITEREAD_IN_GROUP) {
  // write
  double sixone(6.1);
  int seven(7);
  bool truthy(true);
  std::string name("hello");

  {
    Amanzi::FileHDF5 out(Amanzi::getDefaultComm(), "test.h5", FILE_CREATE);
    out.CreateGroup("/group/");
    
    out.WriteAttribute(sixone, "six_one", "/group/");
    out.WriteAttribute(seven, "seven", "/group/");
    out.WriteAttribute(truthy, "truthy", "/group/");
    out.WriteAttribute(name, "name", "/group/");
  }


  {
    Amanzi::FileHDF5 in(Amanzi::getDefaultComm(), "test.h5", FILE_READONLY);
    CHECK_CLOSE(sixone, in.template ReadAttribute<double>("six_one", "/group/"), 1.e-10);
    CHECK_EQUAL(seven, in.template ReadAttribute<int>("seven", "/group/"));
    CHECK_EQUAL(truthy, in.template ReadAttribute<bool>("truthy", "/group/"));
    CHECK_EQUAL(name, in.template ReadAttribute<std::string>("name", "/group/"));
  }
}  


TEST(VECTOR_WRITE_READ_IN_GROUP) {
  auto comm = getDefaultComm();
  int size = comm->getSize();
  int rank = comm->getRank();
  Map_ptr_type map = Teuchos::rcp(new Map_type(size * 3, 3, 0, comm));
  Vector_ptr_type vec = Teuchos::rcp(new Vector_type(map));
  vec->putScalar(rank);

  MultiVector_ptr_type multivec = Teuchos::rcp(new MultiVector_type(map, 2));
  {
    auto mvv = multivec->getLocalViewHost();
    for (int i=0; i!=mvv.extent(0); ++i) {
      mvv(i,0) = i;
      mvv(i,1) = map->getGlobalElement(i);
    }
  }
  
  { // open and write
    Amanzi::FileHDF5 out(Amanzi::getDefaultComm(), "test.h5", FILE_CREATE);
    out.CreateGroup("/group/");
    
    out.WriteVector(*vec, "/group/rank");
    out.WriteMultiVector(*multivec, {"/group/map.0", "/group/map.1"});
  } // close


  Vector_type vec2(map);
  MultiVector_type multivec2(map, 2);
  
  { // open and read
    Amanzi::FileHDF5 in(Amanzi::getDefaultComm(), "test.h5", FILE_READONLY);
    in.ReadVector(vec2, "/group/rank");
    in.ReadMultiVector(multivec2, {"/group/map.0", "/group/map.1"});
  }

  vec2.update(-1, *vec, 1);
  CHECK_CLOSE(0.0, vec2.norm2(), 1.e-10);

  multivec2.update(-1, *multivec, 1);
  Teuchos::Array<double> norms(2, 0.0);
  multivec2.norm2(norms);
  CHECK_CLOSE(0.0, norms[0], 1.e-10);
  CHECK_CLOSE(0.0, norms[1], 1.e-10);

}


TEST(VECTOR_WRITE_READ_IN_GROUP_WITH_NAMES) {
  auto comm = getDefaultComm();
  int size = comm->getSize();
  int rank = comm->getRank();
  Map_ptr_type map = Teuchos::rcp(new Map_type(size * 3, 3, 0, comm));

  std::vector<std::string> names = {"/group/map.local", "/group/map.global"};
  MultiVector_ptr_type multivec = Teuchos::rcp(new MultiVector_type(map, 2));
  {
    auto mvv = multivec->getLocalViewHost();
    for (int i=0; i!=mvv.extent(0); ++i) {
      mvv(i,0) = i;
      mvv(i,1) = map->getGlobalElement(i);
    }
  }
  
  { // open and write
    Amanzi::FileHDF5 out(Amanzi::getDefaultComm(), "test.h5", FILE_CREATE);
    out.CreateGroup("/group/");
    out.WriteMultiVector(*multivec, names);
  } // close


  Vector_type vec2(map);
  MultiVector_type multivec2(map, 2);
  
  { // open and read
    Amanzi::FileHDF5 in(Amanzi::getDefaultComm(), "test.h5", FILE_READONLY);
    in.ReadMultiVector(multivec2, names);
  }

  multivec2.update(-1, *multivec, 1);
  Teuchos::Array<double> norms(2, 0.0);
  multivec2.norm2(norms);
  CHECK_CLOSE(0.0, norms[0], 1.e-10);
  CHECK_CLOSE(0.0, norms[1], 1.e-10);

}


}
