/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

// Tests FileHDF5 writer

#include "UnitTest++.h"

#include <fstream>
#include <iterator>
#include <string>
#include <algorithm>

#include "Teuchos_Array.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "AmanziComm.hh"
#include "AmanziVector.hh"
#include "FileHDF5.hh"
#include "HDF5Reader.hh"

#include "Mesh_MSTK.hh"


SUITE(FILE_HDF5)
{
  using namespace Amanzi;

  TEST(ASCEMIO)
  {
    auto comm = Amanzi::getDefaultComm();
    auto mpi_comm = Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm);
    AMANZI_ASSERT(mpi_comm.get());

    iogroup_conf_t io_config;
    io_config.numIOgroups = 1;
    io_config.commIncoming = *mpi_comm->getRawMpiComm();

    iogroup_t io_group;
    int ierr = parallelIO_IOgroup_init(&io_config, &io_group);
    CHECK(!ierr);

    auto data_file = parallelIO_open_file("ascemio_test.h5", &io_group, FILE_CREATE);
    CHECK(data_file >= 0);

    std::string group_name = "/my_group/";
    ierr = parallelIO_create_dataset_group(group_name.c_str(), data_file, &io_group);
    CHECK(!ierr);

    ierr = parallelIO_close_dataset_group(data_file, &io_group);
    CHECK(!ierr);

    ierr = parallelIO_close_file(data_file, &io_group);
    CHECK(!ierr);

    ierr = parallelIO_IOgroup_cleanup(&io_group);
    CHECK(!ierr);
  }
  
  TEST(ATTR_WRITEREAD)
  {
    // write
    double sixone(6.1);
    int seven(7);
    bool truthy(true);
    std::string name("hello");

    {
      Amanzi::FileHDF5 out(Amanzi::getDefaultComm(), "test.h5", FILE_CREATE);
      out.WriteAttribute("six_one", "/", sixone);
      out.WriteAttribute("seven", "/", seven);
      out.WriteAttribute("truthy", "/", truthy);
      out.WriteAttribute("name", "/", name);
    }
    {
      Amanzi::FileHDF5 in(Amanzi::getDefaultComm(), "test.h5", FILE_READONLY);
      CHECK_CLOSE(
        sixone, in.template ReadAttribute<double>("six_one", "/"), 1.e-10);
      CHECK_EQUAL(seven, in.template ReadAttribute<int>("seven", "/"));
      CHECK_EQUAL(truthy, in.template ReadAttribute<bool>("truthy", "/"));
      CHECK_EQUAL(name, in.template ReadAttribute<std::string>("name", "/"));
    }
  }

  TEST(ATTR_WRITEREAD_IN_GROUP)
  {
    // write
    double sixone(6.1);
    int seven(7);
    bool truthy(true);
    std::string name("hello");

    {
      Amanzi::FileHDF5 out(Amanzi::getDefaultComm(), "test0.h5", FILE_CREATE);
      out.CreateGroup("/group/");

      out.WriteAttribute("six_one", "/group/", sixone);
      out.WriteAttribute("seven", "/group/", seven);
      out.WriteAttribute("truthy", "/group/", truthy);
      out.WriteAttribute("name", "/group/", name);
    }


    {
      Amanzi::FileHDF5 in(Amanzi::getDefaultComm(), "test0.h5", FILE_READONLY);
      CHECK_CLOSE(sixone,
                  in.template ReadAttribute<double>("six_one", "/group/"),
                  1.e-10);
      CHECK_EQUAL(seven, in.template ReadAttribute<int>("seven", "/group/"));
      CHECK_EQUAL(truthy, in.template ReadAttribute<bool>("truthy", "/group/"));
      CHECK_EQUAL(name,
                  in.template ReadAttribute<std::string>("name", "/group/"));
    }
  }


  TEST(VECTOR_WRITE_READ)
  {
    auto comm = getDefaultComm();
    int size = comm->getSize();
    int rank = comm->getRank();
    Map_ptr_type map = Teuchos::rcp(new Map_type(size * 3, 3, 0, comm));

    Vector_ptr_type vec = Teuchos::rcp(new Vector_type(map));
    {
      auto vecv = vec->getLocalViewHost();
      for (int i = 0; i != vecv.extent(0); ++i) { vecv(i, 0) = i + rank * 3; }
    }

    { // open and write
      Amanzi::FileHDF5 out(Amanzi::getDefaultComm(), "test1.h5", FILE_CREATE);
      out.CreateGroup("/group/");

      out.WriteVector("/group/vec", *vec);
    } // close

    // First, read the vec in serial using a different, plain HDF5 reader.
    if (rank == 0) {
      Teuchos::Array<double> vec_g1;

      { // open/RAII
        HDF5Reader reader("test1.h5");
        reader.ReadData("/group/vec", vec_g1);
      } // closes file

      for (int i = 0; i != size * 3; ++i) { CHECK_CLOSE(i, vec_g1[i], 1.e-10); }
    }

    // Next, can we read the vec and compare to what we wrote with?
    Vector_type vec_g2(map);

    { // open and read
      Amanzi::FileHDF5 in(Amanzi::getDefaultComm(), "test1.h5", FILE_READONLY);
      in.ReadVector("/group/vec", vec_g2);
    }

    vec_g2.update(-1, *vec, 1);
    CHECK_CLOSE(0.0, vec_g2.norm2(), 1.e-10);
  }

  TEST(MULTIVECTOR_WRITE_READ)
  {
    auto comm = getDefaultComm();
    int size = comm->getSize();
    int rank = comm->getRank();
    Map_ptr_type map = Teuchos::rcp(new Map_type(size * 3, 3, 0, comm));

    MultiVector_ptr_type multivec = Teuchos::rcp(new MultiVector_type(map, 2));
    {
      auto mvv = multivec->getLocalViewHost();
      for (int i = 0; i != mvv.extent(0); ++i) {
        mvv(i, 0) = (i + rank * 3) * 2;
        mvv(i, 1) = (i + rank * 3) * 2 + 1;
      }
    }

    { // open and write
      Amanzi::FileHDF5 out(Amanzi::getDefaultComm(), "test2.h5", FILE_CREATE);
      out.CreateGroup("/group/");

      out.WriteMultiVector({ "/group/multivec.0", "/group/multivec.1" },
                           *multivec);
    } // close

    // First, read the vec in serial using a different, plain HDF5 reader.
    if (rank == 0) {
      Teuchos::Array<double> multivec_0_g1;
      Teuchos::Array<double> multivec_1_g1;

      { // open/RAII
        HDF5Reader reader("test2.h5");
        reader.ReadData("/group/multivec.0", multivec_0_g1);
        reader.ReadData("/group/multivec.1", multivec_1_g1);
      } // closes file

      for (int i = 0; i != size * 3; ++i) {
        CHECK_CLOSE(i * 2, multivec_0_g1[i], 1.e-10);
        CHECK_CLOSE(i * 2 + 1, multivec_1_g1[i], 1.e-10);
      }
    }

    // Next, can we read the vec and compare to what we wrote with?
    MultiVector_type multivec_g2(map, 2);

    { // open and read
      Amanzi::FileHDF5 in(Amanzi::getDefaultComm(), "test2.h5", FILE_READONLY);
      in.ReadMultiVector({ "/group/multivec.0", "/group/multivec.1" },
                         multivec_g2);
    }

    multivec_g2.update(-1, *multivec, 1);
    Teuchos::Array<double> norms(2, 0.0);
    multivec_g2.norm2(norms);
    CHECK_CLOSE(0.0, norms[0], 1.e-10);
    CHECK_CLOSE(0.0, norms[1], 1.e-10);
  }

  // THIS FAILS due to some issues with ASCEMIO
#if 0

TEST(MULTIVECTOR_BLOCK_WRITE_READ) {
  auto comm = getDefaultComm();
  int size = comm->getSize();
  int rank = comm->getRank();
  Map_ptr_type map = Teuchos::rcp(new Map_type(size * 3, 3, 0, comm));
  
  MultiVector_ptr_type multivec = Teuchos::rcp(new MultiVector_type(map, 2));
  {
    auto mvv = multivec->getLocalViewHost();
    CHECK_EQUAL(3, mvv.extent(0));
    CHECK_EQUAL(2, mvv.extent(1));

    typedef decltype(mvv) view_type;
    static_assert(std::is_same<view_type::array_layout, Kokkos::LayoutLeft>::value, "layout not LayoutLeft");

    for (int i=0; i!=mvv.extent(0); ++i) {
      mvv(i,0) = rank*3*2 + i*2;
      mvv(i,1) = rank*3*2 + i*2 + 1;
    }
  }

  { // test our test to make sure we understand layout
    auto mvv2 = multivec->getLocalViewHost();
    for (int i=0; i!=3; ++i) {
      CHECK_CLOSE((i+rank*3)*2, mvv2(i,0), 1.e-10);
      CHECK_CLOSE((i+rank*3)*2+1, mvv2(i,1), 1.e-10);
    }

    // and in 1d to make sure we understand layout
    // NOTE: this would have to change for LayoutRight
    auto mvv1 = multivec->get1dView();
    for (int i=0; i!=3*2; ++i) {
      std::cout << " " << i << "," << (i%2)*3 + i/2 << "," << mvv1[i] << std::endl;
      CHECK_CLOSE(rank*3*2 + (i%3)*2 + i/3, mvv1[i], 1.e-10);
      if (i > 0) {
        CHECK_EQUAL(&mvv1[i], &mvv1[i-1]+1);
      }
    }
  }
  
  { // open and write
    Amanzi::FileHDF5 out(Amanzi::getDefaultComm(), "test3.h5", FILE_CREATE);
    out.CreateGroup("/group/");
    out.WriteMultiVectorBlock("/group/multivec_block", *multivec);
  } // close

  // First, read the vec in serial using a different, plain HDF5 reader.
  if (rank == 0) {
    Teuchos::SerialDenseMatrix<std::size_t,double> multivec_b_g1;

    { // open/RAII
      HDF5Reader reader("test3.h5");
      reader.ReadMatData("/group/multivec_block", multivec_b_g1);
    } // closes file

    // check contiguous form of data.
    for (int i=0; i!=size*3*2; ++i) {
      int rank = i / (3*2);
      int j = i % (3*2);
      
      CHECK_CLOSE(rank*3*2 + (j%3)*2 + j/3, multivec_b_g1.values()[i], 1.e-10);
    }

    // check non-contiguous form
    for (int i=0; i!=size*3; ++i) {
      // NOTE these indices are flipped because SerialDenseMatrix is LayoutLeft
      CHECK_CLOSE(i*2, multivec_b_g1[0][i], 1.e-10);
      CHECK_CLOSE(i*2+1, multivec_b_g1[1][i], 1.e-10);
    }
  }

  // Next, can we read the vec and compare to what we wrote with?
  MultiVector_type multivec_b_g2(map, 2);
  
  { // open and read
    Amanzi::FileHDF5 in(Amanzi::getDefaultComm(), "test3.h5", FILE_READONLY);
    in.ReadMultiVectorBlock("/group/multivec_block", multivec_b_g2);
  }

  multivec_b_g2.update(-1, *multivec, 1);
  Teuchos::Array<double> norms_b(2, 0.0);
  multivec_b_g2.norm2(norms_b);
  CHECK_CLOSE(0.0, norms_b[0], 1.e-10);
  CHECK_CLOSE(0.0, norms_b[1], 1.e-10);
}

#endif

  TEST(VECTOR_WRITE_READ_IN_GROUP_WITH_NAMES)
  {
    auto comm = getDefaultComm();
    int size = comm->getSize();
    int rank = comm->getRank();
    Map_ptr_type map = Teuchos::rcp(new Map_type(size * 3, 3, 0, comm));

    std::vector<std::string> names = { "/group/map.local",
                                       "/group/map.global" };
    auto multivec = Teuchos::rcp(new MultiVector_type(map, 2));
    {
      auto mvv = multivec->getLocalViewHost();
      for (int i = 0; i != mvv.extent(0); ++i) {
        mvv(i, 0) = i;
        mvv(i, 1) = map->getGlobalElement(i);
      }
    }

    { // open and write
      Amanzi::FileHDF5 out(Amanzi::getDefaultComm(), "test.h5", FILE_CREATE);
      out.CreateGroup("/group/");
      out.WriteMultiVector(names, *multivec);
    } // close


    Vector_type vec2(map);
    MultiVector_type multivec2(map, 2);

    { // open and read
      Amanzi::FileHDF5 in(Amanzi::getDefaultComm(), "test.h5", FILE_READONLY);
      in.ReadMultiVector(names, multivec2);
    }

    multivec2.update(-1, *multivec, 1);
    Teuchos::Array<double> norms(2, 0.0);
    multivec2.norm2(norms);
    CHECK_CLOSE(0.0, norms[0], 1.e-10);
    CHECK_CLOSE(0.0, norms[1], 1.e-10);
  }



}
