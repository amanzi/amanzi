/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
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
      out.writeAttribute("six_one", "/", sixone);
      out.writeAttribute("seven", "/", seven);
      out.writeAttribute("truthy", "/", truthy);
      out.writeAttribute("name", "/", name);
    }
    {
      Amanzi::FileHDF5 in(Amanzi::getDefaultComm(), "test.h5", FILE_READONLY);
      CHECK_CLOSE(sixone, in.template readAttribute<double>("six_one", "/"), 1.e-10);
      CHECK_EQUAL(seven, in.template readAttribute<int>("seven", "/"));
      CHECK_EQUAL(truthy, in.template readAttribute<bool>("truthy", "/"));
      CHECK_EQUAL(name, in.template readAttribute<std::string>("name", "/"));
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
      out.createGroup("/group/");

      out.writeAttribute("six_one", "/group/", sixone);
      out.writeAttribute("seven", "/group/", seven);
      out.writeAttribute("truthy", "/group/", truthy);
      out.writeAttribute("name", "/group/", name);
    }


    {
      Amanzi::FileHDF5 in(Amanzi::getDefaultComm(), "test0.h5", FILE_READONLY);
      CHECK_CLOSE(sixone, in.template readAttribute<double>("six_one", "/group/"), 1.e-10);
      CHECK_EQUAL(seven, in.template readAttribute<int>("seven", "/group/"));
      CHECK_EQUAL(truthy, in.template readAttribute<bool>("truthy", "/group/"));
      CHECK_EQUAL(name, in.template readAttribute<std::string>("name", "/group/"));
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
      auto vecv = vec->getLocalViewHost(Tpetra::Access::ReadWrite);
      for (int i = 0; i != vecv.extent(0); ++i) { vecv(i, 0) = i + rank * 3; }
    }

    { // open and write
      Amanzi::FileHDF5 out(Amanzi::getDefaultComm(), "test1.h5", FILE_CREATE);
      out.createGroup("/group/");

      out.writeVector("/group/vec", *vec);
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
      in.readVector("/group/vec", vec_g2);
    }

    vec_g2.update(-1, *vec, 1);
    CHECK_CLOSE(0.0, vec_g2.norm2(), 1.e-10);
  }

  TEST(MULTIVECTOR_WRITE_READ)
  {
    auto comm = getDefaultComm();
    int size = comm->getSize();
    int rank = comm->getRank();
    int global_size = 0;
    for (int i = 0; i != size; ++i) global_size += (i + 3);
    Map_ptr_type map = Teuchos::rcp(new Map_type(global_size, rank + 3, 0, comm));
    MultiVector_ptr_type multivec = Teuchos::rcp(new MultiVector_type(map, 2));
    {
      auto mvv = multivec->getLocalViewHost(Tpetra::Access::ReadWrite);
      for (int i = 0; i != mvv.extent(0); ++i) {
        mvv(i, 0) = 2 * map->getGlobalElement(i);
        mvv(i, 1) = 2 * map->getGlobalElement(i) + 1;
      }
      std::cout << "extent ( " << rank << ") = " << mvv.extent(0) << "," << mvv.extent(1)
                << std::endl;
      AMANZI_ASSERT(mvv.extent(0) == rank + 3);
      AMANZI_ASSERT(mvv.extent(1) == 2);
    }

    { // open and write
      Amanzi::FileHDF5 out(Amanzi::getDefaultComm(), "test2.h5", FILE_CREATE);
      out.createGroup("/group/");

      //out.writeMultiVector({ "/group/multivec.0", "/group/multivec.1" }, *multivec);
      out.writeMultiVector("/group/multivec", *multivec);
    } // close

    // First, read the vec in serial using a different, plain HDF5 reader.
    // NOTE: this probably SHOULD fail on exactly one of GPU or CPU -- needs a
    // check that this is LayoutLeft (right?) --etc
    // if constexpr(std::is_same<Kokkos::LayoutLeft, typename MultiVector_type::host_view_type::array_layout>::value) {
    if (rank == 0) {
      Teuchos::SerialDenseMatrix<std::size_t, double> multivec_g1;
      { // open/RAII
        HDF5Reader reader("test2.h5");
        reader.ReadMatData("/group/multivec", multivec_g1);
      } // closes file

      CHECK_EQUAL(global_size, multivec_g1.numCols());
      CHECK_EQUAL(2, multivec_g1.numRows());
      for (int i = 0; i != global_size; ++i) {
        CHECK_CLOSE(i * 2, multivec_g1(0, i), 1.e-10);
        CHECK_CLOSE(i * 2 + 1, multivec_g1(1, i), 1.e-10);
      }
    }
    // }
    // END ABOVE FAIL NOTE

    // Next, can we read the vec and compare to what we wrote with?
    MultiVector_type multivec_g2(map, 2);

    { // open and read
      Amanzi::FileHDF5 in(Amanzi::getDefaultComm(), "test2.h5", FILE_READONLY);
      in.readMultiVector("/group/multivec", multivec_g2);
    }

    {
      auto mvv = multivec_g2.getLocalViewHost(Tpetra::Access::ReadOnly);
      for (int i = 0; i != mvv.extent(0); ++i) {
        CHECK_CLOSE(2 * map->getGlobalElement(i), mvv(i, 0), 1.e-10);
        CHECK_CLOSE(2 * map->getGlobalElement(i) + 1, mvv(i, 1), 1.e-10);
      }
    }
    multivec_g2.update(-1, *multivec, 1);
    Teuchos::Array<double> norms(2, 0.0);
    multivec_g2.norm2(norms);
    CHECK_CLOSE(0.0, norms[0], 1.e-10);
    CHECK_CLOSE(0.0, norms[1], 1.e-10);
  }
}
