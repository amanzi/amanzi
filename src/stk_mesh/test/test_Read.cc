// -------------------------------------------------------------
/**
 * @file   test_Read.cc
 * @author William A. Perkins
 * @date Wed May 18 14:15:54 2011
 * 
 * @brief Some unit tests for reading a (serial) Exodus file and
 * building a STK_mesh::Mesh_maps_stk instance.
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 22, 2010 by William A. Perkins
// Last Change: Wed May 18 14:15:54 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <sstream>
#include <UnitTest++.h>
#include <Epetra_MpiComm.h>

#include "Exodus_readers.hh"
#include "Parallel_Exodus_file.hh"
#include "../Mesh_STK_factory.hh"
#include "Auditor.hh"

// -------------------------------------------------------------
// class ReadFixture
// -------------------------------------------------------------
class ReadFixture {
 protected: 

  static const bool verbose;
  virtual void my_read(const std::string& name) = 0;

 public:

  Epetra_MpiComm comm;
  const int nproc;

  Amanzi::AmanziMesh::STK::Mesh_STK_factory mf;
  Amanzi::AmanziMesh::Data::Fields nofields;
  Teuchos::RCP<Amanzi::AmanziMesh::Data::Data> meshdata;
  Amanzi::AmanziMesh::STK::Mesh_STK_Impl_p mesh;
  // Teuchos::RCP<Amanzi::AmanziMesh::STK::Mesh_maps_stk> maps;
  
  /// Default constructor.
  ReadFixture(void)
      : comm(MPI_COMM_WORLD),
        nproc(comm.NumProc()),
        mf(comm.Comm(), 1000)
  { }

  void read(const std::string& name)
  {
    int ierr(0);
    try { 
      my_read(name);
      // maps.reset(new Amanzi::AmanziMesh::STK::Mesh_maps_stk(mesh));
    } catch (const std::exception& e) {
      std::cerr << comm.MyPID() << ": error: " << e.what() << std::endl;
      ierr++;
    }
    int aerr(0);
    comm.SumAll(&ierr, &aerr, 1);

    if (aerr) {
      throw std::runtime_error("Problem in ReadFixture");
    }
  }
  
};

const bool ReadFixture::verbose(false);


// -------------------------------------------------------------
//  class SerialReadFixture
// -------------------------------------------------------------
class SerialReadFixture : public ReadFixture {
 protected:
  
  void my_read(const std::string& name)         
  {
    meshdata.reset(Amanzi::Exodus::read_exodus_file(name.c_str()));
    // meshdata->to_stream(std::cerr, verbose);
    mesh.reset(mf.build_mesh(*meshdata, nofields));
  }

 public:

  /// Default constructor.
  SerialReadFixture(void) 
      : ReadFixture()
  {}

};

// -------------------------------------------------------------
//  class ParallelReadFixture
// -------------------------------------------------------------
class ParallelReadFixture : public ReadFixture {
 protected:
  
  void my_read(const std::string& name) 
  {
    Amanzi::Exodus::Parallel_Exodus_file thefile(comm, name);
    meshdata = thefile.read_mesh();

    for (int p = 0; p < nproc; p++) {
      if (comm.MyPID() == p) {
        std::cerr << std::endl;
        std::cerr << ">>>>>> Process " << p << " Begin <<<<<<" << std::endl;
        meshdata->to_stream(std::cerr, verbose);
        std::cerr << ">>>>>> Process " << p << " End <<<<<<" << std::endl;
        std::cerr << std::endl;
      }
      comm.Barrier();
    }

    Teuchos::RCP<Epetra_Map> cmap(thefile.cellmap());
    Teuchos::RCP<Epetra_Map> vmap(thefile.vertexmap());
    if (verbose) {
      cmap->Print(std::cerr);
      vmap->Print(std::cerr);
    }
    mesh.reset(mf.build_mesh(*meshdata, *cmap, *vmap, nofields));
  }

 public:

  /// Default constructor.
  ParallelReadFixture(void) 
      : ReadFixture()
  {}

};

  

SUITE (Exodus)
{
  TEST (Processes) 
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    const int nproc(comm.NumProc());
    CHECK(nproc >= 1 && nproc <= 4);
  }

  TEST_FIXTURE (SerialReadFixture, SerialReader1)
  {
    if (nproc == 1) {
      read("../exodus/test_files/hex_11x11x11_ss.exo");
      CHECK_EQUAL(mesh->count_entities(stk::mesh::Node, Amanzi::AmanziMesh::OWNED), 11*11*11);
      CHECK_EQUAL(mesh->count_entities(stk::mesh::Face, Amanzi::AmanziMesh::OWNED), 10*10*11*3);
      CHECK_EQUAL(mesh->count_entities(stk::mesh::Element, Amanzi::AmanziMesh::OWNED), 10*10*10);
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Element), 3);
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Node), 20);
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Face), 20);

      stk::mesh::Part *p;
      int count;

      p = mesh->get_set(1, stk::mesh::Face);
      CHECK(p != NULL);
      count = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      CHECK_EQUAL(count, 600);
            
      p = mesh->get_set("element block 20000", stk::mesh::Element);
      CHECK(p != NULL);
      count = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      CHECK_EQUAL(count, 9);
            
      p = mesh->get_set("node set 103", stk::mesh::Node);
      CHECK(p != NULL);
      count = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      CHECK_EQUAL(count, 121);
            
      Auditor audit("hex_11x11x11_ss_", mesh);
      audit();
    }
  }

  TEST_FIXTURE (SerialReadFixture, PrismReader)
  {
    if (nproc == 1) {
      read("../exodus/test_files/prism.exo");
      CHECK_EQUAL(mesh->count_entities(stk::mesh::Node, Amanzi::AmanziMesh::OWNED), 1920);
      CHECK_EQUAL(mesh->count_entities(stk::mesh::Element, Amanzi::AmanziMesh::OWNED), 2634);
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Element), 1);

      Auditor audit("prism_", mesh);
      audit();
    }

  }
        
  TEST_FIXTURE (SerialReadFixture, MixedCoarseReader)
  {
    if (nproc == 1) {
      read("../exodus/test_files/mixed-coarse.exo");
      CHECK_EQUAL(mesh->count_entities(stk::mesh::Node, Amanzi::AmanziMesh::OWNED), 361);
      CHECK_EQUAL(mesh->count_entities(stk::mesh::Element, Amanzi::AmanziMesh::OWNED), 592);
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Element), 5);

      Auditor audit("mixed-coarse_", mesh);
      audit();
    }

  }


  TEST_FIXTURE (SerialReadFixture, MixedReader)
  {
    if (nproc == 1) {
      read("../exodus/test_files/mixed.exo");
      CHECK_EQUAL(mesh->count_entities(stk::mesh::Node, Amanzi::AmanziMesh::OWNED), 6495);
      CHECK_EQUAL(mesh->count_entities(stk::mesh::Element, Amanzi::AmanziMesh::OWNED), 23186);
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Element), 6);

      Auditor audit("mixed_", mesh);
      audit();
    }

  }



  TEST_FIXTURE (SerialReadFixture, SerialReader2)
  {
    if (nproc == 1) {
      read("../exodus/test_files/hex_4x4x4_ss.exo");
      CHECK_EQUAL(mesh->count_entities(stk::mesh::Node, Amanzi::AmanziMesh::OWNED), 4*4*4);
      CHECK_EQUAL(mesh->count_entities(stk::mesh::Face, Amanzi::AmanziMesh::OWNED), 3*3*4*3);
      CHECK_EQUAL(mesh->count_entities(stk::mesh::Element, Amanzi::AmanziMesh::OWNED), 3*3*3);
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Element), 3);
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Node), 21);
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Face), 21);

      stk::mesh::Part *p;
      int count;

      p = mesh->get_set("side set 1", stk::mesh::Face);
      CHECK(p != NULL);
      count = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      CHECK_EQUAL(count, 54);
            
      p = mesh->get_set("element block 20000", stk::mesh::Element);
      CHECK(p != NULL);
      count = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      CHECK_EQUAL(count, 9);
            
      p = mesh->get_set("node set 103", stk::mesh::Node);
      CHECK(p != NULL);
      count = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      CHECK_EQUAL(count, 16);
            
      Auditor audit("hex_4x4x4_ss_", mesh);
      audit();
    }
  }

  TEST_FIXTURE (ParallelReadFixture, ParallelReader1)
  {
    if (nproc > 1 && nproc <= 4) {
      read("../exodus/test_files/split1/hex_11x11x11_ss.par");
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Element), 3);
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Node), 20);
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Face), 20);
      int local, global;
      local = mesh->count_entities(stk::mesh::Node, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 11*11*11);
      local = mesh->count_entities(stk::mesh::Face, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 10*10*11*3);
      local = mesh->count_entities(stk::mesh::Element, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 10*10*10);

      stk::mesh::Part *p;

      p = mesh->get_set("side set 1", stk::mesh::Face);
      CHECK(p != NULL);
      local = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 600);
            
      p = mesh->get_set("element block 20000", stk::mesh::Element);
      CHECK(p != NULL);
      local = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 9);

      p = mesh->get_set("node set 103", stk::mesh::Node);
      CHECK(p != NULL);
      local = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 121);

      Auditor audit("hex_11x11x11_ss.par", mesh);
      audit();
    }
  }            

  TEST_FIXTURE (ParallelReadFixture, ParallelReader2)
  {
    if (nproc > 1 && nproc < 4) {
      read("../exodus/test_files/split1/hex_4x4x4_ss.par");
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Element), 3);
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Node), 21);
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Face), 21);
      int local, global;
      local = mesh->count_entities(stk::mesh::Node, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 4*4*4);
      local = mesh->count_entities(stk::mesh::Face, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 3*3*4*3);
      local = mesh->count_entities(stk::mesh::Element, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 3*3*3);

      stk::mesh::Part *p;

      p = mesh->get_set(1, stk::mesh::Face);
      CHECK(p != NULL);
      local = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 54);
            
      p = mesh->get_set(20000, stk::mesh::Element);
      CHECK(p != NULL);
      local = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 9);

      p = mesh->get_set(103, stk::mesh::Node);
      CHECK(p != NULL);
      local = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 16);

      Auditor audit("hex_4x4x4_ss.par", mesh);
      audit();
    }
  }            


  TEST_FIXTURE (ParallelReadFixture, PrismParallelReader)
  {
    if (nproc == 2 || nproc == 4) {
      read("../exodus/test_files/split1/prism.par");
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Element), 1);
      int local, global;
      local = mesh->count_entities(stk::mesh::Node, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 1920);
      local = mesh->count_entities(stk::mesh::Element, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 2634);

      stk::mesh::Part *p;

      p = mesh->get_set(1, stk::mesh::Element);
      CHECK(p != NULL);
      local = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 2634);

      Auditor audit("prism.par.", mesh);
      audit();
    }
  }            

  TEST_FIXTURE (ParallelReadFixture, MixedCoarseParallelReader)
  {
    if (nproc == 2 || nproc == 4) {
      read("../exodus/test_files/split1/mixed-coarse.par");
      CHECK_EQUAL(mesh->num_sets(stk::mesh::Element), 5);
      int local, global;
      local = mesh->count_entities(stk::mesh::Node, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 361);
      local = mesh->count_entities(stk::mesh::Element, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 592);

      stk::mesh::Part *p;

      p = mesh->get_set(1, stk::mesh::Element);
      CHECK(p != NULL);
      local = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 120);

      p = mesh->get_set(2, stk::mesh::Element);
      CHECK(p != NULL);
      local = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 48);

      p = mesh->get_set(3, stk::mesh::Element);
      CHECK(p != NULL);
      local = mesh->count_entities(*p, Amanzi::AmanziMesh::OWNED);
      comm.SumAll(&local, &global, 1);
      CHECK_EQUAL(global, 48);

      Auditor audit("mixed-coarse.par.", mesh);
      audit();
    }
  }            

  TEST (DirectReader) 
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // FIXME: need to be able to assign path from configuration

    std::string fpath("../exodus/test_files/");
    std::string fname("hex_11x11x11_ss");
    if (comm.NumProc() == 1) {
      fname += ".exo";
    } else {
      fpath += "split1/";
      fname += ".par";
    }
    fname = fpath + fname;

    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> 
      mesh(new Amanzi::AmanziMesh::STK::Mesh_STK(comm, fname.c_str()));

    Auditor audit("stk_mesh_read_", mesh);
    audit();
  }
      
} 

