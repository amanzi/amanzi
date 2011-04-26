// -------------------------------------------------------------
/**
 * @file   Mesh_maps_stk_exodus.cc
 * @author William A. Perkins
 * @date Thu Apr 21 12:58:49 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December 13, 2010 by William A. Perkins
// Last Change: Thu Apr 21 12:58:49 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include "Mesh_maps_stk.hh"
#include "Mesh_factory.hh"
#include "Exodus_readers.hh"
#include "Parallel_Exodus_file.hh"
#include "stk_mesh_error.hh"

namespace STK_mesh {

  // -------------------------------------------------------------
  // Mesh_maps_stk::read_exodus_
  // -------------------------------------------------------------
  void
  Mesh_maps_stk::read_exodus_(const Epetra_MpiComm& comm, 
                              const std::string& fname)
  {
    const int nproc(communicator_->NumProc());
    const int me(communicator_->MyPID());
    int ierr(0);
    
    STK_mesh::Mesh_factory mf(comm.Comm(), 1000);
    Mesh_data::Fields nofields;
    Teuchos::RCP<Mesh_data::Data> meshdata;
    
    try {
      if (nproc == 1) {
        meshdata.reset(ExodusII::read_exodus_file(fname.c_str()));
        mesh_.reset(mf.build_mesh(*meshdata, nofields));
      } else {
        ExodusII::Parallel_Exodus_file thefile(comm, fname);
        meshdata = thefile.read_mesh();
        mesh_.reset(mf.build_mesh(*meshdata, 
                                 *(thefile.cellmap()), 
                                 *(thefile.vertexmap()), 
                                 nofields));
      }
    } catch (const std::exception& e) {
      std::cerr << comm.MyPID() << ": error: " << e.what() << std::endl;
      ierr++;
    }
    int aerr(0);
    comm.SumAll(&ierr, &aerr, 1);
    if (aerr != 0) 
      Exceptions::amanzi_throw( STKMeshError ("Exodus file read error") );
    build_maps_();
  }


  // -------------------------------------------------------------
  // Mesh_maps_stk::Mesh_maps_stk
  // -------------------------------------------------------------
  Mesh_maps_stk::Mesh_maps_stk(const Epetra_MpiComm& comm, 
                               const std::string& fname)
    : mesh_(), 
      entity_map_ (3),          // FIXME: needs to come from the file
      communicator_(comm.Clone())
  {
    read_exodus_(comm, fname);
  }

  Mesh_maps_stk::Mesh_maps_stk(const char *fname, MPI_Comm comm)
    : mesh_(), 
      entity_map_ (3),          // FIXME: needs to come from the file
      communicator_(new Epetra_MpiComm(comm))
  {
    read_exodus_(comm, fname);
  }

  
} // close namespace STK_mesh
