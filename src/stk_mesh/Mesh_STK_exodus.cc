// -------------------------------------------------------------
// file: Mesh_STK_exodus.cc
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created May  9, 2011 by William A. Perkins
// Last Change: Tue May 17 11:45:27 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------


static const char* SCCS_ID = "$Id$ Battelle PNL";

#include "Mesh_STK.hh"
#include "Exodus_readers.hh"
#include "Parallel_Exodus_file.hh"
#include "stk_mesh_error.hh"
#include "Mesh_STK_factory.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace STK {

// -------------------------------------------------------------
//  class Mesh_STK
// -------------------------------------------------------------

// -------------------------------------------------------------
// Mesh_STK::read_exodus_
// -------------------------------------------------------------
void
Mesh_STK::read_exodus_(const std::string& fname)
{
  const int nproc(communicator_->NumProc());
  const int me(communicator_->MyPID());
  int ierr(0);
    
  STK::Mesh_STK_factory mf(communicator_->GetMpiComm(), 1000);
  Data::Fields nofields;
  Teuchos::RCP<Data::Data> meshdata;
    
  try {
    if (nproc == 1) {
      meshdata.reset(Exodus::read_exodus_file(fname.c_str()));
      mesh_.reset(mf.build_mesh(*meshdata, nofields));
    } else {
      Exodus::Parallel_Exodus_file thefile(*communicator_, fname);
      meshdata = thefile.read_mesh();
      mesh_.reset(mf.build_mesh(*meshdata, 
                                *(thefile.cellmap()), 
                                *(thefile.vertexmap()), 
                                nofields));
    }
  } catch (const std::exception& e) {
    std::cerr << communicator_->MyPID() << ": error: " << e.what() << std::endl;
    ierr++;
  }
  int aerr(0);
  communicator_->SumAll(&ierr, &aerr, 1);
  if (aerr != 0) 
    Exceptions::amanzi_throw( STK::Error ("Exodus file read error") );
  build_maps_();
}

// -------------------------------------------------------------
// Mesh_STK::Mesh_STK
// -------------------------------------------------------------
Mesh_STK::Mesh_STK(const Epetra_MpiComm& comm, 
                   const std::string& fname)
    : communicator_(new Epetra_MpiComm(comm)),
      mesh_(), 
      entity_map_(3),            // FIXME: needs to come from the file
      map_owned_(), map_used_()

      
{
  Mesh::set_comm(communicator_->GetMpiComm());
  read_exodus_(fname);
}

Mesh_STK::Mesh_STK(const char *fname, MPI_Comm comm)
    : communicator_(new Epetra_MpiComm(comm)),
      mesh_(), 
      entity_map_(3),           // FIXME: needs to come from the file
      map_owned_(), map_used_()

      
{
  Mesh::set_comm(communicator_->GetMpiComm());
  read_exodus_(fname);
}

} // namespace STK
} // namespace AmanziMesh
} // namespace Amanzi
