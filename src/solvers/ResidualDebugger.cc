/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

License: see $AMANZI_DIR/COPYRIGHT
Author: Ethan Coon

Debugging object for writing vectors to file within an iterative
process for use with vis tools.

------------------------------------------------------------------------- */


#include "Mesh.hh"
#include "hdf5mpi_mesh.hh"

#include "CompositeVector.hh"
#include "TreeVector.hh"
#include "TreeVector_Utils.hh"

#include "ResidualDebugger.hh"


namespace Amanzi {

// Constructor
ResidualDebugger::ResidualDebugger(Teuchos::ParameterList& plist) :
    IOEvent(plist) {

  filebasename_ = plist_.get<std::string>("file name base","amanzi_dbg");
}


void
ResidualDebugger::StartIteration(double time, int cycle, int attempt,
                                 const TreeVectorSpace& space) {
  on_ = DumpRequested(cycle, time);
  time_ = time;
  if (on_) {
    // iterate through the TreeVector finding leaf nodes and write them
    std::vector<Teuchos::RCP<const TreeVectorSpace> > leaves =
      collectTreeVectorLeaves(space);
    vis_.resize(leaves.size());
  
    for (int i=0; i!=leaves.size(); ++i) {
      if (leaves[i]->Data()->HasComponent("cell")) {
        std::string filename = filebasename_ + std::to_string(cycle)
          + "_a" + std::to_string(attempt) + "_v" + std::to_string(i) + ".h5";
        vis_[i] = Teuchos::rcp(new HDF5_MPI(leaves[i]->Data()->Mesh()->get_comm()));
        vis_[i]->createDataFile(filename);
      }
    }
  }
}

  
// Write a vector individually.
void
ResidualDebugger::WriteVector(int iter,
			      const TreeVector& res,
                              const Teuchos::Ptr<const TreeVector>& u,
                              const Teuchos::Ptr<const TreeVector>& du) {

  if (on_) {
    // open files
    for (std::vector<Teuchos::RCP<HDF5_MPI> >::iterator it=vis_.begin();
         it!=vis_.end(); ++it) {
      if (it->get()) {
        (*it)->createTimestep(time_, iter);
        (*it)->open_h5file();
      }
    }

    // write residuals
    std::vector<Teuchos::RCP<const TreeVector> > r_leaves =
      collectTreeVectorLeaves(res);
    for (int i=0; i!=r_leaves.size(); ++i) {
      if (vis_[i]->get()) {
        const Epetra_MultiVector& vec = *r_leaves[i]->Data()
          ->ViewComponent("cell",false);
        for (int j=0; j!=vec.NumVectors(); ++j) {
          std::string my_name = "residual.cell." + std::to_string(j);
          vis_[i]->writeCellDataReal(*vec(j), my_name);
        }
      }
    }

    // write values
    if (u.get()) {
      std::vector<Teuchos::RCP<const TreeVector> > u_leaves =
        collectTreeVectorLeaves(*u);
      for (int i=0; i!=u_leaves.size(); ++i) {
        if (vis_[i]->get()) {
          const Epetra_MultiVector& vec = *u_leaves[i]->Data()
            ->ViewComponent("cell",false);
          for (int j=0; j!=vec.NumVectors(); ++j) {
            std::string my_name = "u.cell." + std::to_string(j);
            vis_[i]->writeCellDataReal(*vec(j), my_name);
          }
        }
      }
    }

    // write corrections
    if (du.get()) {
      std::vector<Teuchos::RCP<const TreeVector> > du_leaves =
        collectTreeVectorLeaves(*du);
      for (int i=0; i!=du_leaves.size(); ++i) {
        if (vis_[i]->get()) {
          const Epetra_MultiVector& vec = *du_leaves[i]->Data()
            ->ViewComponent("cell",false);
          for (int j=0; j!=vec.NumVectors(); ++j) {
            std::string my_name = "du.cell." + std::to_string(j);
            vis_[i]->writeCellDataReal(*vec(j), my_name);
          }
        }
      }
    }

    // close files
    for (std::vector<Teuchos::RCP<HDF5_MPI> >::iterator it=vis_.begin();
         it!=vis_.end(); ++it) {
      if (it->get()) {
        (*it)->close_h5file();
        (*it)->endTimestep();
      }
    }
  }    
}

  
} // namespace Amanzi
