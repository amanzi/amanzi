/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Solvers

  Debugging object for writing vectors to file within an iterative
  process for use with vis tools.
*/

#include "Mesh.hh"
#include "HDF5_MPI.hh"

#include "CompositeVector.hh"
#include "TreeVector.hh"
#include "TreeVector_Utils.hh"

#include "ResidualDebugger.hh"


namespace Amanzi {
namespace AmanziSolvers {

//
// TreeVector does work
// -----------------------------------------------------------------------------
template <>
void
ResidualDebugger::StartIteration<TreeVectorSpace>(double time,
                                                  int cycle,
                                                  int attempt,
                                                  const TreeVectorSpace& space)
{
  on_ = DumpRequested(cycle, time);
  time_ = time;
  if (on_) {
    // iterate through the TreeVector finding leaf nodes and write them
    std::vector<Teuchos::RCP<const TreeVectorSpace>> leaves = collectTreeVectorLeaves_const(space);
    vis_.resize(leaves.size());

    for (int i = 0; i != leaves.size(); ++i) {
      if (leaves[i]->Data()->HasComponent("cell")) {
        std::stringstream filename;
        filename << filebasename_ << cycle << "_a" << attempt << "_v" << i;
        vis_[i] = Teuchos::rcp(new HDF5_MPI(leaves[i]->Data()->Mesh()->get_comm()));
        vis_[i]->setTrackXdmf(true);
        vis_[i]->createMeshFile(leaves[i]->Data()->Mesh(), filename.str() + "_mesh");
        vis_[i]->createDataFile(filename.str());
      }
    }
  }
}


//
// Write a vector individually.
// -----------------------------------------------------------------------------
template <>
void
ResidualDebugger::WriteVector<TreeVector>(int iter,
                                          const TreeVector& res,
                                          const Teuchos::Ptr<const TreeVector>& u,
                                          const Teuchos::Ptr<const TreeVector>& du)
{
  if (on_) {
    // open files
    for (std::vector<Teuchos::RCP<HDF5_MPI>>::iterator it = vis_.begin(); it != vis_.end(); ++it) {
      if (it->get()) {
        (*it)->writeMesh(time_, iter);
        (*it)->createTimestep(time_, iter, "");
        (*it)->open_h5file();
      }
    }

    // write residuals
    std::vector<Teuchos::RCP<const TreeVector>> r_leaves = collectTreeVectorLeaves_const(res);
    for (int i = 0; i != r_leaves.size(); ++i) {
      if (vis_[i].get()) {
        const Epetra_MultiVector& vec = *r_leaves[i]->Data()->ViewComponent("cell", false);
        for (int j = 0; j != vec.NumVectors(); ++j) {
          std::stringstream my_name;
          my_name << "residual.cell." << j;
          vis_[i]->writeCellDataReal(*vec(j), my_name.str());
        }
      }
    }

    // write values
    if (u.get()) {
      std::vector<Teuchos::RCP<const TreeVector>> u_leaves = collectTreeVectorLeaves_const(*u);
      for (int i = 0; i != u_leaves.size(); ++i) {
        if (vis_[i].get()) {
          const Epetra_MultiVector& vec = *u_leaves[i]->Data()->ViewComponent("cell", false);
          for (int j = 0; j != vec.NumVectors(); ++j) {
            std::stringstream my_name;
            my_name << "u.cell." << j;
            vis_[i]->writeCellDataReal(*vec(j), my_name.str());
          }
        }
      }
    }

    // write corrections
    if (du.get()) {
      std::vector<Teuchos::RCP<const TreeVector>> du_leaves = collectTreeVectorLeaves_const(*du);
      for (int i = 0; i != du_leaves.size(); ++i) {
        if (vis_[i].get()) {
          const Epetra_MultiVector& vec = *du_leaves[i]->Data()->ViewComponent("cell", false);
          for (int j = 0; j != vec.NumVectors(); ++j) {
            std::stringstream my_name;
            my_name << "du.cell." << j;
            vis_[i]->writeCellDataReal(*vec(j), my_name.str());
          }
        }
      }
    }

    // close files
    for (std::vector<Teuchos::RCP<HDF5_MPI>>::iterator it = vis_.begin(); it != vis_.end(); ++it) {
      if (it->get()) {
        (*it)->close_h5file();
        (*it)->endTimestep();
      }
    }
  }
}


} // namespace AmanziSolvers
} // namespace Amanzi
