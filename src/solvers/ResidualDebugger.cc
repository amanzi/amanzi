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
#include "CompositeVector.hh"
#include "TreeVector.hh"
#include "DataStructuresHelpers.hh"
#include "OutputFactory.hh"
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
      if (leaves[i]->getData()->hasComponent("cell")) {
        std::stringstream filename;
        filename << filebasename_ << cycle << "_a" << attempt << "_v" << i;
        Teuchos::ParameterList plist;
        plist.set<std::string>("file name base", filename.str());
        vis_[i] = OutputFactory::createForVis(plist,
                AmanziMesh::onMemSpace<MemSpace_kind::HOST>(leaves[i]->getData()->getMesh()));
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
    for (const auto& vis : vis_) {
      if (vis.get()) {
        vis->createTimestep(time_, iter);
      }
    }

    // write residuals
    std::vector<Teuchos::RCP<const TreeVector>> r_leaves = collectTreeVectorLeaves_const(res);
    for (int i = 0; i != r_leaves.size(); ++i) {
      if (vis_[i].get()) {
        Teuchos::ParameterList attrs("residual");
        vis_[i]->write(attrs, *r_leaves[i]->getData());
      }
    }

    // write values
    if (u.get()) {
      std::vector<Teuchos::RCP<const TreeVector>> u_leaves = collectTreeVectorLeaves_const(*u);
      for (int i = 0; i != u_leaves.size(); ++i) {
        if (vis_[i].get()) {
          Teuchos::ParameterList attrs("u");
          vis_[i]->write(attrs, *u_leaves[i]->getData());
        }
      }
    }

    // write corrections
    if (du.get()) {
      std::vector<Teuchos::RCP<const TreeVector>> du_leaves = collectTreeVectorLeaves_const(*du);
      for (int i = 0; i != du_leaves.size(); ++i) {
        if (vis_[i].get()) {
          Teuchos::ParameterList attrs("du");
          vis_[i]->write(attrs, *du_leaves[i]->getData());
        }
      }
    }

    // close files
    for (const auto& vis : vis_) {
      if (vis.get()) {
        vis->finalizeTimestep();
      }
    }
  }
}


} // namespace AmanziSolvers
} // namespace Amanzi
