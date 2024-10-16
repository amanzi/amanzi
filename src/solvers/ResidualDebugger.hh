/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Solvers

*/

/*!

Debugging object for writing vectors to file within an iterative
process for use with vis tools.

.. _residual-debugger-spec:
.. admonition:: residual-debugger-spec

    * `"file name base`" ``[string]`` **amanzi_dbg** Prefix for output filenames.

    INCLUDES:

    - ``[io-event-spec]`` An IOEvent_ spec

*/

#ifndef AMANZI_RESIDUAL_DEBUGGER_HH_
#define AMANZI_RESIDUAL_DEBUGGER_HH_

#include <string>

#include "Teuchos_Ptr.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Tag.hh"
#include "VerboseObject.hh"
#include "IOEvent.hh"
#include "TreeVector.hh"
#include "TreeVector_Utils.hh"
#include "HDF5_MPI.hh"
#include "StateDefs.hh"

namespace Amanzi {

class State;

namespace AmanziSolvers {

class ResidualDebugger : public Utils::IOEvent {
 public:
  // Constructor
  ResidualDebugger(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<State>& S = Teuchos::null,
                   const Tag& tag = Tags::DEFAULT)
    : Utils::IOEvent(plist), S_(S), tag_(tag)
  {
    filebasename_ = plist_.get<std::string>("file name base", "amanzi_dbg");
    additional_vars_ = plist_.get<Teuchos::Array<std::string>>("additional variables", {});
  }

  template <class VectorSpace>
  void StartIteration(int attempt, const VectorSpace& space)
  {}

  template <class Vector>
  void WriteVector(int iter,
                   const Vector& res,
                   const Teuchos::Ptr<const Vector>& u = Teuchos::null,
                   const Teuchos::Ptr<const Vector>& du = Teuchos::null)
  {}


 protected:
  std::string filebasename_;
  bool on_;
  std::vector<Teuchos::RCP<HDF5_MPI>> vis_;

  Teuchos::RCP<State> S_;
  Teuchos::Array<std::string> additional_vars_;
  Tag tag_;
};


template <>
void
ResidualDebugger::StartIteration<TreeVectorSpace>(int attempt, const TreeVectorSpace& space);
template <>
void
ResidualDebugger::WriteVector<TreeVector>(int iter,
                                          const TreeVector& res,
                                          const Teuchos::Ptr<const TreeVector>& u,
                                          const Teuchos::Ptr<const TreeVector>& du);


} // namespace AmanziSolvers
} // namespace Amanzi


#endif
