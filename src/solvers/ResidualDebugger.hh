/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
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

#include "VerboseObject.hh"
#include "IOEvent.hh"
#include "TreeVector.hh"
#include "TreeVector_Utils.hh"
#include "HDF5_MPI.hh"

namespace Amanzi {
namespace AmanziSolvers {

class ResidualDebugger : public IOEvent {

 public:
  // Constructor
  ResidualDebugger(Teuchos::ParameterList& plist) :
      IOEvent(plist) {
    filebasename_ = plist_.get<std::string>("file name base", "amanzi_dbg");
  }

  template <class VectorSpace>
  void StartIteration(double time, int cycle, int attempt,
                      const VectorSpace& space) {}

  template <class Vector>
  void
  WriteVector(int iter,
              const Vector& res,
              const Teuchos::Ptr<const Vector>& u=Teuchos::null,
              const Teuchos::Ptr<const Vector>& du=Teuchos::null) {}

  
 protected:
  std::string filebasename_;
  bool on_;
  double time_;
  std::vector<Teuchos::RCP<HDF5_MPI> > vis_;
};


template <>
void ResidualDebugger::StartIteration<TreeVectorSpace>(double time, int cycle, int attempt,
                                                       const TreeVectorSpace& space);
template <>
void ResidualDebugger::WriteVector<TreeVector>(int iter,
              const TreeVector& res,
              const Teuchos::Ptr<const TreeVector>& u,
              const Teuchos::Ptr<const TreeVector>& du);

  

} // namespace AmanziSolver
} // namespace Amanzi


#endif
