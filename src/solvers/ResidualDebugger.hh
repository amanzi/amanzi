/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

License: see $AMANZI_DIR/COPYRIGHT
Author: Ethan Coon

Debugging object for writing vectors to file within an iterative
process for use with vis tools.

------------------------------------------------------------------------- */

#ifndef AMANZI_RESIDUAL_DEBUGGER_HH_
#define AMANZI_RESIDUAL_DEBUGGER_HH_

#include <string>

#include "Teuchos_Ptr.hpp"
#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"
#include "IOEvent.hh"

namespace Amanzi {

class HDF5_MPI;

namespace AmanziSolvers {

template<class Vector, class VectorSpace>
class ResidualDebugger : public IOEvent {

 public:
  // Constructor
  ResidualDebugger(Teuchos::ParameterList& plist) :
      IOEvent(plist) {
    filebasename_ = plist_.get<std::string>("file name base","amanzi_dbg");
  }
  
  void StartIteration(double time, int cycle, int attempt,
                      const VectorSpace& space);

  void
  WriteVector(int iter,
              const Vector& res,
              const Teuchos::Ptr<const Vector>& u=Teuchos::null,
              const Teuchos::Ptr<const Vector>& du=Teuchos::null);

 protected:
  std::string filebasename_;
  bool on_;
  double time_;
  std::vector<Teuchos::RCP<HDF5_MPI> > vis_;
};


//
// Default doesn't work
// -----------------------------------------------------------------------------
template<class Vector, class VectorSpace>
void
ResidualDebugger<Vector,VectorSpace>::StartIteration(
                                 double time, int cycle, int attempt,
                                 const VectorSpace& space) {}


//  
// Write a vector individually.  Default does not work.
// -----------------------------------------------------------------------------
template<class Vector, class VectorSpace>
void
ResidualDebugger<Vector,VectorSpace>::WriteVector(int iter,
			      const Vector& res,
                              const Teuchos::Ptr<const Vector>& u,
                              const Teuchos::Ptr<const Vector>& du) {}
  

} // namespace AmanziSolver
} // namespace Amanzi

#endif
