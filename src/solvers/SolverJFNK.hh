/*
  Authors: Ethan Coon (ecoon@lanl.gov)

  Interface for using JFNK as a solver with NKA for the linear solve.
*/

#ifndef AMANZI_JFNK_SOLVER_
#define AMANZI_JFNK_SOLVER_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"
#include "SolverNewton.hh"
#include "SolverFnBase.hh"
#include "SolverDefs.hh"

#include "SolverFnBaseJF.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Vector, class VectorSpace>
class SolverJFNK : public SolverNewton<Vector,VectorSpace> {

 public:
  SolverJFNK(Teuchos::ParameterList& plist) :
      SolverNewton<Vector,VectorSpace>(plist) {}

  SolverJFNK(Teuchos::ParameterList& plist,
	     const Teuchos::RCP<SolverFnBase<Vector> >& fn,
	     const VectorSpace& map) :
      SolverNewton<Vector,VectorSpace>(plist) {
    Init(fn, map);
  }

  void Init(const Teuchos::RCP<SolverFnBase<Vector> >& fn,
            const VectorSpace& map);

 protected:
  Teuchos::RCP<SolverFnBaseJF<Vector,VectorSpace> > jf_fnbase_;
  Teuchos::RCP<SolverFnBase<Vector> > wrapped_fnbase_;

};


/* ******************************************************************
* Public Init method.
****************************************************************** */
template<class Vector, class VectorSpace>
void
SolverJFNK<Vector,VectorSpace>::Init(const Teuchos::RCP<SolverFnBase<Vector> >& fn,
        const VectorSpace& map)
{
  wrapped_fnbase_ = fn;

  // wrap the fn in a JF fn
  jf_fnbase_ = Teuchos::rcp(new SolverFnBaseJF<Vector,VectorSpace>(this->plist_,wrapped_fnbase_, map));

  // Init the Newton method with this wrapped SolverFnBase
  SolverNewton<Vector,VectorSpace>::Init(jf_fnbase_, map);
  
  // Change the verbose object header
  this->vo_->set_name("Solver::JFNK");
}

} // namespace
} // namespace

#endif
