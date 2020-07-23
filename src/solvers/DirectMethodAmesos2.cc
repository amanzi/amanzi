/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (amklinv@sandia.gov)
           Ethan Coon (coonet@ornl.gov)
*/
//!  Direct solvers via Trilinos.
/*!

.. warning:: undocumented

*/

#include "Epetra_Vector.h"
#include "Amesos2_Factory.hpp"

#include "Key.hh"
#include "VerboseObject.hh"

#include "DirectMethodAmesos2.hh"

namespace Amanzi {
namespace AmanziSolvers {


/* ******************************************************************
* Initialization from a parameter list.
****************************************************************** */
void DirectMethodAmesos2::set_parameters(Teuchos::ParameterList& plist)
{
  plist_ = plist;

  this->set_name(Keys::cleanPListName(plist.name()));
  solver_name_ = plist.get<std::string>("solver name", "Klu2");

  std::string vo_name = this->name()+" (Amesos " + solver_name_ + ")";
  vo_ = Teuchos::rcp(new VerboseObject(vo_name, plist));

  inited_ = true;
}




/* ******************************************************************
* Update sets symbolic structure
****************************************************************** */
void DirectMethodAmesos2::UpdateInverse() {
  AMANZI_ASSERT(inited_);
  AMANZI_ASSERT(h_.get());

  solver_ = Amesos2::create<Epetra_CrsMatrix, Epetra_MultiVector>(
      solver_name_, h_);
          
  if (!solver_.get()) {
    Errors::Message msg;
    msg << "DirectMethodAmesos2: solver \"" << solver_name_ << "\" is not available";
    Exceptions::amanzi_throw(msg);
  }
  auto plist = Teuchos::rcpFromRef(plist_);
  solver_->setParameters(plist);
  solver_->symbolicFactorization();
  updated_ = true;
}


/* ******************************************************************
* Compute sets symbolic structure
****************************************************************** */
void DirectMethodAmesos2::ComputeInverse() {
  AMANZI_ASSERT(updated_);
  solver_->numericFactorization();
  computed_ = true;
}


int DirectMethodAmesos2::ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const
{
  AMANZI_ASSERT(computed_);

  solver_->setB(Teuchos::rcpFromRef<const Epetra_MultiVector>(v));
  solver_->setX(Teuchos::rcpFromRef<Epetra_MultiVector>(hv));
  
  solver_->solve();
  returned_code_ = 1;
  return returned_code_;
}


}  // namespace AmanziSolvers
}  // namespace Amanzi
