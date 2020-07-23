/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Ifpack.h"

#include "exceptions.hh"
#include "VerboseObject.hh"
#include "PreconditionerIfpack.hh"

namespace Amanzi {
namespace AmanziSolvers {

/* ******************************************************************
* Apply the preconditioner.
* According to IfPack documentation, the error code is set to 0 if 
* the inversion was successful. 
****************************************************************** */
int PreconditionerIfpack::ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const
{
  returned_code_ = IfpILU_->ApplyInverse(v, hv);
  return returned_code_;
}


/* ******************************************************************
 * Initialize the preconditioner.
 ****************************************************************** */
void PreconditionerIfpack::set_parameters(Teuchos::ParameterList& plist)
{
  plist_ = plist;
  std::string vo_name = this->name()+" ("+plist_.get<std::string>("method")+")";
  vo_ = Teuchos::rcp(new VerboseObject(vo_name, plist_));
  initialized_ = true;
}


void PreconditionerIfpack::ComputeInverse()
{
  AMANZI_ASSERT(IfpILU_.get());
  IfpILU_->Compute();
}


/* ******************************************************************
 * Rebuild the preconditioner using the given matrix A.
 ****************************************************************** */
void PreconditionerIfpack::UpdateInverse()
{
  AMANZI_ASSERT(initialized_);
  AMANZI_ASSERT(h_.get());

  Ifpack factory;
  std::string method = plist_.get<std::string>("method");

  // deprecated
  if (method == "block ilu") {
    method = "ILU";
  }
    
  int overlap = plist_.get<int>("overlap", 0);
  // set Amanzi defaults
  if (!plist_.isParameter("schwarz: combine mode"))
    plist_.set<std::string>("schwarz: combine mode", "Add");

  IfpILU_ = Teuchos::rcp(factory.Create(method, &*h_, overlap));
  IfpILU_->SetParameters(plist_);
  IfpILU_->Initialize();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    IfpILU_->Print(*vo_->os());
  }
}

std::string PreconditionerIfpack::returned_code_string() const
{
  switch (returned_code()) {
    case -1 :
      return "Generic Ifpack error.";
    case -2 :
      return "Ifpack says input data not valid.";
    case -3 :
      return "Ifpack says data not correctly preprocessed.";
    case -4 :
      return "Ifpack says problem encountered during algorithm, e.g. divide-by-zero, out-of-bounds, etc.";
    case -5 :
      return "Ifpack out-of-memory";
    case 0:
      return "Ifpack not yet applied.";
    case 1:
      return "success";
  }
  return "unknown error";
}


}  // namespace AmanziSolvers
}  // namespace Amanzi
