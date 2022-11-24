/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov
*/

#ifndef AMANZI_PK_EOS_UTILS_HH_
#define AMANZI_PK_EOS_UTILS_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AmanziComm.hh"
#include "errors.hh"

#include "Factory.hh"

namespace Amanzi {
namespace AmanziEOS {

inline void
ErrorAnalysis(const Comm_ptr_type& comm, int ierr, const std::string& msg)
{
  int ierr_glb;
  comm->MaxAll(&ierr, &ierr_glb, 1);

  if (ierr_glb != 0) {
    Errors::CutTimeStep e;
    e << msg;
    amanzi_throw(e);
  }
}

} // namespace AmanziEOS
} // namespace Amanzi

#endif
