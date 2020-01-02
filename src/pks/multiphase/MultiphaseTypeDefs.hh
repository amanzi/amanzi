/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_MULTIPHASE_TYPEDEFS_HH_
#define AMANZI_MULTIPHASE_TYPEDEFS_HH_

#include "MeshPartition.hh"

#include "WRMmp.hh"

namespace Amanzi {
namespace Multiphase {

typedef std::vector<Teuchos::RCP<WRMmp> > WRMmpList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, WRMmpList> WRMmpPartition;

}  // namespace Flow
}  // namespace Amanzi

#endif

