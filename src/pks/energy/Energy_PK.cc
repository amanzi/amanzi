/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"

#include "GMVMesh.hh"
#include "Mesh.hh"
#include "mfd3d.hh"
#include "OperatorDefs.hh"
#include "State.hh"

#include "Energy_PK.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* default constructor that initializes all pointers to NULL
****************************************************************** */
Energy_PK::Energy_PK() :
    vo_(NULL),
    passwd_("state")
{
}

}  // namespace Energy
}  // namespace Amanzi

