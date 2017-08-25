/*
  Reader for loading data from ExodusII into a vector.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  NOTE: this needs a unit test!
*/


#ifndef AMANZI_FUNCTIONS_EXODUSII_MESH_FUNCTION_HH_
#define AMANZI_FUNCTIONS_EXODUSII_MESH_FUNCTION_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Functions {

void ReadExodusIIMeshFunction(Teuchos::ParameterList& file_list,
        CompositeVector& v);

} // namespace Functions
} // namespace Amanzi
                        
#endif
