/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Collection of non-member functions f2 = f2 @ Map(f1) where 
  Map() connects fields living on different geometric objects.
*/

#ifndef AMANZI_OPERATORS_FIELD_MAPS_HH_
#define AMANZI_OPERATORS_FIELD_MAPS_HH_

#include "Teuchos_RCP.hpp"

namespace Amanzi {

class CompositeVector;

namespace Operators {

int CellToFace_Scale(Teuchos::RCP<CompositeVector> f1, Teuchos::RCP<CompositeVector>& f2);
int CellToFace_ScaleInverse(Teuchos::RCP<const CompositeVector> f1, Teuchos::RCP<CompositeVector>& f2);

}  // namespace Operators
}  // namespace Amanzi


#endif
