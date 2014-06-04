/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "OperatorDefs.hh"
#include "OperatorDiffusionSurface.hh"


namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void OperatorDiffusionSurface::InitDiffusionSurface_(Teuchos::ParameterList& plist)
{
  // factor_ = 1.0;
}

}  // namespace AmanziFlow
}  // namespace Amanzi

