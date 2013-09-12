/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
           Daniil Svyatskiy
*/

#include <strings.h>

#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"
#include "boundary_function.hh"

#include "Flow_typedefs.hh"
#include "RelativePermeability.hh"


namespace Amanzi {
namespace AmanziFlow {

Teuchos::RCP<Matrix> Flow_PK::MatrixFactory(Teuchos::ParameterList& plist)
{
  std::string name = plist.get<std::string>("matrix", "mfd");
  if (name == "mfd") {
    Teuchos::RCP<MatrixMFD> matrix = Teuchos::rcp(new MatrixMFD(mesh_, super_map_));
    matrix->Init(...);  // Local parameter list goes here
    return matrix;
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

