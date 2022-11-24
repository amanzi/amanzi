/*
  Mesh Functions

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Interpolate a depth-based, 1D column of data onto a mesh.  Values are
  prescribed only to cells.  Expected is an HDF5 file in the format:

  Depth coordinates z:
  
  z[:] = (z_0, z_1, ... , z_n)
     z_0 = 0.0
     z_n >= max_z_coordinate of mesh
     z_i > z_(i-1)

  Function values u:
  
  u[:] = (u_0(z_0), u_1(z_1), ..., u_n(z_n))
  
  NOTE: this needs a unit test!
*/

#ifndef AMANZI_FUNCTIONS_COLUMN_MESH_FUNCTION_HH_
#define AMANZI_FUNCTIONS_COLUMN_MESH_FUNCTION_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {

class CompositeVector;
class Function;

namespace Functions {

void
ReadColumnMeshFunction(Teuchos::ParameterList& plist, CompositeVector& v);

void
ReadColumnMeshFunction_ByDepth(const Function& func,
                               const Teuchos::Array<std::string> sidesets,
                               CompositeVector& v);

} // namespace Functions
} // namespace Amanzi

#endif
