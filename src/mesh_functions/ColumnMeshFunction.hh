/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

Interpolate a depth-based, 1D column of data onto a mesh.  Values are
prescribed only to cells.  Expected is an HDF5 or NetCDF file in the format:

Depth coordinates z:

  /z[:] = (z_0, z_1, ... , z_n)
     z_0 = 0.0
     z_n >= max depth of mesh
     z_i > z_(i-1)

Function values u:

  /f[:] = (f_0(z_0), f_1(z_1), ..., f_n(z_n))

.. _column-file-initialization-spec:
.. admonition:: column-file-initialization-spec

   * `"file`" ``[string]`` HDF5 filename
   * `"z header`" ``[string]`` name of the z-coordinate data: `z` above.  Depth
     coordinates (positive downward from the surface), [m]
   * `"f header`" ``[string]`` name of the function data: `f` above.

*/


/*
DEVELOPERS NOTE: this needs a unit test!
*/

#ifndef AMANZI_FUNCTIONS_COLUMN_MESH_FUNCTION_HH_
#define AMANZI_FUNCTIONS_COLUMN_MESH_FUNCTION_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {

class CompositeVector;
class Function;

namespace Functions {

void ReadColumnMeshFunction(Teuchos::ParameterList& plist, CompositeVector& v);

void ReadColumnMeshFunction_ByDepth(const Function& func, CompositeVector& v);

} // namespace Functions
} // namespace Amanzi

#endif
