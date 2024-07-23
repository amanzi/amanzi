/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Mesh Functions

  Interpolate a depth-based, 1D column of data onto a mesh.  Values are
  prescribed only to cells.  Expected is an HDF5 file in the format:

  Depth coordinates z:

  z[:] = (z_0, z_1, ... , z_n)
     z_0 = 0.0
     z_n >= max_z_coordinate of mesh
     z_i > z_(i-1)

  Function values u:

  u[:] = (u_0(z_0), u_1(z_1), ..., u_n(z_n))
*/

#include "dbc.hh"
#include "errors.hh"

#include "CompositeVector.hh"
#include "FunctionFactory.hh"
#include "Function.hh"

#include "ColumnMeshFunction.hh"

namespace Amanzi {
namespace Functions {

void
ReadColumnMeshFunction(Teuchos::ParameterList& plist, CompositeVector& v)
{
  // get filename, data names
  if (!plist.isParameter("file")) {
    Errors::Message message("Missing ColumnMeshFunction parameter \"file\"");
    Exceptions::amanzi_throw(message);
  }
  std::string filename = plist.get<std::string>("file");
  std::string z_str = plist.get<std::string>("z header", "/z");
  std::string f_str = plist.get<std::string>("f header");

  // Create the function
  Teuchos::ParameterList func_list;
  Teuchos::ParameterList& func_sublist = func_list.sublist("function-tabular");
  func_sublist.set("file", filename);
  func_sublist.set("x header", z_str);
  func_sublist.set("y header", f_str);
  FunctionFactory fac;
  auto func = fac.Create(func_list);

  ReadColumnMeshFunction_ByDepth(*func, v);
}


void
ReadColumnMeshFunction_ByDepth(const Function& func, CompositeVector& v)
{
  Epetra_MultiVector& vec = *v.ViewComponent("cell");
  const AmanziMesh::Mesh& mesh = *v.Mesh();
  int d = mesh.getSpaceDimension() - 1;

  int num_columns = mesh.columns.num_columns_owned;
  for (int col = 0; col != num_columns; ++col) {
    std::vector<double> z(1, 0.0);

    const auto& col_faces = mesh.columns.getFaces(col);
    const auto& col_cells = mesh.columns.getCells(col);
    double top_z = mesh.getFaceCentroid(col_faces[0])[d];
    for (int i = 0; i != col_cells.size(); ++i) {
      double cell_z =
        (mesh.getFaceCentroid(col_faces[i])[d] + mesh.getFaceCentroid(col_faces[i + 1])[d]) / 2;
      z[0] = top_z - cell_z;
      vec[0][col_cells[i]] = func(z);
    }
  }
}

} // namespace Functions
} // namespace Amanzi
