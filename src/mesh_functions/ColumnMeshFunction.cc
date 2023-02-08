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

  // starting surface
  Teuchos::Array<std::string> sidesets;
  if (plist.isParameter("surface sideset")) {
    sidesets.push_back(plist.get<std::string>("surface sideset"));
  } else if (plist.isParameter("surface sidesets")) {
    sidesets = plist.get<Teuchos::Array<std::string>>("surface sidesets");
  } else {
    Errors::Message message(
      "Missing InitializeFromColumn parameter \"surface sideset\" or \"surface sidesets\"");
    Exceptions::amanzi_throw(message);
  }

  ReadColumnMeshFunction_ByDepth(*func, sidesets, v);
}


void
ReadColumnMeshFunction_ByDepth(const Function& func,
                               const Teuchos::Array<std::string> sidesets,
                               CompositeVector& v)
{
  Epetra_MultiVector& vec = *v.ViewComponent("cell");

  double z0;
  AmanziMesh::Double_List z(1);

  for (auto setname : sidesets) {
    auto surf_faces = v.Mesh()->getSetEntities(
      setname, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);

    for (auto f : surf_faces) {
      // Collect the reference coordinate z0
      AmanziGeometry::Point x0 = v.Mesh()->getFaceCentroid(f);
      z0 = x0[x0.dim() - 1];

      // Iterate down the column
      AmanziMesh::Entity_ID_View cells;
      cells = v.Mesh()->getFaceCells(f, AmanziMesh::Parallel_kind::OWNED);
      AMANZI_ASSERT(cells.size() == 1);
      AmanziMesh::Entity_ID c = cells[0];

      while (c >= 0) {
        AmanziGeometry::Point x1 = v.Mesh()->getCellCentroid(c);
        z[0] = z0 - x1[x1.dim() - 1];
        vec[0][c] = func(z);
        assert(false); 
        //c = v.Mesh()->cell_get_cell_below(c);
      }
    }
  }
}

} // namespace Functions
} // namespace Amanzi
