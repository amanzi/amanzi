/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iomanip>
#include <iostream>

#include <cstring>

#include "Epetra_MpiComm.h"
#include "Mesh.hh"
#include "MeshInfo.hh"

namespace Amanzi {

void
MeshInfo::WriteMeshCentroids(std::string domain, const AmanziMesh::Mesh& mesh)
{
  std::string filename = plist_.get<std::string>("filename", "meshinfo");

  Teuchos::RCP<Amanzi::HDF5_MPI> output;
  if (single_file_) {
    output = output_.at("domain");
  } else {
    output = output_.at(domain);
  }

  output->createDataFile(filename);
  output->open_h5file();

  int ncells_owned =
    mesh.getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  int n_glob;
  mesh.getComm()->SumAll(&ncells_owned, &n_glob, 1);

  Epetra_BlockMap map(n_glob, ncells_owned, 1, 0, *mesh.getComm());

  // create an auxiliary vector that will hold the centroid and velocity
  int dim = mesh.getSpaceDimension();
  Teuchos::RCP<Epetra_MultiVector> aux = Teuchos::rcp(new Epetra_MultiVector(map, dim));

  std::vector<std::string> name;
  name.resize(0);
  name.push_back("x");
  name.push_back("y");
  if (dim > 2) name.push_back("z");


  for (int n = 0; n < ncells_owned; n++) {
    const AmanziGeometry::Point& xc = mesh.getCellCentroid(n);
    for (int i = 0; i < dim; i++) { (*(*aux)(i))[n] = xc[i]; }
  }

  WriteVector(*aux, name);
  output->close_h5file();
}

} // namespace Amanzi
