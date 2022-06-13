#include <iomanip>
#include <iostream>

#include <cstring>

#include "Epetra_MpiComm.h"
#include "Mesh.hh"
#include "MeshInfo.hh"

namespace Amanzi {

void MeshInfo::WriteMeshCentroids(std::string domain, const AmanziMesh::Mesh& mesh)
{
  std::string filename = plist_.get<std::string>("filename", "meshinfo");
  AMANZI_ASSERT(output_.count(domain));
  output_.at(domain)->createDataFile(filename);
  output_.at(domain)->open_h5file();

  int ncells_owned = mesh.num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int n_glob;
  mesh.get_comm()->SumAll(&ncells_owned, &n_glob, 1);

  Epetra_BlockMap map(n_glob, ncells_owned, 1, 0, *mesh.get_comm());

  // create an auxiliary vector that will hold the centroid and velocity
  int dim = mesh.space_dimension();
  Teuchos::RCP<Epetra_MultiVector> aux = Teuchos::rcp(new Epetra_MultiVector(map, dim));

  std::vector<std::string> name;
  name.resize(0);
  name.push_back("x");
  name.push_back("y");
  if (dim > 2) name.push_back("z");


  for (int n = 0; n < ncells_owned; n++) {
    const AmanziGeometry::Point& xc = mesh.cell_centroid(n);
    for (int i = 0; i < dim; i++) {
      (*(*aux)(i))[n] = xc[i];
    }
  }

  WriteVector(*aux, name);

  output_.at(domain)->close_h5file();
}

} // namespace Amanzi
