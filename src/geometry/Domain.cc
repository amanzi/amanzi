/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Rao Garimella (rao@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "Domain.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class Domain
// -------------------------------------------------------------

// -------------------------------------------------------------
// Domain:: constructors / destructor
// -------------------------------------------------------------

// Constructor
Domain::Domain(const unsigned int dim) : spatial_dimension_(dim)
{
  if (dim < 2u || dim > 3u) {
    std::cerr << "Only 2D and 3D domains are supported" << std::endl;
    throw std::exception();
  }

  GeometricModels.clear();
  FreeRegions.clear();
}


// Copy constructor
Domain::Domain(const Domain& old)
{
  int i, ng, nr;

  spatial_dimension_ = old.spatial_dimension();

  ng = old.Num_Geometric_Models();
  for (i = 0; i < ng; i++) {
    Teuchos::RCP<GeometricModel> g = old.Geometric_Model_i(i);
    GeometricModels.push_back(g);
  }

  nr = old.Num_Free_Regions();
  for (i = 0; i < nr; i++) {
    Teuchos::RCP<Region> r = old.Free_Region_i(i);
    FreeRegions.push_back(r);
  }
}


// Destructor
Domain::~Domain(void)
{
  GeometricModels.clear();
  FreeRegions.clear();
}


// Constructor with lists of geometric models and free regions
Domain::Domain(
  const unsigned int dim,
  const std::vector<Teuchos::RCP<GeometricModel>>& in_geometric_models,
  const std::vector<Teuchos::RCP<Region>>& in_Regions)
  : spatial_dimension_(dim),
    GeometricModels(in_geometric_models),
    FreeRegions(in_Regions)
{
  if (dim != 2 || dim != 3) {
    std::cerr << "Only 2D and 3D domains are supported" << std::endl;
    throw std::exception();
  }
}


// Add a geometric model
void
Domain::Add_Geometric_Model(const Teuchos::RCP<GeometricModel>& gm)
{
  // // Make sure spatial dimension of domain and geometric model are the same

  // if (spatial_dimension_ != gm->dimension()) {
  //   std::cerr << "Spatial dimension of domain and geometric model mismatch"
  //   << std::endl; throw std::exception();
  // }

  GeometricModels.push_back(gm);
}


// Add a Free Region
void
Domain::Add_Free_Region(const Teuchos::RCP<Region>& regptr)
{
  if (spatial_dimension_ < regptr->manifold_dimension()) {
    std::cerr
      << "Spatial dimension of domain is less than that of the free region"
      << std::endl;
    throw std::exception();
  }

  FreeRegions.push_back(regptr);
}


// Number of geometric models
int
Domain::Num_Geometric_Models(void) const
{
  return GeometricModels.size();
}


// Get the i'th Decomposition
Teuchos::RCP<GeometricModel>
Domain::Geometric_Model_i(const int i) const
{
  return GeometricModels[i];
}


// Number of Free Regions
int
Domain::Num_Free_Regions(void) const
{
  return FreeRegions.size();
}


// Get the i'th Free Region
Teuchos::RCP<Region>
Domain::Free_Region_i(const int i) const
{
  return FreeRegions[i];
}

} // namespace AmanziGeometry
} // namespace Amanzi
