/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Rao Garimella (rao@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_DOMAIN_HH_
#define AMANZI_DOMAIN_HH_

#include <vector>

#include "Region.hh"
#include "GeometricModel.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class Domain
// -------------------------------------------------------------
/// A class that represent the entire geometric domain
/**
 * The geometric domain contains one or more objects describing the specifics
 * of how the domain is decomposed into subdomains. In its simplest form,
 * the geometric domain is a list of geometric regions (see
 * AmanziGeometry::Region) that may or may not tile the domain with or
 * without overlap. In its most sophisticated form, the domain has a list
 * of geometric models, each of which is a particular non-overlapping, tiling
 * decomposition of the domain into geometric regions. Other free floating
 * regions may also exist in the domain for purposes like post-processing
 */

class Domain {
 public:
  // Default constructor.
  Domain(const unsigned int dim);

  // Copy constructor
  Domain(const Domain& old);

  // Constructor with Geometric Model List
  Domain(const unsigned int dim,
         const std::vector<Teuchos::RCP<GeometricModel>>& in_geometric_models,
         const std::vector<Teuchos::RCP<Region>>& in_Regions);

  // Destructor
  virtual ~Domain(void);

  inline unsigned int spatial_dimension() const { return spatial_dimension_; }

  // Add a Geometric Model
  void Add_Geometric_Model(const Teuchos::RCP<GeometricModel>& gm);

  // Add a Free Region
  void Add_Free_Region(const Teuchos::RCP<Region>& regptr);

  // Number of Geometric Models
  int Num_Geometric_Models(void) const;

  // Get the i'th Geometric Model
  Teuchos::RCP<GeometricModel> Geometric_Model_i(const int i) const;

  // Number of Free Regions
  int Num_Free_Regions(void) const;

  // Get the i'th Free Region
  Teuchos::RCP<Region> Free_Region_i(const int i) const;

 private:
  // Dimension of domain
  unsigned int spatial_dimension_;

  // List of geometric models in domain
  std::vector<Teuchos::RCP<GeometricModel>> GeometricModels;

  // List of Free Region pointers (regions that are not part of any
  // geometric model)
  std::vector<Teuchos::RCP<Region>> FreeRegions;
};

} // namespace AmanziGeometry
} // namespace Amanzi

#endif
