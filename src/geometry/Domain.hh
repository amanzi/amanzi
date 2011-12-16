/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   Domain.hh
 * @author Rao Garimella
 * @date   Mon Aug  1 09:57:42 2011
 * 
 * @brief  Declaration of the Domain class
 * 
 * 
 */

#ifndef _Domain_hh_
#define _Domain_hh_

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
         const std::vector<GeometricModelPtr>& in_geometric_models, 
         const std::vector<RegionPtr>& in_Regions); 

  // Destructor

  virtual ~Domain(void);


  inline 
  unsigned int spatial_dimension() const
  {
    return spatial_dimension_;
  }


  // Add a Geometric Model

  void Add_Geometric_Model(const GeometricModelPtr& gm);


  // Add a Free Region

  void Add_Free_Region(const RegionPtr& regptr);


  // Number of Geometric Models

  int Num_Geometric_Models(void) const;


  // Get the i'th Geometric Model

  GeometricModelPtr Geometric_Model_i(const int i) const;


  // Number of Free Regions

  int Num_Free_Regions(void) const;


  // Get the i'th Free Region

  RegionPtr Free_Region_i(const int i) const;

private:

  // Dimension of domain

  unsigned int spatial_dimension_;


  // List of geometric models in domain

  std::vector<GeometricModelPtr> GeometricModels;


  // List of Free Region pointers (regions that are not part of any
  // geometric model)

  std::vector<RegionPtr> FreeRegions;

};

} // namespace AmanziGeometry
} // namespace Amanzi

#endif

