/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      William Perkins
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_GEOMETRIC_MODEL_HH_
#define AMANZI_GEOMETRIC_MODEL_HH_

#include <vector>
#include <map>

#include "Teuchos_RCP.hpp"

#include "AmanziTypes.hh"
#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class GeometricModel {
 public:
  // constructor.
  GeometricModel(unsigned int dim);

  // Constructor from parameter list
  GeometricModel(unsigned int dim, Teuchos::ParameterList& gm_param_list,
                 const Comm_type& comm);

  unsigned int dimension() const { return dim_; }

  // Add a Region to a GeometricModel
  void AddRegion(const Teuchos::RCP<Region>& r);

  // Region iterators
  typedef std::vector<Teuchos::RCP<const Region>>::const_iterator
    RegionConstIterator;
  std::size_t RegionSize() const { return regions_.size(); }
  RegionConstIterator RegionBegin() const { return regions_.begin(); }
  RegionConstIterator RegionEnd() const { return regions_.end(); }

  // Get a region by its ID
  Teuchos::RCP<const Region> FindRegion(const int id) const
  {
    return regions_id_.at(id);
  }

  // Get a region by its name
  Teuchos::RCP<const Region> FindRegion(const std::string name) const
  {
    return regions_name_.at(name);
  }

 private:
  // List of regions in this geometric model
  std::vector<Teuchos::RCP<const Region>> regions_;
  std::map<std::string, Teuchos::RCP<const Region>> regions_name_;
  std::map<int, Teuchos::RCP<const Region>> regions_id_;

  unsigned int dim_;
};

} // namespace AmanziGeometry
} // namespace Amanzi

#endif
