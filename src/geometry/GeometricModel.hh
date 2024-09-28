/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins
           Ethan Coon (ecoon@lanl.gov)
*/

/*!

It is not always possible to extract space dimension from provided data.
Therefore, we require the user to provide simple list *domain* with only 
one parameter *spatial dimension*.

.. admonition:: geometric_model-spec

  * `"spatial dimension`" ``[int]`` defined space dimension. The available 
    values are 2 or 3.

.. code-block:: xml

  <ParameterList name="domain">
    <Parameter name="spatial dimension" type="int" value="2"/>
  </ParameterList>

*/

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
  GeometricModel(unsigned int dim, Teuchos::ParameterList& gm_param_list, const Comm_type& comm);

  unsigned int dimension() const { return dim_; }

  // Add a Region to a GeometricModel
  void AddRegion(const Teuchos::RCP<Region>& r);

  bool HasRegion(const std::string& name) const { return regions_name_.count(name); }

  // Region iterators
  typedef std::vector<Teuchos::RCP<const Region>>::const_iterator RegionConstIterator;
  std::size_t size() const { return regions_.size(); }
  RegionConstIterator begin() const { return regions_.begin(); }
  RegionConstIterator end() const { return regions_.end(); }

  // Get a region by its ID
  Teuchos::RCP<const Region> FindRegion(const int id) const
  {
    if (regions_id_.count(id)) return regions_id_.at(id);
    return Teuchos::null;
  }

  // Get a region by its name
  Teuchos::RCP<const Region> FindRegion(const std::string& name) const
  {
    if (regions_name_.count(name)) return regions_name_.at(name);
    return Teuchos::null;
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
