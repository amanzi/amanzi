/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   GeometricModel.hh
 * @author Rao Garimella
 * @date   Sep 15, 2011
 * 
 * @brief  Declaration of the GeometricModel class
 * 
 * 
 */

#ifndef _GeometricModel_hh_
#define _GeometricModel_hh_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"

#include "Region.hh"
#include "RegionFactory.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class GeometricModel
// -------------------------------------------------------------
// A class that represent a geometric model or more specifically, 
// a particular decomposition of the domain into subdomains
/**
 * The geometric model is an object that contains a list of
 * geometric regions that tile the domain (no gaps, no overlaps)
 **/

class GeometricModel {
public:

  // constructor.

  GeometricModel(const unsigned int dim);

  // Copy constructor 

  GeometricModel(const GeometricModel& old);


  // Constructor from parameter list

  GeometricModel(const unsigned int dim, Teuchos::ParameterList gm_param_list,
                 const Epetra_MpiComm *comm);


  // Constructor from a list of regions

  GeometricModel(const unsigned int dim, const std::vector<RegionPtr>& in_Regions); 


  // Destructor

  ~GeometricModel(void);


  // Topological Dimension of geometric model

  inline
  unsigned int dimension() const
  {
    return topo_dimension_;
  }


  // Add a Region to a GeometricModel

  void Add_Region(const RegionPtr& r);


  // Number of Regions

  int Num_Regions(void) const;


  // Get the i'th region of the model

  RegionPtr Region_i(const int i) const;


  // Get a region by its ID
  RegionPtr FindRegion(const int id) const;


  // Get a region by its ID
  RegionPtr FindRegion(const std::string name) const;


  // Check if regions cover the domain extents.  
  // This will work perfectly for domains with rectangular regions
  // but not so for other types of regions

  bool Rough_Check_Tiling(void) const;

private:

  // Topological dimension of the model

  unsigned int topo_dimension_;

  // List of regions in this geometric model

  std::vector<RegionPtr> Regions;

};



  // Smart pointer to an instance of the GeometricModel class
  // RVG: Someone with better C++ knowledge than me could make this work
  //
  // typedef Teuchos::RCP<GeometricModel> GeometricModelPtr;

  typedef GeometricModel *GeometricModelPtr;

} // namespace AmanziGeometry
} // namespace Amanzi

#endif

