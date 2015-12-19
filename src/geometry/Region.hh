/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   Region.hh
 * @author William A. Perkins
 * @date Mon Aug  1 09:57:42 2011
 * 
 * @brief  Declaration of the abstract Region class 
 * 
 * 
 */

#ifndef AMANZI_REGION_HH_
#define AMANZI_REGION_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"

#include "GeometryDefs.hh"
#include "Point.hh"

namespace Amanzi {
namespace AmanziGeometry {


// -------------------------------------------------------------
//  class Region
// -------------------------------------------------------------
// A class to represent a geometric region
/**
 * A Region is just some arbitrary subset of space, that can be
 * specified in a myriad of ways.  At a minimum, there is a need to be
 * able to determine if a point is inside that space.  Other needs to
 * be added later.
 *
 * The region class does not use a constructor based on the XML parameter
 * list because it has to create derived region classes based on the shape 
 * parameter of the region specification.
 *
 * 
 */

class Region {
public:

  /// Default constructor.
  Region(void);

  /// Constructor with name and ID
  Region(const Set_Name& name, const Set_ID id,
         const unsigned int dim=3, const LifeCycleType lifecycle=PERMANENT,
         const VerboseObject *verbobj=NULL);

  /// Copy constructor 
  Region(const Region& old);

  /// Destructor
  virtual ~Region(void);


  /// Set the dimension of the region
  inline
  void set_dimension(const unsigned int dim)
  {
    topo_dimension_ = dim;
  }

  /// Name of the region
  inline
  std::string name() const
  {
    return name_;
  }

  /// Integer identifier of the region
  inline
  Set_ID id() const
  {
    return id_;
  }

  // Topological dimension of region (0 - point, 1 - curve, 2 - surface, 3 - volume)
  inline 
  unsigned int dimension(void) const
  {
    return topo_dimension_;
  }

  // Get the Lifecycle of this region - Do mesh entity sets derived from
  // it have to be kept around or are they temporary and can be destroyed
  // as soon as they are used?
  
  inline
  LifeCycleType lifecycle(void) const 
  {
    return lifecycle_;
  }

  // Get object encoding verbosity of diagnostic messages and output stream

  inline
  const VerboseObject *verbosity_obj(void) const {
    return verbosity_obj_;
  }

  // Type of the region
  virtual RegionType type() const = 0;

  /// Is the specified point inside the Region
  /// Does being on the boundary count as inside or not?
  virtual bool inside(const Point& p) const = 0;


private:

  // Object encoding output stream and verbosity of diagnostics
  const VerboseObject *verbosity_obj_;

  // Lifecycle (Temporary or Permanent)
  LifeCycleType lifecycle_;
  
  // Topological dimension of region (0, 1, 2, 3)
  unsigned int topo_dimension_;

  // Name of identifier
  Set_Name name_;

  // Integer identifier of region
  Set_ID id_;

};

// Useful typedefs
typedef Region* RegionPtr;
typedef std::vector< RegionPtr > RegionVector;

} // namespace AmanziGeometry
} // namespace Amanzi

#endif

