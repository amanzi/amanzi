/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   LogicalRegion.hh
 * @author Rao Garimella
 * @date 
 * 
 * @brief  Declaration of Logical Region class which derives 
 *         an operation on one or two other regions
 * 
 * 
 */

#ifndef _LogicalRegion_hh_
#define _LogicalRegion_hh_

#include "Region.hh"

namespace Amanzi {

  namespace AmanziGeometry {

// -------------------------------------------------------------
//  class LogicalRegion
// -------------------------------------------------------------
/// A region defined by a logical operation on one or two other regions
///
/// Operations supported on a single region are only the NOT operation
///
/// Operations supported on a pair of regions are UNION, SUBTRACT and
/// INTERSECT
///

class LogicalRegion : public Region {
public:

  /// Default constructor 

  LogicalRegion(const Set_Name& name, 
                const Set_ID id, 
                const std::string operation_str,
                const std::vector<std::string> region_names,
                const LifeCycleType lifecycle=PERMANENT,
                const VerboseObject *verbobj=NULL);


  /// Protected copy constructor to avoid unwanted copies.
  LogicalRegion(const LogicalRegion& old);

  /// Destructor
  ~LogicalRegion(void);

  // Type of the region
  inline RegionType type() const { return LOGICAL; }

  // Label in the file
  inline BoolOpType operation() const { return operation_; }

  /// Is the the specified point inside this region
  bool inside(const Point& p) const;

  inline std::vector<std::string> const &component_regions() const 
  { return region_names_; }


protected:  
  BoolOpType operation_; // what logical operation should be performed
  const std::vector<std::string> region_names_; // names of regions in operation
};

/// A smart pointer to LogicalRegion instances
// typedef Teuchos::RCP<LogicalRegion> LogicalRegionPtr;

// RVG: I am not able to correctly code a region factory using smart
// pointers so I will revert to a simpler definition

typedef LogicalRegion *LogicalRegionPtr;

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
