/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  A region defined by a logical operation on one or two other regions
 
  Operations supported on a single region are only the NOT operation
 
  Operations supported on a pair of regions are UNION, SUBTRACT and
  INTERSECT

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

#ifndef AMANZI_REGION_LOGICAL_HH_
#define AMANZI_REGION_LOGICAL_HH_

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionLogical : public Region {
 public:

  // constructor 
  RegionLogical(const std::string& name, 
                const Set_ID id, 
                const std::string& operation_str,
                const std::vector<std::string>& component_regions,
                const LifeCycleType lifecycle=PERMANENT);
  
  // Label in the file
  BoolOpType operation() const { return operation_; }

  /// Is the the specified point inside this region
  bool inside(const Point& p) const;

  std::vector<std::string> component_regions() const {
    return component_regions_; }


protected:  
  BoolOpType operation_; // what logical operation should be performed
  std::vector<std::string> component_regions_;
    // names of regions in operation
};

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
