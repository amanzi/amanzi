/*
  A rectangular region in space, defined by two corner points and 
  normals to sides.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Lipnikov Konstantin (lipnikov@lanl.gov)
           Rao Garimella (rao@lanl.gov)
*/

#ifndef AMANZI_BOX_VOLUME_FRACTIONS_REGION_HH_
#define AMANZI_BOX_VOLUME_FRACTIONS_REGION_HH_

#include <vector>

#include "Point.hh"
#include "Region.hh"
#include "TensorSimple.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionBoxVolumeFractions : public Region {
 public:
  // Default constructor uses two corner points (order not important)
  // and vector of normals (order is important).
  RegionBoxVolumeFractions(const std::string& name,
                           const Set_ID id,
                           const Point& p0, 
                           const Point& p1,
                           const std::vector<Point>& normals,
                           const LifeCycleType lifecycle=PERMANENT);

  // Is the the specified point inside this region
  bool inside(const Point& p) const;

  // Is the box degenerate - zero length in one or more directions and
  // if so in how many directions?
  bool is_degenerate(int *ndeg) const;

 protected:
  const Point p0_, p1_; // two corners of the box
  const std::vector<Point> normals_;

 private:
  TensorSimple N_;
};

}  // namespace AmanziGeometry
}  // namespace Amanzi

#endif
