/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  A region defined by the value of an indicator function in a file

  The region will consist of all mesh elements for which the indicator
  function is a particular value at their centroids

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

#ifndef AMANZI_REGION_COLOR_FUNCTION_HH_
#define AMANZI_REGION_COLOR_FUNCTION_HH_

#include "Epetra_MpiComm.h"

#include "Region.hh"

namespace Amanzi {

class ColorFunction;  
  
namespace AmanziGeometry {

class RegionColorFunction : public Region {
 public:

  // Constructor 
  RegionColorFunction(const std::string& name, 
                      const Set_ID id, 
                      const std::string& file,
                      const int value,
                      const Epetra_MpiComm *comm,
                      const LifeCycleType lifecycle=PERMANENT);

  // Is the the specified point inside this region
  bool inside(const Point& p) const;

protected:  
  std::string file_; // which file are we supposed to read it from
  const int value_;
  Teuchos::RCP<ColorFunction> colorfunc_; // indicator func created from file
};

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
