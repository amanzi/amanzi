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

#include "dbc.hh"
#include "errors.hh"

#include "AmanziTypes.hh"
#include "Point.hh"

#include "FunctionColor.hh"
#include "FunctionColorFactory.hh"

#include "RegionFunctionColor.hh"

namespace Amanzi {
namespace AmanziGeometry {

//
// RegionFunctionColor:: constructor
// -------------------------------------------------------------
RegionFunctionColor::RegionFunctionColor(const std::string& name,
                                         const int id,
                                         const std::string& file,
                                         const int value,
                                         const Comm_type& comm,
                                         const LifeCycleType lifecycle)
  : Region(name, id, true, RegionType::COLORFUNCTION, 0, 0, lifecycle), file_(file), value_(value)
{
  FunctionColorFactory colfunc_factory;
  colorfunc_ = Teuchos::rcp(colfunc_factory.Create(file_, comm));
  set_space_dimension(colorfunc_->getDimension());
}


// -------------------------------------------------------------
// RegionFunctionColor::inside
// -------------------------------------------------------------
bool
RegionFunctionColor::inside(const Point& p) const
{
  int color = (*colorfunc_)(&(p[0]));
  return (color == value_);
}


} // namespace AmanziGeometry
} // namespace Amanzi
