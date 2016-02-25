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
#include "Point.hh"

#include "ColorFunction.hh"
#include "ColorFunctionFactory.hh"

#include "RegionColorFunction.hh"

namespace Amanzi {
namespace AmanziGeometry {

//
// RegionColorFunction:: constructor
// -------------------------------------------------------------
RegionColorFunction::RegionColorFunction(const std::string& name, 
                                         const Set_ID id,
                                         const std::string& file,
                                         const int value,
                                         const Epetra_MpiComm *comm,
                                         const LifeCycleType lifecycle)
  : Region(name, id, true, COLORFUNCTION, 0, 0, lifecycle),
    file_(file),
    value_(value)
{
  // Region dimension is set arbitrarily as 3 since the set of
  // entities in the mesh will determine the dimension

  ColorFunctionFactory colfunc_factory;
  colorfunc_ = Teuchos::rcp(colfunc_factory.Create(file_,*comm));
}


// -------------------------------------------------------------
// RegionColorFunction::inside
// -------------------------------------------------------------
bool
RegionColorFunction::inside(const Point& p) const
{
  int color = (*colorfunc_)(&(p[0]));
  return (color == value_);
}
  

} // namespace AmanziGeometry
} // namespace Amanzi
