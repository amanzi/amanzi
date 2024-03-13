/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

//! RegionFunctionColor: A region defined by the value of an indicator function in a file.
/*!

The list *region: color function* defines a region based a specified integer
color, *value*, in a structured color function file, *file*.  The format of
the color function file is given below in the "Tabulated function file format"
section. As shown in the file, the color values may be specified at the nodes
or cells of the color function grid. A computational cell is assigned the
'color' of the data grid cell containing its cell centroid (cell-based colors)
or the data grid nearest its cell-centroid (node-based colors). Computational
cells sets are then built from all cells with the specified color *Value*.

In order to avoid, gaps and overlaps in specifying materials, it is strongly
recommended that regions be defined using a single color function file.

.. _region-color-function-spec:
.. admonition:: region-color-function-spec

    * `"file`" ``[string]`` File name containing color function.
    * `"value`" ``[int]`` Color that defines the set in the tabulated function file.

Example:

.. code-block:: xml

   <ParameterList name="SOIL_TOP">
     <ParameterList name="region: color function">
       <Parameter name="file" type="string" value="geology_resamp_2D.tf3"/>
       <Parameter name="value" type="int" value="1"/>
     </ParameterList>
   </ParameterList>

*/

#ifndef AMANZI_REGION_COLOR_FUNCTION_HH_
#define AMANZI_REGION_COLOR_FUNCTION_HH_

#include "AmanziTypes.hh"

#include "Region.hh"

namespace Amanzi {

class FunctionColor;

namespace AmanziGeometry {

class RegionFunctionColor : public Region {
 public:
  // Constructor
  RegionFunctionColor(const std::string& name,
                      const int id,
                      const std::string& file,
                      const int value,
                      const Comm_type& comm,
                      const LifeCycleType lifecycle = LifeCycleType::PERMANENT);

  // Is the the specified point inside this region
  bool inside(const Point& p) const;

 protected:
  std::string file_; // which file are we supposed to read it from
  const int value_;
  Teuchos::RCP<FunctionColor> colorfunc_; // indicator func created from file
};

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
