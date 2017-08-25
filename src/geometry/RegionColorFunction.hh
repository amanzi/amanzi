/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! RegionColorFunction: A region defined by the value of an indicator function in a file.

/*
  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

/*!

The list *region: color function* defines a region based a specified
integer color, *value*, in a structured color function file,
*file*. The format of the color function file is given below in
the "Tabulated function file format" section. As
shown in the file, the color values may be specified at the nodes or
cells of the color function grid. A computational cell is assigned
the 'color' of the data grid cell containing its cell centroid
(cell-based colors) or the data grid nearest its cell-centroid
(node-based colors). Computational cells sets are then built from
all cells with the specified color *Value*.

In order to avoid, gaps and overlaps in specifying materials, it is
strongly recommended that regions be defined using a single color
function file. 

* `"file`" ``[string]`` File name.

* `"value`" ``[int]`` Color that defines the set in a tabulated function file.

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
