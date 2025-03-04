/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! RegionLevelSet: A region defined by a level set function.
/*!
List *region: level set* defines a general region using a level set
function f(x) > 0.

.. admonition:: region_level-set

  * `"dimension`"  ``[int]`` region spatial dimension
  * `"formula`" ``[string]`` level set formula

Example:

.. code-block:: xml

   <ParameterList name="TOP_SECTION"> <!-- parent list -->
     <ParameterList name="region: level set">
       <Parameter name="dimension" type="int" value="2" />
       <Parameter name="formula" type="string" value="1 - (x * x + 2 * y * y)" />
       <ParameterList name="expert parameters">
         <Parameter name="tolerance" type="double" value="1.0e-05"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>

*/


#ifndef AMANZI_REGION_LEVEL_SET_HH_
#define AMANZI_REGION_LEVEL_SET_HH_

#include "ExprTK.hh"

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionLevelSet : public Region {
 public:
  // Default constructor uses point and normal
  RegionLevelSet(const std::string& name,
                 const int id,
                 const int dim,
                 const std::string& formula,
                 const LifeCycleType lifecycle = LifeCycleType::PERMANENT);

  // Is the specified point inside this region 
  bool inside(const Point& p) const;

 protected:
  std::shared_ptr<Utils::ExprTK> exprtk_;
};

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
