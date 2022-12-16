/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _SOURCE_H_
#define _SOURCE_H_

#include <string>

#include <Array.H>
#include <FArrayBox.H>
#include <Box.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <LevelBld.H>
#include <AmrLevel.H>
#include <Region.H>

/*
  @ManDoc: Class defining all possible type of sources.
*/
class Source
{
public:

  std::string name;
  /* These specify which variable the source is contributing to.
       var_type: comp or tracer.
       var_id  : name listed in the comp or tracer.
                 ALL: all of the components in var_type.
		 USER: user specified routine in PROB_$D.F
       id: based on var_type, var_id, this is converted to index in State_Type.
  */
  std::string var_type;
  std::string var_id;
  Array<int> id;

  /* This specifies the region in which the source is applied.
     The user give the name but this is internally change to corresponding index
     in region_array.
   */
  int region;

  /* This specifies the distribution of the weight functions
     0: constant value
        param[0]: value
     1: linear/hat function from a point, line or plane
        param[0]: value
	param[1-3]: x1,y1,z1
	param[4-6]: x2,y2,z2
     2: quadratic function from a point, line or plane
        param[0]: value
	param[1-3]: x1,y1,z1
	param[4-6]: x2,y2,z2
     3: exponential function from a point, line or plane
        param[0]: value
	param[1]: exponential coefficient
	param[2-4]: x1,y1,z1
	param[5-7]: x2,y2,z2

   */
  int dist_type;
  Array<Real> val_param;
  Array<Real> dist_param;

  void setVal(FArrayBox& fab,
	       Array<Region*> region_array,
	       const Real* dx);

  void set_constant_val(FArrayBox& fab,
			Array<Region*> region_array,
			const Real* dx);

private:

};
#endif /*_SOURCE_H_*/
