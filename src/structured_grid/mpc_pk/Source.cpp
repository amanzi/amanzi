/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <winstd.H>

#include "Source.H"

void
Source::setVal(FArrayBox& fab,
	       Array<Region*> region_array,
	       const Real* dx)
{
  switch (dist_type)
    {
    case 0:
      set_constant_val(fab,region_array,dx);

    default:
      set_constant_val(fab,region_array,dx);
    }
}

void
Source::set_constant_val(FArrayBox& fab,
			 Array<Region*> region_array,
			 const Real* dx)
{
  for (int i = 0; i<id.size(); i++)
    region_array[region]->setVal(fab,val_param[i],id[i],dx,0);
}
