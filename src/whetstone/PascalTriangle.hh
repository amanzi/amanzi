/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Pascal's triangle.
*/

#ifndef AMANZI_WHETSTONE_PASCAL_TRIANGLE_HH_
#define AMANZI_WHETSTONE_PASCAL_TRIANGLE_HH_

static const int pascal_triangle[55] = {
  1, 
  1, 1,
  1, 2, 1,
  1, 3, 3, 1,
  1, 4, 6, 4, 1,
  1, 5, 10, 10, 5, 1,
  1, 6, 15, 20, 15, 6, 1,
  1, 7, 21, 35, 35, 21, 7, 1,
  1, 8, 28, 56, 70, 56, 28, 8, 1,
  1, 9, 36, 84, 126, 126, 84, 36, 9, 1,
};

#endif
