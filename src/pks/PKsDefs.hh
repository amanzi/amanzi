/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  PKs

*/

#ifndef AMANZI_PKS_DEFS_HH_
#define AMANZI_PKS_DEFS_HH_

namespace Amanzi {

typedef enum {
  NONE,
  SIMPLE,
  SIMPLE_WELL,
  VOLUME,
  VOLUME_FRACTION,
  COUPLING,
  FIELD,
  FIRST_ORDER_EXCHANGE,
  SUBGRID,
  SUBGRID_RETURN,
  WEIGHT,
  WEIGHT_BY_FIELD
} DomainFunctionType;

} // namespace Amanzi

#endif
