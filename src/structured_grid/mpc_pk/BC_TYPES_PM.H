/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _BC_TYPES_PM_H_
#define _BC_TYPES_PM_H_

#if 0
---------------------------------------------------------------------
          CONVERSIONS FROM PHYSICAL TO MATHEMATICAL BCs


PHYS BC        |   Velocity      |    Temperature   |    Scalars     |
----------------------------------------------------------------------
INTERIOR       |                 |                  |                |
   or          |    INT_DIR      |     INT_DIR      |    INT_DIR     |
PERIODIC       |                 |                  |                |
----------------------------------------------------------------------
               |                 |                  |                |
INFLOW         |    EXT_DIR      |     EXT_DIR      |    EXT_DIR     |
               |                 |                  |                |
----------------------------------------------------------------------
               |                 |                  |                |
OUTFLOW        |    FOEXTRAP     |     FOEXTRAP     |    FOEXTRAP    |
               |                 |                  |                |
----------------------------------------------------------------------
NO SLIP WALL   |                 |                  |                |
  with         |    EXT_DIR      |  REFLECT_EVEN    |    HOEXTRAP    |
ADIABATIC TEMP |    U = V = 0    |    dT/dn = 0     |                |
----------------------------------------------------------------------
NO SLIP WALL   |                 |                  |                |
  with         |    EXT_DIR      |     EXT_DIR      |    HOEXTRAP    |
FIXED TEMP     |    U = V = 0    |                  |                |
----------------------------------------------------------------------
SLIP WALL      |                 |                  |                |
  with         |  Un=0 EXT_DIR   |  REFLECT_EVEN    |    HOEXTRAP    |
ADIABATIC TEMP |  Ut   HOEXTRAP  |    dT/dn = 0     |                |
----------------------------------------------------------------------
SLIP WALL      |                 |                  |                |
  with         |  Un=0 EXT_DIR   |     EXT_DIR      |    HOEXTRAP    |
FIXED TEMP     |  Ut   HOEXTRAP  |                  |                |
----------------------------------------------------------------------
               |                 |                  |                |
SYMMETRY       |  Un REFLECT_ODD |  REFLECT_EVEN    | REFLECT_EVEN   |
               |  Ut REFLECT_EVEN|                  |                |
----------------------------------------------------------------------


INT_DIR      : data taken from other grids or interpolated
EXT_DIR      : data specified on EDGE (FACE) of bndry
HOEXTRAP     : higher order extrapolation to EDGE of bndry
FOEXTRAP     : first order extrapolation from last cell in interior
REFLECT_EVEN : F(-n) = F(n)   true reflection from interior cells
REFLECT_ODD  : F(-n) = -F(n)  true reflection from interior cells

The above are defined in amrlib/BC_TYPE.H.  Here, we are only adding
only the seepage boundary.
#endif

#define SEEPAGE 5

#endif
