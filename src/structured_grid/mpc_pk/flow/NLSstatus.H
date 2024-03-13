/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _NLSstatus_H_
#define _NLSstatus_H_

class NLSstatus
{
public:
  enum NLSstatusVal {
    NLS_SUCCESS           = 0,
    NLS_NONLINEAR_FAIL    = 1,
    NLS_LINEAR_FAIL       = 2,
    NLS_CATASTROPHIC_FAIL = 3,
    NLS_INVALID           = 4};

  NLSstatus(NLSstatusVal _val = NLS_INVALID)
    : val(_val)
  {
    BL_ASSERT(val>=0 && val<=4);
  }
  bool operator==(const NLSstatus& rhs) {return val == rhs.val;}
  bool operator<(const NLSstatus& rhs) {return val < rhs.val;}
  bool operator!=(const NLSstatus& rhs) {return val != rhs.val;}
  NLSstatusVal val;
};

#endif
