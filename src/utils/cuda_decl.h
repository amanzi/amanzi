/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#pragma once
// Temporary
// Useful functions to trace NVIDIA calls
#ifdef __CUDACC__
#include <nvToolsExt.h>
#else
inline void nvtxRangePush(const char*){}
inline void nvtxRangePushA(const char*){}
inline void nvtxRangePop(){}
#endif
