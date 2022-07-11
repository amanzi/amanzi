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
