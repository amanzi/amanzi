#ifndef AMANZI_EXECUTOR_HH_
#define AMANZI_EXECUTOR_HH_

#include "Kokkos_Core.hpp"
#include "AmanziTypes.hh"

namespace Amanzi {

// Here's framework function (not for physics developers to write) that
// runs a model on the given device (GPU or CPU).
template<template <typename> class Model, 
         class DeviceType = AmanziDefaultDevice>
void ExecuteModel(const std::string& kernelName,
		  Model<DeviceType>& model,
		  const int beg, const int end,
		  DeviceType /* dev */ = DeviceType ())
{
  using execution_space = typename DeviceType::execution_space;
  using range = Kokkos::RangePolicy<execution_space, int>;
  Kokkos::parallel_for(kernelName, range(beg, end), model);
}

#if 1
template<template <typename> class Model, 
	 class TagType,
         class DeviceType = AmanziDefaultDevice>
void ExecuteModel(const std::string& kernelName,
		  Model<DeviceType>& model,
		  TagType /* tag */,
		  const int beg, const int end,
		  DeviceType /* dev */ = DeviceType ())
{
  using execution_space = typename DeviceType::execution_space;
  using range = Kokkos::RangePolicy<execution_space, TagType, int>;
  std::cout << "launching kernel of type " << kernelName << std::endl;
  Kokkos::parallel_for(kernelName, range(beg, end), model);
}
#endif // 0


} // namespace Amanzi

#endif
