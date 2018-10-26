#ifndef AMANZI_EXECUTOR_HH_
#define AMANZI_EXECUTOR_HH_

namespace Amanzi {

using AmanziDefaultDevice = Kokkos::Device<Kokkos::CudaUVMSpace::execution_space, Kokkos::CudaUVMSpace::memory_space>


// Here's framework function (not for physics developers to write) that
// runs a model on the given device (GPU or CPU).
template<template<class DeviceType> class Model>
void ExecuteModel(const std::string& kernelName,
		  Model<DeviceType>& model,
		  const int beg, const int end)
{
  using execution_space = typename DeviceType::execution_space;
  using range = Kokkos::RangePolicy<execution_space, int>;
  Kokkos::parallel_for(kernelName, range(beg, end), model);
}

} // namespace Amanzi

#endif
