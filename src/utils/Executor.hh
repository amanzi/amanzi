/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#ifndef AMANZI_EXECUTOR_HH_
#define AMANZI_EXECUTOR_HH_

#include "Kokkos_Core.hpp"
#include "AmanziTypes.hh"

namespace Amanzi {

using Entity_ID_View = Kokkos::View<LO*>;

//
// ExecuteModel
// ------------------------------------------------------------
//
// These utility functions execute models over a variety of range policies
//

// Execute over continuous range, no tag.
template <template <typename> class Model,
          class DeviceType = DefaultDevice>
void
ExecuteModel(const std::string& kernelName, Model<DeviceType>& model,
             const int beg, const int end, DeviceType /* dev */ = DeviceType())
{
  using execution_space = typename DeviceType::execution_space;
  using range = Kokkos::RangePolicy<execution_space, int>;
  Kokkos::parallel_for(
    "Evaluator::ExecuteModel",
    kernelName, range(beg, end), model);
}

// Execute over a continuous range, with tag.
template <template <typename> class Model, class TagType,
          class DeviceType = DefaultDevice>
void
ExecuteModel(const std::string& kernelName, Model<DeviceType>& model,
             TagType /* tag */, const int beg, const int end,
             DeviceType /* dev */ = DeviceType())
{
  using execution_space = typename DeviceType::execution_space;
  using range = Kokkos::RangePolicy<execution_space, TagType, int>;
  Kokkos::parallel_for(
    "Evaluator::ExecuteModel",
    kernelName, range(beg, end), model);
}


// Execute over a random range, no tag.
template <template <typename> class Model,
          class DeviceType = DefaultDevice>
void
ExecuteModel(const std::string& kernelName, Model<DeviceType>& model,
             const Entity_ID_View& entities,
             DeviceType /* dev */ = DeviceType())
{
  using execution_space = typename DeviceType::execution_space;
  using range = Kokkos::RangePolicy<execution_space, int>;
  Kokkos::parallel_for(
    "Evaluator::ExecuteModel",
    kernelName,
                       range(0, entities.extent(0)),
                       KOKKOS_LAMBDA(const int i) { model(entities(i)); });
}

// Execute over a random range, with tag.
template <template <typename> class Model, class TagType,
          class DeviceType = DefaultDevice>
void
ExecuteModel(const std::string& kernelName, Model<DeviceType>& model,
             TagType /* tag */, const Entity_ID_View& entities,
             DeviceType /* dev */ = DeviceType())
{
  using execution_space = typename DeviceType::execution_space;
  using range = Kokkos::RangePolicy<execution_space, int>;
  Kokkos::parallel_for(
    "Evaluator::ExecuteModel",
    kernelName, range(0, entities.extent(0)), KOKKOS_LAMBDA(const int i) {
      model(TagType(), entities(i));
    });
}


} // namespace Amanzi

#endif
