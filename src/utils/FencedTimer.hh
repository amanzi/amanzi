/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/
#pragma once

#include "Kokkos_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Amanzi {

struct FencedTimer {
  explicit FencedTimer(Teuchos::Time& t) : mon_(t) {}
  ~FencedTimer() { Kokkos::fence(); }
  Teuchos::TimeMonitor mon_;
};

} // namespace Amanzi

#define AMANZI_TIMER_CAT2(a, b) a##b
#define AMANZI_TIMER_CAT(a, b) AMANZI_TIMER_CAT2(a, b)
#define AMANZI_TIMER(name)                                                        \
  static Teuchos::RCP<Teuchos::Time> AMANZI_TIMER_CAT(amanzi_timer_, __LINE__) = \
    Teuchos::TimeMonitor::getNewCounter(name);                                    \
  Amanzi::FencedTimer AMANZI_TIMER_CAT(amanzi_fenced_, __LINE__)(                 \
    *AMANZI_TIMER_CAT(amanzi_timer_, __LINE__))
