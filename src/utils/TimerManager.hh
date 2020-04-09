/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Nathan Barnett
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_TIMER_MANAGER_HH_
#define AMANZI_TIMER_MANAGER_HH_

// C++ includes
#include <iostream>
#include <map>
#include <memory>

// Our includes
#include "mpi.h"
#include "Timer.hh"

namespace Amanzi {

/**
 * \class  TimerManager
 * \brief  This class provides a method of grouping and managing timers
 *         throughout the code, allowing for a pretty-print off all timers
 *         at the same location/time.
 * \author Nathan Barnett
 * \date   April 13, 2012
 */

class TimerManager {
 public:
  friend std::ostream& operator<<(std::ostream&, TimerManager&);

  TimerManager(){};
  ~TimerManager(){};

  void add(std::string name, Timer::Type type);
  size_t size();
  void start();
  void start(std::string);
  void stop();
  void stop(std::string);
  Timer& operator()(std::string& name);
  void parSync(MPI_Comm comm);
  void print();

 protected:
  std::map<std::string, std::shared_ptr<Timer>> _timer;

 private:
  TimerManager(const TimerManager&);
  TimerManager& operator=(const TimerManager&);
};


// Declare a global instance. Messy, but beats modding all functions to take
// another parameter
extern TimerManager timer_manager;

} // namespace Amanzi

#endif
