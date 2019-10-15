/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

/**
 *  \file   Timer.h
 *  \brief  Timer can be used to return the running time of a program in clock
 * cycles.
 *
 *  \date   Apr 14, 2011
 *  \author Nathan Barnett
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <ctime>
#include <string>
#include <iostream>

#include "mpi.h"

namespace Amanzi {

// Forward declare the manager so that it can be a friend
class TimerManager;

/**
 * \class   Timer
 * \author  Nathan Barnett
 */
class Timer {
  friend std::ostream& operator<<(std::ostream& os, Timer& t);
  friend class TimerManager;

 public:
  enum Type { AVERAGE = 0, ACCUMULATE = 1, ONCE = 2, MAXIMUM = 3, MINIMUM = 4 };

  explicit Timer(std::string name = "Timer_", Type type = ONCE);
  void start();
  void stop();
  clock_t getTicks();
  double getTime();
  ~Timer();
  void parSync(MPI_Comm);

  double getAvgTime() { return _avg_elapsed; }
  double getMaxTime() { return _max_elapsed; }
  double getMinTime() { return _min_elapsed; }

 private:
  std::string
    _name; //!< Name of the timer (defaults to "Timer + _numTimerInstances")
  clock_t _startTime; //!< Time at which the timer was started
  clock_t _stopTime;  //!< Time at which the timer was stopped
  bool _running;      //!< Flag for determining if the timer is running
  static unsigned
    _numTimerInstances;  //!< Keeps track of the number of timers instantiated
  unsigned _id;          //!< Number of this timer
  clock_t _runningTotal; //!< Used to keep track of running totals, etc.
  unsigned
    _numInvocations; //!< Number of times the timer has been started/restarted
  Type _type;        //!< Type of the timer (average, accumulate, etc.)

  double _max_elapsed, _min_elapsed, _avg_elapsed;
};

} // namespace Amanzi

#endif /* TIMER_H_ */
