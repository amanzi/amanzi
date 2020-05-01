/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

/**
 *  \file   Timer.cc
 *  \brief  Timer can be used to return the running time of a program in
 *          clock cycles.
 *
 *  \date   Apr 14, 2011
 *  \author Nathan Barnett
 */

#include "Timer.hh"
#include <sstream>

namespace Amanzi {

unsigned Timer::_numTimerInstances = 0;

/**
 *  \fn         Constructor
 *  \brief      Includes an option for naming the timer - for printout
 *  \param[in]  Name of the timer
 */
Timer::Timer(std::string name, Type type)
  : _startTime(0),
    _stopTime(0),
    _running(false),
    _name(name),
    _type(type),
    _runningTotal(0),
    _numInvocations(1),
    _avg_elapsed(0.0),
    _max_elapsed(0.0),
    _min_elapsed(0.0)
{
  // We increment the number of timer instances and grab the current number
  // for our own
  _id = _numTimerInstances;
  _numTimerInstances++;

  // If we used a generic name, concat the instance number so that the timer
  // names are unique
  if (_name == "Timer_") {
    std::stringstream ss;
    ss << _name << _id;
    _name = ss.str();
  }
}


/**
 *  \fn     Destructor
 */
Timer::~Timer()
{
  // This timer got killed - remove one from the list
  _numTimerInstances--;
}


/**
 *  \fn      start
 *  \brief   Records the number of cpu cycles since the program began as the
 *           start member variable
 */
void
Timer::start()
{
  _running = true;
  _numInvocations++;
  _startTime = clock();
}


/**
 *  \fn     stop
 *  \brief  Records the number of cpu cycles since the program began as the
 *          stop member variable.
 */
void
Timer::stop()
{
  _stopTime = clock();
  _running = false;

  switch (_type) {
  case AVERAGE:
    _runningTotal += (_stopTime - _startTime);
    break;
  case ACCUMULATE:
    _runningTotal += (_stopTime - _startTime);
    break;
  case ONCE:
    _runningTotal = (_stopTime - _startTime);
    break;
  case MAXIMUM:
    _runningTotal = std::max(_runningTotal, (_stopTime - _startTime));
    break;
  case MINIMUM:
    _runningTotal = std::min(_runningTotal, (_stopTime - _startTime));
    break;
  default:
    break;
  }
}


/**
 *  \fn      getTicks
 *  \brief   Returns number of cpu cycles.
 *
 *           If the timer has been stopped, function will return the number
 *           of cpu cycles between the start and stop member variables.
 *           Otherwise, will return the difference between the current cpu
 *           cycle count and the start member variable. Note that this is not
 *           the wall time, but the time the cpu actually uses.
 *  \returns Number of clock cycles (clock_t)
 */
clock_t
Timer::getTicks()
{
  if (_running) {
    switch (_type) {
    case ACCUMULATE:
      return (_runningTotal + (clock() - _startTime));
      break;
    case ONCE:
      return (clock() - _startTime);
      break;
    default:
      throw "That type of timer doesn't make sense to grab while running";
      break;
    }
  } else {
    switch (_type) {
    case AVERAGE:
      return (_runningTotal / _numInvocations);
      break;
    case ACCUMULATE:
      return (_runningTotal);
      break;
    case ONCE:
      return (_runningTotal);
      break;
    case MAXIMUM:
      return (_runningTotal);
      break;
    case MINIMUM:
      return (_runningTotal);
      break;
    default:
      break;
    }
  }
  return 0; 
}


/**
 *  \fn      getTime
 *  \brief   Returns cpu time in seconds.
 *
 *           If the timer has been stopped, function will return the number
 *           of seconds between the start and stop member variables. Otherwise,
 *           will return the number of seconds between the current cpu cycle
 *           count and the start member variable. Note that this is not the
 *           wall time, but the time the cpu actually uses.
 *  \returns Number of seconds (double)
 */
double
Timer::getTime()
{
  if (_running) {
    switch (_type) {
    case ACCUMULATE:
      return (_runningTotal + (clock() - _startTime)) /
             static_cast<double>(CLOCKS_PER_SEC);
      break;
    case ONCE:
      return (clock() - _startTime) / static_cast<double>(CLOCKS_PER_SEC);
      break;
    default:
      throw "That type of timer doesn't make sense to grab while running";
      break;
    }
  } else {
    switch (_type) {
    case AVERAGE:
      return (_runningTotal / _numInvocations) /
             static_cast<double>(CLOCKS_PER_SEC);
      break;
    case ACCUMULATE:
      return (_runningTotal / static_cast<double>(CLOCKS_PER_SEC));
      break;
    case ONCE:
      return (_runningTotal / static_cast<double>(CLOCKS_PER_SEC));
      break;
    case MAXIMUM:
      return (_runningTotal / static_cast<double>(CLOCKS_PER_SEC));
      break;
    case MINIMUM:
      return (_runningTotal / static_cast<double>(CLOCKS_PER_SEC));
      break;
    default:
      break;
    }
  }
  return 0; 
}


void
Timer::parSync(MPI_Comm comm)
{
  double elapsed = getTime();

  MPI_Allreduce(&elapsed, &_max_elapsed, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(&elapsed, &_min_elapsed, 1, MPI_DOUBLE, MPI_MIN, comm);
  MPI_Allreduce(&elapsed, &_avg_elapsed, 1, MPI_DOUBLE, MPI_SUM, comm);
  int numpe(1);
  MPI_Comm_size(comm, &numpe);
  _avg_elapsed /= static_cast<double>(numpe);
}


/**
 *  \fn      friend ostream operator<<
 *  \brief   Outputs text to a ostream
 *  \returns Text string as an output stream
 */
std::ostream&
operator<<(std::ostream& os, Timer& t)
{
  if (t._running)
    os << t._name << " is currently at " << t.getTime() << " seconds.";
  else
    os << t._name << " completed in avg=" << t.getAvgTime()
       << ", min=" << t.getMinTime() << ", max=" << t.getMaxTime()
       << " seconds.";

  return os;
}

} // namespace Amanzi
