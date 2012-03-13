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


namespace Amanzi
{

unsigned Timer::_instances = 0;

/**
 *  \fn         Constructor
 *  \brief      Includes an option for naming the timer - for printout
 *  \param[in]  Name of the timer
 */
Timer::Timer(std::string name) : _start(0), _stop(0), _running(false), _name(name)
{
  // We increment the number of timer instances and grab the current number 
  // for our own
  _num = _instances;
  _instances++;

  // If we used a generic name, concat the instance number so that the timer 
  // names are unique
  if (_name=="Timer_")
  {
    std::stringstream ss;
    ss << _name << _num;
    _name = ss.str();
  }
}


/**
 *  \fn     Destructor
 */
Timer::~Timer()
{
  // This timer got killed - remove one from the list
  _instances--;
}


/**
 *  \fn      start 
 *  \brief   Records the number of cpu cycles since the program began as the 
 *           start member variable
 */
void Timer::start()
{
  _running = true;
  _start   = clock();
}


/**
 *  \fn     stop
 *  \brief  Records the number of cpu cycles since the program began as the 
 *          stop member variable.
 */
void Timer::stop()
{
  _stop    = clock();
  _running = false;
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
clock_t Timer::getTicks()
{
  if ( _running ) 
    return ( clock() - _start );
  else 
    return ( _stop - _start );
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
double Timer::getTime()
{
  if ( _running ) 
    return (double( clock() - _start ) / (double)CLOCKS_PER_SEC);
  else 
    return (double( _stop - _start ) / (double)CLOCKS_PER_SEC);
}

/**
 *  \fn      friend ostream operator<<
 *  \brief   Outputs text to a ostream
 *  \returns Text string as an output stream
 */
std::ostream& operator<<(std::ostream& os, Timer& t)
{
  if (t._running)
    os  << t._name << " is currently at " << t.getTime() << " seconds.";
  else
    os << t._name << " completed in " << t.getTime() << " seconds.";

  return os;
}


} //end of namespace Amanzi
