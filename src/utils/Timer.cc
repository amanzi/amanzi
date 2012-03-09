/**
 *  \file   Timer.cc
 *  \brief  Timer can be used to return the running time of a program in clock cycles.
 *
 *  \date   Apr 14, 2011
 *  \author Nathan Barnett
 */

#include "Timer.h"
#include <sstream>


namespace Amanzi
{


/**
 *  \fn     Constructor
 *  \brief  Includes an option for naming the timer - for printout
 */
Timer::Timer(std::string name="Timer_") : _start(0), _stop(0), _running(false), _name(name)
{
  // If we used a generic name, concat the instance number so that the timer names are unique
  if (_name=="Timer_")
  {
    std::stringstream ss;
    ss << _name << _instance;
    _name = ss.string();
  }
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
  _stop = clock();
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
 */
double Timer::getTime()
{
	if ( _running ) return (( clock() - _start ) / CLOCKS_PER_SEC);
	else 			return (( _stop - _start ) / CLOCKS_PER_SEC);
}

std::ostream& operator<<(std::ostream& os)
{
  os << _name << " completed in " << getTime() << " seconds.";
}


} //end of namespace Amanzi
