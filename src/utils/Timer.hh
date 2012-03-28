/**
 *  \file   Timer.h
 *  \brief  Timer can be used to return the running time of a program in clock cycles.
 *
 *  \date   Apr 14, 2011
 *  \author Nathan Barnett
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <ctime>
#include <string>
#include <iostream>


namespace Amanzi 
{


class Timer
{
  friend std::ostream& operator<<(std::ostream& os, Timer& t);

public:
  explicit Timer(std::string name="Timer_");
  void start();
  void stop();
  clock_t getTicks();
  double getTime();
  ~Timer();

private:
  std::string     _name;      //!< Name of the timer (defaults to "Timer + _instance")
  clock_t 	  _start;     //!< Time at which the timer was started
  clock_t 	  _stop;      //!< Time at which the timer was stopped
  bool 		  _running;   //!< Flag for determining if the timer is running
  static unsigned _instances; //!< Keeps track of the number of timers instantiated
  unsigned        _num;       //!< Number of this timer
};

} //end of namespace Amanzi

#endif /* TIMER_H_ */
