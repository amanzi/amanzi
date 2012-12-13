#ifndef TimerManager_h
#define TimerManager_h

// C includes

// C++ includes
#include <iostream>
#include <map>

// TPL includes
#include <boost/shared_ptr.hpp>

// Our includes
#include "Timer.hh"

#include "mpi.h"

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
  friend std::ostream& operator <<(std::ostream&, TimerManager&);
 
  TimerManager();
  ~TimerManager();

  void    add(std::string name, Timer::Type type);
  size_t  size();
  void    start();
  void    start(std::string);
  void    stop();
  void    stop(std::string);
  Timer&  operator()(std::string name);
  void    parSync(MPI_Comm comm);

protected:
  std::map<std::string, boost::shared_ptr<Timer> > _timer;
  
private:
  TimerManager(const TimerManager&);
  TimerManager& operator=(const TimerManager&);
};

// Declare a global instance. Messy, but beats modding all functions to take another parameter
extern TimerManager timer_manager;

} //end of namespace Amanzi

#endif // TimerManager_h

