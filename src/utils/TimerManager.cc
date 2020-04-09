/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Nathan Barnett
*/

//! <MISSING_ONELINE_DOCSTRING>

// #include <utility>

#include "TimerManager.hh"

namespace Amanzi {

// Instanciate global
TimerManager timer_manager;

/**
 * \fn         add
 * \brief      Add a timer to be managed
 * \param[in]  Timer
 * \returns    void
 * \author     Nathan Barnett
 */
void
TimerManager::add(std::string name, Timer::Type type)
{
  _timer.insert(
    std::make_pair(name, std::shared_ptr<Timer>(new Timer(name, type))));
}


/**
 * \fn         stop
 * \brief      Stops the pecified timer
 * \param[in]  string - name of timer
 * \author     Nathan Barnett
 */
void
TimerManager::stop(std::string timerName)
{
  std::map<std::string, std::shared_ptr<Timer>>::iterator it =
    _timer.find(timerName);
  if (it != _timer.end())
    it->second->stop();
  else
    throw "Unkown timer";
}


/**
 * \fn         start
 * \brief      Stops the specified timer
 * \param[in]  string - name of timer
 * \author     Nathan Barnett
 */
void
TimerManager::start(std::string timerName)
{
  std::map<std::string, std::shared_ptr<Timer>>::iterator it =
    _timer.find(timerName);
  if (it != _timer.end())
    it->second->start();
  else
    throw "Unkown timer";
}


/**
 * \fn         stop
 * \brief      Stops all timers being managed
 * \author     Nathan Barnett
 */
void
TimerManager::stop()
{
  for (std::map<std::string, std::shared_ptr<Timer>>::iterator it =
         _timer.begin();
       it != _timer.end();
       ++it)
    it->second->stop();
}


/**
 * \fn         start
 * \brief      Starts all timers being managed
 * \author     Nathan Barnett
 */
void
TimerManager::start()
{
  for (std::map<std::string, std::shared_ptr<Timer>>::iterator it =
         _timer.begin();
       it != _timer.end();
       ++it)
    it->second->stop();
}


/**
 * \fn         getNumTimers
 * \brief      Returns the number of timers currently managed
 * \returns    size_t Number of Timer objects
 * \author     Nathan Barnett
 */
size_t
TimerManager::size()
{
  return _timer.size();
}


/**
 * \fn         operator()
 * \brief      Returns reference to requested timer
 * \param[in]  string - timer name
 * \returns    Reference to requested timer
 * \author     Nathan Barnett
 */
Timer&
TimerManager::operator()(std::string& timerName)
{
  std::map<std::string, std::shared_ptr<Timer>>::iterator it =
    _timer.find(timerName);
  if (it == _timer.end()) throw "Unkown timer";
  return *(it->second);
}


/**
 * \fn         ostream << operator
 * \brief      friend function for std::ostream for output
 * \param[in]  std::ostream&
 * \param[in]  TimerManager&
 * \returns    std::ostream&
 * \author     Nathan Barnett
 */
std::ostream&
operator<<(std::ostream& os, TimerManager& tm)
{
  os << "**********************************************************\n";
  os << "***                   Timing Summary                   ***\n";
  os << "**********************************************************\n\n";

  // Print info for each of the timers
  for (std::map<std::string, std::shared_ptr<Timer>>::iterator it =
         tm._timer.begin();
       it != tm._timer.end();
       ++it) {
    os.width(30);
    os.fill('.');
    os << *(it->second) << std::endl;
  }

  return os;
}


/**
 * \fn         print
 * \brief      print times to cout
 * \param[in]  std::ostream&
 */
void
TimerManager::print()
{
  for (std::map<std::string, std::shared_ptr<Timer>>::iterator it =
         _timer.begin();
       it != _timer.end();
       ++it) {
    std::cout << "   ";
    std::cout.width(30);
    std::cout.fill('.');
    std::cout << *(it->second) << std::endl;
  }
}


/**
 * \fn         to be written
 */
void
TimerManager::parSync(MPI_Comm comm)
{
  for (std::map<std::string, std::shared_ptr<Timer>>::iterator it =
         _timer.begin();
       it != _timer.end();
       ++it)
    (it->second)->parSync(comm);
}

} // namespace Amanzi
