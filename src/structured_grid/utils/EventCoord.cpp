/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <EventCoord.H>
#include <iostream>

typedef EventCoord::Event Event;
typedef EventCoord::CycleEvent CycleEvent;
typedef EventCoord::TimeEvent TimeEvent;

CycleEvent::CycleEvent(const Array<int>& _cycles)
    : type(CYCLES), cycles(_cycles.size())
{
    for (int i=0; i<cycles.size(); ++i) {
        cycles[i] = _cycles[i];
    }
}

CycleEvent::CycleEvent(const CycleEvent& rhs)
    : start(rhs.start), period(rhs.period), stop(rhs.stop), cycles(rhs.cycles.size()), type(rhs.type)
{
    const Array<int>& v=rhs.cycles;
    for (int i=0; i<v.size(); ++i) {
        cycles[i] = v[i];
    }
}

TimeEvent::TimeEvent(const Array<Real>& _times)
    : type(TIMES), times(_times.size())
{
    for (int i=0; i<times.size(); ++i) {
        times[i] = _times[i];
    }
}

TimeEvent::TimeEvent(const TimeEvent& rhs)
    : start(rhs.start), period(rhs.period), stop(rhs.stop), times(rhs.times.size()), type(rhs.type)
{
    const Array<Real>& v=rhs.times;
    for (int i=0; i<v.size(); ++i) {
        times[i] = v[i];
    }
}

void
EventCoord::Register(const std::string& label,
                     const Event*       event)
{
    // FIXME: Do we have to make the distinction here??
    if (event->IsTime()) {
        std::map<std::string,Event*>::iterator it=timeEvents.find(label);
        if (it!=timeEvents.end()) {
            delete it->second;
        }
        timeEvents[label] = event->Clone();
    }
    else if (event->IsCycle()) {
        std::map<std::string,Event*>::iterator it=cycleEvents.find(label);
        if (it!=cycleEvents.end()) {
            delete it->second;
        }
        cycleEvents[label] = event->Clone();
    }
    else {
        BoxLib::Abort("EventCoord::Invalid event");
    }
}

EventCoord::~EventCoord()
{
    for (std::map<std::string,Event*>::const_iterator it=cycleEvents.begin(); it!=cycleEvents.end(); ++it) {
        delete it->second;
    }
    for (std::map<std::string,Event*>::const_iterator it=timeEvents.begin(); it!=timeEvents.end(); ++it) {
        delete it->second;
    }
}

bool
EventCoord::EventRegistered(const std::string& label) const
{
    if (cycleEvents.find(label) != cycleEvents.end()) {
        return true;
    } else if (timeEvents.find(label) != timeEvents.end()) {
        return true;
    }
    return false;
}

const std::map<std::string,EventCoord::Event*>& EventCoord::CycleEvents() const {return cycleEvents;}
const std::map<std::string,EventCoord::Event*>& EventCoord::TimeEvents() const {return timeEvents;}

static TimeEvent NULLEVENT;

const Event& EventCoord::operator[] (const std::string& label) const
{
    std::map<std::string,Event*>::const_iterator cycle_it = cycleEvents.find(label);
    std::map<std::string,Event*>::const_iterator time_it = timeEvents.find(label);
    if (cycle_it != cycleEvents.end()) {
        return *cycle_it->second;
    } else if (time_it != timeEvents.end()) {
        return *time_it->second;
    }
    BoxLib::Error("EventCoord::operator[] : Event not found");
    return NULLEVENT; // to quiet compiler
}

bool
EventCoord::CycleEvent::ThisEventDue(int cycle, int dCycle) const
{
    int dCycle_red = dCycle;

    if (type == SPS_CYCLES)
    {
        if ( (stop > 0  &&  cycle > stop)  ||  (cycle + dCycle < start) )
        {
            return false;
        }

        if ( ( ( stop < 0 ) || ( cycle <= stop ) )
             && ( cycle + dCycle >= start) )
        {
            int cold_interval = (cycle - start)/period;
            int cnew_interval = (cycle + dCycle - start)/period;
            int cnew_boundary = start + cnew_interval*period;
            int cnew_overage = cycle + dCycle - cnew_boundary;

            // cnew is on boundary
            if ( cnew_overage == 0 ) {
                return true;
            }

            int next_boundary = (cold_interval + 1)*period + start;

            if (next_boundary <= cycle + dCycle) {
                dCycle_red = std::min(dCycle, next_boundary - cycle);
                if (dCycle_red <= dCycle) {
                    BL_ASSERT(dCycle_red > 0);
                    return true;
                }
            }
        }
    }
    else {
        for (int i=0; i<cycles.size(); ++i) {
            if (cycles[i] == cycle+dCycle) {
                return true;
            }
        }
    }
    return false;
}

bool
EventCoord::CycleEvent::ThisEventInit(int cycle) const
{
  if (type == SPS_CYCLES) {
    if ( (stop > 0  &&  cycle > stop)  ||  (cycle < start) ) {
      return false;
    }

    if ( ( ( stop < 0 ) || ( cycle <= stop ) ) && ( cycle >= start) ) {
      return true;
    }
  }
  else {
    for (int i=0; i<cycles.size(); ++i) {
      if (cycle == cycles[i]) {
        return true;
      }
    }
  }
  return false;
}

bool
EventCoord::TimeEvent::ThisEventInit(Real t) const
{
    if (type == SPS_TIMES)
    {
        Real teps = period * 1.e-8;

        if ( (stop > 0  &&  t > stop)  ||  (t < start - teps) )
        {
            return false;
        }

        if ( ( ( stop < 0 ) || ( t - teps <= stop ) )
             &&  t + teps >= start )
        {
          return true;
        }
    }
    else {
      for (int i=0; i<times.size(); ++i) {
        Real teps;
        if (times.size()==1) {
          teps = times[0];
        } else if (i>0) {
          teps = times[i] - times[i-1];
        } else if (i<times.size()-1) {
          teps = times[i+1] - times[i];
        } else {
          BoxLib::Abort("bad logic in time event processing");
        }
        teps *= 1.e-8;

        if ( (t - teps <= times[i]) && (t + teps >= times[i]) ) {
          return true;
        }
      }
    }
    return false;
}


bool
EventCoord::TimeEvent::ThisEventDue(Real t, Real dt, Real& dt_red) const
{
    dt_red = dt;
    if (type == SPS_TIMES)
    {
      Real teps = period * 1.e-8;

        if ( (stop > 0  &&  t > stop)  ||  (t + dt < start - teps) )
        {
            return false;
        }

        if ( ( ( stop < 0 ) || ( t + teps <= stop ) )
             &&  t + dt>= start - teps )
        {
            // Is on interval boundary
            int told_interval = (t - start)/period;
            Real told_boundary = start + told_interval*period;
            Real told_overage = t - told_boundary;

            // told is near boundary
            if ( (std::abs(told_overage - period) < teps)
                 || (std::abs(told_overage) < teps) ) {

                if (dt > period) { // this would step over period, cut back
                    dt_red = period;
                    return true;
                }
            }

            int tnew_interval = (t + dt - start)/period;
            Real tnew_boundary = start + tnew_interval*period;
            Real tnew_overage = t + dt - tnew_boundary;


            // tnew is near boundary
            if ( (std::abs(tnew_overage - period) < teps)
                 || (std::abs(tnew_overage) < teps) ) {
                return true;
            }

            Real next_boundary = (told_interval + 1)*period + start;

            if (stop < 0  ||  next_boundary <= stop) {
                dt_red = std::min(dt, next_boundary - t);
                if (dt_red < dt) {
                    BL_ASSERT(dt_red > 0);
                    return true;
                }
            }
        }
    }
    else {

        for (int i=0; i<times.size(); ++i) {
            Real teps;
            if (times.size()==1) {
                teps = times[0];
            } else if (i>0) {
                teps = times[i] - times[i-1];
            } else if (i<times.size()-1) {
                teps = times[i+1] - times[i];
            } else {
                BoxLib::Abort("bad logic in time event processing");
            }
            teps *= 1.e-8;

            if (times[i] > t  &&  times[i] <= t + dt + teps) {
                dt_red = std::min(dt, times[i]-t);
                return true;
            }

            if (std::abs(times[i] - t) < teps) {

              // In the case where the time is identically on a node, we would
              // double-count, because it would be at the end of one interval
              // and the start of another.  Avoid this by picking the one
              // at the high end.
              //
              // There is a special case though.  If t is at a node, and t+dt is
              // at the next node, then the one at the low side would never be
              // triggered.  I cant think of a way around this, but to check for it
              // explicitly.
              //
              for (int j=0; j<times.size(); ++j) {
                if (i!=j  && std::abs(t + dt - times[j]) < teps) {
                  dt_red = dt;
                  return true;
                }
              }

                Real max_dt = 0;
                for (int j=0; j<times.size(); ++j) {
                    if (t < times[j]) {
                        max_dt = std::max(max_dt,times[j]-t);
                    }
                }

                if (max_dt == 0) {
                    return false;
                }
                else {
                    Real min_dt = max_dt;
                    for (int j=0; j<times.size(); ++j) {
                        if (i!=j  && t < times[j]  &&  t + dt > times[j]) {
                            min_dt = std::min(min_dt, times[j] - t);
                        }
                    }
                    dt_red = std::min(dt, min_dt);
                    return dt_red < dt;
                }
                return true;
            }
        }
    }
    dt_red = -1;
    return false;
}


std::pair<Real,Array<std::string> >
EventCoord::NextEvent(Real t, Real dt, int cycle, int dcycle) const
{
    Array<std::string> events;
    Real delta_t = -1;
    for (std::map<std::string,Event*>::const_iterator it=cycleEvents.begin(); it!=cycleEvents.end(); ++it) {
        const std::string& name = it->first;
        const CycleEvent* event = dynamic_cast<const CycleEvent*>(it->second);
        const CycleEvent::CType& type = event->Type();
        if (event->ThisEventDue(cycle,dcycle)) {
            events.push_back(name);
        }
    }

    for (std::map<std::string,Event*>::const_iterator it=timeEvents.begin(); it!=timeEvents.end(); ++it) {
        const std::string& name = it->first;
        const TimeEvent* event = dynamic_cast<const TimeEvent*>(it->second);
        TimeEvent::TType type = event->Type();
        Real dt_red;
        if (event->ThisEventDue(t,dt,dt_red)) {
            events.push_back(name);
            delta_t = (delta_t < 0  ?  dt_red  :  std::min(delta_t, dt_red));
        }
    }

    return std::pair<Real,Array<std::string> > (delta_t,events);
}

Array<std::string>
EventCoord::InitEvent(Real t, int cycle) const
{
  Array<std::string> events;
  for (std::map<std::string,Event*>::const_iterator it=cycleEvents.begin(); it!=cycleEvents.end(); ++it) {
    const std::string& name = it->first;
    const CycleEvent* event = dynamic_cast<const CycleEvent*>(it->second);
    const CycleEvent::CType& type = event->Type();
    if (event->ThisEventInit(cycle)) {
      events.push_back(name);
    }
  }

  for (std::map<std::string,Event*>::const_iterator it=timeEvents.begin(); it!=timeEvents.end(); ++it) {
    const std::string& name = it->first;
    const TimeEvent* event = dynamic_cast<const TimeEvent*>(it->second);
    TimeEvent::TType type = event->Type();
    if (event->ThisEventInit(t)) {
      events.push_back(name);
    }
  }

  return events;
}

std::ostream& operator<< (std::ostream& os, const CycleEvent& rhs)
{
    os << "CycleEvent::";
    if (rhs.type==CycleEvent::CYCLES) {
        os << "(cycles): ";
        if (rhs.cycles.size()==0) {
            os << "(none)";
        } else {
            for (int i=0; i<rhs.cycles.size(); ++i) {
                os << rhs.cycles[i] << " ";
            }
        }
    } else if (rhs.type==CycleEvent::SPS_CYCLES) {
        os << "(start,period,stop): " << rhs.start << " " << rhs.period << " " << rhs.stop;
    } else {
        os << "<INVALID>";
    }
    return os;
}

std::ostream& operator<< (std::ostream& os, const TimeEvent& rhs)
{
    os << "TimeEvent::";
    if (rhs.type==TimeEvent::TIMES) {
        os << "(times): ";
        if (rhs.times.size()==0) {
            os << "(none)";
        } else {
            for (int i=0; i<rhs.times.size(); ++i) {
                os << rhs.times[i] << " ";
            }
        }
    } else if (rhs.type==TimeEvent::SPS_TIMES) {
        os << "(start,period,stop): " << rhs.start << " " << rhs.period << " " << rhs.stop;
    } else {
        os << "<INVALID>";
    }
    return os;
}

std::ostream& operator<< (std::ostream& os, const EventCoord& rhs)
{
    os << "EventCoord registered events:" << '\n';
    const std::map<std::string,Event*>& cycleEvents = rhs.CycleEvents();
    for (std::map<std::string,Event*>::const_iterator it=cycleEvents.begin(); it!=cycleEvents.end(); ++it) {
        os << "  " << it->first << " "<< dynamic_cast<const CycleEvent&>(*it->second) << std::endl;
    }
    const std::map<std::string,Event*>& timeEvents = rhs.TimeEvents();
    for (std::map<std::string,Event*>::const_iterator it=timeEvents.begin(); it!=timeEvents.end(); ++it) {
        os << "  " << it->first << " "<< dynamic_cast<const TimeEvent&>(*it->second) << std::endl;
    }
    return os;
}
