
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
    if (type == SPS_CYCLES)
    {
        if ( ( ( stop < 0 ) || ( cycle + dCycle < stop ) )
             && ( cycle + dCycle >= start ) )
        {
            return ( (period==1) || ((cycle-start+dCycle)%period == 0) );
        }
    }
    else {
        for (int i=0; i<cycles.size(); ++i) {
            if (cycles[i] == cycle+1) {
                return true;
            }
        }
    }
    return false;
}

bool
EventCoord::TimeEvent::ThisEventDue(Real time, Real dt, Real& dt_red) const
{
    dt_red = dt;
    if (type == SPS_TIMES)
    {
        if ( ( ( stop < 0 ) || ( time < stop ) )
             && ( time >= start ) )
        {
            Real intPart;
            if ( modf( (time-start+dt)/period, &intPart ) == 0 )
            {
                return true;
            }
        }
    }
    else {
        for (int i=0; i<times.size(); ++i) {
            Real eps = 1.e-6*dt;
            if (times[i] > time  &&  times[i] <= time + dt + eps) {
                dt_red = std::min(dt, times[i]-time);
                return true;
            }
            if (times[i] == time) {
                Real max_dt = 0;
                for (int j=0; j<times.size(); ++j) {
                    if (time < times[j]) {
                        max_dt = std::max(max_dt,times[j]-time);
                    }
                }
                if (max_dt == 0) {
                    return false;
                }
                else {
                    Real min_dt = max_dt;
                    for (int j=0; j<times.size(); ++j) {
                        if (i!=j  && time < times[j]  &&  time + dt > times[j]) {
                            min_dt = std::min(min_dt, times[j] - time);
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
EventCoord::NextEvent(Real time, Real dt, int cycle, int dcycle) const
{
    Array<std::string> events;
    Real delta_time = -1;
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
        if (event->ThisEventDue(time,dt,dt_red)) {
            events.push_back(name);
            delta_time = (delta_time < 0  ?  dt_red  :  std::min(delta_time, dt_red));
        }
    }

    return std::pair<Real,Array<std::string> > (delta_time,events);
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
