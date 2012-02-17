
#include <EventCoord.H>
#include <iostream>

bool
EventCoord::CycleEvent::ThisEventDue(int cycle, int dCycle) const
{
    if (type ==SPS)
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
    if (type ==SPS)
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
            if (times[i] > time  &&  times[i] <= time + dt) {
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
                        if (time < times[j]) {
                            min_dt = std::min(min_dt, times[j] - time);
                        }
                    }
                    dt_red = std::min(dt, min_dt);
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
    for (std::map<std::string,CycleEvent>::const_iterator it=cycleEvents.begin(); it!=cycleEvents.end(); ++it) {
        const std::string& name = it->first;
        const CycleEvent& event = it->second;
        CycleEvent::CycleType type = event.Type();
        if (event.ThisEventDue(cycle,dcycle)) {
            events.push_back(name);
        }
    }

    for (std::map<std::string,TimeEvent>::const_iterator it=timeEvents.begin(); it!=timeEvents.end(); ++it) {
        const std::string& name = it->first;
        const TimeEvent& event = it->second;
        TimeEvent::TimeType type = event.Type();
        Real dt_red;
        if (event.ThisEventDue(time,dt,dt_red)) {
            events.push_back(name);
            delta_time = (delta_time < 0  ?  dt_red  :  std::min(delta_time, dt_red));
        }
    }

    return std::pair<Real,Array<std::string> > (delta_time,events);
}

void
EventCoord::InsertCycleEvent(const std::string& label, int start, int period, int stop)
{
    cycleEvents[label] = CycleEvent(start,period,stop);
}

void
EventCoord::InsertCycleEvent(const std::string& label, const Array<int>& cycles)
{
    cycleEvents[label] = CycleEvent(cycles);
}

void
EventCoord::InsertTimeEvent(const std::string& label, Real start, Real period, Real stop)
{
    timeEvents[label] = TimeEvent(start,period,stop);
}

void
EventCoord::InsertTimeEvent(const std::string& label, const Array<Real>& times)
{
    timeEvents[label] = TimeEvent(times);
}

