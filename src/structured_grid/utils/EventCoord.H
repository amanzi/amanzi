/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _EventCoord_H_
#define _EventCoord_H_

#include <string>
#include <map>
#include <cmath>

#include <Array.H>
#include <REAL.H>

#include <iostream>
struct EventCoord
{
    struct Event
    {
        Event() {}
        Event(const Event& rhs) {}
        virtual ~Event() {}
        virtual bool IsCycle() const {return false;}
        virtual bool IsTime() const {return false;}
        virtual bool operator<(const Event& rhs) const
            {return this < &rhs;}
        virtual Event* Clone() const {return new Event(*this);}
    };

    struct CycleEvent
        : public Event
    {
        enum CType {SPS_CYCLES, CYCLES, INVALID};
        CycleEvent() {}
        virtual ~CycleEvent() {}
        CycleEvent(int start, int period, int stop)
            : start(start), period(period), stop(stop), type(SPS_CYCLES) {}
        CycleEvent(const Array<int>& cycles);
        CycleEvent(const CycleEvent& rhs);
        const CType& Type() const {return type;}
        virtual bool IsCycle() const {return true;}
        virtual bool IsTime() const {return false;}
        bool ThisEventDue(int cycle, int dcycle) const;
        bool ThisEventInit(int cycle) const;
        virtual CycleEvent* Clone() const {return new CycleEvent(*this);}

        int start, period, stop;
        Array<int> cycles;
        CType type;
    };

    struct TimeEvent
        : public Event
    {
        enum TType {SPS_TIMES, TIMES, INVALID};
        TimeEvent() {}
        virtual ~TimeEvent() {}
        TimeEvent(Real start, Real period, Real stop)
            : start(start), period(period), stop(stop), type(SPS_TIMES) {}
        TimeEvent(const Array<Real>& times);
        TimeEvent(const TimeEvent& rhs);
        const TType& Type() const {return type;}
        virtual bool IsCycle() const {return false;}
        virtual bool IsTime() const {return true;}
        bool ThisEventDue(Real time, Real dt, Real& dt_red) const;
        bool ThisEventInit(Real time) const;
        virtual TimeEvent* Clone() const {return new TimeEvent(*this);}

        Real start, period, stop;
        Array<Real> times;
        TType type;
    };

    ~EventCoord();
    std::pair<Real,Array<std::string> > NextEvent(Real time, Real dt, int  cycle, int dcycle) const;
    Array<std::string> InitEvent(Real time, int  cycle) const;

    void Register(const std::string& label, const Event*  event);
    bool EventRegistered(const std::string& label) const;
    const Event& operator[] (const std::string& label) const;
    const std::map<std::string,Event*>& CycleEvents() const;
    const std::map<std::string,Event*>& TimeEvents() const;

protected:
    std::map<std::string,Event*> cycleEvents;
    std::map<std::string,Event*> timeEvents;
};

std::ostream& operator<< (std::ostream& os, const EventCoord::CycleEvent& rhs);
std::ostream& operator<< (std::ostream& os, const EventCoord::TimeEvent& rhs);
std::ostream& operator<< (std::ostream& os, const EventCoord& rhs);

#endif
