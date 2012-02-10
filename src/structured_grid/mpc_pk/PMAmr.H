
#ifndef _PMAmr_H_
#define _PMAmr_H_ 

#include <Amr.H>

//
//
// Manages hierarchy of levels for time-dependent AMR computations.
//
// Subclass from Amr to enable specific control of output

struct EventCoord
{
    struct CycleEvent
    {
        enum CycleType {SPS, CYCLES};
        CycleEvent() {}
        CycleEvent(int start, int period, int stop) 
            : start(start), period(period), stop(stop), type(SPS) {}
        CycleEvent(const Array<int>& cycles) 
            : cycles(cycles), type(CYCLES) {}
        CycleType Type() const {return type;}
        bool ThisEventDue(int cycle) const;

    protected:
        int start, period, stop;
        Array<int> cycles;
        CycleType type;
    };

    struct TimeEvent
    {
        enum TimeType {SPS, TIMES};
        TimeEvent() {}
        TimeEvent(Real start, Real period, Real stop) 
            : start(start), period(period), stop(stop), type(SPS) {}
        TimeEvent(const Array<Real>& times) 
            : times(times), type(TIMES) {}
        TimeType Type() const {return type;}
        bool ThisEventDue(Real time, Real dt, Real& dt_red) const;

    protected:
        Real start, period, stop;
        Array<Real> times;
        TimeType type;
    };

    std::pair<Real,Array<std::string> > NextEvent(Real time,
                                                  Real dt,
                                                  int  cycle) const;
    void InsertCycleEvent(const std::string& label, int start, int period, int stop);
    void InsertCycleEvent(const std::string& label, const Array<int>& cycles);

    void InsertTimeEvent(const std::string& label, Real start, Real period, Real stop);
    void InsertTimeEvent(const std::string& label, const Array<Real>& cycles);

protected:
    std::map<std::string,CycleEvent> cycleEvents;
    std::map<std::string,TimeEvent> timeEvents;
};



class PMAmr
    : public Amr
{
public:
    PMAmr() : Amr() {}
    virtual void coarseTimeStep (Real stop_time);
    bool time_to_write_checkpoint(Real cumtime,Real dt_local,int steps);
    static EventCoord& eventCoord() {return event_coord;}

protected:
    static EventCoord event_coord;

private:
    //
    // Disallowed
    //
    PMAmr (const PMAmr& rhs);
    PMAmr& operator= (const PMAmr& rhs);
};

#endif /*_PMAmr_H_*/
