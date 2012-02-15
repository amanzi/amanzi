
#ifndef _PMAmr_H_
#define _PMAmr_H_ 

#include <Amr.H>
#include <EventCoord.H>

class PMAmr
    : public Amr
{
public:
    PMAmr() : Amr() {}
    virtual void coarseTimeStep (Real stop_time);
    virtual void init (Real strt_time,
                       Real stop_time);

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
