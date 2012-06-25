
#ifndef _PMAmr_H_
#define _PMAmr_H_ 

#include <Amr.H>
#include <EventCoord.H>

class Layout; 
class MFTGrowFill; 

#include <Layout.H>

class PMAmr
    : public Amr
{
public:
    PMAmr();
    virtual ~PMAmr();
    virtual void coarseTimeStep (Real stop_time);
    virtual void init (Real strt_time,
                       Real stop_time);

    virtual void initialInit (Real strt_time,
                              Real stop_time);

    static EventCoord& eventCoord() {return event_coord;}

    Real StartTime() const {return strt_time;}

    static Layout& GetLayout() {return layout;}
protected:
    static Layout layout;

protected:
    static EventCoord event_coord;
    static void CleanupStatics ();

private:
    //
    // Disallowed
    //
    PMAmr (const PMAmr& rhs);
    PMAmr& operator= (const PMAmr& rhs);

    Real strt_time;
};

#endif /*_PMAmr_H_*/
