
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

    void pm_timeStep (int  level,
                      Real time,
                      int  iteration,
                      int  niter,
                      Real stop_time);
        

    static EventCoord& eventCoord() {return event_coord;}

    void SetCumTime(Real strt_time) {cumtime = strt_time;}
    Real GetCumTime() const {return cumtime;}
  
    Real StartTime() const {return strt_time;}
    void SetStartTime(Real time) {strt_time = time;}

    static Layout& GetLayout() {return layout;}

    Real Dt0FromPreviousAdvance() const {return dt0_from_previous_advance;}
    Real Dt0BeforeEventCut() const {return dt0_before_event_cut;}

protected:
    static Layout layout;
    static EventCoord event_coord;
    static void CleanupStatics ();

private:
    //
    // Disallowed
    //
    PMAmr (const PMAmr& rhs);
    PMAmr& operator= (const PMAmr& rhs);

    Real strt_time;

    Real dt0_before_event_cut;
    Real dt0_from_previous_advance;
};

#endif /*_PMAmr_H_*/
