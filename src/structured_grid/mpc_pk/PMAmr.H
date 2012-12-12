
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

    virtual void writePlotFile ();
    virtual void checkPoint ();

    void pm_timeStep (int  level,
                      Real time,
                      int  iteration,
                      int  niter,
                      Real stop_time);
        

    static EventCoord& eventCoord() {return event_coord;}

    static Layout& GetLayout() {return layout;}

    Real Dt0FromPreviousAdvance() const {return dt0_from_previous_advance;}
    Real Dt0BeforeEventCut() const {return dt0_before_event_cut;}

protected:
    static Layout layout;
    static EventCoord event_coord;
    static void CleanupStatics ();
    static int plot_file_digits;
    static int chk_file_digits;

private:
    //
    // Disallowed
    //
    PMAmr (const PMAmr& rhs);
    PMAmr& operator= (const PMAmr& rhs);

    Real dt0_before_event_cut;
    Real dt0_from_previous_advance;
};

#endif /*_PMAmr_H_*/
