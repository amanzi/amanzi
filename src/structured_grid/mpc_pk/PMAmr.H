
#ifndef _PMAmr_H_
#define _PMAmr_H_ 

#include <Amr.H>

class Layout; 
class MFTGrowFill; 

#include <Layout.H>
#include <MatIDFiller.H>
#include <Observation.H>
#include <EventCoord.H>
#include <MatFiller.H>

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
                      int  niter);
        

    EventCoord& eventCoord() {return event_coord;}
    Real process_events(bool& write_plotfile_after_step,
                        bool& write_checkpoint_after_step,
                        Array<int>& observations_after_step,
                        bool& begin_tpc,
                        EventCoord& event_coord,
                        Real time, Real dt, int iter, int diter);
    void initial_events(bool& write_plotfile_now,
                        bool& write_checkpoint_now,
                        Array<int>& observations_now,
                        bool& begin_tpc_now,
                        EventCoord& event_coord,
                        Real time, int iter);

    Layout& GetLayout() {return layout;}
    MatFiller* GetMatFiller() {return materialFiller;}
    const Array<Real>& TPCStartTimes() const {return tpc_start_times;}

    Real Dt0FromPreviousAdvance() const {return dt0_from_previous_advance;}
    Real Dt0BeforeEventCut() const {return dt0_before_event_cut;}
    PArray<Observation>& TheObservations() {return observations;}
    Real StopTime() const {return stop_time;}
    int MaxStep() const {return max_step;}
    void RegisterEvent(const std::string& event_label,
                       EventCoord::Event* event); // takes control of pointer

protected:
    void InitializeControlEvents();

    Layout layout;
    PArray<Observation> observations;
    std::string obs_outputfile;

    EventCoord event_coord;
    std::map<std::string,EventCoord::Event*> defined_events; // accumulate all defined, register as needed
    std::string observation_output_file;

    Array<std::string> vis_cycle_macros, vis_time_macros;
    Array<std::string> chk_cycle_macros, chk_time_macros;

    Array<std::string> tpc_labels;
    Array<Real> tpc_start_times;
    Array<Real> tpc_initial_time_steps;
    Array<Real> tpc_initial_time_step_multipliers;
    Array<Real> tpc_maximum_time_steps;

    static void CleanupStatics ();
    int plot_file_digits;
    int chk_file_digits;
    int  max_step;
    Real stop_time;


    void SetUpMaterialServer();
    MatFiller* materialFiller;

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
