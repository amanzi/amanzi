
#ifndef _PMAmr_H_
#define _PMAmr_H_ 

#include <Amr.H>

class Layout; 
class MFTGrowFill; 

#include <Layout.H>
#include <MatFiller.H>
#include <Observation.H>
#include <EventCoord.H>

class PMAmr
    : public Amr
{
public:
    PMAmr();
    virtual ~PMAmr();
    virtual void coarseTimeStep (Real stop_time);
    virtual void init (Real strt_time,
                       Real stop_time);

    void initialInit (Real              strt_time,
		      Real              stop_time,
		      const BoxArray*   lev0_grids = 0,
		      const Array<int>* pmap = 0);
    //
    // First part of initialInit
    //
    void InitializeInit (Real strt_time, Real stop_time,
                         const BoxArray* lev0_grids = 0, const Array<int>* pmap = 0);
    //
    // Second part of initialInit
    //
    void FinalizeInit (Real strt_time, Real stop_time);

    virtual void writePlotFile ();
    virtual void checkPoint ();
    void LinkFinalCheckpoint (int step);

    void pm_timeStep (int  level,
                      Real time,
                      int  iteration,
                      int  niter);
        

    static EventCoord& eventCoord() {return event_coord;}
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
    const Array<Real>& TPCStartTimes() const {return tpc_start_times;}

    Real Dt0FromPreviousAdvance() const {return dt0_from_previous_advance;}
    Real Dt0BeforeEventCut() const {return dt0_before_event_cut;}
    PArray<Observation>& TheObservations() {return observations;}
    Real StopTime() const {return stop_time;}
    int MaxStep() const {return max_step;}
    static void RegisterEvent(const std::string& event_label,
                              EventCoord::Event* event); // takes control of pointer

    static std::pair<Real,std::string> convert_time_units(Real t, const std::string& units);

    void FlushObservations();
    void FlushObservations(std::ostream& out);

    static void Initialize ();
    static void Finalize ();

protected:
    void restart (const std::string& filename); // Note: this should be a virtual, but it is not
    void InitializeControlEvents();
    void UpdateVisitPlotfileList() const;

    Layout layout;
    PArray<Observation> observations;
    std::string obs_outputfile;

    static EventCoord event_coord;
    static std::map<std::string,EventCoord::Event*> defined_events; // accumulate all defined, register as needed
    std::string observation_output_file;

    Array<std::string> vis_cycle_macros, vis_time_macros;
    Array<std::string> chk_cycle_macros, chk_time_macros;

    Array<std::string> tpc_labels;
    Array<Real> tpc_start_times;
    Array<Real> tpc_initial_time_steps;
    Array<Real> tpc_initial_time_step_multipliers;
    Array<Real> tpc_maximum_time_steps;

    static int plot_file_digits;
    static int chk_file_digits;
    int  max_step;
    Real stop_time;

private:
    //
    // Disallowed
    //
    PMAmr (const PMAmr& rhs);
    PMAmr& operator= (const PMAmr& rhs);

    // Private data
    Real dt0_before_event_cut;
    Real dt0_from_previous_advance;
    static bool do_output_time_in_years;
    static bool attempt_to_recover_failed_step;

    // Array of pltfiles written in this run so far
    Array<std::string> plotfiles_written;
};

#endif /*_PMAmr_H_*/
