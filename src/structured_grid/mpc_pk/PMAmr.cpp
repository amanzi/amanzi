#include <PMAmr.H>
#include <PorousMedia.H>
#include <Observation.H>

EventCoord PMAmr::event_coord;
Layout PMAmr::layout;
Real PMAmr::dt_previous = -1;

void
PMAmr::CleanupStatics ()
{
    layout.Clear();
    dt_previous = -1;
}


static bool initialized = false;
PMAmr::PMAmr()
    : Amr()
{
    if (!initialized) {
        BoxLib::ExecOnFinalize(PMAmr::CleanupStatics);
        initialized = true;
    }
    layout.SetParent(this);
}

PMAmr::~PMAmr()
{}

Real
process_events(bool& write_plotfile_after_step,
               bool& write_checkpoint_after_step,
               Array<int>& observations_after_step,
               EventCoord& event_coord,
               Real time, Real dt, int iter, int diter)
{
    write_plotfile_after_step = false;
    write_checkpoint_after_step = false;
    Array<std::string>& vis_cycle_macros = PorousMedia::vis_cycle_macros;
    Array<std::string>& chk_cycle_macros = PorousMedia::chk_cycle_macros;

    Real dt_new = dt;
    std::pair<Real,Array<std::string> > nextEvent = event_coord.NextEvent(time,dt,iter, diter);
    PArray<Observation>& observations = PorousMedia::TheObservationArray();

    if (nextEvent.second.size()) 
    {
        // Process event
        dt_new = nextEvent.first;
        const Array<std::string>& eventList = nextEvent.second;

        for (int j=0; j<eventList.size(); ++j) {

            for (int i=0; i<observations.size(); ++i) {
                Observation& observation = observations[i];
                const std::string& event_label = observation.event_label;
                
                if (eventList[j] == event_label) {
                    observations_after_step.push_back(i);
                }
            }

            for (int k=0; k<vis_cycle_macros.size(); ++k) {
                if (eventList[j] == vis_cycle_macros[k]) {
                    write_plotfile_after_step = true;
                }
            }

            for (int k=0; k<chk_cycle_macros.size(); ++k) {
                if (eventList[j] == chk_cycle_macros[k]) {
                    write_checkpoint_after_step = true;
                }
            }
        }
    }
    return dt_new;
}


void
PMAmr::init (Real strt_time_,
             Real stop_time)
{
    strt_time = strt_time_;

    if (!restart_file.empty() && restart_file != "init")
    {
        restart(restart_file);
    }
    else
    {
        initialInit(strt_time,stop_time);

        bool write_plot, write_check;
        Array<int> initial_observations;
        process_events(write_plot,write_check,initial_observations,event_coord,
                       cumtime, dt_level[0], level_steps[0], 0);

        if (write_plot) {

            writePlotFile(plot_file_root,level_steps[0]);
        }

        if (write_check) {
            checkPoint();
        }
    }
}

void
PMAmr::coarseTimeStep (Real stop_time)
{
    const Real run_strt = ParallelDescriptor::second() ;    
    
    //
    // Compute new dt.
    //
    if (level_steps[0] <= 0 && PorousMedia::DtInit() > 0) 
    {
        int n_factor = 1;
        for (int i = 0; i <= max_level; i++)
        {
            n_factor   *= n_cycle[i];
            dt_level[i] = PorousMedia::DtInit()/( (Real)n_factor );
        }
    }
    else 
    {
        int post_regrid_flag = 0;
        amr_level[0].computeNewDt(finest_level,
                                  sub_cycle,
                                  n_cycle,
                                  ref_ratio,
                                  dt_min,
                                  dt_level,
                                  stop_time,
                                  post_regrid_flag);
        if (dt_previous > 0) {
            dt_level[0] = dt_previous;
            for (int lev = 1; lev <= finest_level; lev++) {
                dt_level[lev] = dt_level[lev-1]/Real(MaxRefRatio(lev-1));
            }
        }
    }


    if (cumtime<strt_time+.001*dt_level[0]  && verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << "\nBEGIN TRANSIENT INTEGRATION, TIME = " << cumtime << '\n' << std::endl;
    }

    bool write_plot, write_check;
    Array<int> observations_to_process;
    Real dt_red = process_events(write_plot,write_check,observations_to_process,event_coord,
                                 cumtime, dt_level[0], level_steps[0], 1);

    if (dt_red > 0  &&  dt_red < dt_level[0]) {
        dt_previous = dt_level[0];
        Array<Real> dt_new(finest_level+1,dt_red);
        for (int lev = 1; lev <= finest_level; lev++) {
            dt_new[lev] = dt_new[lev-1]/Real(MaxRefRatio(lev-1));
        }
        setDtLevel(dt_new);
    }

    // Do time step
    timeStep(0,cumtime,1,1,stop_time);

    cumtime += dt_level[0];

    amr_level[0].postCoarseTimeStep(cumtime);
    
    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_stop = ParallelDescriptor::second() - run_strt;

        ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

        if (verbose>1 && ParallelDescriptor::IOProcessor())
            std::cout << "\nCoarse TimeStep time: " << run_stop << '\n' ;

        long min_fab_bytes = BoxLib::total_bytes_allocated_in_fabs_hwm;
        long max_fab_bytes = BoxLib::total_bytes_allocated_in_fabs_hwm;

        ParallelDescriptor::ReduceLongMin(min_fab_bytes,IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_bytes,IOProc);
        //
        // Reset to zero to calculate high-water-mark for next timestep.
        //
        BoxLib::total_bytes_allocated_in_fabs_hwm = 0;

        if (verbose>1 && ParallelDescriptor::IOProcessor())
            std::cout << "\nFAB byte spread across MPI nodes for timestep: ["
                      << min_fab_bytes << " ... " << max_fab_bytes << "]\n";
    }

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << "\nSTEP = "
                 << level_steps[0]
                  << " TIME = "
                  << cumtime
                  << " DT = "
                  << dt_level[0]
                  << '\n'
                  << std::endl;
    }
    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
        runlog << "STEP = "
               << level_steps[0]
               << " TIME = "
               << cumtime
               << " DT = "
               << dt_level[0]
               << '\n';
    }
    if (record_run_info_terse && ParallelDescriptor::IOProcessor())
        runlog_terse << level_steps[0] << " " << cumtime << " " << dt_level[0] << '\n';

    
    int to_checkpoint = 0;    
    int to_stop       = 0;    
    if (ParallelDescriptor::IOProcessor())
    {
        FILE *fp;
        if ((fp=fopen("dump_and_continue","r")) != 0)
        {
            remove("dump_and_continue");
            to_checkpoint = 1;
            fclose(fp);
        }
        else if ((fp=fopen("stop_run","r")) != 0)
        {
            remove("stop_run");
            to_stop = 1;
            fclose(fp);
        }
        else if ((fp=fopen("dump_and_stop","r")) != 0)
        {
            remove("dump_and_stop");
            to_checkpoint = 1;
            to_stop = 1;
            fclose(fp);
        }
    }
    ParallelDescriptor::Bcast(&to_checkpoint, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&to_stop,       1, ParallelDescriptor::IOProcessorNumber());

    if (write_check || (to_checkpoint==1))
    {
        last_checkpoint = level_steps[0];
        checkPoint();
    }

    if (write_plot)
    {
        last_plotfile = level_steps[0];
        writePlotFile(plot_file_root,level_steps[0]);
    }

    if (to_stop)
    {
        ParallelDescriptor::Barrier();
        if (to_checkpoint == 1)
        {
            BoxLib::Abort("Aborted by user w/checkpoint");
        }
        else {
            BoxLib::Abort("Aborted by user w/o checkpoint");
        }
    }
}

void
PMAmr::initialInit (Real strt_time,
                    Real stop_time)
{
    checkInput();
    //
    // Generate internal values from user-supplied values.
    //
    finest_level = 0;
    //
    // Init problem dependent data.
    //
    int init = true;

    readProbinFile(init);

#ifdef BL_SYNC_RANTABLES
    int iGet(0), iSet(1);
    const int iTableSize(64);
    Real *RanAmpl = new Real[iTableSize];
    Real *RanPhase = new Real[iTableSize];
    FORT_SYNC_RANTABLES(RanPhase, RanAmpl, &iGet);
    ParallelDescriptor::Bcast(RanPhase, iTableSize);
    ParallelDescriptor::Bcast(RanAmpl, iTableSize);
    FORT_SYNC_RANTABLES(RanPhase, RanAmpl, &iSet);
    delete [] RanAmpl;
    delete [] RanPhase;
#endif
    cumtime = strt_time;
    //
    // Define base level grids.
    //
    defBaseLevel(strt_time);
    //
    // Compute dt and set time levels of all grid data.
    //
    amr_level[0].computeInitialDt(finest_level,
                                  sub_cycle,
                                  n_cycle,
                                  ref_ratio,
                                  dt_level,
                                  stop_time);
    //
    // The following was added for multifluid.
    //
    Real dt0   = dt_level[0];
    dt_min[0]  = dt_level[0];
    n_cycle[0] = 1;


    //
    // The following was added to avoid stepping over registered events
    //
    bool write_plot, write_check;
    Array<int> observations_to_process;
    Real dt_red = process_events(write_plot,write_check,observations_to_process,event_coord,
                                 cumtime, dt_level[0], level_steps[0], 1);
    if (dt_red > 0  &&  dt_red < dt0) {
        dt_min[0]  = dt_level[0];
    }


    for (int lev = 1; lev <= max_level; lev++)
    {
        const int fact = sub_cycle ? ref_ratio[lev-1][0] : 1;

        dt0           /= Real(fact);
        dt_level[lev]  = dt0;
        dt_min[lev]    = dt_level[lev];
        n_cycle[lev]   = fact;
    }

    if (max_level > 0)
        bldFineLevels(strt_time);

    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev].setTimeLevel(strt_time,dt_level[lev],dt_level[lev]);

    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev].post_regrid(0,finest_level);
    //
    // Perform any special post_initialization operations.
    //
    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev].post_init(stop_time);

    for (int lev = 0; lev <= finest_level; lev++)
    {
        level_count[lev] = 0;
        level_steps[lev] = 0;
    }

    if (ParallelDescriptor::IOProcessor())
    {
       if (verbose > 1)
       {
           std::cout << "INITIAL GRIDS \n";
           printGridInfo(std::cout,0,finest_level);
       }
       else if (verbose > 0)
       { 
           std::cout << "INITIAL GRIDS \n";
           printGridSummary(std::cout,0,finest_level);
       }
    }

    if (record_grid_info && ParallelDescriptor::IOProcessor())
    {
        gridlog << "INITIAL GRIDS \n";
        printGridInfo(gridlog,0,finest_level);
    }

#ifdef USE_STATIONDATA
    station.init(amr_level, finestLevel());
    station.findGrid(amr_level,geom);
#endif
}

