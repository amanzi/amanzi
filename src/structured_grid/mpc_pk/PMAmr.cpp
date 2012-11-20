#include <PMAmr.H>
#include <PorousMedia.H>
#include <Observation.H>

EventCoord PMAmr::event_coord;
Layout PMAmr::layout;

namespace
{
    //
    // These are all ParmParse'd in.  Set defaults in Initialize block below!!!
    //
    bool pmamr_initialized = false;
    int  regrid_on_restart;
    int  use_efficient_regrid;
    int  compute_new_dt_on_regrid;
    int  plotfile_on_restart;
}


void
PMAmr::CleanupStatics ()
{
    pmamr_initialized = false;
    layout.Clear();
}


static bool initialized = false;
PMAmr::PMAmr()
    : Amr()
{
    if (!pmamr_initialized) {
        regrid_on_restart        = 0;
        use_efficient_regrid     = 0;
        plotfile_on_restart      = 0;
        compute_new_dt_on_regrid = 0;

        BoxLib::ExecOnFinalize(PMAmr::CleanupStatics);
        pmamr_initialized = true;
    }
    layout.SetParent(this);

    dt0_before_event_cut = -1;
    dt0_from_previous_advance = -1;
}

PMAmr::~PMAmr()
{}

static Real
process_events(bool& write_plotfile_after_step,
               bool& write_checkpoint_after_step,
               Array<int>& observations_after_step,
               bool& begin_tpc,
               EventCoord& event_coord,
               Real time, Real dt, int iter, int diter)
{
    write_plotfile_after_step = false;
    write_checkpoint_after_step = false;
    begin_tpc = false;
    Array<std::string>& vis_cycle_macros = PorousMedia::vis_cycle_macros;
    Array<std::string>& vis_time_macros = PorousMedia::vis_time_macros;
    Array<std::string>& chk_cycle_macros = PorousMedia::chk_cycle_macros;
    Array<std::string>& tpc_labels = PorousMedia::tpc_labels;

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

            for (int k=0; k<vis_time_macros.size(); ++k) {
                if (eventList[j] == vis_time_macros[k]) {
                    write_plotfile_after_step = true;
                }
            }

            for (int k=0; k<chk_cycle_macros.size(); ++k) {
                if (eventList[j] == chk_cycle_macros[k]) {
                    write_checkpoint_after_step = true;
                }
            }

            for (int k=0; k<tpc_labels.size(); ++k) {
                if (eventList[j] == tpc_labels[k]) {
                    begin_tpc = true;
                }
            }
        }
    }
    return dt_new;
}


void
PMAmr::init (Real t_start,
             Real t_stop)
{
    if (!restart_file.empty() && restart_file != "init")
    {
        restart(restart_file);
        setStartTime(0); // FIXME: This needs to be written to the checkpoint
    }
    else
    {
        initialInit(t_start,t_stop);

        bool write_plot, write_check, begin_tpc;
        Array<int> initial_observations;

        process_events(write_plot,write_check,initial_observations,begin_tpc,event_coord,
                       cumtime, dt_level[0], level_steps[0], 0);

        if (write_plot) {
            writePlotFile();
        }

        if (write_check) {
            checkPoint();
        }
    }
}

void
PMAmr::pm_timeStep (int  level,
                    Real time,
                    int  iteration,
                    int  niter,
                    Real stop_time)
{
    //
    // Allow regridding of level 0 calculation on restart.
    //
    if (finest_level == 0 && Amr::RegridOnRestart())
    {
        regrid_on_restart = 0;
        //
        // Coarsening before we split the grids ensures that each resulting
        // grid will have an even number of cells in each direction.
        //
        BoxArray lev0(BoxLib::coarsen(geom[0].Domain(),2));
        //
        // Now split up into list of grids within max_grid_size[0] limit.
        //
        lev0.maxSize(max_grid_size[0]/2);
        //
        // Now refine these boxes back to level 0.
        //
        lev0.refine(2);

        //
        // If use_efficient_regrid flag is set, then test to see whether we in fact 
        //    have just changed the level 0 grids. If not, then don't do anything more here.
        //
        if ( !( (use_efficient_regrid == 1) && (lev0 == amr_level[0].boxArray()) ) ) 
        {
            //
            // Construct skeleton of new level.
            //
            AmrLevel* a = (*levelbld)(*this,0,geom[0],lev0,cumtime);

            a->init(amr_level[0]);
            amr_level.clear(0);
            amr_level.set(0,a);

            amr_level[0].post_regrid(0,0);

            if (ParallelDescriptor::IOProcessor())
            {
               if (verbose > 1)
               {
                  printGridInfo(std::cout,0,finest_level);
               }
               else if (verbose > 0)
               {
                  printGridSummary(std::cout,0,finest_level);
               }
            }

            if (record_grid_info && ParallelDescriptor::IOProcessor())
                printGridInfo(gridlog,0,finest_level);
        }
        else
        {
            if (verbose > 0 && ParallelDescriptor::IOProcessor())
                std::cout << "Regridding at level 0 but grids unchanged " << std::endl;
        }
    }
    else
    {
        int lev_top = std::min(finest_level, max_level-1);

        for (int i = level; i <= lev_top; i++)
        {
            const int old_finest = finest_level;

            if (level_count[i] >= regrid_int[i] && amr_level[i].okToRegrid())
            {
                regrid(i,time);

                //
                // Compute new dt after regrid if at level 0 and compute_new_dt_on_regrid.
                //
                if ( compute_new_dt_on_regrid && (i == 0) )
                {
                    int post_regrid_flag = 1;
                    amr_level[0].computeNewDt(finest_level,
                                              sub_cycle,
                                              n_cycle,
                                              ref_ratio,
                                              dt_min,
                                              dt_level,
                                              stop_time, 
                                              post_regrid_flag);
                }

                for (int k = i; k <= finest_level; k++)
                    level_count[k] = 0;

                if (old_finest < finest_level)
                {
                    //
                    // The new levels will not have valid time steps
                    // and iteration counts.
                    //
                    for (int k = old_finest+1; k <= finest_level; k++)
                    {
                        const int fact = sub_cycle ? MaxRefRatio(k-1) : 1;
                        dt_level[k]    = dt_level[k-1]/Real(fact);
                        n_cycle[k]     = fact;
                    }
                }
            }
            if (old_finest > finest_level)
                lev_top = std::min(finest_level, max_level-1);
        }
    }
    //
    // Check to see if should write plotfile.
    // This routine is here so it is done after the restart regrid.
    //
    if (plotfile_on_restart && !(restart_file.empty()) )
    {
	plotfile_on_restart = 0;
	writePlotFile();
    }
    //
    // Advance grids at this level.
    //
    Real dt_taken, dt_suggest;
    bool step_ok = dynamic_cast<PorousMedia&>(amr_level[level]).ml_step_driver(time,dt_level[level],dt_taken,dt_suggest);

    if (step_ok) {
        dt_level[level] = dt_taken;
        int fac = 1;
        for (int lev = level+1; lev<=finest_level; ++lev) {
            fac *= n_cycle[lev];
            dt_level[lev] = dt_level[lev-1] * fac;
        }
        dt0_from_previous_advance = dt_suggest;
    }
    else {
        BoxLib::Abort("Step failed");
    }

    dt_min[level] = iteration == 1 ? dt_suggest : std::min(dt_min[level],dt_suggest);

    level_steps[level]++;
    level_count[level]++;

    if (verbose > 2 && ParallelDescriptor::IOProcessor())
    {
        std::cout << "Advanced "
                  << amr_level[level].countCells()
                  << " cells at level "
                  << level
		  << "    dt0_from_previous_advance is now " << dt0_from_previous_advance
                  << std::endl;
    }

#ifdef USE_STATIONDATA
    station.report(time+dt_level[level],level,amr_level[level]);
#endif

#ifdef USE_SLABSTAT
    AmrLevel::get_slabstat_lst().update(amr_level[level],time,dt_level[level]);
#endif
    //
    // Advance grids at higher level.
    //
    if (level < finest_level)
    {
        const int lev_fine = level+1;

        if (sub_cycle)
        {
            const int ncycle = n_cycle[lev_fine];

            for (int i = 1; i <= ncycle; i++)
                pm_timeStep(lev_fine,time+(i-1)*dt_level[lev_fine],i,ncycle,stop_time);
        }
    }

    amr_level[level].post_timestep(iteration);
}

void
PMAmr::coarseTimeStep (Real stop_time)
{
    const Real run_strt = ParallelDescriptor::second() ;    
    
    int post_regrid_flag = 0;
    amr_level[0].computeNewDt(finest_level,
                              sub_cycle,
                              n_cycle,
                              ref_ratio,
                              dt_min,
                              dt_level,
                              stop_time,
                              post_regrid_flag);

    Real strt_time = startTime();

    if (cumtime<strt_time+.001*dt_level[0]  && verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << "\nBEGIN TRANSIENT INTEGRATION, TIME = " << cumtime << '\n' << std::endl;
    }

    bool write_plot, write_check, begin_tpc;
    Array<int> observations_to_process;

    // NOTE: This is a hack to help keep dt from jumping too much.  Normally, the event processing
    // will only detect an even during the current time step.  If one is detected, the time step
    // is cut to hit the event.  Here, we check the current time interval, but also a guess for
    // the next one.  If one is detected during the next interval, we split the time remaining
    // between the two steps.  Under normal conditions, this will ensure that dt doesn't abruptly
    // change by more than 50%, at least not due to an event detection.
    bool look_ahead_two_steps = true;
    Real dt2_red = -1;
    if (look_ahead_two_steps) {
        dt2_red = process_events(write_plot,write_check,observations_to_process,begin_tpc,event_coord,
                                 cumtime, 2*dt_level[0], level_steps[0] + 1, 1);
        observations_to_process.clear();
        write_check=false;
        write_plot=false;
        begin_tpc=false;
    }

    Real dt_red = process_events(write_plot,write_check,observations_to_process,begin_tpc,event_coord,
                                 cumtime, dt_level[0], level_steps[0], 1);
    
    // Note: if dt_red > 0, then dt_red == dt2_red
    if (dt2_red > 0  &&  dt_red < 0) {
            dt_red = dt2_red / 2; // Nothing in 1 step, but something in 2
    }

    if (dt_red > 0  &&  dt_red < dt_level[0]) {
        
        if (begin_tpc) {
            dt0_before_event_cut = -1; // "forget" current time step, we're headed into a Time Period Control interval
        }
        else if (dt0_before_event_cut < 0) {
            dt0_before_event_cut = dt_level[0];
        }

        Array<Real> dt_new(finest_level+1,dt_red);
        for (int lev = 1; lev <= finest_level; lev++) {
            dt_new[lev] = dt_new[lev-1]/Real(MaxRefRatio(lev-1));
        }
        setDtLevel(dt_new);
    }

    // Do time step
    pm_timeStep(0,cumtime,1,1,stop_time);

    cumtime += dt_level[0];

    amr_level[0].postCoarseTimeStep(cumtime);
    
    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_stop = ParallelDescriptor::second() - run_strt;

        ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

        if (verbose>3 && ParallelDescriptor::IOProcessor())
            std::cout << "\nCoarse TimeStep time: " << run_stop << '\n' ;

        long min_fab_bytes = BoxLib::total_bytes_allocated_in_fabs_hwm;
        long max_fab_bytes = BoxLib::total_bytes_allocated_in_fabs_hwm;

        ParallelDescriptor::ReduceLongMin(min_fab_bytes,IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_bytes,IOProc);
        //
        // Reset to zero to calculate high-water-mark for next timestep.
        //
        BoxLib::total_bytes_allocated_in_fabs_hwm = 0;

        if (verbose>3 && ParallelDescriptor::IOProcessor())
            std::cout << "\nFAB byte spread across MPI nodes for timestep: ["
                      << min_fab_bytes << " ... " << max_fab_bytes << "]\n";
    }

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << "\nSTEP = "
                 << level_steps[0]
                  << " COMPLETE.  TIME = "
                  << cumtime
                  << " DT = "
                  << dt_level[0]
                  << '\n';
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
        checkPoint();
    }

    if (write_plot)
    {
        writePlotFile();
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
PMAmr::initialInit (Real t_start,
                    Real t_stop)
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
    setStartTime(t_start);
    setCumTime(startTime());
    //
    // Define base level grids.
    //
    defBaseLevel(startTime());
    //
    // Compute dt and set time levels of all grid data.
    //
    amr_level[0].computeInitialDt(finest_level,
                                  sub_cycle,
                                  n_cycle,
                                  ref_ratio,
                                  dt_level,
                                  t_stop);
    //
    // The following was added for multifluid.
    //
    Real dt0   = dt_level[0];
    dt_min[0]  = dt_level[0];
    n_cycle[0] = 1;


    //
    // The following was added to avoid stepping over registered events
    //
    bool write_plot, write_check, begin_tpc;
    Array<int> observations_to_process;
    
    Real dt_red = process_events(write_plot,write_check,observations_to_process,begin_tpc,event_coord,
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
        bldFineLevels(startTime());

    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev].setTimeLevel(startTime(),dt_level[lev],dt_level[lev]);

    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev].post_regrid(0,finest_level);
    //
    // Perform any special post_initialization operations.
    //
    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev].post_init(t_stop);

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

