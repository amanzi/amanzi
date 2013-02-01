#include <PMAmr.H>
#include <PorousMedia.H>
#include <Observation.hh>

namespace
{
  //
  // These are all ParmParse'd in, as variables with the same name are
  // not accessible to this derived class...should fix
  // NOTE: Defaults set in Initialize block below
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
}

static bool initialized = false;
PMAmr::PMAmr()
  : Amr(), materialFiller(0)
{
    if (!pmamr_initialized) {
        regrid_on_restart        = 0;
        use_efficient_regrid     = 0;
        plotfile_on_restart      = 0;
        compute_new_dt_on_regrid = 0;
        plot_file_digits         = file_name_digits;
        chk_file_digits          = file_name_digits;

        BoxLib::ExecOnFinalize(PMAmr::CleanupStatics);
        pmamr_initialized = true;
    }

    ParmParse ppa("amr");
    ppa.query("plot_file_digits",plot_file_digits);
    ppa.query("chk_file_digits",chk_file_digits);
    ppa.query("regrid_on_restart",regrid_on_restart);
    ppa.query("use_efficient_regrid",use_efficient_regrid);
    ppa.query("plotfile_on_restart",plotfile_on_restart);
    ppa.query("compute_new_dt_on_regrid",compute_new_dt_on_regrid);

    ParmParse pp;
    pp.query("max_step",max_step);
    pp.query("stop_time",stop_time);

    std::string event_name = "Stop_Time";
    defined_events[event_name] = new EventCoord::TimeEvent(Array<Real>(1,stop_time));
    PMAmr::eventCoord().Register(event_name,defined_events[event_name]);

    layout.SetParent(this);

    dt0_before_event_cut = -1;
    dt0_from_previous_advance = -1;
}

PMAmr::~PMAmr()
{
  layout.Clear();

  for (std::map<std::string,EventCoord::Event*>::iterator it=defined_events.begin(); it!=defined_events.end(); ++it) {
    delete it->second;
  }

  delete materialFiller;
}

void 
PMAmr::writePlotFile ()
{
  int file_name_digits_tmp = file_name_digits;
  file_name_digits = plot_file_digits;
  Amr::writePlotFile();
  file_name_digits = file_name_digits_tmp;
}

void 
PMAmr::checkPoint ()
{
  int file_name_digits_tmp = file_name_digits;
  file_name_digits = chk_file_digits;
  Amr::checkPoint();
  file_name_digits = file_name_digits_tmp;
}

void
PMAmr::RegisterEvent(const std::string& event_label,
                     EventCoord::Event* event)
{
  std::map<std::string,EventCoord::Event*>::iterator it=defined_events.find(event_label);
  if (it!=defined_events.end()) {
    defined_events[event_label] = event;
    eventCoord().Register(event_label,defined_events[event_label]);
  }
  else {
    delete event;
  }
}

Real
PMAmr::process_events(bool& write_plotfile_after_step,
                      bool& write_checkpoint_after_step,
                      Array<int>& observations_after_step,
                      bool& begin_tpc,
                      EventCoord& event_coord,
                      Real time, Real dt, int iter, int diter)
{
  write_plotfile_after_step = false;
  write_checkpoint_after_step = false;
  begin_tpc = false;

  Real dt_new = dt;
  std::pair<Real,Array<std::string> > nextEvent = event_coord.NextEvent(time,dt,iter, diter);

  if (nextEvent.second.size()) {
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
PMAmr::initial_events(bool& write_plotfile_now,
                      bool& write_checkpoint_now,
                      Array<int>& observations_now,
                      bool& begin_tpc_now,
                      EventCoord& event_coord,
                      Real time, int iter)
{
    write_plotfile_now = false;
    write_checkpoint_now = false;
    begin_tpc_now = false;

    Array<std::string> eventList = event_coord.InitEvent(time,iter);

    for (int j=0; j<eventList.size(); ++j) {
      
      for (int i=0; i<observations.size(); ++i) {
        Observation& observation = observations[i];
        const std::string& event_label = observation.event_label;
        
        if (eventList[j] == event_label) {
          observations_now.push_back(i);
        }
      }
      
      for (int k=0; k<vis_cycle_macros.size(); ++k) {
        if (eventList[j] == vis_cycle_macros[k]) {
          write_plotfile_now = true;
        }
      }

      for (int k=0; k<vis_time_macros.size(); ++k) {
        if (eventList[j] == vis_time_macros[k]) {
          write_plotfile_now = true;
        }
      }

      for (int k=0; k<chk_cycle_macros.size(); ++k) {
        if (eventList[j] == chk_cycle_macros[k]) {
          write_checkpoint_now = true;
        }
      }
      
      for (int k=0; k<tpc_labels.size(); ++k) {
        if (eventList[j] == tpc_labels[k]) {
          begin_tpc_now = true;
        }
      }
    }
}


void
PMAmr::init (Real t_start,
             Real t_stop)
{
    SetUpMaterialServer();
    InitializeControlEvents();

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

        initial_events(write_plot,write_check,initial_observations,begin_tpc,event_coord,
                       cumtime, level_steps[0]);

        if (write_plot) {
            int file_name_digits_tmp = file_name_digits;
            file_name_digits = plot_file_digits;
            writePlotFile();
            file_name_digits = file_name_digits_tmp;
        }

        if (write_check) {
            int file_name_digits_tmp = file_name_digits;
            file_name_digits = chk_file_digits;
            checkPoint();
            file_name_digits = file_name_digits_tmp;
        }
    }
}

void
PMAmr::pm_timeStep (int  level,
                    Real time,
                    int  iteration,
                    int  niter)
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
              if (i==0) {
                regrid(i,time);
              }
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
    // Are we in a time period control section?  If so, fix time step
    //
    if (level == 0) {
      int tpc_interval = -1;
      for (int i=0; i<tpc_labels.size(); ++i) {
        if (time >= tpc_start_times[i]) {
          tpc_interval = i;
        }
      }
      Real dt_tpc;
      bool tpc_active = false;
      if (tpc_interval >= 0) {
        if (tpc_initial_time_steps.size()>tpc_interval
            && tpc_initial_time_steps[tpc_interval]>0
            && (std::abs(time - tpc_start_times[tpc_interval]) < tpc_initial_time_steps[tpc_interval])) {
          tpc_active = true;
          dt_tpc = tpc_initial_time_steps[tpc_interval];
        }
        if (tpc_initial_time_step_multipliers.size()>tpc_interval
            && tpc_initial_time_step_multipliers[tpc_interval] > 0) {
          tpc_active = true;
          dt_tpc *= tpc_initial_time_step_multipliers[tpc_interval];
        }
        if (tpc_maximum_time_steps.size()>tpc_interval
            && tpc_maximum_time_steps[tpc_interval] > 0) {
          tpc_active = true;
          dt_tpc = std::min(dt_tpc,tpc_maximum_time_steps[tpc_interval]);
        }
      }
      if (tpc_active && dt_level[level] > dt_tpc) {
        dt_level[level] = dt_tpc;
        for (int lev = level+1; lev<=finest_level; ++lev) {
          dt_level[lev] = dt_level[lev-1] / n_cycle[lev];
        }
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
    PorousMedia& pm = dynamic_cast<PorousMedia&>(amr_level[level]);
    bool step_ok = pm.ml_step_driver(time,iteration,niter,dt_level[level],dt_taken,dt_suggest);

    if (step_ok) {
      if (level==0) {
        dt_level[level] = dt_taken;
        for (int lev = level+1; lev<=finest_level; ++lev) {
          dt_level[lev] = dt_level[lev-1] / n_cycle[lev];
        }
        dt0_from_previous_advance = dt_suggest;
      }
    }
    else {
        BoxLib::Abort("Step failed");
    }

    dt_min[level] = iteration == 1 ? dt_suggest : std::min(dt_min[level],dt_suggest);
    level_steps[level]++;
    level_count[level]++;

#ifdef USE_STATIONDATA
    station.report(time+dt_level[level],level,amr_level[level]);
#endif

#ifdef USE_SLABSTAT
    AmrLevel::get_slabstat_lst().update(amr_level[level],time,dt_level[level]);
#endif
}

void
PMAmr::coarseTimeStep (Real _stop_time)
{
    const Real run_strt = ParallelDescriptor::second();

    // Reset stop time per arg
    if (stop_time != _stop_time) {
      stop_time = _stop_time;
      std::string event_name = "Stop_Time";
      std::map<std::string,EventCoord::Event*>::iterator it=defined_events.find(event_name);
      if (it != defined_events.end()) {
        delete it->second;
      }
      defined_events[event_name] = new EventCoord::TimeEvent(Array<Real>(1,stop_time));
    }
    
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
        std::cout << "\nBEGIN TEMPORAL INTEGRATION, TIME = " << cumtime << '\n' << std::endl;
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
    else
       dt0_before_event_cut = -1;
    // Do time step
    pm_timeStep(0,cumtime,1,1);

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
        int file_name_digits_tmp = file_name_digits;
        file_name_digits = chk_file_digits;
        checkPoint();
        file_name_digits = file_name_digits_tmp;
    }

    if (write_plot)
    {
        int file_name_digits_tmp = file_name_digits;
        file_name_digits = plot_file_digits;
        writePlotFile();
        file_name_digits = file_name_digits_tmp;
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
PMAmr::SetUpMaterialServer()
{
#if 1
  materialFiller = 0;
  return;
#else
  if (materialFiller == 0 || !materialFiller->Initialized()) {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Building MaterialFiller object..." << std::endl;
    }

    int Nlevs = maxLevel() + 1;
    const PArray<Material>& materials = PorousMedia::Materials();
    materialFiller = new MatFiller(geom,refRatio(),materials);
    if (ParallelDescriptor::IOProcessor() && verbose>0) {
      std::cout << "....MaterialFiller object built" << std::endl;
    }
  }
#endif
}


void
PMAmr::initialInit (Real t_start,
                    Real t_stop)
{
    checkInput();

    finest_level = 0;

    int init = true;

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
      // FIXME: Takes only ref[0] here, AND enforces ref-based sub-dt
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
    for (int lev = 0; lev <= finest_level; lev++) {
        amr_level[lev].post_init(t_stop);
    }

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

void PMAmr::InitializeControlEvents()
{
  Real time_eps = 1.e-6; // FIXME: needs to be computed

  // Build time macros
  ParmParse ppa("amr");

  EventCoord& event_coord = PMAmr::eventCoord();

  int n_cmac = ppa.countval("cycle_macros");
  Array<std::string> cmacroNames;
  ppa.getarr("cycle_macros",cmacroNames,0,n_cmac);
  std::map<std::string,int> cmacro_map;
  for (int i=0; i<n_cmac; ++i) {
      std::string prefix = "amr.cycle_macro." + cmacroNames[i];
      ParmParse ppc(prefix);
      std::string type; ppc.get("type",type);
      if (type == "period") {
          int start, period, stop;
          ppc.get("start",start);
          ppc.get("period",period);
          ppc.get("stop",stop);
          defined_events[cmacroNames[i]] = new EventCoord::CycleEvent(start,period,stop);
      }
      else if (type == "cycles" ){
          Array<int> cycles; ppc.getarr("cycles",cycles,0,ppc.countval("cycles"));
          defined_events[cmacroNames[i]] = new EventCoord::CycleEvent(cycles);
      }
      else {
          BoxLib::Abort("Unrecognized cycle macros type");
      }
      cmacro_map[cmacroNames[i]] = i;
  }

  int n_tmac = ppa.countval("time_macros");
  Array<std::string> tmacroNames;
  ppa.getarr("time_macros",tmacroNames,0,n_tmac);
  std::map<std::string,int> tmacro_map;
  for (int i=0; i<n_tmac; ++i) {
      std::string prefix = "amr.time_macro." + tmacroNames[i];
      ParmParse ppt(prefix);
      std::string type; ppt.get("type",type);
      if (type == "period") {
          Real start, period, stop;
          ppt.get("start",start);
          ppt.get("period",period);
          ppt.get("stop",stop);
          defined_events[tmacroNames[i]] = new EventCoord::TimeEvent(start,period,stop);
      }
      else if (type == "times" ){
          Array<Real> times; ppt.getarr("times",times,0,ppt.countval("times"));
          defined_events[tmacroNames[i]] = new EventCoord::TimeEvent(times);
      }
      else {
          BoxLib::Abort("Unrecognized time macros type");
      }
      tmacro_map[tmacroNames[i]] = i;
  }

  ParmParse pp("observation");

  // determine number of observation
  int n_obs = pp.countval("observation");
  std::map<std::string,EventCoord::Event*>::const_iterator eit;

  if (n_obs > 0) {
    observations.resize(n_obs,PArrayManage);
    Array<std::string> obs_names;
    pp.getarr("observation",obs_names,0,n_obs);

    // Get time and cycle macros

    // Get parameters for each observation
    // observation type:0=production,1=mass_fraction,2=mole_fraction,3=saturation
    for (int i=0; i<n_obs; i++) {
      const std::string prefix("observation." + obs_names[i]);
      ParmParse ppr(prefix.c_str());

      std::string obs_type; ppr.get("obs_type",obs_type);
      std::string obs_field; ppr.get("field",obs_field);
      Array<std::string> region_names(1); ppr.get("region",region_names[0]);
      const PArray<Region> obs_regions = PorousMedia::build_region_PArray(region_names);//Should not be a static function

      std::string obs_time_macro, obs_cycle_macro;
      ppr.query("cycle_macro",obs_cycle_macro);
      ppr.query("time_macro",obs_time_macro);

      std::string event_label;
      if (ppr.countval("cycle_macro")>0) {
        eit = defined_events.find(obs_cycle_macro);
        if (eit != defined_events.end()  && eit->second->IsCycle() ) {
          event_label = eit->first;
          event_coord.Register(event_label,eit->second);
        }
        else {
          std::string m = "obs_cycle_macro unrecognized \"" + obs_cycle_macro + "\"";
          BoxLib::Abort(m.c_str());
        }
      }
      else if (ppr.countval("time_macro")>0) {
        eit = defined_events.find(obs_time_macro);
        if (eit != defined_events.end()  && eit->second->IsTime() ) {
          event_label = eit->first;
          event_coord.Register(event_label,eit->second);
        }
        else {
          std::string m = "obs_time_macro unrecognized \"" + obs_time_macro + "\"";
          BoxLib::Abort(m.c_str());
        }
      }
      else {
        std::string m = "Must define either time or cycle macro for observation \"" + obs_names[i] + "\"";
        BoxLib::Abort(m.c_str());
      }

      observations.set(i, new Observation(obs_names[i],obs_field,obs_regions[0],obs_type,event_label));
    }

    // filename for output
    pp.query("output_file",observation_output_file);
  }

  ppa.queryarr("vis_cycle_macros",vis_cycle_macros,0,ppa.countval("vis_cycle_macros"));
  ppa.queryarr("vis_time_macros",vis_time_macros,0,ppa.countval("vis_time_macros"));
  ppa.queryarr("chk_cycle_macros",chk_cycle_macros,0,ppa.countval("chk_cycle_macros"));
  ppa.queryarr("chk_time_macros",chk_time_macros,0,ppa.countval("chk_time_macros"));

  for (int i=0; i<vis_cycle_macros.size(); ++i)
  {
      eit = defined_events.find(vis_cycle_macros[i]);
      if (eit != defined_events.end()  && eit->second->IsCycle() ) {
          event_coord.Register(eit->first,eit->second);
      }
      else {
          std::string m = "vis_cycle_macros contains unrecognized macro name \"" + vis_cycle_macros[i] + "\"";
          BoxLib::Abort(m.c_str());
      }
  }

  for (int i=0; i<vis_time_macros.size(); ++i)
  {
      eit = defined_events.find(vis_time_macros[i]);
      if (eit != defined_events.end()  && eit->second->IsTime() ) {
          event_coord.Register(eit->first,eit->second);
      }
      else {
          std::string m = "vis_time_macros contains unrecognized macro name \"" + vis_time_macros[i] + "\"";
          BoxLib::Abort(m.c_str());
      }
  }
  for (int i=0; i<chk_cycle_macros.size(); ++i)
  {
      eit = defined_events.find(chk_cycle_macros[i]);
      if (eit != defined_events.end()  && eit->second->IsCycle() ) {
          event_coord.Register(eit->first,eit->second);
      }
      else {
          std::string m = "chk_cycle_macros contains unrecognized macro name \"" + chk_cycle_macros[i] + "\"";
          BoxLib::Abort(m.c_str());
      }
  }
  for (int i=0; i<chk_time_macros.size(); ++i)
  {
      eit = defined_events.find(chk_time_macros[i]);
      if (eit != defined_events.end()  && eit->second->IsTime() ) {
          event_coord.Register(eit->first,eit->second);
      }
      else {
          std::string m = "chk_time_macros contains unrecognized macro name \"" + chk_time_macros[i] + "\"";
          BoxLib::Abort(m.c_str());
      }
  }

  //
  // Get run options.
  //
  ParmParse pb("prob");

  int ntps = pb.countval("TPC_Start_Times");
  if (ntps) {
      tpc_start_times.resize(ntps);
      pb.getarr("TPC_Start_Times",tpc_start_times,0,ntps);

      int ndt = pb.countval("TPC_Initial_Time_Step");
      BL_ASSERT(ndt==0 || ndt==ntps);
      tpc_initial_time_steps.resize(ndt,-1);
      if (ndt>0) {
          pb.getarr("TPC_Initial_Time_Step",tpc_initial_time_steps,0,ndt);
      }

      int ndtm = pb.countval("TPC_Initial_Time_Step_Multiplier");
      BL_ASSERT(ndtm==0 || ndtm==ntps);
      tpc_initial_time_step_multipliers.resize(ndt,-1);
      if (ndtm>0) {
          pb.getarr("TPC_Initial_Time_Step_Multiplier",tpc_initial_time_step_multipliers,0,ndtm);
      }

      int ndtmax = pb.countval("TPC_Maximum_Time_Step");
      BL_ASSERT(ndtmax==0 || ndtmax==ntps);
      tpc_initial_time_steps.resize(ndtmax,-1);
      if (ndtmax>0) {
          pb.getarr("TPC_Maximum_Time_Step",tpc_maximum_time_steps,0,ndtmax);
      }
  }

  tpc_labels.resize(ntps);
  for (int i=0; i<ntps; ++i) {
      int ndigits = (int) (std::log10(ntps) + .0001) + 1;
      tpc_labels[i] = BoxLib::Concatenate("Time_Period_Begin_",i,ndigits);
      defined_events[tpc_labels[i]] = new EventCoord::TimeEvent(Array<Real>(1,tpc_start_times[i]));
      PMAmr::eventCoord().Register(tpc_labels[i],defined_events[tpc_labels[i]]);
  }


}
