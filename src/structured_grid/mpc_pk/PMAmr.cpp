/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <ios>
#include <iomanip>
#include <algorithm>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <PMAmr.H>
#include <PorousMedia.H>
#include <Observation.H>
#include <PMAMR_Labels.H>

namespace
{
  const std::string CheckPointVersion("CheckPointVersion_1.0");

  bool initialized = false;

  Real seconds_per_year = 365.25 * 3600 * 24;
}
namespace
{
  std::string visit_plotfile_list_name = "movie.visit";
  //
  // These are all ParmParse'd in, as variables with the same name are
  // not accessible to this derived class...should fix
  // NOTE: Defaults set in Initialize block below
  //
  int  regrid_on_restart;
  int  use_efficient_regrid;
  int  compute_new_dt_on_regrid;
  int  plotfile_on_restart;
  int  mffile_nstreams;
}

std::ostream& operator<<(std::ostream& os, const ExecControl& ec)
{
  os << "{\n";
  os << " label:            " << ec.label << '\n';
  os << " mode:             " << ec.mode  << '\n';
  os << " method:           " << ec.method << '\n';
  os << " max_cycles:       " << ec.max_cycles << '\n';
  os << " start:            " << ec.start << "s (=" << ec.start/(seconds_per_year) << "y)\n";
  os << " end:              " << ec.end << "s (=" << ec.end/(seconds_per_year) << "y)\n";
  os << " init_dt:          ";
  if (ec.init_dt <= 0) {
    os << "<not used>\n";
  }
  else {
    os << ec.init_dt << "s (=" << ec.init_dt/(seconds_per_year) << "y)\n";
  }
  os << " max_dt:           ";
  if (ec.max_dt <= 0) {
    os << "<not used>\n";
  }
  else {
    os << ec.max_dt << "s (=" << ec.max_dt/(seconds_per_year) << "y)\n";
  }
  os << " reduction_factor: " << ec.reduction_factor << '\n';
  os << " increase_factor:  " << ec.increase_factor << '\n';
  os << "}";
  return os;
}

EventCoord PMAmr::event_coord;
std::map<std::string,EventCoord::Event*> PMAmr::defined_events;
bool PMAmr::do_output_time_in_years;
bool PMAmr::attempt_to_recover_failed_step;
int PMAmr::plot_file_digits;
int PMAmr::chk_file_digits;

void
PMAmr::Initialize ()
{
  if (initialized) return;
  //
  // Set all defaults here!!!
  //
  mffile_nstreams          = 1;
  regrid_on_restart        = 0;
  use_efficient_regrid     = 0;
  plotfile_on_restart      = 0;
  compute_new_dt_on_regrid = 0; // FIXME: Avoid cfl calc before velocity field computed on new level, need to revisit
  plot_file_digits         = 5;
  chk_file_digits          = 5;
  do_output_time_in_years  = true;
  attempt_to_recover_failed_step = true;

  BoxLib::ExecOnFinalize(PMAmr::Finalize);

  VisMF::Initialize();

  initialized = true;
}

void
PMAmr::Finalize ()
{
  initialized = false;
}

PMAmr::PMAmr()
  : Amr()
{
    Initialize();
    ParmParse ppa("amr");
    ppa.query("plot_file_digits",plot_file_digits);
    ppa.query("chk_file_digits",chk_file_digits);
    ppa.query("regrid_on_restart",regrid_on_restart);
    ppa.query("use_efficient_regrid",use_efficient_regrid);
    ppa.query("plotfile_on_restart",plotfile_on_restart);
    ppa.query("compute_new_dt_on_regrid",compute_new_dt_on_regrid);
    ppa.query("do_output_time_in_years",do_output_time_in_years);
    ppa.query("attempt_to_recover_failed_step",attempt_to_recover_failed_step);

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
}

static std::string getFilePart(const std::string& path)
{
  std::vector<std::string> parts = BoxLib::Tokenize(path,"/");
  return parts[parts.size()-1];
}

static std::string getDirPart(const std::string& path)
{
  std::vector<std::string> parts = BoxLib::Tokenize(path,"/");
  std::string ret;
  if (path.at(0) == '/') {
    ret = "/";
  }
  else if (parts.size() == 1) {
    ret = "./";
  }

  for (int i=0; i<parts.size()-1; ++i) {
    ret += parts[i];
    if (i!=parts.size()-2) ret += '/';
  }
  return ret;
}

void
PMAmr::UpdateVisitPlotfileList() const
{
  if (ParallelDescriptor::IOProcessor())
  {
    const std::string plot_dir_part = getDirPart(plot_file_root);
    const std::string visit_plotfile_list_fullPath
      = plot_dir_part + "/" + visit_plotfile_list_name;

    std::ofstream ofs(visit_plotfile_list_fullPath.c_str());
    for (int i=0; i<plotfiles_written.size(); ++i) {
      ofs << getFilePart(plotfiles_written[i]) << "/Header" << '\n';
    }
    ofs.close();
  }
}

void
PMAmr::writePlotFile ()
{
  int file_name_digits_tmp = file_name_digits;
  file_name_digits = plot_file_digits;
  Amr::writePlotFile();

  // Not exactly following C++ encapsulation here, but it gets the job done
  plotfiles_written.push_back(BoxLib::Concatenate(plot_file_root,level_steps[0],file_name_digits));
  file_name_digits = file_name_digits_tmp;

  UpdateVisitPlotfileList();
}

void
PMAmr::checkPointObservations () const
{
  if (ParallelDescriptor::IOProcessor())
  {
    std::string chkname = BoxLib::Concatenate(check_file_root,level_steps[0],file_name_digits);
    std::string obsdir = chkname + "/Observations";

    if( ! BoxLib::UtilCreateDirectory(obsdir, 0755)) {
      BoxLib::CreateDirectoryFailed(obsdir);
    }

    std::string obs_ascii = obsdir + "/Header";
    std::string obs_binary = obsdir + "/Data";

    std::ofstream oao;
    oao.open(obs_ascii.c_str(), std::ios::out);
    if (!oao.good())
      BoxLib::FileOpenFailed(obs_ascii);

    std::ofstream obo;
    obo.open(obs_binary.c_str(), std::ios::out|std::ios::app|std::ios::binary);
    if (!obo.good())
      BoxLib::FileOpenFailed(obs_binary);

    oao << observations.size() << '\n';
    Observation::ostream os(oao,obo);
    for (int i=0; i<observations.size(); ++i) {
      observations[i].CheckPoint(os);
    }
  }
}

void
PMAmr::restartObservations (const std::string& chkname)
{
  std::string obsdir = chkname + "/Observations";
  std::string obs_ascii = obsdir + "/Header";
  std::string obs_binary = obsdir + "/Data";

  /*
    The following restart is managed by reading the actually data files on IOProc,
    and broadcasting the ascii and binary streams to all procs so that they can
    all build their own observation set.
   */
  Array<char> aFileCharPtr;
  ParallelDescriptor::ReadAndBcastFile(obs_ascii, aFileCharPtr);
  int Na = aFileCharPtr.size();
  std::string aFileCharPtrString(aFileCharPtr.dataPtr(), Na);
  std::istringstream isa(aFileCharPtrString, std::istringstream::in);

  Array<char> bFileCharPtr;
  ParallelDescriptor::ReadAndBcastFile(obs_binary, bFileCharPtr);
  int Nb = bFileCharPtr.size();
  std::string bFileCharPtrString(bFileCharPtr.dataPtr(), Nb);
  std::istringstream isb(bFileCharPtrString, std::istringstream::in|std::istringstream::binary);

  // Skip file read if the requisite files do not exist in the chkpoint file
  if (isa.good() && isb.good())
  {
    Observation::istream is(isa,isb);
    int n;
    isa >> n; // Get the number of observations as the first entry in the Header file

    BL_ASSERT(n == observations.size() || n==0);      // Assume list of observations has not changed on restart
    for (int i=0; i<n; ++i) {
      const std::string new_event_label = observations[i].event_label;
      observations[i] = Observation(is);             // Replace observation created from input with ones from chkpoint
      observations[i].event_label = new_event_label; // Allow changing time/cycle macros on restart
    }
  }
}

void
PMAmr::checkPoint ()
{
  int file_name_digits_tmp = file_name_digits;
  file_name_digits = chk_file_digits;
  Amr::checkPoint();
  checkPointObservations();
  file_name_digits = file_name_digits_tmp;
}

void
PMAmr::LinkFinalCheckpoint (int step)
{
  if (ParallelDescriptor::IOProcessor()) {
    struct stat statbuff;

    std::string finalCheckpointName = BoxLib::Concatenate(check_file_root,step,chk_file_digits);
    std::string finalCheck_filePart = getFilePart(finalCheckpointName);

    std::string linkName = check_file_root + "_final";

    bool link_exists = ::lstat(linkName.c_str(), &statbuff) != -1;
    if (link_exists) {
      std::cout << "Unlinking \"" << linkName << "\""<< std::endl;
      int ret = unlink(linkName.c_str()); BL_ASSERT(ret == 0);
    }

    std::cout << "Symbolic link, \"" << linkName
	      << "\" created to final checkpoint, \"" << finalCheckpointName
	      << "\""<< std::endl;

    int ret = symlink(finalCheck_filePart.c_str(),linkName.c_str());
    BL_ASSERT(ret == 0);
  }
}

void
PMAmr::RegisterEvent(const std::string& event_label,
                     EventCoord::Event* event = 0)
{
  std::map<std::string,EventCoord::Event*>::const_iterator it=defined_events.find(event_label);
  if (it==defined_events.end()) {
    BL_ASSERT(event!=0);
    defined_events[event_label] = event;
  }
  eventCoord().Register(event_label,defined_events[event_label]);
}

Real
PMAmr::process_events(bool& write_plotfile_after_step,
                      bool& write_checkpoint_after_step,
                      Array<int>& observations_after_step,
                      bool& begin_ecp,
                      EventCoord& event_coord,
                      Real time, Real dt, int iter, int diter)
{
  write_plotfile_after_step = false;
  write_checkpoint_after_step = false;
  begin_ecp = false;

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

      for (int k=0; k<ecs.size(); ++k) {
	if (eventList[j] == ecs[k].label + "_start") {
	  begin_ecp = true;
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
                      bool& begin_ecp_now,
                      EventCoord& event_coord,
                      Real time, int iter)
{
    write_plotfile_now = false;
    write_checkpoint_now = false;
    begin_ecp_now = false;

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

      for (int k=0; k<ecs.size(); ++k) {
	if (eventList[j] == ecs[k].label + "_start") {
	  begin_ecp_now = true;
	}
      }
    }
}

void
PMAmr::init (Real t_start,
             Real t_stop)
{
  InitializeControlEvents();

  if (!restart_chkfile.empty() && restart_chkfile != "init") {

    restart(restart_chkfile);
    setStartTime(0); // FIXME: This needs to be written to the checkpoint

  } else {

    initialInit(t_start, t_stop);

    bool write_plot, write_check, begin_ecp;
    Array<int> initial_observations;

    initial_events(write_plot,write_check,initial_observations,begin_ecp,event_coord,
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

	    layout.SetParent(this);
	    layout.Build(); // Internally destroys itself on rebuild

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
		layout.SetParent(this);
		layout.Build(); // Internally destroys itself on rebuild
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
    // Check to see if should write plotfile.
    // This routine is here so it is done after the restart regrid.
    //
    if (plotfile_on_restart && !(restart_chkfile.empty()) )
    {
	plotfile_on_restart = 0;
        writePlotFile();
    }
    //
    // Advance grids at this level.
    //
    Real dt_taken, dt_suggest;
    PorousMedia& pm = dynamic_cast<PorousMedia&>(amr_level[level]);
    bool step_ok = pm.ml_step_driver(time,iteration,niter,dt_level[level],dt_taken,dt_suggest,attempt_to_recover_failed_step);

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
      BoxLib::Abort("Step failed, and \"attempt_to_recover_failed_step\" is false");
    }

    dt_min[level] = iteration == 1 ? dt_suggest : std::min(dt_min[level],dt_suggest);

    level_steps[level]++;
    level_count[level]++;

    int sum_interval = PorousMedia::SumInterval();
    if (level==0 && sum_interval>0 && level_steps[0]%sum_interval == 0) {
      pm.sum_integrated_quantities();
    }


#ifdef USE_STATIONDATA
    station.report(time+dt_level[level],level,amr_level[level]);
#endif

#ifdef USE_SLABSTAT
    AmrLevel::get_slabstat_lst().update(amr_level[level],time,dt_level[level]);
#endif
}

std::pair<Real,std::string>
PMAmr::convert_time_units(Real t, const std::string& units)
{
  Real t_output;
  std::string units_str;
  std::string units_in = units;
  std::transform(units_in.begin(), units_in.end(), units_in.begin(), toupper);

  if (units_in == "Y") {
    t_output = t/seconds_per_year;
    units_str = "y";
  } else if (units_in == "S") {
    t_output = t;
    units_str = "s";
  }
  else {
    std::cout << "units_str: " << units_in << std::endl;
    BoxLib::Abort();
  }
  return std::pair<Real,std::string>(t_output,units_str);
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

    bool write_plot, write_check, begin_ecp;
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
        dt2_red = process_events(write_plot,write_check,observations_to_process,begin_ecp,event_coord,
                                 cumtime, 2*dt_level[0], level_steps[0] + 1, 1);
        observations_to_process.clear();
        write_check=false;
        write_plot=false;
        begin_ecp=false;
    }

    Real dt_red = process_events(write_plot,write_check,observations_to_process,begin_ecp,event_coord,
                                 cumtime, dt_level[0], level_steps[0], 1);

    // Note: if dt_red > 0, then dt_red == dt2_red
    if (dt2_red > 0  &&  dt_red < 0) {
            dt_red = dt2_red / 2; // Nothing in 1 step, but something in 2
    }

    if (dt_red > 0  &&  dt_red < dt_level[0]) {

        if (begin_ecp) {
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
    else {
       dt0_before_event_cut = -1;
    }

    const ExecControl *ec = GetExecControl(cumtime);
    if (cumtime == ec->start) {
      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Entering Exceution Control Period \""<< ec->label << "\" at time: "
		  << cumtime << "s (=" << cumtime/seconds_per_year << "y)" << std::endl;
	std::cout << *ec << std::endl;
      }
    }

    // Do time step
    pm_timeStep(0,cumtime,1,1);

    cumtime += dt_level[0];

    amr_level[0].postCoarseTimeStep(cumtime);

    static int cnt = 0;
    if (cumtime == ec->end) {
      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Exiting Exceution Control Period \"" << ec->label << "\" at time: "
		  << cumtime << "s (=" << cumtime/seconds_per_year << "y)" << std::endl;
      }
    }

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_stop = ParallelDescriptor::second() - run_strt;

        ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

        if (verbose>3 && ParallelDescriptor::IOProcessor())
            std::cout << "\nCoarse TimeStep time: " << run_stop << '\n' ;

        long min_fab_kilobytes = BoxLib::TotalBytesAllocatedInFabsHWM()/1024;
        long max_fab_kilobytes = min_fab_kilobytes;

        ParallelDescriptor::ReduceLongMin(min_fab_kilobytes,IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_kilobytes,IOProc);
        //
        // Reset to zero to calculate high-water-mark for next timestep.
        //
	BoxLib::ResetTotalBytesAllocatedInFabsHWM();

        if (verbose>3 && ParallelDescriptor::IOProcessor())
            std::cout << "\nFAB kilobyte spread across MPI nodes for timestep: ["
                      << min_fab_kilobytes << " ... " << max_fab_kilobytes << "]\n";
    }

    std::string units_str = do_output_time_in_years ? "Y" : "s";
    std::pair<Real,std::string> t_output = PMAmr::convert_time_units(cumtime,units_str);
    std::pair<Real,std::string> dt_output = PMAmr::convert_time_units(dt_level[0],units_str);
    std::ios_base::fmtflags oldflags = std::cout.flags(); std::cout << std::scientific << std::setprecision(10);
    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << "STEP = "
                 << level_steps[0]
                  << " COMPLETE.  TIME = "
                  << t_output.first << t_output.second
                  << " DT = "
                  << dt_output.first << dt_output.second
                  << "\n\n";
    }
    if (record_run_info && ParallelDescriptor::IOProcessor())
    {
        runlog << "STEP = "
               << level_steps[0]
               << " TIME = "
               << t_output.first << t_output.second
               << " DT = "
               << dt_output.first << dt_output.second
               << "\n\n";
    }
    if (record_run_info_terse && ParallelDescriptor::IOProcessor())
        runlog_terse << level_steps[0] << " " << cumtime << " " << dt_level[0] << '\n';
    std::cout.flags(oldflags);

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
PMAmr::initialInit (Real              strt_time,
		    Real              stop_time,
		    const BoxArray*   lev0_grids,
		    const Array<int>* pmap)
{
  setStartTime(strt_time);
  InitializeInit(strt_time, stop_time, lev0_grids, pmap);

  // This is a subtlety, but in the case where we are initializing the data
  //   from a plotfile, we want to use the time read in from the plotfile as
  //   the start time instead of using "strt_time".
  // The Amr data "cumtime" has been set in InitializeInit; if we are restarting
  //   from a plotfile, then cumtime must be re-defined in that initialization routine.
  //   Thus here we pass "cumtime" rather than "strt_time" to FinalizeInit.
  FinalizeInit  (cumtime, stop_time);
}

void
PMAmr::InitializeInit(Real              strt_time,
		      Real              stop_time,
		      const BoxArray*   lev0_grids,
		      const Array<int>* pmap)
{
  //Amr::InitializeInit(strt_time,stop_time,lev0_grids,pmap);
    BL_COMM_PROFILE_NAMETAG("PMAmr::initialInit TOP");
    checkInput();
    //
    // Generate internal values from user-supplied values.
    //
    finest_level = 0;
    //
    // Init problem dependent data.
    //
    if (!probin_file.empty()) {
      //readProbinFile(init);
    }

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
    // Define base level grids.  Note that if we are restarting from a plotfile, this
    //    routine will call the level 0 AmrLevel initialization which will overwrite cumtime.
    //
    defBaseLevel(strt_time, lev0_grids, pmap);
}

void
PMAmr::FinalizeInit (Real              strt_time,
		     Real              stop_time)
{
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
    bool write_plot, write_check, begin_ecp;
    Array<int> observations_to_process;

    Real dt_red = process_events(write_plot,write_check,observations_to_process,begin_ecp,event_coord,
                                 cumtime, dt_level[0], level_steps[0], 1);
    if (dt_red > 0  &&  dt_red < dt0) {
        dt_min[0]  = dt_level[0];
    }

    for (int lev = 1; lev <= max_level; lev++)
    {
        dt0           /= n_cycle[lev];
        dt_level[lev]  = dt0;
        dt_min[lev]    = dt_level[lev];
    }

    if (max_level > 0)
        bldFineLevels(strt_time);

    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev].setTimeLevel(strt_time,dt_level[lev],dt_level[lev]);

    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev].post_regrid(0,finest_level);

    for (int lev = 0; lev <= finest_level; lev++)
    {
        level_steps[lev] = 0;
        level_count[lev] = 0;
    }

    //
    // Perform any special post_initialization operations.
    //
    for (int lev = 0; lev <= finest_level; lev++)
        amr_level[lev].post_init(stop_time);

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

#if 0
    // Write initial plts and checkpoints no matter what
    {
        int file_name_digits_tmp = file_name_digits;
        file_name_digits = chk_file_digits;
        checkPoint();
        file_name_digits = file_name_digits_tmp;
    }

    {
        int file_name_digits_tmp = file_name_digits;
        file_name_digits = plot_file_digits;
        writePlotFile();
        file_name_digits = file_name_digits_tmp;
    }
#endif

}

// FIXME: HACK  This is copied from Amr.cpp and inserted here.  IT IS NOT A VIRTUAL FUNCTION, SO BEWARE.
// We wanted to avoid reading the probin file, and this is the only way to do everything else but that part.
void
PMAmr::restart (const std::string& filename)
{
    BL_PROFILE("PMAmr::restart()");

    // Just initialize this here for the heck of it
    which_level_being_advanced = -1;

    Real dRestartTime0 = ParallelDescriptor::second();

    DistributionMapping::Initialize();

#if 0
// HACK: Incompatible with slightly older BoxLib version currently distributed with Amanzi....
    if(DistributionMapping::strategy() == DistributionMapping::PFC) {
      Array<IntVect> refRatio;
      Array<BoxArray> allBoxes;
      DistributionMapping::ReadCheckPointHeader(filename, refRatio, allBoxes);
      DistributionMapping::PFCMultiLevelMap(refRatio, allBoxes);
    }
#endif

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "DMCache size = " << DistributionMapping::CacheSize() << std::endl;
      DistributionMapping::CacheStats(std::cout);
    }
    ParallelDescriptor::Barrier();


    VisMF::SetMFFileInStreams(mffile_nstreams);

    int i;

    if (verbose > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "restarting calculation from file: " << filename << std::endl;

    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "RESTART from file = " << filename << '\n';
    //
    // Init problem dependent data.
    //
    // THE WHOLE READON FOR THIS ABOMINATION IS THE FOLLOWING COMMENT
    // readProbinFile(init);
    //
    // Start calculation from given restart file.
    //
    if (record_run_info && ParallelDescriptor::IOProcessor())
        runlog << "RESTART from file = " << filename << '\n';
    //
    // Open the checkpoint header file for reading.
    //
    std::string File = filename;

    File += '/';
    File += "Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    Array<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);
    //
    // Read global data.
    //
    // Attempt to differentiate between old and new CheckPointFiles.
    //
    int         spdim;
    bool        new_checkpoint_format = false;
    std::string first_line;

    std::getline(is,first_line);

    if (first_line == CheckPointVersion)
    {
        new_checkpoint_format = true;
        is >> spdim;
    }
    else
    {
        spdim = atoi(first_line.c_str());
    }

    if (spdim != BL_SPACEDIM)
    {
        std::cerr << "PMAmr::restart(): bad spacedim = " << spdim << '\n';
        BoxLib::Abort();
    }

    is >> cumtime;
    int mx_lev;
    is >> mx_lev;
    is >> finest_level;

    Array<Box> inputs_domain(max_level+1);
    for (int lev = 0; lev <= max_level; lev++)
    {
       Box bx(geom[lev].Domain().smallEnd(),geom[lev].Domain().bigEnd());
       inputs_domain[lev] = bx;
    }

    if (max_level >= mx_lev) {

       for (i = 0; i <= mx_lev; i++) is >> geom[i];
       for (i = 0; i <  mx_lev; i++) is >> ref_ratio[i];
       for (i = 0; i <= mx_lev; i++) is >> dt_level[i];

       if (new_checkpoint_format)
       {
           for (i = 0; i <= mx_lev; i++) is >> dt_min[i];
       }
       else
       {
           for (i = 0; i <= mx_lev; i++) dt_min[i] = dt_level[i];
       }

       Array<int>  n_cycle_in;
       n_cycle_in.resize(mx_lev+1);
       for (i = 0; i <= mx_lev; i++) is >> n_cycle_in[i];
       bool any_changed = false;

       for (i = 0; i <= mx_lev; i++)
           if (n_cycle[i] != n_cycle_in[i])
           {
               any_changed = true;
               if (verbose > 0 && ParallelDescriptor::IOProcessor())
                   std::cout << "Warning: n_cycle has changed at level " << i <<
                                " from " << n_cycle_in[i] << " to " << n_cycle[i] << std::endl;;
           }

       // If we change n_cycle then force a full regrid from level 0 up
       if (max_level > 0 && any_changed)
       {
           level_count[0] = regrid_int[0];
           if ((verbose > 0) && ParallelDescriptor::IOProcessor())
               std::cout << "Warning: This forces a full regrid " << std::endl;
       }


       for (i = 0; i <= mx_lev; i++) is >> level_steps[i];
       for (i = 0; i <= mx_lev; i++) is >> level_count[i];

       //
       // Set bndry conditions.
       //
       if (max_level > mx_lev)
       {
           for (i = mx_lev+1; i <= max_level; i++)
           {
               dt_level[i]    = dt_level[i-1]/n_cycle[i];
               level_steps[i] = n_cycle[i]*level_steps[i-1];
               level_count[i] = 0;
           }

           // This is just an error check
           if (!sub_cycle)
           {
               for (i = 1; i <= finest_level; i++)
               {
                   if (dt_level[i] != dt_level[i-1])
                      BoxLib::Error("restart: must have same dt at all levels if not subcycling");
               }
           }
       }

       if (regrid_on_restart && max_level > 0)
       {
           if (regrid_int[0] > 0)
               level_count[0] = regrid_int[0];
           else
               BoxLib::Error("restart: can't have regrid_on_restart and regrid_int <= 0");
       }

       checkInput();
       //
       // Read levels.
       //
       int lev;
       for (lev = 0; lev <= finest_level; lev++)
       {
           amr_level.set(lev,(*levelbld)());
           amr_level[lev].restart(*this, is);
       }
       //
       // Build any additional data structures.
       //
       for (lev = 0; lev <= finest_level; lev++)
           amr_level[lev].post_restart();

    } else {

       BoxLib::Abort("PMAmr::restart(): max_level is lower on restart than in checkpoint file.  Adjust refinement criteria to achieve desired max_level");

#if 0
       // Remove this robustification since material properties with a different max_level will not be consistent

       if (ParallelDescriptor::IOProcessor())
          BoxLib::Warning("PMAmr::restart(): max_level is lower than before");

       int new_finest_level = std::min(max_level,finest_level);

       finest_level = new_finest_level;

       // These are just used to hold the extra stuff we have to read in.
       Geometry   geom_dummy;
       Real       real_dummy;
       int         int_dummy;
       IntVect intvect_dummy;

       for (i = 0          ; i <= max_level; i++) is >> geom[i];
       for (i = max_level+1; i <= mx_lev   ; i++) is >> geom_dummy;

       for (i = 0        ; i <  max_level; i++) is >> ref_ratio[i];
       for (i = max_level; i <  mx_lev   ; i++) is >> intvect_dummy;

       for (i = 0          ; i <= max_level; i++) is >> dt_level[i];
       for (i = max_level+1; i <= mx_lev   ; i++) is >> real_dummy;

       if (new_checkpoint_format)
       {
           for (i = 0          ; i <= max_level; i++) is >> dt_min[i];
           for (i = max_level+1; i <= mx_lev   ; i++) is >> real_dummy;
       }
       else
       {
           for (i = 0; i <= max_level; i++) dt_min[i] = dt_level[i];
       }

       for (i = 0          ; i <= max_level; i++) is >> n_cycle[i];
       for (i = max_level+1; i <= mx_lev   ; i++) is >> int_dummy;

       for (i = 0          ; i <= max_level; i++) is >> level_steps[i];
       for (i = max_level+1; i <= mx_lev   ; i++) is >> int_dummy;

       for (i = 0          ; i <= max_level; i++) is >> level_count[i];
       for (i = max_level+1; i <= mx_lev   ; i++) is >> int_dummy;

       if (regrid_on_restart && max_level > 0)
       {
           if (regrid_int[0] > 0)
               level_count[0] = regrid_int[0];
           else
               BoxLib::Error("restart: can't have regrid_on_restart and regrid_int <= 0");
       }

       checkInput();

       //
       // Read levels.
       //
       int lev;
       for (lev = 0; lev <= new_finest_level; lev++)
       {
           amr_level.set(lev,(*levelbld)());
           amr_level[lev].restart(*this, is);
       }
       //
       // Build any additional data structures.
       //
       for (lev = 0; lev <= new_finest_level; lev++)
           amr_level[lev].post_restart();
#endif
    }

    for (int lev = 0; lev <= finest_level; lev++)
    {
       Box restart_domain(geom[lev].Domain());
       if (! (inputs_domain[lev] == restart_domain) )
       {
          if (ParallelDescriptor::IOProcessor())
          {
             std::cout << "Problem at level " << lev << '\n';
             std::cout << "Domain according to     inputs file is " <<  inputs_domain[lev] << '\n';
             std::cout << "Domain according to checkpoint file is " << restart_domain      << '\n';
             std::cout << "PMAmr::restart() failed -- box from inputs file does not equal box from restart file" << std::endl;
          }
          BoxLib::Abort();
       }
    }

#ifdef USE_STATIONDATA
    station.init(amr_level, finestLevel());
    station.findGrid(amr_level,geom);
#endif

    // Read in observations
    restartObservations(filename);

    if (verbose > 0)
    {
        Real dRestartTime = ParallelDescriptor::second() - dRestartTime0;

        ParallelDescriptor::ReduceRealMax(dRestartTime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Restart time = " << dRestartTime << " seconds." << '\n';
    }
}

void PMAmr::InitializeControlEvents()
{
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
      else if (type == "cycles" ) {
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
      else if (type == "times" ) {
          Array<Real> times; ppt.getarr("times",times,0,ppt.countval("times"));
          defined_events[tmacroNames[i]] = new EventCoord::TimeEvent(times);
      }
      else {
          BoxLib::Abort("Unrecognized time macros type");
      }
      tmacro_map[tmacroNames[i]] = i;
  }

  ParmParse ppo("observation");

  // determine number of observation
  int n_obs = ppo.countval("observation");
  std::map<std::string,EventCoord::Event*>::const_iterator eit;

  if (n_obs > 0) {
    observations.resize(n_obs,PArrayManage);
    Array<std::string> obs_names;
    ppo.getarr("observation",obs_names,0,n_obs);

    // Get time and cycle macros

    // Get parameters for each observation
    // observation type:0=production,1=mass_fraction,2=mole_fraction,3=saturation
    for (int i=0; i<n_obs; i++) {
      RegionManager* region_manager = PorousMedia::GetRegionManager();
      if (region_manager == 0) {
        BoxLib::Abort("static Region manager must be set up prior to reading observations");
      }

      const std::string prefix("observation." + obs_names[i]);
      ParmParse ppor(prefix.c_str());

      std::string obs_type; ppor.get("obs_type",obs_type);
      std::string obs_field; ppor.get("field",obs_field);

      Array<std::string> obs_time_macros, obs_cycle_macros;
      int ntm = ppor.countval("time_macros");
      if (ntm>0) {
	ppor.getarr("time_macros",obs_time_macros,0,ntm);
      }
      int ncm = ppor.countval("cycle_macros");
      if (ncm>0) {
	ppor.getarr("cycle_macros",obs_cycle_macros,0,ncm);
      }

      std::string event_label;
      if (ncm) {
	for (int j=0; j<obs_cycle_macros.size(); ++j) {
	  eit = defined_events.find(obs_cycle_macros[j]);
	  if (eit != defined_events.end()  && eit->second->IsCycle() ) {
	    event_label = eit->first;
	    RegisterEvent(event_label,eit->second);
	  }
	  else {
	    std::string m = "obs_cycle_macro unrecognized \"" + obs_cycle_macros[j] + "\"";
	    BoxLib::Abort(m.c_str());
	  }
	}
      }
      else if (ntm) {
	for (int j=0; j<obs_time_macros.size(); ++j) {
	  eit = defined_events.find(obs_time_macros[j]);
	  if (eit != defined_events.end()  && eit->second->IsTime() ) {
	    event_label = eit->first;
	    RegisterEvent(event_label,eit->second);
	  }
	  else {
	    std::string m = "obs_time_macro unrecognized \"" + obs_time_macros[j] + "\"";
	    BoxLib::Abort(m.c_str());
	  }
        }
      }
      else {
        std::string m = "Must define either time or cycle macros for observation \"" + obs_names[i] + "\"";
        BoxLib::Abort(m.c_str());
      }

      std::string region_name; ppor.get("region",region_name);
      observations.set(i, new Observation(obs_names[i],obs_field,region_name,obs_type,event_label));
    }

    // filename for output
    ppo.query("output_file",observation_output_file);
  }

  ppa.queryarr("viz_cycle_macros",vis_cycle_macros,0,ppa.countval("viz_cycle_macros"));
  ppa.queryarr("viz_time_macros",vis_time_macros,0,ppa.countval("viz_time_macros"));
  ppa.queryarr("chk_cycle_macros",chk_cycle_macros,0,ppa.countval("chk_cycle_macros"));
  ppa.queryarr("chk_time_macros",chk_time_macros,0,ppa.countval("chk_time_macros"));

  for (int i=0; i<vis_cycle_macros.size(); ++i)
  {
      eit = defined_events.find(vis_cycle_macros[i]);
      if (eit != defined_events.end()  && eit->second->IsCycle() ) {
          RegisterEvent(eit->first,eit->second);
      }
      else {
          std::string m = "viz_cycle_macros contains unrecognized macro name \"" + vis_cycle_macros[i] + "\"";
          BoxLib::Abort(m.c_str());
      }
  }

  for (int i=0; i<vis_time_macros.size(); ++i)
  {
      eit = defined_events.find(vis_time_macros[i]);
      if (eit != defined_events.end()  && eit->second->IsTime() ) {
          RegisterEvent(eit->first,eit->second);
      }
      else {
          std::string m = "viz_time_macros contains unrecognized macro name \"" + vis_time_macros[i] + "\"";
          BoxLib::Abort(m.c_str());
      }
  }
  for (int i=0; i<chk_cycle_macros.size(); ++i)
  {
    eit = defined_events.find(chk_cycle_macros[i]);
    if (eit != defined_events.end()  && eit->second->IsCycle() ) {
          RegisterEvent(eit->first,eit->second);
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
          RegisterEvent(eit->first,eit->second);
      }
      else {
          std::string m = "chk_time_macros contains unrecognized macro name \"" + chk_time_macros[i] + "\"";
          BoxLib::Abort(m.c_str());
      }
  }

  //
  // Get execution control options
  //
  ParmParse pp;
  int necs = pp.countval("exec_controls");
  if (necs) {
    ecs.resize(necs);
    Array<std::string> ec_names(necs);
    pp.getarr("exec_controls",ec_names,0,necs);
    for (int i=0; i<necs; ++i) {
      std::string prefix = "exec_control." + ec_names[i];
      ParmParse ppe(prefix);
      ExecControl& ec = ecs[i];
      ec.label = ec_names[i];
      ppe.get("start",ec.start);
      ppe.get("end",ec.end);
      ppe.get("init_dt",ec.init_dt);
      ppe.get("max_dt",ec.max_dt);
      ppe.get("reduction_factor",ec.reduction_factor);
      ppe.get("increase_factor",ec.increase_factor);
      ppe.get("max_cycles",ec.max_cycles);
      ppe.get("mode",ec.mode);
      ppe.get("method",ec.method);
    }
  }

  for (int i=0; i<necs; ++i) {
    const ExecControl& ec = ecs[i];
    std::string start_label = ec.label + "_start";
    defined_events[start_label] = new EventCoord::TimeEvent(Array<Real>(1,ec.start));
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Registering event: " << start_label << " " << ec.start << std::endl;
    }
    RegisterEvent(start_label,defined_events[start_label]);
  }

  {
    const ExecControl& ec = ecs[necs-1];
    stop_time = ec.end;
    std::string max_time_name = "Stop_Time";
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Registering event: " << max_time_name << " " << stop_time << std::endl;
    }
    defined_events[max_time_name] = new EventCoord::TimeEvent(Array<Real>(1,ec.end));
    RegisterEvent(max_time_name,defined_events[max_time_name]);

    max_step = ec.max_cycles;
    std::string max_cycles_name = "Max_Cycles";
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Registering event: " << max_cycles_name << " " << max_step << std::endl;
    }
    defined_events[max_cycles_name] = new EventCoord::CycleEvent(Array<int>(1,max_step));
    RegisterEvent(max_cycles_name,defined_events[max_cycles_name]);
  }

}

const ExecControl*
PMAmr::GetExecControl(Real time_value) const
{
  for (int i=0; i<ecs.size(); ++i) {
    const ExecControl& ec = ecs[i];
    if (time_value >= ec.start && time_value < ec.end) return &ec;
  }
  return 0;
}

void PMAmr::FlushObservations() {
  if (observation_output_file != "") {
    if (observations.size() && ParallelDescriptor::IOProcessor()) {
      std::cout << "Writing observations to \"" << observation_output_file << "\"" << std::endl;
      std::ofstream ofs;
      ofs.open(observation_output_file.c_str(),std::ios::out);
      FlushObservations(ofs);
      ofs.close();
    }
  }
}

void PMAmr::FlushObservations(std::ostream& out)
{
  // print out observations
  if (observations.size() && ParallelDescriptor::IOProcessor()) {

    const int old_prec = out.precision(16);
    out.setf(std::ios::scientific);

    out << "Observation Name, Region, Functional, Variable, Time, Value\n";
    out << "===========================================================\n";

    for (int i=0; i<observations.size(); ++i) {
      const std::map<int,Real> vals = observations[i].vals;
      for (std::map<int,Real>::const_iterator it=vals.begin();it!=vals.end(); ++it) {
        int j = it->first;
        out << Amanzi::AmanziInput::GlobalData::AMR_to_Amanzi_label_map[observations[i].name]
            << ", " << Amanzi::AmanziInput::GlobalData::AMR_to_Amanzi_label_map[observations[i].region_name]
            << ", " << observations[i].obs_type
            << ", " << observations[i].field
            << ", " << observations[i].times[j]
            << ", " << it->second << std::endl;
      }
    }
    out.precision(old_prec);
  }
}
