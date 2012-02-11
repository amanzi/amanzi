#include <PMAmr.H>
#include <PorousMedia.H>
#include <Observation.H>

EventCoord PMAmr::event_coord;

void
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


    std::pair<Real,Array<std::string> > nextEvent = event_coord.NextEvent(time,dt,iter, diter);
    PArray<Observation>& observations = PorousMedia::TheObservationArray();

    if (nextEvent.second.size()) 
    {
        // Process event
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
}


void
PMAmr::init (Real strt_time,
             Real stop_time)
{
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
                       strt_time, dt_level[0], level_steps[0], 0);

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
    if (level_steps[0] > 0)
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
    }


    bool write_plot, write_check;
    Array<int> observations_to_process;
    process_events(write_plot,write_check,observations_to_process,event_coord,
                   cumtime, dt_level[0], level_steps[0], 1);
    
    // Do time step
    timeStep(0,cumtime,1,1,stop_time);

    cumtime += dt_level[0];

    amr_level[0].postCoarseTimeStep(cumtime);
    
    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_stop = ParallelDescriptor::second() - run_strt;

        ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "\nCoarse TimeStep time: " << run_stop << '\n' ;

        long min_fab_bytes = BoxLib::total_bytes_allocated_in_fabs_hwm;
        long max_fab_bytes = BoxLib::total_bytes_allocated_in_fabs_hwm;

        ParallelDescriptor::ReduceLongMin(min_fab_bytes,IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_bytes,IOProc);
        //
        // Reset to zero to calculate high-water-mark for next timestep.
        //
        BoxLib::total_bytes_allocated_in_fabs_hwm = 0;

        if (ParallelDescriptor::IOProcessor())
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

    
    PArray<Observation>& observations = PorousMedia::TheObservationArray();
    if (observations_to_process.size()) {
        std::cout << " Process observations ";
        for (int i=0; i<observations_to_process.size(); ++i) {
            std::cout << observations[observations_to_process[i]].name << " " << std::endl;
        }
    }

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

bool
EventCoord::CycleEvent::ThisEventDue(int cycle, int dCycle) const
{
    if (type ==SPS)
    {
        if ( ( ( stop < 0 ) || ( cycle < stop ) )
             && ( cycle >= start ) )
        {
            return ( (period==1) || ((cycle-start+dCycle)%period == 0) );
        }
    }
    else {
        for (int i=0; i<cycles.size(); ++i) {
            if (cycles[i] == cycle+1) {
                return true;
            }
        }
    }
    return false;
}

bool
EventCoord::TimeEvent::ThisEventDue(Real time, Real dt, Real& dt_red) const
{
    dt_red = dt;
    if (type ==SPS)
    {
        if ( ( ( stop < 0 ) || ( time < stop ) )
             && ( time >= start ) )
        {
            Real intPart;
            if ( modf( (time-start+dt)/period, &intPart ) == 0 )
            {
                return true;
            }
        }
    }
    else {
        for (int i=0; i<times.size(); ++i) {
            if (times[i] > time  &&  times[i] <= time + dt) {
                dt_red = std::min(dt, times[i]-time);
                return true;
            }
            if (times[i] == time) {
                Real max_dt = 0;
                for (int j=0; j<times.size(); ++j) {
                    if (time < times[j]) {
                        max_dt = std::max(max_dt,times[j]-time);
                    }
                }
                if (max_dt == 0) {
                    return false;
                }
                else {
                    Real min_dt = max_dt;
                    for (int j=0; j<times.size(); ++j) {
                        if (time < times[j]) {
                            min_dt = std::min(min_dt, times[j] - time);
                        }
                    }
                    dt_red = std::min(dt, min_dt);
                }
                return true;
            }
        }
    }
    dt_red = -1;
    return false;
}


std::pair<Real,Array<std::string> >
EventCoord::NextEvent(Real time, Real dt, int cycle, int dcycle) const
{
    Array<std::string> events;
    Real delta_time = -1;
    for (std::map<std::string,CycleEvent>::const_iterator it=cycleEvents.begin(); it!=cycleEvents.end(); ++it) {
        const std::string& name = it->first;
        const CycleEvent& event = it->second;
        CycleEvent::CycleType type = event.Type();
        if (event.ThisEventDue(cycle,dcycle)) {
            events.push_back(name);
        }
    }

    for (std::map<std::string,TimeEvent>::const_iterator it=timeEvents.begin(); it!=timeEvents.end(); ++it) {
        const std::string& name = it->first;
        const TimeEvent& event = it->second;
        TimeEvent::TimeType type = event.Type();
        Real dt_red;
        if (event.ThisEventDue(time,dt,dt_red)) {
            events.push_back(name);
            delta_time = (delta_time < 0  ?  dt_red  :  std::min(delta_time, dt_red));
        }
    }

    return std::pair<Real,Array<std::string> > (delta_time,events);
}

void
EventCoord::InsertCycleEvent(const std::string& label, int start, int period, int stop)
{
    cycleEvents[label] = CycleEvent(start,period,stop);
}

void
EventCoord::InsertCycleEvent(const std::string& label, const Array<int>& cycles)
{
    cycleEvents[label] = CycleEvent(cycles);
}

void
EventCoord::InsertTimeEvent(const std::string& label, Real start, Real period, Real stop)
{
    timeEvents[label] = TimeEvent(start,period,stop);
}

void
EventCoord::InsertTimeEvent(const std::string& label, const Array<Real>& times)
{
    timeEvents[label] = TimeEvent(times);
}

