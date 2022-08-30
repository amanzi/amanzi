#include <Utility.H>
#include <ccse-mpi.H>
#include <EventCoord.H>

// closing DSO objects
#include "VerboseObject_objs.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    EventCoord event_coord;

    std::string cm_1_name = "CM1";
    int cm_1_start = 5;
    int cm_1_period = 7;
    int cm_1_stop = 100;
    EventCoord::CycleEvent ce1(cm_1_start, cm_1_period, cm_1_stop);
    event_coord.Register(cm_1_name, &ce1);

    std::string cm_2_name = "CM2";
    int cm_2_start = 3;
    int cm_2_period = 9;
    int cm_2_stop = -1;
    EventCoord::CycleEvent ce2(cm_2_start, cm_2_period, cm_2_stop);
    event_coord.Register(cm_2_name, &ce2);

    std::string cm_3_name = "CM3";
    Array<int> cm_3_cycles;
    cm_3_cycles.push_back(0);
    cm_3_cycles.push_back(8);
    cm_3_cycles.push_back(9);
    cm_3_cycles.push_back(11);
    cm_3_cycles.push_back(19);
    EventCoord::CycleEvent ce3(cm_3_cycles);
    event_coord.Register(cm_3_name, &ce3);


    std::string tm_1_name = "TM1";
    Real tm_1_start = 0;
    Real tm_1_period = 4;
    Real tm_1_stop = -1;
    EventCoord::TimeEvent te1(tm_1_start, tm_1_period, tm_1_stop);
    event_coord.Register(tm_1_name, &te1);

    std::string tm_2_name = "TM2";
    Array<Real> tm_2_times;
    tm_2_times.push_back(0);
    tm_2_times.push_back(4.5);
    tm_2_times.push_back(4.51);
    tm_2_times.push_back(4.52);
    tm_2_times.push_back(8.9);
    EventCoord::TimeEvent te2(tm_2_times);
    event_coord.Register(tm_2_name, &te2);

    // build test loop
    Real time_new, time_old = 0;
    const Real dt_def = 0.5;
    const Real dt_eps = dt_def*1.e-8;
    int iter_init = 0;
    int diter = 1;
    int iterMAX = 25;
    bool verbose = true;

    Array<std::string> eventList = event_coord.InitEvent(time_old,iter_init);
    if (verbose && eventList.size()>0) {
      std::cout << "tinit = " << time_old;
      for (int j=0; j<eventList.size(); ++j) {
        std::cout << ", event: " << eventList[j] << " ";
      }
      std::cout << std::endl;
    }

    for (int iter=0; iter<iterMAX; iter+=diter) {

        std::pair<Real,Array<std::string> > nextEvent = event_coord.NextEvent(time_old,dt_def,iter,diter);

        if (nextEvent.second.size()) 
        {
            Real dt_new = (nextEvent.first > 0 ? nextEvent.first : dt_def);
            time_new = time_old + dt_new;
            eventList = nextEvent.second;
        }
        else {
            time_new = time_old + dt_def;
            eventList.clear();
        }

        // The time stepper advances from time_old to time_new, arriving at cycle=iter+dcycle
        // If time_old+dt_def would step over an event, we will reduce the interval, and set
        // time_new to when that event(s) would occur.  At this new cycle and time, we would
        // typically then process the list of events that have occurred:
        // 
        if (verbose) {
          std::cout << "told = " << time_old << " tnew = " << time_new << " dt = " << time_new - time_old << " iter=" << iter+diter ;
          for (int j=0; j<eventList.size(); ++j) {
            std::cout << ", event: " << eventList[j] << " ";
          }
          std::cout << std::endl;
        }

        // increment the base time and continue
        time_old = time_new;
    }

    if (std::abs(time_new - 11.4) > dt_eps) {
        std::cerr << "Event Coordinator hosed!" << std::endl;
        throw std::exception();
    }

    BoxLib::Finalize();
    
    return 0;

}
