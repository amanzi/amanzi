#include <Utility.H>
#include <ccse-mpi.H>
#include <EventCoord.H>

#ifdef _OPENMP
#include "omp.h"
#endif

static double epsilon = 1.e-10;

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    EventCoord event_coord;

    std::string cm_1_name = "CM1";
    int cm_1_start = 5;
    int cm_1_period = 7;
    int cm_1_stop = 100;
    event_coord.InsertCycleEvent(cm_1_name, cm_1_start, cm_1_period, cm_1_stop);

    std::string cm_2_name = "CM2";
    int cm_2_start = 3;
    int cm_2_period = 9;
    int cm_2_stop = -1;
    event_coord.InsertCycleEvent(cm_2_name, cm_2_start, cm_2_period, cm_2_stop);

    std::string cm_3_name = "CM3";
    Array<int> cm_3_cycles;
    cm_3_cycles.push_back(8);
    cm_3_cycles.push_back(9);
    cm_3_cycles.push_back(11);
    cm_3_cycles.push_back(19);
    event_coord.InsertCycleEvent(cm_3_name, cm_3_cycles);


    std::string tm_1_name = "TM1";
    Real tm_1_start = 4;
    Real tm_1_period = 4;
    Real tm_1_stop = -1;
    event_coord.InsertTimeEvent(tm_1_name, tm_1_start, tm_1_period, tm_1_stop);

    std::string tm_2_name = "TM2";
    Array<Real> tm_2_times;
    tm_2_times.push_back(4.5);
    tm_2_times.push_back(4.51);
    tm_2_times.push_back(4.52);
    tm_2_times.push_back(8.9);
    event_coord.InsertTimeEvent(tm_2_name, tm_2_times);

    // build test loop
    Real time_new, time_old = 0;
    const Real dt_def = 0.5;
    const Real dt_eps = dt_def*1.e-8;
    int diter = 1;
    int iterMAX = 25;
    Array<std::string> eventList;
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
        //  std::cout << "t=" << time_new << " dt = " << time_new - time_old << " iter=" << iter+diter ;
        //  for (int j=0; j<eventList.size(); ++j) {
        //      std::cout << ", event: " << eventList[j] << " ";
        //}
        //  std::cout << std::endl;

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
