#include "amanzi_structured_grid_simulation_driver.H"
#include "ParmParse.H"
#include "Amr.H"
#include "PorousMedia.H"



Simulator::ReturnType
AmanziStructuredGridSimulationDriver::Run (const MPI_Comm&               mpi_comm,
                                           const Teuchos::ParameterList& input_parameter_list,
                                           ObservationData&              output_observations)
{
    int argc=0;
    char** argv;

    std::string PPfile = Teuchos::getParameter<std::string>(input_parameter_list, "PPfile");

    BoxLib::Initialize(argc,argv,false,mpi_comm);
    ParmParse::Initialize(argc,argv,PPfile.c_str());

    const Real run_strt = ParallelDescriptor::second();

    int  max_step;
    Real strt_time;
    Real stop_time;

    ParmParse pp;

    max_step  = -1;    
    strt_time =  0.0;  
    stop_time = -1.0;  

    pp.query("max_step",max_step);
    pp.query("strt_time",strt_time);
    pp.query("stop_time",stop_time);

    if (strt_time < 0.0)
        BoxLib::Abort("MUST SPECIFY a non-negative strt_time");

    if (max_step < 0 && stop_time < 0.0)
    {
        BoxLib::Abort(
            "Exiting because neither max_step nor stop_time is non-negative.");
    }


    Amr* amrptr = new Amr;

    amrptr->init(strt_time,stop_time);

    while ( amrptr->okToContinue()           &&
           (amrptr->levelSteps(0) < max_step || max_step < 0) &&
           (amrptr->cumTime() < stop_time || stop_time < 0.0) )
    {
        amrptr->coarseTimeStep(stop_time);
    }

    // Process the observations
    const Array<Observation>& observation_array = PorousMedia::TheObservationArray();

    delete amrptr;

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_stop = ParallelDescriptor::second() - run_strt;

    ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

    if (ParallelDescriptor::IOProcessor())
        std::cout << "Run time = " << run_stop << std::endl;

    BoxLib::Finalize(false); // Calling routine responsible for MPI_Finalize call

    return Simulator::SUCCESS;
}
