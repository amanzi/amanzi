#ifndef ENABLE_Structured
#define ENABLE_Structured 1
#endif

#include "amanzi_structured_grid_simulation_driver.H"
#include "ParmParse.H"
#include "InputConverterS.hh"
#include "PMAmr.H"
#include "PMAMR_Labels.H"
#include "PorousMedia.H"

#include "ParmParseHelpers.H"

static std::map<std::string,std::string>& AMR_to_Amanzi_label_map = Amanzi::AmanziInput::AMRToAmanziLabelMap();

// Minimum allowable input file version
static int amanzi_min_major_version = 1;
static int amanzi_min_minor_version = 2;
static int amanzi_min_patch_version = 2;

static
void
EnsureFolderExists(const std::string& full_path)
{
  // Find folder name first, and ensure folder exists
  // FIXME: Will fail on Windows
  const std::vector<std::string>& tokens = BoxLib::Tokenize(full_path,"/");
  std::string dir = (full_path[0] == '/' ? "/" : "");
  for (int i=0; i<tokens.size()-1; ++i) {
    dir += tokens[i];
    if (i<tokens.size()-2) dir += "/";
  }
  
  if(!BoxLib::FileExists(dir)) {
    if ( ! BoxLib::UtilCreateDirectory(dir, 0755)) {
      BoxLib::CreateDirectoryFailed(dir);
    }
  }
}


void
Structured_observations(const PArray<Observation>& observation_array,
      Amanzi::ObservationData& observation_data)
{
  for (int i=0; i<observation_array.size(); ++i)
    {
      std::string label = AMR_to_Amanzi_label_map[observation_array[i].name];
      int ntimes = observation_array[i].times.size();
      std::vector<Amanzi::ObservationData::DataTriple> dt(ntimes);
      const std::map<int, Real> vals = observation_array[i].vals;
      for (std::map<int, Real>::const_iterator it = vals.begin();
           it != vals.end(); ++it) {
        int j = it->first;
        dt[j].value = it->second;
        dt[j].time = observation_array[i].times[j];
        dt[j].is_valid = true;
      }
      observation_data[label] = dt;
    }  
}  

std::pair<bool,std::string>
InputVersionOK(const std::string& version)
{
  bool returnVal = true;
  std::string returnStr;

  std::vector<std::string> version_tokens = BoxLib::Tokenize(version,".");

  if ( version_tokens.size() != 3 ) {
    returnVal = false;
    returnStr = "Input file version must be of the form X.X.X";
  }
  else {
    if (atoi(version_tokens[0].c_str()) < amanzi_min_major_version ) {
      returnVal = false;
      std::stringstream str;
      str << "Input file format must have major version >= ";
      str << amanzi_min_major_version;
      returnStr = str.str();
    }
    else if (atoi(version_tokens[1].c_str()) < amanzi_min_minor_version ) {
      returnVal = false;
      std::stringstream str;
      str << "Input file format must have minor version >= ";
      str << amanzi_min_minor_version;
      returnStr = str.str();
    }
    else if (atoi(version_tokens[2].c_str()) < amanzi_min_patch_version ) {
      returnVal = false;
      std::stringstream str;
      str << "Input file format must have patch version >= ";
      str << amanzi_min_patch_version;
      returnStr = str.str();
    }
  }
  return std::pair<bool,std::string>(returnVal,returnStr);
}

static
bool ConfirmFileExists(const std::string& name) {
  std::ifstream ifs;  
  ifs.open(name.c_str(), std::ios::in);
  bool ok = ifs.good();
  ifs.close();
  return ok;
}

AmanziStructuredGridSimulationDriver::AmanziStructuredGridSimulationDriver(const std::string& input_file):
  input_file_(input_file)
{
}

Amanzi::Simulator::ReturnType
AmanziStructuredGridSimulationDriver::Run(const MPI_Comm&               mpi_comm,
                                          Amanzi::ObservationData&      output_observations)
{
    int argc=0;
    char** argv;

    // Let's start things up.
    BoxLib::Initialize(argc,argv,false,mpi_comm);

    // Parse our input file.
    if (input_file_.find(".xml") != std::string::npos)
    {
      Amanzi::AmanziInput::InputConverterS converter;
      bool valid;
      converter.Init(input_file_, valid);
      if (!valid)
      {
        std::cout << "ERROR: The xml input file " << input_file_ << " that was specified using the command line option --xml_file is not valid." << std::endl;
        throw std::string("amanzi not run");
      }

      // Translate from XML -> ParmParse.
      ParmParse::Initialize(argc,argv,NULL);
      converter.Translate();
    }
    else
    {
      // Try opening it as a ParmParse file.
      ParmParse::Initialize(argc,argv,input_file_.c_str());
    }

    ParmParse pp;

    bool petsc_file_exists = false;
    bool petsc_file_specified = false;

#ifdef BL_USE_PETSC
    std::string petsc_help = "Amanzi-S passthrough access to PETSc help option\n";
    std::string petsc_file_str = "Petsc_Options_File";
    std::string petsc_options_file;

    petsc_file_specified = pp.contains(petsc_file_str.c_str());
    if (petsc_file_specified) {
      pp.get(petsc_file_str.c_str(), petsc_options_file);
      petsc_file_exists = ConfirmFileExists(petsc_options_file);
      if (petsc_file_exists) {
        PetscInitialize(&argc,&argv,petsc_options_file.c_str(),petsc_help.c_str());
      }
      else
      {
        PetscInitializeNoArguments();
      }
    }
    else {
      PetscInitializeNoArguments();
    }
#endif

    // Now its ok to write message
    if (petsc_file_specified) {
      if (petsc_file_exists) {
        if (ParallelDescriptor::IOProcessor()) {
          std::cout << "Initializing PETSc with parameter file: \""
                    << petsc_options_file << "\"" << std::endl;
        }
      }
      else {
        std::cout << "\nWARNING: Couldn't open PETSc parameter file: \""
                  << petsc_options_file << "\" ... continuing anyway\n" << std::endl;
      }
    }

    BL_PROFILE_VAR("main()", pmain);

    // Check version number
    std::string amanzi_version_str = "Amanzi_Input_Format_Version";
    if (!pp.contains(amanzi_version_str.c_str()))
    {
      std::string str = "Must specify a value for the top-level parameter: \"" 
        + amanzi_version_str + "\"";
      BoxLib::Abort(str.c_str());
    }

    std::string version;
    pp.get(amanzi_version_str.c_str(), version);
    std::pair<bool,std::string> status = InputVersionOK(version);
    if (!status.first) {
      if (ParallelDescriptor::IOProcessor()) {
        std::cout << status.second << std::endl;
      }
      ParallelDescriptor::Barrier();
      BoxLib::Abort();
    }

    bool pause_for_debug = false;
    if (pp.contains("Pause_for_Debug"))
      pp.get("Pause_For_Debug", pause_for_debug);
    if ( pause_for_debug && ParallelDescriptor::IOProcessor() ) {
        std::string junk;
        std::cout << "Waiting to attach debugger.  Enter any string to continue ";
        std::cin >> junk;
    }
    ParallelDescriptor::Barrier();

    if ( pause_for_debug && ParallelDescriptor::IOProcessor() ) {
        std::cout << "   continuing run..." << std::endl;
    }

    std::string echo_file;
    if (pp.contains("Echo_Input"))
      pp.get("Echo_Input", echo_file);
    if (!echo_file.empty() && ParallelDescriptor::IOProcessor()) {
      EnsureFolderExists(echo_file);
      std::ofstream ofs; ofs.open(echo_file.c_str());
      bool prettyPrint = false;
      ParmParse::dumpTable(ofs,prettyPrint);
      ofs.close();
      ParallelDescriptor::Barrier();
    }

    const Real run_strt = ParallelDescriptor::second();

    int  max_step;
    Real strt_time;
    Real stop_time;

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

    PMAmr* amrptr = new PMAmr;

    amrptr->init(strt_time,stop_time);
    
    // If we set the regrid_on_restart flag and if we are *not* going to take
    //    a time step then we want to go ahead and regrid here.
    if ( amrptr->RegridOnRestart() && 
         ( (amrptr->levelSteps(0) >= max_step) ||
           (amrptr->cumTime() >= stop_time) ) )
    {
        //
        // Regrid only!
        //
        amrptr->RegridOnly(amrptr->cumTime());
    }
    
    while ( amrptr->okToContinue() )
    {
        amrptr->coarseTimeStep(stop_time);
    }

    // Write final checkpoint and plotfile
    if (amrptr->stepOfLastCheckPoint() < amrptr->levelSteps(0)) {
        amrptr->checkPoint();
    }
    amrptr->LinkFinalCheckpoint(amrptr->stepOfLastCheckPoint());

    if (amrptr->stepOfLastPlotFile() < amrptr->levelSteps(0)) {
        amrptr->writePlotFile();
    }

    // Process the observations
    const PArray<Observation>& observation_array = amrptr->TheObservations();

    Structured_observations(observation_array,output_observations);

    delete amrptr;

    const int IOProc   = ParallelDescriptor::IOProcessorNumber();
    Real      run_stop = ParallelDescriptor::second() - run_strt;

    ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

    if (ParallelDescriptor::IOProcessor())
      {
        std::cout << "Run time = " << run_stop << std::endl;
        std::cout << "SCOMPLETED\n";
      }

    BL_PROFILE_VAR_STOP(pmain);

    BoxLib::Finalize(false); // Calling routine responsible for MPI_Finalize call
#ifdef BL_USE_PETSC
    PetscFinalize();
#endif

    return Amanzi::Simulator::SUCCESS;
}
