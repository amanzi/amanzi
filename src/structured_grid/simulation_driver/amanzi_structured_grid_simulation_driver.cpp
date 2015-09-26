#include "InputTranslator.hh"
#include "amanzi_structured_grid_simulation_driver.H"
#include "ParmParse.H"
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


// v1 spec constructor -- delete when we get rid of v1.2 spec.
AmanziStructuredGridSimulationDriver::AmanziStructuredGridSimulationDriver(const std::string& xmlInFileName)
{
  std::string spec;
  Teuchos::ParameterList driver_parameter_list = Amanzi::AmanziNewInput::translate(xmlInFileName, spec);
  plist_ = new Teuchos::ParameterList(Amanzi::AmanziNewInput::translate(xmlInFileName, spec));
}

AmanziStructuredGridSimulationDriver::AmanziStructuredGridSimulationDriver(xercesc::DOMDocument* input)
{
  plist_ = NULL;
//	ParmParse::Initialize(argc,argv,NULL);
//
//  Amanzi::AmanziInput::InputConverterS converter(input);
//  plist_ = new Teuchos::ParameterList(converter.Translate());
// NOT IN DEFAULT YET.
}

AmanziStructuredGridSimulationDriver::~AmanziStructuredGridSimulationDriver()
{
  if (plist_ != NULL)
    delete plist_;
}

Amanzi::Simulator::ReturnType
AmanziStructuredGridSimulationDriver::Run (const MPI_Comm&               mpi_comm,
                                           Amanzi::ObservationData&      output_observations)
{
    int argc=0;
    char** argv;

    // NOTE: Delay checking input version until we have started parallel so that we can guarantee
    //       a single error message

    bool petsc_file_exists = false;
    bool petsc_file_specified = false;

#ifdef BL_USE_PETSC
    std::string petsc_help = "Amanzi-S passthrough access to PETSc help option\n";
    std::string petsc_file_str = "Petsc Options File";
    std::string petsc_options_file;

    petsc_file_specified = plist_->isParameter(petsc_file_str);
    if (petsc_file_specified) {
      petsc_options_file = Teuchos::getParameter<std::string>(*plist_, petsc_file_str);
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

    BoxLib::Initialize(argc,argv,false,mpi_comm);

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
    std::string amanzi_version_str = "Amanzi Input Format Version";
    if (!plist_->isParameter(amanzi_version_str))
    {
      std::string str = "Must specify a value for the top-level parameter: \"" 
        + amanzi_version_str + "\"";
      BoxLib::Abort(str.c_str());
    }

    std::string version = Teuchos::getParameter<std::string>(*plist_, amanzi_version_str);
    std::pair<bool,std::string> status = InputVersionOK(version);
    if (!status.first) {
      if (ParallelDescriptor::IOProcessor()) {
        std::cout << status.second << std::endl;
      }
      ParallelDescriptor::Barrier();
      BoxLib::Abort();
    }

    bool pause_for_debug = false;
    if (plist_->isParameter("Pause For Debug"))
      {
          pause_for_debug= Teuchos::getParameter<bool>(*plist_, "Pause For Debug");
      }

    if ( pause_for_debug && ParallelDescriptor::IOProcessor() ) {
        std::string junk;
        std::cout << "Waiting to attach debugger.  Enter any string to continue ";
        std::cin >> junk;
    }
    ParallelDescriptor::Barrier();

    if ( pause_for_debug && ParallelDescriptor::IOProcessor() ) {
        std::cout << "   continuing run..." << std::endl;
    }

    if (plist_->isParameter("PPfile"))
      {
	const std::string& PPfile = Teuchos::getParameter<std::string>(*plist_, "PPfile");
	ParmParse::Initialize(argc,argv,PPfile.c_str());
      }

    // Determine whether we need to convert to input file to 
    //native structured format
    bool native = plist_->get<bool>("Native Structured Input",false);
    Teuchos::ParameterList converted_parameter_list;
    if (!native) 
      converted_parameter_list =
	Amanzi::AmanziInput::convert_to_structured(*plist_);
    else
      converted_parameter_list = *plist_;

    if (plist_->isParameter("EchoXMLfile"))
      {
        const std::string& EchoXMLfile = Teuchos::getParameter<std::string>(*plist_, "EchoXMLfile");
        Teuchos::writeParameterListToXmlFile(converted_parameter_list,EchoXMLfile);
      }

    // Stuff away a static copy of the input parameters
    PorousMedia::SetInputParameterList(converted_parameter_list);

    BoxLib::Initialize_ParmParse(converted_parameter_list);

    const Teuchos::ParameterList& echo_list = plist_->sublist("Echo Translated Input");
    if (echo_list.isParameter("Format")) {
      if (echo_list.get<std::string>("Format") == "native") {
	if (ParallelDescriptor::IOProcessor()) {
	  const std::string& pp_file = echo_list.get<std::string>("File Name");
	  EnsureFolderExists(pp_file);
	  std::ofstream ofs; ofs.open(pp_file.c_str());
	  bool prettyPrint = false;
	  ParmParse::dumpTable(ofs,prettyPrint);
	  ofs.close();
	}
      }
      ParallelDescriptor::Barrier();
    }

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
