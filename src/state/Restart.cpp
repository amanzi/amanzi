#include "Restart.hpp"
#include "Epetra_MpiComm.h"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include <iostream>
#include <iomanip>

#include "boost/filesystem.hpp"

Amanzi::Restart::Restart (Teuchos::ParameterList& plist_, Epetra_MpiComm* comm_):
  plist(plist_), disabled(false), comm(comm_), restart_output(NULL)
{
  read_parameters(plist);

  // set the line prefix for output
  this->setLinePrefix("Amanzi::Restart     ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);
  
  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&plist,this);
  restart_output = new Amanzi::HDF5_MPI(*comm);
  restart_output->setTrackXdmf(false);
  
}

// this constructor makes an object that will not create any output
Amanzi::Restart::Restart (): disabled(true), restart_output(NULL)
{
}

Amanzi::Restart::~Restart() 
{
  if (!restart_output) delete restart_output;
}


void Amanzi::Restart::read_parameters(Teuchos::ParameterList& plist)
{
  filebasename = plist.get<string>("File Name Base","amanzi_restart");
  filenamedigits = plist.get<int>("File Name Digits",5);

  boost::filesystem::path fbn(filebasename);
  boost::filesystem::path parent = fbn.parent_path();
  if (!parent.empty()) {
    // we need to check whether the parent path exists
    if (!boost::filesystem::exists(parent)) {
      Errors::Message m("The file path '"+parent.string()+"' used for checkpoints does not exist.");
      Exceptions::amanzi_throw(m);    
    }
    if (!boost::filesystem::is_directory(parent)) {
      Errors::Message m("The file path '"+parent.string()+"' used for checkpoints is not a directory.");
      Exceptions::amanzi_throw(m);    
    }
    if ( (boost::filesystem::status(parent).permissions() & boost::filesystem::owner_write) == 0) {
      Errors::Message m("The directory '"+parent.string()+"' used for checkpoints is not writeable.");
      Exceptions::amanzi_throw(m);       
    }
  }  

  number_of_cycle_intervals = 0;

  // read the cycle data namelist

  if ( plist.isSublist("Cycle Data") ) 
    {
      Teuchos::ParameterList &iclist = plist.sublist("Cycle Data");
      
      interval = iclist.get<int>("Interval",1);
      start = iclist.get<int>("Start",0);
      end = iclist.get<int>("End",-1);
      
      if (iclist.isParameter("Steps"))
	{
	  steps = iclist.get<Teuchos::Array<int> >("Steps");  
	}
    }	
  else
    {
      // error
    }

}


void Amanzi::Restart::dump_state(State& S, bool force)
{
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  if (!disabled) {
    if (force || dump_requested(S.get_cycle())) {
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM,true)) {
	*out << "Writing checkpoint, cycle = " << S.get_cycle() << std::endl;
      }
      
      // create the restart file
      std::stringstream oss; 
      oss.flush();
      oss << filebasename;
      oss.fill('0'); 
      oss.width(filenamedigits);
      oss << std::right << S.get_cycle(); 
      
      restart_output->createDataFile(oss.str());	  
      
      
      // dump all the state vectors into the restart file
      restart_output->writeDataReal(*S.get_pressure(),"pressure");
      restart_output->writeDataReal(*S.get_lambda(),"lambda");
      restart_output->writeDataReal(*S.get_porosity(),"porosity");
      restart_output->writeDataReal(*S.get_particle_density(),"particle density");
      restart_output->writeDataReal(*S.get_water_saturation(),"water saturation");
      restart_output->writeDataReal(*S.get_prev_water_saturation(),"previous water saturation");
      restart_output->writeDataReal(*S.get_water_density(),"water density");
      restart_output->writeDataReal(*S.get_horizontal_permeability(),"horizontal permeability");
      restart_output->writeDataReal(*S.get_vertical_permeability(),"vertical permeability");
      restart_output->writeDataReal(*S.get_material_ids(),"material IDs");
      
      for (int i=0; i<S.get_number_of_components(); i++) {
	std::stringstream tcc_name;
	
	tcc_name << "component " << i;
	
	restart_output->writeDataReal(*(*S.get_total_component_concentration())(i),tcc_name.str());
      }
      
      restart_output->writeDataReal(*S.get_darcy_flux(),"darcy flux");	  
      
      
      restart_output->writeAttrReal(S.get_time(),"time");
      restart_output->writeAttrReal(S.get_last_time(),"last time");
      restart_output->writeAttrInt(S.get_cycle(),"cycle");
      
      restart_output->writeAttrReal((*S.get_gravity())[0],"gravity x");
      restart_output->writeAttrReal((*S.get_gravity())[1],"gravity y");
      restart_output->writeAttrReal((*S.get_gravity())[2],"gravity z");
      
      restart_output->writeAttrReal((*S.get_density()),"density");
      restart_output->writeAttrReal((*S.get_viscosity()),"viscosity");	
      
      restart_output->writeAttrInt(S.get_number_of_components(),"number of components");
      restart_output->writeAttrInt(S.number_of_minerals(),"Number of minerals");
      restart_output->writeAttrInt(S.number_of_ion_exchange_sites(),"Number of ion exchange sites");
      restart_output->writeAttrInt(S.number_of_sorption_sites(),"Number of sorption sites");
      restart_output->writeAttrInt(static_cast<int>(S.use_sorption_isotherms()),
                                   "Use Sorption Isotherms");

      for (int i=0; i<S.get_number_of_components(); i++) {
	std::stringstream name;
	name << "free ion concentrations " << i;
	restart_output->writeDataReal(*(*S.free_ion_concentrations())(i),name.str());
      }
      for (int i=0; i<S.get_number_of_components(); i++) {
	std::stringstream name;
	name << "primary activity coeff " << i;
	restart_output->writeDataReal(*(*S.primary_activity_coeff())(i),name.str());
      }

      if ( ! S.secondary_activity_coeff().is_null() ) {
	for (int i=0; i<S.secondary_activity_coeff()->NumVectors(); i++) {
	  std::stringstream name;
	  name << "secondary activity coeff " << i;
	  restart_output->writeDataReal(*(*S.secondary_activity_coeff())(i),name.str());
	}
      }

      for (int m = 0; m < S.number_of_minerals(); ++m) {
        std::stringstream name;
        name << "mineral volume fractions " << m;
        restart_output->writeDataReal(*(*S.mineral_volume_fractions())(m), name.str());
      }
      for (int m = 0; m < S.number_of_minerals(); ++m) {
	std::stringstream name;
	name << "mineral specific surface area " << m;
        restart_output->writeDataReal(*(*S.mineral_specific_surface_area())(m), name.str());
      }
      if (S.using_sorption()) {
	for (int i=0; i<S.get_number_of_components(); i++) {
	  std::stringstream name;
	  name << "total sorbed " << i;
	  restart_output->writeDataReal(*(*S.total_sorbed())(i),name.str());
	}      
      }
      if (S.use_sorption_isotherms()) {
	for (int i=0; i<S.number_of_sorption_sites(); i++) {
	  std::stringstream name;
	  name << "sorption sites " << i;
	  restart_output->writeDataReal(*(*S.sorption_sites())(i),name.str());
	}
      
	for (int i=0; i<S.number_of_sorption_sites(); i++) {
	  std::stringstream name;
	  name << "surface complex free site conc " << i;
	  restart_output->writeDataReal(*(*S.surface_complex_free_site_conc())(i),name.str());
	}
      }
      for (int i=0; i<S.number_of_ion_exchange_sites(); i++) {
	std::stringstream name;
	name << "ion exchange sites " << i;
	restart_output->writeDataReal(*(*S.ion_exchange_sites())(i),name.str());
      }
      for (int i=0; i<S.number_of_ion_exchange_sites(); i++) {
	std::stringstream name;
	name << "ion exchange ref cation conc " << i;
	restart_output->writeDataReal(*(*S.ion_exchange_ref_cation_conc())(i),name.str());
      }
      if (S.use_sorption_isotherms()) {
	for (int i=0; i<S.get_number_of_components(); i++) {
	  std::stringstream name;
	  name << "isotherm kd " << i;
	  restart_output->writeDataReal(*(*S.isotherm_kd())(i),name.str());
	}      
	for (int i=0; i<S.get_number_of_components(); i++) {
	  std::stringstream name;
	  name << "isotherm freundlich n " << i;
	  restart_output->writeDataReal(*(*S.isotherm_freundlich_n())(i),name.str());
	}
	for (int i=0; i<S.get_number_of_components(); i++) {
	  std::stringstream name;
	  name << "isotherm langmuir b " << i;
	  restart_output->writeDataReal(*(*S.isotherm_langmuir_b())(i),name.str());
	}
      }
     
      restart_output->writeDataReal(*S.get_specific_storage(),"specific storage");
      restart_output->writeDataReal(*S.get_specific_yield(),"specific yield");
      
      
    }
  }
}


void Amanzi::Restart::read_state(State& S, std::string& filename)
{
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab
  
  

  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM,true))	  
    {
      *out << "Reading checkpoint from file " << filename << std::endl;
    }
  
  Amanzi::HDF5_MPI *restart_input = new Amanzi::HDF5_MPI(*comm, filename); 

  int idummy;
  double dummy;

  // first we must read the number of components
  restart_input->readAttrInt(idummy,"number of components");
  S.set_number_of_components(idummy);
  restart_input->readAttrInt(idummy,"Number of minerals");
  S.set_number_of_minerals(idummy);

  // read the attributes
  restart_input->readAttrReal(dummy,"last time");
  S.set_time(dummy);
  restart_input->readAttrReal(dummy,"time");
  S.set_time(dummy);  

  restart_input->readAttrInt(idummy,"cycle");
  S.set_cycle(idummy);

  double g[3];
  restart_input->readAttrReal(g[0],"gravity x");
  restart_input->readAttrReal(g[1],"gravity y");
  restart_input->readAttrReal(g[2],"gravity z");
  S.set_gravity(g);
  
  restart_input->readAttrReal(dummy,"viscosity");
  S.set_viscosity(dummy);

  restart_input->readAttrReal(dummy,"density");
  S.set_water_density(dummy);
  
  // read the vectors
  Epetra_Vector* face_vector = new Epetra_Vector(S.get_mesh().face_epetra_map(false));
  restart_input->readData(*face_vector,"darcy flux");
  S.set_darcy_flux(*face_vector);
  delete face_vector;

  Epetra_Vector* cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_input->readData(*cell_vector,"water saturation");
  S.set_water_saturation(*cell_vector);
  delete cell_vector;
  
  cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_input->readData(*cell_vector,"previous water saturation");
  S.set_prev_water_saturation(*cell_vector);
  delete cell_vector;  

  cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_input->readData(*cell_vector,"water density");
  S.set_water_density(*cell_vector);
  delete cell_vector;

  cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_input->readData(*cell_vector,"pressure");
  S.set_pressure(*cell_vector);
  delete cell_vector;

  face_vector = new Epetra_Vector(S.get_mesh().face_epetra_map(false));
  restart_input->readData(*face_vector,"lambda");
  S.set_lambda(*face_vector);
  delete face_vector;  
  
  cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_input->readData(*cell_vector,"porosity");
  S.set_porosity(*cell_vector);
  delete cell_vector;

  cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_input->readData(*cell_vector,"particle density");
  S.set_particle_density(*cell_vector);
  delete cell_vector;

  cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_input->readData(*cell_vector,"horizontal permeability");
  S.set_horizontal_permeability(*cell_vector);
  delete cell_vector;

  cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_input->readData(*cell_vector,"vertical permeability");
  S.set_vertical_permeability(*cell_vector);
  delete cell_vector;

  cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_input->readData(*cell_vector,"material IDs");
  S.set_material_ids(*cell_vector);
  delete cell_vector;

  Epetra_MultiVector* cell_multivector  = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), S.get_number_of_components());
  for (int i=0; i< S.get_number_of_components(); i++)
    {
      std::stringstream tcc_name;
      tcc_name << "component " << i;

      restart_input->readData(*(*cell_multivector)(i),tcc_name.str()); 
    }
  S.set_total_component_concentration(*cell_multivector);
  delete cell_multivector; 

  


  cell_multivector  = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), 
					     S.get_number_of_components());
  for (int i=0; i< S.get_number_of_components(); i++) {
    std::stringstream name;
    name << "free ion concentrations " << i;
    restart_input->readData(*(*cell_multivector)(i),name.str()); 
  }
  S.set_free_ion_concentrations(*cell_multivector);
  delete cell_multivector;   

  cell_multivector  = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), 
					     S.get_number_of_components());
  for (int i=0; i< S.get_number_of_components(); i++) {
    std::stringstream name;
    name << "primary activity coeff " << i;
    restart_input->readData(*(*cell_multivector)(i),name.str()); 
  }
  S.set_primary_activity_coeff(*cell_multivector);
  delete cell_multivector;

  if ( ! S.secondary_activity_coeff().is_null() ) {
    cell_multivector  = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), 
					       S.secondary_activity_coeff()->NumVectors());
    for (int i=0; i< S.secondary_activity_coeff()->NumVectors(); i++) {
      std::stringstream name;
      name << "secondary activity coeff " << i;
      restart_input->readData(*(*cell_multivector)(i),name.str()); 
    }
    S.set_secondary_activity_coeff(*cell_multivector);
    delete cell_multivector;
  }

  if (S.number_of_minerals() > 0) {
    
    cell_multivector  = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), 
					       S.number_of_minerals());
    for (int i=0; i< S.number_of_minerals(); i++) {
      std::stringstream name;
      name << "mineral volume fractions " << i;
      restart_input->readData(*(*cell_multivector)(i),name.str()); 
    }
    S.set_mineral_volume_fractions(*cell_multivector);
    delete cell_multivector;  

    cell_multivector  = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), 
					       S.number_of_minerals());
    for (int i=0; i< S.number_of_minerals(); i++) {
      std::stringstream name;
      name << "mineral specific surface area " << i;
      restart_input->readData(*(*cell_multivector)(i),name.str()); 
    }
    S.set_mineral_specific_surface_area(*cell_multivector);
    delete cell_multivector;      
  }

  if (S.using_sorption()) {
    cell_multivector  = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), 
					       S.get_number_of_components());
    for (int i=0; i< S.get_number_of_components(); i++) {
      std::stringstream name;
      name << "total sorbed " << i;
      restart_input->readData(*(*cell_multivector)(i),name.str()); 
    }
    S.set_total_sorbed(*cell_multivector);
    delete cell_multivector;      
  }

  if (S.use_sorption_isotherms()) {
    cell_multivector  = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), 
					       S.number_of_sorption_sites());
    for (int i=0; i< S.number_of_sorption_sites(); i++) {
      std::stringstream name;
      name << "sorption sites " << i;
      restart_input->readData(*(*cell_multivector)(i),name.str()); 
    }
    S.set_sorption_sites(*cell_multivector);
    delete cell_multivector;      

    cell_multivector  = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), 
					       S.number_of_sorption_sites());
    for (int i=0; i< S.number_of_sorption_sites(); i++) {
      std::stringstream name;
      name << "surface complex free site conc " << i;
      restart_input->readData(*(*cell_multivector)(i),name.str()); 
    }
    S.set_surface_complex_free_site_conc(*cell_multivector);
    delete cell_multivector;
  }  

  if (S.number_of_ion_exchange_sites() > 0) {
    cell_multivector  = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), 
					       S.number_of_ion_exchange_sites());
    for (int i=0; i< S.number_of_ion_exchange_sites(); i++) {
      std::stringstream name;
      name << "ion exchange sites " << i;
      restart_input->readData(*(*cell_multivector)(i),name.str());
    }
    S.set_ion_exchange_sites(*cell_multivector);
    delete cell_multivector;    
    
    cell_multivector  = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), 
					       S.number_of_ion_exchange_sites());
    for (int i=0; i< S.number_of_ion_exchange_sites(); i++) {
      std::stringstream name;
      name << "ion exchange ref cation conc " << i;
      restart_input->readData(*(*cell_multivector)(i),name.str());
    }
    S.set_ion_exchange_ref_cation_conc(*cell_multivector);
    delete cell_multivector;      
  }
  
  if (S.use_sorption_isotherms()) {
    cell_multivector  = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), 
					       S.get_number_of_components());
    for (int i=0; i< S.get_number_of_components(); i++) {
      std::stringstream name;
      name << "isotherm kd " << i;
      restart_input->readData(*(*cell_multivector)(i),name.str());
    }
    S.set_isotherm_kd(*cell_multivector);
    delete cell_multivector;       

    cell_multivector  = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), 
					       S.get_number_of_components());
    for (int i=0; i< S.get_number_of_components(); i++) {
      std::stringstream name;
      name << "isotherm freundlich n " << i;
      restart_input->readData(*(*cell_multivector)(i),name.str());
    }
    S.set_isotherm_freundlich_n(*cell_multivector);
    delete cell_multivector;       

    cell_multivector  = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), 
					       S.get_number_of_components());
    for (int i=0; i< S.get_number_of_components(); i++) {
      std::stringstream name;
      name << "isotherm langmuir b " << i;
      restart_input->readData(*(*cell_multivector)(i),name.str());
    }
    S.set_isotherm_langmuir_b(*cell_multivector);
    delete cell_multivector;       
  }

  cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_input->readData(*cell_vector,"specific storage");
  S.set_specific_storage(*cell_vector);
  delete cell_vector;    

  cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_input->readData(*cell_vector,"specific yield");
  S.set_specific_yield(*cell_vector);
  delete cell_vector;    
 

  delete restart_input;

}






bool Amanzi::Restart::dump_requested(int cycle)
{
  if (steps.size() > 0) 
    {
      for (int i=0; i<steps.size(); i++) 
	{
	  if (cycle == steps[i])
	    {
	      return true;
	    }
	}
    }
  else if ( (end<0) || (cycle<=end) ) 
    {
      if (start<=cycle)  
	{
	  int cycle_loc = cycle - start;
	  
	  if (cycle_loc % interval == 0) 
	    {
	      return true;
	    }
	  
	}
    }
      
  // if none of the conditions apply we do not dump
  return false;

}


