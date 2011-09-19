#include "Restart.hpp"
#include "Epetra_MpiComm.h"

Amanzi::Restart::Restart (Teuchos::ParameterList& plist_, Epetra_MpiComm* comm_):
  plist(plist_), disabled(false), comm(comm_)
{
  read_parameters(plist);

  restart_output = new Amanzi::HDF5_MPI(*comm);
  restart_output->setTrackXdmf(false);
  
}

// this constructor makes an object that will not create any output
Amanzi::Restart::Restart (): disabled(true)
{
}

Amanzi::Restart::~Restart() 
{
  delete restart_output;
}


void Amanzi::Restart::create_files()
{
  if (!disabled)
    {
      // create file name for the data
      std::stringstream datafilename;  
      datafilename << filebasename << "_data";  
      
      restart_output->createDataFile(datafilename.str());
    }
}

void Amanzi::Restart::read_parameters(Teuchos::ParameterList& plist)
{
  filebasename = plist.get<string>("file base name","amanzi_restart");
  
  number_of_cycle_intervals = 0;

  // read the interval namelists
  for (Teuchos::ParameterList::ConstIterator i = plist.begin(); i != plist.end(); i++)
    {
      // we assume that all sublists will contain one interval sublist
      // and we do not care about the name of the sublist
      if (plist.isSublist(plist.name(i)))
	{      
	  const Teuchos::ParameterList &ilist = plist.sublist( plist.name(i) );

	  if ( ilist.isSublist("cycle range") ) 
	    {
	      const Teuchos::ParameterList &iclist = ilist.sublist("cycle range");
	      
	      cycle_freq.push_back(iclist.get<int>("cycle frequency"));
	      cycle_start.push_back(iclist.get<int>("start cycle"));
	      cycle_end.push_back(iclist.get<int>("end cycle"));
	      
	      number_of_cycle_intervals++;
	    }
	  else
	    {
	      // error
	    }
	  
	}
    }
}


void Amanzi::Restart::dump_state(State& S)
{
  
  if (!disabled) 
    {
      if (dump_requested(S.get_cycle()))
	{
	  // dump all the state vectors into the restart file
	  restart_output->writeCellDataReal(*S.get_pressure(),"pressure");
	  restart_output->writeCellDataReal(*S.get_porosity(),"porosity");
	  restart_output->writeCellDataReal(*S.get_water_saturation(),"water saturation");
	  restart_output->writeCellDataReal(*S.get_water_density(),"water density");
	  restart_output->writeCellDataReal(*S.get_permeability(),"permeability");
	  restart_output->writeCellDataReal(*(*S.get_darcy_velocity())(0),"darcy velocity x");
	  restart_output->writeCellDataReal(*(*S.get_darcy_velocity())(1),"darcy velocity y");
	  restart_output->writeCellDataReal(*(*S.get_darcy_velocity())(2),"darcy velocity z");
	  
	  for (int i=0; i<S.get_number_of_components(); i++)
	    {
	      std::stringstream tcc_name;
	      
	      tcc_name << "component " << i;
	      
	      restart_output->writeCellDataReal(*(*S.get_total_component_concentration())(i),tcc_name.str());
	    }

	  restart_output->writeCellDataReal(*S.get_darcy_flux(),"darcy flux");	  


	  restart_output->writeAttrReal(S.get_time(),"time");
	  restart_output->writeAttrInt(S.get_cycle(),"cycle");

	  restart_output->writeAttrReal((*S.get_gravity())[0],"gravity x");
	  restart_output->writeAttrReal((*S.get_gravity())[1],"gravity y");
	  restart_output->writeAttrReal((*S.get_gravity())[2],"gravity z");
	  
	  restart_output->writeAttrReal((*S.get_density()),"density");
	  restart_output->writeAttrReal((*S.get_viscosity()),"viscosity");	
	  
	  restart_output->writeAttrInt(S.get_number_of_components(),"number of components");
	}
    }
}


void Amanzi::Restart::read_state(State& S)
{
  
  int idummy;
  double dummy;

  // first we must read the number of components
  restart_output->readAttrInt(idummy,"number of components");
  S.set_number_of_components(idummy);

  // now we can create storage
  S.create_storage();
    
  // read the attributes
  restart_output->readAttrReal(dummy,"time");
  S.set_time(dummy);  

  restart_output->readAttrInt(idummy,"cycle");
  S.set_cycle(idummy);

  double g[3];
  restart_output->readAttrReal(g[0],"gravity x");
  restart_output->readAttrReal(g[1],"gravity y");
  restart_output->readAttrReal(g[2],"gravity z");
  S.set_gravity(g);
  
  restart_output->readAttrReal(dummy,"viscosity");
  S.set_viscosity(dummy);

  restart_output->readAttrReal(dummy,"density");
  S.set_water_density(dummy);
  
  // read the vectors
  Epetra_Vector* face_vector = new Epetra_Vector(S.get_mesh().face_epetra_map(false));
  restart_output->readData(*face_vector,"darcy flux");
  S.set_darcy_flux(*face_vector);
  delete face_vector;

  Epetra_Vector* cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_output->readData(*cell_vector,"water saturation");
  S.set_water_saturation(*cell_vector);
  delete cell_vector;
  
  cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_output->readData(*cell_vector,"water density");
  S.set_water_density(*cell_vector);
  delete cell_vector;

  cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_output->readData(*cell_vector,"pressure");
  S.set_pressure(*cell_vector);
  delete cell_vector;
  
  cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_output->readData(*cell_vector,"porosity");
  S.set_porosity(*cell_vector);
  delete cell_vector;

  cell_vector = new Epetra_Vector(S.get_mesh().cell_epetra_map(false));
  restart_output->readData(*cell_vector,"permeability");
  S.set_permeability(*cell_vector);
  delete cell_vector;

  Epetra_MultiVector* cell_multivector = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), 3);
  restart_output->readData(*(*cell_multivector)(0),"darcy velocity x");
  restart_output->readData(*(*cell_multivector)(1),"darcy velocity y");
  restart_output->readData(*(*cell_multivector)(2),"darcy velocity z");
  S.set_darcy_velocity(*cell_multivector);
  delete cell_multivector;        

  cell_multivector  = new Epetra_MultiVector(S.get_mesh().cell_epetra_map(false), S.get_number_of_components());
  for (int i=0; i< S.get_number_of_components(); i++)
    {
      std::stringstream tcc_name;
      tcc_name << "component " << i;

      restart_output->readData(*(*cell_multivector)(i),tcc_name.str()); 
    }
  S.set_total_component_concentration(*cell_multivector);
  delete cell_multivector;        

}






bool Amanzi::Restart::dump_requested(int cycle)
{
  // cycle intervals
  for (int i=0; i<number_of_cycle_intervals; i++)
    {
      if ( (cycle_start[i]<=cycle) && (cycle<=cycle_end[i]) ) 
	{
	  int cycle_loc = cycle - cycle_start[i];
	  
	  // check if a vis dump is requested
	  if (cycle_loc % cycle_freq[i] == 0) 
	    {
	      return true;
	    }
	  
	}
     }

  // if none of the conditions apply we do not dump
  return false;

}


