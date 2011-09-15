#include "Vis.hpp"
#include "Epetra_MpiComm.h"

Amanzi::Vis::Vis (Teuchos::ParameterList& plist_, const Epetra_MpiComm& comm):
  plist(plist_), disabled(false)
{
  read_parameters(plist);

  viz_output = new Amanzi::HDF5_MPI(comm);
  viz_output->setTrackXdmf(true);
  
}

// this constructor makes an object that will not create any output
Amanzi::Vis::Vis (): disabled(true)
{
}

void Amanzi::Vis::set_compnames (std::vector<string>& compnames_)
{
  compnames = compnames_;
}


void Amanzi::Vis::create_files(Amanzi::AmanziMesh::Mesh& mesh)
{
  if (!disabled)
    {

      // create file name for the mesh 
      std::stringstream meshfilename;  
      meshfilename << filebasename << "_mesh";
      
      // create file name for the data
      std::stringstream datafilename;  
      datafilename << filebasename << "_data";  
      
      viz_output->createMeshFile(mesh, meshfilename.str());
      viz_output->createDataFile(datafilename.str());
    }
}

void Amanzi::Vis::read_parameters(Teuchos::ParameterList& plist)
{
  filebasename = plist.get<string>("file base name","amanzi_vis");
  
  number_of_cycle_intervals = 0;
  number_of_time_intervals  = 0;

  // read the interval namelists
  for (Teuchos::ParameterList::ConstIterator i = plist.begin(); i != plist.end(); i++)
    {
      // we assume that all sublists will contain one interval sublist
      // and we do not care about the name of the sublist
      if (plist.isSublist(plist.name(i)))
	{      
	  const Teuchos::ParameterList &ilist = plist.sublist( plist.name(i) );

	  if ( ilist.isSublist("time range") )
	    {
	      const Teuchos::ParameterList &itlist = ilist.sublist("time range");
	      
	      time_freq.push_back(itlist.get<double>("time frequency"));
	      time_start.push_back(itlist.get<double>("start time"));
	      time_end.push_back(itlist.get<double>("end time"));
	      
	      number_of_time_intervals++;	      
	    }
	  else if ( ilist.isSublist("cycle range") ) 
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
};


Amanzi::Vis::~Vis () 
{
  delete viz_output;
}


void Amanzi::Vis::dump_state(double prev_time, double time, int cycle,
		       State& S)
{
  
  if (!disabled) 
    {
      if (dump_requested(prev_time, time, cycle))
	{
	  viz_output->createTimestep(time,cycle);
	  
	  // dump all the state vectors into the vis file
	  viz_output->writeCellDataReal(*S.get_pressure(),"pressure");
	  viz_output->writeCellDataReal(*S.get_porosity(),"porosity");
	  viz_output->writeCellDataReal(*S.get_water_saturation(),"water saturation");
	  viz_output->writeCellDataReal(*S.get_water_density(),"water density");
	  viz_output->writeCellDataReal(*S.get_permeability(),"permeability");
	  viz_output->writeCellDataReal(*(*S.get_darcy_velocity())(0),"darcy velocity x");
	  viz_output->writeCellDataReal(*(*S.get_darcy_velocity())(1),"darcy velocity y");
	  viz_output->writeCellDataReal(*(*S.get_darcy_velocity())(2),"darcy velocity z");
	  
	  for (int i=0; i<S.get_number_of_components(); i++)
	    {
	      std::stringstream tcc_name;
	      tcc_name << "component " << i;
	      
	      viz_output->writeCellDataReal(*(*S.get_total_component_concentration())(i),tcc_name.str());
	    }
	  
	  viz_output->endTimestep();
	  
	}
    }
}


bool Amanzi::Vis::dump_requested(double prev_time, double time, int cycle)
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

   // time intervals 
   for (int i=0; i<number_of_time_intervals; i++)
    {
      if ( (time_start[i]<=time) && (time <= time_end[i]) )
	{
	  double time_loc = time - time_start[i];
	  double prev_time_loc = prev_time - time_start[i];
	  
	  if (  floor(time_loc/time_freq[i]) > floor(prev_time_loc/time_freq[i]) )
	    {
	      return true;
	    }
	}
	
     }

  // if none of the conditions apply we do not dump
  return false;

}
