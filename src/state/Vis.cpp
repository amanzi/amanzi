#include "Vis.hpp"
#include "Epetra_MpiComm.h"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

Amanzi::Vis::Vis (Teuchos::ParameterList& plist_, Epetra_MpiComm* comm_):
  plist(plist_), disabled(false), comm(comm_)
{
  read_parameters(plist);

  // set the line prefix for output
  this->setLinePrefix("Amanzi::Vis         ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);
  
  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&plist,this);
    
  auxnames.resize(0);
  compnames.resize(0);

  viz_output = new Amanzi::HDF5_MPI(*comm);
  viz_output->setTrackXdmf(true);
  
}

// this constructor makes an object that will not create any output
Amanzi::Vis::Vis (): disabled(true)
{
}

void Amanzi::Vis::set_compnames (std::vector<string>& compnames_ )
{
  compnames = compnames_;
}

void Amanzi::Vis::set_auxnames (std::vector<string>& auxnames_ )
{
  auxnames = auxnames_;
}


void Amanzi::Vis::create_files(Amanzi::AmanziMesh::Mesh& mesh)
{
  if (!disabled)
    {
      
      // create file name for the mesh 
      std::stringstream meshfilename;  
      meshfilename.flush();
      meshfilename << filebasename; 
      meshfilename << "_mesh";

      // create file name for the data
      std::stringstream datafilename;  
      datafilename.flush();
      datafilename << filebasename; 
      datafilename << "_data";
      
      viz_output->createMeshFile(mesh, meshfilename.str());
      viz_output->createDataFile(datafilename.str());
    }
}

void Amanzi::Vis::read_parameters(Teuchos::ParameterList& plist)
{
  filebasename = plist.get<string>("File Name Base","amanzi_vis");
  
  enable_gnuplot = plist.get<bool>("enable gnuplot",false);
  if (comm->NumProc() != 1) enable_gnuplot = false;

  if (plist.isParameter("Start_Period_Stop") ) 
    {
      Teuchos::Array<int> sscf = plist.get<Teuchos::Array<int> >("Start_Period_Stop");
      if (sscf.size()==3)
	{
	  start = sscf[0];
	  end = sscf[2];
	  interval = sscf[1];
	}
      else if (sscf.size()==2) 
	{
	  start = sscf[0];
	  interval = sscf[1];
	  end = -1;
	}
    }
  else if (plist.isParameter("Values") )
    {

    }
}


Amanzi::Vis::~Vis () 
{
  delete viz_output;
}


void Amanzi::Vis::dump_state(State& S, Epetra_MultiVector *auxdata)
{
  if (!disabled) 
    {
      if (dump_requested(S.get_cycle()))
	{
	  
	  using Teuchos::OSTab;
	  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
	  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
	  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab
 	  
	  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM,true))	  
	    {
	      *out << "Writing visualization files, cycle = " << S.get_cycle() << std::endl;
	    }
	  


	  viz_output->createTimestep(S.get_time(),S.get_cycle());
	  
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
	      
	      if (compnames.size() == S.get_number_of_components()) 
		{
		  tcc_name << compnames[i]; 
		}
	      else
		tcc_name << "component " << i;
	      
	      viz_output->writeCellDataReal(*(*S.get_total_component_concentration())(i),tcc_name.str());
	    }
	  
	  // write auxillary data
	  if (auxdata != NULL) 
	    {
	      for (int i=0; i<auxdata->NumVectors(); i++)
		{
		  std::stringstream tcc_name;
		  
		  if (auxnames.size() == auxdata->NumVectors()) 
		    {
		      tcc_name << auxnames[i]; 
		    }
		  else
		    tcc_name << "aux " << i;
		  
		  viz_output->writeCellDataReal( *(*auxdata)(i) ,tcc_name.str());
		}	      

	    }



	  viz_output->endTimestep();
	  

	  // now do a gnuplot dump, if gnuplot output is enabled
	  // but only if we are runnung on one processor
	  if (comm->NumProc() == 1  && enable_gnuplot)
	    {
	      write_gnuplot(S.get_cycle(), S, auxdata);
	    }
	}
    }
}


bool Amanzi::Vis::dump_requested(int cycle)
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



void Amanzi::Vis::write_gnuplot(int cycle, State& S, Epetra_MultiVector* auxdata)
{
  // write each variable to a separate file
  // only works on one PE, not on a truly parallel run

  if (cycle==0) 
    {
      std::stringstream fname_p;
      fname_p << "pressure.dat";
      
      std::filebuf fb;
      fb.open (fname_p.str().c_str(), std::ios::out);
      ostream os_p(&fb);
      
      for (int i=0; i< (S.get_pressure())->MyLength(); i++) 
	os_p << (*S.get_pressure())[i] << endl;    
      
      fb.close();
      
      std::stringstream fname_s;
      fname_s << "saturation.dat";
      
      fb.open (fname_s.str().c_str(), std::ios::out);
      ostream os_s(&fb);
      
      for (int i=0; i< (S.get_water_saturation())->MyLength(); i++) 
	os_s << (*S.get_water_saturation())[i] << endl;    
      
      fb.close();
    }

  for (int nc=0; nc<S.get_number_of_components(); nc++) 
    {
      std::stringstream fname;
      fname << "conc_" << nc;
      
      if (compnames.size() == S.get_number_of_components())
	{
	  fname << "_" << compnames[nc];
	}

      fname << "_" << std::setfill('0') << std::setw(5) << cycle << ".dat";
      
      std::filebuf fb;
      fb.open (fname.str().c_str(), std::ios::out);
      ostream os(&fb);

      // now dump the Epetra Vector
      for (int i=0; i< (S.get_total_component_concentration())->MyLength(); i++) 
	os << (*(*S.get_total_component_concentration())(nc))[i] << endl;
      
      fb.close();
    }
  
  // get the auxillary data
  if (auxdata != NULL) 
    {
      int naux = auxdata->NumVectors();
      for (int n=0; n<naux; n++) {
	std::stringstream fname;
	fname << auxnames[n];
	
	fname << "_" << std::setfill('0') << std::setw(5) << cycle << ".dat";
	
	std::filebuf fb;
	fb.open (fname.str().c_str(), std::ios::out);
	ostream os(&fb);      
	
	// now dump the Epetra Vector
	for (int i=0; i< auxdata->MyLength(); i++) 
	  os << (*(*auxdata)(n))[i] << endl;
	fb.close();
      }
    }
}

