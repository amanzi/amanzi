#include "Vis.hh"
#include "Epetra_MpiComm.h"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"


Amanzi::Vis::Vis (Teuchos::ParameterList& plist_, Epetra_MpiComm* comm_):
  plist(plist_), disabled(false), comm(comm_)
{
  read_parameters(plist);

  viz_output = new Amanzi::HDF5_MPI(*comm);
  viz_output->setTrackXdmf(true);
}

// this constructor makes an object that will not create any output
Amanzi::Vis::Vis (): disabled(true)
{
}


void Amanzi::Vis::read_parameters(Teuchos::ParameterList& plist)
{
  filebasename = plist.get<string>("file base name","amanzi_vis");
  
  if ( plist.isSublist("Cycle Data") ) 
    {
      Teuchos::ParameterList &ilist = plist.sublist("Cycle Data");
      
      interval = ilist.get<int>("Interval");
      start = ilist.get<int>("Start");
      end = ilist.get<int>("End");
      
      if (ilist.isParameter("Steps"))
	{
	  steps = ilist.get<Teuchos::Array<int> >("Steps");  
	}
    }  
  else
    {
      Errors::Message m("Amanzi::Vis::read_parameters... Cycle Data sublist does not exist on the Visualization Data list");
      Exceptions::amanzi_throw(m);
    }
}


Amanzi::Vis::~Vis () 
{
  delete viz_output;
}


void Amanzi::Vis::create_files(const Amanzi::AmanziMesh::Mesh& mesh)
{
  if (!is_disabled())
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


void Amanzi::Vis::write_vector(const Epetra_MultiVector& vec, const std::vector<std::string>& names )
{
  if (names.size() < vec.NumVectors()) 
    {
      Errors::Message m("Amanzi::Vis::write_vector... not enough names were specified for the the components of the multi vector");
      Exceptions::amanzi_throw(m);
    }
  
  for (int i=0; i< vec.NumVectors(); i++)
    {
      viz_output->writeCellDataReal( *vec(i), names[i] );  
    }
}


void Amanzi::Vis::write_vector(const Epetra_Vector& vec, std::string name )
{
  viz_output->writeCellDataReal( vec ,name );
}

void Amanzi::Vis::create_timestep(const double& time, const int& cycle)
{
  viz_output->createTimestep(time, cycle);
}


void Amanzi::Vis::finalize_timestep()
{
  viz_output->endTimestep();
}

const bool Amanzi::Vis::is_disabled()
{
  return disabled;
}


const bool Amanzi::Vis::dump_requested(int cycle)
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
      
  // if none of the conditions apply we do not write a visualization dump
  return false;

}

