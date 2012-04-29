#include "Vis.hpp"
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
  filebasename = plist.get<string>("File Name Base","amanzi_vis");
  
  if (plist.isParameter("Visualization Times")) {
    Teuchos::Array<double>vis_times = plist.get<Teuchos::Array<double> >("Visualization Times"); 
    // To improve lookup speed while running, we will put these in a stack 
    for (Teuchos::Array<double>::reverse_iterator rit=vis_times.rbegin(); rit<vis_times.rend(); ++rit)
    {
      visualization_times_.push(*rit);
    }
  }
  
  // Grab the cycle parameter list that we wrote in InputParserIS.cc
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
  } else {
    // Errors::Message m("Amanzi::Vis::read_parameters... Cycle Data sublist does not exist in the Visualization Data list");
    //  Exceptions::amanzi_throw(m);
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


void Amanzi::Vis::write_vector(const Epetra_MultiVector& vec, const std::vector<std::string>& names ) const
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


void Amanzi::Vis::write_vector(const Epetra_Vector& vec, const std::string name ) const
{
  viz_output->writeCellDataReal( vec ,name );  
}


void Amanzi::Vis::create_timestep(const double& time, const int& cycle)
{
  viz_output->createTimestep(time,cycle);
}


void Amanzi::Vis::finalize_timestep() const
{
  viz_output->endTimestep();
}

bool Amanzi::Vis::is_disabled() const
{
  return disabled;
}

/**
 * \fn         dump_requested
 * \brief      Visualization dumps can be specified at time steps or at certain times.
 *             This function checks the current time step/time vs a list specified in 
 *             the input.
 *
 *             The input time defaults to -DBL_MAX, as there were many places in 
 *             the code where the step was available, but not the time (easily).
 * \param[in]  cycle - the time step (e.g. n=16)
 * \param[in]  time - optional (e.g. t=34.65s)
 * \returns    bool - true if the current time/step is a visualization point
 */
bool Amanzi::Vis::dump_requested(const int cycle, const double time)
{

  if (!is_disabled())
  {
    // Test time (e.g. t=34.65s)
    if (!visualization_times_.empty() && time!=-std::numeric_limits<double>::max())
    {
      // If the current timestep is equal to the one on the stack, dump
      if (visualization_times_.top()==time)
      {
        visualization_times_.pop();
        return true;
      }
    }
    // Test time step (e.g. n=16)
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
  }
  // if none of the conditions apply we do not write a visualization dump
  return false;

}

