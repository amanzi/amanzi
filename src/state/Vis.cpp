#include "Vis.hpp"
#include "Epetra_MpiComm.h"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

#include "boost/filesystem.hpp"

Amanzi::Vis::Vis (Teuchos::ParameterList& plist_, Epetra_MpiComm* comm_):
    plist(plist_), disabled(false), comm(comm_), hasCycleData_(false), viz_output(NULL)
{
  read_parameters(plist);

  viz_output = new Amanzi::HDF5_MPI(*comm);
  viz_output->setTrackXdmf(true);
}

// this constructor makes an object that will not create any output
Amanzi::Vis::Vis (): disabled(true), hasCycleData_(false), viz_output(NULL)
{
}


void Amanzi::Vis::read_parameters(Teuchos::ParameterList& plist)
{
  filebasename = plist.get<string>("File Name Base","amanzi_vis");
  
  boost::filesystem::path fbn(filebasename);
  boost::filesystem::path parent = fbn.parent_path();
  if (!parent.empty()) {
    // we need to check whether the parent path exists
    if (!boost::filesystem::exists(parent)) {
      Errors::Message m("The file path '"+parent.string()+"' used for visualization does not exist.");
      Exceptions::amanzi_throw(m);    
    }
    if (!boost::filesystem::is_directory(parent)) {
      Errors::Message m("The file path '"+parent.string()+"' used for visualization is not a directory.");
      Exceptions::amanzi_throw(m);    
    }
    if ( (boost::filesystem::status(parent).permissions() & boost::filesystem::owner_write) == 0) {
      Errors::Message m("The directory '"+parent.string()+"' used for visualization is not writeable.");
      Exceptions::amanzi_throw(m);       
    }
  }
  
  // get the time start period stop data and discrete times
  if (plist.isSublist("time start period stop")) {
    Teuchos::ParameterList& tsps_list = plist.sublist("time start period stop");
    for (Teuchos::ParameterList::ConstIterator it = tsps_list.begin(); it != tsps_list.end(); ++it) {
      std::string name = it->first;
      if (tsps_list.isSublist(name)) {
        Teuchos::ParameterList& itlist = tsps_list.sublist(name);
        if (itlist.isParameter("start period stop")) {
          Teuchos::Array<double> sps = itlist.get<Teuchos::Array<double> >("start period stop");
          time_sps.push_back(sps.toVector());
        }
      }
    }
  }
  if (plist.isParameter("times")) {
    Teuchos::Array<double> vtimes = plist.get<Teuchos::Array<double> >("times");
    times.push_back(vtimes.toVector());
  }
  // Grab the cycle start period stop data
  if ( plist.isSublist("cycle start period stop") ) {
    Teuchos::ParameterList &sps_list = plist.sublist("cycle start period stop");
    // for now we only allow one cycle start period stop entry
    if ( ++ sps_list.begin() == sps_list.end() ) {
      std::string name = sps_list.begin()->first;
      if ( sps_list.isSublist(name)) {
        if (sps_list.sublist(name).isParameter("start period stop")) {
          Teuchos::Array<int> sps = sps_list.sublist(name).get<Teuchos::Array<int> >("start period stop");
          interval = sps[1];
          start = sps[0];
          end = sps[2];
          hasCycleData_ = true;
        }
        if (sps_list.sublist(name).isParameter("cycles")) {
          steps = sps_list.sublist(name).get<Teuchos::Array<int> >("cycles");
          hasCycleData_ = true;
        }
      }
    }
  }
}


Amanzi::Vis::~Vis ()
{
  if (!viz_output) delete viz_output;
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

    Teuchos::Array<std::string> regions;
    if (plist.isParameter("Regions")) {
      regions = plist.get<Teuchos::Array<std::string> >("Regions");
    }
    for (Teuchos::Array<std::string>::const_iterator reg=regions.begin(); reg != regions.end(); ++reg) {
      if (!mesh.valid_set_name(*reg,Amanzi::AmanziMesh::CELL)) {
        Errors::Message m("Amanzi::Vis::create_files... Region \"" + *reg + "\" specified in the Visualization Data list is not a valid region.");
        Exceptions::amanzi_throw(m);
      }
      int local_region_size = mesh.get_set_size(*reg,
                                                Amanzi::AmanziMesh::CELL,
                                                Amanzi::AmanziMesh::OWNED);
      Amanzi::AmanziMesh::Entity_ID_List cell_ids(local_region_size);
      mesh.get_set_entities(*reg, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
                            &cell_ids);

      int global_region_size(0);
      comm->SumAll(&local_region_size,&global_region_size,1);
      Epetra_Map region_map(global_region_size, local_region_size, 0, *comm);
      Teuchos::RCP<Epetra_Vector> mesh_region = Teuchos::rcp(new Epetra_Vector(region_map,false));


      int *indices = new int[local_region_size];
      for (int i=0; i<local_region_size; ++i) indices[i] = i;
      double *values = new double[local_region_size];
      for (int i=0; i<local_region_size; ++i) values[i] = mesh.cell_map(false).GID(cell_ids[i]);

      mesh_region->ReplaceMyValues(local_region_size,values,indices);

      delete [] values;
      delete [] indices;

      viz_output->writeMeshRegion(mesh, *mesh_region, *reg);
    }
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

void Amanzi::Vis::register_with_time_step_manager(TimeStepManager& TSM) {
  if (plist.isParameter("times")) {
    Teuchos::Array<double> times = plist.get<Teuchos::Array<double> >("times");
    TSM.RegisterTimeEvent(times.toVector());
  }
  if (plist.isSublist("time start period stop")) {
    Teuchos::ParameterList &sps_list = plist.sublist("time start period stop");
    for (Teuchos::ParameterList::ConstIterator it = sps_list.begin(); it != sps_list.end(); ++it) {
      std::string name = it->first;
      if ( sps_list.isSublist(name)) {
        Teuchos::ParameterList& itlist = sps_list.sublist(name);
        Teuchos::Array<double> sps = itlist.get<Teuchos::Array<double> >("start period stop");
        TSM.RegisterTimeEvent(sps[0],sps[1],sps[2]);
      }
    }
  }
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
bool Amanzi::Vis::dump_requested(const int cycle, const double time) {
  return dump_requested(cycle) || dump_requested(time);
}


bool Amanzi::Vis::dump_requested(const double time) {
  if (!is_disabled())  {
    // loop over the start period stop triplets
    if (time_sps.size() > 0) {
      for (std::vector<std::vector<double> >::const_iterator it = time_sps.begin(); it != time_sps.end(); ++it) {
        if ( (Amanzi::near_equal(time,(*it)[0]) || time >= (*it)[0]) &&
             (Amanzi::near_equal((*it)[2],-1.0) || time <= (*it)[2] || Amanzi::near_equal(time,(*it)[2]) ) ) {
          if (Amanzi::near_equal(time,(*it)[0])) {
            return true;
          }
          double n_per_tmp = (time - (*it)[0])/(*it)[1];
          double n_periods = floor(n_per_tmp);
          if (Amanzi::near_equal(n_periods+1.0,n_per_tmp)) n_periods += 1.0;
          double tmp = (*it)[0] + n_periods*(*it)[1];
          if (Amanzi::near_equal(time,tmp)) {
            return true;
          }
        }
      }
    }
    // loop over the time arrays
    for (std::vector<std::vector<double> >::const_iterator it = times.begin(); it != times.end(); ++it) {
      for (std::vector<double>::const_iterator j= (*it).begin(); j!=(*it).end(); ++j) {
        if (Amanzi::near_equal(*j,time)) {
          return true;
        }
      }
    }
  }
  return false;
}


bool Amanzi::Vis::dump_requested(const int cycle) {
  if (!is_disabled())  {
    if (hasCycleData_) {
      // Test time step (e.g. n=16)
      if (steps.size() > 0) {
        for (int i=0; i<steps.size(); i++) {
          if (cycle == steps[i]) {
            return true;
          }
        }
      } else if ( (end<0) || (cycle<=end) ) {
        if (start<=cycle) {
          int cycle_loc = cycle - start;

          if (cycle_loc % interval == 0) {
            return true;
          }
        }
      }
    }
  }
  // if none of the conditions apply we do not write a visualization dump
  return false;
}
