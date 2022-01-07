/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

#include <iostream>
#include <ostream>

#include "boost/filesystem.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Epetra_Vector.h"

// Amanzi
#include "CompositeVector.hh"
#include "DomainSet.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "StringExt.hh"

// Amanzi::State
#include "State.hh"
#include "IO.hh"

namespace Amanzi {

// Non-member function for vis.
void WriteVis(Visualization& vis, const State& S)
{
  if (!vis.is_disabled()) {
    // Create the new time step
    vis.CreateTimestep(S.time(), S.cycle(), vis.get_tag());

    // Write all fields to the visualization file, the fields know if they
    // need to be written.
    for (auto r = S.data_begin(); r != S.data_end(); ++r) {
      if (S.GetRecord(r->first).ValidType<CompositeVector>())
        r->second->WriteVis(vis);
    }
    vis.WriteRegions();
    vis.WritePartition();
    vis.FinalizeTimestep();
  }
}


// Non-member function for checkpointing.
void WriteCheckpoint(Checkpoint& chkp, const Comm_ptr_type& comm,
                     const State& S, bool final)
{
  if (!chkp.is_disabled()) {
    // chkp.SetFinal(final);
    chkp.CreateFile(S.cycle());
    for (auto r = S.data_begin(); r != S.data_end(); ++r) {
      r->second->WriteCheckpoint(chkp);
    }
    chkp.Write("mpi_num_procs", comm->NumProc());
    chkp.Finalize();
  }
}


// Non-member function for checkpointing.
void ReadCheckpoint(const Comm_ptr_type& comm, State& S, 
                    const std::string& filename)
{
  Checkpoint chkp(filename, comm);

  // Load the number of processes and ensure they are the same.
  int num_procs(-1);
  chkp.Read("mpi_comm_world_rank", num_procs);
  if (comm->NumProc() != num_procs) {
    std::stringstream ss;
    ss << "Requested checkpoint file " << filename << " was created on " 
       << num_procs << " processes, making it incompatible with this run on "
       << comm->NumProc() << " processes.";
    Errors::Message message(ss.str());
    throw(message);
  }

  // load metadata
  chkp.ReadAttributes(S);

  // load the data
  for (auto data = S.data_begin(); data != S.data_end(); ++data) {
    data->second->ReadCheckpoint(chkp);
  }

  chkp.Finalize();
}


// Non-member function for checkpointing.
double ReadCheckpointInitialTime(const Comm_ptr_type& comm,
                                 std::string filename)
{
  if (!Keys::ends_with(filename, ".h5")) {
    // new style checkpoint
    boost::filesystem::path filepath = boost::filesystem::path(filename) / "domain.h5";
    filename = filepath.string();
  }
  HDF5_MPI checkpoint(comm, filename);

  // load the attributes
  double time(0.);
  checkpoint.open_h5file();
  checkpoint.readAttrReal(time, "time");
  checkpoint.close_h5file();
  return time;
}


// Non-member function for checkpointing position.
int ReadCheckpointPosition(const Comm_ptr_type& comm, std::string filename)
{
  if (!Keys::ends_with(filename, ".h5")) {
    // new style checkpoint
    boost::filesystem::path filepath = boost::filesystem::path(filename) / "domain.h5";
    filename = filepath.string();
  }
  HDF5_MPI checkpoint(comm, filename);

  // load the attributes
  int pos = 0;
  checkpoint.open_h5file();
  checkpoint.readAttrInt(pos, "position");
  checkpoint.close_h5file();
  return pos;
}


// Non-member function for checkpointing observations.
void ReadCheckpointObservations(const Comm_ptr_type& comm,
                                std::string filename,
                                Amanzi::ObservationData& obs_data)
{
  if (!Keys::ends_with(filename, ".h5")) {
    // new style checkpoint
    boost::filesystem::path filepath = boost::filesystem::path(filename) / "domain.h5";
    filename = filepath.string();
  }

  HDF5_MPI checkpoint(comm, filename);
  checkpoint.open_h5file();

  // read observations
  int nlabels, ndata(0), ndata_glb(0);
  int* nobs;
  char** tmp_labels;
  double* tmp_data(NULL);

  checkpoint.readDataString(&tmp_labels, &nlabels, "obs_names");
  if (nlabels > 0) { 
    checkpoint.readAttrInt(&nobs, &nlabels, "obs_numbers");
  }
  for (int i = 0; i < nlabels; ++i) ndata_glb += 2 * nobs[i];
  ndata = (comm->MyPID() == 0) ? ndata_glb : 0;
  checkpoint.readDatasetReal(&tmp_data, ndata, "obs_values");

  checkpoint.close_h5file();

  // populated observations on root
  if (comm->MyPID() == 0) {
    int m(0);
    Amanzi::ObservationData::DataQuadruple data_quad;

    for (int i = 0; i < nlabels; ++i) {
      std::vector<ObservationData::DataQuadruple>& od = obs_data[tmp_labels[i]];
      for (int k = 0; k < nobs[i]; ++k) {
        data_quad.time = tmp_data[m++];
        data_quad.value = tmp_data[m++];
        data_quad.is_valid = true;
        od.push_back(data_quad);
      }
    }
  }

  // clean memory
  for (int i = 0; i < nlabels; i++) free(tmp_labels[i]);
  if (nlabels > 0) {
    free(tmp_labels);
    free(nobs);
    if (tmp_data != NULL) free(tmp_data); 
  }
}


// Non-member function for deforming the mesh after reading a checkpoint file
// that contains the vertex coordinate field (this is written by deformation pks)
//
// FIX ME: Refactor this to make the name more general.  Should align with a
// mesh name prefix or something, and the coordinates should be written by
// state in WriteCheckpoint if mesh IsDeformableMesh() --ETC
void DeformCheckpointMesh(State& S, Key domain)
{
  if (S.HasData("vertex coordinate")) { // only deform mesh if vertex coordinate
                                        // field exists
    AmanziMesh::Mesh* write_access_mesh = const_cast<AmanziMesh::Mesh*>(&*S.GetMesh());

    // get vertex coordinates state
    const CompositeVector& vc = S.Get<CompositeVector>("vertex coordinate");
    vc.ScatterMasterToGhosted("node");
    const Epetra_MultiVector& vc_n = *vc.ViewComponent("node", true);

    int dim = write_access_mesh->space_dimension();
    Amanzi::AmanziMesh::Entity_ID_List nodeids;
    Amanzi::AmanziGeometry::Point new_coords(dim);
    AmanziGeometry::Point_List new_pos, final_pos;

    int nV = vc_n.MyLength();
    for (int n = 0; n != nV; ++n) {
      for (int k = 0; k != dim; ++k)
        new_coords[k] = vc_n[k][n];

      // push back for deform method
      nodeids.push_back(n);
      new_pos.push_back(new_coords);
    }

    // deform the mesh
    if (Keys::starts_with(domain, "column"))
      write_access_mesh->deform(nodeids, new_pos, false, &final_pos);
    else
      write_access_mesh->deform(nodeids, new_pos, true, &final_pos);
  }
}


// Non-member function for statistics
void WriteStateStatistics(const State& S, const VerboseObject& vo)
{
  // sort data in alphabetic order
  std::set<std::string> sorted;
  for (auto it = S.data_begin(); it != S.data_end(); ++it) {
    sorted.insert(it->first);
  }

  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    Teuchos::OSTab tab = vo.getOSTab();
    *vo.os() << "\nField                                    Min/Max/Avg" << std::endl;

    // for (auto it = S.data_begin(); it != S.data_end(); ++it) {
    for (auto name : sorted) {
      if (name.size() > 33) replace_all(name, "temperature", "temp");
      if (name.size() > 33) replace_all(name, "internal_energy", "ie");
      if (name.size() > 33) replace_all(name, "molar", "mol");

      if (S.GetRecord(name).ValidType<CompositeVector>()) {
        std::map<std::string, double> vmin, vmax, vavg;
        S.Get<CompositeVector>(name).MinValue(vmin);
        S.Get<CompositeVector>(name).MaxValue(vmax);
        S.Get<CompositeVector>(name).MeanValue(vavg);

        for (auto c_it = vmin.begin(); c_it != vmin.end(); ++c_it) {
          std::string namedot(name), name_comp(c_it->first);
          if (vmin.size() != 1) namedot.append("." + name_comp);
          namedot.resize(40, '.');
          *vo.os() << namedot << " " << c_it->second << " / " 
                    << vmax[name_comp] << " / " << vavg[name_comp] << std::endl;
        }
      }
      else if (S.GetRecord(name).ValidType<double>()) {
        double vmin = S.Get<double>(name);
        name.resize(40, '.');
        *vo.os() << name << " " << vmin << std::endl;
      }
      else if (S.GetRecord(name).ValidType<AmanziGeometry::Point>()) {
        const auto& p = S.Get<AmanziGeometry::Point>(name);
        name.resize(40, '.');
        *vo.os() << name;
        for (int i = 0; i < p.dim(); ++i) *vo.os() << " " << p[i];
        *vo.os() << std::endl;
      }
    }
  }
}

}  // namespace Amanzi
