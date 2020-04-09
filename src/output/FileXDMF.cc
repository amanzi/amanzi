/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Reads and writes XDMF files.

/*
  Writes XDMF metadata for file i/o.
*/

#include <fstream>
#include "Teuchos_XMLObject.hpp"

#include "FileXDMF.hh"

namespace Amanzi {

// NOTE: this is in the original code, but then not implemented correctly so
// was not written...
// const std::string FileXDMF::xdmfHeader_ = "<?xml version=\"1.0\"
// ?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
const std::string FileXDMF::xdmfHeader_ = "";

const std::string FileXDMF::DatatypeMap<int>::type = "Int";
const std::string FileXDMF::DatatypeMap<float>::type = "Float";
const std::string FileXDMF::DatatypeMap<double>::type = "Float";

FileXDMF::FileXDMF(const std::string& prefix, int n_nodes, int n_elems,
                   int n_conns)
  : filename_prefix_(prefix),
    xmf_data_visit_(GetDataVisit_()),
    xmf_mesh_visit_(GetMeshVisit_()),
    n_nodes_(n_nodes),
    n_elems_(n_elems),
    n_conns_(n_conns),
    cycle_(0),
    last_mesh_cycle_(-1)
{
  WriteFile_(filename_prefix_ + "_mesh.VisIt.xmf", xmf_mesh_visit_);
  WriteFile_(filename_prefix_ + "_data.VisIt.xmf", xmf_data_visit_);

  // strip the prefix for a basename
  // -- for linux/unix/mac directory names
  std::stringstream ss(prefix);
  char delim('/');
  while (std::getline(ss, basename_, delim)) {}
}


void
FileXDMF::CreateTimestep(const double time, const int cycle,
                         const bool write_mesh)
{
  // create the single timestep data xmf
  xmf_data_timestep_ = GetDataTimestep_(time, cycle, write_mesh);
  cycle_ = cycle;
}

void
FileXDMF::CloseTimestep(const double time, const int cycle,
                        const bool write_mesh)
{
  // get the timestep filename
  std::stringstream data_timestep_stream;
  data_timestep_stream << basename_ << "_data.h5." << cycle << ".xmf";

  // add the timestep data filename to the Visit data
  auto node = FindGridNode_(xmf_data_visit_);
  node.addChild(GetGrid_(data_timestep_stream.str()));

  // write the data files
  WriteFile_(filename_prefix_ + "_data.VisIt.xmf", xmf_data_visit_);
  std::stringstream data_full_filename_stream;
  data_full_filename_stream << filename_prefix_ << "_data.h5." << cycle
                            << ".xmf";
  WriteFile_(data_full_filename_stream.str(), xmf_data_timestep_);

  if (write_mesh) {
    // get the mesh timestep
    auto xmf_mesh_timestep = GetMeshTimestep_(time, cycle);
    std::stringstream mesh_timestep_stream;
    mesh_timestep_stream << basename_ << "_mesh.h5." << cycle << ".xmf";

    // add the mesh timestep file to the Visit mesh
    auto node = FindGridNode_(xmf_mesh_visit_);
    node.addChild(GetGrid_(mesh_timestep_stream.str()));

    // write the mesh files
    WriteFile_(filename_prefix_ + "_mesh.VisIt.xmf", xmf_mesh_visit_);
    std::stringstream mesh_full_filename_stream;
    mesh_full_filename_stream << filename_prefix_ << "_mesh.h5." << cycle
                              << ".xmf";
    WriteFile_(mesh_full_filename_stream.str(), xmf_mesh_timestep);
  }

  cycle_ = -1;
}


void
FileXDMF::WriteFile_(const std::string& filename, const Teuchos::XMLObject& xml)
{
  // write xmf
  std::ofstream of(filename.c_str());
  of << FileXDMF::xdmfHeader_ << xml;
  of.close();
}

Teuchos::XMLObject
FileXDMF::GetMeshTimestep_(const double time, const int cycle)
{
  last_mesh_cycle_ = cycle;

  // TODO(barker): add error handling: can't open/write
  Teuchos::XMLObject mesh("Xdmf");

  std::stringstream mesh_name;
  mesh_name << "Mesh " << cycle;

  // build xml object
  mesh.addChild(GetHeaderLocal_(mesh_name.str(), time, cycle, true));
  return mesh;
}


Teuchos::XMLObject
FileXDMF::GetMeshVisit_()
{
  Teuchos::XMLObject xmf("Xdmf");
  xmf.addAttribute("xmlns:xi", "http://www.w3.org/2001/XInclude");
  xmf.addAttribute("Version", "2.0");

  // build xml object
  xmf.addChild(GetHeaderGlobal_());
  return xmf;
}

Teuchos::XMLObject
FileXDMF::GetDataTimestep_(const double time, const int cycle,
                           const bool write_mesh)
{
  Teuchos::XMLObject xmf("Xdmf");
  xmf.addChild(GetHeaderLocal_("Mesh", time, cycle, write_mesh));
  return xmf;
}


Teuchos::XMLObject
FileXDMF::GetDataVisit_()
{
  // TODO(barker): add error handling: can't open/write
  Teuchos::XMLObject xmf("Xdmf");
  xmf.addAttribute("xmlns:xi", "http://www.w3.org/2001/XInclude");
  xmf.addAttribute("Version", "2.0");

  // build xml object
  xmf.addChild(GetHeaderGlobal_());
  return xmf;
}


Teuchos::XMLObject
FileXDMF::GetHeaderGlobal_()
{
  Teuchos::XMLObject domain("Domain");
  Teuchos::XMLObject grid("Grid");
  grid.addAttribute("GridType", "Collection");
  grid.addAttribute("CollectionType", "Temporal");
  domain.addChild(grid);
  return domain;
}


Teuchos::XMLObject
FileXDMF::GetHeaderLocal_(const std::string name, const double time_val,
                          const int cycle, bool write_mesh)
{
  Teuchos::XMLObject domain("Domain");

  Teuchos::XMLObject grid("Grid");
  grid.addAttribute("Name", name);
  domain.addChild(grid);
  grid.addChild(GetTopo_(cycle, write_mesh));
  grid.addChild(GetGeo_(cycle, write_mesh));

  Teuchos::XMLObject time("Time");
  time.addDouble("Value", time_val);
  grid.addChild(time);
  return domain;
}


Teuchos::XMLObject
FileXDMF::GetTopo_(const int cycle, const bool write_mesh)
{
  // NEW MIXED MESH
  Teuchos::XMLObject topo("Topology");
  topo.addAttribute("TopologyType", "Mixed");
  topo.addInt("NumberOfElements", n_elems_);
  topo.addAttribute("Name", "mixedtopo");

  Teuchos::XMLObject DataItem("DataItem");
  DataItem.addAttribute("DataType", "Int");
  DataItem.addInt("Dimensions", n_conns_);
  DataItem.addAttribute("Format", "HDF");

  std::stringstream tmp1;
  if (write_mesh) {
    tmp1 << basename_ << "_mesh.h5"
         << ":/" << cycle << "/Mesh/MixedElements";
  } else {
    tmp1 << basename_ << "_mesh.h5"
         << ":/" << last_mesh_cycle_ << "/Mesh/MixedElements";
  }

  DataItem.addContent(tmp1.str());
  topo.addChild(DataItem);
  return topo;
}


Teuchos::XMLObject
FileXDMF::GetGeo_(const int cycle, const bool write_mesh)
{
  Teuchos::XMLObject geo("Geometry");
  geo.addAttribute("Name", "geo");
  geo.addAttribute("Type", "XYZ");

  Teuchos::XMLObject DataItem("DataItem");
  DataItem.addAttribute("DataType", "Float");

  std::stringstream tmp1;
  tmp1 << n_nodes_ << " "
       << " 3";
  DataItem.addAttribute("Dimensions", tmp1.str());
  DataItem.addAttribute("Format", "HDF");

  std::stringstream tmp;
  if (write_mesh) {
    tmp << basename_ << "_mesh.h5"
        << ":/" << cycle << "/Mesh/Nodes";
  } else {
    tmp << basename_ << "_mesh.h5"
        << ":/" << last_mesh_cycle_ << "/Mesh/Nodes";
  }
  DataItem.addContent(tmp.str());
  geo.addChild(DataItem);
  return geo;
}

Teuchos::XMLObject
FileXDMF::GetGrid_(std::string grid_filename)
{
  // Create xmlObject grid
  Teuchos::XMLObject xi_include("xi:include");
  xi_include.addAttribute("href", grid_filename);
  xi_include.addAttribute("xpointer", "xpointer(//Xdmf/Domain/Grid)");
  return xi_include;
}


Teuchos::XMLObject
FileXDMF::FindGridNode_(Teuchos::XMLObject xmlobject)
{
  Teuchos::XMLObject node, tmp;

  // Step down to child tag==Domain
  for (int i = 0; i < xmlobject.numChildren(); i++) {
    if (xmlobject.getChild(i).getTag() == "Domain") {
      node = xmlobject.getChild(i);
    }
  }

  // Step down to child tag==Grid and Attribute(GridType==Collection)
  for (int i = 0; i < node.numChildren(); i++) {
    tmp = node.getChild(i);
    if (tmp.getTag() == "Grid" && tmp.hasAttribute("GridType")) {
      if (tmp.getAttribute("GridType") == "Collection") { return tmp; }
    }
  }

  // TODO(barker): return some error indicator
  return node;
}


Teuchos::XMLObject
FileXDMF::FindMeshNode_(Teuchos::XMLObject xmlobject)
{
  Teuchos::XMLObject node, tmp;

  // Step down to child tag==Domain
  for (int i = 0; i < xmlobject.numChildren(); i++) {
    if (xmlobject.getChild(i).getTag() == "Domain") {
      node = xmlobject.getChild(i);
    }
  }

  // Step down to child tag==Grid and Attribute(Name==Mesh)
  for (int i = 0; i < node.numChildren(); i++) {
    tmp = node.getChild(i);
    if (tmp.getTag() == "Grid" && tmp.hasAttribute("Name")) {
      if (tmp.getAttribute("Name") == "Mesh") { return tmp; }
    }
  }

  // TODO(barker): return some error indicator
  return node;
}


Teuchos::XMLObject
FileXDMF::GetField_(const std::string& varname, const std::string& location,
                    const std::string& data_type, int length, int cycle)
{
  Teuchos::XMLObject attribute("Attribute");
  attribute.addAttribute("Name", varname);
  attribute.addAttribute("Type", "Scalar");
  attribute.addAttribute("Center", location);

  Teuchos::XMLObject DataItem("DataItem");
  DataItem.addAttribute("Format", "HDF");
  DataItem.addInt("Dimensions", length);

  DataItem.addAttribute("DataType", data_type);
  std::stringstream tmp;
  tmp << basename_ << "_data.h5"
      << ":" << varname << "/" << cycle;
  DataItem.addContent(tmp.str());
  attribute.addChild(DataItem);
  return attribute;
}


} // namespace Amanzi
