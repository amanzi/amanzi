/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_GMV_MESH_HH_
#define AMANZI_GMV_MESH_HH_

#include <string>

extern "C"
{
#include "gmvwrite.h"
}

#include "Mesh.hh"


namespace Amanzi {
namespace GMV {

// Write a GMV file containing only mesh data to be used as a "fromfile".
void
create_mesh_file(const AmanziMesh::Mesh& meshs, std::string filename);

// Opens and initializes a GMV file for writing which references a "fromfile"
// for mesh definition.
void
open_data_file(std::string mesh_fromfile, std::string filename_path,
               unsigned int num_nodes, unsigned int num_cells);

// Opens and initializes a GMV file for writing which references a "fromfile"
// for mesh definition. adds a suffix of the type .000302 with, in this case,
// cycleno=302 and digits= 6
void
open_data_file(std::string mesh_fromfile, std::string filename,
               unsigned int num_nodes, unsigned int num_cells,
               unsigned int cycleno, unsigned int digits);

// Opens and initializes a GMV file which contains mesh data, i.e. doesn't use a
// "fromfile".
void
open_data_file(const AmanziMesh::Mesh& meshs, std::string filename);

// Opens and initializes a GMV file which contains mesh data, i.e. doesn't use a
// "fromfile" adds a suffix of the type .000302 with, in this case, cycleno=302
// and digits= 6
void
open_data_file(const AmanziMesh::Mesh& mesh, std::string filename,
               unsigned int cycleno, unsigned int digits);

// start the variables section (call this after write_cycle or write_time)
void
start_data();

// Writes node data to files which has previously been opened with
// open_data_file.
void
write_node_data(const Epetra_Vector& x, std::string varname);
void
write_node_data(const Epetra_MultiVector& x, const unsigned int component,
                std::string varname);

// Writes cell data to files which has previously been opened with
// open_data_file.
void
write_cell_data(const Epetra_Vector& x, std::string varname);
void
write_cell_data(const Epetra_MultiVector& x, const unsigned int component,
                std::string varname);

// Writes cell data to files which has previously been opened with
// open_data_file.
void
write_face_data(const Epetra_Vector& x, std::string varname);

// Writes the cycle number and time
void
write_cycle(const int cycle);
void
write_time(const double time);

// Finalizes a GMV file which has previously been opened with open_data_file.
void
close_data_file();

// modify suffix string such that it is of form ".0203", where, in this example,
// cycleno=203
void
suffix_no(std::string& suffix, unsigned int cycleno);

} // namespace GMV
} // namespace Amanzi

#endif
