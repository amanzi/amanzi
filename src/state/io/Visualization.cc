/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
      Markus Berndt
*/

//! Visualization writes data and meshes to files.

/*
  Writes to file for visualization using generic Output object.
*/

#include "OutputFactory.hh"
#include "Mesh.hh"
#include "UniqueHelpers.hh"

#include "Visualization.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
Visualization::Visualization(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                             const Teuchos::RCP<AmanziMesh::Mesh>& mesh)
  : IOEvent(plist), output_(std::move(OutputFactory::CreateForVis(*plist, mesh)))
{}

void
Visualization::CreateFile(const double& time, const int& cycle)
{
  output_->CreateFile(time, cycle);
}

void
Visualization::FinalizeFile()
{
  output_->FinalizeFile();
}


// -----------------------------------------------------------------------------
// Write a field with region information
// -----------------------------------------------------------------------------
// void Visualization::WriteRegions() {
// if (regions_.size() > 0) {
//   for (std::map<std::string, Teuchos::Array<std::string>>::const_iterator it
//   =
//            regions_.begin();
//        it != regions_.end(); ++it) {
//     // first make an Epetra_Vector to hold the region information
//     Epetra_Vector reg(mesh_->cell_map(false), true);

//     // loop over the regions and initialize the reg array
//     double reg_index = 1.0;
//     for (Teuchos::Array<std::string>::const_iterator
//              reg_it = (it->second).begin();
//          reg_it != (it->second).end(); ++reg_it, reg_index += 1.0) {
//       // only do something if the user provided a valid region name
//       // for a region that consists of cells
//       if (mesh_->valid_set_name(*reg_it, AmanziMesh::CELL)) {
//         AmanziMesh::Entity_ID_List ids;
//         mesh_->get_set_entities(*reg_it, AmanziMesh::CELL,
//         AmanziMesh::Parallel_type::OWNED,
//                                 &ids);

//         for (AmanziMesh::Entity_ID_List::const_iterator rit = ids.begin();
//              rit != ids.end(); ++rit) {
//           reg[*rit] = reg_index;
//         }
//       }
//     }

//     Write(it->first, reg);
//   }
// }
//}

// -----------------------------------------------------------------------------
// Write a field with region information
// -----------------------------------------------------------------------------
// void Visualization::WritePartition() {
// if (write_partition_) {
//   // first make an Epetra_Vector to hold the partitions information
//   Epetra_Vector reg(mesh_->cell_map(false), false);
//   // loop over the regions and initialize the reg array
//   double part_index = static_cast<double>(mesh_->get_comm()->MyPID());
//   reg.putScalar(part_index);

//   Write("partition", reg);
// }
//}


} // namespace Amanzi
