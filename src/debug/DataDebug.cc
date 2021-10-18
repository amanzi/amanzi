#include "DataDebug.hh"

namespace Amanzi {

DataDebug::DataDebug(Teuchos::RCP<AmanziMesh::Mesh> mesh) :
    mesh_(mesh)
{
}


void DataDebug::write_region_data(std::string& region_name, 
                                  const Epetra_Vector& data, 
                                  std::string& description) {

  if (!mesh_->isValidSetName(region_name, AmanziMesh::Entity_kind::CELL)) {
    throw std::exception();
  }
  unsigned int mesh_block_size = mesh_->getSetSize(region_name,
                                                     AmanziMesh::Entity_kind::CELL,
                                                     AmanziMesh::Parallel_type::OWNED);
  AmanziMesh::Entity_ID_List cell_ids(mesh_block_size);
  cell_ids = mesh_->getSetEntities(region_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
  
  std::cerr << "Printing " << description << " on region " << region_name << std::endl;
  for(AmanziMesh::Entity_ID_List::iterator c = cell_ids.begin(); c != cell_ids.end(); ++c) {
    std::cerr << std::fixed 
              << description << "(" << data.Map().GID(*c) << ") = " << data[*c] << std::endl;
  }
}


void DataDebug::write_region_statistics(std::string& region_name, 
                                        const Epetra_Vector& data, 
                                        std::string& description) {

  if (!mesh_->isValidSetName(region_name, AmanziMesh::Entity_kind::CELL)) {
    throw std::exception();
  }
  unsigned int mesh_block_size = mesh_->getSetSize(region_name,
                                                     AmanziMesh::Entity_kind::CELL,
                                                     AmanziMesh::Parallel_type::OWNED);
  AmanziMesh::Entity_ID_List cell_ids(mesh_block_size);
  cell_ids = mesh_->getSetEntities(region_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED);
  
  // find min and max and their indices
  int max_index(0), min_index(0);
  double max_value(-1e-99), min_value(1e99);
  for(AmanziMesh::Entity_ID_List::iterator c = cell_ids.begin(); c != cell_ids.end(); ++c) {
    if (data[*c] > max_value) {
      max_value = data[*c];
      max_index = data.Map().GID(*c);
    }
    if (data[*c] < min_value) {
      min_value = data[*c];
      min_index = data.Map().GID(*c);
    }      
  }
  
  int num_procs(mesh_->getComm()->NumProc());
  int my_proc(mesh_->getComm()->MyPID());

  double* all_values = new double[num_procs];
  int* all_indices = new int[num_procs];
  mesh_->getComm()->GatherAll(&max_value,all_values,1);
  mesh_->getComm()->GatherAll(&max_index,all_indices,1);
 

  // find global max value and index
  max_value = all_values[0];
  for (int i=1; i<num_procs; ++i) {
    if (all_values[i]>max_value) {
      max_value = all_values[i];
      max_index = all_indices[i];
    }
  }
  
  mesh_->getComm()->GatherAll(&min_value,all_values,1);
  mesh_->getComm()->GatherAll(&min_index,all_indices,1);
  // find global min value and index
  min_value = all_values[0];
  for (int i=1; i<num_procs; ++i) {
    if (all_values[i]<min_value) {
      min_value = all_values[i];
      min_index = all_indices[i];
    }
  }  
  
  delete [] all_indices;
  delete [] all_values;

  // result is the same on all processors, only print once
  if (my_proc == 0) {
    std::cerr << "Printing min/max of " << description << " on region " << region_name << std::endl;   
    std::cerr << std::fixed 
              << description << " min = " << min_value << " in cell " << min_index << std::endl;
    std::cerr << std::fixed 
              << description << " max = " << max_value << " in cell " << max_index << std::endl;
  }
}

}
