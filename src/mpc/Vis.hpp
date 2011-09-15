#ifndef _VIS_HPP_
#define _VIS_HPP_

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Comm.h"
#include "State.hpp"
#include "Mesh.hh"
#include "hdf5mpi_mesh.hh"

namespace Amanzi {

  class Vis {

  public:
    Vis(Teuchos::ParameterList& plist, const Epetra_MpiComm& comm); 
    Vis(); // this object will not create any output 
    ~Vis();
   
    void dump_state (double time, double prev_time, int cycle, State& S);
    void create_files (Amanzi::AmanziMesh::Mesh& mesh);
    void read_parameters(Teuchos::ParameterList& plist);
    bool dump_requested(double time, double prev_time, int cycle);
    void set_compnames(std::vector<std::string>& compnames_);

  private:    
    std::string filebasename; 
    Teuchos::ParameterList plist;
    
    int vizdump_cycle;
    int vizdump_time_count;
    double vizdump_time;
    
    int number_of_cycle_intervals;
    int number_of_time_intervals;

    std::vector<int> cycle_freq;
    std::vector<double> time_freq;

    std::vector<int> cycle_start;
    std::vector<int> cycle_end;
    std::vector<double> time_start;
    std::vector<double> time_end;

    Amanzi::HDF5_MPI *viz_output; 

    // names for the components
    std::vector<std::string> compnames;

    // disable visualization dumps alltogether
    bool disabled;

  };

}
#endif
