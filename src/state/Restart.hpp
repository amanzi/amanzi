#ifndef _RESTART_HPP_
#define _RESTART_HPP_

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Comm.h"
#include "State.hpp"
#include "hdf5mpi_mesh.hh"

namespace Amanzi {

  class Restart {

  public:
    Restart(Teuchos::ParameterList& plist, Epetra_MpiComm *comm); 
    Restart(); // this object will not create any output 
    ~Restart();
   
    void dump_state (State& S);
    void read_state (State& S);
    void create_files ();
    void read_parameters(Teuchos::ParameterList& plist);
    bool dump_requested(int cycle);

  private:    
    std::string filebasename; 
    Teuchos::ParameterList plist;
    
    int restart_cycle;
    int number_of_cycle_intervals;

    std::vector<int> cycle_freq;

    std::vector<int> cycle_start;
    std::vector<int> cycle_end;

    Amanzi::HDF5_MPI *restart_output; 

    // disable restart dumps alltogether
    bool disabled;

    // the Epetra communicator
    Epetra_MpiComm *comm;
  };

}
#endif
