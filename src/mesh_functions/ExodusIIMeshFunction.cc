/*
  Reader for loading data from ExodusII into a vector.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
  
*/
#include "boost/format.hpp"
#include "exodusII.h"

#include "CompositeVector.hh"

namespace Amanzi {
namespace Functions {


void ReadExodusIIMeshFunction(Teuchos::ParameterList& file_list,
                                                     CompositeVector& v)
{
  Epetra_MultiVector& dat_f = *v.ViewComponent("cell", false); 
  int nvectors = dat_f.NumVectors(); 
 
  std::string file_name = file_list.get<std::string>("file"); 
  std::vector<std::string> attributes = 
      file_list.get<Teuchos::Array<std::string> >("attributes").toVector(); 
 
  // open ExodusII file 
  const Epetra_Comm& comm = v.Comm(); 
 
  if (comm.NumProc() > 1) { 
    int ndigits = (int)floor(log10(comm.NumProc())) + 1;
    std::string fmt = boost::str(boost::format("%%s.%%d.%%0%dd") % ndigits);
    file_name = boost::str(boost::format(fmt) % file_name % comm.NumProc() % comm.MyPID());
  } 
 
  int CPU_word_size(8), IO_word_size(0), ierr; 
  float version; 
  int exoid = ex_open(file_name.c_str(), EX_READ, &CPU_word_size, &IO_word_size, &version); 
  if (comm.MyPID() == 0) {
    printf("Trying file: %s ws=%d %d  id=%d\n", file_name.c_str(), CPU_word_size, IO_word_size, exoid); 
  }

  // check if we have to use serial file
  int fail = (exoid < 0) ? 1 : 0;
  int fail_tmp(fail);
  bool distributed_data(true);

  comm.SumAll(&fail_tmp, &fail, 1);
  if (fail == comm.NumProc()) {
    Errors::Message msg("Rao is working on new data layout which we need to proceed.");
    Exceptions::amanzi_throw(msg);

    file_name = file_list.get<std::string>("file"); 
    distributed_data = false;
    if (comm.MyPID() == 0) {
      exoid = ex_open(file_name.c_str(), EX_READ, &CPU_word_size, &IO_word_size, &version); 
      printf("Opening file: %s ws=%d %d  id=%d\n", file_name.c_str(), CPU_word_size, IO_word_size, exoid); 
    }
  } else if (fail > 0) {
    Errors::Message msg("A few parallel Exodus files are missing, but not all.");
    Exceptions::amanzi_throw(msg);
  }
 
  // read database parameters 
  if (comm.MyPID() == 0 || distributed_data) {  
    int dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets; 
    char title[MAX_LINE_LENGTH + 1]; 
    ierr = ex_get_init(exoid, title, &dim, &num_nodes, &num_elem, 
                       &num_elem_blk, &num_node_sets, &num_side_sets); 
 
    int* ids = (int*) calloc(num_elem_blk, sizeof(int)); 
    ierr = ex_get_elem_blk_ids(exoid, ids); 
 
    // read number of variables 
    int num_vars;
    ierr = ex_get_var_param(exoid, "e", &num_vars);
    if (ierr < 0) printf("Exodus file has no variables.\n");

    char* var_names[num_vars];
    for (int i = 0; i < num_vars; i++) {
      var_names[i] = (char*) calloc ((MAX_STR_LENGTH+1), sizeof(char));
    }

    ierr = ex_get_var_names(exoid, "e", num_vars, var_names);
    if (ierr < 0) printf("Exodus file cannot read variable names.\n");

    int var_index(-1), ncells;
    for (int k = 0; k < nvectors; ++k) {
      for (int i = 0; i < num_vars; i++) {
        std::string tmp(var_names[i]);
        if (tmp == attributes[k]) var_index = i + 1;
      }
      if (var_index < 0) printf("Exodus file has no variable \"%s\".\n", attributes[k].c_str());
      printf("Variable \"%s\" has index %d.\n", attributes[k].c_str(), var_index);

      // read variable with the k-th attribute
      int offset = 0; 
      char elem_type[MAX_LINE_LENGTH + 1]; 
      for (int i = 0; i < num_elem_blk; i++) { 
        int num_elem_this_blk, num_attr, num_nodes_elem; 
        ierr = ex_get_elem_block(exoid, ids[i], elem_type, &num_elem_this_blk, 
                                 &num_nodes_elem, &num_attr); 
 
        double* var_values = (double*) calloc(num_elem_this_blk, sizeof(double)); 
        ierr = ex_get_elem_var(exoid, 1, var_index, ids[i], num_elem_this_blk, var_values); 
 
        for (int n = 0; n < num_elem_this_blk; n++) { 
          int c = n + offset; 
          dat_f[k][c] = var_values[n]; 
        } 
        free(var_values); 
        printf("MyPID=%d  ierr=%d  id=%d  ncells=%d\n", comm.MyPID(), ierr, ids[i], num_elem_this_blk); 
 
        offset += num_elem_this_blk; 
      } 
      ncells = offset; 
    }

    for (int i = 0; i < num_vars; i++) {
      free(var_names[i]);
    }
 
    ierr = ex_close(exoid); 
    printf("Closing file: %s ncells=%d error=%d\n", file_name.c_str(), ncells, ierr); 
  }

}



} // namespace Functions
} // namespace Amanzi
                        
