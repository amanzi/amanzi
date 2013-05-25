/*
This is the atate component of the Amanzi code.
 
Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "exodusII.h"

#include "State.hh"

/* ******************************************************************
* We do not verify that the file name matches that in the Mesh
* section of the input Spec.
****************************************************************** */
void State::initialize_from_file_list()
{
  if (! parameter_list.isSublist("File initialization")) return;
  Teuchos::ParameterList& file_list = parameter_list.sublist("File initialization");

  // set permeability, will become separate routines later
  if (! file_list.isSublist("absolute permeability")) return;
  Teuchos::ParameterList& sublist = file_list.sublist("absolute permeability");

  std::string region = sublist.get<std::string>("region");
  std::string file_name = sublist.get<std::string>("file");
  std::string attribute_name = sublist.get<std::string>("attribute");

  // open ExodusII file
  int CPU_word_size(0), IO_word_size(0), ierr; 
  float version;
  int exoid = ex_open(file_name.c_str(), EX_READ, &CPU_word_size, &IO_word_size, &version);

  // read database parameters
  int dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets;
  char title[MAX_LINE_LENGTH + 1];
  ierr = ex_get_init(exoid, title, &dim, &num_nodes, &num_elem, 
                     &num_elem_blk, &num_node_sets, &num_side_sets);

  float attrib[1000];
  int* ids = (int*) calloc(num_elem_blk, sizeof(int));
  ierr = ex_get_elem_blk_ids(exoid, ids);

  for (int i = 0; i < num_elem_blk; i++) {
    ierr = ex_get_elem_attr(exoid, ids[i], attrib);
    cout << ierr << " " << ids[i] << endl;
  }

  ierr = ex_close(exoid); 
  cout << "Result of reading data: " << ierr << endl;
}

