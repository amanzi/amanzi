/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins
*/

#include <iostream>
#include "UnitTest++.h"

#include "Epetra_Map.h"
#include "AmanziComm.hh"

#include "../Parallel_Exodus_file.hh"

// -------------------------------------------------------------
// checkit
// -------------------------------------------------------------
void 
checkit(Amanzi::Exodus::Parallel_Exodus_file & thefile)
{
  Teuchos::RCP<Amanzi::AmanziMesh::Data::Data> themesh(thefile.read_mesh());
  
  int lcell(themesh->parameters().num_elements_), gcell;
  thefile.Comm()->SumAll(&lcell, &gcell, 1);

  Teuchos::RCP<Epetra_Map> cmap(thefile.cellmap());
  CHECK_EQUAL(cmap->NumGlobalElements(), gcell);
  CHECK_EQUAL(cmap->NumMyElements(), lcell);
  CHECK_EQUAL(cmap->MinAllGID(), 1);
  CHECK_EQUAL(cmap->MaxAllGID(), gcell);
  CHECK(cmap->IsOneToOne());
  

  int lvert(themesh->parameters().num_nodes_);
  Teuchos::RCP<Epetra_Map> vmap(thefile.vertexmap());
  CHECK_EQUAL(vmap->NumMyElements(), lvert);
  CHECK(cmap->IsOneToOne()); 

  // FIXME: do some checking on global indexes, whatever that might me
}



