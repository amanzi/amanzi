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
  Teuchos::reduceAll(*thefile.Comm(), Teuchos::REDUCE_SUM, 1, &lcell, &gcell);

  auto cmap = thefile.cellmap();
  CHECK_EQUAL(cmap->getGlobalNumElements(), gcell);
  CHECK_EQUAL(cmap->getNodeNumElements(), lcell);
  CHECK_EQUAL(cmap->getMinGlobalIndex(), 1);
  CHECK_EQUAL(cmap->getMaxGlobalIndex(), gcell);
  CHECK(cmap->isOneToOne());
  

  int lvert(themesh->parameters().num_nodes_);
  auto vmap = thefile.vertexmap();
  CHECK_EQUAL(vmap->getNodeNumElements(), lvert);
  CHECK(cmap->isOneToOne()); 

  // FIXME: do some checking on global indexes, whatever that might me
}



