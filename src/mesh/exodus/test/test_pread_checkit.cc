// -------------------------------------------------------------
/**
 * @file   test_pread_checkit.cc
 * @author William A. Perkins
 * @date Mon May  2 13:06:27 2011
 * 
 * @brief Routines to check ExodusII::Parallel_Exodus_file instances
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 16, 2010 by William A. Perkins
// Last Change: Mon May  2 13:06:27 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <iostream>
#include <UnitTest++.h>

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>

#include "../Parallel_Exodus_file.hh"

// -------------------------------------------------------------
// checkit
// -------------------------------------------------------------
void 
checkit(Amanzi::Exodus::Parallel_Exodus_file & thefile)
{
  Teuchos::RCP<Amanzi::AmanziMesh::Data::Data> themesh(thefile.read_mesh());
  
  int lcell(themesh->parameters().num_elements_), gcell;
  thefile.comm().SumAll(&lcell, &gcell, 1);

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



