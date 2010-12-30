// -------------------------------------------------------------
/**
 * @file   stk_hex_mesh_test.cc
 * @author William A. Perkins
 * @date Wed Dec 29 14:05:33 2010
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December 13, 2010 by William A. Perkins
// Last Change: Wed Dec 29 14:05:33 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------


static const char* SCCS_ID = "$Id$ Battelle PNL";

#include <iostream>

#include <Teuchos_GlobalMPISession.hpp>
#include <Epetra_Vector.h>
#include <Isorropia_EpetraPartitioner.hpp>

#include <boost/format.hpp>

#include "Mesh_maps_stk.hh"
#include "MeshAudit.hh"
#include "cgns_mesh_par.hh"

// -------------------------------------------------------------
// dump_cgns
// -------------------------------------------------------------
void
dump_cgns(const int& me, Mesh_maps_base &maps, const std::string& cgnsout)
{
  CGNS_PAR::create_mesh_file(maps, cgnsout);

  Epetra_Vector part(maps.cell_map(false));
  int nmycell(maps.cell_map(false).NumMyElements());
  std::vector<int> myidx(nmycell, 0);
  for (unsigned int i = 0; i < nmycell; i++) myidx[i] = i;

  std::vector<double> mypart(nmycell, static_cast<double>(me+1.0));
  part.ReplaceMyValues(nmycell, &mypart[0], &myidx[0]);
  
  CGNS_PAR::open_data_file(cgnsout);
  CGNS_PAR::create_timestep(0.0, 0, Mesh_data::CELL);
  CGNS_PAR::write_field_data(part, "Partition");
  CGNS_PAR::close_data_file();
}

// -------------------------------------------------------------
// do_the_audit
// -------------------------------------------------------------
void
do_the_audit(const int& me, Teuchos::RCP<Mesh_maps_base> maps, const std::string& name)
{
  int lresult(0);

  std::ostringstream ofile;
  ofile << name << std::setfill('0') << std::setw(4) << me << ".txt";
  std::ofstream ofs(ofile.str().c_str());
  if (me == 0)
    std::cout << "Writing results to " << ofile.str() << ", etc." << std::endl;
  MeshAudit audit(maps, ofs);
  lresult = audit.Verify();

  int gresult;

  maps->get_comm()->MaxAll(&lresult, &gresult, 1);

  if (me == 0) {
    std::cout << "Mesh \"" << name << "\" " << (gresult ? "has errors" : "OK") << std::endl;
  }

}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
int
main(int argc, char **argv)
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  const int nproc(comm.NumProc());
  const int me(comm.MyPID());
  stk::ParallelMachine pm(comm.Comm());

  const unsigned int mult(nproc);
  const unsigned int isize(mult*4), jsize(4), ksize(4);

  STK_mesh::Mesh_maps_stk *maps_stk = 
    new STK_mesh::Mesh_maps_stk(comm, isize, jsize, ksize);
  Teuchos::RCP<Mesh_maps_base> maps(maps_stk);

  do_the_audit(me, maps, "stk_mesh_test_hex_before");

  dump_cgns(me, *maps, "stk_mesh_test_hex_before.cgns");

  // maps_stk->cellgraph()->Print(std::cout);
  // maps_stk->cell_map(false).Print(std::cout);

  Teuchos::ParameterList params;
  params.set("PARTITIONING_METHOD", "GRAPH");    
  Teuchos::ParameterList& sublist= params.sublist("Zoltan");
  // sublist.set("phg_output_level", "5"); 

  // Teuchos::RCP< const Epetra_CrsGraph> cgraph(maps_stk->cellgraph());
  // Isorropia::Epetra::Partitioner partitioner(cgraph, params, false);
  // partitioner.partition();
  // Teuchos::RCP< Epetra_Map > newcmap(partitioner.createNewMap()); // 0-based

  // newcmap->Print(std::cout);
  // maps_stk->redistribute(*newcmap);

  maps_stk->redistribute(params);

  // maps_stk->cell_map(false).Print(std::cout);

  do_the_audit(me, maps, "stk_mesh_test_hex_after");

  dump_cgns(me, *maps, "stk_mesh_test_hex_after.cgns");

  return 0;
}

