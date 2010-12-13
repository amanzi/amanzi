// -------------------------------------------------------------
/**
 * @file   Mesh_maps_stk_hex.cc
 * @author William A. Perkins
 * @date Mon Dec 13 10:06:17 2010
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December 13, 2010 by William A. Perkins
// Last Change: Mon Dec 13 10:06:17 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include "Mesh_maps_stk.hh"
#include "HexMeshGenerator.hh"
#include "Mesh_factory.hh"

namespace STK_mesh {

  // -------------------------------------------------------------
  // Mesh_maps_stk::Mesh_maps_stk
  // -------------------------------------------------------------
  Mesh_maps_stk::Mesh_maps_stk(const Epetra_MpiComm& comm, 
                               const unsigned int& ni, const unsigned int& nj, const unsigned int& nk,
                               const double& xorigin, 
                               const double& yorigin, 
                               const double& zorigin, 
                               const double& xdelta, 
                               const double& ydelta, 
                               const double& zdelta)
    : mesh_(), 
      entity_map_ (3),
      communicator_(comm.Clone())
  {
    Mesh_data::HexMeshGenerator g(communicator_, 
                                  ni, nj, nk,
                                  xorigin, yorigin, zorigin,
                                  xdelta, ydelta, zdelta);
    Teuchos::RCP<Mesh_data::Data> meshdata(g.generate());
    Teuchos::RCP<Epetra_Map> cmap(g.cellmap(true));
    Teuchos::RCP<Epetra_Map> vmap(g.vertexmap(true));
    stk::ParallelMachine pm(comm.Comm()); // FIXME
    STK_mesh::Mesh_factory mf(pm, 1000);
    Mesh_data::Fields nofields;
    mesh_.reset(mf.build_mesh(*meshdata, *cmap, *vmap, nofields));
    build_maps_();
  }
  
} // close namespace STK_mesh
