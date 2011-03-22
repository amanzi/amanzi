// -------------------------------------------------------------
/**
 * @file   Mesh_maps_stk_hex.cc
 * @author William A. Perkins
 * @date Mon Mar 14 15:34:01 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December 13, 2010 by William A. Perkins
// Last Change: Mon Mar 14 15:34:01 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include "Mesh_maps_stk.hh"
#include "HexMeshGenerator.hh"
#include "Mesh_factory.hh"

namespace STK_mesh {

  // -------------------------------------------------------------
  // Mesh_maps_stk::generate_hex_
  // -------------------------------------------------------------
  void
  Mesh_maps_stk::generate_(const Epetra_MpiComm& comm, 
                           const unsigned int& ni, const unsigned int& nj, const unsigned int& nk,
                           const double& xorigin, 
                           const double& yorigin, 
                           const double& zorigin, 
                           const double& xdelta, 
                           const double& ydelta, 
                           const double& zdelta)
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
    generate_(comm, ni, nj, nk, xorigin, yorigin, zorigin, xdelta, ydelta, zdelta);
  }

  Mesh_maps_stk::Mesh_maps_stk(const double x0, const double y0, const double z0,
                               const double x1, const double y1, const double z1,
                               const int nx, const int ny, const int nz, 
                               Epetra_MpiComm *comm)
    : mesh_(), 
      entity_map_ (3),
      communicator_(comm->Clone())
  {
    double xdelta((x1 - x0)/static_cast<double>(nx));
    double ydelta((y1 - y0)/static_cast<double>(ny));
    double zdelta((z1 - z0)/static_cast<double>(nz));
    generate_(*comm, 
              static_cast<unsigned int>(nx), 
              static_cast<unsigned int>(ny), 
              static_cast<unsigned int>(nz), 
              x0, y0, z0, xdelta, ydelta, zdelta);
  }
  
} // close namespace STK_mesh
