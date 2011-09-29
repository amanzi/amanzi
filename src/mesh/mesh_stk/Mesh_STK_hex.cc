/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
// file: Mesh_STK.cc
// -------------------------------------------------------------
/**
 * @file   Mesh_STK.cc
 * @author William A. Perkins
 * @date Wed Sep 28 10:32:43 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// Created May  2, 2011 by William A. Perkins
// Last Change: Wed Sep 28 10:32:43 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include "Mesh_STK.hh"
#include "Mesh_STK_Impl.hh"
#include "HexMeshGenerator.hh"
#include "Mesh_STK_factory.hh"
#include "RectangularRegion.hh"
#include "GenerationSpec.hh"

namespace Amanzi {
namespace AmanziMesh {

// -------------------------------------------------------------
//  class Mesh_STK
// -------------------------------------------------------------

// -------------------------------------------------------------
// Mesh_maps_stk::generate_
// -------------------------------------------------------------

void
Mesh_STK::generate_(Data::HexMeshGenerator& g)
{
  Teuchos::RCP<Data::Data> meshdata(g.generate());
  Teuchos::RCP<Epetra_Map> cmap(g.cellmap(true));
  Teuchos::RCP<Epetra_Map> vmap(g.vertexmap(true));
  stk::ParallelMachine pm(communicator_->GetMpiComm());
  STK::Mesh_STK_factory mf(pm, 1000);
  Data::Fields nofields;
  mesh_.reset(mf.build_mesh(*meshdata, *cmap, *vmap, nofields));
  build_maps_();
  redistribute();

  // FIXME: this is supposed to be temporary
  fill_setnameid_map_();
}

void
Mesh_STK::generate_(const unsigned int& ni, const unsigned int& nj, const unsigned int& nk,
                    const double& xorigin, 
                    const double& yorigin, 
                    const double& zorigin, 
                    const double& xdelta, 
                    const double& ydelta, 
                    const double& zdelta)
{
  Data::HexMeshGenerator g(communicator_.get(), 
                           ni, nj, nk,
                           xorigin, yorigin, zorigin,
                           xdelta, ydelta, zdelta);
  generate_(g);
}

void
Mesh_STK::generate_(const GenerationSpec& gspec)
{
  int nx(gspec.xcells());
  int ny(gspec.ycells());
  int nz(gspec.zcells());

  AmanziGeometry::Point p0(gspec.domain().point0());
  AmanziGeometry::Point p1(gspec.domain().point1());

  double xdelta((p1.x() - p0.x())/static_cast<double>(nx));
  double ydelta((p1.y() - p0.y())/static_cast<double>(ny));
  double zdelta((p1.z() - p0.z())/static_cast<double>(nz));

  Data::HexMeshGenerator g(communicator_.get(), 
                           nx, ny, nz,
                           p0.x(), p0.y(), p0.z(),
                           xdelta, ydelta, zdelta);

  int nblk(0);

  AmanziGeometry::RegionVector::const_iterator r;
  for (r = gspec.block_begin(); r != gspec.block_end(); ++r, nblk++) {
    std::stringstream s; 
    s << "Mesh block " << nblk + 1;

    g.add_region(nblk+1, s.str(), *r);
  }
  
  generate_(g);
}

// -------------------------------------------------------------
// Mesh_STK:: constructors / destructor
// -------------------------------------------------------------
Mesh_STK::Mesh_STK(const Epetra_MpiComm& comm, 
                        const unsigned int& ni, const unsigned int& nj, const unsigned int& nk,
                        const double& xorigin, 
                        const double& yorigin, 
                        const double& zorigin, 
                        const double& xdelta, 
                        const double& ydelta, 
                        const double& zdelta)
    : communicator_(new Epetra_MpiComm(comm)),
      mesh_(), 
      entity_map_ (3), 
      map_owned_(), map_used_()
      
{
  Mesh::set_comm(communicator_->GetMpiComm());
  generate_(ni, nj, nk, xorigin, yorigin, zorigin, xdelta, ydelta, zdelta);
}

Mesh_STK::Mesh_STK(const double x0, const double y0, const double z0,
                   const double x1, const double y1, const double z1,
                   const int nx, const int ny, const int nz, 
                   Epetra_MpiComm *comm)
    : communicator_(new Epetra_MpiComm(*comm)),
      mesh_(), 
      entity_map_ (3), 
      map_owned_(), map_used_()
      
{
  Mesh::set_comm(communicator_->GetMpiComm());
  double xdelta((x1 - x0)/static_cast<double>(nx));
  double ydelta((y1 - y0)/static_cast<double>(ny));
  double zdelta((z1 - z0)/static_cast<double>(nz));
  generate_(static_cast<unsigned int>(nx), 
            static_cast<unsigned int>(ny), 
            static_cast<unsigned int>(nz), 
            x0, y0, z0, xdelta, ydelta, zdelta);
}

Mesh_STK::Mesh_STK(const GenerationSpec& gspec,
                   Epetra_MpiComm *communicator)
  : communicator_(new Epetra_MpiComm(*communicator)),
    mesh_(), 
    entity_map_ (3), 
    map_owned_(), map_used_()
{
  Mesh::set_comm(communicator_->GetMpiComm());
  generate_(gspec);
}

} // namespace AmanziMesh
} // namespace Amanzi
