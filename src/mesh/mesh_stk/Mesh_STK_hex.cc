/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
// file: Mesh_STK.cc
// -------------------------------------------------------------
/**
 * @file   Mesh_STK.cc
 * @author William A. Perkins
 * @date Mon Aug  1 09:56:58 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// Created May  2, 2011 by William A. Perkins
// Last Change: Mon Aug  1 09:56:58 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include "Mesh_STK.hh"
#include "Mesh_STK_Impl.hh"
#include "HexMeshGenerator.hh"
#include "Mesh_STK_factory.hh"
#include "RectangularRegion.hh"

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
Mesh_STK::generate_(Teuchos::ParameterList &parameter_list)
{
  int nx, ny, nz;
  double x0, y0, z0;
  double x1, y1, z1;

  // read the parameters from the parameter list

  nx = parameter_list.get<int>("Numer of Cells in X");
  ny = parameter_list.get<int>("Numer of Cells in Y");
  nz = parameter_list.get<int>("Numer of Cells in Z");
  
  x0 = parameter_list.get<double>("X_Min", 0);
  x1 = parameter_list.get<double>("X_Max", 1);

  y0 = parameter_list.get<double>("Y_Min", 0);
  y1 = parameter_list.get<double>("Y_Max", 1);

  z0 = parameter_list.get<double>("Z_Min", 0);
  z1 = parameter_list.get<double>("Z_Max", 1);

  double xdelta((x1 - x0)/static_cast<double>(nx));
  double ydelta((y1 - y0)/static_cast<double>(ny));
  double zdelta((z1 - z0)/static_cast<double>(nz));

  Data::HexMeshGenerator g(communicator_.get(), 
                           nx, ny, nz,
                           x0, y0, z0,
                           xdelta, ydelta, zdelta);

  int nblk;

  nblk = parameter_list.get<int>("Number of mesh blocks",0);
  
  if (nblk > 0) {
    for (int nb = 1; nb <= nblk; nb++) {
      std::stringstream s; 
      s << "Mesh block " << nb;

      Teuchos::ParameterList sublist = parameter_list.sublist(s.str());

      // tell the generator about the zone

      AmanziGeometry::Point p0(x0, y0, sublist.get<double>("Z0"));
      AmanziGeometry::Point p1(x1, y1, sublist.get<double>("Z1"));
      AmanziGeometry::RegionPtr r(new AmanziGeometry::RectangularRegion(p0, p1));

      g.add_region(nb, r);
    }
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

Mesh_STK::Mesh_STK(Teuchos::ParameterList &parameter_list,
                   Epetra_MpiComm *communicator)
  : communicator_(new Epetra_MpiComm(*communicator)),
    mesh_(), 
    entity_map_ (3), 
    map_owned_(), map_used_()
{
  Mesh::set_comm(communicator_->GetMpiComm());
  generate_(parameter_list);
}

} // namespace AmanziMesh
} // namespace Amanzi
