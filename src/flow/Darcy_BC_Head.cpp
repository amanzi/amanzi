/*
This is the flow component of the Amanzi code.

 
Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Flow_constants.hpp"
#include "Darcy_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {


/* ******************************************************************
* Process parameter for special treatment of static head b.c.                                           
****************************************************************** */
void Darcy_PK::ProcessShiftWaterTableList()
{
  std::string name("relative position of water table");

  if (dp_list_.isParameter(name)) {
    std::vector<std::string> regions;
    if (dp_list_.isType<Teuchos::Array<std::string> >(name)) {
      regions = dp_list_.get<Teuchos::Array<std::string> >(name).toVector();

      const Epetra_BlockMap& fmap = mesh_->face_map(false);
      shift_water_table_ = Teuchos::rcp(new Epetra_Vector(fmap));
      
      int nregions = regions.size();
      for (int i = 0; i < nregions; i++) {
        CalculateShiftWaterTable(regions[i]);
      }
    }
  }
}


/* ******************************************************************
* Calculate distance to the top of a given surface where the water 
* table is set up. 
* WARNING: works only in 3D.                                            
****************************************************************** */
void Darcy_PK::CalculateShiftWaterTable(const std::string region) 
{   
  double tol = 1e-24;
  Errors::Message msg;

  if (dim == 2) {
    msg << "Darcy PK: This boundary condition is not supported in 2D.";
    Exceptions::amanzi_throw(msg);
  }

  AmanziMesh::Entity_ID_List cells, faces, ss_faces;
  AmanziMesh::Entity_ID_List nodes1, nodes2, common_nodes;

  AmanziGeometry::Point p1(dim), p2(dim), p3(dim);
  std::vector<AmanziGeometry::Point> edges;

  mesh_->get_set_entities(region, AmanziMesh::FACE, AmanziMesh::OWNED, &ss_faces);
  int n = ss_faces.size();

  for (int i = 0; i < n; i++) {
    int f1 = ss_faces[i];
    mesh_->face_get_cells(f1, AmanziMesh::USED, &cells);

    mesh_->face_get_nodes(f1, &nodes1);
    std::sort(nodes1.begin(), nodes1.end());

    int c = cells[0];
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    // find all edges that intersection of boundary faces f1 and f2
    for (int j = 0; j < nfaces; j++) {
      int f2 = faces[j];
      if (f2 != f1) {
        mesh_->face_get_cells(f2, AmanziMesh::USED, &cells);
        int ncells = cells.size();
        if (ncells == 1) {
          mesh_->face_get_nodes(f2, &nodes2);
          std::sort(nodes2.begin(), nodes2.end());
          set_intersection(nodes1, nodes2, &common_nodes);

          int m = common_nodes.size();
          if (m > dim-1) {
            msg << "Darcy PK: Unsupported configuration: two or more common edges.";
            Exceptions::amanzi_throw(msg);
          } else if (m == 1 && dim == 2) {
            int v1 = common_nodes[0];
            mesh_->node_get_coordinates(v1, &p1);
            edges.push_back(p1);
          } else if (m == 2 && dim == 3) {
            int v1 = common_nodes[0], v2 = common_nodes[1];
            mesh_->node_get_coordinates(v1, &p1);
            mesh_->node_get_coordinates(v2, &p2);

            p3 = p1 - p2;
            if (p3[0] * p3[0] + p3[1] * p3[1] > tol * L22(p3)) {  // filter out vertical edges
              edges.push_back(p1);
              edges.push_back(p2);
            }
          }
        }
      }
    }
  }
  int nedges = edges.size();

#ifdef HAVE_MPI
  int gsize;
  const MPI_Comm& comm = mesh_->get_comm()->Comm();
  MPI_Comm_size(comm, &gsize);
  int* edge_counts = new int[gsize];
  MPI_Allgather(&nedges, 1, MPI_INT, edge_counts, 1, MPI_INT, comm);

  // prepare send buffer
  int sendcount = nedges * dim;
  double* sendbuf = NULL;
  if (nedges > 0) sendbuf = new double[sendcount];
  for (int i = 0; i < nedges; i++) {
    for (int k = 0; k < dim; k++) sendbuf[dim * i + k] = edges[i][k];
  } 

  // prepare receive buffer
  for (int i = 0; i < gsize; i++) edge_counts[i] *= dim;
  int recvcount = 0;
  for (int i = 0; i < gsize; i++) recvcount += edge_counts[i];
  double* recvbuf = new double[recvcount];

  int* displs = new int[gsize];
  displs[0] = 0;
  for (int i = 1; i < gsize; i++) displs[i] = edge_counts[i-1] + displs[i-1];

  MPI_Allgatherv(sendbuf, sendcount, MPI_DOUBLE, 
                 recvbuf, edge_counts, displs, MPI_DOUBLE, comm);

  // process receive buffer
  edges.clear();
  nedges = recvcount / dim;
  for (int i = 0; i < nedges; i++) {
    for (int k = 0; k < dim; k++) p1[k] = recvbuf[dim * i + k];
    edges.push_back(p1);
  } 
#endif
  // if (MyPID == 2) for (int i = 0; i < nedges; i++) 
  // printf("i= %5d  x = %12.6f %12.6f  %12.6f\n", i, edges[i][0], edges[i][1], edges[i][2]);

  // calculate head shift
  double edge_length, tol_edge, a, b;
  double rho_g = -rho_ * gravity_[dim - 1];

  for (int i = 0; i < n; i++) {
    int f = ss_faces[i];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    int flag = 0;
    for (int j = 0; j < nedges; j += 2) {
      p1 = edges[j + 1] - edges[j];
      p2 = xf - edges[j];

      edge_length = L22(p1);
      tol_edge = tol * edge_length;

      a = (p1 * p2) / edge_length;
      b = p1[0] * p2[1] - p1[1] * p2[0];
      if (b < tol_edge && a > -tol_edge && a < 1.0 + tol_edge) {
        double z = edges[j][dim - 1] - a * p1[dim - 1];  
        (*shift_water_table_)[f] = z * rho_g;
        if (z > xf[2]) {
          flag = 1;
          break;
        }
      }
    } 
    if (flag == 0) {
      msg << "Darcy PK: The boundary region \"" << region.c_str() << "\" is not piecewise flat.";
      Exceptions::amanzi_throw(msg);
    }
  }

#ifdef HAVE_MPI
  delete [] edge_counts;
  delete [] displs;
  delete [] recvbuf;
  if (sendbuf != NULL) delete [] sendbuf;
#endif

  if (MyPID == 0 && verbosity >= FLOW_VERBOSITY_HIGH) {
    printf("Darcy PK: found %5d boundary edges for region \"%s\"\n", nedges/2, region.c_str());   
  }
}


/* ******************************************************************
* New implementation of the STL function.                                              
****************************************************************** */
void Darcy_PK::set_intersection(const std::vector<int>& v1, 
                                const std::vector<int>& v2, std::vector<int>* vv)
{
  int i(0), j(0), n1, n2;

  n1 = (int)v1.size(); 
  n2 = (int)v2.size();
  vv->clear();

  while (i < n1 && j < n2) {
    if (v1[i] < v2[j]) {
      i++;
    } else if (v2[j] < v1[i]) {
      j++;
    } else {
      vv->push_back(v1[i]);
      i++;
      j++;
    }
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi


