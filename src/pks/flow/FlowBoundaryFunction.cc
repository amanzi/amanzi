/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1)
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#include "FlowBoundaryFunction.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Constructor
****************************************************************** */
FlowBoundaryFunction::FlowBoundaryFunction(const Teuchos::ParameterList& plist)
{
  rainfall_ = false;
  if (plist.isParameter("rainfall")) rainfall_ = plist.get<bool>("rainfall");

  relative_to_top_ = false;
  if (plist.isParameter("relative to top")) relative_to_top_ = plist.get<bool>("relative to top");

  relative_to_bottom_ = false;
  if (plist.isParameter("relative to bottom"))
    relative_to_bottom_ = plist.get<bool>("relative to bottom");

  no_flow_above_water_table_ = false;
  if (plist.isParameter("no flow above water table"))
    no_flow_above_water_table_ = plist.get<bool>("no flow above water table");

  if (plist.isParameter("submodel")) seepage_model_ = plist.get<std::string>("submodel");

  seepage_flux_threshold_ = 0.0;
  if (plist.isParameter("seepage flux threshold"))
    seepage_flux_threshold_ = plist.get<double>("seepage flux threshold");

  if (relative_to_top_ || relative_to_bottom_) {
    regions_ = plist.get<Teuchos::Array<std::string>>("regions").toVector();
    rho_ = plist.sublist("static head").sublist("function-static-head").get<double>("density");
    g_ = plist.sublist("static head").sublist("function-static-head").get<double>("gravity");
  }

  ref_pressure_ = FLOW_PRESSURE_ATMOSPHERIC;
  if (plist.isParameter("reference pressure"))
    ref_pressure_ = plist.get<double>("reference pressure");

  // for screen output
  nedges_ = 0;
}


/* ****************************************************************
* Process additional parameters for BC submodels.
**************************************************************** */
void
FlowBoundaryFunction::ComputeSubmodel(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  int dim = mesh->getSpaceDimension();

  if (rainfall_) {
    for (auto it = begin(); it != end(); ++it) {
      int f = it->first;
      const AmanziGeometry::Point& normal = mesh->getFaceNormal(f);
      it->second[0] *= fabs(normal[dim - 1]) / norm(normal);
    }
  }

  if (relative_to_top_ || relative_to_bottom_) {
    Teuchos::ParameterList vlist;
    vlist.sublist("verbose object");
    Teuchos::RCP<VerboseObject> vo = Teuchos::rcp(new VerboseObject("FlowBoundaryFunction", vlist));
    Teuchos::OSTab tab = vo->getOSTab();

    // lazy allocation of memory
    if (shift_water_table_ == Teuchos::null) {
      const Epetra_BlockMap& fmap = mesh->getMap(AmanziMesh::Entity_kind::FACE, true);
      shift_water_table_ = Teuchos::rcp(new Epetra_Vector(fmap));

      for (int i = 0; i < regions_.size(); ++i) {
        CalculateShiftWaterTable_(mesh, regions_[i]);

        if (vo->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
          *vo->os() << "found " << nedges_ << " top/bottom boundary edges for \"" << regions_[i]
                    << "\"" << std::endl;
        }
      }
    }

    for (auto it = begin(); it != end(); ++it) {
      int f = it->first;
      it->second[0] += (*shift_water_table_)[f];
    }
  }
}


/* ****************************************************************
* Calculate distance to the top of a given surface where the water
* table is set up. We do not distribute computed data. It seems not
* necessary.
*
* WARNING: The implemented algorithm works only in 3D.
**************************************************************** */
void
FlowBoundaryFunction::CalculateShiftWaterTable_(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                                const std::string& region)
{
  double tol = 1e-6;
  Errors::Message msg;

  int dim = mesh->getSpaceDimension();
  if (dim == 2) {
    msg << "Flow PK: \"relative/absolute\" action on static head BC is not supported in 2D.\n";
    Exceptions::amanzi_throw(msg);
  }

  AmanziMesh::Entity_ID_List cells, ss_faces;
  AmanziMesh::Entity_ID_List nodes1, nodes2, common_nodes;

  AmanziGeometry::Point p1(dim), p2(dim), p3(dim);
  std::vector<AmanziGeometry::Point> edges;

  mesh->getSetEntities(
    region, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED, &ss_faces);
  int n = ss_faces.size();

  for (int i = 0; i < n; i++) {
    int f1 = ss_faces[i];
    cells = mesh->getFaceCells(f1, AmanziMesh::Parallel_type::ALL);

    nodes1 = mesh->getFaceNodes(f1);
    std::sort(nodes1.begin(), nodes1.end());

    int c = cells[0];
    const auto& faces = mesh->getCellFaces(c);
    int nfaces = faces.size();

    // find all edges that intersection of boundary faces f1 and f2
    for (int j = 0; j < nfaces; j++) {
      int f2 = faces[j];
      if (f2 != f1) {
        cells = mesh->getFaceCells(f2, AmanziMesh::Parallel_type::ALL);
        int ncells = cells.size();
        if (ncells == 1) {
          nodes2 = mesh->getFaceNodes(f2);
          std::sort(nodes2.begin(), nodes2.end());
          set_intersection_(nodes1, nodes2, &common_nodes);

          int m = common_nodes.size();
          if (m > dim - 1) {
            msg << "Flow PK: Unsupported configuration: two or more common edges.";
            Exceptions::amanzi_throw(msg);
          } else if (m == 1 && dim == 2) {
            int v1 = common_nodes[0];
            p1 = mesh->getNodeCoordinate(v1);
            edges.push_back(p1);
          } else if (m == 2 && dim == 3) {
            int v1 = common_nodes[0], v2 = common_nodes[1];
            p1 = mesh->getNodeCoordinate(v1);
            p2 = mesh->getNodeCoordinate(v2);

            p3 = p1 - p2;
            if (p3[0] * p3[0] + p3[1] * p3[1] > tol * L22(p3)) { // filter out vertical edges
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
  Teuchos::RCP<const Comm_type> comm_p = mesh->getComm();
  Teuchos::RCP<const MpiComm_type> mpi_comm_p =
    Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_p);
  const MPI_Comm& comm = mpi_comm_p->Comm();
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
  for (int i = 1; i < gsize; i++) displs[i] = edge_counts[i - 1] + displs[i - 1];

  MPI_Allgatherv(sendbuf, sendcount, MPI_DOUBLE, recvbuf, edge_counts, displs, MPI_DOUBLE, comm);

  // process receive buffer
  edges.clear();
  nedges = recvcount / dim;
  for (int i = 0; i < nedges; i++) {
    for (int k = 0; k < dim; k++) p1[k] = recvbuf[dim * i + k];
    edges.push_back(p1);
  }
#endif

  // calculate head shift
  double edge_length, tol_edge, a, b;
  double rho_g = rho_ * g_;

  for (int i = 0; i < n; i++) {
    int f = ss_faces[i];
    const AmanziGeometry::Point& xf = mesh->getFaceCentroid(f);

    // fast algorithm for flat faces
    int flag = 0;
    for (int j = 0; j < nedges; j += 2) {
      p1 = edges[j + 1] - edges[j];
      p2 = xf - edges[j];

      edge_length = p1[0] * p1[0] + p1[1] * p1[1];
      tol_edge = tol * edge_length;

      a = (p1[0] * p2[0] + p1[1] * p2[1]) / edge_length;
      b = fabs(p1[0] * p2[1] - p1[1] * p2[0]);

      if (b < tol_edge && a > -0.01 && a < 1.01) {
        double z = edges[j][2] + a * p1[2];
        (*shift_water_table_)[f] = z * rho_g;
        if ((z > xf[2] && relative_to_top_) || (z < xf[2] && relative_to_bottom_)) {
          flag = 1;
          break;
        }
      }
    }

    // slow full search
    if (flag == 0) {
      // msg << "Flow PK: The boundary region \"" << region.c_str() << "\" is not piecewise flat.";
      // Exceptions::amanzi_throw(msg);
      // Instead, we take the closest mid-edge point with a higher z-coordinate.
      double z, d, dmin = 1e+99;
      for (int j = 0; j < nedges; j += 2) {
        p1 = (edges[j] + edges[j + 1]) / 2;
        d = L22(p1 - xf);
        if (((p1[2] > xf[2] && relative_to_top_) || (p1[2] < xf[2] && relative_to_bottom_)) &&
            d < dmin) {
          dmin = d;
          z = p1[2];
        }
      }
      (*shift_water_table_)[f] = z * rho_g;
    }
  }

#ifdef HAVE_MPI
  delete[] edge_counts;
  delete[] displs;
  delete[] recvbuf;
  if (sendbuf != NULL) delete[] sendbuf;
#endif

  nedges_ = nedges / 2;
}


/* ****************************************************************
* New implementation of the STL function.
**************************************************************** */
void
FlowBoundaryFunction::set_intersection_(const std::vector<AmanziMesh::Entity_ID>& v1,
                                        const std::vector<AmanziMesh::Entity_ID>& v2,
                                        std::vector<AmanziMesh::Entity_ID>* vv)
{
  int i(0), j(0), n1, n2;

  n1 = v1.size();
  n2 = v2.size();
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

} // namespace Flow
} // namespace Amanzi
