/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

*/

#include <vector>

// TPLs
#include "Epetra_Vector.h"

// Amanzi
#include "errors.hh"
#include "LinearOperator.hh"
#include "LinearOperatorFactory.hh"
#include "MatrixFE.hh"
#include "MFD3D_CrouzeixRaviart.hh"
#include "MFD3D_Diffusion.hh"
#include "PreconditionerFactory.hh"
#include "SuperMap.hh"
#include "WhetStoneDefs.hh"
#include "Mesh.hh"
#include "GMVMesh.hh"
// Operators
#include "Op.hh"
#include "Op_Cell_Edge.hh"
#include "Op_Cell_Node.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Face_Cell.hh"
#include "Op_SurfaceFace_SurfaceCell.hh"

#include "OperatorDefs.hh"
#include "Operator_Edge.hh"
#include "Operator_FaceCell.hh"
#include "Operator_FaceCellScc.hh"
#include "Operator_FaceCellSff.hh"
#include "Operator_Node.hh"

#include "PDE_DiffusionMFD_Tracer.hh"
//#include "xmof2D.h"

namespace Amanzi {
namespace Operators {

  void PDE_DiffusionMFD_Tracer::InitDiffusion_(Teuchos::ParameterList& plist){

    //PDE_DiffusionMFD::InitDiffusion_(plist);

    cell_surfaces_.resize(ncells_owned);

  }

  void PDE_DiffusionMFD_Tracer::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                               const Teuchos::Ptr<const CompositeVector>& u,
                                               const Teuchos::Ptr<const CompositeVector>& surf_presence,
                                               const Teuchos::Ptr<const CompositeVector>& surface_param){

    // update matrix blocks
    WhetStone::MFD3D_Diffusion mfd(mesh_);
    mfd.ModifyStabilityScalingFactor(factor_);

    AmanziMesh::Entity_ID_List nodes;

    const Epetra_MultiVector& surf_presence_vec = *surf_presence->ViewComponent("cell", true);
    const Epetra_MultiVector& surface_param_vec = *surface_param->ViewComponent("cell", true);

    int num_surfaces = surf_presence_vec.NumVectors();
    
    WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
    Kc(0, 0) = 1.0;

    for (int c = 0; c < ncells_owned; c++) {
    //for (int c = 0; c < 1; c++) {      
      if (K_.get()) Kc = (*K_)[c];
      
      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      WhetStone::DenseMatrix Acell(nnodes, nnodes);

      Acell = 0.;

      for (int s_id=0; s_id<num_surfaces;s_id++){
        if (surf_presence_vec[s_id][c] > 0){
        
          //AmanziGeometry::Point sur_cntr(3);//(0.5, 0.5, 0.5);
          AmanziGeometry::Point sur_norm(3);//(1.0, 0.6, 1.0);
          double surf_d = surface_param_vec[s_id*4+3][c];

          sur_norm[0] = surface_param_vec[s_id*4][c];
          sur_norm[1] = surface_param_vec[s_id*4+1][c];
          sur_norm[2] = surface_param_vec[s_id*4+2][c];


          //sur_norm[0]=1.; sur_norm[1]=0.6; sur_norm[2]=1.;  surf_d=-1.3;
          
          CellSurfaceInterception_(c, s_id, sur_norm, surf_d);
          std::cout<<"CELL: "<<c<<"\n";

          
          if (cell_surfaces_[c].size() > 0){
            if (cell_surfaces_[c][0]->num_faces() > 2){

              cell_surfaces_[c][0]->print();          
              mfd.StiffnessMatrixTracer(c, Kc, sur_norm,
                                        cell_surfaces_[c][0]->surface_pnt(),
                                        cell_surfaces_[c][0]->v_ids(),
                                        cell_surfaces_[c][0]->inter_coef(),
                                        Acell);
            }
          }
          std::cout<<Acell<<"\n";
          //if (c==12) exit(0);
          break;          
        }else{
          //for (int i=0;i<nnodes;i++) Acell(i,i)=1.;
        }
      }

      // std::cout<<"cell "<<c<<" nodes ";
      // for (int i=0;i<nnodes;i++)std::cout<<nodes[i]<<" "; std::cout<<"\n";
      // std::cout<<Acell<<"\n";      
      local_op_->matrices[c] = Acell;     
    }
    
  }

  void PDE_DiffusionMFD_Tracer::ApplyBCs(const Epetra_MultiVector& marker, bool primary, bool eliminate){

     Teuchos::Ptr<const BCs> bc_f, bc_n;
      for (const auto& bc : bcs_trial_){
        if (bc->kind() == AmanziMesh::FACE) {
          bc_f = bc.ptr();
        } else if (bc->kind() == AmanziMesh::NODE) {
          bc_n = bc.ptr();
        }
      }
      ApplyBCs_Nodal_(marker, bc_f.ptr(), bc_n.ptr(), primary, eliminate);
      
  }


/* ******************************************************************
* Apply BCs on nodal operators
****************************************************************** */
  void PDE_DiffusionMFD_Tracer::ApplyBCs_Nodal_(const Epetra_MultiVector& marker,
                                                const Teuchos::Ptr<const BCs>& bc_f,
                                                const Teuchos::Ptr<const BCs>& bc_v,
                                                bool primary, bool eliminate)
{
  AmanziMesh::Entity_ID_List faces, nodes, cells;

  global_op_->rhs()->PutScalarGhosted(0.0);
  Epetra_MultiVector& rhs_node = *global_op_->rhs()->ViewComponent("node", true);

  int nn(0), nm(0);
  for (int c = 0; c != ncells_owned; ++c) {
    bool flag(true);
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];
   
    // process boundary integrals
    // if (bc_f != Teuchos::null) {
    //   const std::vector<int>& bc_model = bc_f->bc_model();
    //   const std::vector<double>& bc_value = bc_f->bc_value();
    //   const std::vector<double>& bc_mixed = bc_f->bc_mixed();

    //   mesh_->cell_get_faces(c, &faces);
    //   int nfaces = faces.size();

    //   for (int n = 0; n != nfaces; ++n) {
    //     int f = faces[n];

    //     if (bc_model[f] == OPERATOR_BC_NEUMANN) {
    //       nn++;
    //       double value = bc_value[f];
    //       double area = mesh_->face_area(f);

    //       mesh_->face_get_nodes(f, &nodes);
    //       int nnodes = nodes.size();

    //       for (int m = 0; m < nnodes; m++) {
    //         int v = nodes[m];
    //         if (bc_v->bc_model()[v] != OPERATOR_BC_DIRICHLET)
    //           rhs_node[0][v] -= value * area / nnodes;
    //       }
    //     } else if (bc_model[f] == OPERATOR_BC_MIXED) {
    //       nm++;
    //       if (flag) {  // make a copy of cell-based matrix
    //         local_op_->matrices_shadow[c] = Acell;
    //         flag = false;
    //       }
    //       double value = bc_value[f];
    //       double area = mesh_->face_area(f);

    //       mesh_->face_get_nodes(f, &nodes);
    //       int nnodes = nodes.size();

    //       for (int m = 0; m < nnodes; m++) {
    //         int v = nodes[m];
    //         if (bc_v->bc_model()[v] != OPERATOR_BC_DIRICHLET)
    //           rhs_node[0][v] -= value * area / nnodes;
    //         Acell(n, n) += bc_mixed[f] * area / nnodes;
    //       }
    //     }
    //   }
    // } 


    
    if (bc_v != Teuchos::null) {
      const std::vector<int>& bc_model = bc_v->bc_model();
      const std::vector<double>& bc_value = bc_v->bc_value();
      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();
      
      std::cout<<"rhs: ";
      for (int n = 0; n != nnodes; ++n) {
        std::cout<<rhs_node[0][nodes[n]]<<" ";
      }
      std::cout<<"\n";
      std::cout<<"original "<<c<<"\n"<<Acell<<"\n";
      
      if (marker[0][c] > 0){
        // essential conditions for test functions
        for (int n = 0; n != nnodes; ++n) {
          int v = nodes[n];
          if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
            if (flag) {  // make a copy of elemental matrix
              local_op_->matrices_shadow[c] = Acell;
              flag = false;
            }
            for (int m = 0; m < nnodes; m++) Acell(n, m) = 0.0;
          }
        }

      

        for (int n = 0; n != nnodes; ++n) {
          int v = nodes[n];
          double value = bc_value[v];

          if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
            if (flag) {  // make a copy of cell-based matrix
              local_op_->matrices_shadow[c] = Acell;
              flag = false;
            }
     
            if (eliminate) {
              for (int m = 0; m < nnodes; m++) {
                rhs_node[0][nodes[m]] -= Acell(m, n) * value;
                Acell(m, n) = 0.0;
              }
            }

            // We take into account multiple contributions to matrix diagonal
            // by dividing by the number of cells attached to a vertex.
            if (primary) {
              mesh_->node_get_cells(v, AmanziMesh::USED, &cells);
              if (v < nnodes_owned) rhs_node[0][v] = value;
              Acell(n, n) = 1.0 / cells.size();
            }
          
          }
        }
      }else{
      
        // essential conditions for test functions
        for (int n = 0; n != nnodes; ++n) {

          int v = nodes[n];
          double value = bc_value[v];
          if (flag) {  // make a copy of elemental matrix
            local_op_->matrices_shadow[c] = Acell;
            flag = false;
          }
          mesh_->node_get_cells(v, AmanziMesh::USED, &cells);

          if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
            if (primary) {
              if (v < nnodes_owned) rhs_node[0][v] = value;
              Acell(n, n) = 1.0 / cells.size();
            }
          }else{
            // bool active(false);
            // for (int j=0; j<cells.size(); j++){
            //   if (marker[0][cells[j]] > 0){
            //     active = true;
            //     break;
            //   }
            // }
            // if (!active){
            //   Acell(n, n) = 1.;
            //   std::cout<<"non-active "<< n<<"\n";
            // }
          }
        }
      }
      std::cout<<"rhs: ";
      for (int n = 0; n != nnodes; ++n) {
        std::cout<<rhs_node[0][nodes[n]]<<" ";
      }
      std::cout<<"\n";
      std::cout<<"cell "<<c<<"\n"<<Acell<<"\n";


      //if (c==12) exit(0);
    }    
  } 

  global_op_->rhs()->GatherGhostedToMaster("node", Add);
}


  void PDE_DiffusionMFD_Tracer::CellSurfaceInterception_(int c, int surf_id,
                                                         const AmanziGeometry::Point& sur_norm, double surf_d){

    AmanziMesh::Entity_ID_List cells, faces, vertices, edges;
    std::vector<int> edges_dirs;
    int dir;
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    Teuchos::RCP<Local_Surface> local_surf = Teuchos::rcp(new Local_Surface(surf_id, sur_norm, surf_d));
    //surface_pnt_[c].resize(nfaces);

    for (int i=0; i<nfaces; i++){
      int f = faces[i];
      AmanziGeometry::Point face_normal = mesh_->face_normal(f, false, c, &dir);
      AmanziGeometry::Point face_centr = mesh_->face_centroid(f);
      AmanziGeometry::Point line_vec = sur_norm^face_normal;
      line_vec *= 1./norm(line_vec);
      int coord0 = 2;
      double eps = 1e-12;
      for (coord0=2; coord0>=0; --coord0){
        if (std::abs(line_vec[coord0]) > eps) break;
      }
      AmanziGeometry::Point line_point(3);

      int coord1, coord2;
      coord1 = (coord0+1)%3;
      coord2 = (coord0+2)%3;
      double a1,a2;
      a1 = sur_norm[coord1]; a2=face_normal[coord1];     
      double b1,b2;
      b1 = sur_norm[coord2]; b2=face_normal[coord2];     
      double d1,d2;
      d1 = surf_d;
      d2 = -face_normal*face_centr;

      double det = a1*b2 - a2*b1;
      line_point[coord0] = 0.;
      line_point[coord1] = (b1*d2 - d1*b2)/ det;
      line_point[coord2] = (a2*d1 - d2*a1)/ det;

      //mesh_->face_get_nodes(f, &vertices);
      mesh_->face_get_edges_and_dirs(f, &edges, &edges_dirs);

      //int nedges = vertices.size();
      int nedges = edges.size();      
      int num_inter = 0;
      AmanziGeometry::Point p_tmp(line_point + line_vec);
      
      std::vector<AmanziGeometry::Point> pnts;
      std::vector<int> vrt_ids;
      std::vector<double> wgt;
      std::vector<int> edges_ids;
      
      AmanziGeometry::Point pa(3), pb(3);
      if (c==715)
        std::cout<<"\n cell "<<c<<"face "<<f<<"\n";
      
      for (int j=0; j<nedges; j++){
        int v_id[2];
        // v_id[0] = vertices[j];
        // v_id[1] = vertices[(j+1)%nedges];
        int e = edges[j];
        mesh_->edge_get_nodes(e, v_id, v_id+1);
        AmanziGeometry::Point vs[2];
        for (int k=0; k<2; k++) mesh_->node_get_coordinates(v_id[k], &vs[k]);
        std::cout << "v_id "<<v_id[0]<<" "<<v_id[1]<<"\n";
        double t1, t2, eps=1e-7;
        bool intersection = LineLineIntersect_(vs[0], vs[1], line_point, p_tmp, eps, pa, pb, &t1, &t2);
        std::cout<<"t1 = "<<t1<<"\n";
        if (intersection){   // if possible to find shortest distance
          if ((t1 > -1e-7)&&(t1<1 + 1e-7)){
            double dist = norm(pa - pb);
            if (dist < eps){
              std::cout << std::setw(5)<< pa.x()<<" "<<std::setw(5)<<pa.y()<<" "<<std::setw(5)<<pa.z()<<"\n";
              pnts.push_back(pa);
              vrt_ids.push_back(v_id[0]);
              vrt_ids.push_back(v_id[1]);
              wgt.push_back(t1);
              edges_ids.push_back(e);
              num_inter++;
            }
          }
        }

      }
      if (num_inter > 2) {
       
        auto it_vrt0 = vrt_ids.begin();
        auto it_wgt0 = wgt.begin();
        auto it_edges_ids0 = edges_ids.begin();
        for (auto it_i=pnts.begin(); it_i != pnts.end(); ++it_i){
          pa = *it_i;
          auto it_vrt1 = it_vrt0 + 2;
          auto it_wgt1 = it_wgt0 + 1;
          auto it_edges_ids1 = it_edges_ids0 + 1; 
          for (auto it_j= it_i + 1; it_j != pnts.end(); ++it_j){
            pb = *it_j;
            if (norm(pa - pb)<1e-6){
              pnts.erase(it_j);
              vrt_ids.erase(it_vrt1);
              vrt_ids.erase(it_vrt1);
              wgt.erase(it_wgt1);
              edges_ids.erase(it_edges_ids1);
            }else{
              it_vrt1 += 2;
              it_wgt1 ++;
              it_edges_ids1 ++;
            }
            if (it_j == pnts.end()) break;
          }
          it_vrt0 += 2;
          it_wgt0 ++;
          it_edges_ids0 ++;
          if (it_i == pnts.end()) break;
        }
        if (pnts.size() > 2){
          Errors::Message msg;
          msg << "Surface intersectes a face in more then 2 points\n";
          Exceptions::amanzi_throw(msg);
        }
      }
      // if (num_inter == 1) {
      //   pnts.erase(it_j);
      //   vrt_ids.erase(it_vrt1);
      //   vrt_ids.erase(it_vrt1+1);
      //   wgt.erase(it_wgt1);
      //   edges_ids.erase(it_edges_ids1);
      //   inter_coef_.erase(it_inter_coef1);
      // } 
      if ( pnts.size() == 2){
        if (norm(pnts[1] - pnts[0]) > 1e-6){
          if (c == 715) std::cout << "\n poly face "<<vrt_ids[0] <<" "<<vrt_ids[1]<<" "<<vrt_ids[2] <<" "<<vrt_ids[3]<<"\n\n";
          local_surf->Add_Face(pnts, vrt_ids, wgt, edges_ids);
        }
      }
      //std::cout<<"\n";
    }

    local_surf->print();

    if (local_surf->num_faces() > 2) 
      cell_surfaces_[c].push_back(local_surf);

    //if (c == 715) exit(0);
      
  }

  /*
   Calculate the shortest segment  between
   two lines P1P2 and P3P4. Calculate also the values of mu_1 and mu_2 where
       int_p1 = P1 + mu_1 (P2 - P1)
       int_p2 = P3 + mu_2 (P4 - P3)
   Return FALSE if no solution exists.
*/

  bool PDE_DiffusionMFD_Tracer::LineLineIntersect_(const AmanziGeometry::Point& p1,
                                                  const AmanziGeometry::Point& p2,
                                                  const AmanziGeometry::Point& p3,
                                                  const AmanziGeometry::Point& p4,
                                                  double EPS,
                                                  AmanziGeometry::Point& int_p1,
                                                  AmanziGeometry::Point& int_p2,
                                                  double *mu_1, double *mu_2){

    AmanziGeometry::Point p13,p43,p21;
    double d1343,d4321,d1321,d4343,d2121;
    double numer,denom;

    p13[0] = p1[0] - p3[0];
    p13[1] = p1[1] - p3[1];
    p13[2] = p1[2] - p3[2];
    p43[0] = p4[0] - p3[0];
    p43[1] = p4[1] - p3[1];
    p43[2] = p4[2] - p3[2];
    if (std::abs(p43[0]) < EPS && std::abs(p43[1]) < EPS && std::abs(p43[2]) < EPS)
      return(false);
    p21[0] = p2[0] - p1[0];
    p21[1] = p2[1] - p1[1];
    p21[2] = p2[2] - p1[2];
    if (std::abs(p21[0]) < EPS && std::abs(p21[1]) < EPS && std::abs(p21[2]) < EPS)
      return(false);

    d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2];
    d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2];
    d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2];
    d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2];
    d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2];

    denom = d2121 * d4343 - d4321 * d4321;
    if (std::abs(denom) < EPS)
      return(false);
    numer = d1343 * d4321 - d1321 * d4343;

    *mu_1 = numer / denom;
    *mu_2 = (d1343 + d4321 * (*mu_1)) / d4343;

    int_p1[0] = p1[0] + *mu_1 * p21[0];
    int_p1[1] = p1[1] + *mu_1 * p21[1];
    int_p1[2] = p1[2] + *mu_1 * p21[2];
    int_p2[0] = p3[0] + *mu_2 * p43[0];
    int_p2[1] = p3[1] + *mu_2 * p43[1];
    int_p2[2] = p3[2] + *mu_2 * p43[2];

    return(true);
    
    
  }

  void PDE_DiffusionMFD_Tracer::OutputGMV_Surface(int surf_id, const Epetra_MultiVector& data, std::string filename){

    gmvwrite_openfile_ir_ascii((char*)filename.c_str(), 4, 8);

    unsigned int num_edges = mesh_->num_entities(AmanziMesh::EDGE, AmanziMesh::OWNED);
    unsigned int num_cells_bulk = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

    std::vector<int> edges2node(num_edges, -1);
    std::vector<AmanziGeometry::Point> nodes;
    std::vector<double> node_data, cell_data;
    std::vector< std::vector<int> > surface_cells;
    int node_id = 0;
    int max_nodes = 0;
    for (int c=0;c<num_cells_bulk; c++){
      int k = -1;
      for (int i=0; i<cell_surfaces_[c].size(); i++){
        // find a local surface with given surf_id
        if (cell_surfaces_[c][i]->surface_id()==surf_id){
          k = i;
          break;
        }
      }
      if (k==-1) continue; // no local surface with given surf_id
      std::vector< std::vector<int> >& v_ids = cell_surfaces_[c][k]->v_ids();
      std::vector< std::vector<int> >& edges_ids = cell_surfaces_[c][k]->edge_ids();
      std::vector< std::vector<double> >& inter_coef = cell_surfaces_[c][k]->inter_coef();
      std::vector< std::vector<AmanziGeometry::Point> >& surf_pnt = cell_surfaces_[c][k]->surface_pnt();
      std::vector<int> poly;

      int num_points = 1;
      int nfaces = inter_coef.size();

      cell_data.push_back(c);
        
      int e1 = edges_ids[0][0];
      if (edges2node[e1]<0) { // add new node to node list
        edges2node[e1] = node_id++;
        nodes.push_back(surf_pnt[0][0]);
        double node_val;
        int v1_id = v_ids[0][0];
        int v2_id = v_ids[0][1];
        node_val = (1 - inter_coef[0][0])*data[surf_id][v1_id] + inter_coef[0][0]*data[surf_id][v2_id];
        node_data.push_back(node_val);
      }
      poly.push_back(edges2node[e1]); // node to polygon list
      AmanziGeometry::Point p1(3),p2(3);
      int e2 = edges_ids[0][1]; // next node
      p2 = surf_pnt[0][1];
      int i_prev = 0;
      
      int l,m;
      double eps = 1e-6;
      while (num_points < nfaces){
        bool find = false;
        for (int i=0; i<nfaces; i++){
          if (i==i_prev) continue;
          if (edges_ids[i][1]==e2){
            l=i;m = 1;
            find = true;
            break;
          }else if (edges_ids[i][0]==e2){
            l=i; m=0;
            find = true;
            break;
          }
        }
        if (!find){
          for (int i=0; i<nfaces; i++){
            if (i==i_prev) continue;
            if (norm(p2 - surf_pnt[i][1]) < eps){
              l=i; m = 1;
              find = true;
              break;
            }else if ( norm(p2 - surf_pnt[i][0]) <eps ){
              l=i; m=0;
              find = true;
              break;
            }
          }
        }

        if (!find){
          AmanziMesh::Entity_ID_List nodes;
          AmanziGeometry::Point vs;
          mesh_->cell_get_nodes(c, &nodes);
          std::cout<<"cell "<<c<<"\n";
          for (int i1=0; i1<nodes.size(); ++i1){
            mesh_->node_get_coordinates(nodes[i1], &vs);
            std::cout<<nodes[i1]<<": "<<vs<<"\n";
          } 
          cell_surfaces_[c][k]->print();
          Errors::Message msg;
          msg << "Can't build polygon in OutputGMV_Surface\n";
          Exceptions::amanzi_throw(msg);
        }
        e1 = edges_ids[l][m];
        if (edges2node[e1]<0) {
          edges2node[e1] = node_id++;
          nodes.push_back(surf_pnt[l][m]);
          double node_val;
          int v1_id = v_ids[l][2*m];
          int v2_id = v_ids[l][2*m+1];
          node_val = (1 - inter_coef[l][m])*data[surf_id][v1_id] + inter_coef[l][m]*data[surf_id][v2_id];
          node_data.push_back(node_val);
        }
        if (e2!=e1) edges2node[e2] = edges2node[e1];
        poly.push_back(edges2node[e1]);
        num_points++;
          
        e2 = edges_ids[l][1-m];
        p2 = surf_pnt[l][1-m];
        i_prev = l;
      }

      max_nodes = std::max(max_nodes, num_points);
      //      std::cout<<"poly cell ";for (int i=0;i<poly.size();i++) std::cout<<poly[i]<<" "; std::cout<<"\n";
      surface_cells.push_back(poly);
    }

    int num_nodes = nodes.size();
    double *x = new double[num_nodes];
    double *y = new double[num_nodes];
    double *z = new double[num_nodes];
    for (int i=0; i<num_nodes; i++){
      x[i] = nodes[i][0];
      y[i] = nodes[i][1];
      z[i] = nodes[i][2];
    }
    gmvwrite_node_data(&num_nodes, x, y, z);
    
    delete [] z;
    delete [] y;
    delete [] x;

    // Write cell info
    unsigned int num_cells = surface_cells.size();
    
    gmvwrite_cell_header(&num_cells);
    unsigned int *xh = new unsigned int[max_nodes];

    for (int i=0; i<num_cells; i++) {
      int nnodes = surface_cells[i].size();
      for (int j=0; j<nnodes; j++) xh[j] = surface_cells[i][j] + 1;
      
      gmvwrite_cell_type((char*) "general 1", nnodes, xh);
    }
    delete [] xh;
    
    gmvwrite_variable_header();
    std::string varname="solution";
    gmvwrite_variable_name_data(1, (char*) varname.c_str(), node_data.data());
    varname="cell_id";
    gmvwrite_variable_name_data(0, (char*) varname.c_str(), cell_data.data());
    gmvwrite_variable_endvars();
    gmvwrite_closefile();

    //GMV::create_mesh_file(*mesh_,  (std::string)"bulk.gmv");

  }

}
}

