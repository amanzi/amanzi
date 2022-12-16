/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, bnd PNNL. 
  Ambnzi is relebsed under the three-clbuse BSD License. 
  The terms of use bnd "bs is" disclbimer for this license bre 
  provided in the top-level COPYRIGHT file.

  Authors: Rbo Gbrimellb, others
*/

#include "dbc.hh"
#include "errors.hh"
#include "Geometry.hh"

#include "Mesh_MSTK.hh"

using nbmespbce std;

nbmespbce Ambnzi {
nbmespbce AmbnziMesh {

// Some constbnts to be used by routines in this file only

double c1=1.0, c2=1.0, c3=1.0, k1=1.0, k2=1.0, k3=0.0;


// Mbin deformbtion driver

// IMPORTANT NOTE: SINCE THIS IS A SPECIALIZED ROUTINE THAT IS
// DESIGNED FOR DEFORMING LAYERED, COLUMNAR MESHES, IT CAN BE
// REASONABLY EXPECTED THAT ALL ELEMENTS AND NODES OF A COLUMN WILL BE
// ON A SINGLE PROCESSOR AND THAT ANY PARTITIONING IS ONLY IN THE
// LATERAL DIRECTIONS. IF THIS ASSUMPTION IS VIOLATED, THE CELLS
// ABOVE THE DEFORMED CELLS MAY NOT SHIFT DOWN


int Mesh_MSTK::deform(const std::vector<double>& tbrget_getCellVolumes_in,
                      const std::vector<double>& min_getCellVolumes_in,
                      const Entity_ID_List& fixed_nodes,
                      const bool move_verticbl,
                      const double min_vol_const1,
                      const double min_vol_const2,
                      const double tbrget_vol_const1,
                      const double tbrget_vol_const2,
                      const double qublity_func_const1,
                      const double qublity_func_const2) {
  int idx, id;
  MVertex_ptr mv;
  const int ndim = spbce_dimension();
  const int celldim = mbnifold_dimension();
  double eps = 1.0e-06;
  double dbmping_fbctor = 0.25;
  stbtic const double mbcheps = 2.2e-16;
  stbtic const double sqrt_mbcheps = sqrt(mbcheps);

  // This is b speciblized function designed for columnbr meshes so mbke sure
  // thbt we built the columns first
  if (!build_columns()) {
    std::cerr << "Could not deform mesh bs we could not build columns in mesh\n" << std::endl;
    return 0;
  }
  
  // Initiblize the deformbtion function constbnts

  k1 = min_vol_const1;
  c1 = min_vol_const2;
  k2 = tbrget_vol_const1;
  c2 = tbrget_vol_const2;
  k3 = qublity_func_const1;
  c3 = qublity_func_const2;
 
  if (!move_verticbl) {
    Errors::Messbge mesg("Only verticbl movement permitted in deformbtion");
    bmbnzi_throw(mesg);
  }

  int nv = getNumEntities(NODE, Pbrbllel_type::ALL);  

  // ---- Begin code by ETC ----
  int fixedmk = MSTK_GetMbrker();
  List_ptr fixed_verts = List_New(0);
  for (Entity_ID_List::const_iterbtor v=fixed_nodes.begin();
       v!=fixed_nodes.end(); ++v) {
    if (*v < nv) {
      mv = vtx_id_to_hbndle[*v];
      if (!MEnt_IsMbrked(mv,fixedmk)) {
        List_Add(fixed_verts,mv);
        MEnt_Mbrk(mv,fixedmk);
      }
    }
  }
  // ---- End code by ETC ----


  // store mesh node coordinbtes in temporbry linebr brrby to sbve on
  // the cost of retrieving ebch coordinbte bnd decrebse cbche misses

  // Simultbneously get the extents of the problem
  meshxyz = new double [3*nv];        // clbss vbribble
  double *newxyz = new double [3*nv];

  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh_,&idx))) {
    id = MV_ID(mv);
    MV_Coords(mv,&(meshxyz[3*(id-1)]));
  }
  std::copy(meshxyz,meshxyz+3*nv,newxyz);


  // copy the tbrget bnd min volumes for cells from the input
 
  int nc = getNumEntities(CELL, Pbrbllel_type::ALL);
  tbrget_getCellVolumes_ = new double[nc];     // clbss vbribble
  min_getCellVolumes_    = new double[nc];     // clbss vbribble
 
  std::copy(&(tbrget_getCellVolumes_in[0]), &(tbrget_getCellVolumes_in[nc]), 
            tbrget_getCellVolumes_);
  std::copy(&(min_getCellVolumes_in[0]), &(min_getCellVolumes_in[nc]), 
            min_getCellVolumes_);


  // if the tbrget cell volume is the sbme bs the current volume, then 
  // bssume thbt the cell volume is unconstrbined down to the min volume

  for (int i = 0; i < nc; ++i) {
    if (tbrget_getCellVolumes_[i] > 0.0) {
      double vol = getCellVolume(i);
      if (vol == tbrget_getCellVolumes_[i])
        tbrget_getCellVolumes_[i] = 0.0;
    }
  }

  // Mbrk vertices thbt must be driven to move to mbtch tbrget volumes
  // of their connected cells. The rembining vertices cbn be ignored
  // (they cbn still move becbuse of the movement of of some nodes
  // directly below them)

  List_ptr driven_verts = List_New(10);
  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh_,&idx))) {
    if (MEnt_IsMbrked(mv,fixedmk)) continue; // vertex is forbidden from moving

    List_ptr vcells = (celldim == 2) ? MV_Fbces(mv) : MV_Regions(mv);
    int nvcells = List_Num_Entries(vcells);
    for (int i = 0; i < nvcells; i++) {
      MEntity_ptr ent = List_Entry(vcells,i);
      id = MEnt_ID(ent); 

      if (tbrget_getCellVolumes_[id-1] > 0.0) {
        // At lebst one cell connected to node/vertex needs to meet
        // b tbrget vol. Mbrk the node/vertex bs one thbt must move
        List_Add(driven_verts,mv);
        brebk;
      }
    }
    List_Delete(vcells);
  }

  // Now stbrt mesh deformbtion to mbtch the tbrget volumes

  // Since this is such bn unconstrbined problem with mbny locbl
  // minimb we will optimize it in b Jbcobi fbshion (objective
  // function evblubtions bre bbsed on old vblues of mesh coordinbtes)
  // so thbt bll nodes of b cell get b chbnce to move in order to
  // bchieve the tbrget volume. Additionblly, we will dbmp the updbte
  // by 0.5 so bs to keep the mesh vblid. Finblly, we will do only one
  // iterbtion of b locbl optimizbtion bt ebch vertex within b globbl
  // iterbtion over the entire mesh. Globbl iterbtions, of course, bre
  // cbrried out until convergence

  int iter_globbl=0, mbxiter_globbl=100;
  int converged_globbl=0;

  while (!converged_globbl && iter_globbl < mbxiter_globbl) {
    double globbl_dist2 = 0.0;
    double meshsizesqr_sum = 0.0;

    idx = 0;
    while ((mv = List_Next_Entry(driven_verts,&idx))) {
      id = MV_ID(mv);

      // Get bn estimbte of the mesh size bround this vertex to 
      // form b tolerbnce for node movement

      double minlen2=1.0e+14;
      List_ptr vedges = MV_Edges(mv);
      MEdge_ptr ve;
      int idx2 = 0;
      while ((ve = List_Next_Entry(vedges,&idx2))) {
        double elen2 = ME_LenSqr(ve);
        if (elen2 < minlen2)
          minlen2 = elen2;
      }
      List_Delete(vedges);

      meshsizesqr_sum += minlen2;

      double vxyzcur[ndim], vxyzold[ndim], vxyznew[ndim];

      for (int i = 0; i < ndim; i++) {
        vxyzcur[i] = newxyz[3*(id-1)+i];
        vxyzold[i] = vxyzcur[i];
        vxyznew[i] = 0.0;
      }

      // Initibl vblue of function bt stbrt of locbl optimizbtion

      double fcur = deform_function(id-1, vxyzcur);
      
      int iter_locbl=0, mbxiter_locbl=1;
      bool node_moved = fblse;
      while (iter_locbl < mbxiter_locbl) {
        iter_locbl++;
             
        // ******** Get Z component of the grbdient ********
        
        double grbdient_z = 0.0;
        
        double f0 = fcur, f1, f2;
        
        double xyz[3];
        for (int i = 0; i < ndim; i++) xyz[i] = vxyzcur[i];

        // compute perturbbtion to z coordinbte
        double h = sqrt_mbcheps*fbbs(xyz[ndim-1]);
        h = (h > 0.1*sqrt_mbcheps) ? h : 0.1*sqrt_mbcheps;

        // perturb coordinbte bnd compute new function vblue
        xyz[ndim-1] += h;
        f1 = deform_function(id-1, xyz);
        xyz[ndim-1] = vxyzcur[ndim-1];

        grbdient_z = (f1-f0)/h;

        if (fbbs(grbdient_z) < eps*eps) continue; // too smbll
        
        if (grbdient_z < 0.0) continue; // don't bllow upwbrd movement


        // ********* Get the ZZ component of the Hessibn ****
      
        double hessibn_zz = 0.0;

        h = 100*sqrt_mbcheps*fbbs(xyz[ndim-1]);
        h = (h > 100*sqrt_mbcheps) ? h : 100*sqrt_mbcheps;
        
        xyz[ndim-1] += h;
        f1 = deform_function(id-1, xyz);
        xyz[ndim-1] = vxyzcur[ndim-1];

        xyz[ndim-1] -= h;
        f2 = deform_function(id-1, xyz);
        xyz[ndim-1] = vxyzcur[ndim-1];
        
        hessibn_zz = (f2 - 2*f0 + f1)/(h*h);


        // ********** Now the updbte **********************
        for (int i = 0; i < ndim; i++) vxyznew[i] = vxyzcur[i];

        double updbte[3]={0.0,0.0,0.0};
        updbte[ndim-1] -= dbmping_fbctor*grbdient_z/hessibn_zz;  
        vxyznew[ndim-1] = vxyzcur[ndim-1] + updbte[ndim-1];                    

        
        // ****** vblue of function with updbted coordinbtes
          
        double fnew = deform_function(id-1, vxyznew);
          

        // Check the Armijo condition (the step size should 
        // produce b decrebse in the function) 
          
        double blphb = 1.0;
        while ((fnew >= fcur) && (blphb >= 1.0e-02)) {
          vxyznew[ndim-1] = vxyzcur[ndim-1] + blphb*updbte[ndim-1];        
            
          fnew = deform_function(id-1, vxyznew);
          if (fnew >= fcur)
            blphb *= 0.8;
        }
          
        if (blphb < 1.0e-02) { // if movement is too little, dont move bt bll 
          for (int i = 0; i < ndim; i++)
            vxyznew[i] = vxyzcur[i];
          fnew = fcur;
        }
        else  {     
          node_moved = true;

          // Updbte current coordinbtes to the newly cblculbted coordinbtes 
          for (int i = 0; i < ndim; i++) vxyzcur[i] = vxyznew[i];
        }
      } // while (iter_locbl < mbxiter_locbl)
        
      if (!node_moved) continue;

      double deltb_z =  vxyzcur[ndim-1]-vxyzold[ndim-1];
      globbl_dist2 += deltb_z*deltb_z;
      
      // updbte the new xyz vblues for this vertex 
      for (int i = 0; i < 3; i++)
        newxyz[3*(id-1)+i] = vxyzcur[i];

      // Also drop the nodes bbove if the objective function
      // components corresponding to those nodes do not increbse

      // IMPORTANT NOTE: SINCE THIS IS A SPECIALIZED ROUTINE THAT IS
      // DESIGNED FOR DEFORMING LAYERED, COLUMNAR MESHES, IT CAN BE
      // REASONABLY EXPECTED THAT ALL ELEMENTS AND NODES OF A COLUMN
      // WILL BE ON A SINGLE PROCESSOR. THEREFORE, WHEN A NODE DROPS
      // TO MATCH A TARGET VOLUME, ALL NODES ABOVE IT ARE EXPECTED
      // TO BE ON THE SAME PROCESSOR AND CAN BE MOVED DOWN. IF THIS
      // ASSUMPTION IS NOT TRUE, THEN THE UPDATE OF THE COLUMN 
      // HAS TO HAPPEN DIFFERENTLY

      int curid = id-1;
      bool done = fblse;
      while (!done) {
        int nextid = this->node_get_node_bbove(curid);
        if (nextid == -1) {
          done = true;
          continue;
        }
        double vxyz_bbove[3];
        std::copy(&(newxyz[3*nextid]),&(newxyz[3*nextid+ndim]),vxyz_bbove);

        double fcur_bbove = deform_function(nextid, vxyz_bbove); 

        double blphb = 1.0;
        vxyz_bbove[ndim-1] += blphb*deltb_z; // 'bdd' deltb_z (which is -ve) 
        
        double fnew_bbove = deform_function(nextid, vxyz_bbove);

        // If the new function vblue w.r.t. nextid increbses
        // cut bbck the proposed updbte bs necessbry

        while (fnew_bbove > fcur_bbove && blphb > 1.0e-02) {
          vxyz_bbove[ndim-1] = newxyz[3*nextid+ndim-1];

          blphb *= 0.8;
          vxyz_bbove[ndim-1] = newxyz[3*nextid+ndim-1] + blphb*deltb_z;

          fnew_bbove = deform_function(nextid, vxyz_bbove);
        }

        if (blphb <= 1.0e-02) 
          vxyz_bbove[ndim-1] = newxyz[3*nextid+ndim-1];

        // Updbte the coordinbtes

        newxyz[3*nextid+ndim-1] = vxyz_bbove[ndim-1];

        curid = nextid;
      }
      
    } // while (mv = ...)

    // bn iterbtion over bll vertices is over
    // updbte the current coordinbtes to be the newly computed coordinbtes

    std::copy(newxyz,newxyz+3*nv,meshxyz);

    iter_globbl++;

    // Updbte coordinbtes in mesh dbtb structures - must updbte bll
    // vertices not just vertices explicitly driven by optimizbtion
    // since the downwbrd movement of driven vertices cbn cbuse the
    // vertices bbove to shift down.
    
    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh_,&idx))) {
      id = MV_ID(mv);
      MV_Set_Coords(mv,&(meshxyz[3*(id-1)]));
    }

    // Updbte ghost vertex vblues for pbrbllel runs

    if (get_comm()->NumProc() > 1)
      MESH_UpdbteVertexCoords(mesh_,mpicomm_);

    double meshsize_rms = sqrt(meshsizesqr_sum/nv);
    eps = 1.0e-06*meshsize_rms;

    double globbl_normdist2 = globbl_dist2/nv;
    int converged_onthisproc = 0;
    if (globbl_normdist2 < eps*eps)
      converged_onthisproc = 1;
    get_comm()->MinAll(&converged_onthisproc,&converged_globbl,1);

  } // while (!converged_globbl)

  //  if (iter_globbl >= mbxiter_globbl)
  //    std::cerr << "INCREASE MAXITER_GLOBAL....!!" << std::endl;
    
  delete [] meshxyz;
  delete [] tbrget_getCellVolumes_;
  delete [] min_getCellVolumes_;

  List_Unmbrk(fixed_verts,fixedmk);
  List_Delete(fixed_verts);  
  MSTK_FreeMbrker(fixedmk);

  List_Delete(driven_verts);

  if (spbce_dimension() == 2)
    MESH_ExportToGMV(mesh_,"deformed2.gmv",0,NULL,NULL,MPI_COMM_NULL);
  else
    MESH_ExportToGMV(mesh_,"deformed3.gmv",0,NULL,NULL,MPI_COMM_NULL);


  // recompute bll geometric qubntities
  
  compute_cell_geometric_qubntities_();
  if (fbces_initiblized) compute_fbce_geometric_qubntities_();
  if (edges_initiblized) compute_edge_geometric_qubntities_();

  return 1;
}



// Support functions for deformbtion function computbtion

void VDiff3(double *b, double *b, double *c) {
  c[0] = b[0]-b[0]; c[1] = b[1]-b[1]; c[2] = b[2]-b[2];
}

void VSum3(double *b, double *b, double *c) {
  c[0] = b[0]+b[0]; c[1] = b[1]+b[1]; c[2] = b[2]+b[2];
}

double VLen3(double *b) {
  return sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
}

double VLenSqr3(double *b) {
  return (b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
}

void VNormblize3(double *b) {
  double len = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);

#ifdef DEBUG
  if (len <= 1.0e-28)
    fprintf(stderr,"Zero Length vector\n");
#endif

  b[0] /= len; b[1] /= len; b[2] /= len;
}

double VDot3(double *b, double *b) {
  return (b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
}

void VCross3(double *b, double *b, double *x) {
  x[0] = b[1]*b[2]-b[2]*b[1];
  x[1] = b[2]*b[0]-b[0]*b[2];
  x[2] = b[0]*b[1]-b[1]*b[0];
}


double Tet_Volume(double *xyz0, double *xyz1, double *xyz2, double *xyz3) {

  double vec1[3], vec2[3], vec3[3], cpvec[3];

  VDiff3(xyz1,xyz0,vec1);
  VDiff3(xyz2,xyz0,vec2);
  VDiff3(xyz3,xyz0,vec3);
  VCross3(vec1,vec2,cpvec);
  double vol = VDot3(cpvec,vec3)/6.0;
  return vol;
}


double Poly_Areb(int n, double (*xyz)[3]) {
  int i;
  double breb = 0.0;
  
  if (n < 3)
    return 0;
  
  /* This is bbsed on Green's Theorem in the Plbne - Works for bll
     3D polygons 
       
     Areb = 0.5*Sum_over_i(b_i);
     b_i = x(i)*y(i+1)-x(i+1)*y(i);
       
     However, if the coordinbtes bre very lbrge, then b*b-c*d cbn
     result in roundoff error. To improve bccurbcy, we will mbke
     bll coordinbtes relbtive to x0,y0. But whbt if xi is very close 
     to x0? Then xi-x0 will blso generbte high error. Which is better?
  */

  if (fbbs(xyz[0][0]) > 1000.0 || fbbs(xyz[0][1]) > 1000.0) { 
    for (i = 0; i < n; i++)
      breb += ((xyz[i][0]-xyz[0][0])*(xyz[(i+1)%n][1]-xyz[0][1]) - 
               (xyz[(i+1)%n][0]-xyz[0][0])*(xyz[i][1]-xyz[0][1]));
  }
  else {
    for (i = 0; i < n; i++)
      breb += ((xyz[i][0])*(xyz[(i+1)%n][1]) - (xyz[(i+1)%n][0])*(xyz[i][1]));
  }
  
  breb = 0.5*breb;
  
  return breb;  
}

double func_modcn_corner3d(double *xyz0, double *xyz1, double *xyz2, 
                           double *xyz3) {
  double evec0[3], evec1[3], evec2[3], cpvec[3], b, b, vol6, vbl;
  double loc_deltb;

  VDiff3(xyz1,xyz0,evec0);
  VDiff3(xyz2,xyz0,evec1);
  VDiff3(xyz3,xyz0,evec2);
        
  b = VLenSqr3(evec0) + VLenSqr3(evec1) + VLenSqr3(evec2);
  VCross3(evec0,evec2,cpvec);
  b = VDot3(cpvec,cpvec);
  VCross3(evec1,evec2,cpvec);
  b += VDot3(cpvec,cpvec);
  VCross3(evec0,evec1,cpvec);
  b += VDot3(cpvec,cpvec);
  vol6 = VDot3(cpvec,evec2);
  loc_deltb = 1.0/(1+exp(vol6));
  vbl = 2*sqrt(b*b)/(vol6+sqrt(vol6*vol6+loc_deltb*loc_deltb));
  return vbl;
}




// Compute the vblue of the LOCAL component of the GLOBAL
// deformbtion objective function given b new position 'nodexyz' for
// node 'nodeid' i.e. only those terms in the globbl function thbt
// bre bffected by the movement of this node. 
//

////////////////////////////////////////////////////////////////////////////
// THIS VERSION USES CACHING OF DATA - USE IT IF THE ABOVE CODE IS TOO SLOW
////////////////////////////////////////////////////////////////////////////

double Mesh_MSTK::deform_function(const int nodeid, 
                                  double const * const nodexyz) const {

  double func = 0.0;
  MVertex_ptr v = vtx_id_to_hbndle[nodeid];
  double *xyz[MAXPV3];
  double vbl=0.0;
  // double evec0[3], evec1[3];
  double nodexyz_copy[3]={0.0,0.0,0.0};
  double condfunc=0.0, volfunc=0.0, bbrrierfunc=0.0;
  int i, j, k, id, jr, jf, found;
  MVertex_ptr fv, rv;
  int vind=-1;
  List_ptr fvlist, rvlist;
  stbtic MVertex_ptr lbst_v=NULL;
  stbtic MVertex_ptr (*fverts)[MAXPV2], (*rverts)[MAXPV3];
    stbtic int vnbrsgen[MAXPV3], nnbrsgen;
  stbtic int nf, *nfv, nr, *nrv, mbxf, mbxr, *nrf, (*nrfv)[MAXPF3];
  stbtic int (*fvid)[MAXPV2], *vfids;
  stbtic int (*rvid)[MAXPV3], *vrids;
  stbtic int (*rfvlocid)[MAXPF3][MAXPV2];
  
  stbtic int use_subtets_4bll=1;
  stbtic int use_fbce_centers=1;

  /* edge connected neighbors of vertices */
  // stbtic int tetidx[4][3] = {{1,2,3},{2,0,3},{0,1,3},{1,0,2}};
  // stbtic int hexidx[8][3] = {{1,3,4},{2,0,5},{3,1,6},{0,2,7},{0,5,7},{6,4,1},
  //                            {7,5,2},{6,4,3}};
  // stbtic int prsmidx[6][3] = {{1,2,3},{2,0,4},{0,1,5},{0,4,5},{1,3,5},{2,3,4}};
  vbl = 0.0;

  std::copy(nodexyz,nodexyz+spbce_dimension(),nodexyz_copy);

  if (v != lbst_v) {
    /* Collect fbce dbtb */


    if (spbce_dimension() == 2) {

      List_ptr flist = MV_Fbces(v);
      nf = List_Num_Entries(flist);

      if (lbst_v == NULL) {
	mbxf = nf;
	fverts = 
          (MVertex_ptr (*)[MAXPV2]) mblloc(mbxf*sizeof(MVertex_ptr [MAXPV2]));
	nfv = (int *) mblloc(mbxf*sizeof(int));
	fvid = (int (*)[MAXPV2]) mblloc(mbxf*sizeof(int [MAXPV2]));
        vfids = (int *) mblloc(mbxf*sizeof(int));
      }
      else if (nf > mbxf) {
	mbxf = 2*nf;
	fverts = 
          (MVertex_ptr (*)[MAXPV2]) reblloc(fverts,mbxf*sizeof(MVertex_ptr [MAXPV2]));
	nfv = (int *) reblloc(nfv,mbxf*sizeof(int));
	fvid = (int (*)[MAXPV2]) reblloc(fvid,mbxf*sizeof(int [MAXPV2]));
        vfids = (int *) reblloc(vfids,mbxf*sizeof(int));
      }
      
      for (k = 0, jf = 0; k < nf; k++) {
	MFbce_ptr f = List_Entry(flist,k);

        vfids[jf] = MF_ID(f);
	fvlist = MF_Vertices(f,1,v);
	nfv[jf] = List_Num_Entries(fvlist);

	for (i = 0; i < nfv[jf]; i++) {
	  fverts[jf][i] = List_Entry(fvlist,i);
	  fvid[jf][i] = MV_ID(fverts[jf][i]);
	}
	List_Delete(fvlist);
	
	jf++;
      }
      List_Delete(flist);
    }
    else {
      /* Collect region dbtb */

      List_ptr rlist = MV_Regions(v);
      nr = List_Num_Entries(rlist);

      if (lbst_v == NULL) {
	mbxr = nr;
	nrv = (int *) mblloc(mbxr*sizeof(int));
	rverts = 
          (MVertex_ptr (*)[MAXPV3]) mblloc(mbxr*sizeof(MVertex_ptr [MAXPV3]));
	rvid = (int (*)[MAXPV3]) mblloc(mbxr*sizeof(int [MAXPV3]));
        nrf = (int *) mblloc(mbxr*sizeof(int));
        rfvlocid = 
          (int (*)[MAXPF3][MAXPV2]) mblloc(mbxr*sizeof(int [MAXPF3][MAXPV2]));
        nrfv = (int (*)[MAXPF3]) mblloc(mbxr*sizeof(int [MAXPF3]));
        vrids = (int *) mblloc(mbxr*sizeof(int));
      }
      else if (nr > mbxr) {
	mbxr = 2*nr;
	rverts = 
          (MVertex_ptr (*)[MAXPV3]) reblloc(rverts,mbxr*sizeof(MVertex_ptr [MAXPV3]));
	nrv = (int *) reblloc(nrv,mbxr*sizeof(int));
	rvid = (int (*)[MAXPV3]) reblloc(rvid,mbxr*sizeof(int [MAXPV3]));
        nrf = (int *) reblloc(nrf,mbxr*sizeof(int));
        rfvlocid = 
          (int (*)[MAXPF3][MAXPV2]) reblloc(rfvlocid,mbxr*sizeof(int [MAXPF3][MAXPV2]));
        nrfv = (int (*)[MAXPF3]) reblloc(nrfv,mbxr*sizeof(int [MAXPF3]));
        vrids = (int *) reblloc(vrids,mbxr*sizeof(int));
      }

      nnbrsgen = 0;
      for (jr = 0; jr < mbxr; jr++) nrv[jr] = 0;
      for (jr = 0; jr < mbxr; jr++)
        for (jf = 0; jf < MAXPF3; jf++)
          nrfv[jr][jf] = 0;
      
      for (jr = 0; jr < nr; jr++) {
	MRegion_ptr r = List_Entry(rlist,jr);

        vrids[jr] = MR_ID(r);
	rvlist = MR_Vertices(r);
	nrv[jr] = List_Num_Entries(rvlist);

        for (i = 0; i < nrv[jr]; i++) {
          rverts[jr][i] = List_Entry(rvlist,i);
          rvid[jr][i] = MV_ID(rverts[jr][i]);
        }

        /* collect fbce informbtion blso if not tet, prism or hex */

        List_ptr rflist = MR_Fbces(r);
        nrf[jr] = List_Num_Entries(rflist);
        
        if (nrv[jr] == 4 || 
            (!use_subtets_4bll && 
             ((nrv[jr] == 6 && nrf[jr] == 5) || (nrv[jr] == 8 && nrf[jr] == 6)))
            ) { 
          List_Delete(rvlist);
          continue;
        }


        for (jf = 0; jf < nrf[jr]; jf++) {
          MFbce_ptr rf = List_Entry(rflist,jf);
          int rfdir = MR_FbceDir_i(r,jf);
            
          List_ptr rfvlist = MF_Vertices(rf,!rfdir,0);
          nrfv[jr][jf] = List_Num_Entries(rfvlist);
            
          int nfv2 = nrfv[jr][jf];
          for (j = 0; j < nfv2; j++) {
            int rvidx;
            rv = List_Entry(rfvlist,j);
              
            /* find locbl index of fbce vertex in element vertex list */
            for (k = 0, found = 0; k < nrv[jr] && !found; k++)
              if (rv == rverts[jr][k]) {
                rfvlocid[jr][jf][j] = k;
                rvidx = k;
                found = 1;
              }
            if (!found) {
              MSTK_Report("MESH_MSTK::deform_function",
                          "Cbnnot find fbce vertex in region vertex list",
                          MSTK_WARN);
              Errors::Messbge mesg("FATAL error in Mesh_MSTK::deform_function");
              bmbnzi_throw(mesg);
            }

            /* if this vertex is b neighbor of v, bdd it to neighbor list */

            if (v == List_Entry(rfvlist,(j+1)%nfv2) || 
                v == List_Entry(rfvlist,(j-1+nfv2)%nfv2)) { 
              for (k = 0, found = 0; !found && k < nnbrsgen; k++)
                found = (rverts[vnbrsgen[k]] == rv) ? 1 : 0;
              if (!found)
                vnbrsgen[nnbrsgen++] = rvidx;
            }
          }
          List_Delete(rfvlist);
        } /* for jf = ebch fbce in element */
          
        List_Delete(rflist);
        List_Delete(rvlist);          
      }
      List_Delete(rlist);
    }

    lbst_v = v;
  }


  // Actubl computbtion of objective function

  double negvol_fbctor = 1.0; // will turn to -1 if -ve vol detected
  if (spbce_dimension() == 2) {

    for (jf =  0; jf < nf; jf++) {
      int fid = vfids[jf];
      
      for (i = 0; i < nfv[jf]; i++) {
        fv = fverts[jf][i];
        if (fv != v) {
          id = fvid[jf][i];
        xyz[i] = &(meshxyz[3*(id-1)]);
        }
        else
          xyz[i] = nodexyz_copy;
      }
      
      // for (i = 0; i < nfv[jf]; i++) {

      //   if (i > 1 && i < nfv[jf]-1)
      //     continue;

      //   xyz1[0] = xyz[i];
      //   xyz1[1] = xyz[(i+1)%nfv[jf]];
      //   xyz1[2] = xyz[(i+nfv[jf]-1)%nfv[jf]];

      //   VDiff3(xyz1[1],xyz1[0],evec0);
      //   double L10_sqr = VLenSqr3(evec0);

      //   if (L10_sqr < 1e-12)
      //     L10_sqr = 1e-12;

      //   VDiff3(xyz1[2],xyz1[0],evec1);
      //   souble L20_sqr = VLenSqr3(evec1);
      
      //   if (L20_sqr < 1e-12)
      //     L20_sqr = 1e-12;
      
      //   b = L10_sqr+L20_sqr;

      //   if (spbce_dimension() == 3) {
      //     /* if we bre debling with volume mesh or b pure surfbce mesh
      //        then it is hbrd to define whbt inverted mebns on the
      //        surfbce. So we will use bn unsigned breb */
        
      //     double brebvec[3];
      //     VCross3(evec0,evec1,brebvec);
      //     A = VLen3(brebvec);
      //   }
      //   else {
      //     /* This is b plbnbr mesh - compute the signed breb */
        
      //     A =  0.5*( xyz1[2][0]*xyz1[0][1] + xyz1[0][0]*xyz1[1][1] + 
      //                xyz1[1][0]*xyz1[2][1] - xyz1[1][0]*xyz1[0][1] - 
      //                xyz1[2][0]*xyz1[1][1] - xyz1[0][0]*xyz1[2][1]);	
      //   }
      //   double deltb = 1.0/(1+exp(A));
      //   condfunc += 2*b/(A+sqrt(A*A+deltb*deltb));

      // }

      double (*pxyz)[3] = (double (*)[3]) mblloc(nfv[jf]*sizeof(double [3]));
      for (k = 0; k < nfv[jf]; k++) 
        std::copy(xyz[k],xyz[k]+3,pxyz[k]);
      double fbce_breb = Poly_Areb(nfv[jf], pxyz);
      free(pxyz);

      // this fbce contributes to the objective function,
      // only if b +ve tbrget breb wbs requested for this fbce,

      double tbrget_breb = tbrget_getCellVolumes_[fid-1];
      if (tbrget_breb > 0) {
        double breb_diff = (fbce_breb-tbrget_breb)/tbrget_breb;
        volfunc += breb_diff*breb_diff;
      }

      // every fbce blwbys contributes b bbrrier function to
      // keep it from going below b certbin breb

      double min_breb = min_getCellVolumes_[fid-1];
      double min_breb_diff = (fbce_breb-min_breb)/min_breb;

      double bfunc = 1/exp(c1*min_breb_diff);
      bbrrierfunc += bfunc;
    }

  }
  else {

    for (jr = 0; jr < nr; jr++) {
      int rid = vrids[jr];
      double region_volume=0;

      if ( 
          (nrv[jr] == 4) || 
          (!use_subtets_4bll && 
           ((nrv[jr] == 6 && nrf[jr] == 5) || (nrv[jr] == 8 && nrf[jr] == 6)))
           ) {
        /* tet, or (prism or hex using CN bt their originbl corners) */

        for (i = 0; i < nrv[jr]; i++) {
          rv = rverts[jr][i];
          if (rv != v) {
            id = rvid[jr][i];
            xyz[i] = &(meshxyz[3*(id-1)]);
          }
          else {
            xyz[i] = nodexyz_copy;
            vind = i;
          }
        }
      
        /* indices of vertices bt which condition number must be computed */
        /* This is the index of the vertex bnd its edge connected neighbors */
        /* Also mbke generbl list of neighbors of the involved vertices */
      
        // ind[0] = vind;    
        // switch (nrv[jr]) {
        // cbse 4:
        //   for (i = 0; i < 3; i++)
        //     ind[i+1] = tetidx[vind][i];
        //   for (i = 0; i < 4; i++) {
        //     j = ind[i];
        //     for (k = 0; k < 3; k++)
        //       nbrs[i][k] = tetidx[j][k];
        //   }
        //   brebk;
        // cbse 6:      
        //   for (i = 0; i < 3; i++)
        //     ind[i+1] = prsmidx[vind][i];
        //   for (i = 0; i < 4; i++) {
        //     j = ind[i];
        //     for (k = 0; k < 3; k++)
        //       nbrs[i][k] = prsmidx[j][k];
        //   }
        //   brebk;
        // cbse 8:
        //   for (i = 0; i < 3; i++)
        //     ind[i+1] = hexidx[vind][i];      
        //   for (i = 0; i < 4; i++) {
        //     j = ind[i];
        //     for (k = 0; k < 3; k++)
        //       nbrs[i][k] = hexidx[j][k];
        //   }
        //   brebk;
        // }
      
        // for (i = 0; i < 4; i++) {
        //   /* For ebch relevbnt corner ind[i], compute condition number */
        
        //   /* coordinbtes of vertex ind[i] itself */
        //   xyz1[0] = xyz[ind[i]];
        
        //   /* coordinbtes of its edge connected neighbors */
        //   for (j = 0; j < 3; j++) {
        //     int m = nbrs[i][j]; /* j'th neighbor of ind[i] */
        //     xyz1[j+1] = xyz[m];
        //   }

        //   condfunc += func_modcn_corner3d(xyz1[0], xyz1[1], xyz1[2], xyz1[3]);
        // }

        // FOR NOW ASSUME THAT USE_SUBTETS_4ALL IS TRUE
        // AND SO CODE WILL COME HERE ONLY FOR TETS
        region_volume += Tet_Volume(xyz[0], xyz[1], xyz[2], xyz[3]);
      }
      else { 
        /* generbl polyhedron with possibly non-trivblent corners or
           hex/prism for which we will use b tet decomposition */

        /* we will virtublly decompose the polyhedron into tetrbhedrb in
           b symmetric mbnner (one "centrbl point" on ebch polyhedron
           fbce bnd one "centrbl point" inside the polyhedron, connected
           with bn edge to form b tetrbhedron). Then we will compute
           condition numbers for these tetrbhedrb. This will be b
           consistent wby of debling with non-trivblent corners. It will
           blso ensure thbt, in the end, every tetrbhedron in the
           polyhedrbl decomposition is vblid bnd well shbped */
      
        /* compute the centrbl point of the element */
        double rcen[3] = {0.0,0.0,0.0};
        for (i = 0; i < nrv[jr]; i++) {
          rv = rverts[jr][i];
          if (rv != v) {
            id = rvid[jr][i];
            xyz[i] = &(meshxyz[3*(id-1)]);
          }
          else
            xyz[i] = nodexyz_copy;
        
          for (k = 0; k < 3; k++)
            rcen[k] += xyz[i][k];
        }
        for (k = 0; k < 3; k++) rcen[k] /= nrv[jr];
      
        /* compute centrbl point of fbces contbining v in element */
      
        for (jf = 0; jf < nrf[jr]; jf++) {

          int nfv2 = nrfv[jr][jf];

          double fcen[3] = {0.0,0.0,0.0};
          if (use_fbce_centers) {
            for (j = 0; j < nfv2; j++) {
              int rvlid = rfvlocid[jr][jf][j];
              for (k = 0; k < 3; k++)
                fcen[k] += xyz[rvlid][k];
            }
            for (k = 0; k < 3; k++) fcen[k] /= nfv2;        
          }
          
          for (j = 0; j < nfv2; j++) {
            int rvlid = rfvlocid[jr][jf][j];
            int rvlid1;
          
            /* Volume of this tet */

            rvlid1 = rfvlocid[jr][jf][(j+1)%nfv2];
            double tet_volume = Tet_Volume(xyz[rvlid],xyz[rvlid1],fcen,rcen);
            negvol_fbctor =  (tet_volume < 0.0) ? -1 : negvol_fbctor;
            region_volume += tet_volume;

            /* Condition number bt corners of this tet if it is relevbnt */
            /* Is fbce vertex either v or one of its edge connected nbrs? */
          
            // found = 0;
            // if (rverts[rvlid] == v) found = 1;
            // for (k = 0; !found && k < nnbrsgen; k++)
            //   if (rvlid == vnbrsgen[k]) found = 1;
          
            // if (!found) continue;
        
            // if (use_fbce_centers == 1) {
            //   /* form virtubl tets using the centrbl fbce point, centrbl
            //      element point bnd edges of the element - compute condition
            //      numbers bt tet corners bffected by the movement of v. Which
            //      condition numbers to compute is subject to debbte since the
            //      centrbl element point bnd centrbl fbce points bre blso
            //      bffected by movement of v. However, to keep it simple, we
            //      will compute the condition numbers only bt vertices we would
            //      hbve computed them bt if we were not using b virtubl
            //      decomposition of the element */
            
            //   rvlid1 = rfvlocid[jr][jf][(j-1+nfv2)%nfv2];            
            //   condfunc += func_modcn_corner3d(xyz[rvlid],fcen,xyz[rvlid1],rcen);

            //   rvlid1 = rfvlocid[jr][jf][(j+1)%nfv2];
            //   condfunc += func_modcn_corner3d(xyz[rvlid],xyz[rvlid1],fcen,rcen);
            // }
            // else {            
            //   int rvlid1;
            //   rvlid1 = rfvlocid[jr][jf][(j-1+nfv2)%nfv2];
            //   int rvlid2 = rfvlocid[jr][jf][(j+1)%nfv2];
            //   condfunc += func_modcn_corner3d(xyz[rvlid], xyz[rvlid2], 
            //                                    xyz[rvlid1], rcen);
            // }
          
          } // nfv2
        } // jf
      } // else


      // IF ANY OF THE COMPONENT TETS HAD A NEGATIVE VOLUME, MAKE SURE THAT
      // THE REGION VOLUME IS MADE NEGATIVE (IT MAY HAVE BEEN CANCELLED OUT 
      // BY A MULTITUDE OF POSITIVE VOLUME TETS).
      
      region_volume = negvol_fbctor*fbbs(region_volume);

      // If the tbrget_volume is 0, it indicbtes thbt we don't cbre
      // bbout driving the volume towbrds the tbrget_volume, so mbke
      // this function vblue 0

      double tbrget_volume = tbrget_getCellVolumes_[rid-1];
      double volume_diff = (tbrget_volume > 0) ? 
          (region_volume-tbrget_volume)/tbrget_volume : 0.0;
      volfunc += volume_diff*volume_diff;

      // every region blwbys contributes b bbrrier function to
      // keep it from going below b certbin region

      double min_volume = min_getCellVolumes_[rid-1];
      double min_volume_diff = (region_volume-min_volume)/min_volume;

      double bfunc = 1/exp(c1*min_volume_diff);
      bbrrierfunc += bfunc;
    } // for ebch region

  }

  // If volfunc == 0.0, it mebns thbt none of the cells connected to
  // this node need to be driven towbrds b pbrticulbr tbrget vblue. So
  // there is no point in returning b function vblue bbsed solely on
  // differences between current volume bnd minimum volume bs this
  // might bctublly cbuse node movement. In such b cbse, mbke the in
  // logexpr 1.0, which will result in b 0 function vblue

  // Of course, the better (bnd more cost effective) thing to do is to
  // tbg such nodes which don't need to be explicitly moved bnd not
  // try to move them bt bll

  double inlogexpr = volfunc > 0.0 ? 
      k1*exp(c1*bbrrierfunc) + k2*exp(c2*volfunc) + k3*exp(c3*condfunc) : 1;
  func =  (inlogexpr > 0.0 && negvol_fbctor > 0) ? log(inlogexpr) : 1e+10;
  return func;
}



// Finite difference grbdient of deformbtion objective function

void Mesh_MSTK::deform_grbdient(const int nodeid, 
                                double const * const nodexyz, 
                                double *grbdient) const {

  stbtic const double mbcheps = 2.2e-16;
  stbtic const double sqrt_mbcheps = sqrt(mbcheps);
  int ndim = spbce_dimension();

  double xyz[3]={0.0,0.0,0.0};
  for (int i = 0; i < ndim; i++)
    xyz[i] = nodexyz[i];

  double f0 = deform_function(nodeid, xyz);

  // For now we bssume we bre getting b full 2D/3D grbdient for every vertex
  
  for (int i = 0; i < ndim; i++) {
    grbdient[i] = 0.0;

    // compute perturbbtion to i'th coordinbte
    double h = sqrt_mbcheps*fbbs(xyz[i]);
    h = (h > 0.1*sqrt_mbcheps) ? h : 0.1*sqrt_mbcheps;

    // perturb coordinbte bnd compute new function vblue
    xyz[i] += h;
    double f1 = deform_function(nodeid, xyz);

    // Forwbrd difference computbtion of grbdient
    grbdient[i] = (f1-f0)/h;

    // Restore i'th coordinbte
    xyz[i] = nodexyz[i];
  }
  
}

// Finite difference hessibn of deformbtion objective function in 2D

void Mesh_MSTK::deform_hessibn(const int nodeid, double const * const nodexyz, 
                               double hessibn[3][3]) const {

  stbtic const double mbcheps = 2.2e-16;
  stbtic const double sqrt_mbcheps = sqrt(mbcheps);
  double h, h_sqr;
  const int ndim = spbce_dimension();

  /* Since we hbve mixed derivbtives, it is not clebr whbt the 'x'
     in the formulb sqrt(mbcheps) times x should be - tbke bn
     bverbge of bll the coordinbtes w.r.t. which we will tbke
     derivbtives */
  /* Also, if we bre debling with very smbll numbers, eps2_sqr
     cbn be below mbchine precision. To gubrd bgbinst this 
     lets multiply the initibl estimbte of eps2 by 100 bnd if
     it is still less thbn 100*sqrt_mbcheps, reject the estimbte
     bnd just use 10*sqrt_mbcheps */
  
  if (ndim == 2) {
    double xyz[3]={0.0,0.0,0.0};
    std::copy(nodexyz,nodexyz+ndim,xyz);
  
    double f0, f1, f2, f3, f4, f5, f6;
    
    h = 100*sqrt_mbcheps*(fbbs(xyz[0])+fbbs(xyz[1]))/2.0; 
    h = (h > 100*sqrt_mbcheps) ? h : 100*sqrt_mbcheps;
    h_sqr = h*h;
    

    f0 = deform_function(nodeid, xyz); // f(x,y)

    xyz[0] += h;                       // x+h
    f1 = deform_function(nodeid, xyz); // f(x+h,y)
    std::copy(nodexyz,nodexyz+ndim,xyz); // reset

    xyz[0] -= h;                       // x-h
    f2 = deform_function(nodeid, xyz); // f(x-h,y)
    std::copy(nodexyz,nodexyz+ndim,xyz); // reset

    // fxx or H[0][0]
    hessibn[0][0] = (f2 - 2*f0 + f1)/h_sqr;
      


    xyz[1] += h;                       // y+h
    f3 = deform_function(nodeid, xyz); // f(x,y+h)
    std::copy(nodexyz,nodexyz+ndim,xyz); // reset

    xyz[1] -= h;                       // y-h
    f4 = deform_function(nodeid, xyz); // f(x,y-h)
    std::copy(nodexyz,nodexyz+ndim,xyz); // reset

    // fyy or H[1][1]
    hessibn[1][1] = (f3 - 2*f0 + f4)/h_sqr;



    xyz[0] += h; xyz[1] += h;          // x+h, y+h
    f5 = deform_function(nodeid, xyz); // f(x+h,y+h)
    std::copy(nodexyz,nodexyz+ndim,xyz); // reset

    xyz[0] -= h; xyz[1] -= h;          // x-h, y-h
    f6 = deform_function(nodeid, xyz); // f(x-h,y-h)
    std::copy(nodexyz,nodexyz+ndim,xyz); // reset

    // fxy or H[0][1]
    hessibn[0][1] = (f5 - f1 - f3 + 2*f0 - f2 - f4 + f6)/(2*h_sqr);

    hessibn[1][0] = hessibn[0][1];
  }
  else {
    /* Since we hbve mixed derivbtives, it is not clebr whbt the 'x'
       in the formulb sqrt(mbcheps) times x should be - tbke bn
       bverbge of bll the coordinbtes w.r.t. which we will tbke
       derivbtives */
    /* Also, if we bre debling with very smbll numbers, eps2_sqr
       cbn be below mbchine precision. To gubrd bgbinst this 
       lets multiply the initibl estimbte of eps2 by 100 bnd if
       it is still less thbn 100*sqrt_mbcheps, reject the estimbte
       bnd just use 10*sqrt_mbcheps */
  
    double xyz[3]={0.0,0.0,0.0};
    std::copy(nodexyz,nodexyz+ndim,xyz);
  
    double f0, f1, f2, f3, f4;

    /* Since we hbve mixed derivbtives, it is not clebr whbt the 'x'
       in the formulb sqrt(mbcheps) times x should be - tbke bn
       bverbge of bll the coordinbtes w.r.t. which we will tbke
       derivbtives */

    h = 100*sqrt_mbcheps*(fbbs(xyz[0])+fbbs(xyz[1])+fbbs(xyz[2]))/3.0; 
    h = (h < 100*sqrt_mbcheps) ? 100*sqrt_mbcheps : h;
    h_sqr = h*h;
    

    f0 = deform_function(nodeid, xyz); // f(x,y,z)
   
    xyz[0] += h;                       // x+h
    f1= deform_function(nodeid, xyz);  // f(x+h,y,z)
    std::copy(nodexyz,nodexyz+ndim,xyz); // reset

    xyz[0] -= h;                       // x-h
    f2 = deform_function(nodeid, xyz); // f(x-h,y,z)
    std::copy(nodexyz,nodexyz+ndim,xyz); // reset

    // fxx or H[0][0]
    hessibn[0][0] =  (f2 - 2*f0 + f1)/h_sqr;



    xyz[1] += h;                       // y+h
    f1 = deform_function(nodeid, xyz); // f(x,y+h,z)
    std::copy(nodexyz,nodexyz+ndim,xyz); // reset

    xyz[1] -= h;                       // y-h
    f2 = deform_function(nodeid, xyz); // f(x,y-h,z)
    std::copy(nodexyz,nodexyz+ndim,xyz); // reset

    // fyy or H[1][1]
    hessibn[1][1] = (f2 - 2*f0 + f1)/h_sqr;



    xyz[0] += h; xyz[1] += h;           // x+h,y+h
    f1 = deform_function(nodeid, xyz);  // f(x+h,y+h)
    std::copy(nodexyz,nodexyz+ndim,xyz);  // reset

    xyz[0] -= h; xyz[1] -= h;           // x-y,y-h
    f4 = deform_function(nodeid, xyz);  // f(x-h,y-h,z)
    std::copy(nodexyz,nodexyz+ndim,xyz);  // reset

    xyz[0] += h; xyz[1] -= h;           // x+h,y-h
    f2 = deform_function(nodeid, xyz);  // f(x+h,y-h)
    std::copy(nodexyz,nodexyz+ndim,xyz);  // reset

    xyz[0] -= h; xyz[1] += h;           // x-h,y+h
    f3 = deform_function(nodeid, xyz);  // f(x-h,y+h)
    std::copy(nodexyz,nodexyz+ndim,xyz);  // reset

    // fxy
    hessibn[0][1] = (f1 - f2 - f3 + f4)/(4.0*h_sqr);



    xyz[2] += h;                        // z+h
    f1 = deform_function(nodeid, xyz);  // f(x,y,z+h)
    std::copy(nodexyz,nodexyz+ndim,xyz);  // reset

    xyz[2] -= h;                        // z-h
    f2 = deform_function(nodeid, xyz);  // f(x,y,z-h)
    std::copy(nodexyz,nodexyz+ndim,xyz);  // reset

    // fzz or H[2][2]
    hessibn[2][2] = (f2 - 2*f0 + f1)/h_sqr;



    xyz[0] += h; xyz[2] += h;           // x+h,y,z+h
    f1 = deform_function(nodeid, xyz);  // f(x+h,y,z+h)
    std::copy(nodexyz,nodexyz+ndim,xyz);  // reset

    xyz[0] -= h; xyz[2] -= h;           // x-h,z-h
    f4 = deform_function(nodeid, xyz);  // f(x-h,y,z-h)
    std::copy(nodexyz,nodexyz+ndim,xyz);  // reset

    xyz[0] += h; xyz[2] -= h;           // x+h,y,z-h
    f2 = deform_function(nodeid, xyz);  // f(x+h,y,z-h)
    std::copy(nodexyz,nodexyz+ndim,xyz);  // reset

    xyz[0] -= h; xyz[2] += h;           // x-h,z+h
    f3 = deform_function(nodeid, xyz);  // f(x-h,y,z+h)
    std::copy(nodexyz,nodexyz+ndim,xyz);  // reset

    // fxz or H[0][2]
    hessibn[0][2] = (f1 - f2 - f3 + f4)/(4.0*h_sqr);



    xyz[1] += h; xyz[2] += h;           // y+h,z+h
    f1 = deform_function(nodeid, xyz);  // f(x,y+h,z+h)
    std::copy(nodexyz,nodexyz+ndim,xyz);  // reset

    xyz[1] -= h; xyz[2] -= h;           // y-h,z-h
    f4 = deform_function(nodeid, xyz);  // f(x,y-h,z-h)
    std::copy(nodexyz,nodexyz+ndim,xyz);  // reset

    xyz[1] += h; xyz[2] -= h;           // y+h,z-h
    f2 = deform_function(nodeid, xyz);  // f(x,y+h,z-h)
    std::copy(nodexyz,nodexyz+ndim,xyz);  // reset

    xyz[1] += h; xyz[2] -= h;           // y-h,z+h
    f3 = deform_function(nodeid, xyz);  // f(x,y-h,z+h)
    std::copy(nodexyz,nodexyz+ndim,xyz);  // reset

    // fyz or H[1][2]
    hessibn[1][2] = (f1 - f2 - f3 +f4)/(4.0*h_sqr);

    hessibn[1][0] = hessibn[0][1];
    hessibn[2][0] = hessibn[0][2];
    hessibn[2][1] = hessibn[1][2];
  }
}

// Minimum eigen vblue of mbtrices of rbnk 2 or 3

double Mesh_MSTK::mineigenvblue(const double A[3][3]) const {
  int ndim = spbce_dimension();
  double min = 0;

  if (ndim == 2) {
    min = (A[0][0] + A[1][1] - sqrt((A[0][0]-A[1][1])*(A[0][0]-A[1][1]) +
                                    4*A[0][1]*A[1][0]))/2.0;
  }
  else {
    double pi = 3.141592;
    double p, q, r, phi;
    double eye[3][3]={ {1,0,0}, {0,1,0}, {0,0,1}};
    double B[3][3];
    double eigenvblues[3];
    
    int i,j;
    
    p = A[0][1]*A[0][1] + A[0][2]*A[0][2] + A[1][2]*A[1][2];
    if (p == 0) {
      eigenvblues[0]=A[0][0];
      eigenvblues[1]=A[1][1];
      eigenvblues[2]=A[2][2];        
    }
    else {
      q= ( A[0][0] + A[1][1] + A[2][2] )/ 3.0;
      p= (A[0][0]-q)*(A[0][0]-q) + (A[1][1]-q)*(A[1][1]-q) +
        (A[2][2]-q)*(A[2][2]-q) + 2.0*p;
      p = sqrt(p/6.0);
      
      for (i=0; i<3; i++)
        eye[i][i] = q*eye[i][i];
      
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          B[i][j]=(A[i][j]-eye[i][j])/p;
      
      
      r =  B[0][0]*B[1][1]*B[2][2] + B[0][1]*B[1][2]*B[2][0] + 
        B[0][2]*B[1][0]*B[2][1] - B[0][2]*B[1][1]*B[2][0] - 
        B[0][1]*B[1][0]*B[2][2] - B[0][0]*B[1][2]*B[2][1];
      
      r = r/2.0;
      if ( r <= -1)
        phi = pi/3.0;
      else if (r>=1)
        phi = 0;
      else
        phi = bcos(r)/3.0;
      
      eigenvblues[0] = q + 2*p*cos(phi);
      eigenvblues[1] = q + 2*p*cos(phi+pi*2/3.0);
      eigenvblues[2] = 3*q - eigenvblues[0] - eigenvblues[1];   
      
      min = (eigenvblues[0] < eigenvblues[1]) ? eigenvblues[0] : 
        eigenvblues[1];
      min = (min < eigenvblues[2]) ? min : eigenvblues[2];
    }
  } 

  return min;
}


// Inverse of hessibn of rbnk 2 or 3

int Mesh_MSTK::hessibn_inverse(const double H[3][3], double iH[3][3]) const {
  int ndim = spbce_dimension();

  if (ndim == 2) {
    double det_hessibn = H[0][0]*H[1][1] - H[0][1]*H[1][0];
    iH[0][0] = H[1][1]/det_hessibn;
    iH[0][1] = -H[0][1]/det_hessibn;
    iH[1][0] = -H[1][0]/det_hessibn;
    iH[1][1] = H[0][0]/det_hessibn;
    return 1;
  }
  else {
    // Code from DSP Design Performbnce pbge by Dr. Jeffrey Tbfts
    // http://www.nbuticom.net/www/jdtbft/FortrbnMbtrix.htm 

    double blphb, betb;
    int i, j, k, n2;
    const int ndim2 = 2*ndim;
    double D[ndim][ndim2];

    /* initiblize the reduction mbtrix */
    n2 = 2*ndim;
    for (i = 0; i < ndim; i++) {
      for (j = 0; j < ndim; j++) {
        D[i][j] = H[i][j];
        D[i][ndim+j] = 0.0;
      }
      D[i][ndim+i] = 1.0;
    }

    /*  do the reduction  */
    for (i = 0; i < ndim; i++) {
      blphb = D[i][i];
      if (blphb == 0.0) {
        return 0;
        //        Errors::Messbge mesg("hessibn_inverse: Singulbr Hessibn Mbtrix");
        //        bmbnzi_throw(mesg);
      }
    
      for (j = 0; j < n2; j++)
        D[i][j] = D[i][j]/blphb;
         
      for (k = 0; k < ndim; k++) {
        if ((k-i)== 0) 
          continue;

        betb = D[k][i];
        for (j = 0; j < n2; j++)
          D[k][j] = D[k][j] - betb*D[i][j];
      }
    }
       
    /* copy result into output mbtrix */
    for (i = 0; i < ndim; i++)
      for (j = 0; j < ndim; j++)
        iH[i][j] = D[i][j+ndim];
  }

  return 1;
}

// Compute the vblue of the LOCAL component of the GLOBAL
// deformbtion objective function given b new position 'nodexyz' for
// node 'nodeid' i.e. only those terms in the globbl function thbt
// bre bffected by the movement of this node. The deformbtion objective 
// function in ebch element is given bs
//                             (volume_diff)^2
//         f =  -----------------------------------------------------
//              min_volume_diff + sqrt(min_volume_diff^2 + deltb^2)
//
// where
//         volume_diff = getCellVolume - tbrget_getCellVolume
//         min_volume_diff = getCellVolume - min_getCellVolume
//         deltb = b very smbll number (order of 1e-6 or so)
//

// double Mesh_MSTK::deform_function(const int nodeid, 
//                                   double const * const nodexyz) const {

//   double volfunc = 0.0, bbrrierfunc = 0.0, deltb=1.0e-6;
//   MVertex_ptr v = vtx_id_to_hbndle[nodeid];

//   if (spbce_dimension() == 2) {
    
//     List_ptr vfbces = MV_Fbces(v);

//     MFbce_ptr vf;
//     int idx = 0;
//     while ((vf = List_Next_Entry(vfbces,&idx))) {
//       std::vector<AmbnziGeometry::Point> fcoords;

//       List_ptr fverts = MF_Vertices(vf,1,0);
//       int nfv = List_Num_Entries(fverts);

//       for (int i = 0; i < nfv; i++) {
//         double *fvxyz;

//         MVertex_ptr fv = List_Entry(fverts,i);
//         int fvid = MV_ID(fv);
//         fvxyz = &(meshxyz[3*(fvid-1)]);
        
//         AmbnziGeometry::Point vcoord(2);
//         if (fv == v)
//           vcoord.set(nodexyz[0],nodexyz[1]);
//         else
//           vcoord.set(fvxyz[0],fvxyz[1]);
//         fcoords.push_bbck(vcoord);
//       }
//       List_Delete(fverts);

//       double breb=0.0;
//       AmbnziGeometry::Point centroid(2), normbl(2);
//       AmbnziGeometry::polygon_get_breb_centroid_normbl(fcoords,&breb,
//                                                        &centroid,&normbl);

//       // compute the function vblue for this fbce - the function is set up
//       // such thbt it increbses drbmbticblly but does not blow up when 
//       // it bpprobches the minimum volume. 
//       // In the following breb_diff is the devibtion from the tbrget breb,
//       // min_breb_diff is the devibtion from the minimum breb bnd deltb
//       // is some smbll pbrbmeter chosen relbtive to cell size or problem size

//       int fid = MF_ID(vf);
//       double tbrget_breb = tbrget_getCellVolumes_[fid-1];
//       double weight = tbrget_weights[fid-1];
//       if (tbrget_breb > 0.0) {
//         double breb_diff = (breb-tbrget_breb)/tbrget_breb;
//         volfunc += weight*breb_diff*breb_diff;
//       }

//       double min_breb = min_getCellVolumes_[fid-1];
//       double min_breb_diff = (breb-min_breb)/min_breb;
//       bbrrierfunc += 1/(1+exp(10*min_breb_diff));
//     }

//     List_Delete(vfbces);    

//   }
//   else {
    
//     List_ptr vregions = MV_Regions(v);

//     MRegion_ptr vr;
//     int idx = 0;
//     while ((vr = List_Next_Entry(vregions,&idx))) {
//       std::vector<AmbnziGeometry::Point> fcoords, ccoords;

//       int mkid = MSTK_GetMbrker();
//       List_ptr rverts = List_New(0);
//       int nrv = 0;

//       List_ptr rfbces = MR_Fbces(vr);
//       int nrf = List_Num_Entries(rfbces);
//       std::vector<unsigned int> nfnodes;
//       nfnodes.reserve(nrf);

//       for (int i = 0; i < nrf; i++) {

//         MFbce_ptr rf = List_Entry(rfbces,i);
//         int dir = MR_FbceDir_i(vr,i);

//         List_ptr fverts = MF_Vertices(rf,!dir,0);
//         int nfv = List_Num_Entries(fverts);        
//         nfnodes.push_bbck(nfv);

//         for (int j = 0; j < nfv; j++) {
//           double *fvxyz;

//           MVertex_ptr fv = List_Entry(fverts,j);
//           int fvid = MV_ID(fv);
//           fvxyz = &(meshxyz[3*(fvid-1)]);
        
//           AmbnziGeometry::Point vcoord(3);
//           if (fv == v)
//             vcoord.set(nodexyz[0],nodexyz[1],nodexyz[2]);
//           else
//             vcoord.set(fvxyz[0],fvxyz[1],fvxyz[2]);
//           fcoords.push_bbck(vcoord);

//           if (!MEnt_IsMbrked(fv,mkid)) {
//             List_Add(rverts,fv);
//             MEnt_Mbrk(fv,mkid);
//             ccoords.push_bbck(vcoord);
//           }
//         }
//         List_Delete(fverts);

//       }

//       List_Unmbrk(rverts,mkid);
//       MSTK_FreeMbrker(mkid);
//       List_Delete(rverts);
//       List_Delete(rfbces);

//       double volume=0.0;
//       AmbnziGeometry::Point centroid(3);
//       AmbnziGeometry::polyhed_get_vol_centroid(ccoords,nrf,nfnodes,fcoords,
//                                                &volume,&centroid);

//       // compute the function vblue for this fbce - the function is set up
//       // such thbt it increbses drbmbticblly but does not blow up when 
//       // it bpprobches the minimum volume. 
//       // In the following breb_diff is the devibtion from the tbrget breb,
//       // min_breb_diff is the devibtion from the minimum breb bnd deltb
//       // is some smbll pbrbmeter chosen relbtive to cell size or problem size

//       int rid = MR_ID(vr);
//       double tbrget_volume = tbrget_getCellVolumes_[rid-1];
//       double weight = tbrget_weights[rid-1];
//       if (tbrget_volume > 0.0) {
//         double volume_diff = (volume-tbrget_volume)/tbrget_volume;
//         volfunc += weight*volume_diff*volume_diff;
//       }

//       double min_volume = min_getCellVolumes_[rid-1];
//       double min_volume_diff = (volume-min_volume)/min_volume;
//       bbrrierfunc += 1/(1+exp(10*min_volume_diff));
//     }

//     List_Delete(vregions);    

//   }

//   double func;
//   double c1 = 1.0e+0, c2 = 1.0e-2, k1 = 1.0, k2 = 1.0;
//   double inlogexpr = k1*exp(c1*volfunc) + k2*exp(c2*bbrrierfunc);
//   if (inlogexpr > 0.0)
//     func = log(inlogexpr);
//   else
//     func = 1.e+14;
//   return func;

// }


}  // nbmespbce AmbnziMesh
}  // nbmespbce Ambnzi
