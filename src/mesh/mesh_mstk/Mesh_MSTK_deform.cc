/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella, others
*/

#include "dbc.hh"
#include "errors.hh"
#include "Geometry.hh"

#include "Mesh_MSTK.hh"

using namespace std;

namespace Amanzi {
namespace AmanziMesh {

// Some constants to be used by routines in this file only

double c1 = 1.0, c2 = 1.0, c3 = 1.0, k1 = 1.0, k2 = 1.0, k3 = 0.0;


// Main deformation driver

// IMPORTANT NOTE: SINCE THIS IS A SPECIALIZED ROUTINE THAT IS
// DESIGNED FOR DEFORMING LAYERED, COLUMNAR MESHES, IT CAN BE
// REASONABLY EXPECTED THAT ALL ELEMENTS AND NODES OF A COLUMN WILL BE
// ON A SINGLE PROCESSOR AND THAT ANY PARTITIONING IS ONLY IN THE
// LATERAL DIRECTIONS. IF THIS ASSUMPTION IS VIOLATED, THE CELLS
// ABOVE THE DEFORMED CELLS MAY NOT SHIFT DOWN


int
Mesh_MSTK::deform(const std::vector<double>& target_cell_volumes_in,
                  const std::vector<double>& min_cell_volumes_in,
                  const Entity_ID_List& fixed_nodes,
                  const bool move_vertical,
                  const double min_vol_const1,
                  const double min_vol_const2,
                  const double target_vol_const1,
                  const double target_vol_const2,
                  const double quality_func_const1,
                  const double quality_func_const2)
{
  int idx, id;
  MVertex_ptr mv;
  const int ndim = space_dimension();
  const int celldim = manifold_dimension();
  double eps = 1.0e-06;
  double damping_factor = 0.25;
  static const double macheps = 2.2e-16;
  static const double sqrt_macheps = sqrt(macheps);

  // This is a specialized function designed for columnar meshes so make sure
  // that we built the columns first
  if (!build_columns()) {
    std::cerr << "Could not deform mesh as we could not build columns in mesh\n" << std::endl;
    return 0;
  }

  // Initialize the deformation function constants

  k1 = min_vol_const1;
  c1 = min_vol_const2;
  k2 = target_vol_const1;
  c2 = target_vol_const2;
  k3 = quality_func_const1;
  c3 = quality_func_const2;

  if (!move_vertical) {
    Errors::Message mesg("Only vertical movement permitted in deformation");
    amanzi_throw(mesg);
  }

  int nv = num_entities(NODE, Parallel_type::ALL);

  // ---- Begin code by ETC ----
  int fixedmk = MSTK_GetMarker();
  List_ptr fixed_verts = List_New(0);
  for (Entity_ID_List::const_iterator v = fixed_nodes.begin(); v != fixed_nodes.end(); ++v) {
    if (*v < nv) {
      mv = vtx_id_to_handle[*v];
      if (!MEnt_IsMarked(mv, fixedmk)) {
        List_Add(fixed_verts, mv);
        MEnt_Mark(mv, fixedmk);
      }
    }
  }
  // ---- End code by ETC ----


  // store mesh node coordinates in temporary linear array to save on
  // the cost of retrieving each coordinate and decrease cache misses

  // Simultaneously get the extents of the problem
  meshxyz = new double[3 * nv]; // class variable
  double* newxyz = new double[3 * nv];

  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh_, &idx))) {
    id = MV_ID(mv);
    MV_Coords(mv, &(meshxyz[3 * (id - 1)]));
  }
  std::copy(meshxyz, meshxyz + 3 * nv, newxyz);


  // copy the target and min volumes for cells from the input

  int nc = num_entities(CELL, Parallel_type::ALL);
  target_cell_volumes_ = new double[nc]; // class variable
  min_cell_volumes_ = new double[nc];    // class variable

  std::copy(&(target_cell_volumes_in[0]), &(target_cell_volumes_in[nc]), target_cell_volumes_);
  std::copy(&(min_cell_volumes_in[0]), &(min_cell_volumes_in[nc]), min_cell_volumes_);


  // if the target cell volume is the same as the current volume, then
  // assume that the cell volume is unconstrained down to the min volume

  for (int i = 0; i < nc; ++i) {
    if (target_cell_volumes_[i] > 0.0) {
      double vol = cell_volume(i);
      if (vol == target_cell_volumes_[i]) target_cell_volumes_[i] = 0.0;
    }
  }

  // Mark vertices that must be driven to move to match target volumes
  // of their connected cells. The remaining vertices can be ignored
  // (they can still move because of the movement of of some nodes
  // directly below them)

  List_ptr driven_verts = List_New(10);
  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh_, &idx))) {
    if (MEnt_IsMarked(mv, fixedmk)) continue; // vertex is forbidden from moving

    List_ptr vcells = (celldim == 2) ? MV_Faces(mv) : MV_Regions(mv);
    int nvcells = List_Num_Entries(vcells);
    for (int i = 0; i < nvcells; i++) {
      MEntity_ptr ent = List_Entry(vcells, i);
      id = MEnt_ID(ent);

      if (target_cell_volumes_[id - 1] > 0.0) {
        // At least one cell connected to node/vertex needs to meet
        // a target vol. Mark the node/vertex as one that must move
        List_Add(driven_verts, mv);
        break;
      }
    }
    List_Delete(vcells);
  }

  // Now start mesh deformation to match the target volumes

  // Since this is such an unconstrained problem with many local
  // minima we will optimize it in a Jacobi fashion (objective
  // function evaluations are based on old values of mesh coordinates)
  // so that all nodes of a cell get a chance to move in order to
  // achieve the target volume. Additionally, we will damp the update
  // by 0.5 so as to keep the mesh valid. Finally, we will do only one
  // iteration of a local optimization at each vertex within a global
  // iteration over the entire mesh. Global iterations, of course, are
  // carried out until convergence

  int iter_global = 0, maxiter_global = 100;
  int converged_global = 0;

  while (!converged_global && iter_global < maxiter_global) {
    double global_dist2 = 0.0;
    double meshsizesqr_sum = 0.0;

    idx = 0;
    while ((mv = List_Next_Entry(driven_verts, &idx))) {
      id = MV_ID(mv);

      // Get an estimate of the mesh size around this vertex to
      // form a tolerance for node movement

      double minlen2 = 1.0e+14;
      List_ptr vedges = MV_Edges(mv);
      MEdge_ptr ve;
      int idx2 = 0;
      while ((ve = List_Next_Entry(vedges, &idx2))) {
        double elen2 = ME_LenSqr(ve);
        if (elen2 < minlen2) minlen2 = elen2;
      }
      List_Delete(vedges);

      meshsizesqr_sum += minlen2;

      double vxyzcur[ndim], vxyzold[ndim], vxyznew[ndim];

      for (int i = 0; i < ndim; i++) {
        vxyzcur[i] = newxyz[3 * (id - 1) + i];
        vxyzold[i] = vxyzcur[i];
        vxyznew[i] = 0.0;
      }

      // Initial value of function at start of local optimization

      double fcur = deform_function(id - 1, vxyzcur);

      int iter_local = 0, maxiter_local = 1;
      bool node_moved = false;
      while (iter_local < maxiter_local) {
        iter_local++;

        // ******** Get Z component of the gradient ********

        double gradient_z = 0.0;

        double f0 = fcur, f1, f2;

        double xyz[3];
        for (int i = 0; i < ndim; i++) xyz[i] = vxyzcur[i];

        // compute perturbation to z coordinate
        double h = sqrt_macheps * fabs(xyz[ndim - 1]);
        h = (h > 0.1 * sqrt_macheps) ? h : 0.1 * sqrt_macheps;

        // perturb coordinate and compute new function value
        xyz[ndim - 1] += h;
        f1 = deform_function(id - 1, xyz);
        xyz[ndim - 1] = vxyzcur[ndim - 1];

        gradient_z = (f1 - f0) / h;

        if (fabs(gradient_z) < eps * eps) continue; // too small

        if (gradient_z < 0.0) continue; // don't allow upward movement


        // ********* Get the ZZ component of the Hessian ****

        double hessian_zz = 0.0;

        h = 100 * sqrt_macheps * fabs(xyz[ndim - 1]);
        h = (h > 100 * sqrt_macheps) ? h : 100 * sqrt_macheps;

        xyz[ndim - 1] += h;
        f1 = deform_function(id - 1, xyz);
        xyz[ndim - 1] = vxyzcur[ndim - 1];

        xyz[ndim - 1] -= h;
        f2 = deform_function(id - 1, xyz);
        xyz[ndim - 1] = vxyzcur[ndim - 1];

        hessian_zz = (f2 - 2 * f0 + f1) / (h * h);


        // ********** Now the update **********************
        for (int i = 0; i < ndim; i++) vxyznew[i] = vxyzcur[i];

        double update[3] = { 0.0, 0.0, 0.0 };
        update[ndim - 1] -= damping_factor * gradient_z / hessian_zz;
        vxyznew[ndim - 1] = vxyzcur[ndim - 1] + update[ndim - 1];


        // ****** value of function with updated coordinates

        double fnew = deform_function(id - 1, vxyznew);


        // Check the Armijo condition (the step size should
        // produce a decrease in the function)

        double alpha = 1.0;
        while ((fnew >= fcur) && (alpha >= 1.0e-02)) {
          vxyznew[ndim - 1] = vxyzcur[ndim - 1] + alpha * update[ndim - 1];

          fnew = deform_function(id - 1, vxyznew);
          if (fnew >= fcur) alpha *= 0.8;
        }

        if (alpha < 1.0e-02) { // if movement is too little, dont move at all
          for (int i = 0; i < ndim; i++) vxyznew[i] = vxyzcur[i];
          fnew = fcur;
        } else {
          node_moved = true;

          // Update current coordinates to the newly calculated coordinates
          for (int i = 0; i < ndim; i++) vxyzcur[i] = vxyznew[i];
        }
      } // while (iter_local < maxiter_local)

      if (!node_moved) continue;

      double delta_z = vxyzcur[ndim - 1] - vxyzold[ndim - 1];
      global_dist2 += delta_z * delta_z;

      // update the new xyz values for this vertex
      for (int i = 0; i < 3; i++) newxyz[3 * (id - 1) + i] = vxyzcur[i];

      // Also drop the nodes above if the objective function
      // components corresponding to those nodes do not increase

      // IMPORTANT NOTE: SINCE THIS IS A SPECIALIZED ROUTINE THAT IS
      // DESIGNED FOR DEFORMING LAYERED, COLUMNAR MESHES, IT CAN BE
      // REASONABLY EXPECTED THAT ALL ELEMENTS AND NODES OF A COLUMN
      // WILL BE ON A SINGLE PROCESSOR. THEREFORE, WHEN A NODE DROPS
      // TO MATCH A TARGET VOLUME, ALL NODES ABOVE IT ARE EXPECTED
      // TO BE ON THE SAME PROCESSOR AND CAN BE MOVED DOWN. IF THIS
      // ASSUMPTION IS NOT TRUE, THEN THE UPDATE OF THE COLUMN
      // HAS TO HAPPEN DIFFERENTLY

      int curid = id - 1;
      bool done = false;
      while (!done) {
        int nextid = this->node_get_node_above(curid);
        if (nextid == -1) {
          done = true;
          continue;
        }
        double vxyz_above[3];
        std::copy(&(newxyz[3 * nextid]), &(newxyz[3 * nextid + ndim]), vxyz_above);

        double fcur_above = deform_function(nextid, vxyz_above);

        double alpha = 1.0;
        vxyz_above[ndim - 1] += alpha * delta_z; // 'add' delta_z (which is -ve)

        double fnew_above = deform_function(nextid, vxyz_above);

        // If the new function value w.r.t. nextid increases
        // cut back the proposed update as necessary

        while (fnew_above > fcur_above && alpha > 1.0e-02) {
          vxyz_above[ndim - 1] = newxyz[3 * nextid + ndim - 1];

          alpha *= 0.8;
          vxyz_above[ndim - 1] = newxyz[3 * nextid + ndim - 1] + alpha * delta_z;

          fnew_above = deform_function(nextid, vxyz_above);
        }

        if (alpha <= 1.0e-02) vxyz_above[ndim - 1] = newxyz[3 * nextid + ndim - 1];

        // Update the coordinates

        newxyz[3 * nextid + ndim - 1] = vxyz_above[ndim - 1];

        curid = nextid;
      }

    } // while (mv = ...)

    // an iteration over all vertices is over
    // update the current coordinates to be the newly computed coordinates

    std::copy(newxyz, newxyz + 3 * nv, meshxyz);

    iter_global++;

    // Update coordinates in mesh data structures - must update all
    // vertices not just vertices explicitly driven by optimization
    // since the downward movement of driven vertices can cause the
    // vertices above to shift down.

    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh_, &idx))) {
      id = MV_ID(mv);
      MV_Set_Coords(mv, &(meshxyz[3 * (id - 1)]));
    }

    // Update ghost vertex values for parallel runs

    if (get_comm()->NumProc() > 1) MESH_UpdateVertexCoords(mesh_, mpicomm_);

    double meshsize_rms = sqrt(meshsizesqr_sum / nv);
    eps = 1.0e-06 * meshsize_rms;

    double global_normdist2 = global_dist2 / nv;
    int converged_onthisproc = 0;
    if (global_normdist2 < eps * eps) converged_onthisproc = 1;
    get_comm()->MinAll(&converged_onthisproc, &converged_global, 1);

  } // while (!converged_global)

  //  if (iter_global >= maxiter_global)
  //    std::cerr << "INCREASE MAXITER_GLOBAL....!!" << std::endl;

  delete[] meshxyz;
  delete[] target_cell_volumes_;
  delete[] min_cell_volumes_;

  List_Unmark(fixed_verts, fixedmk);
  List_Delete(fixed_verts);
  MSTK_FreeMarker(fixedmk);

  List_Delete(driven_verts);

  if (space_dimension() == 2)
    MESH_ExportToGMV(mesh_, "deformed2.gmv", 0, NULL, NULL, MPI_COMM_NULL);
  else
    MESH_ExportToGMV(mesh_, "deformed3.gmv", 0, NULL, NULL, MPI_COMM_NULL);


  // recompute all geometric quantities

  compute_cell_geometric_quantities_();
  if (faces_initialized) compute_face_geometric_quantities_();
  if (edges_initialized) compute_edge_geometric_quantities_();

  return 1;
}


// Support functions for deformation function computation

void
VDiff3(double* a, double* b, double* c)
{
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
}

void
VSum3(double* a, double* b, double* c)
{
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
}

double
VLen3(double* a)
{
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

double
VLenSqr3(double* a)
{
  return (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

void
VNormalize3(double* a)
{
  double len = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

#ifdef DEBUG
  if (len <= 1.0e-28) fprintf(stderr, "Zero Length vector\n");
#endif

  a[0] /= len;
  a[1] /= len;
  a[2] /= len;
}

double
VDot3(double* a, double* b)
{
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

void
VCross3(double* a, double* b, double* x)
{
  x[0] = a[1] * b[2] - a[2] * b[1];
  x[1] = a[2] * b[0] - a[0] * b[2];
  x[2] = a[0] * b[1] - a[1] * b[0];
}


double
Tet_Volume(double* xyz0, double* xyz1, double* xyz2, double* xyz3)
{
  double vec1[3], vec2[3], vec3[3], cpvec[3];

  VDiff3(xyz1, xyz0, vec1);
  VDiff3(xyz2, xyz0, vec2);
  VDiff3(xyz3, xyz0, vec3);
  VCross3(vec1, vec2, cpvec);
  double vol = VDot3(cpvec, vec3) / 6.0;
  return vol;
}


double
Poly_Area(int n, double (*xyz)[3])
{
  int i;
  double area = 0.0;

  if (n < 3) return 0;

  /* This is based on Green's Theorem in the Plane - Works for all
     3D polygons

     Area = 0.5*Sum_over_i(a_i);
     a_i = x(i)*y(i+1)-x(i+1)*y(i);

     However, if the coordinates are very large, then a*b-c*d can
     result in roundoff error. To improve accuracy, we will make
     all coordinates relative to x0,y0. But what if xi is very close
     to x0? Then xi-x0 will also generate high error. Which is better?
  */

  if (fabs(xyz[0][0]) > 1000.0 || fabs(xyz[0][1]) > 1000.0) {
    for (i = 0; i < n; i++)
      area += ((xyz[i][0] - xyz[0][0]) * (xyz[(i + 1) % n][1] - xyz[0][1]) -
               (xyz[(i + 1) % n][0] - xyz[0][0]) * (xyz[i][1] - xyz[0][1]));
  } else {
    for (i = 0; i < n; i++)
      area += ((xyz[i][0]) * (xyz[(i + 1) % n][1]) - (xyz[(i + 1) % n][0]) * (xyz[i][1]));
  }

  area = 0.5 * area;

  return area;
}

double
func_modcn_corner3d(double* xyz0, double* xyz1, double* xyz2, double* xyz3)
{
  double evec0[3], evec1[3], evec2[3], cpvec[3], a, b, vol6, val;
  double loc_delta;

  VDiff3(xyz1, xyz0, evec0);
  VDiff3(xyz2, xyz0, evec1);
  VDiff3(xyz3, xyz0, evec2);

  a = VLenSqr3(evec0) + VLenSqr3(evec1) + VLenSqr3(evec2);
  VCross3(evec0, evec2, cpvec);
  b = VDot3(cpvec, cpvec);
  VCross3(evec1, evec2, cpvec);
  b += VDot3(cpvec, cpvec);
  VCross3(evec0, evec1, cpvec);
  b += VDot3(cpvec, cpvec);
  vol6 = VDot3(cpvec, evec2);
  loc_delta = 1.0 / (1 + exp(vol6));
  val = 2 * sqrt(a * b) / (vol6 + sqrt(vol6 * vol6 + loc_delta * loc_delta));
  return val;
}


// Compute the value of the LOCAL component of the GLOBAL
// deformation objective function given a new position 'nodexyz' for
// node 'nodeid' i.e. only those terms in the global function that
// are affected by the movement of this node.
//

////////////////////////////////////////////////////////////////////////////
// THIS VERSION USES CACHING OF DATA - USE IT IF THE ABOVE CODE IS TOO SLOW
////////////////////////////////////////////////////////////////////////////

double
Mesh_MSTK::deform_function(const int nodeid, double const* const nodexyz) const
{
  double func = 0.0;
  MVertex_ptr v = vtx_id_to_handle[nodeid];
  double* xyz[MAXPV3];
  // double evec0[3], evec1[3];
  double nodexyz_copy[3] = { 0.0, 0.0, 0.0 };
  double condfunc = 0.0, volfunc = 0.0, barrierfunc = 0.0;
  int i, j, k, id, jr, jf, found;
  MVertex_ptr fv, rv;
  int vind = -1;
  List_ptr fvlist, rvlist;
  static MVertex_ptr last_v = NULL;
  static MVertex_ptr(*fverts)[MAXPV2], (*rverts)[MAXPV3];
  static int vnbrsgen[MAXPV3], nnbrsgen;
  static int nf, *nfv, nr, *nrv, maxf, maxr, *nrf, (*nrfv)[MAXPF3];
  static int(*fvid)[MAXPV2], *vfids;
  static int(*rvid)[MAXPV3], *vrids;
  static int(*rfvlocid)[MAXPF3][MAXPV2];

  static int use_subtets_4all = 1;
  static int use_face_centers = 1;

  /* edge connected neighbors of vertices */
  // static int tetidx[4][3] = {{1,2,3},{2,0,3},{0,1,3},{1,0,2}};
  // static int hexidx[8][3] = {{1,3,4},{2,0,5},{3,1,6},{0,2,7},{0,5,7},{6,4,1},
  //                            {7,5,2},{6,4,3}};
  // static int prsmidx[6][3] = {{1,2,3},{2,0,4},{0,1,5},{0,4,5},{1,3,5},{2,3,4}};

  std::copy(nodexyz, nodexyz + space_dimension(), nodexyz_copy);

  if (v != last_v) {
    /* Collect face data */


    if (space_dimension() == 2) {
      List_ptr flist = MV_Faces(v);
      nf = List_Num_Entries(flist);

      if (last_v == NULL) {
        maxf = nf;
        fverts = (MVertex_ptr(*)[MAXPV2])malloc(maxf * sizeof(MVertex_ptr[MAXPV2]));
        nfv = (int*)malloc(maxf * sizeof(int));
        fvid = (int(*)[MAXPV2])malloc(maxf * sizeof(int[MAXPV2]));
        vfids = (int*)malloc(maxf * sizeof(int));
      } else if (nf > maxf) {
        maxf = 2 * nf;
        fverts = (MVertex_ptr(*)[MAXPV2])realloc(fverts, maxf * sizeof(MVertex_ptr[MAXPV2]));
        nfv = (int*)realloc(nfv, maxf * sizeof(int));
        fvid = (int(*)[MAXPV2])realloc(fvid, maxf * sizeof(int[MAXPV2]));
        vfids = (int*)realloc(vfids, maxf * sizeof(int));
      }

      for (k = 0, jf = 0; k < nf; k++) {
        MFace_ptr f = List_Entry(flist, k);

        vfids[jf] = MF_ID(f);
        fvlist = MF_Vertices(f, 1, v);
        nfv[jf] = List_Num_Entries(fvlist);

        for (i = 0; i < nfv[jf]; i++) {
          fverts[jf][i] = List_Entry(fvlist, i);
          fvid[jf][i] = MV_ID(fverts[jf][i]);
        }
        List_Delete(fvlist);

        jf++;
      }
      List_Delete(flist);
    } else {
      /* Collect region data */

      List_ptr rlist = MV_Regions(v);
      nr = List_Num_Entries(rlist);

      if (last_v == NULL) {
        maxr = nr;
        nrv = (int*)malloc(maxr * sizeof(int));
        rverts = (MVertex_ptr(*)[MAXPV3])malloc(maxr * sizeof(MVertex_ptr[MAXPV3]));
        rvid = (int(*)[MAXPV3])malloc(maxr * sizeof(int[MAXPV3]));
        nrf = (int*)malloc(maxr * sizeof(int));
        rfvlocid = (int(*)[MAXPF3][MAXPV2])malloc(maxr * sizeof(int[MAXPF3][MAXPV2]));
        nrfv = (int(*)[MAXPF3])malloc(maxr * sizeof(int[MAXPF3]));
        vrids = (int*)malloc(maxr * sizeof(int));
      } else if (nr > maxr) {
        maxr = 2 * nr;
        rverts = (MVertex_ptr(*)[MAXPV3])realloc(rverts, maxr * sizeof(MVertex_ptr[MAXPV3]));
        nrv = (int*)realloc(nrv, maxr * sizeof(int));
        rvid = (int(*)[MAXPV3])realloc(rvid, maxr * sizeof(int[MAXPV3]));
        nrf = (int*)realloc(nrf, maxr * sizeof(int));
        rfvlocid = (int(*)[MAXPF3][MAXPV2])realloc(rfvlocid, maxr * sizeof(int[MAXPF3][MAXPV2]));
        nrfv = (int(*)[MAXPF3])realloc(nrfv, maxr * sizeof(int[MAXPF3]));
        vrids = (int*)realloc(vrids, maxr * sizeof(int));
      }

      nnbrsgen = 0;
      for (jr = 0; jr < maxr; jr++) nrv[jr] = 0;
      for (jr = 0; jr < maxr; jr++)
        for (jf = 0; jf < MAXPF3; jf++) nrfv[jr][jf] = 0;

      for (jr = 0; jr < nr; jr++) {
        MRegion_ptr r = List_Entry(rlist, jr);

        vrids[jr] = MR_ID(r);
        rvlist = MR_Vertices(r);
        nrv[jr] = List_Num_Entries(rvlist);

        for (i = 0; i < nrv[jr]; i++) {
          rverts[jr][i] = List_Entry(rvlist, i);
          rvid[jr][i] = MV_ID(rverts[jr][i]);
        }

        /* collect face information also if not tet, prism or hex */

        List_ptr rflist = MR_Faces(r);
        nrf[jr] = List_Num_Entries(rflist);

        if (nrv[jr] == 4 || (!use_subtets_4all &&
                             ((nrv[jr] == 6 && nrf[jr] == 5) || (nrv[jr] == 8 && nrf[jr] == 6)))) {
          List_Delete(rvlist);
          continue;
        }


        for (jf = 0; jf < nrf[jr]; jf++) {
          MFace_ptr rf = List_Entry(rflist, jf);
          int rfdir = MR_FaceDir_i(r, jf);

          List_ptr rfvlist = MF_Vertices(rf, !rfdir, 0);
          nrfv[jr][jf] = List_Num_Entries(rfvlist);

          int nfv2 = nrfv[jr][jf];
          for (j = 0; j < nfv2; j++) {
            int rvidx;
            rv = List_Entry(rfvlist, j);

            /* find local index of face vertex in element vertex list */
            for (k = 0, found = 0; k < nrv[jr] && !found; k++)
              if (rv == rverts[jr][k]) {
                rfvlocid[jr][jf][j] = k;
                rvidx = k;
                found = 1;
              }
            if (!found) {
              MSTK_Report("MESH_MSTK::deform_function",
                          "Cannot find face vertex in region vertex list",
                          MSTK_WARN);
              Errors::Message mesg("FATAL error in Mesh_MSTK::deform_function");
              amanzi_throw(mesg);
            }

            /* if this vertex is a neighbor of v, add it to neighbor list */

            if (v == List_Entry(rfvlist, (j + 1) % nfv2) ||
                v == List_Entry(rfvlist, (j - 1 + nfv2) % nfv2)) {
              for (k = 0, found = 0; !found && k < nnbrsgen; k++)
                found = (rverts[vnbrsgen[k]] == rv) ? 1 : 0;
              if (!found) vnbrsgen[nnbrsgen++] = rvidx;
            }
          }
          List_Delete(rfvlist);
        } /* for jf = each face in element */

        List_Delete(rflist);
        List_Delete(rvlist);
      }
      List_Delete(rlist);
    }

    last_v = v;
  }


  // Actual computation of objective function

  double negvol_factor = 1.0; // will turn to -1 if -ve vol detected
  if (space_dimension() == 2) {
    for (jf = 0; jf < nf; jf++) {
      int fid = vfids[jf];

      for (i = 0; i < nfv[jf]; i++) {
        fv = fverts[jf][i];
        if (fv != v) {
          id = fvid[jf][i];
          xyz[i] = &(meshxyz[3 * (id - 1)]);
        } else
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

      //   a = L10_sqr+L20_sqr;

      //   if (space_dimension() == 3) {
      //     /* if we are dealing with volume mesh or a pure surface mesh
      //        then it is hard to define what inverted means on the
      //        surface. So we will use an unsigned area */

      //     double areavec[3];
      //     VCross3(evec0,evec1,areavec);
      //     A = VLen3(areavec);
      //   }
      //   else {
      //     /* This is a planar mesh - compute the signed area */

      //     A =  0.5*( xyz1[2][0]*xyz1[0][1] + xyz1[0][0]*xyz1[1][1] +
      //                xyz1[1][0]*xyz1[2][1] - xyz1[1][0]*xyz1[0][1] -
      //                xyz1[2][0]*xyz1[1][1] - xyz1[0][0]*xyz1[2][1]);
      //   }
      //   double delta = 1.0/(1+exp(A));
      //   condfunc += 2*a/(A+sqrt(A*A+delta*delta));

      // }

      double(*pxyz)[3] = (double(*)[3])malloc(nfv[jf] * sizeof(double[3]));
      for (k = 0; k < nfv[jf]; k++) std::copy(xyz[k], xyz[k] + 3, pxyz[k]);
      double face_area = Poly_Area(nfv[jf], pxyz);
      free(pxyz);

      // this face contributes to the objective function,
      // only if a +ve target area was requested for this face,

      double target_area = target_cell_volumes_[fid - 1];
      if (target_area > 0) {
        double area_diff = (face_area - target_area) / target_area;
        volfunc += area_diff * area_diff;
      }

      // every face always contributes a barrier function to
      // keep it from going below a certain area

      double min_area = min_cell_volumes_[fid - 1];
      double min_area_diff = (face_area - min_area) / min_area;

      double bfunc = 1 / exp(c1 * min_area_diff);
      barrierfunc += bfunc;
    }

  } else {
    for (jr = 0; jr < nr; jr++) {
      int rid = vrids[jr];
      double region_volume = 0;

      if ((nrv[jr] == 4) || (!use_subtets_4all &&
                             ((nrv[jr] == 6 && nrf[jr] == 5) || (nrv[jr] == 8 && nrf[jr] == 6)))) {
        /* tet, or (prism or hex using CN at their original corners) */

        for (i = 0; i < nrv[jr]; i++) {
          rv = rverts[jr][i];
          if (rv != v) {
            id = rvid[jr][i];
            xyz[i] = &(meshxyz[3 * (id - 1)]);
          } else {
            xyz[i] = nodexyz_copy;
            vind = i;
          }
        }

        /* indices of vertices at which condition number must be computed */
        /* This is the index of the vertex and its edge connected neighbors */
        /* Also make general list of neighbors of the involved vertices */

        // ind[0] = vind;
        // switch (nrv[jr]) {
        // case 4:
        //   for (i = 0; i < 3; i++)
        //     ind[i+1] = tetidx[vind][i];
        //   for (i = 0; i < 4; i++) {
        //     j = ind[i];
        //     for (k = 0; k < 3; k++)
        //       nbrs[i][k] = tetidx[j][k];
        //   }
        //   break;
        // case 6:
        //   for (i = 0; i < 3; i++)
        //     ind[i+1] = prsmidx[vind][i];
        //   for (i = 0; i < 4; i++) {
        //     j = ind[i];
        //     for (k = 0; k < 3; k++)
        //       nbrs[i][k] = prsmidx[j][k];
        //   }
        //   break;
        // case 8:
        //   for (i = 0; i < 3; i++)
        //     ind[i+1] = hexidx[vind][i];
        //   for (i = 0; i < 4; i++) {
        //     j = ind[i];
        //     for (k = 0; k < 3; k++)
        //       nbrs[i][k] = hexidx[j][k];
        //   }
        //   break;
        // }

        // for (i = 0; i < 4; i++) {
        //   /* For each relevant corner ind[i], compute condition number */

        //   /* coordinates of vertex ind[i] itself */
        //   xyz1[0] = xyz[ind[i]];

        //   /* coordinates of its edge connected neighbors */
        //   for (j = 0; j < 3; j++) {
        //     int m = nbrs[i][j]; /* j'th neighbor of ind[i] */
        //     xyz1[j+1] = xyz[m];
        //   }

        //   condfunc += func_modcn_corner3d(xyz1[0], xyz1[1], xyz1[2], xyz1[3]);
        // }

        // FOR NOW ASSUME THAT USE_SUBTETS_4ALL IS TRUE
        // AND SO CODE WILL COME HERE ONLY FOR TETS
        region_volume += Tet_Volume(xyz[0], xyz[1], xyz[2], xyz[3]);
      } else {
        /* general polyhedron with possibly non-trivalent corners or
           hex/prism for which we will use a tet decomposition */

        /* we will virtually decompose the polyhedron into tetrahedra in
           a symmetric manner (one "central point" on each polyhedron
           face and one "central point" inside the polyhedron, connected
           with an edge to form a tetrahedron). Then we will compute
           condition numbers for these tetrahedra. This will be a
           consistent way of dealing with non-trivalent corners. It will
           also ensure that, in the end, every tetrahedron in the
           polyhedral decomposition is valid and well shaped */

        /* compute the central point of the element */
        double rcen[3] = { 0.0, 0.0, 0.0 };
        for (i = 0; i < nrv[jr]; i++) {
          rv = rverts[jr][i];
          if (rv != v) {
            id = rvid[jr][i];
            xyz[i] = &(meshxyz[3 * (id - 1)]);
          } else
            xyz[i] = nodexyz_copy;

          for (k = 0; k < 3; k++) rcen[k] += xyz[i][k];
        }
        for (k = 0; k < 3; k++) rcen[k] /= nrv[jr];

        /* compute central point of faces containing v in element */

        for (jf = 0; jf < nrf[jr]; jf++) {
          int nfv2 = nrfv[jr][jf];

          double fcen[3] = { 0.0, 0.0, 0.0 };
          if (use_face_centers) {
            for (j = 0; j < nfv2; j++) {
              int rvlid = rfvlocid[jr][jf][j];
              for (k = 0; k < 3; k++) fcen[k] += xyz[rvlid][k];
            }
            for (k = 0; k < 3; k++) fcen[k] /= nfv2;
          }

          for (j = 0; j < nfv2; j++) {
            int rvlid = rfvlocid[jr][jf][j];
            int rvlid1;

            /* Volume of this tet */

            rvlid1 = rfvlocid[jr][jf][(j + 1) % nfv2];
            double tet_volume = Tet_Volume(xyz[rvlid], xyz[rvlid1], fcen, rcen);
            negvol_factor = (tet_volume < 0.0) ? -1 : negvol_factor;
            region_volume += tet_volume;

            /* Condition number at corners of this tet if it is relevant */
            /* Is face vertex either v or one of its edge connected nbrs? */

            // found = 0;
            // if (rverts[rvlid] == v) found = 1;
            // for (k = 0; !found && k < nnbrsgen; k++)
            //   if (rvlid == vnbrsgen[k]) found = 1;

            // if (!found) continue;

            // if (use_face_centers == 1) {
            //   /* form virtual tets using the central face point, central
            //      element point and edges of the element - compute condition
            //      numbers at tet corners affected by the movement of v. Which
            //      condition numbers to compute is subject to debate since the
            //      central element point and central face points are also
            //      affected by movement of v. However, to keep it simple, we
            //      will compute the condition numbers only at vertices we would
            //      have computed them at if we were not using a virtual
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
        }   // jf
      }     // else


      // IF ANY OF THE COMPONENT TETS HAD A NEGATIVE VOLUME, MAKE SURE THAT
      // THE REGION VOLUME IS MADE NEGATIVE (IT MAY HAVE BEEN CANCELLED OUT
      // BY A MULTITUDE OF POSITIVE VOLUME TETS).

      region_volume = negvol_factor * fabs(region_volume);

      // If the target_volume is 0, it indicates that we don't care
      // about driving the volume towards the target_volume, so make
      // this function value 0

      double target_volume = target_cell_volumes_[rid - 1];
      double volume_diff =
        (target_volume > 0) ? (region_volume - target_volume) / target_volume : 0.0;
      volfunc += volume_diff * volume_diff;

      // every region always contributes a barrier function to
      // keep it from going below a certain region

      double min_volume = min_cell_volumes_[rid - 1];
      double min_volume_diff = (region_volume - min_volume) / min_volume;

      double bfunc = 1 / exp(c1 * min_volume_diff);
      barrierfunc += bfunc;
    } // for each region
  }

  // If volfunc == 0.0, it means that none of the cells connected to
  // this node need to be driven towards a particular target value. So
  // there is no point in returning a function value based solely on
  // differences between current volume and minimum volume as this
  // might actually cause node movement. In such a case, make the in
  // logexpr 1.0, which will result in a 0 function value

  // Of course, the better (and more cost effective) thing to do is to
  // tag such nodes which don't need to be explicitly moved and not
  // try to move them at all

  double inlogexpr =
    volfunc > 0.0 ? k1 * exp(c1 * barrierfunc) + k2 * exp(c2 * volfunc) + k3 * exp(c3 * condfunc) :
                    1;
  func = (inlogexpr > 0.0 && negvol_factor > 0) ? log(inlogexpr) : 1e+10;
  return func;
}


// Finite difference gradient of deformation objective function

void
Mesh_MSTK::deform_gradient(const int nodeid, double const* const nodexyz, double* gradient) const
{
  static const double macheps = 2.2e-16;
  static const double sqrt_macheps = sqrt(macheps);
  int ndim = space_dimension();

  double xyz[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < ndim; i++) xyz[i] = nodexyz[i];

  double f0 = deform_function(nodeid, xyz);

  // For now we assume we are getting a full 2D/3D gradient for every vertex

  for (int i = 0; i < ndim; i++) {
    gradient[i] = 0.0;

    // compute perturbation to i'th coordinate
    double h = sqrt_macheps * fabs(xyz[i]);
    h = (h > 0.1 * sqrt_macheps) ? h : 0.1 * sqrt_macheps;

    // perturb coordinate and compute new function value
    xyz[i] += h;
    double f1 = deform_function(nodeid, xyz);

    // Forward difference computation of gradient
    gradient[i] = (f1 - f0) / h;

    // Restore i'th coordinate
    xyz[i] = nodexyz[i];
  }
}

// Finite difference hessian of deformation objective function in 2D

void
Mesh_MSTK::deform_hessian(const int nodeid, double const* const nodexyz, double hessian[3][3]) const
{
  static const double macheps = 2.2e-16;
  static const double sqrt_macheps = sqrt(macheps);
  double h, h_sqr;
  const int ndim = space_dimension();

  /* Since we have mixed derivatives, it is not clear what the 'x'
     in the formula sqrt(macheps) times x should be - take an
     average of all the coordinates w.r.t. which we will take
     derivatives */
  /* Also, if we are dealing with very small numbers, eps2_sqr
     can be below machine precision. To guard against this
     lets multiply the initial estimate of eps2 by 100 and if
     it is still less than 100*sqrt_macheps, reject the estimate
     and just use 10*sqrt_macheps */

  if (ndim == 2) {
    double xyz[3] = { 0.0, 0.0, 0.0 };
    std::copy(nodexyz, nodexyz + ndim, xyz);

    double f0, f1, f2, f3, f4, f5, f6;

    h = 100 * sqrt_macheps * (fabs(xyz[0]) + fabs(xyz[1])) / 2.0;
    h = (h > 100 * sqrt_macheps) ? h : 100 * sqrt_macheps;
    h_sqr = h * h;


    f0 = deform_function(nodeid, xyz); // f(x,y)

    xyz[0] += h;                             // x+h
    f1 = deform_function(nodeid, xyz);       // f(x+h,y)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    xyz[0] -= h;                             // x-h
    f2 = deform_function(nodeid, xyz);       // f(x-h,y)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    // fxx or H[0][0]
    hessian[0][0] = (f2 - 2 * f0 + f1) / h_sqr;


    xyz[1] += h;                             // y+h
    f3 = deform_function(nodeid, xyz);       // f(x,y+h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    xyz[1] -= h;                             // y-h
    f4 = deform_function(nodeid, xyz);       // f(x,y-h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    // fyy or H[1][1]
    hessian[1][1] = (f3 - 2 * f0 + f4) / h_sqr;


    xyz[0] += h;
    xyz[1] += h;                             // x+h, y+h
    f5 = deform_function(nodeid, xyz);       // f(x+h,y+h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    xyz[0] -= h;
    xyz[1] -= h;                             // x-h, y-h
    f6 = deform_function(nodeid, xyz);       // f(x-h,y-h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    // fxy or H[0][1]
    hessian[0][1] = (f5 - f1 - f3 + 2 * f0 - f2 - f4 + f6) / (2 * h_sqr);

    hessian[1][0] = hessian[0][1];
  } else {
    /* Since we have mixed derivatives, it is not clear what the 'x'
       in the formula sqrt(macheps) times x should be - take an
       average of all the coordinates w.r.t. which we will take
       derivatives */
    /* Also, if we are dealing with very small numbers, eps2_sqr
       can be below machine precision. To guard against this
       lets multiply the initial estimate of eps2 by 100 and if
       it is still less than 100*sqrt_macheps, reject the estimate
       and just use 10*sqrt_macheps */

    double xyz[3] = { 0.0, 0.0, 0.0 };
    std::copy(nodexyz, nodexyz + ndim, xyz);

    double f0, f1, f2, f3, f4;

    /* Since we have mixed derivatives, it is not clear what the 'x'
       in the formula sqrt(macheps) times x should be - take an
       average of all the coordinates w.r.t. which we will take
       derivatives */

    h = 100 * sqrt_macheps * (fabs(xyz[0]) + fabs(xyz[1]) + fabs(xyz[2])) / 3.0;
    h = (h < 100 * sqrt_macheps) ? 100 * sqrt_macheps : h;
    h_sqr = h * h;


    f0 = deform_function(nodeid, xyz); // f(x,y,z)

    xyz[0] += h;                             // x+h
    f1 = deform_function(nodeid, xyz);       // f(x+h,y,z)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    xyz[0] -= h;                             // x-h
    f2 = deform_function(nodeid, xyz);       // f(x-h,y,z)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    // fxx or H[0][0]
    hessian[0][0] = (f2 - 2 * f0 + f1) / h_sqr;


    xyz[1] += h;                             // y+h
    f1 = deform_function(nodeid, xyz);       // f(x,y+h,z)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    xyz[1] -= h;                             // y-h
    f2 = deform_function(nodeid, xyz);       // f(x,y-h,z)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    // fyy or H[1][1]
    hessian[1][1] = (f2 - 2 * f0 + f1) / h_sqr;


    xyz[0] += h;
    xyz[1] += h;                             // x+h,y+h
    f1 = deform_function(nodeid, xyz);       // f(x+h,y+h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    xyz[0] -= h;
    xyz[1] -= h;                             // x-y,y-h
    f4 = deform_function(nodeid, xyz);       // f(x-h,y-h,z)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    xyz[0] += h;
    xyz[1] -= h;                             // x+h,y-h
    f2 = deform_function(nodeid, xyz);       // f(x+h,y-h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    xyz[0] -= h;
    xyz[1] += h;                             // x-h,y+h
    f3 = deform_function(nodeid, xyz);       // f(x-h,y+h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    // fxy
    hessian[0][1] = (f1 - f2 - f3 + f4) / (4.0 * h_sqr);


    xyz[2] += h;                             // z+h
    f1 = deform_function(nodeid, xyz);       // f(x,y,z+h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    xyz[2] -= h;                             // z-h
    f2 = deform_function(nodeid, xyz);       // f(x,y,z-h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    // fzz or H[2][2]
    hessian[2][2] = (f2 - 2 * f0 + f1) / h_sqr;


    xyz[0] += h;
    xyz[2] += h;                             // x+h,y,z+h
    f1 = deform_function(nodeid, xyz);       // f(x+h,y,z+h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    xyz[0] -= h;
    xyz[2] -= h;                             // x-h,z-h
    f4 = deform_function(nodeid, xyz);       // f(x-h,y,z-h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    xyz[0] += h;
    xyz[2] -= h;                             // x+h,y,z-h
    f2 = deform_function(nodeid, xyz);       // f(x+h,y,z-h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    xyz[0] -= h;
    xyz[2] += h;                             // x-h,z+h
    f3 = deform_function(nodeid, xyz);       // f(x-h,y,z+h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    // fxz or H[0][2]
    hessian[0][2] = (f1 - f2 - f3 + f4) / (4.0 * h_sqr);


    xyz[1] += h;
    xyz[2] += h;                             // y+h,z+h
    f1 = deform_function(nodeid, xyz);       // f(x,y+h,z+h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    xyz[1] -= h;
    xyz[2] -= h;                             // y-h,z-h
    f4 = deform_function(nodeid, xyz);       // f(x,y-h,z-h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    xyz[1] += h;
    xyz[2] -= h;                             // y+h,z-h
    f2 = deform_function(nodeid, xyz);       // f(x,y+h,z-h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    xyz[1] += h;
    xyz[2] -= h;                             // y-h,z+h
    f3 = deform_function(nodeid, xyz);       // f(x,y-h,z+h)
    std::copy(nodexyz, nodexyz + ndim, xyz); // reset

    // fyz or H[1][2]
    hessian[1][2] = (f1 - f2 - f3 + f4) / (4.0 * h_sqr);

    hessian[1][0] = hessian[0][1];
    hessian[2][0] = hessian[0][2];
    hessian[2][1] = hessian[1][2];
  }
}

// Minimum eigen value of matrices of rank 2 or 3

double
Mesh_MSTK::mineigenvalue(const double A[3][3]) const
{
  int ndim = space_dimension();
  double min = 0;

  if (ndim == 2) {
    min = (A[0][0] + A[1][1] -
           sqrt((A[0][0] - A[1][1]) * (A[0][0] - A[1][1]) + 4 * A[0][1] * A[1][0])) /
          2.0;
  } else {
    double pi = 3.141592;
    double p, q, r, phi;
    double eye[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
    double B[3][3];
    double eigenvalues[3];

    int i, j;

    p = A[0][1] * A[0][1] + A[0][2] * A[0][2] + A[1][2] * A[1][2];
    if (p == 0) {
      eigenvalues[0] = A[0][0];
      eigenvalues[1] = A[1][1];
      eigenvalues[2] = A[2][2];
    } else {
      q = (A[0][0] + A[1][1] + A[2][2]) / 3.0;
      p = (A[0][0] - q) * (A[0][0] - q) + (A[1][1] - q) * (A[1][1] - q) +
          (A[2][2] - q) * (A[2][2] - q) + 2.0 * p;
      p = sqrt(p / 6.0);

      for (i = 0; i < 3; i++) eye[i][i] = q * eye[i][i];

      for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) B[i][j] = (A[i][j] - eye[i][j]) / p;


      r = B[0][0] * B[1][1] * B[2][2] + B[0][1] * B[1][2] * B[2][0] + B[0][2] * B[1][0] * B[2][1] -
          B[0][2] * B[1][1] * B[2][0] - B[0][1] * B[1][0] * B[2][2] - B[0][0] * B[1][2] * B[2][1];

      r = r / 2.0;
      if (r <= -1)
        phi = pi / 3.0;
      else if (r >= 1)
        phi = 0;
      else
        phi = acos(r) / 3.0;

      eigenvalues[0] = q + 2 * p * cos(phi);
      eigenvalues[1] = q + 2 * p * cos(phi + pi * 2 / 3.0);
      eigenvalues[2] = 3 * q - eigenvalues[0] - eigenvalues[1];

      min = (eigenvalues[0] < eigenvalues[1]) ? eigenvalues[0] : eigenvalues[1];
      min = (min < eigenvalues[2]) ? min : eigenvalues[2];
    }
  }

  return min;
}


// Inverse of hessian of rank 2 or 3

int
Mesh_MSTK::hessian_inverse(const double H[3][3], double iH[3][3]) const
{
  int ndim = space_dimension();

  if (ndim == 2) {
    double det_hessian = H[0][0] * H[1][1] - H[0][1] * H[1][0];
    iH[0][0] = H[1][1] / det_hessian;
    iH[0][1] = -H[0][1] / det_hessian;
    iH[1][0] = -H[1][0] / det_hessian;
    iH[1][1] = H[0][0] / det_hessian;
    return 1;
  } else {
    // Code from DSP Design Performance page by Dr. Jeffrey Tafts
    // http://www.nauticom.net/www/jdtaft/FortranMatrix.htm

    double alpha, beta;
    int i, j, k, n2;
    const int ndim2 = 2 * ndim;
    double D[ndim][ndim2];

    /* initialize the reduction matrix */
    n2 = 2 * ndim;
    for (i = 0; i < ndim; i++) {
      for (j = 0; j < ndim; j++) {
        D[i][j] = H[i][j];
        D[i][ndim + j] = 0.0;
      }
      D[i][ndim + i] = 1.0;
    }

    /*  do the reduction  */
    for (i = 0; i < ndim; i++) {
      alpha = D[i][i];
      if (alpha == 0.0) {
        return 0;
        //        Errors::Message mesg("hessian_inverse: Singular Hessian Matrix");
        //        amanzi_throw(mesg);
      }

      for (j = 0; j < n2; j++) D[i][j] = D[i][j] / alpha;

      for (k = 0; k < ndim; k++) {
        if ((k - i) == 0) continue;

        beta = D[k][i];
        for (j = 0; j < n2; j++) D[k][j] = D[k][j] - beta * D[i][j];
      }
    }

    /* copy result into output matrix */
    for (i = 0; i < ndim; i++)
      for (j = 0; j < ndim; j++) iH[i][j] = D[i][j + ndim];
  }

  return 1;
}

// Compute the value of the LOCAL component of the GLOBAL
// deformation objective function given a new position 'nodexyz' for
// node 'nodeid' i.e. only those terms in the global function that
// are affected by the movement of this node. The deformation objective
// function in each element is given as
//                             (volume_diff)^2
//         f =  -----------------------------------------------------
//              min_volume_diff + sqrt(min_volume_diff^2 + delta^2)
//
// where
//         volume_diff = cell_volume - target_cell_volume
//         min_volume_diff = cell_volume - min_cell_volume
//         delta = a very small number (order of 1e-6 or so)
//

// double Mesh_MSTK::deform_function(const int nodeid,
//                                   double const * const nodexyz) const {

//   double volfunc = 0.0, barrierfunc = 0.0, delta=1.0e-6;
//   MVertex_ptr v = vtx_id_to_handle[nodeid];

//   if (space_dimension() == 2) {

//     List_ptr vfaces = MV_Faces(v);

//     MFace_ptr vf;
//     int idx = 0;
//     while ((vf = List_Next_Entry(vfaces,&idx))) {
//       std::vector<AmanziGeometry::Point> fcoords;

//       List_ptr fverts = MF_Vertices(vf,1,0);
//       int nfv = List_Num_Entries(fverts);

//       for (int i = 0; i < nfv; i++) {
//         double *fvxyz;

//         MVertex_ptr fv = List_Entry(fverts,i);
//         int fvid = MV_ID(fv);
//         fvxyz = &(meshxyz[3*(fvid-1)]);

//         AmanziGeometry::Point vcoord(2);
//         if (fv == v)
//           vcoord.set(nodexyz[0],nodexyz[1]);
//         else
//           vcoord.set(fvxyz[0],fvxyz[1]);
//         fcoords.push_back(vcoord);
//       }
//       List_Delete(fverts);

//       double area=0.0;
//       AmanziGeometry::Point centroid(2), normal(2);
//       AmanziGeometry::polygon_get_area_centroid_normal(fcoords,&area,
//                                                        &centroid,&normal);

//       // compute the function value for this face - the function is set up
//       // such that it increases dramatically but does not blow up when
//       // it approaches the minimum volume.
//       // In the following area_diff is the deviation from the target area,
//       // min_area_diff is the deviation from the minimum area and delta
//       // is some small parameter chosen relative to cell size or problem size

//       int fid = MF_ID(vf);
//       double target_area = target_cell_volumes_[fid-1];
//       double weight = target_weights[fid-1];
//       if (target_area > 0.0) {
//         double area_diff = (area-target_area)/target_area;
//         volfunc += weight*area_diff*area_diff;
//       }

//       double min_area = min_cell_volumes_[fid-1];
//       double min_area_diff = (area-min_area)/min_area;
//       barrierfunc += 1/(1+exp(10*min_area_diff));
//     }

//     List_Delete(vfaces);

//   }
//   else {

//     List_ptr vregions = MV_Regions(v);

//     MRegion_ptr vr;
//     int idx = 0;
//     while ((vr = List_Next_Entry(vregions,&idx))) {
//       std::vector<AmanziGeometry::Point> fcoords, ccoords;

//       int mkid = MSTK_GetMarker();
//       List_ptr rverts = List_New(0);
//       int nrv = 0;

//       List_ptr rfaces = MR_Faces(vr);
//       int nrf = List_Num_Entries(rfaces);
//       std::vector<unsigned int> nfnodes;
//       nfnodes.reserve(nrf);

//       for (int i = 0; i < nrf; i++) {

//         MFace_ptr rf = List_Entry(rfaces,i);
//         int dir = MR_FaceDir_i(vr,i);

//         List_ptr fverts = MF_Vertices(rf,!dir,0);
//         int nfv = List_Num_Entries(fverts);
//         nfnodes.push_back(nfv);

//         for (int j = 0; j < nfv; j++) {
//           double *fvxyz;

//           MVertex_ptr fv = List_Entry(fverts,j);
//           int fvid = MV_ID(fv);
//           fvxyz = &(meshxyz[3*(fvid-1)]);

//           AmanziGeometry::Point vcoord(3);
//           if (fv == v)
//             vcoord.set(nodexyz[0],nodexyz[1],nodexyz[2]);
//           else
//             vcoord.set(fvxyz[0],fvxyz[1],fvxyz[2]);
//           fcoords.push_back(vcoord);

//           if (!MEnt_IsMarked(fv,mkid)) {
//             List_Add(rverts,fv);
//             MEnt_Mark(fv,mkid);
//             ccoords.push_back(vcoord);
//           }
//         }
//         List_Delete(fverts);

//       }

//       List_Unmark(rverts,mkid);
//       MSTK_FreeMarker(mkid);
//       List_Delete(rverts);
//       List_Delete(rfaces);

//       double volume=0.0;
//       AmanziGeometry::Point centroid(3);
//       AmanziGeometry::polyhed_get_vol_centroid(ccoords,nrf,nfnodes,fcoords,
//                                                &volume,&centroid);

//       // compute the function value for this face - the function is set up
//       // such that it increases dramatically but does not blow up when
//       // it approaches the minimum volume.
//       // In the following area_diff is the deviation from the target area,
//       // min_area_diff is the deviation from the minimum area and delta
//       // is some small parameter chosen relative to cell size or problem size

//       int rid = MR_ID(vr);
//       double target_volume = target_cell_volumes_[rid-1];
//       double weight = target_weights[rid-1];
//       if (target_volume > 0.0) {
//         double volume_diff = (volume-target_volume)/target_volume;
//         volfunc += weight*volume_diff*volume_diff;
//       }

//       double min_volume = min_cell_volumes_[rid-1];
//       double min_volume_diff = (volume-min_volume)/min_volume;
//       barrierfunc += 1/(1+exp(10*min_volume_diff));
//     }

//     List_Delete(vregions);

//   }

//   double func;
//   double c1 = 1.0e+0, c2 = 1.0e-2, k1 = 1.0, k2 = 1.0;
//   double inlogexpr = k1*exp(c1*volfunc) + k2*exp(c2*barrierfunc);
//   if (inlogexpr > 0.0)
//     func = log(inlogexpr);
//   else
//     func = 1.e+14;
//   return func;

// }


} // namespace AmanziMesh
} // namespace Amanzi
