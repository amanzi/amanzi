#include "Geometry.hh"
#include "dbc.hh"
#include "errors.hh"

#include "Mesh.hh"

using namespace std;

namespace Amanzi
{

namespace AmanziMesh
{

int Mesh::compute_geometric_quantities() const {

  int ncells = num_entities(CELL,USED);

  for (int i = 0; i < ncells; i++) {
    double volume;
    AmanziGeometry::Point centroid(spacedim);

    compute_cell_geometry(i,&volume,&centroid);

    cell_volumes.push_back(volume);
    cell_centroids.push_back(centroid);
  }
  
  
  int nfaces = num_entities(FACE,USED);

  for (int i = 0; i < nfaces; i++) {
    double area;
    AmanziGeometry::Point centroid(spacedim), normal0(spacedim), 
      normal1(spacedim);
    
    // normal0 and normal1 are outward normals of the face with
    // respect to the cell0 and cell1 of the face. The natural normal
    // of the face points out of cell0 and into cell1. If one of these
    // cells do not exist, then the normal is the null vector.

    compute_face_geometry(i,&area,&centroid,&normal0,&normal1);
  
    face_areas.push_back(area);
    face_centroids.push_back(centroid);
    face_normal0.push_back(normal0);  
    face_normal1.push_back(normal1);  
  }

  geometry_precomputed = true;

  return 1;

} // Mesh::compute_geometric_quantities



int Mesh::compute_cell_geometry(const Entity_ID cellid, double *volume, AmanziGeometry::Point *centroid) const {


  if (celldim == 3) {

    // 3D Elements with possibly curved faces
    // We have to build a description of the element topology
    // and send it into the polyhedron volume and centroid
    // calculation routine
    
    Entity_ID_List faces;
    std::vector<unsigned int> nfnodes;
    std::vector<int> fdirs;
    std::vector<AmanziGeometry::Point> ccoords, cfcoords, fcoords;
    
    cell_get_faces_and_dirs(cellid,&faces,&fdirs);
  
    int nf = faces.size();
  
    for (int j = 0; j < nf; j++) {
      
      face_get_coordinates(faces[j],&fcoords);
      nfnodes.push_back(fcoords.size());
      
      if (fdirs[j] == 1) {
        for (int k = 0; k < nfnodes[j]; k++)
          cfcoords.push_back(fcoords[k]);
      }
      else {
        for (int k = nfnodes[j]-1; k >=0; k--)
          cfcoords.push_back(fcoords[k]);
      }
    }
    
    cell_get_coordinates(cellid,&ccoords);
    
    AmanziGeometry::polyhed_get_vol_centroid(ccoords,nf,nfnodes,
                                             cfcoords,volume,
                                             centroid);
    return 1;
  }
  else if (celldim == 2) {

    std::vector<AmanziGeometry::Point> ccoords;
    
    cell_get_coordinates(cellid,&ccoords);

    AmanziGeometry::Point normal(spacedim);
    
    AmanziGeometry::polygon_get_area_centroid_normal(ccoords,volume,centroid,&normal);

    return 1;
  }

  return 0;
} // Mesh::compute_cell_geometry


  int Mesh::compute_face_geometry(const Entity_ID faceid, double *area, 
                                  AmanziGeometry::Point *centroid, 
                                  AmanziGeometry::Point *normal0, 
                                  AmanziGeometry::Point *normal1) const {

  AmanziGeometry::Point_List fcoords;

  (*normal0).set(0.0L);
  (*normal1).set(0.0L);

  if (celldim == 3) {

    // 3D Elements with possibly curved faces
    // We have to build a description of the element topology
    // and send it into the polyhedron volume and centroid
    // calculation routine

    face_get_coordinates(faceid,&fcoords);
      
    AmanziGeometry::Point normal(3);
    AmanziGeometry::polygon_get_area_centroid_normal(fcoords,area,centroid,&normal);      

    Entity_ID_List cellids;    
    face_get_cells(faceid, USED, &cellids);
    
    for (int i = 0; i < cellids.size(); i++) {
      Entity_ID_List cellfaceids;
      std::vector<int> cellfacedirs;
      int dir = 1;
      
      cell_get_faces_and_dirs(cellids[i], &cellfaceids, &cellfacedirs);

      bool found = false;
      for (int j = 0; j < cellfaceids.size(); j++) {
        if (cellfaceids[j] == faceid) {
          found = true;
          dir = cellfacedirs[j];
          break;
        }
      }

      ASSERT(found);

      if (dir == 1)
        *normal0 = normal;
      else
        *normal1 = -normal;
    }

    return 1;
  }
  else if (celldim == 2) {

    if (spacedim == 2) {   // 2D mesh

      face_get_coordinates(faceid,&fcoords);

      AmanziGeometry::Point evec = fcoords[1]-fcoords[0];
      *area = sqrt(evec*evec);
      
      *centroid = 0.5*(fcoords[0]+fcoords[1]);
      
      AmanziGeometry::Point normal = AmanziGeometry::Point(evec[1],-evec[0]);

      Entity_ID_List cellids;
      face_get_cells(faceid, USED, &cellids);
      
      for (int i = 0; i < cellids.size(); i++) {
        Entity_ID_List cellfaceids;
        std::vector<int> cellfacedirs;
        int dir = 1;
        
        cell_get_faces_and_dirs(cellids[i], &cellfaceids, &cellfacedirs);
        
        bool found = false;
        for (int j = 0; j < cellfaceids.size(); j++) {
          if (cellfaceids[j] == faceid) {
            found = true;
            dir = cellfacedirs[j];
            break;
          }
        }
        
        ASSERT(found);
        
        if (dir == 1)
          *normal0 = normal;
        else
          *normal1 = -normal;
      }
      
       return 1;
    }
    else {  // Surface mesh - cells are 2D, coordinates are 3D

      // edge normals are ambiguous for surface mesh
      // So we won't compute them

      face_get_coordinates(faceid,&fcoords);

      AmanziGeometry::Point evec = fcoords[1]-fcoords[0];
      *area = sqrt(evec*evec);

      *centroid = 0.5*(fcoords[0]+fcoords[1]);

      Entity_ID_List cellids;
      face_get_cells(faceid, USED, &cellids);
      
      for (int i = 0; i < cellids.size(); i++) {
        Entity_ID_List cellfaceids;
        std::vector<int> cellfacedirs;
        int dir = 1;
        
        cell_get_faces_and_dirs(cellids[i], &cellfaceids, &cellfacedirs);
        
        bool found = false;
        for (int j = 0; j < cellfaceids.size(); j++) {
          if (cellfaceids[j] == faceid) {
            found = true;
            dir = cellfacedirs[j];
            break;
          }
        }
        
        ASSERT(found);

        AmanziGeometry::Point cvec = fcoords[0]-cell_centroids[cellids[i]];
        AmanziGeometry::Point trinormal = cvec^evec;

        AmanziGeometry::Point normal = evec^trinormal;
        
        double len = norm(normal);
        normal /= len;
        normal *= *area;
 
        if (dir == 1)
          *normal0 = normal;
        else
          *normal1 = normal; // Note that we are not flipping the sign here
      }      
 
      return 1;
    }

  }

  return 0;

} // Mesh::compute_face_geometry



// Volume/Area of cell

double Mesh::cell_volume (const Entity_ID cellid, const bool recompute) const {

  if (!geometry_precomputed) {
    compute_geometric_quantities();
    return cell_volumes[cellid];
  }
  else {
    if (recompute) {
      double volume;
      AmanziGeometry::Point centroid(spacedim);
      compute_cell_geometry(cellid, &volume, &centroid);
      return volume;
    }
    else
      return cell_volumes[cellid];
  }
}
  
// Area/length of face
  
double Mesh::face_area(const Entity_ID faceid, const bool recompute) const {

  if (!geometry_precomputed) {
    compute_geometric_quantities();
    return face_areas[faceid];
  }
  else {
    if (recompute) {
      double area;
      AmanziGeometry::Point centroid(spacedim);
      AmanziGeometry::Point normal0(spacedim), normal1(spacedim);
      compute_face_geometry(faceid, &area, &centroid, &normal0, &normal1);
      return area;
    }
    else
      return face_areas[faceid];
  }
}
  

// Centroid of cell

AmanziGeometry::Point Mesh::cell_centroid (const Entity_ID cellid, const bool recompute) const {

  if (!geometry_precomputed) {
    compute_geometric_quantities();
    return cell_centroids[cellid];
  }
  else {
    if (recompute) {
      double volume;
      AmanziGeometry::Point centroid(spacedim);
      compute_cell_geometry(cellid, &volume, &centroid);
      return centroid;
    }
    else
      return cell_centroids[cellid];
  }
  
}

// Centroid of face

AmanziGeometry::Point Mesh::face_centroid (const Entity_ID faceid, const bool recompute) const {

  if (!geometry_precomputed) {
    compute_geometric_quantities();
    return face_centroids[faceid];
  }
  else {
    if (recompute) {
      double area;
      AmanziGeometry::Point centroid(spacedim);
      AmanziGeometry::Point normal0(spacedim), normal1(spacedim);
      compute_face_geometry(faceid, &area, &centroid, &normal0, &normal1);
      return centroid;
    }
    else
      return face_centroids[faceid];
  }
  
}

// Normal to face
// The vector is normalized and then weighted by the area of the face
//
// If recompute is TRUE, then the normal is recalculated using current
// face coordinates but not stored. (If the recomputed normal must be
// stored, then call recompute_geometric_quantities). 
//
// If cellid is not specified, the normal is the natural normal of the face
// If cellid is specified, the normal is the outward normal with respect
// to the cell. In planar and solid meshes, the normal with respect to 
// the cell on one side of the face is just the negative of the normal 
// with respect to the cell on the other side. In general surfaces meshes,
// this will not be true at C1 discontinuities


AmanziGeometry::Point Mesh::face_normal (const Entity_ID faceid, const bool recompute, const Entity_ID cellid) const {

  AmanziGeometry::Point normal0(spacedim), normal1(spacedim);
    
  if (!geometry_precomputed) {
    compute_geometric_quantities();

    normal0 = face_normal0[faceid];
    normal1 = face_normal1[faceid];
  }
  else {
    if (recompute) {
      double area;
      AmanziGeometry::Point centroid(spacedim);
      compute_face_geometry(faceid, &area, &centroid, &normal0, &normal1);
    }
    else {      
      normal0 = face_normal0[faceid];
      normal1 = face_normal1[faceid];
    }
  }

  if (cellid == -1) {
    // Just the natural normal of the face
    // Since normal0 and normal1 are outward facing normals with respect
    // to their respective cells, we can return normal0 as is but have
    // to negate normal1.

    if (L22(normal0) != 0.0)
      return normal0;
    else {
      ASSERT(L22(normal1) != 0.0);
      return -normal1;
    }
  }
  else {
    Entity_ID_List faceids;
    std::vector<int> face_dirs;

    cell_get_faces_and_dirs(cellid, &faceids, &face_dirs);

    int nf = faceids.size();
    bool found = false;
    int dir = 1;
    for (int i = 0; i < nf; i++)
      if (faceids[i] == faceid) {
        dir = face_dirs[i];
        found = true;
        break;
      }
    
    ASSERT(found);
    
    if (dir == 1) {
      ASSERT(L22(normal0) != 0.0);
      return normal0;
    }
    else {
      ASSERT(L22(normal1) != 0.0);
      return normal1;
    }
  }

  return normal0;
}



// Get set ID given the name of the set - return 0 if no match is found

Set_ID Mesh::set_id_from_name(const std::string setname) const
{
  if (!geometric_model_) return 0;

  unsigned int ngr = geometric_model_->Num_Regions();
  for (int i = 0; i < ngr; i++) {
    AmanziGeometry::RegionPtr rgn = geometric_model_->Region_i(i);
    
    if (rgn->name() == setname)
      return rgn->id();
  }

  return 0;
}

// Get set name given the ID of the set - return 0 if no match is found

std::string Mesh::set_name_from_id(const int setid) const
{
  std::string nullname("");
 
  if (!geometric_model_) return nullname;

  unsigned int ngr = geometric_model_->Num_Regions();
  for (int i = 0; i < ngr; i++) {
    AmanziGeometry::RegionPtr rgn = geometric_model_->Region_i(i);
    
    if (rgn->id() == setid)
      return rgn->name();
  }

  return 0;
}



// Is there a set with this id and entity type

bool Mesh::valid_set_id(Set_ID id, Entity_kind kind) const
{

  if (!geometric_model_) return false;

  unsigned int gdim = geometric_model_->dimension();    
  
  unsigned int ngr = geometric_model_->Num_Regions();
  for (int i = 0; i < ngr; i++) {
    AmanziGeometry::RegionPtr rgn = geometric_model_->Region_i(i);
    
    unsigned int rdim = rgn->dimension();
    
    if (rgn->id() == id) {
      
      // For regions of type Labeled Set and Color Function, the
      // dimension parameter is not guaranteed to be correct

      if (rgn->type() == AmanziGeometry::LABELEDSET ||
          rgn->type() == AmanziGeometry::COLORFUNCTION) return true;

      // If we are looking for a cell set the region has to be 
      // of the same topological dimension as the cells
      
      if (kind == CELL && rdim == celldim) return true;
      
      // If we are looking for a side set, the region has to be 
      // one topological dimension less than the cells
      
      if (kind == FACE && rdim == celldim-1) return true;
      
      // If we are looking for a node set, the region can be of any
      // dimension upto the spatial dimension of the domain
      
      if (kind == NODE) return true;
      
    }
  }

  return false;
}



// Is there a set with this name and entity type

bool Mesh::valid_set_name(std::string name, Entity_kind kind) const
{

  if (!geometric_model_) return false;

  unsigned int gdim = geometric_model_->dimension();    
  
  unsigned int ngr = geometric_model_->Num_Regions();
  for (int i = 0; i < ngr; i++) {
    AmanziGeometry::RegionPtr rgn = geometric_model_->Region_i(i);
    
    unsigned int rdim = rgn->dimension();
    
    if (rgn->name() == name) {

      // For regions of type Labeled Set and Color Function, the dimension 
      // parameter is not guaranteed to be correct

      if (rgn->type() == AmanziGeometry::LABELEDSET ||
          rgn->type() == AmanziGeometry::COLORFUNCTION) return true;

      // If we are looking for a cell set the region has to be 
      // of the same topological dimension as the cells or it
      // has to be a point region
      
      if (kind == CELL && (rdim >= celldim || rdim == 0)) return true;
      
      // If we are looking for a side set, the region has to be 
      // one topological dimension less than the cells
      
      if (kind == FACE && rdim >= celldim-1) return true;
      
      // If we are looking for a node set, the region can be of any
      // dimension upto the spatial dimension of the domain
      
      if (kind == NODE) return true;
      
    }
  }

  return false;
}








bool Mesh::point_in_cell(const AmanziGeometry::Point &p, const Entity_ID cellid) const
{
  std::vector<AmanziGeometry::Point> ccoords;

  if (celldim == 3) {

    // 3D Elements with possibly curved faces
    // We have to build a description of the element topology
    // and send it into the polyhedron volume and centroid
    // calculation routine
    
    int nf;
    Entity_ID_List faces;
    std::vector<unsigned int> nfnodes;
    std::vector<int> fdirs;
    std::vector<AmanziGeometry::Point> cfcoords;
    
    cell_get_faces_and_dirs(cellid,&faces,&fdirs);
    
    nf = faces.size();
    
    for (int j = 0; j < nf; j++) {
      std::vector<AmanziGeometry::Point> fcoords;

      face_get_coordinates(faces[j],&fcoords);
      nfnodes.push_back(fcoords.size());
      
      if (fdirs[j] == 1) {
        for (int k = 0; k < nfnodes[j]; k++)
          cfcoords.push_back(fcoords[k]);
      }
      else {
        for (int k = nfnodes[j]-1; k >=0; k--)
          cfcoords.push_back(fcoords[k]);
      }
    }
    
    cell_get_coordinates(cellid,&ccoords);
    
    return AmanziGeometry::point_in_polyhed(p,ccoords,nf,nfnodes,cfcoords);

  }
  else if (celldim == 2) {
    
    cell_get_coordinates(cellid,&ccoords);
    
    return AmanziGeometry::point_in_polygon(p,ccoords);
    
  }

  return false;
}


// Deform the mesh according to a given set of new node positions
// If keep_valid is true, the routine will cut back node displacement
// if the cells connected to a moved node become invalid

int Mesh::deform (const Entity_ID_List nodeids,
                  const AmanziGeometry::Point_List new_positions,
                  const bool keep_valid,
                  AmanziGeometry::Point_List *final_positions) {

  int status = 1;

  ASSERT(nodeids.size() == new_positions.size());
  ASSERT(final_positions != NULL);

  final_positions->clear();


  // Once we start moving nodes around, the precomputed/cached
  // geometric quantities are no longer valid. So any geometric calls
  // must use the "recompute=true" option until the end of this routine
  // where we once again call compute_geometric_quantities 

  int nn = nodeids.size();

  bool done_outer = false; 
  int iter = 0, maxiter = 5;

  while (!done_outer) {
    double totdisp2 = 0.0;

    for (int j = 0; j < nn; j++) {
      Entity_ID node = nodeids[j];
      
      AmanziGeometry::Point oldcoords, newcoords, dispvec;
      Entity_ID_List cells;
      
      node_get_coordinates(node,&oldcoords);
      dispvec = new_positions[j]-oldcoords;
      
      node_get_cells(node,USED,&cells);
      int nc = cells.size();
      
      double mult = 1.0;
      bool done = false;
      bool allvalid = true;
      
      while (!done) {
        
        newcoords = oldcoords + mult*dispvec;
        
        node_set_coordinates(node,newcoords);
        
        if (keep_valid) { // check if the cells remain valid
          allvalid = true;
          for (int k = 0; k < nc; k++)
            if (cell_volume(cells[k],true) < 0.0) {
              allvalid = false;
              break;
          }
        }
        
        if (allvalid)
          done = true;
        else {
          if (mult < 1.0e-5)
            done = true;
          mult = mult/2.0;
        }
      } // while (!done)
      
      if (!allvalid) { // could not move the node even a bit
        status = 0;    // perhaps the mesh was invalid to start?
        node_set_coordinates(node,oldcoords);
      }
      else {
        AmanziGeometry::Point actual_dispvec = newcoords - oldcoords;
        totdisp2 += L22(actual_dispvec);
      }
    } // for (j = 0; j < nn; j++) 

    if (totdisp2 < 1.0e-12)
      done_outer = 1;

    if (++iter == maxiter)
      break;
  } // while (!done_outer)


  for (int j = 0; j < nn; j++) {
    Entity_ID node = nodeids[j];
    
    AmanziGeometry::Point newcoords;

    node_get_coordinates(node,&newcoords);
    final_positions->push_back(newcoords);
  }


  // recompute all geometric quantities

  compute_geometric_quantities();

  return status;
}


// Figure out columns of cells in a structured mesh and cache the
// information for later.

// The code currently makes the assumption that the "bottom" of the is
// a flat surface in the XY plane. It then builds up information about
// the cell above and cell below for each cell based on the
// orientation of the face normals w.r.t the z direction. If the mesh
// is highly warped, this could lead to ambiguities. Also,
// intersecting columns in an unstructured mesh will lead to an
// exception being thrown. These data structures are never populated
// if these operators are never called. The above and below cells are
// computed for all cells the first time one of the routines,
// cell_get_cellabove or cell_get_cellbelow are called and then this
// info is cached



int Mesh::build_columns() const {
  int status = 1;

  // Find the faces at the bottom of the domain. For now we assume
  // that these are all the boundary faces whose normal points in the
  // negative z-direction

  int nf = num_entities(FACE,USED);
  int nc = num_entities(CELL,USED);

  // Initialize cell_below and cell_above so that we can assign using
  // cell_below[i] = j type operations below

  cell_cellbelow.resize(nc);
  cell_cellbelow.assign(nc,-1);
  cell_cellabove.resize(nc);
  cell_cellabove.assign(nc,-1);
  

  for (int i = 0; i < nf; i++) {
    
    Entity_ID_List fcells;
    face_get_cells(i,USED,&fcells); // Should this be USED or OWNED?

    if (fcells.size() != 1) continue;

    AmanziGeometry::Point normal = face_normal(i,false,fcells[0]);
    normal /= norm(normal);
      
    AmanziGeometry::Point negzvec(spacedim);
    if (spacedim == 2)
      negzvec.set(0.0,-1.0);
    else if (spacedim == 3)
      negzvec.set(0.0,0.0,-1.0);
    
    double dp = negzvec*normal;
    
    if (fabs(dp-1.0) > 1.0e-06) continue;

    // found a boundary face with a downward facing normal

    // Just to make sure we are not making a mistake, lets check that
    // the centroid of the cell is above the centroid of the face

    AmanziGeometry::Point ccen(spacedim),fcen(spacedim);
    ccen = cell_centroid(fcells[0]);
    fcen = face_centroid(i);

    AmanziGeometry::Point cfvec = fcen-ccen;
    cfvec /= norm(cfvec);

    dp = negzvec*cfvec;

    if (fabs(dp-1.0) > 1.0e-06) continue;

    // Now we are quite sure that this is a face at the bottom of the
    // mesh/domain

    // Walk through the cells until we get to the top of the domain

    Entity_ID cur_cell = fcells[0];
    Entity_ID bot_face = i;
    Entity_ID top_face = -1;
    Entity_ID_List fcells2, cfaces;
    std::vector<int> cfdirs;

    bool done = false;
    while (!done) {
      
      // Faces of current cell

      cell_get_faces_and_dirs(cur_cell,&cfaces,&cfdirs);

      // Find the top face of the cell as the face whose outward
      // normal from the current cell is most aligned with the Z
      // direction

      double mindp = 999.0;
      top_face = -1;
      for (int j = 0; j < cfaces.size(); j++) {
        normal = face_normal(cfaces[j]);
        if (cfdirs[j] == -1) normal *= -1;
        normal /= norm(normal);
        
        dp = normal*negzvec;
        if (dp < mindp) {
          mindp = dp;
          top_face = cfaces[j];
        }
      }


      assert(top_face != bot_face);

      face_get_cells(top_face,USED,&fcells2);
      if (fcells2.size() == 2) {

        if (cell_cellabove[cur_cell] != -1) {  // intersecting column of cells 
          status = 0;
          Errors::Message mesg("Intersecting column of cells");
          Exceptions::amanzi_throw(mesg);
        }

        if (fcells2[0] == cur_cell) {
          cell_cellabove[cur_cell] = fcells2[1];
          cell_cellbelow[fcells2[1]] = cur_cell;
          cur_cell = fcells2[1];
        }
        else if (fcells2[1] == cur_cell) {
          cell_cellabove[cur_cell] = fcells2[0];
          cell_cellbelow[fcells2[0]] = cur_cell;
          cur_cell = fcells2[0];
        }
        else {
          status = 0;
          Errors::Message mesg("Unlikely problem in face to cell connectivity");
          Exceptions::amanzi_throw(mesg);
        }

        bot_face = top_face;
      }
      else 
        done = true;

    } // while (!done)
      
  }

  columns_built = true;
  return status;
}


Entity_ID Mesh::cell_get_cell_above(const Entity_ID cellid) const {

  if (!columns_built)
    build_columns();

  return cell_cellabove[cellid];
}

Entity_ID Mesh::cell_get_cell_below(const Entity_ID cellid) const {
  
  if (!columns_built)
    build_columns();

  return cell_cellbelow[cellid];
}

} // close namespace AmanziMesh
} // close namespace Amanzi
