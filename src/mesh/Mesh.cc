#include "Geometry.hh"
#include "dbc.hh"

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
    AmanziGeometry::Point centroid(spacedim), normal(spacedim);
    
    compute_face_geometry(i,&area,&centroid,&normal);
  
    face_areas.push_back(area);
    face_centroids.push_back(centroid);
    face_normals.push_back(normal);
  }

  geometry_precomputed = true;

} // Mesh::precompute_geometric_quantities



int Mesh::compute_cell_geometry(const Entity_ID cellid, double *volume, AmanziGeometry::Point *centroid) const {


  if (celldim == 3) {

    // 3D Elements with possibly curved faces
    // We have to build a description of the element topology
    // and send it into the polyhedron volume and centroid
    // calculation routine
    
    std::vector<Entity_ID> faces;
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
    
    AmanziGeometry::polygon_get_area_centroid(ccoords,volume,centroid);

    return 1;
  }

  return 0;
} // Mesh::compute_cell_geometry


int Mesh::compute_face_geometry(const Entity_ID faceid, double *area, AmanziGeometry::Point *centroid, AmanziGeometry::Point *normal) const {

  AmanziGeometry::Point_List fcoords;

  if (celldim == 3) {

    // 3D Elements with possibly curved faces
    // We have to build a description of the element topology
    // and send it into the polyhedron volume and centroid
    // calculation routine

    face_get_coordinates(faceid,&fcoords);
      
    AmanziGeometry::polygon_get_area_centroid(fcoords,area,centroid);
      
    *normal = AmanziGeometry::polygon_get_normal(fcoords);
    *normal *= *area;

    return 1;
  }
  else if (celldim == 2) {

    if (spacedim == 2) {   // 2D mesh

      face_get_coordinates(faceid,&fcoords);

      AmanziGeometry::Point evec = fcoords[1]-fcoords[0];
      *area = sqrt(evec*evec);
      
      *centroid = 0.5*(fcoords[0]+fcoords[1]);
      
      *normal = AmanziGeometry::Point(evec[1],-evec[0]);

      return 1;
    }
    else {  // Surface mesh - cells are 2D, coordinates are 3D

      // edge normals are ambiguous for surface mesh
      // So we won't compute them

      face_get_coordinates(faceid,&fcoords);

      AmanziGeometry::Point evec = fcoords[1]-fcoords[0];
      *area = sqrt(evec*evec);

      *centroid = 0.5*(fcoords[0]+fcoords[1]);

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
      AmanziGeometry::Point centroid(spacedim), normal(spacedim);
      compute_face_geometry(faceid, &area, &centroid, &normal);
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
      AmanziGeometry::Point centroid(spacedim), normal(spacedim);
      compute_face_geometry(faceid, &area, &centroid, &normal);
      return area;
    }
    else
      return face_centroids[faceid];
  }
  
}

// Normal to face
// The vector is normalized and then weighted by the area of the face

AmanziGeometry::Point Mesh::face_normal (const Entity_ID faceid, const bool recompute) const {

  if (!geometry_precomputed) {
    compute_geometric_quantities();
    return face_normals[faceid];
  }
  else {
    if (recompute) {
      double area;
      AmanziGeometry::Point centroid(spacedim), normal(spacedim);
      compute_face_geometry(faceid, &area, &centroid, &normal);
      return area;
    }
    else
      return face_normals[faceid];
  }
  
}



// Get set ID given the name of the set - return 0 if no match is found

unsigned int Mesh::set_id_from_name(const std::string setname) const
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

bool Mesh::valid_set_id(unsigned int id, Entity_kind kind) const
{

  if (!geometric_model_) return false;

  unsigned int gdim = geometric_model_->dimension();    
  
  unsigned int ngr = geometric_model_->Num_Regions();
  for (int i = 0; i < ngr; i++) {
    AmanziGeometry::RegionPtr rgn = geometric_model_->Region_i(i);
    
    unsigned int rdim = rgn->dimension();
    
    if (rgn->id() == id) {
      
      // For regions of type Labeled Set, the dimension parameter is
      // not guaranteed to be correct

      if (rgn->type() == AmanziGeometry::LABELEDSET) return true;

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

      // For regions of type Labeled Set, the dimension parameter is
      // not guaranteed to be correct

      if (rgn->type() == AmanziGeometry::LABELEDSET) return true;

      // If we are looking for a cell set the region has to be 
      // of the same topological dimension as the cells or it
      // has to be a point region
      
      if (kind == CELL && (rdim == celldim || rdim == 0)) return true;
      
      // If we are looking for a side set, the region has to be 
      // one topological dimension less than the cells
      
      if (kind == FACE && rdim == celldim-1) return true;
      
      // If we are looking for a node set, the region can be of any
      // dimension upto the spatial dimension of the domain
      
      if (kind == NODE) return true;
      
    }
  }
  
}










// Temporary routines for backward compatibility
//
// ----- ATTENTION - ATTENTION - ATTENTION -----
// NEW CODE SHOULD NOT USE THESE FUNCTIONS SINCE THEY WILL GO AWAY
//



void Mesh::cell_to_faces (unsigned int cell,
                          std::vector<unsigned int>::iterator begin,
                          std::vector<unsigned int>::iterator end) const
{
  cell_to_faces (cell, &(*begin), &(*end));
};

void Mesh::cell_to_faces (unsigned int cell,
                          unsigned int* begin, unsigned int *end) const
{
  Entity_ID_List cfaces;

  cell_get_faces ((Entity_ID)cell, &cfaces);

  assert (cfaces.size() <= (unsigned int) (end-begin));
  std::copy (cfaces.begin(), cfaces.end(), begin);
};


void Mesh::cell_to_face_dirs (unsigned int cell,
                              std::vector<int>::iterator begin,
                              std::vector<int>::iterator end) const
{
  cell_to_face_dirs (cell, &(*begin), &(*end));
};
void Mesh::cell_to_face_dirs (unsigned int cell,
                              int * begin, int * end) const
{
  vector<int> cfdirs;

  cell_get_face_dirs ((Entity_ID)cell, &cfdirs);

  assert (cfdirs.size() <= (unsigned int) (end-begin));
  std::copy (cfdirs.begin(), cfdirs.end(), begin);
};



void Mesh::cell_to_nodes (unsigned int cell,
                          std::vector<unsigned int>::iterator begin,
                          std::vector<unsigned int>::iterator end) const
{
  cell_to_nodes (cell, &(*begin), &(*end));
};
void Mesh::cell_to_nodes (unsigned int cell,
                          unsigned int * begin, unsigned int * end) const
{
  Entity_ID_List cnodes;

  cell_get_nodes ((Entity_ID)cell, &cnodes);

  assert (cnodes.size() <= (unsigned int) (end-begin));
  std::copy (cnodes.begin(), cnodes.end(), begin);
};




void Mesh::face_to_nodes (unsigned int face,
                          std::vector<unsigned int>::iterator begin,
                          std::vector<unsigned int>::iterator end) const
{
  face_to_nodes(face, &(*begin), &(*end));
};
void Mesh::face_to_nodes (unsigned int face,
                          unsigned int * begin, unsigned int * end) const
{
  Entity_ID_List fnodes;

  face_get_nodes ((Entity_ID)face, &fnodes);

  assert (fnodes.size() <= (unsigned int) (end-begin));
  std::copy (fnodes.begin(), fnodes.end(), begin);
};



void Mesh::node_to_coordinates (unsigned int node,
                                std::vector<double>::iterator begin,
                                std::vector<double>::iterator end) const
{
  node_to_coordinates (node, &(*begin), &(*end));
};
void Mesh::node_to_coordinates (unsigned int node,
                                double * begin,
                                double * end) const
{
  AmanziGeometry::Point xyz;

  node_get_coordinates ((Entity_ID)node, &xyz);

  assert (xyz.dim() == (end-begin));
  for (int i = 0; i < xyz.dim(); i++) begin[i] = xyz[i];
};

void Mesh::face_to_coordinates (unsigned int face,
                                std::vector<double>::iterator begin,
                                std::vector<double>::iterator end) const
{
  face_to_coordinates (face, &(*begin), &(*end));
};
void Mesh::face_to_coordinates (unsigned int face,
                                double * begin,
                                double * end) const
{
  vector<AmanziGeometry::Point> fxyz;

  face_get_coordinates ((Entity_ID)face, &fxyz);

  int dim = fxyz[0].dim();
  int nfn = fxyz.size();
  assert ((nfn*dim) <= (end-begin));

  for (int i = 0; i < nfn; i++)
    for (int j = 0; j < dim; j++)
      begin[i*dim+j] = fxyz[i][j];
};

void Mesh::cell_to_coordinates (unsigned int cell,
                                std::vector<double>::iterator begin,
                                std::vector<double>::iterator end) const
{
  cell_to_coordinates (cell, &(*begin), &(*end));
};
void Mesh::cell_to_coordinates (unsigned int cell,
                                double * begin,
                                double * end) const
{
  vector<AmanziGeometry::Point> cxyz;

  cell_get_coordinates ((Entity_ID)cell, &cxyz);

  int dim = cxyz[0].dim();
  int ncn = cxyz.size();
  assert ((ncn*dim) <= (end-begin));

  for (int i = 0; i < ncn; i++)
    for (int j = 0; j < dim; j++)
      begin[i*dim+j] = cxyz[i][j];
};


unsigned int Mesh::count_entities (Entity_kind kind,
                                   Parallel_type ptype) const
{
  return num_entities(kind, ptype);
};



  // Number of sets containing entities of type 'kind' in mesh


  // Id numbers
void Mesh::get_set_ids (Entity_kind kind,
                        std::vector<unsigned int>::iterator begin,
                        std::vector<unsigned int>::iterator end) const
{
  get_set_ids (kind, &(*begin), &(*end));
};
void Mesh::get_set_ids (Entity_kind kind,
                        unsigned int * begin,
                        unsigned int * end) const
{
  Set_ID_List setids;

  get_set_ids (kind, &setids);

  assert (setids.size() <= (end-begin));
  std::copy (setids.begin(), setids.end(), begin);
};


void Mesh::get_set (unsigned int set_id, Entity_kind kind,
                    Parallel_type ptype,
                    std::vector<unsigned int>::iterator begin,
                    std::vector<unsigned int>::iterator end) const
{
  get_set (set_id, kind, ptype, &(*begin), &(*end));
};
void Mesh::get_set (unsigned int set_id, Entity_kind kind,
                    Parallel_type ptype,
                    unsigned int * begin,
                    unsigned int * end) const
{
  Entity_ID_List setents;

  get_set_entities ((Set_ID)set_id, kind, ptype, &setents);

  assert (setents.size() <= (end-begin));
  std::copy (setents.begin(), setents.end(), begin);
};



  // Number of sets containing entities of type 'kind' in mesh
  // 
  // DEPRECATED due to ambiguity in determining what types of sets
  // some regions are supposed to create (a planar region can 
  // result in sidesets or nodesets

unsigned int Mesh::num_sets(const Entity_kind kind) const
{
  int nsets = 0;

  std::cerr << "THIS ROUTINE (num_sets) IS DEPRECATED" << std::endl;
  std::cerr << "It might work but there is no guarantee" << std::endl;
  std::cerr << "that it will work in general situations" << std::endl;
  
  if (!geometric_model_) return 0;

  unsigned int gdim = geometric_model_->dimension();    
  
  unsigned int ngr = geometric_model_->Num_Regions();
  for (int i = 0; i < ngr; i++) {
    AmanziGeometry::RegionPtr rgn = geometric_model_->Region_i(i);
    
    unsigned int rdim = rgn->dimension();
    
    // If we are looking for a cell set the region has to be 
    // of the same topological dimension as the cells
    
    if (kind == CELL && rdim == celldim) 
      {
        nsets++;
      }
    
    // If we are looking for a side set, the region has to be 
    // one topological dimension less than the cells
    
    else if (kind == FACE && rdim == celldim-1) 
      {
        nsets++;
      }
    
    // If we are looking for a node set, the region can be of any
    // dimension upto the cell dimension
    
    else if (kind == NODE && rdim <= celldim)
      {
        nsets++;
      }
  }

  return nsets;
  
}



  // Ids of sets containing entities of 'kind'
  // 
  // DEPRECATED due to ambiguity in determining what types of sets
  // some regions are supposed to create (a planar region can 
  // result in sidesets or nodesets

void Mesh::get_set_ids (const Entity_kind kind, std::vector<Set_ID> *setids) const 
{
  int i, nsets=0;

  std::cerr << "THIS ROUTINE (get_set_ids_) IS DEPRECATED" << std::endl;
  std::cerr << "It might work but there is no guarantee" << std::endl;
  std::cerr << "that it will work in general situations" << std::endl;
  
  assert(setids != NULL);

  setids->clear();

  if (!geometric_model_) return;
  
  int ngr = geometric_model_->Num_Regions();
  for (int i = 0; i < ngr; i++) {
    AmanziGeometry::RegionPtr rgn = geometric_model_->Region_i(i);
    
    unsigned int rdim = rgn->dimension();
    
    // If we are looking for a cell set the region has to be 
    // of the same topological dimension as the cells
    
    if (kind == CELL && rdim == celldim)
      {
        setids->push_back(rgn->id());
      }

      // If we are looking for a side set, the region has to be 
      // one topological dimension less than the cells
      
    else if (kind == FACE && rdim == celldim-1)
      {
        setids->push_back(rgn->id());
      }

      // If we are looking for a node set, the region can be of any
      // dimension upto the cell dimension
      
    else if (kind == NODE && rdim <= celldim) 
      {
        setids->push_back(rgn->id());
      }
    
  }
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
    std::vector<Entity_ID> faces;
    std::vector<unsigned int> nfnodes;
    std::vector<int> fdirs;
    std::vector<AmanziGeometry::Point> cfcoords;
    
    cell_get_faces(cellid,&faces);
    cell_get_face_dirs(cellid,&fdirs);
    
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

} // close namespace AmanziMesh
} // close namespace Amanzi
