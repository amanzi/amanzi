#include "Mesh.hh"
#include <Geometry.hh>

using namespace std;


namespace Amanzi
{

  using namespace AmanziGeometry;

  
  namespace AmanziMesh
  {

    int Mesh::precompute_geometric_quantities() {
      int ncells;
      double volume, area, len;
  
      ncells = num_entities(CELL,USED);
    
      if (celldim == 3) {
	Point centroid(3), normal(3);
	std::vector<Point> ccoords;

	for (int i = 0; i < ncells; i++) {

	  if (cell_get_type(i) == TET) {
 
	    cell_get_coordinates(i,&ccoords);
	    AmanziGeometry::tet_get_vol_centroid(ccoords,&volume,&centroid);

	  }
	  else {
	    // 3D Elements with possibly curved faces
	    // We have to build a description of the element topology
	    // and send it into the polyhedron volume and centroid 
	    // calculation routine

	    int nf;
	    std::vector<Entity_ID> faces;
	    std::vector<unsigned int> nfnodes;
	    std::vector<int> fdirs;
	    std::vector<Point> ccoords, cfcoords, fcoords;

	    cell_get_faces(i,&faces);
	    cell_get_face_dirs(i,&fdirs);

	    nf = faces.size();
	    
	    for (int j = 0; j < nf; j++) {

	      face_get_coordinates(faces[j],&fcoords);
	      nfnodes[j] = fcoords.size();

	      if (fdirs[j] == 1) {
		for (int k = 0; k < nfnodes[j]; k++)
		  cfcoords.push_back(fcoords[k]);
	      }
	      else {
		for (int k = nfnodes[j]-1; k >=0; k--)
		  cfcoords.push_back(fcoords[k]);
	      }
	    }
	    
	    cell_get_coordinates(i,&ccoords);	   

	    AmanziGeometry::polyhed_get_vol_centroid(ccoords,nf,nfnodes,
						     cfcoords,&volume,
						     &centroid);
	  }

	  cell_volumes.push_back(volume);
	  cell_centroids.push_back(centroid);
	}

	int nfaces = num_entities(FACE,USED);
	std::vector<Point> fcoords;

	for (int i = 0; i < nfaces; i++) {
	  face_get_coordinates(i,&fcoords);

	  AmanziGeometry::polygon_get_area_centroid(fcoords,&area,&centroid);
	  face_areas.push_back(area);

	  normal = AmanziGeometry::polygon_get_normal(fcoords);
	  normal *= area;
	  face_normals.push_back(normal);
	}
      }
      else if (celldim == 2) {	
	
	if (spacedim == 2) {   // 2D mesh
	  
	  Point centroid(2);
	  std::vector<Point> ccoords;	  

	  for (int i = 0; i < ncells; i++) {
	    cell_get_coordinates(i,&ccoords);
	    
	    AmanziGeometry::polygon_get_area_centroid(ccoords,&volume,&centroid);
	    cell_volumes.push_back(volume); // yes, area stored as a volume
	    cell_centroids.push_back(centroid);
	  }

	  int nfaces = num_entities(FACE,USED);
	  std::vector<Point> fcoords;
	  Point evec(2), normal(2);

	  for (int i = 0; i < nfaces; i++) {
	    face_get_coordinates(i,&fcoords);

	    evec = fcoords[1]-fcoords[0];
	    area = sqrt(evec*evec);
	    face_areas.push_back(area);
	    
	    normal = Point(evec[1],-evec[0]);
	    face_normals.push_back(normal);
	  }

	}
	else {  // Surface mesh - cells are 2D, coordinates are 3D

	  Point centroid(3);
	  int ncells = num_entities(CELL,USED);
	  std::vector<Point> ccoords;	  
    
	  for (int i = 0; i < ncells; i++) {
	    cell_get_coordinates(i,&ccoords);
	    
	    AmanziGeometry::polygon_get_area_centroid(ccoords,&volume,&centroid);
	    cell_volumes.push_back(volume); // yes, area stored as a volume
	    cell_centroids.push_back(centroid);
	  }

	  // edge normals are ambiguous for surface mesh
	  // So we won't compute them

	}

      }

    } // Mesh::precompute_geometric_quantities()



    // Mesh entity geometry
    //--------------
    //
    
    
    // Volume/Area of cell

    double Mesh::cell_volume (const Entity_ID cellid)
    {
      if (!geometry_precomputed) {
	precompute_geometric_quantities();
	geometry_precomputed = true;
      }

      return cell_volumes[cellid];
    }
    
    // Area/length of face

    double Mesh::face_area(const Entity_ID faceid)
    {
      if (!geometry_precomputed) {
	precompute_geometric_quantities();
	geometry_precomputed = true;
      }

      return face_areas[faceid];
    }
    
    
    // Centroid of cell
    AmanziGeometry::Point Mesh::cell_centroid (const Entity_ID cellid)
    {
      if (!geometry_precomputed) {
	precompute_geometric_quantities();
	geometry_precomputed = true;
      }

      return cell_centroids[cellid];
    }
    
    // Normal to face
    // The vector is normalized and then weighted by the area of the face

    AmanziGeometry::Point Mesh::face_normal (const Entity_ID faceid)
    {
      if (!geometry_precomputed) {
	precompute_geometric_quantities();
	geometry_precomputed = true;
      }

      return face_normals[faceid];
    }
    
  } // end namespace AmanziMesh

}  // end namespace Amanzi  
