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
	Point centroid(spacedim), normal(spacedim);
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



    // Temporary routines for backward compatibility
    //
    // ----- ATTENTION - ATTENTION - ATTENTION -----
    // NEW CODE SHOULD NOT USE THESE FUNCTIONS SINCE THEY WILL GO AWAY
    //
    


    void Mesh::cell_to_faces (unsigned int cell, 
			      std::vector<unsigned int>::iterator begin, 
			      std::vector<unsigned int>::iterator end) 
    {
      cell_to_faces (cell, &(*begin), &(*end));
    };
    
    void Mesh::cell_to_faces (unsigned int cell, 
			      unsigned int* begin, unsigned int *end) 
    {
      Entity_ID_List cfaces;
      
      cell_get_faces ((Entity_ID)cell, &cfaces);

      assert (cfaces.size() <= (unsigned int) (end-begin));
      std::copy (cfaces.begin(), cfaces.end(), begin);
    };


    void Mesh::cell_to_face_dirs (unsigned int cell, 
				  std::vector<int>::iterator begin, 
			    std::vector<int>::iterator end) 
    {
      cell_to_face_dirs (cell, &(*begin), &(*end));
    };
    void Mesh::cell_to_face_dirs (unsigned int cell, 
				  int * begin, int * end) 
    {
      vector<int> cfdirs;

      cell_get_face_dirs ((Entity_ID)cell, &cfdirs);

      assert (cfdirs.size() <= (unsigned int) (end-begin));
      std::copy (cfdirs.begin(), cfdirs.end(), begin);
    };
  
    

    void Mesh::cell_to_nodes (unsigned int cell, 
			      std::vector<unsigned int>::iterator begin, 
			      std::vector<unsigned int>::iterator end) 
    {
      cell_to_nodes (cell, &(*begin), &(*end));
    };
    void Mesh::cell_to_nodes (unsigned int cell, 
			      unsigned int * begin, unsigned int * end) 
    {
      Entity_ID_List cnodes;

      cell_get_nodes ((Entity_ID)cell, &cnodes);

      assert (cnodes.size() <= (unsigned int) (end-begin));
      std::copy (cnodes.begin(), cnodes.end(), begin);
    };
      



    void Mesh::face_to_nodes (unsigned int face, 
			      std::vector<unsigned int>::iterator begin, 
			      std::vector<unsigned int>::iterator end) 
    {
      face_to_nodes(face, &(*begin), &(*end));
    };
    void Mesh::face_to_nodes (unsigned int face, 
			      unsigned int * begin, unsigned int * end) 
    {
      Entity_ID_List fnodes;

      face_get_nodes ((Entity_ID)face, &fnodes);

      assert (fnodes.size() <= (unsigned int) (end-begin));
      std::copy (fnodes.begin(), fnodes.end(), begin);
    };



    void Mesh::node_to_coordinates (unsigned int node, 
				    std::vector<double>::iterator begin, 
				    std::vector<double>::iterator end) 
    {
      node_to_coordinates (node, &(*begin), &(*end));
    };
    void Mesh::node_to_coordinates (unsigned int node, 
				    double * begin, 
				    double * end) 
    {
      Point xyz;
      
      node_get_coordinates ((Entity_ID)node, &xyz);
      
      assert (xyz.dim() == (end-begin));
      for (int i = 0; i < xyz.dim(); i++) begin[i] = xyz[i];
    };

    void Mesh::face_to_coordinates (unsigned int face, 
				    std::vector<double>::iterator begin, 
				    std::vector<double>::iterator end)
    {
      face_to_coordinates (face, &(*begin), &(*end));
    };
    void Mesh::face_to_coordinates (unsigned int face, 
				    double * begin, 
				    double * end)
    {
      vector<Point> fxyz;

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
				    std::vector<double>::iterator end) 
    {
      cell_to_coordinates (cell, &(*begin), &(*end));
    };
    void Mesh::cell_to_coordinates (unsigned int cell, 
				    double * begin,
				    double * end) 
    {
      vector<Point> cxyz;
      
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


    // Unchanged in new interface
    // unsigned int num_sets(Entity_kind kind) const {};
      
    // Unchanged in new interface
    // unsigned int get_set_size (unsigned int set_id, 
    //			       Entity_kind kind,
    //			       Parallel_type ptype) const {};

    // Id numbers
    void Mesh::get_set_ids (Entity_kind kind, 
			    std::vector<unsigned int>::iterator begin, 
			    std::vector<unsigned int>::iterator end)
    {
      get_set_ids (kind, &(*begin), &(*end));
    };
    void Mesh::get_set_ids (Entity_kind kind, 
			    unsigned int * begin, 
			    unsigned int * end)
    {
      Set_ID_List setids;

      get_set_ids (kind, &setids);

      assert (setids.size() <= (end-begin));
      std::copy (setids.begin(), setids.end(), begin);
    };

    // Unchanged in new interface
    // bool valid_set_id (unsigned int id, Entity_kind kind) const {};
      

    void Mesh::get_set (unsigned int set_id, Entity_kind kind, 
			Parallel_type ptype, 
			std::vector<unsigned int>::iterator begin, 
			std::vector<unsigned int>::iterator end)
    {
      get_set (set_id, kind, ptype, &(*begin), &(*end));
    };
    void Mesh::get_set (unsigned int set_id, Entity_kind kind, 
			Parallel_type ptype, 
			unsigned int * begin, 
			unsigned int * end)
    {
      Entity_ID_List setents;

      get_set_entities ((Set_ID)set_id, kind, ptype, &setents);

      assert (setents.size() <= (end-begin));
      std::copy (setents.begin(), setents.end(), begin);
    };
    
  } // end namespace AmanziMesh

}  // end namespace Amanzi  
