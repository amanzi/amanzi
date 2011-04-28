#include "Geometry.hh"

using namespace std;


namespace Amanzi
{

  namespace AmanziGeometry
  {

    // Return the volume and centroid of a general polyhedron
    //
    // ccoords  - vertices of the polyhedron (in no particular order)
    // nf       - number of faces of polyhedron
    // nfnodes  - number of nodes for each face
    // fcoords  - linear array of face coordinates in in ccw manner
    //            assuming normal of face is pointing out (
    //
    // So if the polyhedron has 5 faces with 5,3,3,3 and 3 nodes each
    // then entries 1-5 of fcoords describes face 1, entries 6-9
    // describes face 2 and so on
    //
    // So much common work has to be done for computing the centroid
    // and volume calculations that they have been combined into one
    //
    // The volume of all polyhedra except tets is computed as a sum of
    // volumes of tets created by connecting the polyhedron center to
    // a face center and an edge of the face      


    void polyhed_get_vol_centroid(const std::vector<Point> ccoords, 
				  const unsigned int nf, 
				  const std::vector<unsigned int> nfnodes,
				  const std::vector<Point> fcoords,
				  double *volume,
				  Point *centroid)
    {

      // Initialize to sane values

      centroid->set(0.0);
      (*volume) = 0.0;
	
      // Compute the geometric center of all face nodes

      Point center(3);

      int np = ccoords.size();
      if (np < 4) {
	cout << "Not a polyhedron" << std::endl;
	return;
      }

      if (np == 4) { // is a tetrahedron

	tet_get_vol_centroid(ccoords,volume,centroid);

      }
      else { // if (np > 4), polyhedron with possibly curved faces
	
	for (int i = 0; i < np; i++)
	  center += ccoords[i];
	center /= np;

	int offset = 0;
	for (int i = 0; i < nf; i++) {

	  // geometric center of all face nodes
	  
	  Point fcenter(3);
	  for (int j = 0; j < nfnodes[i]; j++)
	    fcenter += fcoords[offset+j];
	  fcenter /= nfnodes[i];
	  
	  for (int j = 0; j < nfnodes[i]; j++) { // for each edge of face
	    
	    // form tet from edge of face, face center and cell center
	    
	    Point tcentroid(3);
	    double tvolume;
	    int k, kp1;
	    Point v1(3), v2(3), v3(3);
	    
	    k = offset+j;
	    kp1 = offset+(j+1)%nfnodes[i];
	    
	    tcentroid = (center+fcenter+fcoords[k]+fcoords[kp1])/4.0;
	    v1 = fcoords[kp1]-center;
	    v2 = fcoords[k]-center;
	    v3 = fcenter-center;
	    tvolume = (v1^v2)*v3;
	    
	    (*centroid) += tvolume*tcentroid;      // sum up 1st moment
	    (*volume) += tvolume;                  // sum up 0th moment
	    
	  } // for each edge of face

	  offset += nfnodes[i];

	} // for each face


	(*centroid) /= (*volume);  // centroid = 1st moment / 0th moment
      } // end if (np > 4)
 
    } // polyhed_get_vol_centroid





    // Special function for tet volume and centroid

    void tet_get_vol_centroid(const std::vector<Point> ccoords, 
			      double *volume,
			      Point *centroid)
    {
      // Initialize to sane values
      
      centroid->set(0.0);
      (*volume) = 0.0;
      
      int np = ccoords.size();
      if (np < 4) {
	cout << "Less than 4 points supplied for tetrahedron" << std::endl;
	return;
      }
      else if (np > 4) 
	cout << "More than 4 points supplied for tetrahedron" << std::endl;

      (*centroid) = (ccoords[0]+ccoords[1]+ccoords[2]+ccoords[3])/4.0;

      Point v1 = ccoords[1]-ccoords[0];
      Point v2 = ccoords[2]-ccoords[0];
      Point v3 = ccoords[3]-ccoords[0];
      (*volume) = ((v1^v2)*v3)/6.0;
 
    } // tet_get_vol_centroid






    // Compute area of polygon 

    // In 2D, the area and centroid are computed by a contour integral
    // around the perimeter.
    //
    // In 3D, the area and centroid computed by connecting a "center"
    // point of the polygon to the edges of the polygon and summing
    // the moments of the resulting triangles

    void polygon_get_area_centroid(const std::vector<Point> coords,
				   double *area, Point *centroid) {


      (*area) = 0;
	centroid->set(0.0);

      unsigned int np = coords.size();

      if (np < 3) {
	cout << "Degenerate polygon - area is zero" << std::endl;
	return;
      }


      int dim = coords[0].dim();

      switch (dim) {
      case 2: {
	// We may have to adjust coordinates relative to one
	// corner of the polygon if the coordinates are very
	// large and are generating large roundoff errors.

	double moment_1x=0, moment_1y=0;

	for (int i = 0; i < np; i++) {
	  (*area) += ((coords[i].x())*(coords[(i+1)%np].y()) - 
		      (coords[(i+1)%np].x())*(coords[i].y()));

	  moment_1x += ( ((coords[i].x())*(coords[(i+1)%np].y()) - 
			  (coords[(i+1)%np].x())*(coords[i].y())) *
			 (coords[i].x()+coords[i+1].x()) );

	  moment_1y += ( ((coords[i].x())*(coords[(i+1)%np].y()) - 
			  (coords[(i+1)%np].x())*(coords[i].y())) *
			 (coords[i].y()+coords[i+1].y()) );
	}
	    
	(*area) /= 2;

	moment_1x /= 6;
	moment_1y /= 6;
       
	(*centroid).set( moment_1x/(*area), moment_1y/(*area) );

	break;
      }
      case 3: {

	Point center(3);
	
	// Compute a center point 
	
	for (int i = 0; i < np; i++)
	  center += coords[i];
	center /= np;
	  
	if (coords.size() == 3) { // triangle - straightforward
	  Point v1 = coords[2]-coords[1];
	  Point v2 = coords[0]-coords[1];

	  Point v3 = v1^v2;

	  *area = sqrt(v1*v1);
	  (*centroid) = center;   
	}
	else {
	  // Compute the area of each triangle formed by
	  // the center point and each polygon edge
	  
	  *area = 0.0;
	  for (int i = 0; i < np; i++) {
	    Point v1 = coords[i]-center;
	    Point v2 = coords[(i+1)%np]-center;
	    
	    Point v3 = v1^v2;
	    
	    (*area) += 0.5*sqrt(v3*v3);
	    (*centroid) += (coords[i]+coords[(i+1)%np]+center)/3.0;
	  }

	  (*centroid) /= (*area);
	}
	break;
      }
      default:
	cout << "Invalid coordinate dimensions" << std::endl;
      } // end switch(dim)
      
    } // polygon_get_area_centroid
  




    // Get normalized normal of polygon
    // In 2D, the normal is unambiguous - the normal is evaluated at one corner
    // In 3D, the procedure evaluates the normal at each corner and averages it

    Point polygon_get_normal(const std::vector<Point> coords)
    {
      Point normal;

      if (coords.size() < 3) {
	cout << "Degenerate polygon - Cannot compute normal" << std::endl;	
	Point dummynormal(0,0,0);
	return dummynormal;
      }

      int dim = coords[0].dim();

      switch (dim) {
      case 2: {
	// sufficient to evaluate normal at one corner

	Point v1 = coords[2]-coords[1];
	Point v2 = coords[0]-coords[1];

	normal = v1^v2;
	break;
      }
      case 3: {
	// evaluate normal at all corners and average

	unsigned int np = coords.size();
	for (int i = 0; i < np; i++) {
	  Point v1 = coords[(i+1)%np]-coords[i];
	  Point v2 = coords[(i-1+np)%np]-coords[i];
	  normal += v1^v2;
	}
	break;
      }
      default:
	cout << "Invalid coordinate dimensions" << std::endl;
      }

      normal /= norm(normal);

      return normal;

    } // polygon_get_normal



  } // end namespace AmanziGeometry

} // end namespace Amanzi


