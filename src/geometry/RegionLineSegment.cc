/*
  Geometry

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)

  A region defined by a line segment.
*/

#include <vector>

#include "dbc.hh"
#include "errors.hh"

#include "Point.hh"
#include "Geometry.hh"
#include "RegionLineSegment.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------------
// Constructor
// -------------------------------------------------------------------
RegionLineSegment::RegionLineSegment(
    const std::string& name, const int id,
    const Point& p0, const Point& p1,
    const LifeCycleType lifecycle)
  : Region(name, id, true, RegionType::LINE_SEGMENT, p0.dim(), p0.dim(), lifecycle),
    p0_(p0),
    p1_(p1)
{
  set_manifold_dimension(3);
  //line_points_.clear();
  //line_frac_.clear();
  
  Errors::Message msg;
  if (p0_.dim() != p1_.dim()) {
    msg << "Mismatch in dimensions of end points of RegionLineSegment \""
        << Region::get_name() << "\"";
    Exceptions::amanzi_throw(msg);
  }

  double eps = 1e-15;
  if (norm(p0_ - p1_) < eps) {
    msg <<" Zero length line segment \""<< Region::get_name() <<"\" is NOT allowed."; 
    Exceptions::amanzi_throw(msg);
  }
}


// -------------------------------------------------------------------
// Implementation of a virtual member function.
// -------------------------------------------------------------------
bool RegionLineSegment::inside(const Point& p) const
{
  Errors::Message mesg("In/out check not implemented for line segment sets");
  Exceptions::amanzi_throw(mesg);
  return false;
}


// -------------------------------------------------------------------
// Implementation of a virtual member function.
// We have to analyze 
// -------------------------------------------------------------------
double RegionLineSegment::intersect(
    const std::vector<Point>& polytope,
    const std::vector<std::vector<int> >& faces) const
{
  int mdim, sdim;

  mdim = get_manifold_dimension();
  sdim = polytope[0].dim();

  if ((mdim == 3)&&(sdim==3)) {
    std::vector<Point> line(2);
    bool intersct = false;
    double eps = 1e-12;
    line[0] = p0_;
    line[1] = p1_;
    std::vector<Point> inter_pnts(2);
    Point v1(p0_.dim()), vp(p0_.dim());  
    double tau[2]={0., 0.};
    int num_int = 0;
    int num_line_int = 0;
    for (int i=0; i<faces.size(); i++) {
      intersct = false;
      std::vector<Point> plane(faces[i].size());
      for (int j=0;j<plane.size();j++) plane[j] = polytope[faces[i][j]];
      double t = PlaneLineIntersection(plane, line);

      v1 = p0_ + t*(p1_ - p0_);
      double diff_x = 0., diff_y = 0.;
      for (int j=0; j<plane.size(); j++) {
        diff_x += fabs(plane[j].x() - v1.x());
        diff_y += fabs(plane[j].y() - v1.y());
      }
      if (diff_x < eps) { // projection to plane y-z
        vp[0] = v1[1]; vp[1] = v1[2];
        for (int j=0; j<plane.size(); j++) {
          plane[j][0] = plane[j][1];
          plane[j][1] = plane[j][2];
        }
      } else if (diff_y < eps) { //projection to plane x-z
        vp[0] = v1[0]; vp[1] = v1[2];
        for (int j=0; j<plane.size(); j++) {         
          plane[j][1] = plane[j][2];
        }
      } else {
        vp = v1;
      }
          
      intersct = point_in_polygon(vp, plane);

      if (intersct) {
        if (std::fabs(t) < eps) t = 0;
        tau[num_line_int] = t;
        num_line_int++;
        if ((t>=0)&&(t<=1)) {
          inter_pnts[num_int] = v1;
          num_int++;
        }
      }
    }
    double len;
    if (num_int == 0) {
      if (num_line_int==2) {
        if (tau[0]*tau[1]<0) {
          len = norm(p1_ - p0_);
          return len;
        }       
      }
      return 0.;
    }
    else if (num_int == 1) {
      len = 0.;
      if ((tau[0]<0)||(tau[1]<0)) {
        // v1 = 0.5*(p0_ + inter_pnts[0]);
        len = norm(p0_ - inter_pnts[0]);        
      } else if ((tau[0]>=0) && (tau[1]>=0)) {
        //v1 = 0.5*(p1_ + inter_pnts[0]);
        len = norm(p1_ - inter_pnts[0]);
      }
      return len;
    } else if (num_int == 2) {
      len = norm(inter_pnts[1] - inter_pnts[0]);
      return len;
    } else {
      Errors::Message mesg("More than two intersection points in RegionLineSegment intersect function"); 
      Exceptions::amanzi_throw(mesg);
    }
  }

  return 0.;
}


void RegionLineSegment::ComputeInterLinePoints(const std::vector<Point>& polytope,
                                               const std::vector<std::vector<int> >& faces,
                                               Point& res_point) const
{
  int mdim, sdim;
  double eps = 1e-12;

  mdim = get_manifold_dimension();
  sdim = polytope[0].dim();

  if ((mdim == 3)&&(sdim==3)) {
    std::vector<Point> line(2);
    bool intersct = false;
    line[0] = p0_;
    line[1] = p1_;
    std::vector<Point> inter_pnts(2);
    Point v1(p0_.dim()), vp(p0_.dim());    
    int num_int = 0;
    double tau[2]={0., 0.};
    int num_line_int = 0;    


    for (int i=0; i<faces.size(); i++) {
      intersct = false;
      std::vector<Point> plane(faces[i].size());
      for (int j=0;j<plane.size();j++) plane[j] = polytope[faces[i][j]];
      double t = PlaneLineIntersection(plane, line);
      //if ((t<0)||(t>1)) continue;

      v1 = p0_ + t*(p1_ - p0_);
      double diff_x = 0., diff_y = 0.;
      for (int j=0; j<plane.size(); j++) {
        diff_x += fabs(plane[j].x() - v1.x());
        diff_y += fabs(plane[j].y() - v1.y());
      }
      if (diff_x < eps) { // projection to plane y-z
        vp[0] = v1[1]; vp[1] = v1[2];
        for (int j=0; j<plane.size(); j++) {
          plane[j][0] = plane[j][1];
          plane[j][1] = plane[j][2];
        }
      } else if (diff_y < eps) { //projection to plane x-z
        vp[0] = v1[0]; vp[1] = v1[2];
        for (int j=0; j<plane.size(); j++) {         
          plane[j][1] = plane[j][2];
        }
      } else {
        vp = v1;
      }

      intersct = point_in_polygon(vp, plane);

      if (intersct) {
        if (std::fabs(t) < eps) t = 0;
        tau[num_line_int] = t;
        num_line_int++;
        if ((t>=0)&&(t<=1)) {
          inter_pnts[num_int] = v1;
          num_int++;
        }
      }
    }

    if (num_int == 0) {
      res_point = 0.5*(p1_ + p0_);
      return;
      //Errors::Message mesg("No intersection points in RegionLineSegment"); 
      //Exceptions::amanzi_throw(mesg);
    }
    else if (num_int == 1) {
      if ((tau[0]<0)||(tau[1]<0)) {
        res_point = 0.5*(p0_ + inter_pnts[0]);
      } else if ((tau[0]>=0) && (tau[1]>=0)) {
        res_point = 0.5*(p1_ + inter_pnts[0]);
      }
      return;
    } else if (num_int == 2) {
      res_point =0.5*(inter_pnts[1] + inter_pnts[0]); 
      return;
    }
    else {
      Errors::Message mesg("More than two intersection points in RegionLineSegment intersect function"); 
      Exceptions::amanzi_throw(mesg);
    }
  }
}


/* 
   Compute intersection of line and plane which defined 
   by 3 points of a face.
*/
double PlaneLineIntersection(const std::vector<Point>& plane,
                             const std::vector<Point>& line)
{
  std::vector<double> row1(4, 1.), row2(4,1);
  row2[3] = 0.;

  std::vector<double> smatr1(12), smatr2(12);

  int k=0;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      smatr1[k] = plane[i][j];
      smatr2[k] = plane[i][j];
      k++;
    }
  }

  for (int i=0;i<3;i++) {
    smatr1[k] = line[0][i];
    smatr2[k] = line[1][i] - line[0][i];
    k++;
  }

  
  double t1 = det_aux(row1, smatr1);
  double t2 = det_aux(row2, smatr2);

  // std::cout<< "PlaneLineIntersection: either plane or line is degenerate\n";
  // std::cout.precision(16);
  // std::cout<<"line "<<line[0]<<" : "<<line[1]<<"\n";
  // std::cout<<t1<<" "<<t2<<"\n\n";
  // for (int i=0; i<4; i++) std::cout<<row1[i]<<" "; std::cout<<"\n";
  // for (int i=0; i<3; i++) {
  //   for (int j=0;j<4;j++) std::cout<<smatr1[i + j*3]<<" ";std::cout<<"\n";
  // }
  // std::cout<<"\n";
  // for (int i=0; i<4; i++) std::cout<<row2[i]<<" "; std::cout<<"\n";
  // for (int i=0; i<3; i++) {
  //   for (int j=0;j<4;j++) std::cout<<smatr2[i + j*3]<<" ";std::cout<<"\n";
  // }

  if (fabs(t2) < 1e-10) {
    return nan("");
  }

  return -t1/t2;
}


/*
    Computes determinant of matrix 4x4.
    x x x x      - first_row
    x x x x
    x x x x      - submatr
    x x x x
*/
double det_aux(const std::vector<double>& first_row,
               const std::vector<double>& submatr)
{
  if ((first_row.size()<4)||(submatr.size()<12)) {
    Errors::Message mesg("Incorrect parameters in det_aux"); 
    Exceptions::amanzi_throw(mesg);
  }

  double res = 0.;

  res += first_row[0]*(submatr[3]*submatr[7]*submatr[11] +
                       submatr[6]*submatr[10]*submatr[5] +
                       submatr[4]*submatr[8]*submatr[9] -
                       submatr[5]*submatr[7]*submatr[9] -
                       submatr[4]*submatr[6]*submatr[11] -
                       submatr[3]*submatr[8]*submatr[10]);

  res -= first_row[1]*(submatr[0]*submatr[7]*submatr[11] +
                       submatr[6]*submatr[10]*submatr[2] +
                       submatr[1]*submatr[8]*submatr[9] -
                       submatr[2]*submatr[7]*submatr[9] -
                       submatr[1]*submatr[6]*submatr[11] -
                       submatr[0]*submatr[8]*submatr[10]);

  res += first_row[2]*(submatr[0]*submatr[4]*submatr[11] +
                       submatr[3]*submatr[10]*submatr[2] +
                       submatr[1]*submatr[5]*submatr[9] -
                       submatr[2]*submatr[4]*submatr[9] -
                       submatr[1]*submatr[3]*submatr[11] -
                       submatr[0]*submatr[5]*submatr[10]);

  res -= first_row[3]*(submatr[0]*submatr[4]*submatr[8] +
                       submatr[3]*submatr[7]*submatr[2] +
                       submatr[1]*submatr[5]*submatr[6] -
                       submatr[2]*submatr[4]*submatr[6] -
                       submatr[1]*submatr[3]*submatr[8] -
                       submatr[0]*submatr[5]*submatr[7]);

  return res;
}

} // namespace AmanziGeometry
} // namespace Amanzi
