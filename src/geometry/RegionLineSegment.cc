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
    const std::string& name, const Set_ID id,
    const Point& p0, const Point& p1,
    const LifeCycleType lifecycle)
  : Region(name, id, true, LINE_SEGMENT, p0.dim(), p0.dim(), lifecycle),
    p0_(p0),
    p1_(p1)
{

  line_points_.clear();
  line_frac_.clear();
  
  Errors::Message msg;
  if (p0_.dim() != p1_.dim()) {
    msg << "Mismatch in dimensions of end points of RegionLineSegment \""
        << Region::name() << "\"";
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

  mdim = manifold_dimension();
  sdim = polytope[0].dim();

  if ((mdim == 3)&&(sdim==3)){
    std::vector<Point> line(2);
    bool intersct = false;
    line[0] = p0_;
    line[1] = p1_;
    std::vector<Point> inter_pnts(2);
    Point v1(p0_.dim());    
    int num_int = 0;
    for (int i=0; i<faces.size(); i++){
      intersct = false;
      std::vector<Point> plane(faces[i].size());
      for (int j=0;j<plane.size();j++) plane[j] = polytope[faces[i][j]];
      double t = PlaneLineIntersection(plane, line);
      if ((t<0)||(t>1)) continue;
      std::cout<<"t= "<<t<<"\n";    

      v1 = p0_ + t*(p1_ - p0_);
      //std::cout << "intersection "<<v1<<"\n";
      intersct = point_in_polygon(v1, plane);
      if (intersct){
        inter_pnts[num_int] = v1;
        num_int++;
      }
    }
    if (num_int == 0) return 0.;
    else if (num_int == 1){
      double t = norm(inter_pnts[0] - p0_)/norm(p1_ - p0_);
      double len;
      if (t>=0.5){
        v1 = 0.5*(p1_ + inter_pnts[0]);
        len = norm(p1_ - inter_pnts[0]);
      }else{
        v1 = 0.5*(p0_ + inter_pnts[0]);
        len = norm(p0_ - inter_pnts[0]);
      }
      // line_points_.push_back(v1);
      // line_frac_.push_back(len);
      return 1;
    }else if (num_int == 2){
      v1 =0.5*(inter_pnts[1] + inter_pnts[0]); 
      // line_points_.push_back(v1);
      // line_frac_.push_back(norm(inter_pnts[1] - inter_pnts[0]));
      return 1;
    }
    else {
      Errors::Message mesg("More than two intersection points in RegionLineSegment intersect function"); 
      Exceptions::amanzi_throw(mesg);
    }
  }



  return 0.;

}

  /* 
     Compute intersection of line and plane which defined 
     by 3 points of a face.
  */


double PlaneLineIntersection(const std::vector<Point>& plane,
                             const std::vector<Point>& line){

  std::vector<double> row1(4, 1.), row2(4,1);
  row2[3] = 0.;

  std::vector<double> smatr1(12), smatr2(12);

  int k=0;
  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      smatr1[k] = plane[i][j];
      smatr2[k] = plane[i][j];
      k++;
    }
  }

  for (int i=0;i<3;i++){
    smatr1[k] = line[0][i];
    smatr2[k] = line[1][i] - line[0][i];
    k++;
  }

  
  double t1 = det_aux(row1, smatr1);
  double t2 = det_aux(row2, smatr2);

    //    std::cout<< "PlaneLineIntersection: either plane or line is degenerate\n";
  // std::cout.precision(16);
  // std::cout<<"line "<<line[0]<<" : "<<line[1]<<"\n";
  // std::cout<<t1<<" "<<t2<<"\n\n";
  // for(int i=0; i<4; i++) std::cout<<row1[i]<<" "; std::cout<<"\n";
  // for(int i=0; i<3; i++){
  //   for (int j=0;j<4;j++) std::cout<<smatr1[i + j*3]<<" ";std::cout<<"\n";
  // }
  // std::cout<<"\n";
  // for(int i=0; i<4; i++) std::cout<<row2[i]<<" "; std::cout<<"\n";
  // for(int i=0; i<3; i++){
  //   for (int j=0;j<4;j++) std::cout<<smatr2[i + j*3]<<" ";std::cout<<"\n";
  // }



  if (abs(t2) < 1e-10){
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
               const std::vector<double>& submatr){


  if ((first_row.size()<4)||(submatr.size()<12)){
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
