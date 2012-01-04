/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   GenerationSpec.cc
 * @author William A. Perkins
 * @date Wed Sep 28 10:20:33 2011
 * 
 * @brief  
 * 
 * 
 */


#include <Teuchos_ParameterList.hpp>

#include "GenerationSpec.hh"

namespace Amanzi {
namespace AmanziMesh {

// -------------------------------------------------------------
//  class GenerationSpec
// -------------------------------------------------------------

// -------------------------------------------------------------
// GenerationSpec:: constructors / destructor
// -------------------------------------------------------------
GenerationSpec::GenerationSpec(const Teuchos::ParameterList& parameter_list)
  : domain_(NULL), nx_(0), ny_(0), nz_(0), blocks_()
{
  parse_(parameter_list);
}

GenerationSpec::~GenerationSpec(void)
{
  // empty
}

// -------------------------------------------------------------
// GenerationSpec::parse_
// -------------------------------------------------------------
void
GenerationSpec::parse_(const Teuchos::ParameterList& parameter_list)
{
  double x0, y0, z0;
  double x1, y1, z1;

  // read the parameters from the parameter list

  nx_ = parameter_list.get<int>("Number of Cells in X");
  ny_ = parameter_list.get<int>("Number of Cells in Y");
  nz_ = parameter_list.get<int>("Number of Cells in Z");
  
  x0 = parameter_list.get<double>("X_Min");
  x1 = parameter_list.get<double>("X_Max");

  y0 = parameter_list.get<double>("Y_Min");
  y1 = parameter_list.get<double>("Y_Max");

  z0 = parameter_list.get<double>("Z_Min");
  z1 = parameter_list.get<double>("Z_Max");


  AmanziGeometry::Point p0(x0, y0, z0);
  AmanziGeometry::Point p1(x1, y1, z1);
  domain_ = new AmanziGeometry::BoxRegion("GenDomain", 0, p0, p1);
  

  // This part is encapsulated in geometric model
  // The mesh specific mesh generation procedures will
  // read the geometric model and create the necessary regions

  //  int nblk(0);
  //  try {
  //    nblk = parameter_list.get<int>("Number of mesh blocks");
  //  } catch (const Teuchos::Exceptions::InvalidParameterName& e) {
  //    // this is OK, just eat the exception
  // }
  
  // if (nblk > 0) {
  //   for (int nb = 1; nb <= nblk; nb++) {
  //     std::stringstream s; 
  //     s << "Mesh block " << nb;

  //     Teuchos::ParameterList sublist = parameter_list.sublist(s.str());

  //     // tell the generator about the zone

  //     AmanziGeometry::Point p0(x0, y0, sublist.get<double>("Z0"));
  //     AmanziGeometry::Point p1(x1, y1, sublist.get<double>("Z1"));
  //     AmanziGeometry::RegionPtr r(new AmanziGeometry::RectangularRegion(p0, p1));

  //     blocks_.push_back(r);
  //   }
  // }
}

} // end namespace AmanziMesh
} // end namespace Amanzi
