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

  Teuchos::Array<int> ncells = parameter_list.get< Teuchos::Array<int> >("Number of Cells");
  Teuchos::Array<double> low_corner = parameter_list.get< Teuchos::Array<double> >("Domain Low Corner");
  Teuchos::Array<double> high_corner = parameter_list.get< Teuchos::Array<double> >("Domain High Corner");

  unsigned int dimension = ncells.size();

  nx_ = ncells[0];
  ny_ = ncells[1];
  if (dimension == 3) nz_ = ncells[2];

  AmanziGeometry::Point p0(dimension), p1(dimension);

  p0.set(&(low_corner[0]));
  p1.set(&(high_corner[0]));


  domain_ = new AmanziGeometry::BoxRegion("GenDomain", 0, p0, p1);

}

} // end namespace AmanziMesh
} // end namespace Amanzi
