/**
 * @file   ColorFunctionRegion.cc
 * @author Rao Garimella
 * @date 
 * 
 * @brief  Implementation of Indicator Region class which derives its
 *         definition from the value of an indicator function in a file
 * 
 * 
 */

#include "color-function-factory.hh"

#include "ColorFunctionRegion.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class ColorFunctionRegion
// -------------------------------------------------------------

// -------------------------------------------------------------
// ColorFunctionRegion:: constructors / destructor
// -------------------------------------------------------------
ColorFunctionRegion::ColorFunctionRegion(const std::string name, 
                                         const unsigned int id,
                                         const std::string file,
                                         const int value,
                                         const Epetra_MpiComm *comm)
  : Region(name,id,3),
    file_(file), value_(value)
{
  // Region dimension is set arbitrarily as 3 since the set of
  // entities in the mesh will determine the dimension

  ColorFunctionFactory colfunc_factory;
  try {
    colorfunc_ = colfunc_factory.Create(file_,*comm);
  }
  catch (Errors::Message &msg) {
    Exceptions::amanzi_throw(msg);
  }

}

ColorFunctionRegion::ColorFunctionRegion(const char *name, 
                                         const unsigned int id,
                                         const char *file,
                                         const int value,
                                         const Epetra_MpiComm *comm)
  : Region(name,id,3),
    file_(file), value_(value)
{
  // Region dimension is set arbitrarily as 3 since the set of
  // entities in the mesh will determine the dimension

  ColorFunctionFactory colfunc_factory;
  try {
    colorfunc_ = colfunc_factory.Create(file_,*comm);
  }
  catch (Errors::Message &msg) {
    Exceptions::amanzi_throw(msg);
  }
}

ColorFunctionRegion::ColorFunctionRegion(const ColorFunctionRegion& old)
  : Region(old),file_(old.file_),value_(old.value_)
{
}

ColorFunctionRegion::~ColorFunctionRegion(void)
{
  delete colorfunc_;
}


// -------------------------------------------------------------
// ColorFunctionRegion::inside
// -------------------------------------------------------------
bool
ColorFunctionRegion::inside(const Point& p) const
{
  try {
    int color = (*colorfunc_)(&(p[0]));
    return (color == value_);
  }
  catch (Errors::Message &msg) {
    std::cerr << msg.what() << std::endl;
    return false;
  }
}
  

} // namespace AmanziGeometry
} // namespace Amanzi
