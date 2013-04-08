/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   ColorFunctionRegion.hh
 * @author Rao Garimella
 * @date 
 * 
 * @brief  Declaration of ColorFunction Region class which derives 
 *         its definition from specified value of an color function
 *         i.e. the region exists wherever the value of the color
 *         function matches a specified value
 * 
 * 
 */

#ifndef _ColorFunctionRegion_hh_
#define _ColorFunctionRegion_hh_

#include "Epetra_MpiComm.h"

#include "color-function.hh"

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class ColorFunctionRegion
// -------------------------------------------------------------
/// A region defined by the value of an indicator function in a file
///
/// The region will consist of all mesh elements for which the indicator
/// function is a particular value at their centroids

class ColorFunctionRegion : public Region {
public:

  /// Default constructor 

  ColorFunctionRegion(const std::string name, 
                      const unsigned int id, 
                      const std::string file,
                      const int value,
                      const Epetra_MpiComm *comm);

  ColorFunctionRegion(const char *name, 
                      const unsigned int id, 
                      const char *file,
                      const int value,
                      const Epetra_MpiComm *comm);


  /// Protected copy constructor to avoid unwanted copies.
  ColorFunctionRegion(const ColorFunctionRegion& old);

  /// Destructor
  ~ColorFunctionRegion(void);

  // Type of the region
  inline RegionType type() const { return COLORFUNCTION; }

  /// Is the the specified point inside this region
  bool inside(const Point& p) const;

protected:  
  std::string file_; // which file are we supposed to read it from
  const int value_;
  const ColorFunction *colorfunc_; // indicator func created from file
};

/// A smart pointer to ColorFunctionRegion instances
// typedef Teuchos::RCP<ColorFunctionRegion> ColorFunctionRegionPtr;

// RVG: I am not able to correctly code a region factory using smart
// pointers so I will revert to a simpler definition

typedef ColorFunctionRegion *ColorFunctionRegionPtr;

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
