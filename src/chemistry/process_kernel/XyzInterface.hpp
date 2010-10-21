/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef __XyzInterface_hpp__
#define __XyzInterface_hpp__

/*
** Description:  class for interfacing XYZ to the ascem/amanzi geochemistry code
*/

#include <string>
#include <vector>

#include "Geochemistry.hpp"

class XyzInterface : Geochemistry {
 public:
  virtual ~XyzInterface();

  virtual void Advance();
  virtual void CommitState();

 protected:
  XyzInterface();

 private:
};

#endif  // __XyzInterface_hpp__
