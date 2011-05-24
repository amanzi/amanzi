// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: stk_mesh_error.hh
// -------------------------------------------------------------
/**
 * @file   stk_mesh_error.hh
 * @author William A. Perkins
 * @date Mon May  2 14:29:37 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// Created April 21, 2011 by William A. Perkins
// Last Change: Mon May  2 14:29:37 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _stk_mesh_error_hh_
#define _stk_mesh_error_hh_

#include "errors.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace STK {

class Error : public Errors::Message 
{
 public:
  explicit Error(void) : Errors::Message() {};
  explicit Error(const char* message) : Errors::Message(message) {};
  ~Error(void) throw() {};
};

} // close namespace STK 
} // close namespace Mesh 
} // close namespace Amanzi 
#endif
