/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for liquid water. See the permafrost physical properties notes for
  references and documentation of this EOS at:
*/

#include "EOS_Water.hh"

namespace Amanzi {
namespace EOS {

// registry of method
Utils::RegisteredFactory<EOS, EOS_Water> EOS_Water::factory_("liquid water");

}  // namespace EOS
}  // namespace Amanzi
