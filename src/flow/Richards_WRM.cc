/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Neil Carlson (nnc@lanl.gov), 
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "Richards_PK.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Convertion p -> s.                                               
****************************************************************** */
void Richards_PK::DeriveSaturationFromPressure(const Epetra_Vector& p, Epetra_Vector& s)
{
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm->WRM();  

  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[*i];
      s[*i] = WRM[mb]->saturation(pc);
    }
  }
}


/* ******************************************************************
* Convertion s -> p.                                               
****************************************************************** */
void Richards_PK::DerivePressureFromSaturation(const Epetra_Vector& s, Epetra_Vector& p)
{
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm->WRM();  

  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = WRM[mb]->capillaryPressure(s[*i]);
      p[*i] = atm_pressure - pc;
    }
  }
}


/* ******************************************************************
* Clip pressure using pressure threshold.
****************************************************************** */
void Richards_PK::ClipHydrostaticPressure(double pmin, Epetra_MultiVector& p)
{
  for (int c = 0; c < ncells_owned; c++) p[0][c] = std::max(p[0][c], pmin);
}


/* ******************************************************************
* Clip pressure using constant saturation.
****************************************************************** */
void Richards_PK::ClipHydrostaticPressure(double pmin, double s0, Epetra_MultiVector& p)
{
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm->WRM();  

  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      if (p[0][*i] < pmin) {
        double pc = WRM[mb]->capillaryPressure(s0);
        p[0][*i] = atm_pressure - pc;
      }
    }
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi



