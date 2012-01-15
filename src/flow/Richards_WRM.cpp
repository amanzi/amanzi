/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (nnc@lanl.gov), 
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Richards_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Defines relative permeability ONLY for cells.                                               
****************************************************************** */
void Richards_PK::calculateRelativePermeability(const Epetra_Vector& p)
{
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();

    //AmanziMesh::Set_ID_List block;
    std::vector<AmanziMesh::Set_ID> block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) (*Krel_cells)[*i] = WRM[mb]->k_relative(p[*i]);
  }
}


/* ******************************************************************
* Defines relative permeabilities ONLY for faces. 
* WARNING; we assume that K is a diagonal tensor. (lipnikov@lanl.gov)                                              
****************************************************************** */
void Richards_PK::calculateRelativePermeabilityUpwindGravity(const Epetra_Vector& p)
{
  calculateRelativePermeability(p);  // populates cell-based permeabilities
  FS->copyMasterCell2GhostCell(*Krel_cells);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c=0; c<number_owned_cells; c++) {
    mesh_->cell_get_face_dirs(c, &dirs);
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n=0; n<nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f); 
      if ((normal * gravity) * dirs[n] >= 0.0) (*Krel_faces)[f] = (*Krel_cells)[c];
      else if (bc_markers[f] != FLOW_BC_FACE_NULL) (*Krel_faces)[f] = (*Krel_cells)[c];
    }
  }
}


/* ******************************************************************
* Use analytical formula for derivative dS/dP.                                               
****************************************************************** */
void Richards_PK::derivedSdP(const Epetra_Vector& p, Epetra_Vector& ds)
{
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);

    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
      
    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) ds[*i] = WRM[mb]->d_saturation(p[*i]);
  }
}


/* ******************************************************************
* Convertion p -> s.                                               
****************************************************************** */
void Richards_PK::deriveSaturationFromPressure(const Epetra_Vector& p, Epetra_Vector& s)
{
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);

    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
      
    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) s[*i] = WRM[mb]->saturation(p[*i]);
  }
}


/* ******************************************************************
* Convertion s -> p.                                               
****************************************************************** */
void Richards_PK::derivePressureFromSaturation(double s, Epetra_Vector& p)
{
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);

    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
      
    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) p[*i] = WRM[mb]->pressure(s);
  } 
 
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);

    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
      
    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) p[*i] = WRM[mb]->pressure(s);
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi



