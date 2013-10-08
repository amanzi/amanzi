/*
  This is the transport component of Amanzi. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Teuchos_RCP.hpp"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

#include "mfd3d_diffusion.hh"
#include "nlfv.hh"
#include "tensor.hh"
#include "PreconditionerFactory.hh"

#include "TransportDefs.hh"
#include "Dispersion.hh"

namespace Amanzi {
namespace AmanziTransport {

/* *******************************************************************
* Initialization of a class.
******************************************************************* */
void Dispersion::Init()
{
  ncells_owned  = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  nfaces_owned  = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  dim = mesh_->space_dimension();

  D.resize(ncells_owned);
  for (int c = 0; c < ncells_owned; c++) D[c].init(dim, 2);
}


/* *******************************************************************
* Initialization of a class.
******************************************************************* */
void Dispersion::Init(
    std::vector<Teuchos::RCP<DispersionModel> >* specs,
    Teuchos::RCP<const AmanziMesh::Mesh> mesh, Teuchos::RCP<Transport_State> TS) 
{
  specs_ = specs;
  mesh_ = mesh; 
  TS_ = TS;
}


/* *******************************************************************
* Initialization of a preconditioner
******************************************************************* */
void Dispersion::InitPreconditioner(
    const std::string& prec_name, const Teuchos::ParameterList& prec_list)
{
  AmanziPreconditioners::PreconditionerFactory factory;
  preconditioner_ = factory.Create(prec_name, prec_list);
}


/* *******************************************************************
 * Calculate a dispersive tensor the from Darcy fluxes. The flux is
 * assumed to be scaled by face area.
 ****************************************************************** */
void Dispersion::CalculateDispersionTensor(
    const Epetra_Vector& darcy_flux, 
    const Epetra_Vector& porosity, const Epetra_Vector& saturation)
{
  for (int mb = 0; mb < specs_->size(); mb++) {
    Teuchos::RCP<DispersionModel> spec = (*specs_)[mb]; 

    std::vector<AmanziMesh::Entity_ID> block;
    for (int r = 0; r < (spec->regions).size(); r++) {
      std::string region = (spec->regions)[r];
      mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

      AmanziMesh::Entity_ID_List::iterator c;
      for (c = block.begin(); c != block.end(); c++) {
        D[*c].PutScalar(0.0); 
        if (spec->model == TRANSPORT_DISPERSIVITY_MODEL_ISOTROPIC) {
          for (int i = 0; i < dim; i++) {
            D[*c](i, i) = spec->alphaL + spec->D * spec->tau;
          }
          D[*c] *= porosity[*c] * saturation[*c];
        } else {
          WhetStone::MFD3D_Diffusion mfd3d(mesh_);

          AmanziMesh::Entity_ID_List faces;
          std::vector<int> dirs;
          AmanziGeometry::Point velocity(dim);

          mesh_->cell_get_faces_and_dirs(*c, &faces, &dirs);
          int nfaces = faces.size();

          std::vector<double> flux(nfaces);
          for (int n = 0; n < nfaces; n++) flux[n] = darcy_flux[faces[n]];
          mfd3d.RecoverGradient_MassMatrix(*c, flux, velocity);
          velocity /= porosity[*c];  // pore velocity

          double velocity_value = norm(velocity);
          double anisotropy = spec->alphaL - spec->alphaT;

          for (int i = 0; i < dim; i++) {
            D[*c](i, i) = spec->D * spec->tau + spec->alphaT * velocity_value;
            for (int j = i; j < dim; j++) {
              double s = anisotropy * velocity[i] * velocity[j];
              if (velocity_value) s /= velocity_value;
              D[*c](j, i) = D[*c](i, j) += s;
            }
          }

          D[*c] *= porosity[*c] * saturation[*c];
        }
      }
    }
  }
}


/* ******************************************************************
* Adds time derivative to the cell-based part of MFD algebraic system.
****************************************************************** */
void Dispersion::AddTimeDerivative(
    double dT, const Epetra_Vector& porosity, const Epetra_Vector& saturation)
{
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);

  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = volume * porosity[c] * saturation[c] / dT;

    int c_GID = cmap_wghost.GID(c);
    App_->SumIntoGlobalValues(1, &c_GID, &factor);
  }
}


}  // namespace AmanziTransport
}  // namespace Amanzi



