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
#include "Transport_PK.hh"

namespace Amanzi {
namespace Transport {

/* *******************************************************************
* Calculate a dispersive tensor the from Darcy fluxes. The flux is
* assumed to be scaled by face area.
******************************************************************* */
void Transport_PK::CalculateDispersionTensor_(
    const Epetra_MultiVector& darcy_flux, 
    const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation)
{
  D.resize(ncells_owned);
  for (int c = 0; c < ncells_owned; c++) D[c].init(dim, 2);

  for (int mb = 0; mb < dispersion_models_.size(); mb++) {
    Teuchos::RCP<DispersionModel> spec = dispersion_models_[mb]; 

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
          D[*c] *= porosity[0][*c] * saturation[0][*c];
        } else {
          WhetStone::MFD3D_Diffusion mfd3d(mesh_);

          AmanziMesh::Entity_ID_List faces;
          AmanziGeometry::Point velocity(dim);

          mesh_->cell_get_faces(*c, &faces);
          int nfaces = faces.size();

          std::vector<double> flux(nfaces);
          for (int n = 0; n < nfaces; n++) flux[n] = darcy_flux[0][faces[n]];
          mfd3d.RecoverGradient_MassMatrix(*c, flux, velocity);
          velocity /= porosity[0][*c];  // pore velocity

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

          D[*c] *= porosity[0][*c] * saturation[0][*c];
        }
      }
    }
  }
}


/* *******************************************************************
* Calculate diffusion tensor if no dispersion is given.
******************************************************************* */
int Transport_PK::CalculateDiffusionTensor_(
    const std::string component_name,
    const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation)
{
  if (diffusion_models_ == Teuchos::null) return -1;

  double md = diffusion_models_->FindComponentValue(component_name);
  if (md == 0.0) return -1;

  D.resize(ncells_owned);
  for (int c = 0; c < ncells_owned; c++) { 
    D[c].init(dim, 1);
    D[c](0, 0) = md * porosity[0][c] * saturation[0][c];
  }
}


/* *******************************************************************
* Add molecular diffusion to the existing dispersive tensor.
******************************************************************* */
void Transport_PK::AddMolecularDiffusion_(
    const std::string component_name,
    const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation)
{
  if (diffusion_models_ == Teuchos::null) return;

  double md = diffusion_models_->FindComponentValue(component_name);
  if (md == 0.0) return;

  for (int mb = 0; mb < dispersion_models_.size(); mb++) {
    Teuchos::RCP<DispersionModel> spec = dispersion_models_[mb]; 

    std::vector<AmanziMesh::Entity_ID> block;
    for (int r = 0; r < (spec->regions).size(); r++) {
      std::string region = (spec->regions)[r];
      mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

      AmanziMesh::Entity_ID_List::iterator c;
      for (c = block.begin(); c != block.end(); c++) {
        D[*c] += md * spec->tau * porosity[0][*c] * saturation[0][*c];
      }
    }
  }
}


/* ******************************************************************
* Adds time derivative to the cell-based part of MFD algebraic system.
****************************************************************** */
/*
void Dispersion::AddTimeDerivative(
    double dT, const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation)
{
  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);

  for (int c = 0; c < ncells_owned; c++) {
    double volume = mesh_->cell_volume(c);
    double factor = volume * porosity[0][c] * saturation[0][c] / dT;

    int c_GID = cmap_wghost.GID(c);
    App_->SumIntoGlobalValues(1, &c_GID, &factor);
  }
}
*/


}  // namespace Transport
}  // namespace Amanzi



