/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $AMANZI_DIR/COPYRIGHT
Author: ??

Effectively stolen from Amanzi, with few modifications.

Keys required by Monty:


Fetch
- temperature
- saturation_ice
- saturation_liquid




Put
(timestep)
- <energy flux to surface>
- <water flux to surface>
- <water source in subsurface>  <-- roots


(initialization)
- temperature
- saturation_ice
- saturation_liquid
- saturation_gas   (1 - s_i - s_l)


------------------------------------------------------------------------- */

#ifndef ATS_CLM_DRIVER_HH
#define ATS_CLM_DRIVER_HH

#include <Domain.hh>
#include <GeometricModel.hh>
#include <State.hh>
#include <coordinator.hh>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Simulator.hh"

namespace Amanzi {

struct ATSCLMDriver
  : public Teuchos::VerboseObject<ATSCLMDriver>
{
  int32_t Initialize(const MPI_Comm& mpi_comm,
                  int * col_types,
                  int num_cols,
                  int num_types);

  int32_t Finalize();

  int32_t Advance(double dt, bool force_vis=false);

  int32_t SetInitCLMData(double* T, double* sl, double* si);
  int32_t SetCLMData(double* e_flux, double* w_flux);
  int32_t GetCLMData(double* T, double* sl, double* si);

 protected:
  // size of data
  int ncells_surf_;
  int ncells_sub_;

  // debug
  bool coord_setup_, coord_init_;

  // ATS internals
  Teuchos::RCP<AmanziGeometry::Domain> simdomain_;
  Teuchos::RCP<AmanziGeometry::GeometricModel> geom_model_;
  Teuchos::RCP<State> S_;
  Teuchos::RCP<State> S_next_;
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<Coordinator> coordinator_;

  // Maps from CLM to ATS
  Teuchos::RCP<Epetra_Map> surf_clm_map_;
  Teuchos::RCP<Epetra_Map> sub_clm_map_;
  Teuchos::RCP<Epetra_Import> surf_importer_;
  Teuchos::RCP<Epetra_Import> sub_importer_;

  Teuchos::EVerbosityLevel verbosity_;

  // list of region names
  std::vector<std::string> clm_type_region_names;

 private:
  int32_t SetData_(std::string key, double* data, int length);
  int32_t GetData_(std::string key, double* data, int length);
  int32_t SetSurfaceData_(std::string key, double* data, int length);
  int32_t GetSurfaceData_(std::string key, double* data, int length);

};

} // namespace Amanzi

#endif // ATS_CLM_DRIVER_HH
