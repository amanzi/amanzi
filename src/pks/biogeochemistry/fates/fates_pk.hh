/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Daniil Svyatsky (dasvyat@lanl.gov), Xu Chonggang
----------------------------------------------------------------------------- */

#ifndef PK_FATES_HH_
#define PK_FATES_HH_

#include "PK_Factory.hh"
#include "pk_physical_default.hh"
#include "ISO_Fortran_binding.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseVector.h"

#include "VerboseObject.hh"
#include "TreeVector.hh"

#include <string.h>

namespace Amanzi {
namespace Vegetation {

  typedef struct {  
    int nlevbed;
    int nlevdecomp;
    int patchno;
    int altmax_lastyear_indx_col;
    double temp_veg24_patch;
    double latdeg, londeg;
  } site_info;

#ifdef __cplusplus
  extern "C"
  {
#endif
    void init_ats_fates(int*, site_info*);
    void init_soil_depths(int*, int*, site_info*, double*, double*, double*, double*);
    void init_coldstart(int* );
    void fatesreadparameters(int*, CFI_cdesc_t * clm, CFI_cdesc_t * fates);
    void fatesreadpfts(int*, CFI_cdesc_t * clm, CFI_cdesc_t * fates);
    void dynamics_driv_per_site(int*, int*, site_info*, double*,
                                double*, double*, double*, double*,
                                double*, double*, double*, double*); 
  
#ifdef __cplusplus
  } // extern "C"
#endif 

  class FATES_PK : public PK_Physical_Default  {

  public:
    FATES_PK (Teuchos::ParameterList& FElist,
              const Teuchos::RCP<Teuchos::ParameterList>& plist,
              const Teuchos::RCP<State>& S,
              const Teuchos::RCP<TreeVector>& solution);
  
    // Virtual destructor
    virtual ~FATES_PK() {}

    // main methods
    // -- Initialize owned (dependent) variables.
    //virtual void setup(const Teuchos::Ptr<State>& S);
    virtual void Setup(const Teuchos::Ptr<State>& S);

    // -- Initialize owned (dependent) variables.
    //virtual void initialize(const Teuchos::Ptr<State>& S);
    virtual void Initialize(const Teuchos::Ptr<State>& S);


    // -- Update diagnostics for vis.
    virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S){}; 

    virtual bool AdvanceStep(double t_old, double t_new, bool reinit);
    
    // -- Commit any secondary (dependent) variables.
    virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);

      // -- provide a timestep size
    virtual double get_dt() {
      return dt_;
    }
    virtual void set_dt(double dt) {
      dt_ = dt;
    }

  protected:


    // void FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
    //                     Teuchos::Ptr<Epetra_SerialDenseVector> col_vec, bool copy);
    // void ColDepthDz_(AmanziMesh::Entity_ID col,
    //                  Teuchos::Ptr<Epetra_SerialDenseVector> depth,
    //                  Teuchos::Ptr<Epetra_SerialDenseVector> dz);

    
    double dt_;
    Teuchos::RCP<const AmanziMesh::Mesh> mesh_surf_;
    Key domain_surf_;
    Key trans_key_;

    int patchno_, nlevdecomp_;
    int ncells_owned_, ncells_per_col_;
    std::vector<site_info> site_;

  // factory registration
  static RegisteredPKFactory<FATES_PK> reg_;
};

}  // namespace Vegetation
}  // namespace Amanzi

#endif
