/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  PK for coupling of surface and subsurface transport PKs
*/

#ifndef ATS_AMANZI_COUPLEDTRANSPORT_PK_HH_
#define ATS_AMANZI_COUPLEDTRANSPORT_PK_HH_

#include "Teuchos_RCP.hpp"

#include "pk_mpcsubcycled.hh"
#include "PK.hh"

namespace Amanzi {

  class CoupledTransport_PK: public PK_MPCSubcycled_ATS {

  public: 
    CoupledTransport_PK(Teuchos::ParameterList& pk_tree_or_fe_list,
                     const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                     const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<TreeVector>& soln);
    ~CoupledTransport_PK(){}

    // PK methods
    // -- dt is the minimum of the sub pks
    virtual double get_dt();
    virtual void set_dt(double dt);

    // -- advance each sub pk from t_old to t_new.
    virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);

    virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);

    std::string name() { return "surface subsurface transport";} 

  private:

    void InterpolateCellVector(const Epetra_MultiVector& v0, const Epetra_MultiVector& v1, 
                               double dt_int, double dt, Epetra_MultiVector& v_int) ;


    // factory registration
    static RegisteredPKFactory<CoupledTransport_PK> reg_;
};

}  // namespace Amanzi
#endif
