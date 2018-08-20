/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  PK for coupling of surface and subsurface transport PKs
*/

#ifndef ATS_AMANZI_MORPHOLOGY_PK_HH_
#define ATS_AMANZI_MORPHOLOGY_PK_PK_HH_

#include "Teuchos_RCP.hpp"

#include "weak_mpc.hh"
#include "pk_physical_bdf_default.hh"
#include "PK.hh"
#include "Debugger.hh"

namespace Amanzi {

  class Morphology_PK: public WeakMPC{

  public: 
    Morphology_PK(Teuchos::ParameterList& pk_tree_or_fe_list,
                     const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                     const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<TreeVector>& soln);
    ~Morphology_PK(){}

    // PK methods
    // -- dt is the minimum of the sub pks
    virtual double get_dt();
    //virtual void set_dt(double dt);
    virtual void Setup(const Teuchos::Ptr<State>& S);
    virtual void Initialize(const Teuchos::Ptr<State>& S);

    // -- advance each sub pk from t_old to t_new.
    virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);

    //virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);

    std::string name() { return name_;} 

  protected:

    Teuchos::RCP<PK_PhysicalBDF_Default> flow_pk_;
    Teuchos::RCP<PK_PhysicalBDF_Default> sed_transport_pk_;

    // debugger for dumping vectors
    Teuchos::RCP<Debugger> flow_db_;
    Teuchos::RCP<Debugger> trans_db_;
    
    // factory registration
    static RegisteredPKFactory<Morphology_PK> reg_;
};

}  // namespace Amanzi
#endif
