/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------


Author: Daniil Svyatsky

Default base with a few methods implemented for ATS
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_PHYSICAL_BDF_ATS_HH_
#define AMANZI_PK_PHYSICAL_BDF_ATS_HH_

#include "errors.hh"
#include "PK_Default.hh"
#include "PK_PhysicalBDF.hh"
#include "Operator.hh"
#include "Teuchos_TimeMonitor.hpp"
#include "BDF1_TI.hh"
#include "BDFFnBase.hh"

namespace Amanzi {

class PK_PhysicalBDF_ATS : public PK_PhysicalBDF {

public:
  PK_PhysicalBDF_ATS(){};

  PK_PhysicalBDF_ATS(Teuchos::ParameterList& pk_tree,
                     const Teuchos::RCP<Teuchos::ParameterList>& glist,
                     const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<TreeVector>& soln);

  PK_PhysicalBDF_ATS(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                     Teuchos::ParameterList& FElist,
                     const Teuchos::RCP<TreeVector>& solution);

// Virtual destructor
  virtual ~PK_PhysicalBDF_ATS(){};

  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);
  virtual void CommitStep (double t_old, double t_new, const Teuchos::RCP<State>& S);
  virtual double get_dt() { return dt_; }
  // -- transfer operators
  virtual void set_dt(double dt){ dt_ = dt;}
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u){};
  virtual bool ModifyCorrection(double h, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u){};
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up){};
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult ModifyCorrection(double h, 
                                                                             Teuchos::RCP<const TreeVector> res,
                                                                             Teuchos::RCP<const TreeVector> u,
                                                                             Teuchos::RCP<TreeVector> du) {};
    
  // Default preconditioner is Picard
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
    *Pu = *u;
    return 0;
  }

  // updates the preconditioner, default does nothing
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {}

  // Default implementations of BDFFnBase methods.
  // -- Compute a norm on u-du and return the result.
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du);


  // -- Experimental approach -- calling this indicates that the time
  //    integration scheme is changing the value of the solution in
  //    state.
  virtual void ChangedSolution();
  virtual double BoundaryValue(const Teuchos::RCP<const Amanzi::CompositeVector>& solution, int face_id);
  virtual void ApplyBoundaryConditions_(const Teuchos::Ptr<CompositeVector>& u);
  virtual void UpdateContinuationParameter(double lambda);

  // PC operator access
  Teuchos::RCP<Operators::Operator> preconditioner() { return preconditioner_; }

 // BC access
  std::vector<int>& bc_markers() { return bc_markers_; }
  std::vector<double>& bc_values() { return bc_values_; }
  Teuchos::RCP<Operators::BCs> BCs() { return bc_; }

protected: 

  void DeriveFaceValuesFromCellValues_(const Teuchos::Ptr<CompositeVector>& cv);

// error criteria
  Key conserved_key_;
  Key cell_vol_key_;
  double atol_, rtol_, fluxtol_;

 // PC
  Teuchos::RCP<Operators::Operator> preconditioner_;

  // BCs
  std::vector<int> bc_markers_;
  std::vector<double> bc_values_;
  Teuchos::RCP<Operators::BCs> bc_;

  // name of domain, associated mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::string domain_;

  // solution and evaluator
  std::string key_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> solution_evaluator_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> db_;

  //ENORM struct
  typedef struct ENorm_t {
    double value;
    int gid;
  } ENorm_t;

  // preconditioner assembly control
  bool assemble_preconditioner_;

  // timestep control
  double dt_;
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > time_stepper_;

  // timing
  Teuchos::RCP<Teuchos::Time> step_walltime_;

};

} // namespace

#endif

