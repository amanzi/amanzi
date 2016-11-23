/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for a StrongMPC which uses a preconditioner in which the
block-diagonal cell-local matrix is dense.  If the system looks something
like:

A( y1, y2, x, t ) = 0
B( y1, y2, x, t ) = 0

where y1,y2 are spatially varying unknowns that are discretized using the MFD
method (and therefore have both cell and face unknowns), an approximation to
the Jacobian is written as

[  dA_c/dy1_c  dA_c/dy1_f   dA_c/dy2_c       0      ]
[  dA_f/dy1_c  dA_f/dy1_f      0              0      ]
[  dB_c/dy1_c     0          dB_c/dy2_c  dB_c/dy2_f ]
[      0           0          dB_f/dy2_c  dB_f/dy2_f ]


Note that the upper left block is the standard preconditioner for the A
system, and the lower right block is the standard precon for the B system,
and we have simply added cell-based couplings, dA_c/dy2_c and dB_c/dy1_c.

In the temperature/pressure system, these correspond to d_water_content /
d_temperature and d_energy / d_pressure.

------------------------------------------------------------------------- */

#ifndef MPC_COUPLED_CELLS_HH_
#define MPC_COUPLED_CELLS_HH_

#include "pk_physical_bdf_default.hh"
#include "strong_mpc.hh"

namespace Amanzi {

namespace Operators { class TreeOperator; }

class MPCCoupledCells : public StrongMPC<PK_PhysicalBDF_Default> {
 public:

  MPCCoupledCells(Teuchos::ParameterList& FElist,
                  const Teuchos::RCP<Teuchos::ParameterList>& plist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& solution):
    PK(FElist, plist, S, solution),
    StrongMPC<PK_PhysicalBDF_Default>(FElist, plist, S, solution),
    decoupled_(false) {}

  virtual void Setup(const Teuchos::Ptr<State>& S);

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

 protected:
  Key dA_dy2_key_;
  Key dB_dy1_key_;
  Key A_key_;
  Key y2_key_;
  Key B_key_;
  Key y1_key_;
  Teuchos::RCP<Operators::TreeOperator> preconditioner_;
  Teuchos::RCP<Operators::TreeOperator> linsolve_preconditioner_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  bool decoupled_;

  // cruft for easier global debugging
  Teuchos::RCP<Debugger> db_;
};


} //  namespace
#endif
