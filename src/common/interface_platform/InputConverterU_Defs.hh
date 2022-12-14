/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (previos version, ParserIS)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Input Converter

*/

#include <sstream>
#ifndef INPUT_CONVERTER_U_DEFS_HH__
#  define INPUT_CONVERTER_U_DEFS_HH__

// universal physical constants
#  define ATMOSPHERIC_PRESSURE 101325.0
#  define GRAVITY_MAGNITUDE 9.80665

// default material constants
#  define PARTICLE_DENSITY 1.0

#  define ELL_MUALEM 0.5
#  define ELL_BURDINE 2.0

// other default constats
#  define MAXIMUM_TIMESTEP 4.3234e+17
#  define MINIMUM_TIMESTEP 1.0e-6
#  define RESTART_TIMESTEP 1.0

#  define VERBOSITY_DEFAULT "Medium"

// all linear solvers
#  define LINEAR_SOLVER_TOL 1.0e-15
#  define LINEAR_SOLVER_MAXITER 100
#  define LINEAR_SOLVER_METHOD "gmres"
#  define LINEAR_SOLVER_PC "Hypre AMG"

// preconditioners: trilinos ml
#  define TRILINOS_ML_SMOOTHER "Jacobi"
#  define TRILINOS_ML_AGG_THR 0.0
#  define TRILINOS_ML_NSMOOTH 3
#  define TRILINOS_ML_NCYC 2
#  define TRILINOS_ML_OUTPUT 0
#  define TRILINOS_ML_MAXLVLS 40
#  define TRILINOS_ML_PRECTYPE "MGV"
#  define TRILINOS_ML_AGGTYPE "Uncoupled-MIS"
#  define TRILINOS_ML_AGGDAMP 1.33333
#  define TRILINOS_ML_EIGENANAL_TYPE "cg"
#  define TRILINOS_ML_EIGENANAL_ITERS 10
#  define TRILINOS_ML_SMOOTH_DAMP 1.0
#  define TRILINOS_ML_SMOOTH_PRE_POST "both"
#  define TRILINOS_ML_SMOOTH_DAMP 1.0
#  define TRILINOS_ML_CSOLVE_TYPE "Amesos-KLU"
#  define TRILINOS_ML_CSOLVE_MAX_SIZE 256

// preconditioners: hypre amg
#  define HYPRE_AMG_TOL 0.0
#  define HYPRE_AMG_NCYC 5
#  define HYPRE_AMG_NSMOOTH 3
#  define HYPRE_AMG_STR_THR 0.5

// preconditioners: trilinos ilu
#  define TRILINOS_ILU_OLV 0
#  define TRILINOS_ILU_RLXVAL 1.0
#  define TRILINOS_ILU_RELTHR 1.0
#  define TRILINOS_ILU_ABSTHR 0.0
#  define TRILINOS_ILU_LVLFILL 0

// all nolninear solvers
#  define NONLINEAR_TOLERANCE 1.0e-5
#  define INC_DIVERG_FACTOR 1000.0
#  define MAX_DIVERG_ITERATIONS 3

// nonlinear solver: NKA
#  define NKA_DIVERG_TOL 1.0e10
#  define NKA_LIMIT_ITERATIONS 20
#  define NKA_NUM_VECTORS 10

// nonlinear solver: Picard
#  define PICARD_SOLVER "AztecOO"
#  define PICARD_TOLERANCE 1.0e-8
#  define PICARD_MAX_ITERATIONS 400
#  define PIC_PRECONDITIONER "Hypre AMG"

#  define LIN_SOLVE_TOL 1.0e-15
#  define LIN_SOLVE_MAXITER 100
// time integrator
#  define TI_SOLVER "AztecOO"
#  define TI_PRECONDITIONER "Hypre AMG"
#  define TI_PLAMBDA_SOLVER "GMRES with Hypre AMG"
#  define TI_TIMESTEP_CONTROLLER "standard"
#  define TI_MIN_ITERATIONS 10
#  define TI_MAX_ITERATIONS 15
#  define TI_TS_INCREASE_FACTOR 1.2
#  define TI_TS_REDUCTION_FACTOR 0.8
#  define TI_MAX_PC_LAG 5
#  define TI_TOL_RELAX_FACTOR 1.0
#  define TI_TOL_RELAX_FACTOR_DAMPING 1.0

// other constants
#  define FLOW_STEADY_REGIME 1
#  define FLOW_TRANSIENT_REGIME 2
#  define FLOW_BOTH_REGIMES 3
#  define TRANSPORT_SUBCYCLING true

#  define MOLAR_MASS_WATER 0.0180153333333

#endif
