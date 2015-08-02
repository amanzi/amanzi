/*
  This is the input component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <sstream>
#ifndef INPUT_CONVERTER_U_DEFS_HH__ 
#define INPUT_CONVERTER_U_DEFS_HH__

#define ATMOSPHERIC_PRESSURE            101325.0
#define GRAVITY_MAGNITUDE               9.80665
#define PARTICLE_DENSITY                1.0

#define VERBOSITY_DEFAULT               "Medium"

// all linear solvers
#define LINEAR_SOLVER_TOL               1.0e-16
#define LINEAR_SOLVER_MAXITER           100
#define LINEAR_SOLVER_METHOD            "gmres"
#define LINEAR_SOLVER_PC                "Hypre AMG"

// preconditioners: trilinos ml
#define TRILINOS_ML_SMOOTHER            "Jacobi"
#define TRILINOS_ML_AGG_THR             0.0
#define TRILINOS_ML_NSMOOTH             3
#define TRILINOS_ML_NCYC                2
#define TRILINOS_ML_OUTPUT              0
#define TRILINOS_ML_MAXLVLS             40
#define TRILINOS_ML_PRECTYPE            "MGV"
#define TRILINOS_ML_AGGTYPE             "Uncoupled-MIS"
#define TRILINOS_ML_AGGDAMP             1.33333
#define TRILINOS_ML_EIGENANAL_TYPE      "cg"
#define TRILINOS_ML_EIGENANAL_ITERS     10
#define TRILINOS_ML_SMOOTH_DAMP         1.0
#define TRILINOS_ML_SMOOTH_PRE_POST     "both"
#define TRILINOS_ML_SMOOTH_DAMP         1.0
#define TRILINOS_ML_CSOLVE_TYPE         "Amesos-KLU"
#define TRILINOS_ML_CSOLVE_MAX_SIZE     256

// preconditioners: hypre amg
#define HYPRE_AMG_TOL                   0.0
#define HYPRE_AMG_NCYC                  5
#define HYPRE_AMG_NSMOOTH               3
#define HYPRE_AMG_STR_THR               0.5

// preconditioners: trilinos ilu
#define TRILINOS_ILU_OLV                0
#define TRILINOS_ILU_RLXVAL             1.0
#define TRILINOS_ILU_RELTHR             1.0
#define TRILINOS_ILU_ABSTHR             0.0
#define TRILINOS_ILU_LVLFILL            0

#endif
