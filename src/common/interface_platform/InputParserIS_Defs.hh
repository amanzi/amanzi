#ifndef INPUT_PARSERIS_DEFS_HH__ 
#define INPUT_PARSERIS_DEFS_HH__

namespace Amanzi {
namespace AmanziInput {

#define AMANZI_OLD_INPUT_VERSION_MAJOR 1
#define AMANZI_OLD_INPUT_VERSION_MINOR 2
#define AMANZI_OLD_INPUT_VERSION_MICRO 2

}  // namespace AmanziInput
}  // namespace Amanzi


// define defaults that are used in the input parser here
#define ATMOSPHERIC_PRESSURE            101325.0


#define GRAVITY_MAGNITUDE               9.80665
#define PARTICLE_DENSITY                1.0
#define CHEM_TRANS_DT_RATIO             1.0
#define MAX_DIVERGENT_ITERATIONS        3      // default for native spec

#define VERBOSITY_DEFAULT               "Medium"

#define USE_PICARD                      false

#define TRANSPORT_SUBCYCLING            true

#define ELL_MUALEM                      0.5
#define ELL_BURDINE                     2.0

#define ST_MAX_ITER                     15
#define ST_MIN_ITER                     10
#define ST_LIMIT_ITER                   20
#define STEADY_NONLINEAR_TOLERANCE      1.0e-5
#define ST_NONLIN_DAMP                  1.0
#define ST_TS_RED_FACTOR                0.8
#define ST_TS_INC_FACTOR                1.2
#define ST_MAX_TS                       1.0e+10
#define ST_MIN_TS                       1.0e-20
#define ST_MAX_PREC_LAG                 5
#define ST_ERROR_REL_TOL                0.0
#define ST_ERROR_ABS_TOL                1.0
#define ST_MAX_DIVERGENT_ITERATIONS     3
#define ST_NONLIN_INIT_GUESS_EXTR_ORD   1
#define ST_NONLIN_INIT_TS_FACTOR        1.0
#define ST_NONLIN_INIT_TS_FACTOR_DAMP   1.0
#define ST_PRECOND                      "Hypre AMG"
#define ST_SOLVER                       "AztecOO"
#define ST_INIT_SOLVER                  "AztecOO"
#define ST_PLAMB_SOLVER                 "AztecOO"
#define ST_INIT_DARCY_BOOL              true
#define ST_DIVERG_FACT                  1000.0
#define ST_SP_DT_INCR_FACTOR            1.0    // this is the dt increase factor for single phase
#define ST_CLIP_SAT                     0.6
#define ST_NKA_DIVGD_TOL                1.0e10
#define ST_NKA_NUMVEC                   10
#define ST_TS_STRATEGY                  "standard"
#define ST_TS_CONTROLLER                "standard"


#define TR_MAX_ITER                     15
#define TR_MIN_ITER                     10
#define TR_LIMIT_ITER                   20
#define TRANSIENT_NONLINEAR_TOLERANCE   1.0e-5
#define TR_NONLIN_DAMP                  1.0
#define TR_TS_RED_FACTOR                0.8
#define TR_TS_INC_FACTOR                1.2
#define TR_MAX_TS                       1.0e+8
#define TR_MIN_TS                       1.0e-20
#define TR_MAX_PREC_LAG                 5
#define TR_ERROR_REL_TOL                0.0
#define TR_ERROR_ABS_TOL                1.0
#define TR_MAX_DIVERGENT_ITERATIONS     3
#define TR_NONLIN_INIT_GUESS_EXTR_ORD   1
#define TR_NONLIN_INIT_TS_FACTOR        1.0
#define TR_NONLIN_INIT_TS_FACTOR_DAMP   1.0
#define TR_PRECOND                      "Hypre AMG" 
#define TR_SOLVER                       "AztecOO"
#define TR_SOLVER_DARCY                 "AztecOO"
#define TR_INIT_SOLVER                  "AztecOO"
#define TR_PLAMB_SOLVER                 "AztecOO"
#define TR_DIVERG_FACT                  1000.0
#define TR_SP_DT_INCR_FACTOR            1.0    // this is the dt increase factor for single phase
#define TR_CLIP_SAT                     0.6
#define TR_NKA_DIVGD_TOL                1.0e10
#define TR_NKA_NUMVEC                   10
#define TR_TS_STRATEGY                  "standard"
#define TR_TS_CONTROLLER                "standard"
#define TR_INIT_DARCY_BOOL              false

 
#define PIC_INIT_DARCY                  true
#define PIC_CLIP_SAT                    0.9
#define PIC_CLIP_PRESSURE               50000.0
#define PICARD_TOLERANCE                1.0e-8
#define PIC_MAX_ITER                    400
#define PIC_PRECOND                     "Hypre AMG"
#define PIC_SOLVE                       "AztecOO"
#define PIC_ERROR_METHOD                "pressure"
#define PIC_METHOD                      "Picard"

#define LIN_SOLVE_TOL                   1.0e-16
#define LIN_SOLVE_MAXITER               100
#define LIN_SOLVE_METHOD                "gmres"
#define LIN_SOLVE_PREC                  "Hypre AMG"

#define ML_SMOOTHER                     "Jacobi"
#define ML_AGG_THR                      0.0
#define ML_NSMOOTH                      3
#define ML_NCYC                         2
#define ML_OUTPUT                       0
#define ML_MAXLVLS                      40
#define ML_PRECTYPE                     "MGV"
#define ML_AGGTYPE                      "Uncoupled-MIS"
#define ML_AGGDAMP                      1.33333
#define ML_EIGENANAL_TYPE               "cg"
#define ML_EIGENANAL_ITERS              10
#define ML_SMOOTH_DAMP                  1.0
#define ML_SMOOTH_PRE_POST              "both"
#define ML_SMOOTH_DAMP                  1.0
#define ML_CSOLVE_TYPE                  "Amesos-KLU"
#define ML_CSOLVE_MAX_SIZE              256

#define AMG_TOL                         0.0
#define AMG_NCYC                        5
#define AMG_NSMOOTH                     3
#define AMG_STR_THR                     0.5

#define ILU_OLV                         0
#define ILU_RLXVAL                      1.0
#define ILU_RELTHR                      1.0
#define ILU_ABSTHR                      0.0
#define ILU_LVLFILL                     0

#define BCHYDRST_COORD                  "Absolute"

#define MAXIMUM_TIME_STEP               4.3234e+17
#define RESTART_TIME_STEP               1.0
#define TI_RESCUE_REDUCTION_FACTOR      0.5

#define STEADY_REGIME                   1
#define TRANSIENT_REGIME                2
#define BOTH_REGIMES                    3

#endif
