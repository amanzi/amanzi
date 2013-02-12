#ifndef INPUTPARSERISDEFAULTS__ 
#define INPUTPARSERISDEFAULTS__

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
#define ST_MAX_PREC_LAG                 5
#define ST_ERROR_REL_TOL                0.0
#define ST_ERROR_ABS_TOL                1.0
#define ST_MAX_DIVERGENT_ITERATIONS     3
#define ST_NONLIN_INIT_GUESS_EXTR_ORD   1
#define ST_NONLIN_INIT_TS_FACTOR        1.0
#define ST_NONLIN_INIT_TS_FACTOR_DAMP   1.0
#define ST_PRECOND                      "Hypre AMG"
#define ST_SOLVER                       "AztecOO"
#define ST_INIT_DARCY_BOOL              true
#define ST_DIVERG_FACT                  1000.0


#define TR_MAX_ITER                     15
#define TR_MIN_ITER                     10
#define TR_LIMIT_ITER                   20
#define TRANSIENT_NONLINEAR_TOLERANCE   1.0e-5
#define TR_NONLIN_DAMP                  1.0
#define TR_TS_RED_FACTOR                0.8
#define TR_TS_INC_FACTOR                1.2
#define TR_MAX_TS                       1.0e+8
#define TR_MAX_PREC_LAG                 5
#define TR_ERROR_REL_TOL                0.0
#define TR_ERROR_ABS_TOL                1.0
#define TR_MAX_DIVERGENT_ITERATIONS     3
#define TR_NONLIN_INIT_GUESS_EXTR_ORD   1
#define TR_NONLIN_INIT_TS_FACTOR        1.0
#define TR_NONLIN_INIT_TS_FACTOR_DAMP   1.0
#define TR_PRECOND                      "Hypre AMG"
#define TR_SOLVER                       "AztecOO"
#define TR_DIVERG_FACT                  1000.0

#define PIC_INIT_DARCY                  true
#define PIC_CLIP_SAT                    0.9
#define PICARD_TOLERANCE                1.0e-8
#define PIC_MAX_ITER                    400
#define PIC_PRECOND                     "Hypre AMG"
#define PIC_SOLVE                       "AztecOO"
#define PIC_ERROR_METHOD                "pressure"
#define PIC_METHOD                      "Picard"

#define LIN_SOLVE_TOL                   1.0e-16
#define LIN_SOLVE_MAXITER               100
#define LIN_SOLVE_METHOD                "GMRES"
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

#endif
