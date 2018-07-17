# -*- mode: cmake -*-

project (BUTLER_POD_2D)

add_amanzi_comparison_test(butler_pod_2d
                           INPUT amanzi_butler_pod_2d-u.xml
                           FILES mesh_cylinder.exo
                           TOLERANCE 2e-3
                           OBSERVATION Drawdown
                           REFERENCE ${BUTLER_POD_2D_SOURCE_DIR}/golden_output/observations.out)
