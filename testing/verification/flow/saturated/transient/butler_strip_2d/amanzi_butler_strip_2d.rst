Transient 2D Heterogeneous Strip Confined Flow to a Pumping Well
================================================================

Capabilties Tested
------------------


Background
----------

This file documents Amanzi code verification of flow simulations for
three selected cases: (1) an infinite strip aquifer embedded between
two other aquifers, (2) a heterogeneous cylinder (or disc in 2D case)
embedded in a background aquifer, and (3) uniform aquifers in bounded
domains with fixed heads at the upper- and down-stream directions and
no-flow at lateral boundaries.  The original FORTRAN codes that solve
the analytical solutions for cases (1) and (2) were obtained from Greg
Ruskauf. Two modifications were made for compilation on Linux
systems: (a) The keyword ‘mode’ in OPEN statements was removed because
it is not standard, and (b) variables in COMMON blocks were
re-arranged to eliminate a warning message (padding of additional 4
bytes).  These codes were compiled using the gfortran compiler on the
Linux environment. Note that these analytical solutions were derived
for two-dimensional problems only, and therefore in numerical
simulations, we used three-dimensional domains with unit thickness for
comparison.

Model
-----


Problem Specification
---------------------

Schematic
~~~~~~~~~

.. figure:: schematic/butler_strip_schematic.png
    :figclass: align-center
    :width: 600 px

    ** Schematic of the Butler and Liu Linear Strip verification problem **

Variables
~~~~~~~~~


Results and Comparison
----------------------


References
----------


About
-----

* Directory: testing/verification/flow/transient/butler_strip_2d

* Authors:  Dylan Harp, Zhiming Lu

* Maintainer(s): 

* Input Files:

Status
~~~~~~

Add notes here about the status of the test.  

.. todo:: 

  * Documentation:
