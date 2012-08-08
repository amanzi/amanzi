Test Suite 1
~~~~~~~~~~~~~~

This is documentation for the Amanzi run that is intended to compare
the results from Amanzi against the semi-analytical results documented
in "A set of Analytical Benchmarks to Test Numerical Models of Flow
and Transport in Soils." by J. Vanderborght,
et. al. http://vzj.geoscienceworld.org/content/4/1/206.abstract

This page is a summary of the first test case noted as the first line
in Table 3 of that paper.  It is three cases that are all
"Steady-state flux in layered soil profiles".  The difference between
the three cases is the fact the first case a is loam over sand, b is
sand over loam and c is clay over sand.  The labeling of the case is
based on the Figure 2. results so files and data associated with case
2a should match Fig 2. a in the Vanderborght paper.

Geometry 
------------

.. image:: figs/Geometry.png

Stated Initial Conditions
------------------------------

Pressure :math:`\psi` for all x (depths) and t=0: -200cm (assume cm Hg) = 81747 Pa

Stated Boundary Conditions
------------------------------

Flow at surface for all time:  .5 cm/d = 5.78703704E-8 m/s 
:math:`\partial \psi / \partial x` at 200cm depth=0

Case 2a
^^^^^^^^

Case 2b
^^^^^^^^

Cace 2c
^^^^^^^^

