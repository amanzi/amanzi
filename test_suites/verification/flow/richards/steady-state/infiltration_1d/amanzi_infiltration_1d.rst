Infiltration 1D
===============

Capabilities Tested
-------------------

For details on this test, see :ref:`about_sand_clay`.


Background
----------

Verification problems from the literature have been identified to test
isothermal, single-phase, variably saturated flow.  We initially focus
on test problems that address the two most widely used k-s-p
functions, Mualem-van Genuchten :cite:`scinfil-Mualem_1976` :cite:`scinfil-vanGenuchten_1980` and Brooks-Corey :cite:`scinfil-brooks1964hydraulic`.  These include
steady-state and transient tests with Dirichlet and Neumann boundary
conditions.

This is documentation for the Amanzi run that is intended to compare
the results from Amanzi against the semi-analytical results documented
in "A set of Analytical Benchmarks to Test Numerical Models of Flow
and Transport in Soils." by J. Vanderborght,
et. al. http://vzj.geoscienceworld.org/content/4/1/206.abstract :cite:`scinfil-vanderborght2005set`.

This page is a summary of the first test case noted as the first line
in Table 3 of that paper.  It is three cases that are all
"Steady-state flux in layered soil profiles".  The difference between
the three cases is as follows:
case #1 is 0.5 m of clay and 1.5 m of sand;
case #2 is 0.5 m of loam and 1.5 m of sand;
case #3 is 1.5 m of loam and 0.5 m of sand.
The results should match that in the Vanderborght paper.


Model
-----

Initial condition.
Pressure :math:`\psi` for all x (depths) and t=0: -200cm (assume cm Hg) = 81747 Pa

Boundary conditions. 
The pressure at the bottom z=0m is 99630.6336 Pa.
Flow at the opposite end, z=2m is fixed at 0.5 cm/d = 5.78703704E-8 m/s.

.. image:: geometry.png
  :align: center
  :width: 200px


Problem Specification
---------------------

We use a box domain with hight 2 m. The other dimenstions are equal to 1 m.


Schematic
~~~~~~~~~


Mesh
~~~~

We consider a column mesh with 200 cells in the vertical direction.


Variables
~~~~~~~~~


Case #1: Sand Clay Layers
-------------------------

The steady-state solution is shown below.


Results and Comparison
~~~~~~~~~~~~~~~~~~~~~~

.. plot:: verification/flow/richards/steady-steate/infiltration_clay_sand_1d/amanzi_infiltration_1d-c.py
   :align: center


Case #2 Loam Sand Layers
------------------------

The steady-state solution is shown below.


Results and Comparison
~~~~~~~~~~~~~~~~~~~~~~

.. plot:: verification/flow/richards/steady-steate/infiltration_loam_sand_1d/amanzi_infiltration_1d-a.py
   :align: center


Case #3: Sand Loam Layers
-------------------------

The steady-state solution is shown below.


Results and Comparison
~~~~~~~~~~~~~~~~~~~~~~

.. plot:: verification/flow/richards/steady-steate/infiltration_loam_sand_1d/amanzi_infiltration_1d-b.py
   :align: center


References
----------

.. bibliography:: /bib/ascem.bib
   :filter: docname in docnames
   :style:  alpha
   :keyprefix: scinfil-

.. _about_sand_clay:


About
-----

* Directory:  testing/verification/flow/richards/steady-state/infiltration_1d

* Author:  

* Maintainer:  David Moulton (moulton@lanl.gov)

* Input Files:

  * amanzi_infiltration_clay_sand_1d-u.xml
  * amanzi_infiltration_loam_sand_1d-u.xml
  * amanzi_infiltration_sand_loam_1d-u.xml

    * Spec Version 2.3, unstructured mesh framework
    * mesh:  generated internally 


