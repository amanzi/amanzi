2D Dispersive Transport in a Uniform Flow Field (diagonal)
==========================================================

Introduction
~~~~~~~~~~~~
This test problem is a slight modification of the previous test,
the direction of the velocity field is rotated by :math:`45^\circ` 
counter clockwise with respect to the mesh.


Defining Variables
~~~~~~~~~~~~~~~~~~
The background mesh consists of square cells with size :math:`H=15` m.
It has 83 grid cells in the x-direction and 37 grid cells in the y-direction. 
The mesh is gradually refined toward the source such that the well is
represented by a square cell of :math:`h=0.46875` [m] (:math:`h = H/32`).
The mesh refinement adds 8.4% more cells.

.. figure:: figures/mesh.png 
    :figclass: align-center

    **Computational mesh with 7468 cells.**


Results and Comparison
~~~~~~~~~~~~~~~~~~~~~~
The plume structure is characterized by three line cuts.
The first cut is given by line :math:`y=0` that goes through the well.
The two other cuts are given by lines :math:`x=0` and :math:`x=424`.

.. plot:: verification/transport/dispersion_45_point_2d/amanzi_dispersion_45_point_2d-a.py
   :align: center

.. include:: table_centerline.txt

The analytic data were computed with the AT123DAT software package.
A difference observed near downstream boundary for the Amanzi with the 
first-order transport scheme (Amanzi(1st)).
The second-order transport scheme provides excellent match.

.. plot:: verification/transport/dispersion_45_point_2d/amanzi_dispersion_45_point_2d-b.py
   :align: center


