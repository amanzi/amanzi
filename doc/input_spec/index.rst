.. Amanzi documentation master file, created by
   sphinx-quickstart on Thu Mar 12 14:16:05 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Amanzi
$$$$$$

.. image:: images/amanzi.png

Amanzi provides a flexible and extensible parallel flow and reactive transport simulation 
capability for environmental applications. 
It includes general polyhedral mesh infrastructure, various discretizations of process models, 
including traditional finite volume schemes, mimetic finite differences, and nonlinear finite volumes. 
In addition, it provides advanced nonlinear solvers, such as Nonlinear Krylov Acceleration (NKA) 
and Anderson Acceleration, and leverages Trilinos-ML and Hypre Algebraic Multigrid for scalable 
solvers. 
The reaction of contaminants and minerals carried by the flow through the surrounding rock and 
soil is modeled through a general software interface called Alquimia that allows Amanzi to 
interface with a variety of powerful geochemistry engines including PFLOTRAN and CrunchFlow. 
The code is parallel and leverages open-source parallel frameworks such as Trilinos and PETSc. 
Amanzi has been used to model contaminant migration at various DOE waste sites (e.g., Nevada 
National Security Site, and Savannah River), and is generally applicable to groundwater 
contaminant migration under partially saturated, nonisothermal conditions in fractured media.

.. toctree::
   :maxdepth: 2
   :caption: User Guide:

   input_specs

..   
.. toctree::
   :maxdepth: 2
   :caption: Developers Guide:
             

Indices and tables
%%%%%%%%%%%%%%%%%%

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
