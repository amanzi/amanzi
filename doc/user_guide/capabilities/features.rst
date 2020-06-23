Amanzi Features
=================================

*Amanzi* is
a full featured multiprocess simulator that can run efficiently on a wide variety
of platforms from a laptop to high performance computers utilizing
thousands of cores. It is designed to allow analysts to perform quick
field calculations or high fidelity massively parallel simulations as
the situation dictates. *Amanzi* is built in a modular manner using
the best numerical analysis tools such as state-of-the art solvers
from Sandia and Livermore National labs, mesh infrastructure from Los
Alamos, Sandia, Berkeley and Argonne National Laboratories as well as
many other established open source software components. This allows
*Amanzi* to be robust, maintanable and extensible. *Amanzi* is
extensively tested using an agile unit testing system as well
benchmarked against established subsurface flow and transport codes
from around the DOE complex.

Some of the salient features of *Amanzi* are listed below:

1. *Amanzi* can handle structured meshes (with adaptive mesh
   refinement) or unstructured meshes (including polyhedral elements).
2. It can internally generate a mesh or read a previously generated
   unstructured mesh in the `Exodus II <http://sourceforge.net/projects/exodusii/>`_ format.
3. Simulations can be performed representing the mesh cells in 2 and 3 dimensions.
4. *Amanzi* can run on a single processor on a laptop or thousands of cores on
   large supercomputers
5. Users can specify simulation input to *Amanzi* in a simple
   human-readable XML file that can be hand generated or created 
   using the ASCEM user platform, Akuna_.
6. Boundary conditions and materials can be specified using geometric
   regions defined on the computational domain. Regions may be
   combined together using Boolean operations to create other complex
   regions.
7. *Amanzi* supports checkpointing and restart to allow for simulations
   to be paused and resumed.
8. It supports output in an efficient format called
   `XDMF <http://www.xdmf.org>`_ that is recognized by two open-source
   high performance visualization tools,
   `ViSiT <http://wci.llnl.gov/codes/visit>`_ and
   `Paraview <http://www.paraview.org>`_.
9. *Amanzi* also allows observations at certain points or along
   specified planes to be written out to data files in support of
   uncertainty quantification and sensitivity analysis.
10. *Amanzi* supports multiple process models and more process
    capabilities are being added constantly. The process capabilities
    of *Amanzi* are detailed in the section called
    :ref:`Amanzi Process Capabilities <Amanzi-Process-Models>`.

.. _Akuna : http://akuna.labworks.org
