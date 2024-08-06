.. _Amanzi Input:

Amanzi Input Description
========================

To run *Amanzi*, the minimum input required is a description of 
the simulation domain and parameters specified in an XML file 
using the schema described in the 
:ref:`Amanzi XML Schema <Amanzi XML Schema>`
page. Additionally, *Amanzi* might need an input mesh (for 
unstructured mesh based runs where the mesh is not internally 
generated) as well as any other auxilliary files (like color function 
file for specifying a region of the domain - see 
:ref:`Color Function Region  <Color Function Region>`).

.. toctree::
   :maxdepth: 2
   
   units.rst
   input_schema.rst
   mesh.rst
   chemistry.rst
