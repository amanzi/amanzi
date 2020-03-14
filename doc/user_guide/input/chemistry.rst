Chemistry in Amanzi with Alquimia
----------------------------------------------

Amanzi makes available through the Alquimia interface geochemical
capabilities available in widely used geochemical codes. 

The geochemical capabilities of the reactive transport code PFLOTRAN
(http://ees.lanl.gov/pflotran/ and https://bitbucket.org/pflotran/pflotran-dev/wiki/Home) 
are available within Amanzi through the Alquimia interface. The specifications required 
for a PFLOTRAN input file are available in the PFLOTRAN User Manual 
(http://www.pflotran.org/docs/user_manual.pdf).

The geochemical capabilities of the reactive transport code CRUNCHFLOW 
(http://www.csteefel.com/CrunchFlowIntroduction.html) are available through the 
Alquimia interface. The User Manual for CrunchFlow can be found here: 
http://www.csteefel.com/CrunchFlowManual.pdf.

For convenience, the user is kindly referred to the examples available in the Amanzi documentation 
(see :doc:`../benchmarking/chemistry/index`), and the related input files available in the repository,
for additional examples of input files for the different geochemical engines.

About Alquimia:
+++++++++++++++

Alquimia is a biogeochemistry API and wrapper library that provides a 
unified interface to the biogeochemistry routines from high quality 
reactive transport simulators such as PFLOTRAN or CrunchFlow, 
allowing any subsurface flow and transport simulator to access a range of functionality.

The official repository for Alquimia API documentation and reference API is 
https://bitbucket.org/berkeleylab/alquimia. The alquimia API is intended to be open and implemented 
by anyone. The reference implementation developed by Lawrence Berkeley National Laboratory is open
source, licensed under a BSD License.
