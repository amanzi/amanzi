Multiphysics Models
-------------------

*Amanzi* allows the user to couple physical models in a number of ways. 
For a single continuum, the following multiphysics models can be created:

  * Flow and transport 
  * Flow and reactive transport
  * Thermal flow 

Two spatial continua are supported by *Amanzi*:

  * Fracture network
  * Porous matrix

A fracture network could be meshed by the user and given to *Amanzi* via an 
input Exodus file. Alternatevely, a fracture network can be automatically 
extracted from a volumetric mesh using labeled mesh sets, see the Exodus manual
for more detail.

*Amanzi* allows the user to couple physical processes defineid on two spatial continua.
Example include:

  * Flow and transport in a coupled matrix-fracture network system
  * Flow and reactive transport in a coupled matrix-fracture network system


