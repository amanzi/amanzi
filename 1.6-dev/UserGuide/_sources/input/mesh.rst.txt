Simulation Mesh
=================

When running *Amanzi* in the structured mode, the mesh is generated internally
and the user is not required to specify it explicitly.

When running *Amanzi* in the unstructured mode, it is possible to ask
*Amanzi* to internally generate the mesh or read a mesh that is generated
by an external tool like `LaGriT <http://lagrit.lanl.gov>`_ and written out
in the `Exodus II <http://sourceforge.net/projects/exodusii/>`_ format. 

If the mesh is internally generated and the simulation is being run on
multiple cores, the mesh will automatically be distributed over the
multiple processors.

If the mesh is externally generated and the simulation is on multiple
cores, then the user has two options for distributing the mesh across
processors.  In the first option, the complete mesh exists in a single
*Exodus II* file which is read in, partitioned and distributed to the
various processors. The input file must have an extension of *.exo*
for this option. The second option is to prepartition the mesh on disk
using the *splitexo* script under the *bin* directory of the *Amanzi*
tree into files and letting *Amanzi* read each one in one processor
(*Nemesis I* files are *Exodus II* files with extra information
allowing them to be used in a parallel setting). To specify a
prepartitioned file in the *Nemesis I* format, users must specify the
filename with a *.exo* or *.par* extension. No processor or partition
number needs to be specified even though the files on disk have a
processor suffix. Thus, a serial mesh might be specified as
*mymesh.exo* while a mesh partitioned into four parts, *mymesh.par.0*,
*mymesh.par.1*, *mymesh.par.2*, *mymesh.par.2* would be specified as
*mymesh.exo* or *mymesh.par*.
