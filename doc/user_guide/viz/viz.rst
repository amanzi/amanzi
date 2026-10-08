.. viz:

============
Introduction
============

Amanzi writes out visualization information using HDF5 and Xdmf.  HDF5 is a platform independent binary file format utilized by many high performance computing applications.  Xdmf uses XML to describe data, store light data internally, and locate heavy data in a separate HDF5 file.  Xdmf provides a standard interface for exchanging information between applications.  

Mesh information is written out to an HDF5 binary file.  Field data is written to a separate HDF5 file.  At each timestep that visualization data is written out, a corresponding Xdmf file is written linking to the mesh and field data available for that timestep.  A single Xdmf file ending with ``VisIt.xmf`` collects all of the individual timestep Xdmf files.  VisIt and ParaView are third party visualization programs that have built in data readers for Xdmf files that can read the individual timestep files or the time series collection in the ``VisIt.xmf`` file.

Following are descriptions for installing and using both VisIt and ParaView.
