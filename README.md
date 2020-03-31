# Amanzi

Master:
[![Build Status](https://travis-ci.org/amanzi/amanzi.svg?branch=master)](https://travis-ci.org/amanzi/amanzi)
Amanzi-0.88:
[![Build Status](https://travis-ci.org/amanzi/amanzi.svg?branch=amanzi-0.88)](https://travis-ci.org/amanzi/amanzi)


### BACKGROUND

Amanzi provides a flexible and extensible parallel flow and reactive transport simulation capability for environmental applications. It includes general polyhedral mesh infrastructure, which leverages MSTK, advanced discretizations of process models, including traditional finite volume schemes, mimetic finite differences, and nonlinear finite volumes. In addition, it provides advanced nonlinear solvers, such as Nonlinear Krylov Acceleration (NKA) and Anderson Acceleration, and leverages Trilinos-ML and Hypre Algebraic Multigrid for scalable solvers. The reaction of contaminants and minerals carried by the flow through the surrounding rock and soil is modeled through a general software interface called Alquimia that allows Amanzi to interface with a variety of powerful geochemistry engines including PFLOTRAN and CrunchFlow. The code is parallel and leverages open-source parallel frameworks such as Trilinos, PETSc. Amanzi has been used to model contaminant migration at various DOE waste sites (e.g., Nevada National Security Site, and Savannah River), and is generally applicable to groundwater contaminant migration under partially saturated, nonisothermal conditions and its interaction with surface water. 

The multiphysics framework in Amanzi is called Arcos, and it provides modelers with the flexibility they need to creatively decompose complex problems and explore a variety of mixed-dimensional model configurations to develop understanding and make predictions of environmental systems. In particular, Arcos provides flexibility for hierarchical weak and strong coupling of processes with subcycling of mixed dimensions. This capability in conjuction with Amanzi's powerful mesh infrastructure, which supports the splitting and subsetting of meshes,  enables creative conceptual modeling. Applications in Amanzi include, coupling flow and transport on discrete-fracture-networks (DFNs) and the background matrix, while applications in the [Advanced Terrestrial Simulator (ATS)](https://amanzi.github.io/ats) include, integrated hydrology coupling surface and subsurface processes; and an intermediate scale thermal hydrology model of polygonal tundra based on one-dimensional columns coupled to the two-dimensional surface.

### BUILDING/INSTALLING AND RUNNING AMANZI

See the [INSTALL.md](INSTALL.md) file in this directory.


### COPYRIGHT AND LICENSE

The copyright and license are contained in the top-level file
[COPYRIGHT.md](COPYRIGHT.md).  The copyright is held jointly by the DOE laboratories that
are contributing code to the development of Amanzi. 


