=======================================
XML Input Specification (Current)
=======================================

The ''current'' XML Input Specification describes the XML input file that should work with
a clean clone and build of Amanzi right now.  However, to get started I just (4/21/11) copied the one specification from Phase 1, so HPC Toolset developers need to review their sections and update as necessary.

.. contents:: **Table of Contents**


Creating a Parameter List
=======================================

The XML input file must be framed at the beginning and end by the following statements:


.. code-block:: xml

  <ParameterList name="Main">

  </ParameterList>

where the name can be anything. Here I have named the main list "Main".

We have several sublists that are read and interpreted in various parts of Amanzi. Currently, the sublists must be named 

* Mesh
* MPC
* State
* Chemistry
* Flow
* Transport

The the following sections we first give an example of each of the sublists, and explain the options therein.

1. Mesh
=======================================

Example:

.. code-block:: xml

  <ParameterList name="Mesh">
    <Parameter name="Mesh Class" type="string" value="Simple"/>
    <ParameterList name="Simple Mesh Parameters">
      <Parameter name="X_Min" type="double" value="0.0"/>
      <Parameter name="X_Max" type="double" value="1.0"/>
      <Parameter name="Y_Min" type="double" value="0.0"/>
      <Parameter name="Y_Max" type="double" value="1.0"/>
      <Parameter name="Z_Min" type="double" value="0.0"/>
      <Parameter name="Z_Max" type="double" value="1.0"/>
      <Parameter name="Numer of Cells in X" type="int" value="10"/> 
      <Parameter name="Numer of Cells in Y" type="int" value="10"/> 
      <Parameter name="Numer of Cells in Z" type="int" value="10"/>
      <Parameter name="Number of mesh blocks" type="int" value="2"/>
      <ParameterList name="Mesh block 1">
	<Parameter name="Z0" type="double" value="-0.1"/>
	<Parameter name="Z1" type="double" value="0.3"/>
      </ParameterList>
      <ParameterList name="Mesh block 2">
	<Parameter name="Z0" type="double" value="0.3"/>
	<Parameter name="Z1" type="double" value="1.1"/>	
      </ParameterList>
    </ParameterList>
    <ParameterList name="MOAB Mesh Parameters">
      <Parameter name="Exodus file name" type="string" value="test/moab_mesh.exo"/>
    </ParameterList>   
  </ParameterList>

The Mesh parameter list has one parameter named "Mesh Class" that selects which of two mesh classes is being used to create the mesh. The two options for the value of the parameter "Mesh Class" are "Simple" or "MOAB". 

1.1. Simple Mesh
---------------------------------------

In the case where value="Simple", the sublist named "Simple Mesh Parameters" triggers the creation of a simple serial 3D rectangular mesh. The user must specify the dimensions of the mesh, the number of mesh blocks which are layers in the mesh and are defined by each mesh block's minimum and maximum z-coordinates (the parameters named "Z0" and "Z1", respectively. Make sure that there are no gaps or overlaps in the definition of mesh blocks, or not all cells will be part of one unique mesh block.

The sublists that define the individual mesh blocks must be named "Mesh block 1", "Mesh block 2", and so forth. The mesh block IDs that are associated with these mesh blocks are numbered subsequently starting at zero, so that the ID associated with "Mesh block 1" is 0, the ID associated with "Mesh block 2" is 1, and so forth.

There will be six side sets in the mesh numbered 0 through 5 and are as follows:

 * Side set 0 is incident with the yz plane at x=X_Max.
 * Side set 1 is incident with the xz plane at y=Y_Max.
 * Side set 2 is incident with the yz plane at x=X_Min.
 * Side set 3 is incident with the xz plane at y=Y_Min.
 * Side set 4 is incident with the xy plane at z=Z_Min.
 * Side set 5 is incident with the xy plane at z=Z_Max.


1.2. MOAB Mesh
---------------------------------------

In the case where value="MOAB" there is just a single parameter, in the sublist named "MOAB Mesh Parameters". This parameter is named "Exodus file name" and is used to specify the Exodus file name that will be read to initialize the MOAB mesh.

2. MPC
=======================================

Example:

.. code-block:: xml

  <ParameterList name="MPC">
    <Parameter name="Start Time" type="double" value="0.0"/>
    <Parameter name="End Time" type="double" value="0.1"/>
    <Parameter name="End Cycle" type="int" value="10000"/>
    <Parameter name="Flow model" type="string" value="Darcy"/>
    <Parameter name="disable Flow_PK" type="string" value="no"/>
    <Parameter name="disable Transport_PK" type="string" value="no"/>
    <Parameter name="disable Chemistry_PK" type="string" value="yes"/>
    <Parameter name="Viz dump cycle frequency" type="int" value="10"/>
    <Parameter name="Viz dump time frequency" type="double" value="0.05"/>
    <ParameterList name="CGNS">
      <Parameter name="File name" type="string" value="test1.cgns"/>
    </ParameterList>
  </ParameterList> 

In the MPC parameter list, the user specifies the following parameters:

 * "Start Time" the start time of the simulation
 * "End Time" the end time of the simulation
 * "End Cycle" the end cycle of the simulation 
 * "Flow model" specifies the choice of flow model.  The choices are currently `"Darcy"` for saturated flow and  `"Richards"` for unsaturated flow.
 * "disable Flow_PK" is used to disable flow in the simulation. In this case the user should specify a mesh block wise constant darcy flow field in the State namelist.
 * "disable Transport_PK" is used to disable transport in the simulation.
 * "disable Chemistry_PK" is used to disable chemistry in the simulation.
 * "Viz dump cycle frequency" is used to generate visualization dumps every so many cycles.
 * "Viz dump time frequency" is used to generate visualization dumps every so many time increments.

The sublist named "CGNS" is used to specify the filename for the CGNS visualization dumps. 


3. State
=======================================

Example:

.. code-block:: xml

  <ParameterList name="State">
    <Parameter name="Number of mesh blocks" type="int" value="2"/>
    <Parameter name="Number of component concentrations" type="int" value="3"/>
    <Parameter name="Constant water saturation" type="double" value="1.0"/>
    <Parameter name="Constant water density" type="double" value="1000.0"/>
    <Parameter name="Constant viscosity" type="double" value="1.0"/>
    <Parameter name="Gravity x" type="double" value="0.0"/>
    <Parameter name="Gravity y" type="double" value="0.0"/>
    <Parameter name="Gravity z" type="double" value="0.0"/>

    <ParameterList name="Mesh block 1"> 
      <Parameter name="Mesh block ID" type="int" value="0"/>
      <Parameter name="Constant porosity" type="double" value="0.2"/>
      <Parameter name="Constant permeability" type="double" value="10.0"/>
      <Parameter name="Constant component concentration 0" type="double" value="0.0"/>
      <Parameter name="Constant component concentration 1" type="double" value="0.1"/>
      <Parameter name="Constant component concentration 2" type="double" value="0.2"/>
      <Parameter name="Constant Darcy flux x" type="double" value="1.0"/>
      <Parameter name="Constant Darcy flux y" type="double" value="0"/>
      <Parameter name="Constant Darcy flux z" type="double" value="0"/> -->
    </ParameterList>

    <ParameterList name="Mesh block 2"> 
      <Parameter name="Mesh block ID" type="int" value="1"/>
      <Parameter name="Constant porosity" type="double" value="0.2"/>
      <Parameter name="Constant permeability" type="double" value="1.0"/>
      <Parameter name="Constant component concentration 0" type="double" value="0.1"/>
      <Parameter name="Constant component concentration 1" type="double" value="0.3"/>
      <Parameter name="Constant component concentration 2" type="double" value="0.5"/>
      <Parameter name="Constant Darcy flux x" type="double" value="1.0"/>
      <Parameter name="Constant Darcy flux y" type="double" value="0"/> 
      <Parameter name="Constant Darcy flux z" type="double" value="0"/> 
    </ParameterList>    

  </ParameterList>


The following parameters can be set in the State parameter list:

 * "Number of mesh blocks" defines the number of mesh blocks in the mesh. This number should be consistent with the definition of the mesh.
 * "Number of component concentrations" defines the number of total component concentrations that the state object will keep track of.
 * "Constant water saturation" sets the water saturation to a constant on the entire mesh.
 * "Constant water density" sets the water density to a constant on the entire mesh.
 * "Constant viscosity" sets the viscosity to a constant on the entire mesh.
 * "Gravity x", "Gravity y", and "Gravity z" are used to specify the gravity vector.

Additionally, some values can be initialized by mesh block. The parameter lists that specify these values per mesh block must be named "Mesh block 1", "Mesh block 2" and so forth. The number of the mesh block sub parameter list must equal the parameter "Number of mesh blocks" above.

In each mesh block sub parameter list the user can specify the following parameters:

 * "Mesh block ID" is the ID number of the specified mesh block. This must equal a mesh block ID specified in the mesh  file  
 * "Constant porosity" defines the porosity on the entire mesh block.
 * "Constant permeability" defines the permeability on the entire mesh block.
 * "Constant component concentration 0", "Constant component concentration 1" and so forth, define the total component concentration for the mesh block. There must be as many definitions for the total component concentration as there are components (see "Number of component concentrations" above).
 * "Constant Darcy flux x", "Constant Darcy flux y", and "Constant Darcy flux z" define a constant Darcy flux on the mesh block. These parameters are optional and will not be used if flow is on.

4. Transport
=======================================

Example:


.. code-block:: xml

  <ParameterList name="Transport">
    <Parameter name="CFL" type="double" value="0.5"/>   
    <!-- debug and developers options -->
    <Parameter name="enable internal tests" type="string" value="no"/>   
    <Parameter name="internal tests tolerance" type="double" value="1e-6"/>   
    <Parameter name="verbosity level" type="int" value="0"/>  
    <Parameter name="maximal time step" type="double" value="10"/>
    <!-- end of debug and developers options -->
    
    <ParameterList name="Transport BCs">
      <Parameter name="number of BCs" type="int" value="1"/>
      <ParameterList name="BC 0">
	<Parameter name="Side set ID" type="int" value="3"/>
	<Parameter name="Type" type="string" value="Constant"/>
	<Parameter name="Component 0" type="double" value="1.0"/>
	<Parameter name="Component 1" type="double" value="0.6"/>
	<Parameter name="Component 2" type="double" value="0.2"/>
      </ParameterList>  
      <ParameterList name="BC 1"> 
        <Parameter name="Side set ID" type="int" value="50001"/> 
        <Parameter name="Type" type="string" value="Constant"/> 
        <Parameter name="Component 0" type="double" value="1.0"/> 
        <Parameter name="Component 2" type="double" value="1.0e-4"/>     
      </ParameterList>
    </ParameterList>
  </ParameterList>

In the Transport sublist the user specifies the following parameters: 
 
 * "CFL" is the Courant–Friedrichs–Lewy number. It must be strictly bigger than 0 and less or equal to 1. Default value is 1. 
 * "enable internal tests" turns on/off build-in tests. This option is useful for code development; therefore its default value is "no". 
 * "verbosity level" sets up the volume of information printed out by the transport. It must be any non-negative integer. This option is useful for code development; therefore, its default value is 0.
 * "internal tests tolerance" is the relative tolerance for internal tests. This is the developers option. Default value is 1e-6.
 * "maximal time step" overwrites the calculated time step. This is the developers option.  
 	 
The boundary conditions sublist consists of a few similar sublists related to boundary side sets. The number of these sublists can be both bigger or smaller than the number of defined side sets. Each of the sublists may contain only a few components. The other components will be automatically set to zero. Note that the boundary conditions have to be set up mathematically only on influx boundaries. If it is not done, the default boundary condition is always zero.   

 * "number of BCs" is the total number of boundary conditions (i.e. subsequent sublists). 
 * "Side set ID" is the side set id in the mesh model. 
 * "Type" specifies the boundary condition. At the moment only constant boundary conditions are available. Put a ticket if you need a different type of boundary condition. 
 * "Component X" specified the value of component X on this boundary. 


5. Chemistry
=======================================

Example:

.. code-block:: xml

  <ParameterList name="Chemistry">
    <Parameter name="Thermodynamic Database Format" type="string" value="simple" />
    <Parameter name="Thermodynamic Database File" type="string" value="fbasin-uo2-5-component.bgd" />
    <Parameter name="Verbosity" type="int" value="1" />
    <Parameter name="Activity Model" type="string" value="debye-huckel" />
    <Parameter name="Tolerance" type="double" value="1.5e-12"/>
    <Parameter name="Max Time Step (s)" type="double" value="86400.0"/>
    <Parameter name="Maximum Newton Iterations" type="int" value="150"/>
    <Parameter name="Using sorption" type="string" value="yes"/>
    <Parameter name="Free ion concentrations provided" type="string" value="yes"/>
    <ParameterList name="Initial Conditions">
      <Parameter name="Number of minerals" type="int" value="3"/>
      <Parameter name="Number of ion exchange sites" type="int" value="0"/>
      <Parameter name="Number of sorption sites" type="int" value="0"/>
      <Parameter name="Number of mesh blocks" type="int" value="1"/>
      <ParameterList name="Mesh block 1"> 
        <Parameter name="Mesh block ID" type="int" value="0"/>
        <ParameterList name="Free Ion Species">
	  <Parameter name="Free Ion Species 0" type="double" value="4.36476e-16"/>  <!-- Al+++ -->
	  <Parameter name="Free Ion Species 1" type="double" value="3.16772e-08"/>  <!-- H+ -->
	  <Parameter name="Free Ion Species 2" type="double" value="1.00000e-06"/>  <!-- HPO4-2 -->
	  <Parameter name="Free Ion Species 3" type="double" value="1.87000e-04"/>  <!-- SiO2(aq) -->
	  <Parameter name="Free Ion Species 4" type="double" value="1.84374e-20"/>  <!-- UO2++ -->
        </ParameterList>
        <ParameterList name="Minerals">
          <Parameter name="Mineral 0" type="double" value="0.15"/>  <!-- Kaolinite -->
          <Parameter name="Mineral 1" type="double" value="0.21"/>  <!-- Quartz -->
          <Parameter name="Mineral 2" type="double" value="0.0"/>   <!-- (UO2)3(PO4)2.4H2O -->
        </ParameterList>
      </ParameterList>
    </ParameterList>
  </ParameterList>


'''Note: all chemistry names and values are case sensitive.'''

+------------------------------------+---------------+------------------+-----------------------------+----------------------------------------------------------------------------------------+
|  Parameter name                    |  Type         |  Default Value   | Optional Values             | Purpose                                                                                |
+====================================+===============+==================+=============================+========================================================================================+
| `"Thermodynamic Database Format"`  | string        | `"simple`"       | `"simple"`                  | set the format of the database                                                         |
+------------------------------------+---------------+------------------+-----------------------------+----------------------------------------------------------------------------------------+
| `"Thermodynamic Database File"`    | string        | `"dummy.dbs"`    |  ---                        | path name to the chemistry database file, relative to the program execution directory. |
+------------------------------------+---------------+------------------+-----------------------------+----------------------------------------------------------------------------------------+

The following parameters are optional in the Chemistry parameter list:

+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+
|  Parameter name                       | Type          | Default Value    | Optional Values             | Purpose                                                                             |
+=======================================+===============+==================+=============================+=====================================================================================+
| `"Verbosity"`                         | int           | 0                | 0, 1, 2, 3, 4, 5, 6, ...    | set the verbosity level of chemistry: 0=silent, 1=terse warnings, 2=verbose details,|
|                                       |               |                  |                             |  3=debug, 4=debug beaker, 5=debug database file, ....                               | 
+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+
| `"Activity Model"`                    | string        | `"unit`"         | `"unit"`, `"debye-huckel"`  | set the model used for activity corrections                                         |
+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+
| `"Tolerance"`                         | double        | 1.0e-12          |  ---                        | set the tolerance for newton iterations within chemistry                            |
+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+
| `"Maximum Newton Iterations"`         | int           | 200              | ---                         | set the maximum number of newton iterations for chemistry.                          |
+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+
| `"Max Time Step (s)"`                 | double        | 9.9e9            | ---                         | set the maximum time step allowed for chemistry.                                    |
+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+
| `"Using sorption"`                    | string        | `"no"`           | `"yes"`                     | Tells the chemistry module whether to allocate memory for sorption.                 |
+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+
| `"Free ion concentrations provided"`  | string        | `"no"`           | `"yes"`                     | Tells chemistry that in initial guess for free ion concentrations is provided in    |
|                                       |               |                  |                             | the xml file.                                                                       |
+---------------------------------------+---------------+------------------+-----------------------------+-------------------------------------------------------------------------------------+

The initial condition list must have a `"Mesh Block"` parameter list for each mesh block, mesh block numbering should correspond to the other mesh block lists. Each mesh block list will have parameter lists for the non-zero elements of the chemistry. Valid parameter list names are: `"Free Ion Species"` `"Minerals"` `"Ion Exchange Sites"` `"Sorption Sites"`.

Each initial condition list should contain a parameter name constructed like `"Type #"` where `"Type"` is `"Mineral"`, `"Free Ion Species"`, `"Ion Exchange Site"` `"Sorption Site"` and `"#"` in the integer identifier, starting with zero.  

'''Note: it is recommended that you include an xml comment with the species or mineral name after each initial condition. The xml parser expects every instance of `"--"` to mark a comment, so species names with negative charges should be written as `"SO4-2"` rather than `"SO4--"`.'''

5.1 Chemistry Thermodynamic data file formats 
-------------------------------------------------

5.1.1 XML Database format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

not yet implemented

5.1.2 Simple Database format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Importing thermodynamic data into the chemistry module using the `"simple"` (file extension `"bgd"`) format requires the user to explicitly specify all the species and reactions for the problem. There is no basis switching or automatic species and reaction selection. Below is an example of a `"simple"` database file for a five component uranium problem with mineral dissolution and surface complexation:

::

 <Primary Species
 # name               ; debye-huckel a0 ; charge ; GMW     

 Al+++                ;   9.0 ;   3.0 ;  26.9815
 H+                   ;   9.0 ;   1.0 ;   1.0079
 HPO4--               ;   4.0 ;  -2.0 ;  95.9793
 SiO2(aq)             ;   3.0 ;   0.0 ;  60.0843
 UO2++                ;   4.5 ;   2.0 ;  270.028

 <Aqueous Equilibrium Complexes
 # name               =  coeff primary_name  coeff primary_name  ; log10(Keq) 25C ; debye-huckel a0 ; charge ; GMW      

 OH-                  =  1.0 H2O  -1.0 H+                ;    13.9951 ;   3.5 ;  -1.0 ;  17.0073 
 AlOH++               =  1.0 H2O  1.0 Al+++  -1.0 H+     ;     4.9571 ;   4.5 ;   2.0 ;  43.9889 
 Al(OH)2+             =  2.0 H2O  1.0 Al+++  -2.0 H+     ;    10.5945 ;   4.0 ;   1.0 ;  60.9962 
 Al(OH)3(aq)          =  3.0 H2O  1.0 Al+++  -3.0 H+     ;    16.1577 ;   3.0 ;   0.0 ;  78.0034 
 Al(OH)4-             =  4.0 H2O  1.0 Al+++  -4.0 H+     ;    22.8833 ;   4.0 ;  -1.0 ;  95.0107 
 UO2OH+               =  1.0 H2O  -1.0 H+  1.0 UO2++     ;     5.2073 ;   4.0 ;   1.0 ;  287.035 
 UO2(OH)2(aq)         =  2.0 H2O  -2.0 H+  1.0 UO2++     ;    10.3146 ;   3.0 ;   0.0 ;  304.042 
 UO2(OH)3-            =  3.0 H2O  -3.0 H+  1.0 UO2++     ;    19.2218 ;   4.0 ;  -1.0 ;   321.05 
 UO2(OH)4--           =  4.0 H2O  -4.0 H+  1.0 UO2++     ;    33.0291 ;   4.0 ;  -2.0 ;  338.057 
 (UO2)2OH+++          =  1.0 H2O  -1.0 H+  2.0 UO2++     ;     2.7072 ;   5.0 ;   3.0 ;  557.063 
 (UO2)2(OH)2++        =  2.0 H2O  -2.0 H+  2.0 UO2++     ;     5.6346 ;   4.5 ;   2.0 ;   574.07 
 (UO2)3(OH)4++        =  4.0 H2O  -4.0 H+  3.0 UO2++     ;     11.929 ;   4.5 ;   2.0 ;  878.112 
 (UO2)3(OH)5+         =  5.0 H2O  -5.0 H+  3.0 UO2++     ;    15.5862 ;   4.0 ;   1.0 ;   895.12 
 (UO2)3(OH)7-         =  7.0 H2O  -7.0 H+  3.0 UO2++     ;    31.0508 ;   4.0 ;  -1.0 ;  929.135 
 (UO2)4(OH)7+         =  7.0 H2O  -7.0 H+  4.0 UO2++     ;    21.9508 ;   4.0 ;   1.0 ;  1199.16 
 UO2(H2PO4)(H3PO4)+   =  3.0 H+  2.0 HPO4--  1.0 UO2++   ;   -22.7537 ;   4.0 ;   1.0 ;   465.01 
 UO2(H2PO4)2(aq)      =  2.0 H+  2.0 HPO4--  1.0 UO2++   ;   -21.7437 ;   3.0 ;   0.0 ;  464.002 
 UO2HPO4(aq)          =  1.0 HPO4--  1.0 UO2++           ;    -8.4398 ;   3.0 ;   0.0 ;  366.007 
 UO2H2PO4+            =  1.0 H+  1.0 HPO4--  1.0 UO2++   ;   -11.6719 ;   4.0 ;   1.0 ;  367.015 
 UO2H3PO4++           =  2.0 H+  1.0 HPO4--  1.0 UO2++   ;   -11.3119 ;   4.5 ;   2.0 ;  368.023 
 UO2PO4-              =  -1.0 H+  1.0 HPO4--  1.0 UO2++  ;    -2.0798 ;   4.0 ;  -1.0 ;  364.999 

 <Minerals
 # name               =  coeff primary_name  coeff primary_name  ; log10(Keq) 25C ; GMW      ; molar volume [cm^2/mol] ; SSA [m^2/g] 

 Kaolinite            =  5.00 H2O  2.00 Al+++  -6.00 H+  2.00 SiO2(aq)  ;     6.8101 ;   258.16 ;    99.52 ;   1.0 
 Quartz               =  1.00 SiO2(aq)  ;    -3.9993 ;  60.0843 ;   22.688 ;   1.0 
 (UO2)3(PO4)2.4H2O    =  4.00 H2O  -2.00 H+  2.00 HPO4--  3.00 UO2++  ;   -27.0349 ;  1072.09 ;    500.0 ;   1.0 

 <Mineral Kinetics
 # name               ; TST ; log10_rate_constant double     moles_m2_sec 

 Kaolinite            ; TST ; log10_rate_constant    -16.699 moles_m2_sec 
 Quartz               ; TST ; log10_rate_constant      -18.0 moles_m2_sec 
 (UO2)3(PO4)2.4H2O    ; TST ; log10_rate_constant      -10.0 moles_m2_sec 

 <Surface Complex Sites
 # name               ; surface_density

 >FeOH                ; 6.3600E-03
 >AlOH                ; 6.3600E-03
 >SiOH                ; 6.3600E-03

 <Surface Complexes
 # name               =  coeff surface site  coeff primary_name  ; log10(Keq) 25C ; charge 

 >SiOUO3H3++          =  1.0 >SiOH  1.0 H2O  1.0 UO2++  ;       5.18 ;   2.0 
 >SiOUO3H2+           =  1.0 >SiOH  1.0 H2O  -1.0 H+  1.0 UO2++  ;       5.18 ;   1.0 
 >SiOUO3H             =  1.0 >SiOH  1.0 H2O  -2.0 H+  1.0 UO2++  ;       5.18 ;   0.0 
 >SiOUO3-             =  1.0 >SiOH  1.0 H2O  -3.0 H+  1.0 UO2++  ;      12.35 ;  -1.0 
 >SiOUO2(OH)2-        =  1.0 >SiOH  2.0 H2O  -3.0 H+  1.0 UO2++  ;      12.35 ;  -1.0 
 >FeOHUO3             =  1.0 >FeOH  1.0 H2O  -2.0 H+  1.0 UO2++  ;       3.05 ;   0.0 
 >FeOHUO2++           =  1.0 >FeOH  1.0 UO2++  ;      -6.63 ;   2.0 
 >AlOUO2+             =  1.0 >AlOH  -1.0 H+  1.0 UO2++  ;      -3.13 ;   1.0 

Note the following about this format:

 * Any line starting with a `"#"` or space character is a comment. 

 * Data is separated into sections, where each section of the file is starts with a line containing `"<Section Name"`. The valid section names are: `"Primary Species"`, `"Aqueous Equilibrium Complexes"`, `"Minerals"`, `"Mineral Kinetics"`, `"General Kinetics"`, `"Ion Exchange Sites"`, `"Ion Exchange Complexes"`, `"Surface Complex Sites"`, `"Surface Complexes"`. The less than character, `"<"`, should be the first character on the line and there is no space between the character and the section name.

 * Sections should be ordered so that the primary species, minerals, and exchange sites come before any reactions using those species.

 * Each line within a section is a semi-colon delimited

 * A primary species line contains the primary species must contain:

   ::

     # name               ; debye-huckel a0 ; charge ; GMW [grams/mole]    
     Al+++                ;   9.0 ;   3.0 ;  26.9815

 * An aqueous equilibrium complex line contains a reaction and data for the reaction on a single line:

   ::

     # name               =  coeff primary_name  coeff primary_name ... ; log10(Keq) 25C ; debye-huckel a0 ; charge ; GMW [grams/mole]     
     OH-                  =  1.0 H2O  -1.0 H+  ;    13.9951 ;   3.5 ;  -1.0 ;  17.0073 

   The reaction is written as product species = reactants.... The coefficient of the product aqueous complex is assumed to be 1.0, and the reactants must be primary species. The equilibrium constant is for a fixed temperature of 25C.

 * Minerals and other complexes follow the same convention as aqueous equilibrium complexes, with additional data as needed.

   ::

     <Minerals
     # name               =  coeff primary_name  coeff primary_name ... ; log10(Keq) 25C ; GMW      ; molar volume [cm^2/mol] ; SSA [m^2/g] 

     <Surface Complexes
     # name               =  coeff surface site  coeff primary_name ... ; log10(Keq) 25C ; charge 

     These are all minerals present in the system during the simulation, including those that may precipitate later. They are used for calculating saturation states, but not equilibrium or kinetic calculations.

 * The mineral kinetics section lists the name of a mineral found in the mineral section, the type of rate law, and rate parameters for that law.

   :: 

     <Mineral Kinetics
     # name               ; TST ; log10_rate_constant double     moles_m2_sec ; primary_name coeff ....
 
   Currently only the `"TST"` rate law is implemented. The keywords "log10_rate_constant" and "moles_m2_sec" must be present in the line, but no unit conversion are currently preformed. The    modifying primary species terms follow the rate constant, along with their exponent coefficients.

 * Surface complex sites are listed by name and surface density:

   ::

     <Surface Complex Sites 
     # name               ; surface_density [moles sites / m^2 mineral]


6. Flow
=======================================

Example:

.. code-block:: xml

  <ParameterList name="Flow">

    <Parameter name="Max Iterations" type="int" value="100"/>
    <Parameter name="Error Tolerance" type="double" value="1.0e-13"/>
    <Parameter name="Nonlinear Solver" type="string" value="NLK"/>
    <Parameter name="Preconditioner Update Frequency" type="int" value="1"/>
    
    <ParameterList name="Flow BC">
      <Parameter name="number of BCs" type="int" value="2"/>
      <ParameterList name="BC00">
	<Parameter name="Type" type="string" value="Darcy Constant"/>
	<Parameter name="BC value" type="double" value="-1.0" />
	<Parameter name="Side set ID" type="int" value="3" />
      </ParameterList>  
      <ParameterList name="BC01">
	<Parameter name="Type" type="string" value="Pressure Constant"/>
	<Parameter name="BC value" type="double" value="0.0" />
	<Parameter name="Side set ID" type="int" value="1" />
      </ParameterList>  
     </ParameterList>
  </ParameterList>


The following parameters must be specified in the  Flow parameter list:

 * "Max Iterations" defines the maximum number of iterations the flow solver is allowed to take.
 * "Error Tolerance" defines the error tolerance to which the flow solver will attempt to solve the flow equation.

When the Richards flow model is used the following parameters must also be specified:

 * "Nonlinear Solver" defines the choice of nonlinear solver.  The valid choices are `"JFNK"` for Jacobian-free Newton-Krylov and `"NLK"` for the Nonlinear Krylov method.
 * "Preconditioner Update Frequency" sets how frequently the preconditioner will be recomputed during the iterative nonlinear solve.  With the value 1 it will be recomputed every iteration, with 2 every other iteration, and so forth.

The sub list named "Flow BC" contains the definition of the boundary conditions for the flow equations. The number of these conditions that are specified is defined by the parameter named "number of BCs". The boundary condition sublists must be named "BC00", "BC01" and so forth. Each of these boundary condition sublists must contain the following paramters:

 * "Type" defines the boundary condition type, allowed values are "Darcy Constant", "Pressure Constant", "Static Head", or "No Flow".
 * "BC value" is the value that should be applied, its interpretation depends on the parameter "Type" above.
 * "Side set ID" is the ID number of the side set in the mesh where the boundary condition should be applied.

The default boundary condition is "No Flow". It is applied to all boundary faces that are in side sets that do not have a corresponding boundary condition sublist.
