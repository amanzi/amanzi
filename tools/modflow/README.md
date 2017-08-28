## MODFLOW to Amanzi ##
### *A conversion script using PyLaGriT* ###
------------------------------------

In its current form, this program converts Modflow materials and elevation files into an AVS-UCD mesh, with facesets in the binary Exodus file format.

### 0. Requirements ###
It is necessary that PyLaGriT, along with the Exodus library, be installed in order for this script to work.
LaGriT with Exodus can be built by following the instructions on the [LaGriT GitHub repository](https://github.com/lanl/LaGriT).

To configure PyLaGriT, visit the [PyLaGriT documentation](https://lanl.github.io/LaGriT/gettingstarted.html).

### 1. Usage ###
The script is intended to be both easy to use and completely configurable. The core of user interaction lies in the `PARAMS` file - this file contains all of the different variables that one should need in order to get the desired output. It should not be necessary to edit the source to accomplish what you want - ideally, all mesh configurations can be set up through the parameters file.

You may use and edit the provided `PARAMS` file or create your own. To pass in a parameter file, run the line

    python MODFLOW_to_Amanzi.py -p name_of_param_file

If no parameter file is passed, the script will attempt to locate and load the file `PARAMS`.
If ParaView is installed (and configured within PyLaGriT), you may view the final mesh by passing `-v` or `--view` as an argument.

### 2. Configuring parameters ###
The included `PARAMS` file lists all of the variables required for this version of the script. If a custom parameter file is passed in and some variables are missing, the script will prompt you to enter the values. Entering nothing into the prompt and pressing `ENTER` will cause the script to use the default values for those variables. The default values can be viewed/edited by looking in the source at the function `verify_params()`.

### 3. Contact Information ###
This script is still under active development and has not been tested for the general case. To report issues or suggestions, please email [livingston@lanl.gov](mailto:livingston@lanl.gov).

### 4. Supported variables ###

| Variable 		   | Type   | Description     |
|-------------------|--------|-----------------|
| `INPUT_BNDS`  		| `str`  | Materials file input  |
| `INPUT_ELEV`  		| `str`  | *(Optional)* Elevation file  |
| `EXPORT_NAME`	 	| `str`  | AVS outfile name  |
| `EXO_FILE`		 	| `str`  | Exodus   |
| `NCOLS` 			| `int`  | Number of columns  |
| `NROWS`  			| `int`  | Number of rows  |
| `DX`  				| `real` | Cell spacing X  |
| `DY`  				| `real` | Cell spacing Y  |
| `HEIGHT`  			| `real` | Z-thickness of final mesh |
| `MINX`  			| `real` | Lower-left hand corner X  |
| `MINY`  			| `real` | Lower-left hand corner Y  |
| `RMMAT_EDGE` 		| `int`  | IMAT material to remove |
| `MAXMAT` 			| `int`  | Max IBND material |
| `IBND_NOFLOW` 		| `int`  | Materials value from `INPUT_BNDS` to declare as the noflow region |
| `IBND_HEAD` 		| `int`  | Materials value from `INPUT_BNDS` to declare as the head region |
| `IBND_EDGE` 		| `int`  | Materials value from `INPUT_BNDS` to declare as the edge region |
| `IBND_HALITE` 		| `int`  | Materials value from `INPUT_BNDS` to declare as the halite region |
| `IMAT_ACTIVE` 		| `int`  | Positive integer corresponding to `IBND_ACTIVE` value |
| `IMAT_NOFLOW` 		| `int`  | Positive integer corresponding to `IBND_NOFLOW` value |
| `IMAT_HEAD` 		| `int`  | Positive integer corresponding to `IBND_HEAD` value |
| `IMAT_EDGE` 		| `int`  | Positive integer corresponding to `IBND_EDGE` value |
| `IMAT_HALITE` 		| `int`  | Positive integer corresponding to `IBND_HALITE` value |
| `FS_BOTTOM` 		| `int`  | set the bottom faces as the `n`-th faceset  |
| `FS_TOP` 			| `int`  | set the top faces as the `n`-th faceset  |
| `FS_EAST` 			| `int`  | set the east faces as the `n`-th faceset  |
| `FS_NORTH` 			| `int`  | set the north faces as the `n`-th faceset  |
| `FS_WEST` 			| `int`  | set the west faces as the `n`-th faceset  |
| `FS_SOUTH` 			| `int`  | set the south faces as the `n`-th faceset  |
| `FS_HALITE` 		| `int`  | set the halite faces as the `n`-th faceset  |
| `FS_HEAD` 			| `int`  | set the head faces as the `n`-th faceset  |
| `FS_NOFLOW` 		| `int`  | set the noflow boundary faces as the `n`-th faceset  |