# MPC pre-processor instructions for Abaqus Analyses

This pre-processor package creates an "interface_properties.txt" file that is used by the uexternaldb user subroutine in Abaqus.
The interface properties file contains two key pieces of information for each MPC:
1. The node-to-element map that links the shell nodes on an interface with their attached elements.
1. The cross-section dimensions: d (or h), bf, tf, tw; 
where d is the total depth, bf is the flange width, tf is the flange thickness, and tw is the web thickness.

The function to use is:
```Python
run_mpc_preprocessor(input_file, sections_file)
```
 where `input_file` is the path to an input file, and `sections_file` contains the dimensions for each interface.
 This function `run_mpc_preprocessor` outputs a file "interface_props.txt" that needs to be input to the uexternaldb subroutine.
 
 ## Usage
 
 ### Input file
 The input file corresponds to the Abaqus job of interest.
 The user MPC needs to be specified in the standard format (e.g., with the jtype corresponding to the beam node, and all the nodes specified).
 This format is described in the MPC report.
 
 ### Sections file
 The sections file contains the cross-section dimensions corresponding to each interface.
 Examples are provided in the testing/ directory.
 A sample sections file is:
 ```
*Isection
2245,2247
500.,300.,28.,14.5

*Isection
22459,22479
2500.,2300.,128.,114.5
```
 
- The `*Isection` specifies a section start
- The `2245,2247` are beam nodes corresponding to MPCs with the following dimensions
- d = 500., bf = 300., tf = 28., tw = 14.5
- The rest of the lines are similar but for different nodes and dimensions.
