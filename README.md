# BF-MPC: Best-fit Constraint Equations for Coupling Mixed-dimension Simulation Models with Wide Flange Cross Sections

Component macro models are a method of reducing the computational effort of full continuum component finite element models, while maintaining their solution fidelity.
Macro models contain 1D beam-column (or beam) element domains and 2D/3D continuum elements.
Some coupling method is required to join these domains of differing kinematics - a best-fit multipoint constraint method (BF-MPC) is proposed in this repository for this purpose.
Details of the method can be found in reference [1] at the end of this README.

This repository mainly contains two things: (1) the BF-MPC user subroutine for Abaqus, and (2) a Python pre-processor for Abaqus models.
The user subroutine can be used in a valid Abaqus model to couple beam-column element and continuum element domains (i.e., shell or solid elements).
The pre-processor is useful to generate a file that is used by the user MPC subroutine.


**IMPORTANT:** The code in the main branch is currently undergoing modification to make it compatible with recent versions of Abaqus (2019+) and the new Intel Fortran compilers. A new release with the finalized code and example analyses will be tagged when the modifications are complete.

Notes:
- The source code for the BF-MPC method is provided for clarity regarding the method and its implementation.
The current implementation is not considered to be optimal in terms of minimizing the number of operations carried out within the routine nor in minimizing the memory used by the subroutine.
- The remainder of this README file provides instructions for how to obtain, install, and use the BF-MPC method in Abaqus analyses.
The method is general, and could be used in other finite element simulation platforms, however, this functionality has not yet been implemented.
- The BF-MPC method is currently only implemented for wide flange (i.e., I-shaped) cross sections.
The functionality can be extended for other functions given that the warping function can be computed at each continuum node.
- See reference [1] at the end for the theory of the coupling method and some applications.

## Installation

Command-line instructions are provided to obtain and use the BF-MPC subroutine and the Python pre-processor.

1. Clone the BF-MPC repository

With `git` installed, open your command-line tool of choice (e.g., Command Prompt, Terminal, etc.) and enter:
```
git clone git@github.com:ahartloper/bf-mpc.git
```

2. Locate the BF-MPC subroutine

The Abaqus MPC user subroutine for the BF-MPC method is contained within the file `bf-mpc/bf-mpc_subroutine.for`.
There are three supporting files that must be in the same directory: `jacobi_eigenvalue.f90` (used to solve eigenvalue problems), `mpc_modules.f90` (contains functions used by the subroutine) and `uexternaldb_subr.f90` (used to load the data file).

## Usage

### BF-MPC user subroutine

The following are general instructions to use the BF-MPC in an Abaqus model.
Examples are provided in the `examples` directory.

1. Specify the user subroutine.

The BF-MPC user subroutine needs to be included with the analysis.
In Abaqus/CAE:
```
Double-click on the associated Job > General tab > Locate the bf-mpc_subroutine.for file
```
It is recommended that the user subroutine and the supporting files are copied to the work directory of the model.

2. Add the `*MPC` keyword definitions

See the Abaqus keyword manual for more information on the `*MPC` keyword.
For each coupling, the following keyword needs to be specified:
```
*MPC, MODE=NODE, USER
<jtype>, <beam node id>, <continuum node id 1>, ..., <continuum node id 14>,
0, <continuum node id 15>, ...
0, <continuum node id 31>, ...
```
where `jtype` is the mode specification of the BF-MPC subroutine, `continuum node id` is the id number of a node in the continuum domain on the interface, and `beam node id` is the id number of the beam node on the interface.
The "0" at the start of the second and further lines is necessary because the Abaqus input file processor does not read the first entry.
**Note** `jtype` should be the `beam node id` as it is used to retrieve the section properties information.
Multiple lines may be necessary depending on the number of continuum nodes because the Abaqus input file processor allows a maximum of 16 entries per line (14 continuum nodes on the first line and 15 on each subsequent line).
This keyword needs to be repeated for each coupling in the model.
For example:
```
*MPC, user, mode=node
2245,2245,   2,   3,   6,   8,   9,  12,  62,  63,  64,  65,  66, 175, 176, 177,
0,    178, 179, 229, 230, 231, 232, 233, 342, 343, 344, 345, 346, 347, 348, 349,
0,    350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364
```

## Contributing

Bug fixes and contributions can be raised by opening a new issue in the BF-MPC repository.

## Authors

Code written and maintained by Alex Hartloper (alexander.hartloper@epfl.ch).

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.
The `jacobi_eigenvalue.f90` used in the project is provided under the GNU LGPL library by John Burkardt.

## Acknowledgments

- Dimitrios Lignos and Albano de Castro e Sousa for their invaluable input.

## References
[1] Hartloper A.R., de Castro e Sousa A., and Lignos D.G. (2022). "Best-fit constraint equations for coupling mixed-dimension simulation models with wide flange cross sections", Finite Elements in Analysis and Design. DOI: doi.org/10.1016/j.finel.2022.103782.
