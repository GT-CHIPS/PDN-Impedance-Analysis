PDN Impedance Analysis Tool for Heteregeneous Integration 
--------------------------------------------------------------------------------------------------------------------
Copyright (c) 2018 Hakki Mert Torun
Power Delivery Network (PDN) Impedance Analysis for Heteregenous Integration. 
Components in PDN analysis include:
 - PCB P/G plane,
 - C4/u-bump bump array with customizable ground/signal/power pattern
 - TSV array with customizable ground/signal/power pattern
 - via array with customizable ground/signal/power pattern
 
All components are parameterized with respect to their material and geometrical dimensions. 
The code is tested on Matlab R2017b.
The papers from which this tool is based on are provided inside the scripts in appropriate places.
Example output of the code is provided as a plot.
This material is based on work supported by DARPA CHIPS project under award N00014-17-1-2950.
For questions and queries, please contact: htorun3@gatech.edu

--------------------------------------------------------------------------------------------------------------------

RUNNING THE CODE:
--------------------------------------------------------------------------------------------------------------------
The project contains 6 MATLAB scripts:
--"calc_PDN_impedance.m"
This the main code where parameterization is included.
T-Matrices of individual components are calculated and cascaded to each other.
Proper matrix reductions are performed for parallel excitation of elements in TSV/bump/via arrays.

--"calc_T_TSV.m"
Multiconductor RLGC matrices of coupled TSV arrays with depletion capacitance.
Please see the comments inside script for parameterization and ground/signal/power assignment.

--"calc_T_bump.m"
T-Matrix calculation for coupled C4/via/u-bumps. Bumps are approximated as cylindrical structures.
Please see the comments inside script for parameterization and ground/signal/power assignment.

--"calc_Z_plane.m"
Z-Matrix of a power/ground plane pair. Arbitrary port definitions can be assigned to any location
on the planes for detailed analysis. Please see the comments inside the script for details.

--"calc_T_grid.m"
T-Matrix of a power/ground grid on interposer. Arbitrary port definitions can be assigned to any location
on the grid for detailed analysis. Please see the comments inside the script for details.
