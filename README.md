FPVCentralFoam-3.0.1
==================

## Description

This code realizes a steady laminar flamelet approach for turbulent non-premixed combustion.
The solver is based on ''rhoCentralFoam'', i.e. it is density-based (For high speed flow, Mach > 1), compressible and runs with both LES and RAS turbulence.
 
The theory is mainly taken from the work of N. Peters and is based on the view of a turbulent flame as an ensemble of laminar flamelets.
The calculation of these flamelets is a one-dimensional problem and can be done in a pre-processing step.
Integration using a presumed beta-Probability Density Function (PDF) accounts for the interaction between turbulent fluctuations and chemistry.
The results of the pre-processing procedure are stored in tables which are accessed during the simulation.
Values of interest, such as species mass fraction or enthalpy, are looked-up and connected to the flow using three parameters - the mixture fraction, its variance and the progress parameter.
In doing so, the expensive solution of chemical mechanisms during run-time can be avoided and the run-time thus reduces significantly.

## Installation

This version works with OpenFOAM-3.0.1

* Prepare a directory on your system, e.g.:  

  `mkdir ~/OpenFOAM/FPVCentralFoam/`

* Download FPVCentralFoam using git:

  `git clone https://github.com/weixian001/FPVCentralFoam-v3.0.1.git ~/OpenFOAM/FPVCentralFoam/`

* Set an environment variable to the FPVCentralFoam src folder:

  `export LIB_FPVCentralFoam_SRC=~/OpenFOAM/FPVCentralFoam/src/`

* Execute `./Allwmake`

## Notes

This solver was based on the work done by Prof Pfitzner, FlameletFoam created at the Universität der Bundeswehr München, Thermodynamics Institute (Prof. Pfitzner). http://sourceforge.net/projects/openfoam-extend/files/OpenFOAM_Workshops/OFW8_2013_Jeju/Fri/Track3/HagenMuller-OFW8.tar/download

The FPVFoam was developed by a NanyangCFD team at Nanyang Technological University, Singapore lead by Prof Chan. Main contributor is Wei Xian Lim (weixian001@e.ntu.edu.sg).

# FPVFoam-v3.0.1
