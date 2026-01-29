# DTLreactingFoam-12

## General Information
A package for high fidelity simulations of laminar reacting flows in OpenFOAM-12 with low computational cost, incorporating both the detailed transport model (DTM) and the polynomial fit transport model (FTM) based on the principle of kinetic gas theory [1]. To enhance computational efficiency, it was integrated with the time-correlated thermophysical property calculation (coTHERM) method. This technique can significantly reduce the computational cost of numerical simulations using DTM/FTM in OpenFOAM-12 while preserving accuracy. Readers are referred to our paper for all validation data. Readers are also referred to https://github.com/danhnam11/DTLreactingFoam-10 and https://github.com/danhnam11/DTLreactingFoam-8 for DTLreactingFoam in OpenFOAM-10 and OpenFOAM-8, respectively.

## Installation
- The complete installation of the OpenFOAM-12 framework in a Linux operating system is required before installing this package, as it is designed for the Linux-based OpenFOAM-12 version. 
- Prepare a directory on your system, for example, _yourDirectory_:

		mkdir ~/OpenFOAM/yourDirectory/
		cd ~/OpenFOAM/yourDirectory/	
- Download source files using git: 

		git clone https://github.com/danhnam11/DTLreactingFoam-12.git

- Specify the path of the _src_ directory of this package to an environment variable named _LIB_DTL12_SRC_. Suppose the _DTLreactingFoam-12_ have downloaded into _yourDirectory_. Then the following commands should be executed to specify the path of the _src_:

		echo "export LIB_DTL12_SRC=~/OpenFOAM/yourDirectory/DTLreactingFoam-12/src/" >> ~/.bashrc
		source ~/.bashrc

- To compile the necessary libraries and solver, go to _DTLreactingFoam-12_ directory and run the _Allwmake_ script (it may take one hour to be finished):

		cd ~/OpenFOAM/yourDirectory/DTLreactingFoam-12/
		./Allwmake

- After successful compilation, the following libraries are saved at _$FOAM_USER_LIBBIN_ :

		libspecie.so
		libthermophysicalProperties.so
		libfluidThermophysicalModels.so
		libmulticomponentThermophysicalModels.so
		libchemistryModel.so
		libsolidThermo.so
		libthermophysicalTransportModel.so
		libfluidThermophysicalTransportModel.so
		libsolidThermophysicalTransportModels.so
		libfluidThermoThermophysicalTransportModels.so
		libfluidMulticomponentThermophysicalTransportModels.so
		libphaseFluidMulticomponentThermophysicalTransportModels.so
		libphaseFluidThermophysicalTransportModels.so
		libcoupledThermophysicalTransportModels.so
		libphaseSolidThermophysicalTransportModels.so
		libradiationModels.so
		libcombustionModels.so
		libfvConstraints.so
		libfvModels.so
		libinterRegionFvModels.so
		librotorDisk.so
		libspecieTransfer.so
		libmulticomponentFluid.so
		libDTLreactingFoam.so
  
- and the following executable programs are saved at _$FOAM_USER_APPBIN_ :

  		chemkinToFoam
  		DTMchemkinToFoam
		FTMchemkinToFoam


- These newly compiled libraries, modular solvers, and utilities are now ready for use.
- It is important to note that if a different solver/modular solver (i.e., program) relies on any of the aforementioned compiled libraries, its corresponding _options_ file, located in the _Make_ directory, must be updated accordingly. The solver/modular sovler should then be recompiled to prevent potential conflicts, such as segmentation faults. For reference, the _Make_ directory of the _DTLreactingFoam_ modular solver included in this package provides a convenient example.

- To remove all compiled libraries and solvers, go to _DTLreactingFoam-12_ directory and run the _Allwclean_ script:

		cd ~/OpenFOAM/yourDirectory/DTLreactingFoam-12/
		./Allwclean

## Using this package 
Upon completing the compilation process, the _DTLreactingFoam_ modular solver can be utilized by specifying its name (i.e., _DTLreactingFoam_) in the solver keyword of controlDict dicitonary. All important instructions for case setting are provided in _documentations_ directory.

## Tutorials
Several test cases (e.g., canonical cases in combustion) are available in the _tutorials_ directory.

	cd ~/OpenFOAM/yourDirectory/DTLreactingFoam-12/tutorials/

## Authors 
This package was developed at the Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering, Ulsan National Institute of Science and Technology (UNIST), Korea (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/). If you publish results obtained by using this package, please cite our paper as follows:
- D. N. Nguyen, J. H. Lee, C. S. Yoo, DTLreactingFoam: An efficient CFD tool for laminar reacting flow simulations using detailed chemistry and transport with time-correlated thermophysical properties, Computer Physics Communications 322 (2026) 110052 (https://doi.org/10.1016/j.cpc.2026.110052).

If you need help with installation or have any questions, feel free to reach out: 
- danhnam11@gmail.com or nam.nguyendanh@hust.edu.vn 

## Reference
[1] R. J. Kee, F. M. Rupley, E. Meeks, J. A. Miller, CHEMKIN-III: A FORTRAN chemical kinetics package for the analysis of gas-phase chemical and plasma kinetics, SAND96-8216 (1996).

