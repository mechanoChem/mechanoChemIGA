<B>mechanoChemIGA: Modeling solid mechanics, chemistry, and their interactions</B> <br>
=======================================================================

Developed by the Computational Physics Group at the University of Michigan. 
http://www.umich.edu/~compphys/index.html <br>

<B>List of contributors:</B> <br>
Greg Teichert (Lead Developer) <br>
Shiva Rudraraju <br>
Koki Sagiyama <br>
Krishna Garikipati <br>

<B>[Code documentation](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChem/master/doxygen/html/index.html)</B><br>

<B>Overview</B> <br>
=======================================================================
The mechanoChemIGA code is an isogeometric analysis based code used to solve the partial differential equations describing solid mechanics (including gradient elasticity) and chemistry (including the Cahn-Hilliard phase field model). It is built on the PetIGA [https://bitbucket.org/dalcinl/petiga/] and PETSc [https://www.mcs.anl.gov/petsc/] libraries, and it uses the automatic differentiation capabilities of the Sacado package from the Trilino library [https://trilinos.org/packages/sacado/]. <br>


<B>Version information</B>
=======================================================================
This is version 0.2.1. <br>


<B>License</B>
=======================================================================
GNU Lesser General Public License (LGPL). Please see the file LICENSE for details. Note that the functions IGAElementNextFormFunction and IGAComputeProjectionFunction in the file src/output.cc, as well as all functions 
in the file src/petigasnes_mod.h were derived from the PetIGA/src/petigasnes.c source code
in the PetIGA library [https://bitbucket.org/dalcinl/petiga/]. Accordingly,
we include the license/copyright notice for the PetIGA library here in the file LICENSE_PetIGA 
to apply to the above functions.<br>


<B>Acknowledgements</B>
=======================================================================
This code has been developed under the support of the following: <br>

Toyota Research Institute, Award #849910 "Computational framework for data-driven, predictive, multi-scale and multi-physics modeling of battery materials" <br>
NSF DMREF grant: DMR1436154 "DMREF: Integrated Computational Framework for Designing Dynamically Controlled Alloy-Oxide Heterostructures" <br>
NSF CDI Type I grant: CHE1027729 "Meta-Codes for Computational Kinetics" <br>
DOE BES, Division of Materials Sciences and Engineering: Award #DE-SC0008637 that funds the PRedictive Integrated Structural Materials Science (PRISMS) Center at University of Michigan <br>


<B>Referencing this code</B>
=======================================================================
If you write a paper using results obtained with the help of this code,  please consider citing one or more of the following: <br>

"A variational treatment of material configurations with application to interface motion and microstructural evolution" (Journal of the Mechanics and Physics of Solids) <br>
G. Teichert, S. Rudraraju, K. Garikipati  <br>

<pre>
\@article{Teichert2016a,
	title	= "A variational treatment of material configurations with application to interface motion and microstructural evolution ",
	journal = "Journal of the Mechanics and Physics of Solids ",
	volume 	= "99",
	pages 	= "338 - 356",
	year 	= "2017",
	issn 	= "0022-5096",
	doi 	= "https://doi.org/10.1016/j.jmps.2016.11.008",
	url 	= "http://www.sciencedirect.com/science/article/pii/S0022509616305221",
	author 	= "Gregory H. Teichert and Shiva Rudraraju and Krishna Garikipati",
}
</pre>

"A comparison of Redlich-Kister polynomial and cubic spline representations of the chemical potential in phase field computations" (Computational Materials Science) <br>
G. Teichert, H. Gunda, S. Rudraraju, A. Natarajan, B. Puchala, K. Garikipati, A. Van der Ven <br>

<pre>
\@article{Teichert2016b,
	title	= "A comparison of Redlich-Kister polynomial and cubic spline representations of the chemical potential in phase field computations ",
	journal = "Computational Materials Science ",
	volume 	= "128",
	pages 	= "127 - 139",
	year 	= "2017",
	issn 	= "0927-0256",
	doi 	= "https://doi.org/10.1016/j.commatsci.2016.11.024",
	url 	= "http://www.sciencedirect.com/science/article/pii/S0927025616305754",
	author 	= "Gregory H. Teichert and N.S. Harsha Gunda and Shiva Rudraraju and Anirudh Raju Natarajan and Brian Puchala and Krishna Garikipati and Anton Van der Ven",
}
</pre>

"Mechano-chemical spinodal decomposition: A phenomenological theory of phase transformations in multi-component, crystalline solids" (Nature npj Computational Materials) <br>
S. Rudraraju, A. Van der Ven, K. Garikipati  <br>

<pre>
\@article{Rudraraju2016,
  Title                    = {Phenomenological treatment of chemo-mechanical spinodal decomposition},
  Author                   = {S. Rudraraju and A. Van der Ven and K. Garikipati},
  Journal                  = {npj Computational Materials},
  Year                     = {2016},
  Volume                   = {2},
  Doi                      = {10.1038/npjcompumats.2016.12}
}
</pre>


<B>Installation</B>
=======================================================================
1) Install PETSc: <br>

-Download and extract PETSc source code. <br>
-Quick installation as follows (the symbol $ denotes the command prompt): <br>

	$ ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-fblaslapack --download-mpich
	$ make all test

-Set the appropriate PETSC_DIR and PETSC_ARCH environment variables, e.g. <br>

	$ export PETSC_DIR=/path/to/petsc-3.6.4
	$ export PETSC_ARch=arch-linux2-c-debug

Note that this also installs mpich. It is possible to use an existing version of mpi by including a flag to its directory. If you will be using the local PETSc installation of mpich, set the following:

	$ alias mpirun=$PETSC_DIR/$PETSC_ARCH/bin/mpirun

If you would like to use a parallel direct solver, we recommend using superlu_dist. Add the following flags when configuring PETSc: --download-metis --download-parmetis --download-superlu_dist <br>

Download: http://www.mcs.anl.gov/petsc/download/index.html <br>
Installation instructions: http://www.mcs.anl.gov/petsc/documentation/installation.html <br>


2) Install PetIGA: <br>

-Download and extract or clone the PetIGA source code. <br>

-Enter the PetIGA top directory and install using make: <br>

	$ make all
	$ make test

-Export path to PetIGA directory: <br>

	$ export PETIGA_DIR=/path/to/petiga/

Download: https://bitbucket.org/dalcinl/petiga/downloads (or clone the bitbucket repository) <br>
Installation instructions: https://bitbucket.org/dalcinl/petiga/ <br>


3) Install CMake: <br>

Download: https://cmake.org/download/ <br>


4) Install the Sacado package from Trilinos (version 11.10.2 recommended) <br>

-Download and extract or clone the Trilinos source code. <br>
-Simple installation of Sacado as follows (from the Trilinos top directory), with the desire /path/to/trilinos/installation/: <br>

	$ mkdir build
	$ cd build
	$ cmake -DTrilinos_ENABLE_Sacado=ON -DTrilinos_ENABLE_Teuchos=OFF -DCMAKE_INSTALL_PREFIX=/path/to/trilinos/installation/ ../
	$ make install

-Export path to Trilinos: <br>

	$ export TRILINOS_DIR=/path/to/trilinos/installation/

Download: http://trilinos.org/oldsite/download/ <br>
Installation instructions: https://trilinos.org/docs/files/TrilinosBuildReference.html <br>

5) Install igakit (required to convert binary output files to .vtk files)

Download and installation instructions: https://bitbucket.org/dalcinl/igakit

<B>Usage</B>
=======================================================================
To run the example initial boundary value problems, navigate to the desired example folder. Create a makefile using cmake: <br>

	$ cmake CMakeLists.txt

To compile the code (default is debug mode): <br>

	$ make

To switch the compilation to release mode: <br>

	$ make release

To switch the compilation to debug mode: <br>

	$ make debug

To run the code (replace "nproc" with the number of processors to be use; superlu_dist is used by default if available, gmres is used otherwise): <br>

	$ mpirun -np nproc ./main

To use use gmres instead of superlu_dist: <br>

	$ mpirun -np nproc ./main -ksp_type gmres -pc_type none

The ouput files created by the code are .dat binary files and a fieldInfo.txt file. These files can be converted to .vtk files and visualized using tools such as VisIt or ParaView using the igakit package (see step 5 of the installation instructions). To do this file conversion, run the initBounValProbs/writeVTKFile.py script from the directory containing the output files. For example, if the .dat and fieldInfo.txt ouput files were located in the initBounValProb/nonGradientMechanics/3D folder, the following commands would create the .vtk files: <br>

	$ cd initBounValProb/nonGradientMechanics/3D
	$ python ../../writeVTKFiles.py


