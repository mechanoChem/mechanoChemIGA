<B>mechanoChem: Modeling solid mechanics, chemistry, and their interactions</B> <br>
=======================================================================

Developed by the Computational Physics Group at the University of Michigan. 
http://www.umich.edu/~compphys/index.html <br>

<B>List of contributors:</B> <br>
Greg Teichert (Lead Developer) <br>
Shiva Rudraraju <br>
Krishna Garikipati <br>

<B>Code documentation:</B> https://goo.gl/VNAt2Y <br>

<B>Overview</B> <br>
=======================================================================
The mechanoChem code is an isogeometric analysis based code used to solve the partial differential equations describing solid mechanics (including gradient elasticity) and chemistry (including the Cahn-Hilliard phase field model). It is built on the PetIGA [https://bitbucket.org/dalcinl/petiga/] and PETSc [https://www.mcs.anl.gov/petsc/] libraries, and it uses the automatic differentiation capabilities of the Sacado package from the Trilino library [https://trilinos.org/packages/sacado/]. <br>


<B>Version information</B>
=======================================================================
This is version 0.1, the intial release of the code. <br>


<B>License</B>
=======================================================================
GNU Lesser General Public License (LGPL). Please see the file LICENSE for details. <br>


<B>Acknowledgements</B>
=======================================================================
This code has been developed under the support of the following: <br>

NSF DMREF grant: DMR1436154 "DMREF: Integrated Computational Framework for Designing Dynamically Controlled Alloy-Oxide Heterostructures" <br>
NSF CDI Type I grant: CHE1027729 "Meta-Codes for Computational Kinetics" <br>
DOE BES, Division of Materials Sciences and Engineering: Award #DE-SC0008637 that funds the PRedictive Integrated Structural Materials Science (PRISMS) Center at University of Michigan <br>


<B>Referencing this code</B>
=======================================================================
If you write a paper using results obtained with the help of this code,  please consider citing one or more of the following: <br>

"A variational treatment of material configurations with application to interface motion and microstructural evolution" (under review) <br>
G. Teichert, S. Rudraraju, K. Garikipati  <br>

<pre>
\@article{Teichert2016a,
  Title                    = {A variational treatment of material configurations with application to interface motion and microstructural evolution},
  Author                   = {G. Teichert and S. Rudraraju and K. Garikipati},
  Journal                  = {ArXiv e-prints},
  archivePrefix            = "arXiv",
  eprint                   = {1608.05355},
  Year                     = {2016},
  Month                    = {Jul},
  Url                      = {https://arxiv.org/abs/1608.05355}
}
</pre>

"A comparison of Redlich-Kister polynomial and cubic spline representations of the chemical potential in phase field computations" (under review) <br>
G. Teichert, H. Gunda, S. Rudraraju, A. Natarajan, B. Puchala, K. Garikipati, A. Van der Ven <br>

<pre>
\@article{Teichert2016b,
  Title                    = {A comparison of Redlich-Kister polynomial and cubic spline representations of the chemical potential in phase field computations},
  Author                   = {G. Teichert and N. S. H. Gunda and S. Rudraraju and A. R. Natarajan and B. Puchala and A. Van der Ven and K. Garikipati},
  Journal                  = {ArXiv e-prints},
  archivePrefix            = "arXiv",
  eprint                   = {1609.00704},
  Year                     = {2016},
  Month                    = {Aug},
  Url                      = {https://arxiv.org/abs/1609.00704}
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

To run the code (replace "nproc" with the number of processors to be use), with some recommended flags: <br>

	$ mpirun -np nproc ./main -ts_monitor -snes_type newtontr -ksp_type fgmres
