#!/bin/sh
####  PBS preamble

#PBS -N configForce

#PBS -M greght@umich.edu
#PBS -m abe

# Change the number of cores (ppn=1), amount of memory, and walltime:
#PBS -l nodes=10:ppn=16:prisms:sandybridge,pmem=3800mb,walltime=48:00:00
#PBS -j oe
#PBS -V

# Change "example_flux" to the name of your Flux allocation:
#PBS -A prismsproject_fluxoe
#PBS -q fluxoe
#PBS -l qos=flux
#PBS -p -1000

####  End PBS preamble

#  Show list of CPUs you ran on, if you're running under PBS
if [ -n "$PBS_NODEFILE" ]; then cat $PBS_NODEFILE; fi

#  Change to the directory you submitted from
if [ -n "$PBS_O_WORKDIR" ]; then cd $PBS_O_WORKDIR; fi

#  Put your job commands here:

#make
fileName=main

rm *.o $fileName

#compile
#Sacado
mpicxx -c -O3 -fPIC  -Wall -Wwrite-strings -I/$PETSC_DIR/$PETSC_ARCH/include -I/$PETSC_DIR/include -I$PETIGA_DIR/$PETSC_ARCH/include -I$PETIGA_DIR/include -I$TRILINOS_DIR/include/ $fileName.c 
mpicxx -o $fileName $fileName.o -Wl,-rpath,$PETIGA_DIR/$PETSC_ARCH/lib -L$PETIGA_DIR/$PETSC_ARCH/lib -lpetiga -Wl,-rpath,$PETSC_DIR/$PETSC_ARCH/lib -L$PETSC_DIR/$PETSC_ARCH/lib  -lpetsc -Wl,-rpath,$PETSC_DIR/$PETSC_ARCH/lib $TRILINOS_DIR/lib/libteuchoscore.so

#run
module load mkl


#Run
#mpirun -np 63 ./$fileName -ts_monitor -snes_monitor -snes_converged_reason -ksp_converged_reason -ts_adapt_type none -snes_type newtontr -ts_max_snes_failures 500 -snes_max_it 100 -snes_max_funcs 50000 -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps
#mpirun -np 16 ./$fileName -ts_monitor -snes_monitor -snes_converged_reason -ksp_converged_reason -ts_adapt_type none -snes_type newtontr -ts_max_snes_failures 500 -snes_max_it 100 -snes_max_funcs 500000 -ksp_type fgmres -pc_type ksp -ksp_ksp_type gmres -ksp_pc_type jacobi -snes_mf -ksp_gmres_restart 1200 -ksp_gmres_preallocate -ksp_gmres_modifiedgramschmidt -ksp_max_it 5001
#mpirun -np 32 ./$fileName -ts_monitor -snes_monitor -snes_converged_reason -ksp_converged_reason -ts_adapt_type none -snes_type newtontr -ts_max_snes_failures 500 -snes_max_it 100 -snes_max_funcs 50000 -ksp_type fgmres -pc_type ksp -ksp_ksp_type gmres -ksp_pc_type jacobi -ksp_gmres_restart 1200 -ksp_gmres_preallocate -ksp_gmres_modifiedgramschmidt -ksp_max_it 3001

mpirun -np 160 ./$fileName -ts_monitor -snes_monitor -snes_converged_reason -ksp_converged_reason -ts_adapt_type none -snes_type newtontr -ts_max_snes_failures 500 -snes_max_it 100 -snes_max_funcs 500000 -ksp_max_it 10001