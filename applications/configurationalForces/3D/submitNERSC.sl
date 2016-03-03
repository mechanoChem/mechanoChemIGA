#!/bin/bash -l
#SBATCH -p regular
#SBATCH -N 10
#SBATCH -t 36:00:00
#SBATCH -J greght_config
#SBATCH -o greght_config.o%j
 

#Cori has 32 cores per compute node
fileName=main

rm *.o

#Run
srun -n 320 ./$fileName -ts_monitor -snes_monitor -snes_converged_reason -ksp_converged_reason -ts_adapt_type none -snes_type newtontr -ts_max_snes_failures 500 -snes_max_it 20 -snes_max_funcs 50000 -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package superlu_dist

