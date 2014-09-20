module load mkl
fileName=mechanoChemoND
#fileName=nonConvexMechanicsND

rm *.o $fileName
#Compile according to AD type
#Sacado
mpicxx -c -O3 -fPIC  -Wall -Wwrite-strings -Wno-unused-variable -Wno-unused-value  -Wno-uninitialized -Wno-strict-aliasing -Wno-unknown-pragmas -I/$PETSC_DIR/$PETSC_ARCH/include -I/$PETSC_DIR/include -I$PETIGA_DIR/$PETSC_ARCH/include -I$PETIGA_DIR/include -I$TRILINOS_DIR/include/ $fileName.c 
mpicxx -o $fileName $fileName.o -Wl,-rpath,$PETIGA_DIR/$PETSC_ARCH/lib -L$PETIGA_DIR/$PETSC_ARCH/lib -lpetiga -Wl,-rpath,$PETSC_DIR/$PETSC_ARCH/lib -L$PETSC_DIR/$PETSC_ARCH/lib  -lpetsc -Wl,-rpath,$PETSC_DIR/$PETSC_ARCH/lib /opt/software/numerics/trilinos-11.4.3-Source/install/lib//libteuchoscomm.dylib	 /opt/software/numerics/trilinos-11.4.3-Source/install/lib//libteuchoscore.dylib /opt/software/numerics/trilinos-11.4.3-Source/install/lib//libteuchosnumerics.dylib /opt/software/numerics/trilinos-11.4.3-Source/install/lib//libteuchosparameterlist.dylib /opt/software/numerics/trilinos-11.4.3-Source/install/lib//libteuchosremainder.dylib

#Adept
#mpicxx -c -O3 -fPIC -DADEPT_RECORDING_PAUSABLE -Wall -Wwrite-strings -Wno-unused-variable -Wno-unused-value  -Wno-uninitialized -Wno-strict-aliasing -Wno-unknown-pragmas -pipe -I/$PETSC_DIR/$PETSC_ARCH/include -I/$PETSC_DIR/include -I$PETIGA_DIR/$PETSC_ARCH/include -I$PETIGA_DIR/include -I/opt/software/numerics/adept-1.0/include $fileName.c
#mpicxx -o $fileName $fileName.o -Wl,-rpath,$PETIGA_DIR/$PETSC_ARCH/lib -L$PETIGA_DIR/$PETSC_ARCH/lib -lpetiga -Wl,-rpath,$PETSC_DIR/$PETSC_ARCH/lib -L$PETSC_DIR/$PETSC_ARCH/lib  -lpetsc -Wl,-rpath,$PETSC_DIR/$PETSC_ARCH/lib /opt/software/numerics/adept-1.0/lib/libadept.a

time mpiexec -np 4 ./$fileName -file_prefix "test" -N 80 -dt 2.5e-5 -ch_monitor -ts_monitor -snes_monitor -snes_converged_reason -log_summary -ts_max_snes_failures 200 -snes_max_it 200  -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps 
#-snes_type newtontr
#-pc_type svd -pc_svd_monitor
#-ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps
#-snes_linesearch_monitor

#time mpiexec -np 4 ./$fileName -file_prefix "test" -N 10 -ch_output -ch_monitor -ts_monitor -snes_monitor -snes_converged_reason -dt 0.01 -ts_max_snes_failures 200 -snes_max_it 200 -snes_type newtontr 

#mpiexec -np 8 ./nonConvexMechanics2D -file_prefix "test" -N 400 -ch_output -ch_monitor -ts_monitor -snes_monitor -snes_converged_reason -dt 0.01 -ts_max_snes_failures 200 -snes_max_it 200 -snes_type newtontr -ksp_type preonly -pc_type ilu -pc_factor_mat_solver_package mumps

#mpiexec -np 1 ./$fileName -file_prefix "test" -N 50 -ch_output -ch_monitor -ts_monitor -snes_monitor -snes_converged_reason -ksp_type gmres -ksp_gmres_restart 1200 -pc_type hypre -log_summary -ksp_converged_reason 
#-ts_max_snes_failures -1
#-snes_linesearch_type basic
#-ksp_gmres_restart 300 -pc_type hypre
