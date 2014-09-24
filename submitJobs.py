import shutil
import subprocess
from popen2 import popen2
import time, datetime, os
from subprocess import call

FS=[1]
Explicit=[1]
GRID=[100,200,400]
DT=[1.0e-5,1.0e-6,1.0e-7]
BC={0:"SHEAR",1:"FREE",2:"FIXED"}
FLUX={0:"ALL",1:"TOPBOTTOM",2:"SIDES",3:"QUENCH"}
LAMBDA=[1.0,10.0,100.0]
#parameters
El=1.0
compileDefs=["-D DIM=2", "-D EXPLICIT", "-D finiteStrain", "-D ADSacado", "-D numVars=27"]

#create today directory
today = datetime.date.today()
todaystr = today.isoformat()
dirPath0=os.path.join(os.getcwd(),"results")
dirPath=os.path.join(os.getcwd(),"results",todaystr)
if not os.path.exists(dirPath0):
    os.mkdir(dirPath0)
if not os.path.exists(dirPath):
    os.mkdir(dirPath)

#iterate
for dt, N, bc, flux, lam in [(dt, N, bc, flux, lam) for dt in DT for N in GRID for bc in BC.keys() for flux in FLUX.keys() for lam in LAMBDA]:
    fileName="Flux"+FLUX[flux]+"-Eg"+str(lam)+"-bc"+BC[bc]+'-N'+str(N)+'-dt'+str(dt)
    shutil.copyfile ("/scratch/prismsproject_flux/rudraa/mechanochemo/mechanochemocode/mechanoChemoND.c", fileName+".c")
    print fileName
    varDefs=["-D dtVAL="+str(dt),"-D NVAL="+str(N),"-D bcVAL="+str(bc),"-D FLUX="+str(flux),"-D EgVAL="+str(lam)]
    
    #compile code
    cmd = ["mpicxx"] + compileDefs + varDefs + ["-c", "-O3", "-fPIC","-I/nfs/mcfs_home/rudraa/Public/petsc/petsc/shared-optimized/include", "-I/nfs/mcfs_home/rudraa/Public/petsc/petsc/include", "-I/nfs/mcfs_home/rudraa/Public/petiga/PetIGA/include", "-I/nfs/mcfs_home/rudraa/Public/trilinos/ver10.12.2/install/include", fileName+".c"];  
    p = subprocess.Popen(cmd);  
    p.wait();  

    #link code
    cmd = ["mpicxx", "-o", fileName, fileName+".o","-L/nfs/mcfs_home/rudraa/Public/petiga/PetIGA/shared-optimized/lib", "-lpetiga", "-L/nfs/mcfs_home/rudraa/Public/petsc/petsc/shared-optimized/lib", "-lpetsc","/nfs/mcfs_home/rudraa/Public/trilinos/ver10.12.2/install/lib/libteuchos.so"];  
    p = subprocess.Popen(cmd);  
    p.wait();  

    #copy code,exe to dirPath
    dirPath2=os.path.join(dirPath, fileName)
    if not os.path.exists(dirPath2):
        os.mkdir(dirPath2)
    shutil.copy (fileName+".c", dirPath2+"/"+fileName+".c")
    shutil.copy (fileName, dirPath2+"/"+fileName)
    
    #write options to options file
    f = open(dirPath2+"/options.log", 'wt')
    f.write(str(compileDefs))
    f.write("\n")
    f.write(str(varDefs))
    f.write("\nGrid:"+str(N)+" dt:"+str(dt)+" bc:"+str(bc))
    f.close()

    #PBS
    # Open a pipe to the qsub command.
    output, input = popen2('qsub')    
    
    # Customize your options here
    job_name = fileName
    walltime = "24:00:00"
    if (dt<1.0e-6):
         walltime = "48:00:00"
    processors = "nodes=1:ppn=16,pmem=3800mb,qos=flux"
    command = "time mpiexec -np 16 "+dirPath2+"/"+fileName+" -ts_monitor -snes_monitor -snes_converged_reason -log_summary -ts_max_snes_failures 200 -snes_max_it 200  -ksp_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -snes_linesearch_type basic"
    job_string = """#!/bin/bash
#PBS -N %s
#PBS -l walltime=%s
#PBS -l %s
cd %s
%s""" % (job_name, walltime, processors, dirPath2, command)
    # Send job_string to qsub
    input.write(job_string)
    input.close()
    # Print your job and the system response to the screen
    print job_string
    print output.read()
    time.sleep(0.1)


'''
