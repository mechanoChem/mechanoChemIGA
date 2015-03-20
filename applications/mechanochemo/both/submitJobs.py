import shutil
import subprocess
import time, datetime, os

#
cluster_name="prismsproject_fluxoe"
#
GRID=[100,200]
DT=[1.0e-6, 1.0e-7]
#BC={0:"SHEAR",1:"FREE",2:"FIXED"}
BC={1:"FREE",2:"FIXED"}
#FLUX={0:"ALL",1:"TOPBOTTOM",2:"SIDES",3:"QUENCH"}
#FLUX={1:"TOPBOTTOM"}
FLUX={3:"QUENCH"}
LAMBDA=[1.0e-5, 1.0e-6, 1.0e-7]

#gmres
#runOptions="-ts_monitor -snes_monitor -snes_converged_reason -pc_type hypre -ksp_type gmres -ksp_gmres_restart 600 -ts_adapt_type none -ts_max_snes_failures 100000 -snes_max_it 200 -snes_type newtontr -snes_stol 1.0e-10 -snes_trtol 0.0"

#mumps
runOptions="-ts_monitor -snes_monitor -snes_converged_reason -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -ts_adapt_type none -ts_max_snes_failures 100000 -snes_max_it 200 -snes_type newtontr -snes_stol 1.0e-10 -snes_trtol 0.0"

#create today directory
today = datetime.date.today()
todaystr = today.isoformat()
dirPath0=os.path.join(os.getcwd(),"resultsQuenchRes2")
dirPath=os.path.join(os.getcwd(),"resultsQuenchRes2",todaystr)
if not os.path.exists(dirPath0):
    os.mkdir(dirPath0)
if not os.path.exists(dirPath):
    os.mkdir(dirPath)

#iterate
for dt, N, bc, flux, lam in [(dt, N, bc, flux, lam) for dt in DT for N in GRID for bc in BC.keys() for flux in FLUX.keys() for lam in LAMBDA]:
    fileName='N'+str(N)+"l"+str(lam)+"t"+str(dt)+"F"+FLUX[flux]+"B"+BC[bc]
    shutil.copyfile ("./main.c", fileName+".c")
    print fileName
    varDefs=["-D dtVal="+str(dt),"-D NVal="+str(N),"-D bcVAL="+str(bc),"-D FLUX="+str(flux),"-D El="+str(lam)]
    
    #compile code
    cmd = ["mpicxx"] + varDefs + ["-c", "-O3", "-fPIC","-I/nfs/mcfs_home/rudraa/Public/petsc/petsc/shared-optimizedDevDec182014/include", "-I/nfs/mcfs_home/rudraa/Public/petsc/petsc/include",  "-I/nfs/mcfs_home/rudraa/Public/petiga/PetIGA/shared-optimizedDevDec182014/include", "-I/nfs/mcfs_home/rudraa/Public/petiga/PetIGA/include", "-I/nfs/mcfs_home/rudraa/Public/trilinos/ver11.12.1/install/include/", fileName+".c"];  
    p = subprocess.Popen(cmd);  
    p.wait();  

    #link code
    cmd = ["mpicxx", "-o", fileName, fileName+".o","-L/nfs/mcfs_home/rudraa/Public/petiga/PetIGA/shared-optimizedDevDec182014/lib", "-lpetiga", "-L/nfs/mcfs_home/rudraa/Public/petsc/petsc/shared-optimizedDevDec182014/lib", "-lpetsc","/nfs/mcfs_home/rudraa/Public/trilinos/ver11.12.1/install/lib/libteuchoscore.so"];  
    p = subprocess.Popen(cmd);  
    p.wait();  

    #copy code,exe to dirPath
    dirPath2=os.path.join(dirPath, fileName)
    if not os.path.exists(dirPath2):
        os.mkdir(dirPath2)
    try:
        shutil.copy (fileName+".c", dirPath2+"/"+fileName+".c")
        shutil.copy (fileName, dirPath2+"/"+fileName)
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)
        continue
 #write options to options file
    f = open(dirPath2+"/options.log", 'wt')
    #f.write(str(compileDefs))
    #f.write("\n")
    f.write(str(varDefs))
    f.write("\nGrid:"+str(N)+" dt:"+str(dt)+" bc:"+str(bc))

    #PBS
    # Open a pipe to the qsub command.
    #output, input = popen2('qsub')    
    
    # Customize your options here
    job_name = fileName
    walltime = "48:00:00"
    if (dt<=1.0e-6):
         walltime = "48:00:00"
    processors = "nodes=1:ppn=16:prisms:sandybridge,pmem=3800mb,qos=flux"
    command = "time mpirun -np 8 "+dirPath2+"/"+fileName+" "+ runOptions
    if (N>=300):
        processors = "nodes=4:ppn=16:prisms:sandybridge,pmem=3800mb,qos=flux"
        command = "time mpirun -np 64 "+dirPath2+"/"+fileName+" "+ runOptions
    elif (N>=200):
        command = "time mpirun -np 16 "+dirPath2+"/"+fileName+" "+ runOptions
    print command
    job_string = """#!/bin/bash
#PBS -A %s
#PBS -N %s
#PBS -l walltime=%s
#PBS -l %s
#PBS -q fluxoe
#PBS -n
#PBS -V
#PBS -m n
module load intel-comp
module load mkl
cd %s
%s""" % (cluster_name, job_name, walltime, processors, dirPath2, command)
    # Send job_string to qsub
    fpbs = open(dirPath2+"/pbsJob", 'wt')
    fpbs.write(job_string)
    fpbs.close()
    cmd = ["qsub","-k","oe",dirPath2+"/pbsJob"];  
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE);  
    pbsOut=p.communicate()[0]
    p.wait();
    #write to options log
    f.write("\n\n")
    f.write(job_string)
    f.write("\n\n")
    f.write(pbsOut)
    f.close()
    print pbsOut
    time.sleep(0.1)
