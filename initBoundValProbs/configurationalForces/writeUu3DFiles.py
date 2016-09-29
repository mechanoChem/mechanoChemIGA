import glob
from igakit.io import PetIGA,VTK
from igakit.nurbs import NURBS
import numpy as np
import itertools
import scipy.io
import os, sys

data={};
cwdir=os.getcwd();
geom = PetIGA().read('mesh.dat')
inc=1;
for filename1, filename2 in zip(glob.glob('outU*.dat'), glob.glob('outE*.dat')):
    #print filename1, filename2
    filename2= "outE" + filename1.split(".")[0][4:] + ".dat"
    outname = "out" + filename1.split(".")[0][4:] + ".vtk"
    sol2 = PetIGA().read_vec(filename1,geom)
    sol1 = PetIGA().read_vec(filename2,geom)
    sol=np.concatenate((sol1,sol2),axis=3);
    nrb=NURBS(geom.knots, (geom.points,geom.weights),sol)
    VTK().write(outname,nrb,scalars=dict(e2=0, e3=1, well=2, dist=3),vectors=dict(Displacement=[6,7,8],displacement=[9,10,11]))
