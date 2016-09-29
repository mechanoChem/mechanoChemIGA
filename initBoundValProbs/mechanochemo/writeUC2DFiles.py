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
    sol=np.concatenate((sol1,sol2),axis=2);
    nrb=NURBS(geom.knots, (geom.points,geom.weights),sol)
    VTK().write(outname,nrb,scalars=dict(c=5, e2=0, e6=1, c2=2),vectors=dict(displacement=[3,4]))
