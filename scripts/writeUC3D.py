import glob
from igakit.io import PetIGA,VTK
from igakit.nurbs import NURBS
import numpy as np
import itertools
import scipy.io

writeMat=0;

geom = PetIGA().read('mesh.dat')
data={};
inc=1;
for filename1, filename2 in zip(glob.glob('outU*.dat'), glob.glob('outE*.dat')):
    #print filename1, filename2
    outname = "out" + filename1.split(".")[0][4:] + ".vtk"
    sol1 = PetIGA().read_vec(filename1,geom)
    sol2 = PetIGA().read_vec(filename2,geom)
    sol=np.concatenate((sol1,sol2),axis=3);
    nrb=NURBS(geom.knots, (geom.points,geom.weights),sol)
    VTK().write(outname,nrb,scalars=dict(c=3, e2=4, e3=5, well=6, depth=7),vectors=dict(displacement=[0,1,2]))

