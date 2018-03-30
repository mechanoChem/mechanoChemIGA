import glob
from igakit.io import PetIGA,VTK
from igakit.nurbs import NURBS
import numpy as np
import itertools
import scipy.io
import os, sys

scalars = dict()
vectors = dict()

fin = open('fieldInfo.txt','r')
data = fin.read().splitlines()
dim = int(data[0])

l = data[3].split()
nPrjtnScalars = int(l[0])
for i in range(nPrjtnScalars):
    scalars[l[i+1]] = i

l = data[4].split()
nPrjtnVectors = int(l[0])
for i in range(nPrjtnVectors):
    j = nPrjtnScalars+i*dim
    vectors[l[i+1]] = range(j,j+dim)

l = data[1].split()
nSolScalars = int(l[0])
for i in range(nSolScalars):
    scalars[l[i+1]] = nPrjtnScalars+dim*nPrjtnVectors+i

l = data[2].split()
nSolVectors = int(l[0])
for i in range(nSolVectors):
    j = nPrjtnScalars+dim*nPrjtnVectors+nSolScalars+i*dim
    vectors[l[i+1]] = range(j,j+dim)

data={};
cwdir=os.getcwd();
geom = PetIGA().read('mesh.dat')
inc=1;
if (nPrjtnScalars + nPrjtnVectors > 0):
    for filename1, filename2 in zip(glob.glob('outU*.dat'), glob.glob('outE*.dat')):

        filename2= "outE" + filename1.split(".")[0][4:] + ".dat"
        outname = "out" + filename1.split(".")[0][4:] + ".vtk"
        sol2 = PetIGA().read_vec(filename1,geom)
        sol1 = PetIGA().read_vec(filename2,geom)
        if (sol1.ndim == dim):
            sol1 = np.expand_dims(sol1, axis=dim)
        if (sol2.ndim == dim):
            sol2 = np.expand_dims(sol2, axis=dim)
        sol=np.concatenate((sol1,sol2),axis=dim);
        nrb=NURBS(geom.knots, (geom.points,geom.weights),sol)
        VTK().write(outname,nrb,scalars=scalars,vectors=vectors)
else:
    for filename in glob.glob('outU*.dat'):

        outname = "out" + filename.split(".")[0][4:] + ".vtk"
        sol = PetIGA().read_vec(filename,geom)
        nrb=NURBS(geom.knots, (geom.points,geom.weights),sol)
        VTK().write(outname,nrb,scalars=scalars,vectors=vectors)
