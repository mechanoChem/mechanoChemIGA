import glob
from igakit.io import PetIGA,VTK
from igakit.nurbs import NURBS
import numpy as np
import itertools
import scipy.io
import os, sys

writeMat=0;

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
    #sol=np.concatenate((sol1,sol2),axis=2);
    nrb=NURBS(geom.knots, (geom.points,geom.weights),sol)
    #VTK().write(outname,nrb,scalars=dict(e2=0, e3=1),vectors=dict(Displacement=[2,3]))
    VTK().write(outname,nrb,scalars=dict(e2=0, e3=1, well=2, dist=3),vectors=dict(Displacement=[6,7,8],displacement=[9,10,11]))
    #VTK().write(outname,nrb,scalars=dict(e2=0, e3=1, well=2, dist=3),vectors=dict(Displacement=[4,5],displacement=[6,7]))
    if (writeMat):
        #write mesh into X
        if not(data):
            X = geom.points
            W = geom.weights
            val=np.zeros((X.shape[0]*X.shape[1], 2))
            temp=-1
            for i, j in itertools.product(range(X.shape[0]), range(X.shape[1])):
                temp=temp+1
                for l in range(2):
                    val[temp,l]=X[i][j][l]
            data['X'] = val
            print "wrote mesh";
        #write solution fields
        val=np.zeros((sol1.shape[0]*sol1.shape[1], 2))
        temp=-1
        for i, j in itertools.product(range(sol1.shape[0]), range(sol1.shape[1])):
            temp=temp+1
            val[temp,0]=sol1[i][j][0] #e2
            val[temp,1]=sol1[i][j][1] #e3
        data['T'+filename1[4:-4]] = val
        print "wrote T" + filename1[4:-4]; 
        data['numIncs']=inc; inc=inc+1;
if writeMat:
    scipy.io.savemat('values2DCE2.mat', data)
