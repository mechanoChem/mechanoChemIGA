import glob
from igakit.io import PetIGA,VTK
from igakit.nurbs import NURBS
import numpy as np
import itertools
import scipy.io

writeMat=1;

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
    VTK().write(outname,nrb,scalars=dict(c=3,e2=4, e3=5, dist=7, well=6),vectors=dict(displacement=[0,1,2]))
    if (writeMat):
        if data:
            X = geom.points
            W = geom.weights
            val=np.zeros((X.shape[0]*X.shape[1]*X.shape[2], 3))
            temp=-1
            for i, j, k in itertools.product(range(X.shape[0]), range(X.shape[1]), range(X.shape[2])):
                temp=temp+1
                for l in range(3):
                    val[temp,l]=X[i][j][k][l]
            data['X'] = val
        val=np.zeros((sol1.shape[0]*sol1.shape[1]*sol1.shape[2], 3))
        temp=-1
        for i, j, k in itertools.product(range(sol1.shape[0]), range(sol1.shape[1]), range(sol1.shape[2])):
            temp=temp+1
            val[temp,0]=sol1[i][j][k][3]
            val[temp,1]=sol2[i][j][k][0]
            val[temp,2]=sol2[i][j][k][1]
        data['T'+filename1[4:-4]] = val
        data['numIncs']=inc; inc=inc+1;
if writeMat:
    scipy.io.savemat('values3DE23.mat', data)

