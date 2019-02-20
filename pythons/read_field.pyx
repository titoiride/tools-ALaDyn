#!/usr/bin/python

import os
import struct
import numpy as np

from cpython cimport array
from ctypes import *
import array

def read_ALaDyn_bin(dir_path,file_name,grid_no_grid):
    
    # - #

    cdef int i, j, k, nx, ny, nz, nproc_y, nproc_z, N_param, offset
    cdef int npx, npy, npz, offsety, offsetz
    cdef int counter_z,counter_y
    path     = os.path.join(os.path.join(dir_path,file_name+'.bin'))
    f        = open(str(path),'rb')
    path     = c_char_p(os.path.join(os.path.join(dir_path,file_name+'.bin')))
    cpath=os.path.join(os.environ["HOME"],'Pic/forked_tools-ALaDyn/pythons')
    read_bin= CDLL(cpath+'/read_binary.so')
    #- vector length -#
    struct.unpack('i', f.read(4))
    N_param = struct.unpack('i', f.read(4))[0]

    struct.unpack('i', f.read(4))
    struct.unpack('i', f.read(4))
    int_param=[]
    for i in range(0,N_param):
        int_param.append( struct.unpack('i', f.read(4)) )
    struct.unpack('i', f.read(4))
    nx= int(int_param[3][0])
    ny= int(int_param[4][0])
    nz= int(int_param[6][0])
    nproc_y = int(int_param[0][0])
    nproc_z = int(int_param[1][0])
    ndimension=int(int_param[14][0])
    struct.unpack('i', f.read(4))
    for i in range(0,N_param):
        struct.unpack('f', f.read(4))
    struct.unpack('i', f.read(4))

    #---***---#
    totlen=nx*ny*nz
    fieldarr=(c_float*totlen)    
    gridx=(c_float*nx)
    gridy=(c_float*ny)
    gridz=(c_float*nz)

    print('total grid size: n=(',nx,ny,nz,')')
    print('number of Np processors: Np_y=',nproc_y,'Np_z=', nproc_z)
    r=fieldarr()    
    x=gridx()
    y=gridy()
    z=gridz()


    read_bin.read_binary(byref(r), byref(x), byref(y), byref(z), path)
    r=np.ndarray(buffer=r,dtype=np.float32,shape=(nx,ny,nz),order='F')
    x=np.ndarray(buffer=x,dtype=np.float32,shape=(nx))
    y=np.ndarray(buffer=y,dtype=np.float32,shape=(ny))
    z=np.ndarray(buffer=z,dtype=np.float32,shape=(nz))

    if(ndimension==2):
        r=r[:,:,0]

    r=np.flip(r,1)
    if(ndimension==3):
        r=np.flip(r,3)

    if grid_no_grid == 'nogrid':
        return r
    if (ndimension==2):
        return (r,x,y)

    return (r,x,y,z)