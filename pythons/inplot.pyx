#!/usr/bin/python

import sys,os
import struct
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
from scipy import integrate
from cpython cimport array
from ctypes import *
import array


def get_path(directory):
    path=os.path.join(os.getcwd(),directory)

    return path

def export_parameters(dir_path,file_name):
    path     = os.path.join(os.path.join(dir_path,file_name))
    file_read=open(path+'.bin','rb')
    
    file_read.seek(4)
    Nparam=struct.unpack('i',file_read.read(4))[0]
    file_read.seek(16)
    integerdata_temp=struct.unpack(Nparam*'i',file_read.read(Nparam*4))
    file_read.seek(104)
    realdata_temp=struct.unpack(Nparam*'f',file_read.read(Nparam*4))

    file_read.close()

    integerdata=np.array(integerdata_temp[:],dtype=[('proc_y','i4'),('proc_z','i4'),('proc_x','i4'),('nx_tot','i4'),('ny_tot','i4'),('ny_per_proc','i4'),('nz_tot','i4'),('nz_per_proc','i4'),('boundary_x','i4'),('boundary_y','i4'),('boundary_z','i4'),('env_non_env','i4'),('unknown_param2','i4'),('unknown_param3','i4'),('ndim','i4'),('unknown_param5','i4'),('unknown_param6','i4'),('unknown_param7','i4'),('unknown_param8','i4'),('unknown_param9','i4')])

    realdata=np.array(realdata_temp[:],dtype=[('time','f4'),('x_min','f4'),('x_max','f4'),('y_min','f4'),('y_max','f4'),('z_min','f4'),('z_max','f4'),('pulse_duration','f4'),('waist','f4'),('n_over_nc','f4'),('a0','f4'),('lambda_0','f4'),('E0','f4'),('unknown','f4'),('np_per_cell','f4'),('unknown_param5','f4'),('unknown_param6','f4'),('unknown_param7','f4'),('unknown_param8','f4'),('unknown_param9','f4')])

    return (integerdata,realdata)

def read_ALaDyn_bin(dir_path,file_name,grid_no_grid):
    
    # - #

    cdef int i, j, k, nx, ny, nz, nproc_y, nproc_z, N_param, offset
    cdef int npx, npy, npz, offsety, offsetz
    cdef int counter_z,counter_y
    path     = os.path.join(os.path.join(dir_path,file_name+'.bin'))
    f        = open(str(path),'rb')
    path     = c_char_p(os.path.join(os.path.join(dir_path,file_name+'.bin')))
    cpath='/home/HPC/terzani/Pic/tools-ALaDyn/pythons'
    read_bin= CDLL(cpath+'/read_binary.so')
    #- vector length -#
    struct.unpack('i', f.read(4))
    N_param = struct.unpack('i', f.read(4))[0]
    print(N_param)
    struct.unpack('i', f.read(4))
    struct.unpack('i', f.read(4))
    int_param=[]
    for i in range(0,N_param):
        int_param.append( struct.unpack('i', f.read(4)) )
    struct.unpack('i', f.read(4))
    #nx= c_uint(int_param[3][0])
    #ny= c_uint(int_param[4][0])
    #nz= c_uint(int_param[6][0])
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
    #nx=int(nx.value)
    #ny=int(ny.value)
    #nz=int(nz.value)
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

    if grid_no_grid == 'nogrid':
        return r
    if (ndimension==2):
        return (r,x,y)

    return (r,x,y,z)

def read_particle_phasespace_bycomponent(dir_path,file_name,component):


    path     = os.path.join(dir_path,file_name)
    f        = open(path+'.bin','rb')
    cmp = []
    while True:
        try:
            Np=struct.unpack('i', f.read(4))
        except struct.error:
            if(component == 'X'):  print component+'component has been read'
            if(component == 'Y'):  print component+'component has been read'
            if(component == 'Z'):  print component+'component has been read'
            if(component == 'Px'): print component+'component has been read'
            if(component == 'Py'): print component+'component has been read'
            if(component == 'Pz'): print component+'component has been read'
            if(component == 'W'):  print component+'component has been read'
            if(component == 'Q'):  print component+'component has been read'
            break
        jump= 8 #2*ndim+2
        vars=[]
        try:
            
            for i in range(0,Np[0]):
                vars=struct.unpack(jump*'f', f.read(4*jump))
                if(component == 'X'): 	cmp.append(vars[0])
                if(component == 'Y'): 	cmp.append(vars[1])
                if(component == 'Z'): 	cmp.append(vars[2])
                if(component == 'Px'): 	cmp.append(vars[3])
                if(component == 'Py'): 	cmp.append(vars[4])
                if(component == 'Pz'): 	cmp.append(vars[5])
                if(component == 'W'): 	cmp.append(vars[6])
                if(component == 'Q'): 	cmp.append(vars[7])
        except:
            pass
    f.close()
    cmp=np.asarray(cmp)
    return np.array(cmp)

def read_particle_phasespace_bycomponent_2d(dir_path,file_name,component):


    path     = os.path.join(dir_path,file_name)  
    f        = open(path+'.bin','rb')
    cmp = []
    while True:
        try:			
            Np=struct.unpack('i', f.read(4))
        except struct.error:
            if(component == 'X'):  print component+'component has been read'
            if(component == 'Y'):  print component+'component has been read'
            if(component == 'Px'): print component+'component has been read'
            if(component == 'Py'): print component+'component has been read'
            if(component == 'W'):  print component+'component has been read'
            if(component == 'Q'):  print component+'component has been read'
            break
        jump=6 #2*ndim+2
        vars=[]
        for i in range(0,Np[0]):
            vars=struct.unpack('f'*jump, f.read(4*jump))
            if(component == 'X'): 	cmp.append(vars[0])
            if(component == 'Y'): 	cmp.append(vars[1])
            if(component == 'Px'): 	cmp.append(vars[2])
            if(component == 'Py'): 	cmp.append(vars[3])
            if(component == 'W'): 	cmp.append(vars[4])
            if(component == 'Q'): 	cmp.append(vars[5])
    f.close()
    cmp=np.asarray(cmp)
    return np.array(cmp)

def get_a(path,a_field_number,grid_no_grid):

    integerparam,realparam=export_parameters(path,'Renvout'+a_field_number)
    
    ndimension=int(integerparam['ndim'])
    a_real=read_ALaDyn_bin(path,'Renvout'+a_field_number,'nogrid')
    Envelope=np.zeros(a_real.shape[0])
    if grid_no_grid=='grid':
        if(ndimension==3):
            a_im,x,y,z=read_ALaDyn_bin(path,'Ienvout'+a_field_number,'grid')
            Envelope=np.sqrt(a_real**2+a_im**2)
            return (Envelope,x,y,z)
        if(ndimension==2):
            a_im,x,y=read_ALaDyn_bin(path,'Ienvout'+a_field_number,'grid')
            Envelope=np.sqrt(a_real**2+a_im**2)
            return (Envelope,x,y)
    a_im=read_ALaDyn_bin(path,'Ienvout'+a_field_number,'nogrid')
    Envelope=np.sqrt(a_real**2+a_im**2)
    return Envelope

def convert_a_in_e(a_field_number,path):

    from scipy.interpolate import interp1d
    from scipy.signal import hilbert
    import scipy.integrate as integrate
    cdef int i, j, k
    integerparam,realparam=export_parameters(path,'Renvout'+a_field_number)
    

    time=realparam['time']
    nx=integerparam['nx_tot']
    ny=integerparam['ny_tot']
    nz=integerparam['nz_tot']
    lam_0=realparam['lambda_0']
    speed_of_light=1.0
    omega_0=2*np.pi*speed_of_light/lam_0
    k_0=2*np.pi/lam_0
    e_field=np.zeros((nx,ny,nz))

    ndimension=int(integerparam['ndim'])



    if(ndimension==2):
        a_real,x,y=read_ALaDyn_bin(path,'Renvout'+a_field_number,'grid')
        a_im=read_ALaDyn_bin(path,'Ienvout'+a_field_number,'nogrid')
  
        phi=k_0*(x-time)
        deltax=(x[-1]-x[0])/float(nx)
        a_oscillating=a_real
        for j in range(ny):
            a_oscillating[:,j]=a_real[:,j]*np.cos(phi)-a_im[:,j]*np.sin(phi)

    if(ndimension==3):
        a_real,x,y,z=read_ALaDyn_bin(path,'Renvout'+a_field_number,'grid')
        a_im=read_ALaDyn_bin(path,'Ienvout'+a_field_number,'nogrid')
        phi=k_0*(x-time)
        deltax=(x[-1]-x[0])/float(nx)
        a_oscillating=a_real
        for k in range(nz): 
            for j in range(ny):
                a_oscillating[:,j,k]=a_real[:,j,k]*np.cos(phi)-a_im[:,j,k]*np.sin(phi)
    
    a_diff=np.diff(a_oscillating,axis=0)/deltax
    if(ndimension==2):
        a_diff=np.append(a_diff,[a_diff[-1,:]],axis=0)
    if(ndimension==3):
        a_diff=np.append(a_diff,[a_diff[-1,:,:]],axis=0)

    e_field=a_diff
#   a_real_diff=np.diff(a_real,axis=0)/deltax
#   a_im_diff=np.diff(a_im,axis=0)/deltax
#   a_real_diff=np.append(a_real_diff,[a_real_diff[-1,:,:]],axis=0)
#   a_im_diff=np.append(a_im_diff,[a_im_diff[-1,:,:]],axis=0)

#   for k in range(nz):
#       for j in range(ny):
#           discrete_cos=np.zeros(nx)
#           discrete_sin=np.zeros(nx)
#           for i in range(nx):
#               discrete_cos[i]=(speed_of_light*a_real_diff[i,j,k]+omega_0*a_im[i,j,k])*np.cos(k_0*x[i])
#               discrete_sin[i]=-(speed_of_light*a_im_diff[i,j,k]+omega_0*a_real[i,j,k])*np.sin(k_0*x[i])
#           envelope=np.abs(hilbert(discrete_cos+discrete_sin))
#           e_field[:,j,k]=envelope
    if(ndimension==2):
        return (e_field,x,y)

    return (e_field,x,y,z)

def tracking_parameters(dir_path,file_number):

    path=os.path.join(dir_path,'El_track_out'+file_number+'.dat')

    file_dat=open(path,'r')

    lines=file_dat.readlines()
    tracking_parameters=np.zeros(1,dtype=[('time','f4'),('dt','f4'),('tot_nproc','i4'),('phasespace_dimensions','i4'),('spatial_dimensions','i4'),('number_timesteps','i4'),('jumped_timestep','i4'),('current_tot_particles','i4'),('initial_tot_particles','i4')])[0]
    tracking_parameters[0]=float(lines[2])
    tracking_parameters[1]=float(lines[4])
    tracking_parameters[2]=int(lines[7])
    tracking_parameters[3]=int(lines[9])
    tracking_parameters[4]=int(lines[11])
    tracking_parameters[5]=int(lines[13])
    tracking_parameters[6]=int(lines[15])
    tracking_parameters[7]=int(lines[17].split()[0])
    tracking_parameters[8]=int(lines[17].split()[1])

    return tracking_parameters


def track_data(dir_path,file_number):

    path=os.path.join(dir_path,'El_track_out'+file_number)

    file_dat=open(path+'.dat','r')
    lines=file_dat.readlines()
    time=float(lines[2])
    dt=float(lines[4])
    tot_nproc=int(lines[7])
    phasespace_dimensions=int(lines[9])
    spatial_dimensions=int(lines[11])
    number_timesteps=int(lines[13])
    timestep_jump=int(lines[15])
    current_tot_particles=int(lines[17].split()[0])
    initial_tot_particles=int(lines[17].split()[1])
    if (phasespace_dimensions==8):
         track=np.zeros((number_timesteps,initial_tot_particles),dtype=[('x','d'),('y','d'),('z','d'),('vx','d'),('vy','d'),('vz','d'),('gamma','d'),('nindex','i')])
    else:
         track=np.zeros((number_timesteps,initial_tot_particles),dtype=[('x','d'),('y','d'),('vx','d'),('vy','d'),('gamma','d'),('nindex','i')])

    xp=np.zeros((number_timesteps,initial_tot_particles))
    yp=np.zeros((number_timesteps,initial_tot_particles))
    vxp=np.zeros((number_timesteps,initial_tot_particles))
    vyp=np.zeros((number_timesteps,initial_tot_particles))
    gamma=np.zeros((number_timesteps,initial_tot_particles))
    part_index=np.zeros((number_timesteps,initial_tot_particles),dtype=np.int)
    check_index=np.full(initial_tot_particles,1)


    if (phasespace_dimensions==8):
        zp=np.zeros((number_timesteps,initial_tot_particles))
        vzp=np.zeros((number_timesteps,initial_tot_particles))

    jump=phasespace_dimensions*number_timesteps

    file_bin=open(path+'.bin','rb')

    if(phasespace_dimensions==6):
        while True:
            try:			
                Np=struct.unpack('i', file_bin.read(4))
            except struct.error:
                break
            for k in range(0,Np[0]):
                vars=struct.unpack(jump*'f', file_bin.read(4*jump))
                for j in range(0,number_timesteps):
                    kk=int(vars[j*phasespace_dimensions+5])-1
                    xp[j,kk]=float(vars[j*phasespace_dimensions])
                    yp[j,kk]=float(vars[j*phasespace_dimensions+1])
                    vxp[j,kk]=float(vars[j*phasespace_dimensions+2])
                    vyp[j,kk]=float(vars[j*phasespace_dimensions+3])
                    gamma[j,kk]=float(vars[j*phasespace_dimensions+4])
                    part_index[j,kk]=kk+1
        
    if(phasespace_dimensions==7):
        while True:
            try:			
                Np=struct.unpack('i', file_bin.read(4))
            except struct.error:
                break
            for k in range(0,Np[0]):
                vars=struct.unpack(jump*'f', file_bin.read(4*jump))
                for j in range(0,number_timesteps):
                    kk=int(vars[j*phasespace_dimensions+7])-1
                    xp[j,kk]=float(vars[j*phasespace_dimensions])
                    yp[j,kk]=float(vars[j*phasespace_dimensions+1])
                    zp[j,kk]=float(vars[j*phasespace_dimensions+2])
                    vxp[j,kk]=float(vars[j*phasespace_dimensions+3])
                    vyp[j,kk]=float(vars[j*phasespace_dimensions+4])
                    vzp[j,kk]=float(vars[j*phasespace_dimensions+5])
                    gamma[j,kk]=float(vars[j*phasespace_dimensions+6])
                    part_index[j,kk]=kk+1

    track['x']=xp
    track['y']=yp
    track['vx']=vxp
    track['vy']=vyp
    track['gamma']=gamma
    track['nindex']=part_index
    if(phasespace_dimensions==7):
        track['z']=zp
        track['vz']=vzp
   
    return track

def nearest_index(track,component,number):
        
        ind=[]
        ind.append(track['nindex'][0,0])
        current_diff=np.abs(track[component][0,0]-number)
        for i in range(1,track[component].shape[1]):
            if (np.abs(track[component][0,i]-number)<current_diff):
                ind=[]
                current_diff=np.abs(track[component][0,i]-number)
                ind.extend([track['nindex'][0,i],i])
            elif (np.abs(track[component][0,i]-number)==current_diff):
                ind.extend([track['nindex'][0,i],i])
        print 'The particle(s) with component', component,'closest to', number,'has(have) index', ind[::2]
        print 'The distance is', current_diff           
        print 'Its(Their) initial coordinates are:'
       
        for i in range(1,len(ind)/2+1):
            print 'x=',track['x'][0,ind[2*i-1]],' y=',track['y'][0,ind[2*i-1]]

        print 'Its(Their) position in the array is (are):'
        
        for i in range(1,len(ind)/2+1):
            print 'Array index to be used',ind[2*i-1]


def define_tracking_time_array(path,file_number):

        track_param=tracking_parameters(path,file_number)
        steps=track_param['number_timesteps']
        dt=track_param['dt']
        time=np.arange(steps)
        time=time*dt

        return time

def diff_functions(case_number,file_number, *args, **kwargs): 

        dir_path=os.getcwd()
        pathenv=os.path.join(dir_path,'env'+case_number+'/00'+file_number)

        pathpic=os.path.join(dir_path,'pic'+case_number+'/00'+file_number)
        _,yp,_,vyp,_=track_data(pathpic,file_number)

        _,ye,_,vye,_=track_data(pathenv,file_number)

        interpenv=interp.interp1d(ye[0,:],vye[-1,:],kind='cubic')

        interppic=interp.interp1d(yp[0,:],vyp[-1,:],kind='cubic')

        ymin=max(np.amin(ye[0,:]),np.amin(yp[0,:]))
        ymax=min(np.amax(ye[0,:]),np.amax(yp[0,:]))
        diff=integrate.quad(lambda x: np.abs(interpenv(x)-interppic(x)),ymin,ymax)

        print 'Difference is ',diff[0],'with an error of ', diff[1]

        return (interpenv,interppic)


def plot_env_pic(**kwargs):

        plot_number=0
        plot_label='a0'
        folder_name=4
        component='transverse'
        leg_pos=0
        case_number=[]
        pic_or_not=True
        if ('first' in kwargs):
                case_number.append(kwargs['first'])
                plot_number+=1
        if ('second' in kwargs):
                case_number.append(kwargs['second'])
                plot_number+=1
        if ('third' in kwargs):
                case_number.append(kwargs['third'])
                plot_number+=1
        if ('fourth' in kwargs):
                case_number.append(kwargs['fourth'])
                plot_number+=1
        if ('fifth' in kwargs):
                case_number.append(kwargs['fifth'])
                plot_number+=1
        if ('time' in kwargs):
                folder_name=kwargs['time']
        if ('pic' in kwargs):
                pic_or_not=kwargs['pic']
        if ('component' in kwargs):
                component=kwargs['component']
	
        dir_path=os.getcwd()

        if('label' in kwargs):
                plot_label=kwargs['label']
        labenv=np.array([''],dtype=[('a0','U30'),('tau','U30'),('w0','U30'),('lam0','U30')])
        plt.ion()
        f=plt.figure(figsize=(13,7))
        ax=f.add_subplot(111)
        colors=["#FF0066","#00FF00","#FFDE00","#996600","#FF66FF"]
        colors2=["#800033","#008000","#806F00","#332200","#B300B3"]
        if(pic_or_not):
		#colors2=['b','b','b','b','b']
                transparency=0.4
                markenv='*'
        else:
		#colors2=["#FF0066","#00FF00","#FFDE00","#996600","#FF66FF"]
                transparency=1
                markenv='x'
        for i in range(plot_number):

                pathenv=os.path.join(dir_path,'env'+str(case_number[i]),str(folder_name).zfill(4))
                pathpic=os.path.join(dir_path,'pic'+str(case_number[i]),str(folder_name).zfill(4))
                if(pic_or_not):
                         trackpic=track_data(pathpic,str(folder_name).zfill(2))
                trackenv=track_data(pathenv,str(folder_name).zfill(2))
	
                intdata,realdata=export_parameters(pathenv,'Exfout'+str(folder_name).zfill(2))
	
                a0=round(realdata['a0'],1)
                tau=round(realdata['pulse_duration']*2.*np.arccos(1/2**0.25)/np.pi/0.3,1)
                waist=realdata['waist']
                lam0=round(realdata['lambda_0'],1)
                if not (pic_or_not):
                        labenv['a0']='$a_0= $'+str(a0)
                        labenv['tau']='$L_{FWHM}= $'+str(tau)+'$fs$'
                        labenv['w0']='$w_0= $'+str(waist)+'$\mu m$'
                        labenv['lam0']='$\lambda_0= $'+str(lam0)+'$\mu m$'

                if (component=='transverse'):
                        ax.set_ylabel('$p_x/mc$',fontsize=18)
                        if (plot_label=='a0'):
				
                                if(pic_or_not):
                                        ax.scatter(trackpic['y'][0,:],trackpic['vy'][-1,:],c=colors[i],s=6,label='$a_0= $'+str(a0))
                                ax.scatter(trackenv['y'][0,:],trackenv['vy'][-1,:],marker=markenv,s=1,c=colors2[i],alpha=transparency,label=labenv['a0'][0])
			
                        if (plot_label=='tau'):

                                if(pic_or_not):
                                        ax.scatter(trackpic['y'][0,:],trackpic['vy'][-1,:],c=colors[i],s=6,label='$L_{FWHM}= $'+str(tau)+'$fs$')
                                ax.scatter(trackenv['y'][0,:],trackenv['vy'][-1,:],marker=markenv,s=1,c=colors2[i],alpha=transparency,label=labenv['tau'][0])
			
                        if (plot_label=='w0'):

                                if(pic_or_not):
                                        ax.scatter(trackpic['y'][0,:],trackpic['vy'][-1,:],c=colors[i],s=6,label='$w_0= $'+str(waist)+'$\mu m$')
                                ax.scatter(trackenv['y'][0,:],trackenv['vy'][-1,:],marker=markenv,s=1,c=colors2[i],alpha=transparency,label=labenv['w0'][0])

                        if (plot_label=='lambda'):

                                if(pic_or_not):
                                        ax.scatter(trackpic['y'][0,:],trackpic['vy'][-1,:],c=colors[i],s=6,label='$\lambda_0= $'+str(lam0)+'$\mu m$')
                                ax.scatter(trackenv['y'][0,:],trackenv['vy'][-1,:],marker=markenv,s=1,c=colors2[i],alpha=transparency,label=labenv['lam0'][0])
                
                if (component=='longitudinal'):
                        ax.set_ylabel('$p_z/mc$',fontsize=18)
                        if (plot_label=='a0'):
                                if(pic_or_not):
                                        ax.scatter(trackpic['y'][0,:],trackpic['vx'][-1,:],c=colors[i],s=6,label='$a_0= $'+str(a0))
                                ax.scatter(trackenv['y'][0,:],trackenv['vx'][-1,:],marker=markenv,s=1,c=colors2[i],alpha=transparency,label=labenv['a0'][0])

                        if (plot_label=='tau'):

                                if(pic_or_not):
                                        ax.scatter(trackpic['y'][0,:],trackpic['vx'][-1,:],c=colors[i],s=6,label='$L_{FWHM}= $'+str(tau)+'$fs$')
                                ax.scatter(trackenv['y'][0,:],trackenv['vx'][-1,:],marker=markenv,s=1,c=colors2[i],alpha=transparency,label=labenv['tau'][0])

                        if (plot_label=='w0'):

                                if(pic_or_not):
                                        ax.scatter(trackpic['y'][0,:],trackpic['vx'][-1,:],c=colors[i],s=6,label='$w_0= $'+str(waist)+'$\mu m$')
                                ax.scatter(trackenv['y'][0,:],trackenv['vx'][-1,:],marker=markenv,s=1,c=colors2[i],alpha=transparency,label=labenv['w0'][0])
                        
                        if (plot_label=='lambda'):

                                if(pic_or_not):
                                        ax.scatter(trackpic['y'][0,:],trackpic['vx'][-1,:],c=colors[i],s=6,label='$\lambda_0= $'+str(lam0)+'$\mu m$')
                                ax.scatter(trackenv['y'][0,:],trackenv['vx'][-1,:],marker=markenv,s=1,c=colors2[i],alpha=transparency,label=labenv['lam0'][0])

        if ('legend_position' in kwargs):
                 leg_pos=kwargs['legend_position']
	
        ax.legend(loc=leg_pos,fontsize=18)
        ax.set_xlabel('x $(\mu m)$',fontsize=18)

def data(path,label,n_dimensions):
   
    import os.path as os

    if(label=='env'):
         global Exenv
         Exenv=[]
         global Eyenv
         Eyenv=[]
         global rhoenv
         rhoenv=[]
         global A
         A=[]
         global Bzenv
         Bzenv=[]
         global enerenv
         enerenv=[]
         global xenv
         xenv=[]
         global yenv
         yenv=[]
         global zenv
         zenv=[]
         global rhofluid
         rhofluid=[]

    if(label=='pic'):
        global Expic
        Expic=[]
        global Eypic
        Eypic=[]
        global rhopic
        rhopic=[]
        global Bzpic
        Bzpic=[]
        global enerpic
        enerpic=[]
        global xpic
        xpic=[]
        global ypic
        ypic=[]
        global zpic
        zpic=[]

    if(n_dimensions==2):
        for i in path:
            if(label=='env'):
                if(os.isfile(os.join(i,'Exfout'+i[-2:]+'.bin'))):
                    print 'Now reading Exfout',i[-2:]
                    temp,tempx,tempy=read_ALaDyn_bin(i,'Exfout'+i[-2:],'grid')
                    Exenv.append(temp)
                    xenv.append(tempx)
                    yenv.append(tempy)
                if(os.isfile(os.join(i,'Aenvout'+i[-2:]+'.bin'))):
                    print 'Now reading Aenvout',i[-2:]
                    A.append(read_ALaDyn_bin(i,'Aenvout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Renvout'+i[-2:]+'.bin'))):
                    print 'Now reading Aenvout',i[-2:]
                    A.append(get_a(i,i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Eyfout'+i[-2:]+'.bin'))):
                    print 'Now reading Eyfout',i[-2:]
                    Eyenv.append(read_ALaDyn_bin(i,'Eyfout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Edenout'+i[-2:]+'.bin'))):
                    print 'Now reading Edenout',i[-2:]
                    rhoenv.append(read_ALaDyn_bin(i,'Edenout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Fdenout'+i[-2:]+'.bin'))):
                    print 'Now reading Fdenout',i[-2:]
                    rhofluid.append(read_ALaDyn_bin(i,'Edenout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Bzfout'+i[-2:]+'.bin'))):
                    print 'Now reading Bzfout',i[-2:]
                    Bzenv.append(read_ALaDyn_bin(i,'Bzfout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Elenout'+i[-2:]+'.bin'))):
                    print 'Now reading Elenout',i[-2:]
                    enerenv.append(read_ALaDyn_bin(i,'Elenout'+i[-2:],'nogrid'))
            if(label=='pic'):
                if(os.isfile(os.join(i,'Exfout'+i[-2:]+'.bin'))):
                    print 'Now reading Exfout',i[-2:]
                    temp,tempx,tempy=read_ALaDyn_bin(i,'Exfout'+i[-2:],'grid')
                    Expic.append(temp)
                    xpic.append(tempx)
                    ypic.append(tempy)
                if(os.isfile(os.join(i,'Eyfout'+i[-2:]+'.bin'))):
                    print 'Now reading Eyfout',i[-2:]
                    Eypic.append(read_ALaDyn_bin(i,'Eyfout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Edenout'+i[-2:]+'.bin'))):
                    print 'Now reading Edenout',i[-2:]
                    rhopic.append(read_ALaDyn_bin(i,'Edenout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Bzfout'+i[-2:]+'.bin'))):
                    print 'Now reading Bzfout',i[-2:]
                    Bzpic.append(read_ALaDyn_bin(i,'Bzfout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Elenout'+i[-2:]+'.bin'))):
                    print 'Now reading Elenout',i[-2:]
                    enerpic.append(read_ALaDyn_bin(i,'Elenout'+i[-2:],'nogrid'))
    if(n_dimensions==3):
        for i in path:
            if(label=='env'):
                if(os.isfile(os.join(i,'Exfout'+i[-2:]+'.bin'))):
                    print 'Now reading Exfout',i[-2:]
                    temp,tempx,tempy,tempz=read_ALaDyn_bin(i,'Exfout'+i[-2:],'grid')
                    Exenv.append(temp)
                    xenv.append(tempx)
                    yenv.append(tempy)
                    zenv.append(tempz)
                if(os.isfile(os.join(i,'Aenvout'+i[-2:]+'.bin'))):
                    print 'Now reading Aenvout',i[-2:]
                    A.append(read_ALaDyn_bin(i,'Aenvout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Renvout'+i[-2:]+'.bin'))):
                    print 'Now reading Aenvout',i[-2:]
                    A.append(get_a(i,i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Eyfout'+i[-2:]+'.bin'))):
                    print 'Now reading Eyfout',i[-2:]
                    Eyenv.append(read_ALaDyn_bin(i,'Eyfout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Edenout'+i[-2:]+'.bin'))):
                    print 'Now reading Edenout',i[-2:]
                    rhoenv.append(read_ALaDyn_bin(i,'Edenout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Fdenout'+i[-2:]+'.bin'))):
                    print 'Now reading Fdenout',i[-2:]
                    rhofluid.append(read_ALaDyn_bin(i,'Fdenout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Bzfout'+i[-2:]+'.bin'))):
                    print 'Now reading Bzfout',i[-2:]
                    Bzenv.append(read_ALaDyn_bin(i,'Bzfout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Elenout'+i[-2:]+'.bin'))):
                    print 'Now reading Elenout',i[-2:]
                    enerenv.append(read_ALaDyn_bin(i,'Elenout'+i[-2:],'nogrid'))
            if(label=='pic'):
                if(os.isfile(os.join(i,'Exfout'+i[-2:]+'.bin'))):
                    print 'Now reading Exfout',i[-2:]
                    temp,tempx,tempy,tempz=read_ALaDyn_bin(i,'Exfout'+i[-2:],'grid')
                    Expic.append(temp)
                    xpic.append(tempx)
                    ypic.append(tempy)
                    zpic.append(tempz)
                if(os.isfile(os.join(i,'Eyfout'+i[-2:]+'.bin'))):
                    print 'Now reading Eyfout',i[-2:]
                    Eypic.append(read_ALaDyn_bin(i,'Eyfout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Edenout'+i[-2:]+'.bin'))):
                    print 'Now reading Edenout',i[-2:]
                    rhopic.append(read_ALaDyn_bin(i,'Edenout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Bzfout'+i[-2:]+'.bin'))):
                    print 'Now reading Bzfout',i[-2:]
                    Bzpic.append(read_ALaDyn_bin(i,'Bzfout'+i[-2:],'nogrid'))
                if(os.isfile(os.join(i,'Elenout'+i[-2:]+'.bin'))):
                    print 'Now reading Elenout',i[-2:]
                    enerpic.append(read_ALaDyn_bin(i,'Elenout'+i[-2:],'nogrid'))


def temporal_average(gamma,time,lambda_0,n):

    omega_0=2*np.pi/lambda_0
    t=lambda_0
    delta_t=(time[-1]-time[0])/time.shape[0]
    gammaint=np.full(gamma.shape[0],1./delta_t)
    gammaroot_mean=np.full(gamma.shape[0],1./delta_t)
    timesteps_in_period=t/delta_t
    a=n*int(timesteps_in_period/2.)
    for i in range(a,gamma.shape[0]-a):
       gammaint[i]=np.sum(gamma[i-a:i+a])
    gammaint=gammaint*delta_t/(n*t)
    for i in range(a,gamma.shape[0]-a):
       gammaroot_mean[i]=np.sum(gamma[i-a:i+a]**2)
    gammaroot_mean=np.sqrt(gammaroot_mean*delta_t/(n*t))

    
    return gammaint,gammaroot_mean
    
