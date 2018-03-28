#!/usr/bin/python

import sys,os
import struct
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
from scipy import integrate


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

    integerdata=np.array(integerdata_temp[:],dtype=[('proc_y','i4'),('proc_z','i4'),('proc_x','i4'),('nx_tot','i4'),('ny_tot','i4'),('ny_per_proc','i4'),('nz_tot','i4'),('nz_per_proc','i4'),('boundary_x','i4'),('boundary_y','i4'),('boundary_z','i4'),('env_non_env','i4'),('unknown_param2','i4'),('unknown_param3','i4'),('unknown_param4','i4'),('unknown_param5','i4'),('unknown_param6','i4'),('unknown_param7','i4'),('unknown_param8','i4'),('unknown_param9','i4')])

    realdata=np.array(realdata_temp[:],dtype=[('time','f4'),('x_min','f4'),('x_max','f4'),('y_min','f4'),('y_max','f4'),('z_min','f4'),('z_max','f4'),('pulse_duration','f4'),('waist','f4'),('n_over_nc','f4'),('a0','f4'),('lambda_0','f4'),('E0','f4'),('unknown','f4'),('np_per_cell','f4'),('unknown_param5','f4'),('unknown_param6','f4'),('unknown_param7','f4'),('unknown_param8','f4'),('unknown_param9','f4')])

    return (integerdata,realdata)

def read_ALaDyn_bin(dir_path,file_name,grid_no_grid):
    
    # - #
    path     = os.path.join(os.path.join(dir_path,file_name+'.bin'))
    f        = open(path,'rb')

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
    nx= int_param[3][0]
    ny= int_param[4][0]
    nz= int_param[6][0]
    nproc_y = int_param[0][0]
    nproc_z = int_param[1][0]
    struct.unpack('i', f.read(4))
    for i in range(0,N_param):
        struct.unpack('f', f.read(4))
    struct.unpack('i', f.read(4))

    #---***---#
    r = np.zeros((nx,ny,nz))
    print('total grid size: n=(',nx,ny,nz,')')
    rr=[]
    print('number of Np processors: Np_y=',nproc_y,'Np_z=', nproc_z)

    offsetz = 0
    for counter_z in range(0,nproc_z):
        offsety = 0

        for counter_y in range(0,nproc_y):
            struct.unpack('i', f.read(4))
            npx= struct.unpack('i', f.read(4))[0]
            npy= struct.unpack('i', f.read(4))[0]
            npz= struct.unpack('i', f.read(4))[0]
            struct.unpack('i', f.read(4))
            
            struct.unpack('i', f.read(4))
            for k in range(0,npz):
                for j in range(0,npy):
                    for i in range(0,npx):
                        r[i,j+offsety,k+offsetz] = struct.unpack('f', f.read(4))[0]
            struct.unpack('i', f.read(4))
            offsety += npy
        offsetz += npz

    
    if grid_no_grid == 'nogrid':
    
        return r
    
    #--- * --- * --- * --- * --- * ---#
    #- reading grid -#
    struct.unpack('i', f.read(4))
    X=struct.unpack(nx*'f', f.read(nx*4))
    struct.unpack('i', f.read(4))
    
    struct.unpack('i', f.read(4))
    Y=struct.unpack(ny*'f', f.read(ny*4))
    struct.unpack('i', f.read(4))
    
    struct.unpack('i', f.read(4))
    Z=struct.unpack(nz*'f', f.read(nz*4))
    struct.unpack('i', f.read(4))
    X=np.asarray(X)
    Y=np.asarray(Y)
    Z=np.asarray(Z)
    
    x=X; y=Y; z=Z
    
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

    a_real=read_ALaDyn_bin(path,'Renvout'+a_field_number,'nogrid')
    Envelope=np.zeros(a_real.shape[0])
    if grid_no_grid=='grid':
        a_im,x,y,z=read_ALaDyn_bin(path,'Ienvout'+a_field_number,'grid')
        Envelope=np.sqrt(a_real**2+a_im**2)
        return (Envelope,x,y,z)
    a_im=read_ALaDyn_bin(path,'Ienvout'+a_field_number,'nogrid')
    Envelope=np.sqrt(a_real**2+a_im**2)
    return Envelope

def convert_a_in_e(a_field_number,path):

    from scipy.interpolate import interp1d
    from scipy.signal import hilbert
    import scipy.integrate as integrate

    integerparam,realparam=export_parameters(path,'Renvout'+a_field_number)

    nx=integerparam['nx_tot']
    ny=integerparam['ny_tot']
    nz=integerparam['nz_tot']
    lam_0=realparam['lambda_0']
    speed_of_light=0.3
    omega_0=2*np.pi*speed_of_light/lam_0
    k_0=2*np.pi/lam_0
    e_field=np.zeros((nx,ny,nz))


    a_real,x,y,z=read_ALaDyn_bin(path,'Renvout'+a_field_number,'grid')
    a_im=read_ALaDyn_bin(path,'Ienvout'+a_field_number,'nogrid')

    deltax=(x[-1]-x[0])/float(nx)

    a_real_diff=np.diff(a_real,axis=0)/deltax
    a_im_diff=np.diff(a_im,axis=0)/deltax
    a_real_diff=np.append(a_real_diff,[a_real_diff[-1,:,:]],axis=0)
    a_im_diff=np.append(a_im_diff,[a_im_diff[-1,:,:]],axis=0)

    for k in range(nz):
        for j in range(ny):
            discrete_cos=np.zeros(nx)
            discrete_sin=np.zeros(nx)
            for i in range(nx):
                discrete_cos[i]=(speed_of_light*a_real_diff[i,j,k]+omega_0*a_im[i,j,k])*np.cos(k_0*x[i])
                discrete_sin[i]=-(speed_of_light*a_im_diff[i,j,k]+omega_0*a_real[i,j,k])*np.sin(k_0*x[i])
            envelope=np.abs(hilbert(discrete_cos+discrete_sin))
            e_field[:,j,k]=envelope

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

    xp=np.zeros((number_timesteps,initial_tot_particles))
    yp=np.zeros((number_timesteps,initial_tot_particles))
    vxp=np.zeros((number_timesteps,initial_tot_particles))
    vyp=np.zeros((number_timesteps,initial_tot_particles))
    part_index=np.zeros((number_timesteps,initial_tot_particles),dtype=np.int)
    check_index=np.full(initial_tot_particles,1)


    if (phasespace_dimensions==7):
        zp=np.zeros((number_timesteps,initial_tot_particles))
        vzp=np.zeros((number_timesteps,initial_tot_particles))

    jump=phasespace_dimensions*number_timesteps

    file_bin=open(path+'.bin','rb')

    if(phasespace_dimensions==5):
        while True:
            try:			
                Np=struct.unpack('i', file_bin.read(4))
            except struct.error:
                break
            for k in range(0,Np[0]):
                vars=struct.unpack(jump*'f', file_bin.read(4*jump))
                for j in range(0,number_timesteps):
                    kk=int(vars[j*phasespace_dimensions+4])-1
                    xp[j,kk]=float(vars[j*phasespace_dimensions])
                    yp[j,kk]=float(vars[j*phasespace_dimensions+1])
                    vxp[j,kk]=float(vars[j*phasespace_dimensions+2])
                    vyp[j,kk]=float(vars[j*phasespace_dimensions+3])
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
                    kk=int(vars[j*phasespace_dimensions+6])-1
                    xp[j,kk]=float(vars[j*phasespace_dimensions])
                    yp[j,kk]=float(vars[j*phasespace_dimensions+1])
                    zp[j,kk]=float(vars[j*phasespace_dimensions+2])
                    vxp[j,kk]=float(vars[j*phasespace_dimensions+3])
                    vyp[j,kk]=float(vars[j*phasespace_dimensions+4])
                    vzp[j,kk]=float(vars[j*phasespace_dimensions+5])
                    part_index[j,kk]=kk+1


    # Removing "sick" particles 
    check_index[part_index[0,:]==0]=0
   # xp=xp[:,check_index!=0]
   # yp=yp[:,check_index!=0]
   # vxp=vxp[:,check_index!=0]
   # vyp=vyp[:,check_index!=0]
   # part_index=part_index[:,check_index!=0]
   # if(phasespace_dimensions==7):
   #     zp=zp[:,check_index!=0]
   #     vzp=vzp[:,check_index!=0]

   # final_tot_particles=part_index.shape[1]

   # print initial_tot_particles-final_tot_particles,'particles have been removed because they were wrongly interpreted'

    # Masking particle phasespace after they escaped the box
    #for i in range(final_tot_particles):
    #    for j in range(number_timesteps):
    #        if(part_index[j,i]==0):
    #            part_index[j:,i]=0
    #            break
    #xp=np.ma.masked_array(xp,mask=part_index==0)
    #yp=np.ma.masked_array(yp,mask=part_index==0)
    #vxp=np.ma.masked_array(vxp,mask=part_index==0)
    #vyp=np.ma.masked_array(vyp,mask=part_index==0)
    #part_index=np.ma.masked_array(part_index,mask=part_index==0)
    #if(phasespace_dimensions==7):
    #    zp=np.ma.masked_array(zp,mask=part_index==0)
    #    vzp=np.ma.masked_array(vzp,mask=part_index==0)

    #    return (xp,yp,zp,vxp,vyp,vzp,part_index)

    return (xp,yp,vxp,vyp,part_index)



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
	if ('first' in kwargs):
		case_number.append(kwargs['first'])
		plot_number+=1
	if ('second' in kwargs):
		case_number.append(kwargs['second'])
		plot_number+=1
	if ('third' in kwargs):
		case_number.append(kwargs['third'])
		plot_number+=1
	if ('time' in kwargs):
		folder_name=kwargs['time']

	if ('component' in kwargs):
		component=kwargs['component']
	
	dir_path=os.getcwd()

	if('label' in kwargs):
		plot_label=kwargs['label']

	plt.ion()
	f=plt.figure(figsize=(16,9))
	ax=f.add_subplot(111)
	colors=["#FF0066","#00FF00","#FFDE00"]
	for i in range(plot_number):

		pathenv=os.path.join(dir_path,'env'+str(case_number[i]),str(folder_name).zfill(4))
		pathpic=os.path.join(dir_path,'pic'+str(case_number[i]),str(folder_name).zfill(4))
		_,yp,vxp,vyp,_=track_data(pathpic,str(folder_name).zfill(2))
		_,ye,vxe,vye,_=track_data(pathenv,str(folder_name).zfill(2))
	
		intdata,realdata=export_parameters(pathpic,'Exfout'+str(folder_name).zfill(2))
	
		a0=round(realdata['a0'],1)
		tau=round(realdata['pulse_duration']*2.*np.arccos(1/2**0.25)/np.pi/0.3,1)
		waist=realdata['waist']
		
		if (component=='transverse'):
			ax.set_ylabel('$p_x/mc$',fontsize=18)
			if (plot_label=='a0'):
				ax.scatter(yp[0,:],vyp[-1,:],c=colors[i],s=6,label='$a_0= $'+str(a0))
				ax.scatter(ye[0,:],vye[-1,:],c='b',marker='x',s=1,alpha=0.4)
			
			if (plot_label=='tau'):

				ax.scatter(yp[0,:],vyp[-1,:],c=colors[i],s=6,label='$L_{FWHM}= $'+str(tau)+'$fs$')
				ax.scatter(ye[0,:],vye[-1,:],c='b',marker='x',s=1,alpha=0.4)
			
			if (plot_label=='wy'):

				ax.scatter(yp[0,:],vyp[-1,:],c=colors[i],s=6,label='$w_0= $'+str(waist)+'$\mu m$')
				ax.scatter(ye[0,:],vye[-1,:],c='b',marker='x',s=1,alpha=0.4)

                if (component=='longitudinal'):
			ax.set_ylabel('$p_z/mc$',fontsize=18)
                        if (plot_label=='a0'):
                                ax.scatter(yp[0,:],vxp[-1,:],c=colors[i],s=6,label='$a_0= $'+str(a0))
                                ax.scatter(ye[0,:],vxe[-1,:],c='b',marker='x',s=1,alpha=0.4)

                        if (plot_label=='tau'):

                                ax.scatter(yp[0,:],vxp[-1,:],c=colors[i],s=6,label='$L_{FWHM}= $'+str(tau)+'$fs$')
                                ax.scatter(ye[0,:],vxe[-1,:],c='b',marker='x',s=1,alpha=0.4)

                        if (plot_label=='wy'):

                                ax.scatter(yp[0,:],vxp[-1,:],c=colors[i],s=6,label='$w_0= $'+str(waist)+'$\mu m$')
                                ax.scatter(ye[0,:],vxe[-1,:],c='b',marker='x',s=1,alpha=0.4)
	if ('legend_position' in kwargs):
		leg_pos=kwargs['legend_position']
	
	ax.legend(loc=leg_pos,fontsize=18)
	ax.set_xlabel('x $(\mu m)$',fontsize=18)
