#ipython
# plot_MCGO_nc.py 
#YuP 2023

# Reading the *.nc file specified below by file=...

from numpy import *
from mpl_toolkits.mplot3d import Axes3D

from pylab import *
from matplotlib import rc 
from matplotlib.pyplot import cm,figure,axes,plot,xlabel,ylabel,title,savefig,show

import os
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import time
import pylab as pylab
import scipy.io.netcdf as nc
#matplotlib.interactive(True) # no plots on screen
matplotlib.interactive(False) # with plots on screen


#-----------------------------------------------
# NetCDF issues: machine-dependent
# Try netcdf=4 to envoke netCDF4,
# Or try netcdf=2 to work with older netCDF.
netcdf=4
#-----------------------------------------------
if netcdf==4: from netCDF4 import Dataset # YuP
#-----------------------------------------------

file='NSTX_ZOW_Maxw2_RF000_dt05_taufrac1m4_10K_t10ms_r39m.nc'



fnt  = 14   # font size for axis numbers (see 'param=' below) 
linw = 1.0  # LineWidth for linear plots
linwc= 1.0  # LineWidth for contour plots
dot_size=2. # Marker size in plots of particles as dots [pt]
Ncont= 30   # Number of contour levels
DCmn=0.0
DCmx=0.0
i_R_start=1 # Range of flux surfaces where to plot f(u,pitch)
i_R_stop =20

isave_eps=1 # To save eps format files (png are saved in any case)
iplot_ptcl_list=0 # Particles as dots in (R,Z) and in (Vpar,Vper)
            #(but only if saved, save_ptcl_list='enabled' in mcgoinput)


c=2.99792458e8 # speed of light [m/s]

#------------------------------------------------------------
#set fonts and line thicknesses
params = {
    'axes.linewidth': linw,
    'lines.linewidth': linw,
    'axes.labelsize': fnt+4,
    'text.fontsize': fnt+4,
    'legend.fontsize': fnt,
    'xtick.labelsize':fnt,
    'ytick.labelsize':fnt,
    'xtick.linewidth':linw,
    'ytick.linewidth':linw,
    'font.weight'  : 'regular',
    'format' : '%0.1e'
}
##pylab.rcParams.update(params)
#rc.defaults() #to restore defaults
mpl.rcParams['font.size']=fnt+2  # set font size for text in mesh-plots

#Input netcdf file into a structure:
if netcdf==2: 
    s_file=nc.netcdf_file(file,'r')
    if igzip==1:
        os.popen('gzip file')

if netcdf==4: 
    s_file= Dataset(file, 'r', format='NETCDF4')

print('The input file, ',file,', contains:')
print('========================================')
print("The global attributes: ",s_file.dimensions.keys())     
print("File contains variables: ",s_file.variables.keys())
print('========================================')

raxis=s_file.variables['raxis'].getValue()  #getValue() for scalar
zaxis=s_file.variables['zaxis'].getValue()  #getValue() for scalar
raxis=np.asscalar(raxis) # [m]
zaxis=np.asscalar(zaxis) # [m]

xlimiter= s_file.variables['xlimiter'] #[m] R-coord of limiter(wall)
ylimiter= s_file.variables['ylimiter'] #[m] Z-coord of limiter(wall)
xlimiter=np.asarray(xlimiter)
ylimiter=np.asarray(ylimiter)

irbnd=s_file.variables['irbnd'].getValue()  #getValue() for scalar
irbnd=np.asscalar(irbnd) #Length of R grid for distr func (at bin bndries)
ivbnd=s_file.variables['ivbnd'].getValue() 
ivbnd=np.asscalar(ivbnd) #Length of vel. grid for distr func
iptchbnd=s_file.variables['iptchbnd'].getValue() 
iptchbnd=np.asscalar(iptchbnd) #Length of pitch-angle grid for distr func

radbnd= s_file.variables['radbnd']
radbnd=np.asarray(radbnd)   #Rad. grid [m] for distr. func.
vbnd= s_file.variables['vbnd']
vbnd=np.asarray(vbnd)       #Vel. grid [m/s] for distr. func.
ptchbnd= s_file.variables['ptchbnd']
ptchbnd=np.asarray(ptchbnd) #Pitch angle grid [rad] for distr. func.

psibin= s_file.variables['psibin']
psibin=np.asarray(psibin)   #PSI pol.flux [Wb] at radbnd grid points
rho_sqpolflx= s_file.variables['rho_sqpolflx']
rho_sqpolflx=np.asarray(rho_sqpolflx) #Norm-ed rho~sqrt(pol.flux) at radbnd

vdstb= s_file.variables['vdstb']
vdstb=np.asarray(vdstb) #Midplane distr.func aver over [tim_fdist_1;tim_fdist_2]
#print 'shape of vdstb', shape(vdstb)

cosy=np.asmatrix(np.cos(ptchbnd)) #make a matrix (1,iptchbnd) {not same as vector}
cosy=cosy.transpose()          # transpose to (iptchbnd,1) shape
siny=np.asmatrix(np.sin(ptchbnd)) #make a matrix (1,iptchbnd) {not same as vector}
siny=siny.transpose()          # transpose to (iptchbnd,1) shape
xx=np.asmatrix(vbnd) # make a matrix (1,ivbnd)   [m/s]
X=np.dot(cosy,xx)/c    # (iptchbnd,ivbnd) matrix
Y=np.dot(siny,xx)/c    # (iptchbnd,ivbnd) matrix

# For evaluating density from f():
siny_dy_2pi= np.sin(ptchbnd)*np.gradient(ptchbnd)*6.28 # 2pi*sin(theta)*dtheta
vv_dv= vbnd*vbnd*np.gradient(vbnd)  # v^2 dv
siny_dy_2pi=np.asmatrix(siny_dy_2pi)
vv_dv=np.asmatrix(vv_dv)
#print 'shapes  vv_dv, siny_dy_2pi:', shape(vv_dv), shape(siny_dy_2pi)



if iplot_ptcl_list==1:
    rend= s_file.variables['rend'] #[m] R-coord at t=tend
    zend= s_file.variables['zend'] #[m] Z-coord at t=tend
    rend=np.asarray(rend)
    zend=np.asarray(zend)
    vparend= s_file.variables['vparend'] #[m/s] Vpar at t=tend
    vperend= s_file.variables['vperend'] #[m/s] Vper at t=tend
    vparend=np.asarray(vparend)
    vperend=np.asarray(vperend)


#--------------------------------------------------------------------------
if iplot_ptcl_list==1:
    fig10=plt.figure() # Particle positions over (R,Z) plane, at t=tend
    ax=plt.subplot(111) #-------------------------
    ax.set_aspect(1.0)
    #plt.axis([Zmin,Zmax,-Rmax,Rmax])
    plt.hold(True)
    plt.grid(True)
    plt.minorticks_on() # To add minor ticks
    plt.tick_params(which='both',  width=1)
    plt.tick_params(which='major', length=7)
    plt.tick_params(which='minor', length=4, color='k')
    plt.xticks(rotation=90)
    plt.xlabel('$R$  $(m)$')
    plt.ylabel('$Z$  $(m)$')
    plt.plot(xlimiter,ylimiter,'k',linewidth=linw*2)    
    plt.plot(rend,zend,'k.',markersize=dot_size)  #Large arrays; consider stride
    plt.plot(raxis,zaxis,'r+',markersize=dot_size*5)    
    savefig('mcgo_particles_RZ.png')
    if isave_eps==1: # To save eps format files (png are saved in any case)
        savefig('mcgo_particles_RZ.eps')
    show() 


#--------------------------------------------------------------------------
if iplot_ptcl_list==1:
    fig20=plt.figure() # Particles in (Vpar,Vperp), at t=tend
    ax=plt.subplot(111) #-------------------------
    ax.set_aspect(1.0)
    plt.hold(True)
    plt.grid(True)
    plt.minorticks_on() # To add minor ticks
    plt.tick_params(which='both',  width=1)
    plt.tick_params(which='major', length=7)
    plt.tick_params(which='minor', length=4, color='k')
    #plt.xticks(rotation=90)
    plt.xlabel(r'$u_{||}/c$')
    plt.ylabel(r'$u_{\perp}/c$')
    plt.plot(vparend/c,vperend/c,'k.',markersize=dot_size)  #Large arrays; consider stride
    savefig('mcgo_particles_uu.png')
    if isave_eps==1: # To save eps format files (png are saved in any case)
        savefig('mcgo_particles_uu.eps')
    show() 
    

#Make plots of distr.function at each radial point (i_R)
i_R_stop=min(i_R_stop,irbnd)
#print 'irbnd=',irbnd, shape(rho_sqpolflx),shape(vdstb)

for i_R in range(i_R_start,i_R_stop+1,1):

    # Form suffix for filename/plots :
    if i_R<10:
        i_R_index = '00'+str(i_R)
    elif i_R<100:
        i_R_index = '0'+str(i_R) 
    else:
        i_R_index = str(i_R)
    #----------------------------
    
    rho=rho_sqpolflx[i_R-1] # To add in title of each plot
    R=radbnd[i_R-1] # [m]
    
    fdist= vdstb[i_R-1,:,:]
    DDD=np.nan_to_num(np.log10(fdist))
    DDD=DDD.transpose()
    #print 'shape of DDD',shape(DDD) # Should be same as for X or Y arrays
    
    #Evaluate density [m^-3]
    #print shape(fdist), shape(vv_dv)
    dens= (vv_dv)*fdist*transpose(siny_dy_2pi)
    dens=np.asscalar(dens)
    print('i_R=',i_R,'   rho=',rho,'   R[m]=',R,'  n [m^-3]=', dens)
    
    #if DCmn==0. and DCmx==0.: 
    #    DCmin=np.min(DDD)
    #    DCmax=np.max(DDD)
    #else:
    #    DCmin = DCmn
    #    DCmax = DCmx
    DCmax=np.max(DDD)
    DCmin= DCmax-8 # limits for plots of log10(f)
    
    fig30=plt.figure()  
    ax=plt.subplot(111)
    ax.set_aspect(1.0)
    ##plt.hold(True)
    plt.grid(True)
    plt.minorticks_on() # To add minor ticks
    plt.tick_params(which='both',  width=1)
    plt.tick_params(which='major', length=7)
    plt.tick_params(which='minor', length=4, color='k')
    
    xlabel(r"$u_{||}/c$",fontsize=fnt+4)
    txt= '$log10(f)$ $at$ ' + r"$\rho=$"+r"%1.3f" %(rho) \
        + r"  $R=$"+r"%1.3f" %(R) +"$m$"    
    title(txt, fontsize=fnt+2, y=1.05)
    
    levels=np.arange(DCmin,DCmax,(DCmax-DCmin)/(Ncont-1))
    CS=plt.contour(X,Y,DDD,levels,linewidths=linwc,cmap=plt.cm.brg)
    l,b,w,h = plt.gca().get_position().bounds
    CB=plt.colorbar(orientation='vertical',shrink=0.46,\
    ticks=[-8,-6,-4,-2,0,2,4,6,8,9,10,11,12,13,14,15,16,17,18,19,20],format='%1.0f') 
    ll,bb,ww,hh = CB.ax.get_position().bounds
    CB.ax.set_position([l+w*0.04, bb, ww, hh])    
    ##CB.lines.set_linewidth(3)
    
    savefig('f'+'_r'+str(i_R_index)+'.png',dpi=600)
    