# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 18:18:56 2020

@author: Jon Staggs
"""

import numpy as np
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

pi = np.pi

nspace = 120
ivxmin = -10
ivxmax = 10
ivymin = -10
ivymax = 10
ivzmin = -10
ivzmax = 10

#Make sure stability is ensured via CFL number 
alpha_x = .05 #change in x
beta_v = .6 #max velocity magnitude
delta_t = .05 #change in t

vxmax = beta_v*float(ivxmax)
vxmin = beta_v*float(ivxmin)
vymax = beta_v*float(ivymax)
vymin = beta_v*float(ivymin)
vzmax = beta_v*float(ivzmax)
vzmin = beta_v*float(ivzmin)

ndf = 1.0 #initial gas density
u_f = -1.0 #initial gas velocity
T_f = 1.0 #initial gas temp 

nspace = 120

relax_t = 1.0 #scaled relaxation time
e_sense = 1.0 #sensible energy, 1 is placeholder 
gamma = 1.0
R = 1.0 #gas constant


Kn = 1.0 #Knusden number of gas flow, 1 is placeholder
Ma = 1.0 #mach number


#1 D bulk velocities
y_bulk =0.0
z_bulk =0.0


#CFL number for stability 
cfl = (beta_v*delta_t)/(alpha_x)

dens = np.array([0.0 for ns in range(1,122)])

#Need to make 4d array for phi, phi_1 and 3d for phi_f

phii = np.array([[[[0.0 for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)])
phi_1 = np.array([[[[0.0 for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)])
phi_2 = np.array([0.0 for ns in range(1,122)])
phi_f = np.array([[[[float(ns) for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)])
#phi_t = np.array([[[[[t for t in range(0,500)]for ns in range(0,121)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)])
del_phi = [[[[0.0 for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)]


tauxx = np.array([0.0 for ns in range(0,nspace+1)])
tauxy = np.array([0.0 for ns in range(0,nspace+1)])

qx = np.array([0.0 for ns in range(0,nspace+1)])


Tcal = np.array([0.0 for ns in range(0,nspace+1)])

e_s = np.array([0.0 for ns in range(1,122)])

x_flux = np.array([[[[0.0 for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)])
y_flux = np.array([[[[0.0 for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)])
z_flux = np.array([[[[0.0 for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)])


T_denom = .5/T_f
norm = ndf/((2.0*pi*T_f)**(3/2))


def gen_np_array(num):
    return np.array([0.0 for ns in range(0, num + 1)])

def gen_array(num):
    return [0.0 for i in range(1, num)]


def maxwell():
    phi_f = np.array([[[[0.0 for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)])
    norm = ndf/((2.0*pi*T_f)**(3/2))
    for ns in range(1,121):
        for k in range(0,21):
            for j in range(0,21):
                for i in range(0,21):
                    phi_f[i][j][k][ns] = norm*np.exp(-(((beta_v*float(i-10)))**2+(beta_v*float(j-10))**2+(beta_v*float(k-10))**2)*T_denom)
    return(phi_f)



def convect1(phi):
    #first and last points in space are determined by boundary conditions
    for ns in range(2,120): #1 to 119
        for k in range(0,21): #-10 to 10
            for j in range(0,21): #-10 to 10
                for i in range(0,10): #-10 to -1
                    del_phi[i][j][k][ns] = cfl*float(i-10)*(phi[i][j][k][ns+1]-phi[i][j][k][ns])
                    phi_1[i][j][k][ns]=phi[i][j][k][ns]-del_phi[i][j][k][ns]
                    #No convection for =0
                phi_1[10][j][k][ns] = phi[10][j][k][ns]
                #for > 0 upwind is on the right
                for i in range(11,21): #1 to 10
                    del_phi[i][j][k][ns] = cfl*float(i-10)*(phi[i][j][k][nspace]-phi[i][j][k][nspace-1])
                    phi_1[i][j][k][ns]=phi[i][j][k][ns]-del_phi[i][j][k][ns]
    #left side boundary condition - for i >0 specular reflection of incoming (i<0) phi 
    #Normal convective update for i<0
    #we overwrite phi values at first space point for i>0 using -i values
    for k in range(0,21): #-10 to 10
        for j in range(0,21): #-10 to 10
            for i in range(0,11): #-10 to 0
                #for i<0 upwind is on the right (ns+1)
                del_phi[i][j][k][1] = cfl*float(i-10)*(phi[i][j][k][2]-phi[i][j][k][1])
                phi_1[i][j][k][1]= phi[i][j][k][1]-del_phi[i][j][k][1]
            for i in range(11,21): #0 to 10
                phi_1[i][j][k][1] = phi[i-(2*(i-10))][j][k][1]
    #Right side boundary condition - freestream inflow for i<0
    #Normal convection update for i >0
    for k in range(0,21): #-10 to 10
        for j in range(0,21): #-10 to 10
            for i in range(0,10): #-10 to -1
                phi_1[i][j][k][nspace] = phi_f[i][j][k][ns]
            #for i > 0 upwind is on the left (ns-1)
            for i in range(11,21): #1 to 10
                del_phi[i][j][k][nspace] =cfl*float(i-10)*(phi[i][j][k][nspace]-phi[i][j][k][nspace-1])
                phi_1[i][j][k][nspace] = phi[i][j][k][nspace]-del_phi[i][j][k][nspace]
    for ns in range(1,nspace+1):
        for k in range(0,21): #-10 to 10
            for j in range(0,21): #-10 to 10
                for i in range(0,10): #-10 to -1 
                    phi[i][j][k][ns]=phi_1[i][j][k][ns]
                #No update for i =0 since there is no convection
                for i in range(11,21): #1 to 10
                    phi[i][j][k][ns] = phi_1[i][j][k][ns]
    return np.copy(phi)

# phi_0 = maxwell()
# phi_t=[]

# for t in range(0,50):
#     phi_0 = convect1(phi_0)
#     print("is it working?")
#     ax = plt.axes(projection="3d")
#     x = [i for i in range(-10,11)]
#     y = [j for j in range(-10,11)]
#     X,Y = np.meshgrid(x,y)
#     phitest = []
#     for i in range(0,21):
#         for j in range(0,21):
#             phitest.append(phi_0[10][i][j][60])
#     print("is it working?")
#     phitest2 = np.reshape(phitest,(21,21))
    
#     fig = plt.figure()
#     fig.tight_layout() 
#     ax.plot_wireframe(X,Y,phitest2)
#     ax.set_xlabel('v_y')
#     ax.set_ylabel('v_z')
#     ax.set_zlabel('phi({})'.format(t))
#     plt.tight_layout(pad=40)
#     plt.show()    
    
#     phi2test =[]
#     for ns in range(0,120):
#         phi2test.append(phi_0[10][10][10][ns])
    
#     plt.plot(phi2test)
#     plt.ylabel("phi({}+1)".format(t))
#     plt.xlabel("v_x")
#     plt.ylim(bottom = 0)
#     plt.ylim(top= 0.09)
#     plt.show()    



def halfconvect(phi):
    #first and last points in space are determined by boundary conditions
    for ns in range(2,120): #1 to 119
        for k in range(0,21): #-10 to 10
            for j in range(0,21): #-10 to 10
                for i in range(0,10): #-10 to -1
                    del_phi[i][j][k][ns] = (cfl * 0.5) *float(i-10)*(phi[i][j][k][ns+1]-phi[i][j][k][ns])
                    phi_1[i][j][k][ns]=phi[i][j][k][ns]-del_phi[i][j][k][ns]
                    #No convection for =0
                phi_1[10][j][k][ns] = phi[10][j][k][ns]
                #for > 0 upwind is on the right
                for i in range(11,21): #1 to 10
                    del_phi[i][j][k][ns] = (cfl * 0.5) *float(i-10)*(phi[i][j][k][nspace]-phi[i][j][k][nspace-1])
                    phi_1[i][j][k][ns]=phi[i][j][k][ns]-del_phi[i][j][k][ns]
    #left side boundary condition - for i >0 specular reflection of incoming (i<0) phi 
    #Normal convective update for i<0
    #we overwrite phi values at first space point for i>0 using -i values
    for k in range(0,21): #-10 to 10
        for j in range(0,21): #-10 to 10
            for i in range(0,11): #-10 to 0
                #for i<0 upwind is on the right (ns+1)
                del_phi[i][j][k][1] = (cfl * 0.5)*float(i-10)*(phi[i][j][k][2]-phi[i][j][k][1])
                phi_1[i][j][k][1]= phi[i][j][k][1]-del_phi[i][j][k][1]
            for i in range(11,21): #0 to 10
                phi_1[i][j][k][1] = phi[i-(2*(i-10))][j][k][1]
    #Right side boundary condition - freestream inflow for i<0
    #Normal convection update for i >0
    for k in range(0,21): #-10 to 10
        for j in range(0,21): #-10 to 10
            for i in range(0,10): #-10 to -1
                phi_1[i][j][k][nspace] = phi_f[i][j][k][ns]
            #for i > 0 upwind is on the left (ns-1)
            for i in range(11,21): #1 to 10
                del_phi[i][j][k][nspace] = (cfl * 0.5)*float(i-10)*(phi[i][j][k][nspace]-phi[i][j][k][nspace-1])
                phi_1[i][j][k][nspace] = phi[i][j][k][nspace]-del_phi[i][j][k][nspace]
    for ns in range(1,nspace+1):
        for k in range(0,21): #-10 to 10
            for j in range(0,21): #-10 to 10
                for i in range(0,10): #-10 to -1 
                    phi[i][j][k][ns]=phi_1[i][j][k][ns]
                #No update for i =0 since there is no convection
                for i in range(11,21): #1 to 10
                    phi[i][j][k][ns] = phi_1[i][j][k][ns]          
    return np.copy(phi)

                    
    
def f_density(phi):
    dens = gen_np_array(120)
    for ns in range(1,121):
        for k in range(0,21):
            for j in range(0,21):
                for i in range(0,21):
                    dens[ns] = dens[ns]+ phi[i][j][k][ns]
        dens[ns]=dens[ns]*beta_v**3
    return(dens)



    
def f_x_velocity(phi,density):
    x_vel = gen_np_array(120)
    for ns in range(1,121):
        for k in range(0,21):
            for j in range(0,21):
                for i in range(0,21):
                    x_vel[ns] += phi[i][j][k][ns]*float(beta_v)*float(i-10)
        x_vel[ns]=x_vel[ns]/density[ns]*beta_v**3 
    return(x_vel)
    


def f_y_velocity(phi,density):
    y_vel = gen_np_array(120)
    for ns in range(1,121):
        for k in range(0,21):
            for j in range(0,21):
                for i in range(0,21):
                    y_vel[ns] += phi[i][j][k][ns]*float(beta_v)*float(j-10)
        y_vel[ns]=y_vel[ns]/density[ns]*beta_v**3
    return(y_vel)
    
    
    
def f_z_velocity(phi,density):
    z_vel = gen_np_array(120)
    for ns in range(1,121):
        for k in range(0,21):
            for j in range(0,21):
                for i in range(0,21):
                    z_vel[ns] += phi[i][j][k][ns]*float(beta_v)*float(k-10)
        z_vel[ns]=z_vel[ns]/density[ns]*beta_v**3
    return(z_vel)
    

    
    
def f_visc_stress_xx(phi,beta_v,x_vel,y_vel,z_vel):
    tauxx = gen_np_array(120)
    for ns in range(1,121):
        for k in range(0,21):
            for j in range(0,21):
                for i in range(0,21):
                    Cx = (beta_v*float(i-10)-x_vel[ns])
                    Cy = (float(j-10)*beta_v-y_vel[ns])
                    Cz = (float(k-10)*beta_v-z_vel[ns])
                    Csq = Cx**2 + Cy**2 + Cz**2
                    tauxx[ns] += phi[i][j][k][ns]*(-Cx*Cx+Csq/3.0)
        tauxx[ns] = tauxx[ns]*(beta_v**3)
    return(tauxx)


                                 
def f_visc_stress_xy(phi,x_vel,y_vel):
    tauxy = gen_np_array(120)
    for ns in range(1,121):
        for k in range(0,21):
            for j in range(0,21):
                for i in range(0,21):
                    Cx = (float(i-10)*beta_v-x_vel[ns])
                    Cy = (float(j-10)*beta_v-y_vel[ns])
                    tauxy[ns]+=phi[i][j][k][ns]*(-Cx*Cy)
        tauxy[ns] =tauxy[ns]*(beta_v**3)
    return(tauxy)

                 
    
def f_temp(phi,x_velocity,density):
    Tcal = gen_np_array(120)
    for ns in range(1,121):
        for k in range(0,21):
            for j in range(0,21):
                for i in range(0,21):
                    Cxsq = (beta_v*float(i-10)-x_velocity[ns])**2
                    Cysq = float(j-10)*float(j-10)*beta_v**2
                    Czsq = float(k-10)*float(k-10)*beta_v**2
                    Csq = Cxsq + Cysq + Czsq
                    Tcal[ns]+=phi[i][j][k][ns]*Csq
        Tcal[ns]=Tcal[ns]*(beta_v**3)/(3.0*dens[ns])
    return(Tcal)

    
def f_temp_11(phi,x_velocity,dens):
    Tcal = gen_np_array(120)
    for ns in range(1,121):
        for k in range(0,21):
            for j in range(0,21):
                for i in range(0,21):
                    Cxsq = (beta_v*float(i-10)-x_velocity[ns])**2
                    Tcal[ns] += phi[i][j][k][ns]*Cxsq
        Tcal[ns]=Tcal[ns]*(beta_v**3)/(dens[ns])
    return(Tcal)
    
        
def f_temp_22(phi,dens):
    Tcal = gen_np_array(120)
    for ns in range(1,121):
        for k in range(0,21):
            for i in range(0,21):
                for j in range(0,21):
                    Cysq = (beta_v*float(j-10))**2
                    Tcal[ns] += phi[i][j][k][ns]*Cysq
        Tcal[ns]=Tcal[ns]*(beta_v**3)/(dens[ns])
    return(Tcal)
    

    
def f_temp_33(phi,dens):
    Tcal = gen_np_array(120)
    for ns in range(1,121):
        for i in range(0,21):
            for j in range(0,21):
                for k in range(0,21):
                    Czsq = (beta_v*float(k-10))**2
                    Tcal[ns] += phi[i][j][k][ns]*Czsq
        Tcal[ns]=Tcal[ns]*(beta_v**3)/(dens[ns])
    return(Tcal)
        


def f_heat_flux(phi,x_vel,y_vel,z_vel):
    qx = gen_np_array(120)
    for ns in range(0,121):
        for k in range(0,21):
            for j in range(0,21):
                for i in range(0,21):
                    Cx = (beta_v*float(i-10)-x_vel[ns])
                    Cy = (beta_v*float(j-10)-y_vel[ns])
                    Cz = (beta_v*float(k-10)-z_vel[ns])
                    Csq = Cx**2+Cy**2+Cz**2
                    qx[ns]+= phi[i][j][k][ns]*Csq
        qx[ns]=qx[ns]*0.5*beta_v**3
    return(qx)

    
def f_del_v_x(phi,x_velocity):
    del_x_flux = np.array([[[[0.0 for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)])
    for ns in range(1,121):
        for k in range(0,21):
            for j in range(0,21):
                #Boundary conditions in x-velocity
                for i in range(1,20):
                    del_x_flux[i][j][k][ns] = ((float(i-10)*beta_v - x_velocity[ns])*(phi[i+1][j][k][ns]-phi[i-1][j][k][ns])/(2*beta_v)+phi[i][j][k][ns])*delta_t
                #Left Side Boundary Conditions
                del_x_flux[0][j][k][ns] = ((float(i-10)*beta_v -x_velocity[ns])*(phi[1][j][k][ns])/(2*beta_v) +phi[0][j][k][ns])*delta_t
                #Right Side Boundary Conditions
                del_x_flux[20][j][k][ns] = ((float(i-10)*beta_v-x_velocity[ns])*((-1)*phi[20][j][k][ns])/(2*beta_v) + phi[20][j][k][ns])*delta_t
    return(del_x_flux)

def f_del_v_y(phi,y_velocity):
    del_y_flux = np.array([[[[0.0 for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)])
    for ns in range(1,121):
        for k in range(0,21):
            for i in range(0,21):
                #Boundary conditions in x-velocity
                for j in range(1,20):
                    del_y_flux[i][j][k][ns] = ((float(j-10))*beta_v*(phi[i][j+1][k][ns]-phi[i][j-1][k][ns])/(2*beta_v)+phi[i][j][k][ns])*delta_t
                #Left Side Boundary Conditions 
                del_y_flux[i][0][k][ns] = ((float(j-10))*beta_v*(phi[i][1][k][ns])/(beta_v*2) +phi[i][0][k][ns])*delta_t
                #Right Side Boundary Conditions
                del_y_flux[i][20][k][ns] = (float(j-10))*beta_v*(((-1)*phi[i][19][k][ns])/(beta_v*2) + phi[i][20][k][ns])*delta_t
    return(del_y_flux)

    
def f_del_v_z(phi,z_velocity):
    del_z_flux = np.array([[[[0.0 for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)])
    for ns in range(1,121):
        for i in range(0,21):
            for j in range(0,21):
                #Boundary conditions in x-velocity
                for k in range(1,20):
                    del_z_flux[i][j][k][ns] = ((float(k-10))*beta_v*(phi[i][j][k+1][ns]-phi[i][j][k-1][ns])/(beta_v*2)+phi[i][j][k][ns])*delta_t
                #Left Side Boundary Conditions
                del_z_flux[i][j][0][ns] = ((float(k-10))*beta_v*(phi[i][j][1][ns])/(beta_v*2) +phi[i][j][0][ns])*delta_t
                #Right Side Boundary Conditions
                del_z_flux[i][j][20][ns] = ((float(k-10))*beta_v*((-1)*phi[i][j][19][ns])/(beta_v*2) + phi[i][j][20][ns])*delta_t
    return(del_z_flux)

    
def f_v_flux(del_v_x,del_v_y,del_v_z):
    v_flux_1 = np.array([[[[0.0 for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)])
    for ns in range(1,121):
        for k in range(0,21):
            for j in range(0,21):
                for i in range(0,21):
                    v_flux_1[i][j][k][ns] = (del_v_x[i][j][k][ns]+del_v_y[i][j][k][ns]+del_v_z[i][j][k][ns])
    return(v_flux_1)

    
    
def f_double_del_v_x(phi):
    double_del_phi_v_x_1 = np.array([[[[0.0 for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)])
    for ns in range(1,nspace+1):
        for k in range(0,21):
            for j in range(0,21):
                for i in range(1,20):
                    double_del_phi_v_x_1[i][j][k][ns]=((phi[i+1][j][k][ns]-2*phi[i][j][k][ns]+phi[i-1][j][k][ns])/(beta_v**2))*delta_t
                #Leftside boundary conditions
                double_del_phi_v_x_1[0][j][k][ns] = ((phi[1][j][k][ns]-2*phi[0][j][k][ns])/(beta_v**2))*delta_t
                #Rightside boundary conditions
                double_del_phi_v_x_1[20][j][k][ns] = (((-2)*phi[20][j][k][ns]+phi[19][j][k][ns])/(beta_v**2))*delta_t
    return(double_del_phi_v_x_1)

    
def f_double_del_v_y(phi):
    double_del_phi_v_y_1 = np.array([[[[0.0 for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)])
    for ns in range(1,nspace+1):
        for k in range(0,21):
            for i in range(0,21):
                for j in range(1,20):
                    double_del_phi_v_y_1[i][j][k][ns]=((phi[i][j+1][k][ns]-2*phi[i][j][k][ns]+phi[i][j-1][k][ns])/(beta_v**2))*delta_t
                #Leftside boundary conditions
                double_del_phi_v_y_1[i][0][k][ns] = ((phi[i][1][k][ns]-2*phi[i][0][k][ns])/(beta_v**2))*delta_t
                #Rightside boundary conditions
                double_del_phi_v_y_1[i][20][k][ns] = (((-2)*phi[i][20][k][ns]+phi[i][19][k][ns])/(beta_v**2))*delta_t
    return(double_del_phi_v_y_1)
    
    
def f_double_del_v_z(phi):
    double_del_phi_v_z_1 = np.array([[[[0.0 for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)])
    for ns in range(1,nspace+1):
        for i in range(0,21):
            for j in range(0,21):
                for k in range(1,20):
                    double_del_phi_v_z_1[i][j][k][ns]=((phi[i][j][k+1][ns]-2*phi[i][j][k][ns]+phi[i][j][k-1][ns])/(beta_v**2))*delta_t
                #Leftside boundary conditions
                double_del_phi_v_z_1[i][j][0][ns] = ((phi[i][j][1][ns]-2*phi[i][j][0][ns])/(beta_v**2))*delta_t
                #Rightside boundary conditions
                double_del_phi_v_z_1[i][j][20][ns] = (((-2)*phi[i][j][20][ns]+phi[i][j][19][ns])/(beta_v**2))*delta_t
    return(double_del_phi_v_z_1)

    
    
def f_laplacian_v(phi,temp_11,temp_22,temp_33,double_del_v_x,double_del_v_y,double_del_v_z):
    laplacian_v_1 = [[[[0.0 for ns in range(1,122)]for k in range(0,21)]for j in range(0,21)]for i in range(0,21)]
    for ns in range(1,nspace+1):
        for k in range(0,21):
            for j in range(0,21):
                for i in range(0,21): #should these be + ?
                    laplacian_v_1[i][j][k][ns] = (temp_11[ns]*double_del_v_x[i][j][k][ns])+(temp_22[ns]*double_del_v_y[i][j][k][ns])+(temp_33[ns]*double_del_v_z[i][j][k][ns])
    return(laplacian_v_1)
    
def f_collision(phi,v_flux,laplacian):
    collision_1 = [[[[0.0 for ns in range(1,122)] for k in range(0,21)] for j in range(0,21)] for i in range(0,21)]
    for ns in range(1,121):
        for k in range(0,21):
            for j in range(0,21):
                for i in range(0,21):
                    collision_1[i][j][k][ns]= (1/relax_t)*(1/Kn)*(1/Ma)*(v_flux[i][j][k][ns]+laplacian[i][j][k][ns])
                    phi[i][j][k][ns] += collision_1[i][j][k][ns] 
    return np.copy(phi)
    
       

def main():
    phi = maxwell()
    
    #dog test confirms convect1 is working roughly as intended 
    # dog = convect1(phi)
    
    # dogtest = []
    # for ns in range(0,21):
    #     dogtest.append(dog[ns][10][10][60])

    # plt.plot(dogtest)
    # plt.ylabel("convect(phi)")
    # plt.xlabel("v_x")
    # plt.ylim(bottom = 0)
    # plt.ylim(top= 0.09)
    # plt.show()     
    
    
    ### Plotting ####
   
    phitest3 =[]
    for ns in range(0,21):
        phitest3.append(phi[ns][10][10][60])
    
    plt.plot(phitest3)
    plt.ylabel("phi_0")
    plt.xlabel("v_x")
    plt.ylim(bottom = 0)
    plt.ylim(top= 0.09)
    plt.show()        
    
    ##################
    
    for t in range(0, 100):
        halfconvect(phi)
        halfphi = halfconvect(phi)
        
        ### Plotting ####
        
        phitest4 =[]
        for ns in range(0,21):
            phitest4.append(halfphi[ns][10][10][60])
        
        plt.plot(phitest4)
        plt.ylabel("phi({})".format(t))
        plt.xlabel("v_x")
        plt.ylim(bottom = 0)
        plt.ylim(top= 0.09)
        plt.show()          
        
        ################
        
        ### Printing test values ####
  
        density = f_density(halfconvect(phi))
        print(t, "is it working?")
        x_vel = f_x_velocity(halfconvect(phi), density)
        y_vel = f_y_velocity(halfconvect(phi), density)
        z_vel = f_z_velocity(halfconvect(phi), density)
        temp_11 = f_temp_11(halfconvect(phi), x_vel, density)
        temp_22 = f_temp_22(halfconvect(phi), density)
        temp_33 = f_temp_33(halfconvect(phi), density)
        dd_x = f_double_del_v_x(halfconvect(phi))
        dd_y = f_double_del_v_y(halfconvect(phi))
        dd_z = f_double_del_v_z(halfconvect(phi))
        flux = f_v_flux(f_del_v_x(halfconvect(phi), x_vel), f_del_v_y(halfconvect(phi), y_vel), f_del_v_z(halfconvect(phi),z_vel))
        laplace = f_laplacian_v(halfconvect(phi), temp_11, temp_22, temp_33, dd_x, dd_y, dd_z)
        collide = f_collision(halfconvect(phi), flux, laplace)

        density = f_density(collide)
        print(density, "big this is density at t=",t)
        x_vel = f_x_velocity(collide, density)
        print(x_vel,"big this is  vel at t=",t)
        y_vel = f_y_velocity(collide, density)
        z_vel = f_z_velocity(collide, density)
        temp_11 = f_temp_11(collide, x_vel, density)
        print(temp_11,"big this is temp11 t=",t)
        temp_22 = f_temp_22(collide, density)
        print(temp_22,"big this is temp22 t=",t)
        temp_33 = f_temp_33(collide, density)  
        print(temp_33,"big this is temp33 t=",t)
        
        
        ### Plotting ####        
        phitest5 =[]
        for ns in range(0,21):
            phitest5.append(collide[ns][10][10][60])
        
        plt.plot(phitest5)
        plt.ylabel("phi_1*")
        plt.xlabel("v_x")
        plt.ylim(bottom = 0)
        plt.ylim(top= 0.09)
        plt.show()  
        
        ###############
          
        for mini_t in range(0, 8):
            phi_star = convect1(collide)
           
            ### Plotting ####
            
            phi_star_test = []
            for ns in range(0,21):
                phi_star_test.append(phi_star[ns][10][10][60])
             
            plt.plot(phi_star_test)
            plt.ylabel("phi 1.5*({})".format(mini_t))
            plt.xlabel("v_x")
            plt.ylim(bottom = 0)
            plt.ylim(top= 0.09)
            plt.show()  
          
            ###############  

            ### Printing test values #### 
            
            density = f_density(phi_star)
            print(density, "mini this is density at t=",mini_t)
            x_vel = f_x_velocity(phi_star, density)
            print(x_vel,"mini this is  vel at t=",mini_t)
            y_vel = f_y_velocity(phi_star, density)
            z_vel = f_z_velocity(phi_star, density)
            temp_11 = f_temp_11(phi_star, x_vel, density)
            print(temp_11,"mini this is temp11 t=",mini_t)
            temp_22 = f_temp_22(phi_star, density)
            print(temp_22,"mini this is temp22 t=",mini_t)
            temp_33 = f_temp_33(phi_star, density)
            print(temp_33,"mini this is temp33 t=",mini_t)
            
            ############################
            
            dd_x = f_double_del_v_x(phi_star)
            dd_y = f_double_del_v_y(phi_star)
            dd_z = f_double_del_v_z(phi_star)
            flux = f_v_flux(f_del_v_x(phi_star, x_vel), f_del_v_y(phi_star, y_vel), f_del_v_z(phi_star,z_vel))
            laplace = f_laplacian_v(phi_star, temp_11, temp_22, temp_33, dd_x, dd_y, dd_z)
            collide = f_collision(phi_star, flux, laplace)
            
            ### Plotting ####
            
            colltest =[]
            for ns in range(0,21):
                colltest.append(collide[ns][10][10][60])
            print("phi(",mini_t+2.5,"0) does it look whack?")
            plt.plot(colltest)
            plt.ylabel("phi({}".format(mini_t+2.5))
            plt.xlabel("v_x")
            plt.ylim(bottom = 0)
            plt.ylim(top= 0.09)
            plt.show()           
            
            ################
            
            
        phi = halfconvect(collide)
        x_vel = f_x_velocity(phi, density)
        y_vel = f_y_velocity(phi, density)
        z_vel = f_z_velocity(phi, density)
        temp_11 = f_temp_11(phi, x_vel, density)
        temp_22 = f_temp_22(phi, density)
        temp_33 = f_temp_33(phi, density)
        
        
        ### Plotting ####

        ax = plt.axes(projection="3d")
        x = [i for i in range(-10,11)]
        y = [j for j in range(-10,11)]
        X,Y = np.meshgrid(x,y)
        conv2test = []
        for i in range(0,21):
            for j in range(0,21):
                conv2test.append(phi[10][i][j][60])
        
        conv2test2 = np.reshape(conv2test,(21,21))
        
        fig = plt.figure()
        fig.tight_layout() 
        ax.plot_wireframe(X,Y,conv2test2)
        ax.set_xlabel('v_y')
        ax.set_ylabel('v_z')
        ax.set_zlabel('phi({}+1)'.format(time))
        plt.tight_layout(pad=40)
        plt.show()    
        
        conv2t =[]
        for ns in range(0,21):
            conv2t.append(phi[ns][10][10][60])
        
        plt.plot(conv2t)
        plt.ylabel("phi({}+1)".format(time))
        plt.xlabel("v_x")
        plt.ylim(bottom = 0)
        plt.ylim(top= 0.09)
        plt.show()    
        
        ####################

main()