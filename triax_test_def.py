'''
Created on Jul 22, 2016

@author: benni
'''

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import os

from IPsimulator import umat_setup
from IPsimulator import loading
from IPsimulator import constants

from objective_functions_tsaiwu import strains_lat

vec_strains_lat = np.vectorize(strains_lat)

#----------------------------------------------------------------------------------
# DEFINITION AND INITIALIZATION OF VARIABLES
#----------------------------------------------------------------------------------

umat_aux_vars = umat_setup.set_variables(mat_file_name='test-material-tsaiwu.inp')    # load material data and initialize auxiliary variables for UMAT

output_dict = {}                        # empty dictionary for matlab output file
np.set_printoptions(linewidth=100)      # reset linewidth of numpy print command
size_pres = (9, 5)                      # width, height in inches suitable for pp slides
plt.figure(1, size_pres)                # create first figure for plots



# #---------------------------------------------------------------
# # uniaxial compressive test
# #---------------------------------------------------------------
#         
# # number of increments
# n_inc = 300
# sig_r = 0.
# delta_eps_a = -5.0e-5
#   
# beta = np.pi/2.                           # beta describes the angle between the global x1-axis (axial loading direction) and the direction of the foliation planes    
# rot_angles = [0.0, beta, 0.0]       # definition of rotation angles alpha, beta, gamma of material orientation
#        
# sig1 = np.empty((n_inc+1,6))
# eps1 = np.empty_like(sig1)
# statev1 = np.empty((n_inc+1,umat_aux_vars['nstatv']))
# outer_it1 = np.empty((n_inc+1,))
# inner_it1 = np.empty((n_inc+1,constants.max_it))
#        
# sig1[:,:], eps1[:,:], statev1[:,:], outer_it1[:], inner_it1[:,:] = loading.triax_ct_simulation(n_inc, sig_r, delta_eps_a,
#                                                                                              umat_aux_vars, rot_angles)
# 
# # label
# lab = '$\sigma_r=0$ MPa'
# plt.plot(eps1[:,0], sig1[:,0], color='b', marker='', label=lab)     # plot axial strain vs axial stress
# plt.plot(eps1[:,1], sig1[:,0], color='b', marker='')                           # plot radial strain vs axial stress curve 
# 
# # print outer_it1
# # print outer_it1[n_inc]
# # print inner_it1[n_inc,:]
# # print eps1   



#---------------------------------------------------------------
# triax compression test (sig_r=-12.5 MPa)
#---------------------------------------------------------------


# HORIZONTAL FOLIATION PLANES:
# ----------------------------
# number of increments
# n_inc = 100
# sig_r = -25.0
# delta_eps_a = -5.0e-5
n_inc = 45
sig_r = -12.5
delta_eps_a = -1e-4
   
beta = np.pi/180. * 0.                 # beta describes the angle between the global x1-axis (axial loading direction) and the direction of the foliation planes
rot_angles = [0.0, beta, 0.0]        # definition of rotation angles alpha, beta, gamma of material orientation

# create variables for results        
sig3 = np.empty((n_inc+1,6))
eps3 = np.empty_like(sig3)
statev3 = np.empty((n_inc+1,umat_aux_vars['nstatv']))
outer_it3 = np.empty((n_inc+1,))
inner_it3 = np.empty((n_inc+1,constants.max_it))

# start simulation        
sig3[:,:], eps3[:,:], statev3[:,:], outer_it3[:], inner_it3[:,:] = loading.triax_ct_simulation(n_inc, sig_r, delta_eps_a,
                                                                                             umat_aux_vars, rot_angles)

# compute lateral (circumferential) strains
eps_lat3 = vec_strains_lat(eps3[:,0], eps3[:,1], 35.)
 
# label
lab = r'$\sigma_r=-12.5$ MPa, $\alpha=0$'

# plot simulation results
plt.plot(eps3[1:,2]-eps3[1,2], sig3[1:,2], color='g', marker='', label=lab)     # plot axial strain vs axial stress
# plt.plot(eps3[1:,1]-eps3[1,1], sig3[1:,2], color='g')                           # plot radial strain vs axial stress curve 
plt.plot(eps_lat3[1:]-eps_lat3[1], sig3[1:,2], color='g')                           # plot radial strain vs axial stress curve 
 
# # load and plot experimental data
# expData = np.loadtxt('parameter_identification/results/triax_res_exp06p_hardening.txt')
# plt.plot(expData[:,2], expData[:,0]-25., color='g', linestyle='--')
# plt.plot(expData[:,3], expData[:,0]-25., color='g', linestyle='--')
        
# print outer_it3
# print outer_it3[n_inc]
# print inner_it3[n_inc,:]
# print statev3[1,:]



# INCLINED FOLIATION PLANES:
# ----------------------------

n_inc = 20

# redefine beta and start simulation again
beta = np.pi/180. * 60. 
rot_angles = [0.0, beta, 0.0]
 
# create variables for results        
sig4 = np.empty((n_inc+1,6))
eps4 = np.empty_like(sig4)
statev4 = np.empty((n_inc+1,umat_aux_vars['nstatv']))
outer_it4 = np.empty((n_inc+1,))
inner_it4 = np.empty((n_inc+1,constants.max_it))
 
# start simulation        
sig4[:,:], eps4[:,:], statev4[:,:], outer_it4[:], inner_it4[:,:] = loading.triax_ct_simulation(n_inc, sig_r, delta_eps_a,
                                                                                             umat_aux_vars, rot_angles)

# compute lateral (circumferential) strains
eps_lat4 = vec_strains_lat(eps4[:,0], eps4[:,1], 35.)
 
# label
lab = r'$\sigma_r=-12.5$ MPa, $\alpha=60$'

# plot simulation results
plt.plot(eps4[1:,2]-eps4[1,2], sig4[1:,2], color='r', marker='', label=lab)     # plot axial strain vs axial stress
# plt.plot(eps4[1:,1]-eps4[1,1], sig4[1:,2], color='r')                           # plot radial strain vs axial stress curve 
plt.plot(eps_lat4[1:]-eps_lat4[1], sig4[1:,2], color='r')                     # plot radial strain vs axial stress curve 
 
 
 
# VERTICAL FOLIATION PLANES:
# ----------------------------

n_inc = 27

# redefine beta and start simulation again
beta = np.pi/180. * 90. 
rot_angles = [0.0, beta, 0.0]
 
# create variables for results        
sig5 = np.empty((n_inc+1,6))
eps5 = np.empty_like(sig5)
statev5 = np.empty((n_inc+1,umat_aux_vars['nstatv']))
outer_it5 = np.empty((n_inc+1,))
inner_it5 = np.empty((n_inc+1,constants.max_it))
 
# start simulation        
sig5[:,:], eps5[:,:], statev5[:,:], outer_it5[:], inner_it5[:,:] = loading.triax_ct_simulation(n_inc, sig_r, delta_eps_a,
                                                                                             umat_aux_vars, rot_angles)
 

# compute lateral (circumferential) strains
eps_lat5 = vec_strains_lat(eps5[:,0], eps5[:,1], 35.)

# label
lab = r'$\sigma_r=-12.5$ MPa, $\alpha=90$'

# plot simulation results
plt.plot(eps5[1:,2]-eps5[1,2], sig5[1:,2], color='b', marker='', label=lab)     # plot axial strain vs axial stress
# plt.plot(eps5[1:,1]-eps5[1,1], sig5[1:,2], color='b')                           # plot radial strain vs axial stress curve 
plt.plot(eps_lat5[1:]-eps_lat5[1], sig5[1:,2], color='b')                           # plot radial strain vs axial stress curve 


#---------------------------------------------------------------
# load experimental data and plot it
#---------------------------------------------------------------

# horizontal identification data
exp04 = np.loadtxt('triax_res_exp04p_hardening.txt')
plt.plot(exp04[:,0], exp04[:,2], color='g', linestyle='--')
plt.plot(exp04[:,1], exp04[:,2], color='g', linestyle='--')

# inclined identification data
exp30 = np.loadtxt('triax_res_exp30p_hardening_first_peak.txt')
plt.plot(exp30[:,0], exp30[:,2], color='r', linestyle='--')
plt.plot(exp30[:,1], exp30[:,2], color='r', linestyle='--')

# vertical identification data
exp15 = np.loadtxt('triax_res_exp15p_hardening.txt')
plt.plot(exp15[:,0], exp15[:,2], color='b', linestyle='--')
plt.plot(exp15[:,1], exp15[:,2], color='b', linestyle='--')


#---------------------------------------------------------------
# plot options for strain stress curve in figure(1)
#---------------------------------------------------------------
title = 'BBT-QuartzPhyllite'
plt.ticklabel_format(style='sci', axis='x', scilimits=(-2,2))
plt.title(title, fontsize=14)
plt.xlabel('$\epsilon_{a}, \epsilon_{r}$ [-]')
plt.ylabel('$\sigma_{a}$ [N/mm$^2$]')
plt.tick_params(labelsize=12)
plt.legend()
plt.grid(True)
plt.tight_layout()

# fig_name = '/home/benni/Documents/Triax-Simulations/As-%.1f_mg0-%.1f.pdf' % (mat_data[17], mat_data[11])
# plt.savefig(fig_name)


#---------------------------------------------------------------
# save simulation data
#---------------------------------------------------------------
np.savetxt('triax-sim-res-hori-125-tsaiwu.txt', np.column_stack((eps3[1:,2]-eps3[1,2], 
                                                                 eps_lat3[1:]-eps_lat3[1], 
                                                                 sig3[1:,2])))
np.savetxt('triax-sim-res-incl-125-tsaiwu.txt', np.column_stack((eps4[1:,2]-eps4[1,2], 
                                                                 eps_lat4[1:]-eps_lat4[1], 
                                                                 sig4[1:,2])))
np.savetxt('triax-sim-res-vert-125-tsaiwu.txt', np.column_stack((eps5[1:,2]-eps5[1,2], 
                                                                 eps_lat5[1:]-eps_lat5[1], 
                                                                 sig5[1:,2])))


plt.show()