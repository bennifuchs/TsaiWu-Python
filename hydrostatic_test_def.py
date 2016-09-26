'''
Created on Jul 25, 2016

@author: benni
'''
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

import os

from IPsimulator import constants
from IPsimulator import umat_setup
from IPsimulator import loading


#----------------------------------------------------------------------------------
# DEFINITION AND INITIALIZATION OF VARIABLES
#----------------------------------------------------------------------------------
 
umat_aux_vars = umat_setup.set_variables(mat_file_name='test-material.inp')    # load material data and initialize auxiliary variables for UMAT 

print 'number of state variables: ', umat_aux_vars['nstatv']
 
output_dict = {}                        # empty dictionary for matlab output file
np.set_printoptions(linewidth=100)      # reset linewidth of numpy print command
plt.figure(1)                           # create first figure for plots
 
#---------------------------------------------------------------
# general_load_steps_simulation (strain driven hydrostatic compressive loading)
#---------------------------------------------------------------

# angle beta between horizontal plane and bedding plane of rock
beta = 0.
 
rot_angles = [0.0, beta, 0.0]        # definition of rotation angles list with alpha, beta, gamma of material orientation
 
# step and number of increments definition
n_steps = 1                                 # number of steps
increments = np.zeros((n_steps,))           # initialize increments per step array
increments[:] = [50]                       # define increments per step
  
# active/inactive load components (delta epsilon/sigma) definition
delta_sig_def = np.ndarray((n_steps,6), dtype=bool)     # boolean array for definition of active/inactive stress increment components in each step
delta_sig_def[:,0:3] = False                            # all steps: delta_sig_11, delta_sig_22, delta_sig_33 inactive (delta_eps11, delta_eps_22, delta_eps_33 active)
delta_sig_def[:,3:] = True                              # all steps: delta_sig_12 to delta_sig_23 active (shear stresses)
  
  
# load components definition
delta_sig_in = np.zeros((n_steps,6))                    # zeros array for definition of stress increments in each step, all active stress components set to zero (shear stresses)
delta_eps_in = np.zeros_like(delta_sig_in)              # zeros array for definition of strain increments in each step
  
# set delta_eps components
delta_eps_in[:,0:3] = -8.0e-5               # assign delta_eps_11, delta_eps_22, delta_eps33 values for all steps

# create list with all load definition arrays  
steps = (increments, delta_sig_in, delta_eps_in, delta_sig_def)

# define and initialize arrays for return data of general load step simulation
n_incs_all = sum(increments)+1                              # number of all increments
sig1 = np.empty((n_incs_all,6))                             # stress output array
eps1 = np.empty_like(sig1)                                  # strain output array
statev1 = np.empty((n_incs_all,umat_aux_vars['nstatv']))    # state vars output array
outer_it1 = np.empty((n_incs_all,))                         # outer iteration numbers output array
inner_it1 = np.empty((n_incs_all,constants.max_it))         # inner iteration numbers output array

# call general load step simulation
sig1[:,:], eps1[:,:], statev1[:,:], outer_it1[:], inner_it1[:,:] = loading.general_load_steps_simulation(steps, umat_aux_vars, rot_angles)
 
output_dict.update({'sig_hydrostatic_compression': sig1, 'eps_hydrostatic_compression': eps1})  # append sigma and epsilon arrays to output dictionary
plt.plot(eps1[:,0], sig1[:,0], color='b', marker='', label='hydrostatic compressive loading' )                        # plot strain stress curve

print 'number of outer iterations in last load increment: '  
print outer_it1[increments.sum()]
print 'number of inner iterations per outer iterations in last load increment: '
print inner_it1[increments.sum()]
print 'number of active yield functions in each load increment: '
print statev1[:,0]
print 'internal hardening variable after each load increment: '
print statev1[:,15]


#---------------------------------------------------------------
# plot options for strain stress curve in figure(2)
#---------------------------------------------------------------
title = 'triaxial compression test (strain driven)'
plt.ticklabel_format(style='sci', axis='x', scilimits=(-2,2))
plt.title(title)
plt.xlabel(r'$\varepsilon_{11}$ [-]')
plt.ylabel('$\sigma_{11}$ [N/mm$^2$]')
plt.legend(loc=2)
plt.grid(True)

# os.chdir('./output-files')
# sio.savemat('results-hydrostatic-loading-straindriven', output_dict, oned_as='column')
# os.chdir('./..')

plt.show()