'''
Created on Jan 22, 2014

@author: benni
'''
import numpy as np
# import modleon as md
import tsaiwu as md
from error_functions import UmatError

def call_umat(sigma_old, state_vars_old,
              epsilon_old, delta_epsilon,
              umat_aux_vars):
    
    # copy values of input arrays to working variables
    # that's because the fortran-umat changes its input values
    # which in this case will only result in a change of the working variables,
    # so that the input variables are not changed
    sig_umat = np.array(sigma_old, order='F')
    statev_umat = np.array(state_vars_old, order='F')
    
    # perform stress update with umat
    return_val = md.stress_update.umat_tsaiwu(sig_umat, statev_umat, umat_aux_vars['sse'], umat_aux_vars['spd'], umat_aux_vars['scd'],
                                       epsilon_old, delta_epsilon,
                                       umat_aux_vars['time'], umat_aux_vars['dtime'],
                                       umat_aux_vars['temp'], umat_aux_vars['dtemp'],
                                       umat_aux_vars['predef'], umat_aux_vars['dpred'],
                                       umat_aux_vars['cmname'],
                                       umat_aux_vars['ndi'], umat_aux_vars['nshr'], umat_aux_vars['ntens'],
                                       umat_aux_vars['nstatv'],
                                       umat_aux_vars['props'], umat_aux_vars['nprops'],
                                       umat_aux_vars['coords'], umat_aux_vars['drot'], umat_aux_vars['celent'],
                                       umat_aux_vars['dfgrd0'], umat_aux_vars['dfgrd1'],
                                       umat_aux_vars['noel'], umat_aux_vars['npt'], umat_aux_vars['layer'],
                                       umat_aux_vars['kspt'], umat_aux_vars['kstep'], umat_aux_vars['kinc'])

    # unpack results variables and update temporary variables
    sig_umat[:] = return_val[0]
    statev_umat[:] = return_val[1]
    umat_status = statev_umat[umat_aux_vars['nstatv']-1]
    C_ep_umat = np.array(return_val[2], order='F')

    # check return status of umat
    if umat_status != 0:
        # an error occured in umat
        print 'error in Umat'
        print 'error code: ', umat_status
        raise UmatError
    
    # return working variables
    return (sig_umat, statev_umat, C_ep_umat)