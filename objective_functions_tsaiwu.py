'''
Created on Jul 6, 2016

@author: benni
'''

import numpy as np
from YieldFunctions import yfTsaiWuLists
import shelve
import scipy.optimize as scop
import matplotlib.pyplot as plt
import scipy.special as special



def residuals_yf(params, sig_ax, sig_r, beta, weights):
    
    # extract single structural vector and tensor parameters from parameter list
    F2, F3, F12, F23, F44 = params
    
    # compute dependent structural tensor components
    F22 = -F12-F23
    F33 = -2.*F23
    
    # temporarily overwrite of F2 and F3
#     F2 = 3.
#     F3 =10.
    
    # create and initialize structural vector
    strucVec = np.array([F2, F2, F3, 0., 0., 0.])      # transverse isotropy
    
    # create and initialize structural tesnor
    strucTens = np.array([[F22, F12, F23, 0., 0., 0.],
                          [F12, F22, F23, 0., 0., 0.],
                          [F23, F23, F33, 0., 0., 0.],
                          [0., 0., 0., F44, 0., 0.],
                          [0., 0., 0., 0., F44, 0.],
                          [0., 0., 0., 0., 0., 2.*(F22-F12)]])  # transverse isotropy
    
    # create yfTsaiWuLists object
    yf = yfTsaiWuLists([sig_r, sig_r, sig_ax,
                        np.zeros_like(sig_r), np.zeros_like(sig_r), np.zeros_like(sig_r)], 'global',
                       strucVec, 'local',
                       strucTens, 'local', 
                       beta)
    
    # return array with absolute yf values multiplied by the weight factors in the weights array    
    return np.absolute(yf.yf)*weights
    

    
    
def objective_function_yf(params, sig_ax, sig_r, beta, weights):
    
    # compute residuals array
    residuals = residuals_yf(params, sig_ax, sig_r, beta, weights)
    
    # return squared L2-norm of the residuals array (scalar values, which equals the least squares sum)
    return np.linalg.norm(residuals)**2


def objective_function_yf_penalized(params, sig_ax, sig_r, beta, weights, pen_fac):
    
    # compute objective function
    obf = objective_function_yf(params, sig_ax, sig_r, beta, weights)

    # compute penalty term
    penalty = 2*params[0] + params[1]
        
    pen = min(penalty,0.)
    
    obf = obf + pen_fac*0.5*pen**2
        
    return obf



def strains_lat(eps_22, eps_33, r0):
    # computes lateral strains for a 0 deg triax test (vertical bedding planes),
    # based on strains in 22- and 33-direction
    
    # compute initial circumference
    u0 = 2.*np.pi*r0
    
    # compute main radii of ellipsis, based on the strains
    r_22 = r0*(1.+eps_22)
    r_33 = r0*(1.+eps_33)
    
    # check which radius is larger and assign it to a and the smaller one to b
    if r_22 >= r_33:
        a = r_22
        b = r_33
    else:
        a = r_33
        b = r_22
    
    # compute elliptical circumference of deformed section
    aux_e = 1.-b**2/a**2
    u_ellipsis = 4.*a*special.ellipe(aux_e)
    
    # compute lateral strains, based on the elliptical circumference and return it
    return u_ellipsis/u0 - 1.


    
if __name__ == '__main__':
    
    # open binary data file with experimental data needed for identification
    db = shelve.open('triax-data-graz-identification')
    
    # read entries
    for key in db.keys():
        print key + ':'
        print db[key]
    
    # extract needed entries from file
    fcu_inp = db['fcu_meas_mean']
    sig_ax_inp = db['sig_ax']
    
    # close binary file
    db.close()
    
    # create array with axial peak stress values
#     sig_ax = np.array([-fcu_inp[0], -fcu_inp[1], -fcu_inp[2],
#                        sig_ax_inp[0], sig_ax_inp[3], sig_ax_inp[5]])
    sig_ax = np.array([-fcu_inp[0], -fcu_inp[1], -fcu_inp[2],
                       fcu_inp[0]*0.1, fcu_inp[1]*0.1, fcu_inp[2]*0.1,
                       sig_ax_inp[0], sig_ax_inp[3], sig_ax_inp[5]])
    
    # create array with corresponding confing pressures
#     sig_r = np.array([0., 0., 0., -12.5, -12.5, -12.5])
    sig_r = np.array([0., 0., 0.,
                      0., 0., 0.,
                      -12.5, -12.5, -12.5])
    
    # create array with corresponding loading angles in rad
#     beta = np.pi/180.*np.array([0., 60., 90., 0., 60., 90.])
    beta = np.pi/180.*np.array([0., 60., 90., 
                                0., 60., 90.,
                                0., 60., 90.])
    
    # create weights array to gain a better fit for uniaxial stress states
    weights = np.ones_like(sig_ax)
#     weights[0:3] = 4.
    weights[0:3] = 1.
    weights[3:6] = 4.
    
    print residuals_yf([-3., -10., -5., -6., 0.], sig_ax, sig_r, beta, weights), "\n"
    
    # perform leastsquare optimization    
    res = scop.leastsq(residuals_yf, [3., 10., -5., -6., 0.], (sig_ax, sig_r, beta, weights), full_output=False)

    params = res[0]
    print 'identified leastsquares Tsai-Wu parameters: '
    print 'F2: ', params[0]
    print 'F3: ', params[1]
    print 'F12: ', params[2]
    print 'F23: ', params[3]
    print 'F44: ', params[4]
    print 'F22: ', -params[2]-params[3]
    print 'F33: ', -2.*params[3], "\n"

    res_tmp = scop.minimize(objective_function_yf_penalized, [-3., -10., 5., 6., -2.], (sig_ax, sig_r, beta, weights, 0.), 
                         method='L-BFGS-B', bounds=((0.,None), (0.,None), (None, 0.), (None, 0.), (0., None)))
    
    print 'identified penalized and bounded (L-BFGS-B optimization) Tsai-Wu parameters: '
    print 'F2: ', res_tmp.x[0]
    print 'F3: ', res_tmp.x[1]
    print 'F12: ', res_tmp.x[2]
    print 'F23: ', res_tmp.x[3]
    print 'F44: ', res_tmp.x[4]
    print 'F22: ', -res_tmp.x[2]-res_tmp.x[3]
    print 'F33: ', -2.*res_tmp.x[3], "\n"
    
#     for i in range(200, 600, 100):
#          
#         res2 = scop.minimize(objective_function_yf_penalized, res_tmp.x, (sig_ax, sig_r, beta, weights, float(i)), 
#                              method='L-BFGS-B', bounds=((0.,None), (0.,None), (None, 0.), (None, 0.), (0., None)))
#         
#         res_tmp = res2
#         
#         
#     print res2 
     
    # define resolution and ranges for plotting
    n = 200
    range_max = 200.
    range_min = -500.
    range1 = np.linspace(range_max, range_min, n)
     
    # create meshgrid objects
    sig0_mesh, sig33_mesh = np.meshgrid(range1, range1)
     
    # loading angle in rad
    beta = np.pi/180.*0.
     
    # structural vector
    F2 = params[0]
    F3 = params[1]
    strucVec = np.array([F2, F2, F3, 0., 0., 0.])      # transverse isotropy
     
    # structural tensor
    F12 = params[2]
    F23 = params[3]
    F22 = -F12-F23
    F33 = -2.*F23
    F44 = params[4]
    strucTens = np.array([[F22, F12, F23, 0., 0., 0.],
                          [F12, F22, F23, 0., 0., 0.],
                          [F23, F23, F33, 0., 0., 0.],
                          [0., 0., 0., F44, 0., 0.],
                          [0., 0., 0., 0., F44, 0.],
                          [0., 0., 0., 0., 0., 2.*(F22-F12)]])  # transverse isotropy
     
    # create yield function object with identified parameters
    yfIdentified = yfTsaiWuLists([sig0_mesh, sig0_mesh, sig33_mesh, 
                                  np.zeros_like(sig0_mesh), np.zeros_like(sig0_mesh), np.zeros_like(sig0_mesh)], 'global', 
                                 strucVec, 'local',
                                 strucTens, 'local',
                                 beta, internal_csys='global')
     
    yfIdentified.plot_yf_contour(0, 2, scale1= np.sqrt(2.), scaled=True, color='blue')
    
    yfIdentified.update(beta=np.pi/180.*60.)
    yfIdentified.plot_yf_contour(0, 2, scale1= np.sqrt(2.), scaled=True, color='red')

    yfIdentified.update(beta=np.pi/180.*90.)
    yfIdentified.plot_yf_contour(0, 2, scale1= np.sqrt(2.), scaled=True, color='green')
    
    plt.plot(sig_r*np.sqrt(2.), sig_ax, 'o')
    
    plt.show()
