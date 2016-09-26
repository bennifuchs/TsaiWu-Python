'''
Created on Feb 5, 2014

@author: benni
'''
import numpy as np
import single_element_simulator as element_sim
import constants
from error_functions import OuterNewtonError, OuterNewtonConvergenceError, UmatError, OuterNewtonSubsteppingError
import time

def triax_ct_simulation(n_inc, stress_radial, strain_inc_axial,
                        umat_aux_vars, rot_angles):
    # driver function for strain driven triaxial compression test simulation:
    #    1) hydrostatic loading, 
    #    2) axial strains increased and radial stresses kept constant
    #    axial loading is associated with the global x1 direction
    
    # declaration and initialization of stress increment array           
    delta_sig_in = np.zeros((n_inc,6))                  # initialization of input stress increment array, all components set to zero
    
    delta_sig_def = np.ndarray((n_inc,6), dtype=bool)   # initialization of boolean array for definition of defined/undefined stress increment components
    delta_sig_def[:,2] = False                          # delta_s_33 component set to undefined
    delta_sig_def[:,[0,1,3,4,5]] = True                 # remaining stress increment components set to defined

    # declaration and initialization of strain increment array                 
    delta_eps_in = np.zeros_like(delta_sig_in)          # initialization of input strain increment array, all components set to zero
    delta_eps_in[:,2] = strain_inc_axial                # delta_eps_33 component set to axial strain increment
    
    delta_eps_def = np.ndarray((n_inc,6), dtype=bool)   # initialization of boolean array for definition of defined/undefined strain increment components
    delta_eps_def[:,2] = True                           # delta_eps_33 component set to defined
    delta_eps_def[:,[0,1,3,4,5]] = False                # remaining strain increment components set to undefined
      
    # set first load increment to stress driven hydrostatic loading
    delta_sig_in[0,:] = [stress_radial, stress_radial, stress_radial, 0., 0., 0.]
    delta_sig_def[0,:] = [True, True, True, True, True, True]
    delta_eps_in[0,:] = [0., 0., 0., 0., 0., 0.]
    delta_eps_def[0,:] = [False, False, False, False, False, False]
      
    # output arrays
    sig_out = np.zeros((n_inc+1,6))         # output stresses array
    eps_out = np.zeros_like(sig_out)        # output strains array
    statev_out = np.zeros((n_inc+1,umat_aux_vars['nstatv']))        # output state variables array
    n_outer_it_out = np.zeros((n_inc+1,))                 # outer iteration numbers array
    n_inner_it_out = np.zeros((n_inc+1,constants.max_it))        # inner iteration numbers array
      
   
    time_umat = 0                       # initialize umat time with zero
    time_overall1 = time.clock()        # set overall time stamp

    # start loading:
    print 'triax compression test simulation started'
    for i in range(n_inc):               # loop over all increments
        
        try:
            # call mixed strain-stress update routine in the single_element_simulator module,
            # and store return values in output arrays
            #print 'increment: ', i
            sig_out[i+1,:], eps_out[i+1,:], \
            statev_out[i+1,:], \
            n_outer_it_out[i+1], n_inner_it_out[i+1,:], time_umat1 = element_sim.strain_stress_update(sig_out[i,:], delta_sig_in[i,:],
                                                                                                eps_out[i,:], delta_eps_in[i,:],
                                                                                                delta_sig_def[i,:], delta_eps_def[i,:],
                                                                                                statev_out[i,:],
                                                                                                umat_aux_vars, rot_angles)
            time_umat += time_umat1       # update umat time
        except (OuterNewtonConvergenceError, UmatError):
            try:
                            sig_out[i+1,:], eps_out[i+1,:], \
                            statev_out[i+1,:], \
                            n_outer_it_out[i+1], n_inner_it_out[i+1,:], time_umat1 = global_substepper(sig_out[i,:], delta_sig_in[i,:],
                                                                                                eps_out[i,:], delta_eps_in[i,:],
                                                                                                delta_sig_def[i,:], delta_eps_def[i,:],
                                                                                                statev_out[i,:],
                                                                                                umat_aux_vars, rot_angles)
                            time_umat += time_umat1
                            print 'substepping in increment: ', i
            except OuterNewtonSubsteppingError:
                print 'increment with Substepping error: ', i
                raise
        except OuterNewtonError:
            print 'increment with outer Newton error: ', i
            raise
#         except UmatError:
#             print 'increment with Umat error: ', i
#             raise            
    
    
    time_overall = time.clock()-time_overall1       # compute overall time
    
    print 'triax compression test simulation finished'
    print 'overall time spent: ', time_overall
    print 'time spent in outer newton: ', time_overall-time_umat
    print 'time spent in umat: ', time_umat, '\n'
           
    return (sig_out, eps_out, statev_out,
            n_outer_it_out, n_inner_it_out)

    

def general_load_steps_simulation(LoadSteps,
                         umat_aux_vars, rot_angles):
    # driver function for general load steps simulation
    
#     n_steps = 3
#     increments = np.zeros((n_steps,))
#     increments[:] = [200, 400, 600]
# 
#     delta_sig_in = np.zeros((n_steps,6))
#     delta_eps_in = np.zeros_like(delta_sig_in)
#     delta_sig_in[0,:] = [-4.3, 0., 0., 0., 0., 0.]
#     delta_eps_in[1,:] = [-1.0e-5, 0., 0., 0., 0., 0.]
#     delta_eps_in[2,:] = [-1.0e-5, 0., 0., 0., 0., 0.]
#     
#     delta_sig_def = np.ndarray((n_steps,6), dtype=bool)
#     delta_eps_def = np.ndarray((n_steps,6), dtype=bool)
#     delta_sig_def[0,:] = [True for i in range(6)]
#     delta_sig_def[1,:] = [True, False, False, False, False, False]
#     delta_sig_def[2,:] = [True, False, False, False, False, False]
#     delta_eps_def[:,:] = np.logical_not(delta_sig_def)
#     
#     steps = (increments, delta_sig_in, delta_eps_in, delta_sig_def, delta_eps_def)
    
    
    
    local_incs, local_delta_sig_in, local_delta_eps_in, local_delta_sig_def = LoadSteps      # unpack LoadSteps list
    n_steps = len(local_incs)           # number of steps
    n_all_incs = local_incs.sum()       # number of all increments

        
    # input arrays
    all_delta_sig_in = np.zeros((n_all_incs,6))                  # initialization of global input stress increment array, all components set to zero
    all_delta_eps_in = np.zeros_like(all_delta_sig_in)           # initialization of global input strain increment array, all components set to zero
    all_delta_sig_def = np.ndarray((n_all_incs,6), dtype=bool)   # initialization of global boolean array for definition of defined/undefined stress increment components
    all_delta_eps_def = np.ndarray((n_all_incs,6), dtype=bool)   # initialization of global boolean array for definition of defined/undefined strain increment components
    
    i=0             # initial global start increment of each step
    step_no=0       # initial step number
    for step_incs in local_incs:                                            # loop over all steps
        
        all_delta_sig_in[i:i+step_incs,:] = local_delta_sig_in[step_no,:]            # assign local stress increments of each step to global stress increment array
        all_delta_eps_in[i:i+step_incs,:] = local_delta_eps_in[step_no,:]            # assign local strain increments of each step to global strain increment array
        all_delta_sig_def[i:i+step_incs,:] = local_delta_sig_def[step_no,:]       # assign local stress definition of defined/undefined stress increment components of each step to global boolean stress increment array
        #all_delta_eps_def[i:i+step_incs,:] = delta_eps_def[step_no,:]      # assign local stress definition of defined/undefined stress increment components of each step to global boolean stress increment array
        
        step_no += 1
        i += step_incs
        
    all_delta_eps_def[:,:] = np.logical_not(all_delta_sig_def)      # assign local stress definition of defined/undefined stress increment components of each step to global boolean stress increment array (negotian of boolean stress array)

        
    # output arrays
    sig_out = np.zeros((n_all_incs+1,6))                                # output stresses array
    eps_out = np.zeros_like(sig_out)                                    # output strains array
    statev_out = np.zeros((n_all_incs+1,umat_aux_vars['nstatv']))       # output state variables array
    n_outer_it_out = np.zeros((n_all_incs+1,))                          # outer iteration numbers array
    n_inner_it_out = np.zeros((n_all_incs+1,constants.max_it))                      # inner iteration numbers array


    time_umat = 0                       # initialize umat time with zero
    time_overall1 = time.clock()        # set overall time stamp


    # start loading:
    print 'general load step simulation started'
    for i in range(int(n_all_incs)):               # loop over all increments
        
        try:
            # call mixed strain-stress update routine in the single_element_simulator module,
            # and store return values in output arrays
            #print 'increment: ', i
            sig_out[i+1,:], eps_out[i+1,:], \
            statev_out[i+1,:], \
            n_outer_it_out[i+1], n_inner_it_out[i+1,:], time_umat1 = element_sim.strain_stress_update(sig_out[i,:], all_delta_sig_in[i,:],
                                                                                                eps_out[i,:], all_delta_eps_in[i,:],
                                                                                                all_delta_sig_def[i,:], all_delta_eps_def[i,:],
                                                                                                statev_out[i,:],
                                                                                                umat_aux_vars, rot_angles)
            time_umat += time_umat1       # update umat time
        except OuterNewtonConvergenceError:
            try:
                            sig_out[i+1,:], eps_out[i+1,:], \
                            statev_out[i+1,:], \
                            n_outer_it_out[i+1], n_inner_it_out[i+1,:], time_umat1 = global_substepper(sig_out[i,:], all_delta_sig_in[i,:],
                                                                                                eps_out[i,:], all_delta_eps_in[i,:],
                                                                                                all_delta_sig_def[i,:], all_delta_eps_def[i,:],
                                                                                                statev_out[i,:],
                                                                                                umat_aux_vars, rot_angles)
                            time_umat += time_umat1
            except OuterNewtonSubsteppingError:
                print 'increment with Substepping error: ', i
                raise
        except OuterNewtonError:
            print 'increment with outer Newton error: ', i
            raise
        except UmatError:
            print 'increment with Umat error: ', i
            raise
            
    
    time_overall = time.clock()-time_overall1       # compute overall time
    
    print 'general load step simulation finished'
    print 'overall time spent: ', time_overall
    print 'time spent in outer newton: ', time_overall-time_umat
    print 'time spent in umat: ', time_umat, '\n'
           
    return (sig_out, eps_out, statev_out,
            n_outer_it_out, n_inner_it_out)
    

def global_substepper(sigma_old, delta_sigma,
                         epsilon_old, delta_epsilon,
                         delta_sigma_def, delta_epsilon_def,
                         state_vars_old,
                         umat_aux_vars, rot_angles):
        
    delta_eps_sub = np.array(delta_epsilon)                     # strain increments array for substeps
    delta_sig_sub = np.array(delta_sigma)                       # stress increments array for substeps    
    delta_eps_sub[:] = 1./float(constants.n_substeps) * delta_eps_sub     # fill substepping strain increments array with the values of delta_epsilon divided by no of substeps
    delta_sig_sub[:] = 1./float(constants.n_substeps) * delta_sig_sub     # fill substepping stress increments array with the values of delta_epsilon divided by no of substeps
    
    eps_out = np.array(epsilon_old)             # strain output array for substepping, initialized with old strain values
    sig_out = np.array(sigma_old)               # stress output array for substepping, initialized with old stress values
    state_vars_out = np.array(state_vars_old)   # state vars output array, initialized with old state variables
    
    n_inner_it = np.zeros((constants.max_it))   # array for number of inner iterations for every outer iteration per substep
    
    time_umat = 0.                              # set umat time counter to zero
    # loop over all substeps and perform outer newton iteration
    for i in range(constants.n_substeps):
        
        try:
            sig_out[:], eps_out[:], \
            state_vars_out[:], \
            n_outer_it, n_inner_it[:], time_umat1 = element_sim.strain_stress_update(sig_out, delta_sig_sub,
                                                        eps_out, delta_eps_sub,
                                                        delta_sigma_def, delta_epsilon_def,
                                                        state_vars_out,
                                                        umat_aux_vars, rot_angles)
            time_umat += time_umat1     # update time counter
            
        except OuterNewtonConvergenceError:
            print 'global substepping increment with outer Newton convergence error: ', i
            raise OuterNewtonSubsteppingError     
        except OuterNewtonError:
            print 'global substepping increment with outer Newton error: ', i
            raise OuterNewtonSubsteppingError
        except UmatError:
            print 'global substepping increment with Umat error: ', i
            raise OuterNewtonSubsteppingError
            
        
    return (sig_out, eps_out, state_vars_out, n_outer_it, n_inner_it, time_umat)
    
    