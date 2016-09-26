'''
Created on Dec 13, 2013

@author: benni
'''
import numpy as np
import constants
# import modleon as md
import tsaiwu as md
from error_functions import OuterNewtonError, OuterNewtonConvergenceError
import umat_driver
import time

# mixed stress/strain driven
def strain_stress_update(sigma_old, delta_sigma,
                         epsilon_old, delta_epsilon,
                         delta_sigma_def, delta_epsilon_def,
                         state_vars_old,
                         umat_aux_vars, rot_angles):
    # performs outer newton iteration
    # i.e. computes strain vector at the end of a load/stress increment
    #
    # input values:
    # ------------
    # sigma_old      stress vector at the beginning of the load increment
    # delta_sigma    stress increment
    # epsilon_old    strain vector at the beginning of the load increment
    # delta_epsilon  strain increment
    # delta_sigma_def    boolean list with active components of stress increment vector (active/defined means True, inactive/undefined means False)
    # delta_epsilon_def  boolean list with active components of strain increment vector (active/defined means True, inactive/undefined means False)
    # state_vars_old     state variables array at the beginning of the step
    # umat_aux_vars  dictionary with auxiliary unimportant variables needed by umat
    # rot_angles     tuple with rotation angles which define the directions of the material axes of the rock within the global coordinate system
    #                    alpha: rotation angle, around z-axis, eastwards positive (zyz convention)
    #                    beta: rotation angle, around rotated y-axis, downwards positive (zyz convention)
    #                    gamma: third rotation angle, around rotated z-axis, eastwards positive (zyz convention)
    #
    # return values:
    # --------------
    # epsilon        strain vector at the end of the load increment
    # state_vars     state variables array at the end of the load increment
    
    # -----------------------------------------------------
    # DATA CHECK
    # -----------------------------------------------------    
    # list for number of inner iterations performed within each outer iteration
    no_inner_it = [-1 for i in range(constants.max_it)]
    
    # unpacking of rotation angles
    alpha, beta, gamma = rot_angles
    
    # rotate old stress and old strain vector into material direction
    sigma_old_rot = md.stress_routines.rotate_sigma(sigma_old, alpha, beta, gamma)          # rotated old stress vector of last step
    epsilon_old_rot = md.stress_routines.rotate_epsilon(epsilon_old, alpha, beta, gamma)    # rotated old strain vector of last step

    # check delta_sigma_def and delta_epsilon_def
    for i in range(6):        
        
        # check if ith entry of delta_sigma_def equals ith entry of delta_epsilon_def
        if delta_sigma_def[i] == delta_epsilon_def[i]:
            
            # entries are the same, raise exception
            raise OuterNewtonError
    

    # -----------------------------------------------------
    # DATA PREPARATION
    # -----------------------------------------------------
    # number of unknown stress components
    n_sigma_undef = sum(delta_sigma_def == False )
    #print "number of unknown stress components: ", n_sigma_undef
    # number of unknown strain components
    #n_epsilon_undef = 6-n_sigma_undef
    #print "number of unknown strain components: ", n_epsilon_undef, "\n"
    
    # create list with indices of known delta sigma components
    # and unknown delta sigma components
    ind_dsig_def = [i for i in range(6) if delta_sigma_def[i]==True]
    ind_dsig_undef = [i for i in range(6) if delta_sigma_def[i]==False]
    #print "indices of defined stress components: ", ind_dsig_def
    #print "indices of undefined stress components: ", ind_dsig_undef, "\n"
    
    # create list with indices of known delta epsilon components
    # and unknown delta epsilon components
    ind_deps_def = ind_dsig_undef
    ind_deps_undef = ind_dsig_def
    #print "indices of defined strain components: ", ind_deps_def
    #print "indices of undefined strain components: ", ind_deps_undef, "\n"



    # -----------------------------------------------------
    # INITIAL COMPUTATIONS
    # -----------------------------------------------------    
    # define vector of unknown quantities and initialize it with converged values from the last step
    y = np.empty((6,))
    y[ind_dsig_undef] = sigma_old[ind_dsig_undef]
    y[ind_deps_undef] = epsilon_old[ind_deps_undef]
    
    # define vector for increment of unknown quantities (are computed in newton iteration)
    delta_y = np.zeros_like(y)
    
    # define first approximation of updated strain vector
    # undefined components are set to converged values of last step
    # defined components are set to defined value
    eps_tmp = np.array(epsilon_old)
    eps_tmp[ind_deps_def] = eps_tmp[ind_deps_def] + delta_epsilon[ind_deps_def]
    #print "first approximation of updated strain vector: ", eps_tmp, "\n"

    # define first approximation of updated stess vector
    # undefined components are set to converged values of last step
    # defined components are set to defined value
    sigma_tmp = np.array(sigma_old)
    sigma_tmp[ind_dsig_def] = sigma_tmp[ind_dsig_def] + delta_sigma[ind_dsig_def]
    #print "first approximation of updated stress vector: ", sigma_tmp, "\n"
    
    # initialize temporary elastoplastic tangent as empty array
    C_ep_tmp = np.empty((6,6), order='F')
    
    # initialize jacobian as empty array
    jac = np.empty_like(C_ep_tmp)

# -----------------------------------------------------
# BEGIN NEWTON ITERATION
# -----------------------------------------------------
    
    # define temporary variables
    sig_umat =np.empty_like(sigma_old)          # temporary stress vector from umat
    statev_tmp = np.empty_like(state_vars_old)  # temporary state variables from umat
    delta_eps_tmp = np.empty_like(epsilon_old)  # temporary strain increment vector as input to umat

    sig_umat_rot = np.empty_like(sig_umat)                  # rotated temporary stress vector from umat
    C_ep_tmp_rot = np.empty_like(C_ep_tmp)                  # rotated temporary elastoplastic tangent
    delta_eps_tmp_rot = np.empty_like(delta_eps_tmp)        # temporary rotated strain increment vector as input to umat
    
    # define residuum vector R
    R = np.empty((6,))
    
    k = 1       # counter for number of iterations
    t_all = 0   # set time counter to zero
    while True:
        
        if k > constants.max_it:
            # solution is not converging
            raise OuterNewtonConvergenceError            
        
        # compute temporary strain increment vector
        delta_eps_tmp[:] = eps_tmp - epsilon_old
        #print '\t delta_eps_tmp: ', delta_eps_tmp

        # ------------------------------------------------------------------        
        # rotate input variables for umat in direction of rock orientation
        delta_eps_tmp_rot[:] = md.stress_routines.rotate_epsilon(delta_eps_tmp, alpha, beta, gamma) 
        t1 = time.clock()       # set first time stamp   
        # perform return mapping step by calling the umat through the special calling function
        # to obtain updated stress vector, state variables vector and elastic-plastic tangent matrix
        sig_umat_rot[:], statev_tmp[:], C_ep_tmp_rot[:,:] = umat_driver.call_umat(sigma_old_rot, state_vars_old, 
                                                                                       epsilon_old_rot, delta_eps_tmp_rot, 
                                                                                       umat_aux_vars)        
        t2 = time.clock()       # set second time stamp
        t_all += (t2-t1)        # update time counter
        no_inner_it[k-1] = statev_tmp[umat_aux_vars['nstatv']-2] # store number of inner iterations
        
#         print 'C_ep_tmp out of Umat:'
#         print C_ep_tmp_rot
        
        # rotate output variables from umat back into global coordinate system
        sig_umat[:] = md.stress_routines.rotate_sigma(sig_umat_rot, -alpha, -beta, -gamma)
        C_ep_tmp[:,:] = md.stress_routines.rotate_material_matrix(C_ep_tmp_rot, -alpha, -beta, -gamma)
#         print 'C_ep_tmp rotated back:'
#         print C_ep_tmp
        # ------------------------------------------------------------------        
        
        # debug information
        #print '\t n_inner_it: ', statev_tmp[umat_aux_vars['nstatv']-2]
        #print '\t alpha: ', statev_tmp[15]
        #print '\t alpha_d: ', statev_tmp[25]
        #print '\t omega: ', statev_tmp[32]
        #print '\t C_ep_tmp: '
        #print C_ep_tmp
        
        # -----------------------------------------------------
        # COMPUTATION OF RESIDUUM VECTOR
        # -----------------------------------------------------    
        R[:] = sig_umat - sigma_tmp
        
        # debug information
        #print '\t R: ', R        
        
        # compute norm of residuum vector
        R_norm = np.linalg.norm(R)
        
        # debug information
#         print '\t R_norm: ', R_norm
        
        # check if convergence criterion is met
        if R_norm < constants.tolerance:           
            
            # norm is small than tolerance,
            # exit loop
            break
        
        # -----------------------------------------------------
        # COMPUTATION OF DERIVATIVES OF RESIDUUM VECTOR (dR/dy - jacobian)
        # -----------------------------------------------------
        # set jacobian to zero
        jac[:,:] = 0.
        
        # assignment of dR_i/dsig_undef_j = -1.0
        # to the right positions of the jacobian
        jac[ind_dsig_undef, ind_dsig_undef] = -1.0
        
        # assignment of dR_i/deps_undef_j = Cep_ij
        # to the right positions of the jacobian
        jac[:,ind_deps_undef] = C_ep_tmp[:,ind_deps_undef]
            
#         print 'jac:'
#         print jac
        
        # -----------------------------------------------------
        # SOLVE SYSTEM OF LINEAR EQUATIONS
        # -----------------------------------------------------
        delta_y[:] = np.linalg.solve(jac, -R)
        
        # update y
        y[:] = y + delta_y
 
        
        # -----------------------------------------------------
        # UPDATE STRESS AND STRAIN VECTOR WITH VALUES FROM y
        # AND COMPUTE NEW RESIDUUM
        # -----------------------------------------------------        
        # assign stress values from y to sigma_tmp
        sigma_tmp[ind_dsig_undef] = y[ind_dsig_undef]
        #print '\t sigma_tmp: ', sigma_tmp
        
        # assign strain values from y to eps_tmp
        eps_tmp[ind_deps_undef] = y[ind_deps_undef]
        #print '\t eps_tmp: ', eps_tmp
        
        # increase counter
        k+=1

    
#     # derive number of inner newton iterations
#     no_inner_it = statev_tmp[umat_aux_vars['nstatv']-2]
    
    # return epsilon, sigma, C_ep_tmp, state variables, number of outer newton iterations
    # and array with number of inner newton iterations per outer iteration
    return (sigma_tmp, eps_tmp, statev_tmp, k, no_inner_it, t_all)


# # mixed stress/strain driven
# def strain_stress_update_anisotropic(sigma_old, delta_sigma,
#                                      epsilon_old, delta_epsilon,
#                                      delta_sigma_def, delta_epsilon_def,
#                                      state_vars_old,
#                                      umat_aux_vars, rot_angles):
#     # performs outer newton iteration for anisotropic elastic-plastic material
#     # by rotating boundary conditions towards material orientation
#     # and performing outer newton iteration in material orientation coordinate system
#     # residuum vector R is defined in material orientation coordinate system
#     # i.e. computes strain vector at the end of a load/stress increment
#     #
#     # input values:
#     # ------------
#     # sigma_old      stress vector at the beginning of the load increment
#     # delta_sigma    stress increment
#     # epsilon_old    strain vector at the beginning of the load increment
#     # delta_epsilon  strain increment
#     # delta_sigma_def    boolean list with active components of stress increment vector (active/defined means True, inactive/undefined means False)
#     # delta_epsilon_def  boolean list with active components of strain increment vector (active/defined means True, inactive/undefined means False)
#     # state_vars_old     state variables array at the beginning of the step
#     # umat_aux_vars  dictionary with auxiliary unimportant variables needed by umat
#     # rot_angles     tuple with rotation angles which define the directions of the material axes of the rock within the global coordinate system
#     #                    alpha: rotation angle, around z-axis, eastwards positive (zyz convention)
#     #                    beta: rotation angle, around rotated y-axis, downwards positive (zyz convention)
#     #                    gamma: third rotation angle, around rotated z-axis, eastwards positive (zyz convention)
#     #
#     # return values:
#     # --------------
#     # epsilon        strain vector at the end of the load increment
#     # state_vars     state variables array at the end of the load increment
#     
#     # -----------------------------------------------------
#     # DATA CHECK
#     # -----------------------------------------------------    
#     # list for number of inner iterations performed within each outer iteration
#     no_inner_it = [-1 for i in range(constants.max_it)]
#     
#     # unpacking of rotation angles
#     alpha, beta, gamma = rot_angles
#     
#     # compute rotation matrix for stress and strain tensor
#     T_eps = md.stress_routines.calc_rot_matrix_epsilon(alpha, beta, gamma)      # rotation matrix for strains
#     T_sig = md.stress_routines.calc_rot_matrix_sigma(alpha, beta, gamma)        # rotation matrix for stresses
#     
#     # rotate old stress and old strain vector into material direction
#     sigma_old_rot = md.stress_routines.rotate_sigma(sigma_old, alpha, beta, gamma)          # rotated old stress vector of last step
#     epsilon_old_rot = md.stress_routines.rotate_epsilon(epsilon_old, alpha, beta, gamma)    # rotated old strain vector of last step
# #     print 'sigma_old_rot: '
# #     print sigma_old_rot
# #     print 'epsilon_old_rot: '
# #     print epsilon_old_rot
# 
#     # check delta_sigma_def and delta_epsilon_def
#     for i in range(6):        
#         
#         # check if ith entry of delta_sigma_def equals ith entry of delta_epsilon_def
#         if delta_sigma_def[i] == delta_epsilon_def[i]:
#             
#             # entries are the same, raise exception
#             raise OuterNewtonError
#     
# 
#     # -----------------------------------------------------
#     # DATA PREPARATION
#     # -----------------------------------------------------
#     # number of unknown stress components
#     n_sigma_undef = sum(delta_sigma_def == False )
#     #print "number of unknown stress components: ", n_sigma_undef
#     # number of unknown strain components
#     #n_epsilon_undef = 6-n_sigma_undef
#     #print "number of unknown strain components: ", n_epsilon_undef, "\n"
#     
#     # create list with indices of known delta sigma components
#     # and unknown delta sigma components
#     ind_dsig_def = [i for i in range(6) if delta_sigma_def[i]==True]
#     ind_dsig_undef = [i for i in range(6) if delta_sigma_def[i]==False]
#     #print "indices of defined stress components: ", ind_dsig_def
#     #print "indices of undefined stress components: ", ind_dsig_undef, "\n"
#     
#     # create list with indices of known delta epsilon components
#     # and unknown delta epsilon components
#     ind_deps_def = ind_dsig_undef
#     ind_deps_undef = ind_dsig_def
#     #print "indices of defined strain components: ", ind_deps_def
#     #print "indices of undefined strain components: ", ind_deps_undef, "\n"
# 
# 
# 
#     # -----------------------------------------------------
#     # INITIAL COMPUTATIONS
#     # -----------------------------------------------------    
#     # define vector of unknown quantities and initialize it with converged values from the last step
#     y = np.empty((6,))
#     y[ind_dsig_undef] = sigma_old[ind_dsig_undef]
#     y[ind_deps_undef] = epsilon_old[ind_deps_undef]
#     
#     # define vector for increment of unknown quantities (are computed in newton iteration)
#     delta_y = np.zeros_like(y)
#     
#     # define first approximation of updated strain vector
#     # undefined components are set to converged values of last step
#     # defined components are set to defined value
#     eps_tmp = np.array(epsilon_old)
#     eps_tmp[ind_deps_def] = eps_tmp[ind_deps_def] + delta_epsilon[ind_deps_def]
#     #print "first approximation of updated strain vector: ", eps_tmp, "\n"
# 
#     # define first approximation of updated stess vector
#     # undefined components are set to converged values of last step
#     # defined components are set to defined value
#     sigma_tmp = np.array(sigma_old)
#     sigma_tmp[ind_dsig_def] = sigma_tmp[ind_dsig_def] + delta_sigma[ind_dsig_def]
#     #print "first approximation of updated stress vector: ", sigma_tmp, "\n"
#     
#     # initialize temporary rotated elastoplastic tangent as empty array
#     C_ep_tmp_rot = np.empty((6,6), order='F')
#     
#     # initialize jacobian as empty array
#     jac = np.empty_like(C_ep_tmp_rot)
# 
# # -----------------------------------------------------
# # BEGIN NEWTON ITERATION
# # -----------------------------------------------------
#     
#     # define temporary variables
#     sig_umat_rot =np.empty_like(sigma_old)      # temporary rotated stress vector from umat
#     statev_tmp = np.empty_like(state_vars_old)  # temporary state variables from umat
#     delta_eps_tmp = np.empty_like(epsilon_old)  # temporary strain increment vector
#     delta_eps_tmp_rot = np.empty_like(delta_eps_tmp)        # temporary rotated strain increment vector as input to umat
# 
#     sigma_tmp_rot = np.empty_like(sigma_tmp)                # rotated approximation of updated stress vector
#     T_eps_u = np.zeros_like(T_eps)                          # rotation matrix for unknown strain components
#     T_sig_u = np.zeros_like(T_sig)                          # rotation matrix for unknown stress components    
#     T_eps_u[:,ind_deps_undef] = T_eps[:,ind_deps_undef]     # fill columns with indices of unknown strain components with corresponding columns of strain transformation matrix
#     T_sig_u[:,ind_dsig_undef] = T_sig[:,ind_dsig_undef]     # fill columns with indices of unknown stress components with corresponding columns of stress transformation matrix
#     
#     # define residuum vector R
#     R = np.empty((6,))
#     
#     k = 1       # counter for number of iterations
#     t_all = 0   # set time counter to zero
#     while True:
#         
#         if k > constants.max_it:
#             # solution is not converging
#             raise OuterNewtonConvergenceError            
#         
#         # compute temporary strain increment vector
#         delta_eps_tmp[:] = eps_tmp - epsilon_old
# #         print '\t delta_eps_tmp: ', delta_eps_tmp
# 
#         # ------------------------------------------------------------------        
#         # rotate input variables for umat in direction of rock orientation
#         delta_eps_tmp_rot[:] = md.stress_routines.rotate_epsilon(delta_eps_tmp, alpha, beta, gamma) 
# #         print '\tdelta_eps_tmp_rot: ', delta_eps_tmp_rot
#         t1 = time.clock()       # set first time stamp   
#         # perform return mapping step by calling the umat through the special calling function
#         # to obtain updated stress vector, state variables vector and elastic-plastic tangent matrix
#         sig_umat_rot[:], statev_tmp[:], C_ep_tmp_rot[:,:] = umat_driver.call_umat(sigma_old_rot, state_vars_old, 
#                                                                                        epsilon_old_rot, delta_eps_tmp_rot, 
#                                                                                        umat_aux_vars)        
#         t2 = time.clock()       # set second time stamp
#         t_all += (t2-t1)        # update time counter
#         no_inner_it[k-1] = statev_tmp[umat_aux_vars['nstatv']-2] # store number of inner iterations
# #         print 'yield flags: ', statev_tmp[16:18]
# #         print 'no inner-it: ', no_inner_it[k-1]
# #         print 'alpha_p: ', statev_tmp[15]
#         # ------------------------------------------------------------------        
#         
#         # debug information
# #         print '\t n_inner_it: ', statev_tmp[umat_aux_vars['nstatv']-2]
# #         print '\t alpha: ', statev_tmp[15]
# #         print '\t alpha_d: ', statev_tmp[25]
# #         print '\t delta_eps_p_vol: ', sum(statev_tmp[26:29])
# #         print '\t omega: ', statev_tmp[32]
# #         print '\tsig_umat_rot: ', sig_umat_rot
# #         print '\tsig_umat rotated back', md.stress_routines.rotate_sigma(sig_umat_rot, -alpha, -beta, -gamma)
# #         print '\t C_ep_tmp_rot: '
# #         print C_ep_tmp_rot
#         
#         # -----------------------------------------------------
#         # COMPUTATION OF RESIDUUM VECTOR
#         # -----------------------------------------------------    
#         # rotate sigma_tmp
#         sigma_tmp_rot[:] = md.stress_routines.rotate_sigma(sigma_tmp, alpha, beta, gamma)
# #         print '\tsigma_tmp_rot: ', sigma_tmp_rot
#         
#         # compute residuum vector
#         R[:] = sig_umat_rot - sigma_tmp_rot
#         
#         # debug information
# #         print '\t R: ', R        
#         
#         # compute norm of residuum vector
#         R_norm = np.linalg.norm(R)
#         
#         # debug information
# #         print '\t R_norm: ', R_norm, '\n'
#         
#         # check if convergence criterion is met
#         if R_norm < constants.tolerance:           
#             
#             # norm is small than tolerance,
#             # exit loop
#             break
#         
#         # -----------------------------------------------------
#         # COMPUTATION OF DERIVATIVES OF RESIDUUM VECTOR (dR/dy - jacobian)
#         # -----------------------------------------------------
#         # set jacobian to zero
#         jac[:,:] = 0.
#         
#         # computation of jacobian
#         jac[:,:] = np.dot(C_ep_tmp_rot, T_eps_u) - T_sig_u
#         
# #         print 'jac:'
# #         print jac
#         
#         # -----------------------------------------------------
#         # SOLVE SYSTEM OF LINEAR EQUATIONS
#         # -----------------------------------------------------
#         delta_y[:] = np.linalg.solve(jac, -R)
#         
#         # update y
#         y[:] = y + delta_y
#         
#         # debug information
# #         print '\tdelta_y: ', delta_y
# #         print '\ty: ', y
#         
#         # -----------------------------------------------------
#         # UPDATE STRESS AND STRAIN VECTOR WITH VALUES FROM y
#         # AND COMPUTE NEW RESIDUUM
#         # -----------------------------------------------------        
#         # assign stress values from y to sigma_tmp
#         sigma_tmp[ind_dsig_undef] = y[ind_dsig_undef]
# #         print '\t sigma_tmp: ', sigma_tmp
#         
#         # assign strain values from y to eps_tmp
#         eps_tmp[ind_deps_undef] = y[ind_deps_undef]
# #         print '\t eps_tmp: ', eps_tmp
#         
#         # increase counter
#         k+=1
# 
#     
# #    # derive number of inner newton iterations
# #     no_inner_it = statev_tmp[umat_aux_vars['nstatv']-2]
#         
#     # return epsilon, sigma, C_ep_tmp, state variables, number of outer newton iterations
#     # and array with number of inner newton iterations per outer iteration
#     return (sigma_tmp, eps_tmp, statev_tmp, k, no_inner_it, t_all)



# # mixed stress/strain driven
# def strain_stress_update_anisotropic2(sigma_old, delta_sigma,
#                                          epsilon_old, delta_epsilon,
#                                          delta_sigma_def, delta_epsilon_def,
#                                          state_vars_old,
#                                          umat_aux_vars, rot_angles):
#     # performs outer newton iteration for anisotropic elastic-plastic material
#     # by rotating stress vector and elastoplastic tangent from material orientation coordinate system towards global coordinate system
#     # and performing outer newton iteration in global coordinate system
#     # residuum vector R is defined in global coordinate system
#     # i.e. computes strain vector at the end of a load/stress increment
#     #
#     # input values:
#     # ------------
#     # sigma_old      stress vector at the beginning of the load increment
#     # delta_sigma    stress increment
#     # epsilon_old    strain vector at the beginning of the load increment
#     # delta_epsilon  strain increment
#     # delta_sigma_def    boolean list with active components of stress increment vector (active/defined means True, inactive/undefined means False)
#     # delta_epsilon_def  boolean list with active components of strain increment vector (active/defined means True, inactive/undefined means False)
#     # state_vars_old     state variables array at the beginning of the step
#     # umat_aux_vars  dictionary with auxiliary unimportant variables needed by umat
#     # rot_angles     tuple with rotation angles which define the directions of the material axes of the rock within the global coordinate system
#     #                    alpha: rotation angle, around z-axis, eastwards positive (zyz convention)
#     #                    beta: rotation angle, around rotated y-axis, downwards positive (zyz convention)
#     #                    gamma: third rotation angle, around rotated z-axis, eastwards positive (zyz convention)
#     #
#     # return values:
#     # --------------
#     # epsilon        strain vector at the end of the load increment
#     # state_vars     state variables array at the end of the load increment
#     
#     # -----------------------------------------------------
#     # DATA CHECK
#     # -----------------------------------------------------
# 
#     # list for number of inner iterations performed within each outer iteration
#     no_inner_it = [-1 for i in range(constants.max_it)]
#     
#     # unpacking of rotation angles
#     alpha, beta, gamma = rot_angles
#         
#     # rotate old stress and old strain vector into material direction
#     sigma_old_rot = md.stress_routines.rotate_sigma(sigma_old, alpha, beta, gamma)          # rotated old stress vector of last step
#     epsilon_old_rot = md.stress_routines.rotate_epsilon(epsilon_old, alpha, beta, gamma)    # rotated old strain vector of last step
# 
#     # check delta_sigma_def and delta_epsilon_def
#     for i in range(6):        
#         
#         # check if ith entry of delta_sigma_def equals ith entry of delta_epsilon_def
#         if delta_sigma_def[i] == delta_epsilon_def[i]:
#             
#             # entries are the same, raise exception
#             raise OuterNewtonError
#     
# 
#     # -----------------------------------------------------
#     # DATA PREPARATION
#     # -----------------------------------------------------
#     # number of unknown stress components
#     n_sigma_undef = sum(delta_sigma_def == False )
#     #print "number of unknown stress components: ", n_sigma_undef
#     # number of unknown strain components
#     #n_epsilon_undef = 6-n_sigma_undef
#     #print "number of unknown strain components: ", n_epsilon_undef, "\n"
#     
#     # create list with indices of known delta sigma components
#     # and unknown delta sigma components
#     ind_dsig_def = [i for i in range(6) if delta_sigma_def[i]==True]
#     ind_dsig_undef = [i for i in range(6) if delta_sigma_def[i]==False]
#     #print "indices of defined stress components: ", ind_dsig_def
#     #print "indices of undefined stress components: ", ind_dsig_undef, "\n"
#     
#     # create list with indices of known delta epsilon components
#     # and unknown delta epsilon components
#     ind_deps_def = ind_dsig_undef
#     ind_deps_undef = ind_dsig_def
#     #print "indices of defined strain components: ", ind_deps_def
#     #print "indices of undefined strain components: ", ind_deps_undef, "\n"
# 
# 
# 
#     # -----------------------------------------------------
#     # INITIAL COMPUTATIONS
#     # -----------------------------------------------------    
#     # define vector of unknown quantities and initialize it with converged values from the last step
#     y = np.empty((6,))
#     y[ind_dsig_undef] = sigma_old[ind_dsig_undef]
#     y[ind_deps_undef] = epsilon_old[ind_deps_undef]
#     
#     # define vector for increment of unknown quantities (are computed in newton iteration)
#     delta_y = np.zeros_like(y)
#     
#     # define first approximation of updated strain vector
#     # undefined components are set to converged values of last step
#     # defined components are set to defined value
#     eps_tmp = np.array(epsilon_old)
#     eps_tmp[ind_deps_def] = eps_tmp[ind_deps_def] + delta_epsilon[ind_deps_def]
#     #print "first approximation of updated strain vector: ", eps_tmp, "\n"
# 
#     # define first approximation of updated stess vector
#     # undefined components are set to converged values of last step
#     # defined components are set to defined value
#     sigma_tmp = np.array(sigma_old)
#     sigma_tmp[ind_dsig_def] = sigma_tmp[ind_dsig_def] + delta_sigma[ind_dsig_def]
#     #print "first approximation of updated stress vector: ", sigma_tmp, "\n"
#     
#     # initialize temporary elastoplastic tangent as empty array
#     C_ep_tmp = np.empty((6,6), order='F')
#     
#     # initialize jacobian as empty array
#     jac = np.empty_like(C_ep_tmp)
# 
# # -----------------------------------------------------
# # BEGIN NEWTON ITERATION
# # -----------------------------------------------------
#     
#     # define temporary variables
#     sig_umat =np.empty_like(sigma_old)          # temporary stress vector from umat
#     statev_tmp = np.empty_like(state_vars_old)  # temporary state variables from umat
#     delta_eps_tmp = np.empty_like(epsilon_old)  # temporary strain increment vector as input to umat
# 
#     sig_umat_rot = np.empty_like(sig_umat)                  # rotated temporary stress vector from umat
#     C_ep_tmp_rot = np.empty_like(C_ep_tmp)                  # rotated temporary elastoplastic tangent
#     delta_eps_tmp_rot = np.empty_like(delta_eps_tmp)        # temporary rotated strain increment vector as input to umat
#     
#     # define residuum vector R
#     R = np.empty((6,))
#     
#     k = 1       # counter for number of iterations
#     t_all = 0   # set time counter to zero
#     while True:
#         
#         if k > constants.max_it:
#             # solution is not converging
#             raise OuterNewtonConvergenceError            
#         
#         # compute temporary strain increment vector
#         delta_eps_tmp[:] = eps_tmp - epsilon_old
# #         print '\t delta_eps_tmp: ', delta_eps_tmp
# 
#         # ------------------------------------------------------------------        
#         # rotate input variables for umat in direction of rock orientation
#         delta_eps_tmp_rot[:] = md.stress_routines.rotate_epsilon(delta_eps_tmp, alpha, beta, gamma) 
# #         print '\tdelta_eps_tmp_rot: ', delta_eps_tmp_rot
#         t1 = time.clock()       # set first time stamp   
#         # perform return mapping step by calling the umat through the special calling function
#         # to obtain updated stress vector, state variables vector and elastic-plastic tangent matrix
#         sig_umat_rot[:], statev_tmp[:], C_ep_tmp_rot[:,:] = umat_driver.call_umat(sigma_old_rot, state_vars_old, 
#                                                                                        epsilon_old_rot, delta_eps_tmp_rot, 
#                                                                                        umat_aux_vars)        
#         t2 = time.clock()       # set second time stamp
#         t_all += (t2-t1)        # update time counter
#         no_inner_it[k-1] = statev_tmp[umat_aux_vars['nstatv']-2] # store number of inner iterations
#         # ------------------------------------------------------------------        
#         
#         # debug information
# #         print '\t n_inner_it: ', statev_tmp[umat_aux_vars['nstatv']-2]
# #         print '\t alpha: ', statev_tmp[15]
# #         print '\t alpha_d: ', statev_tmp[25]
# #         print '\t omega: ', statev_tmp[32]
# #         print '\tsig_umat_rot: ', sig_umat_rot
# #         print '\tsig_umat rotated back', md.stress_routines.rotate_sigma(sig_umat_rot, -alpha, -beta, -gamma)
# #         print '\t C_ep_tmp_rot: '
# #         print C_ep_tmp_rot
# 
# #         eps_test = np.linalg.solve(C_ep_tmp_rot, sig_umat_rot)
# #         print '\teps_test: '
# #         print eps_test
# #         print '\tC_ep_tmp_rot_inv: '
# #         print np.linalg.inv(C_ep_tmp_rot)
#         
#         # -----------------------------------------------------
#         # COMPUTATION OF RESIDUUM VECTOR
#         # -----------------------------------------------------    
#         # rotate sigma_umat and C_ep back
#         sig_umat[:] = md.stress_routines.rotate_sigma(sig_umat_rot, -alpha, -beta, -gamma)
#         C_ep_tmp[:,:] = md.stress_routines.rotate_material_matrix(C_ep_tmp_rot, -alpha, -beta, -gamma)
# #         print '\t C_ep_tmp: '
# #         print C_ep_tmp
# #         print '\tdifference: '
# #         print C_ep_tmp_rot-C_ep_tmp
# #         print '\tsig_umat: ', sig_umat
#         
#         # compute residuum vector
#         R[:] = sig_umat - sigma_tmp
#         
# #         # debug information
# #         print '\t R: ', R        
#         
#         # compute norm of residuum vector
#         R_norm = np.linalg.norm(R)
#         
#         # debug information
# #         print '\t R_norm: ', R_norm, '\n'
#         
#         # check if convergence criterion is met
#         if R_norm < constants.tolerance:           
#             
#             # norm is small than tolerance,
#             # exit loop
#             break
#         
#         # -----------------------------------------------------
#         # COMPUTATION OF DERIVATIVES OF RESIDUUM VECTOR (dR/dy - jacobian)
#         # -----------------------------------------------------
#         # set jacobian to zero
#         jac[:,:] = 0.
#         
#         # computation of jacobian
#         jac[:,ind_deps_undef] = C_ep_tmp[:,ind_deps_undef]
#         jac[ind_dsig_undef,ind_dsig_undef] = -1.
#         
# #         print 'jac:'
# #         print jac
#         
#         # -----------------------------------------------------
#         # SOLVE SYSTEM OF LINEAR EQUATIONS
#         # -----------------------------------------------------
#         delta_y[:] = np.linalg.solve(jac, -R)
#         
#         # update y
#         y[:] = y + delta_y
#         
# #         debug information
# #         print '\tdelta_y: ', delta_y
# #         print '\ty: ', y
#         
#         # -----------------------------------------------------
#         # UPDATE STRESS AND STRAIN VECTOR WITH VALUES FROM y
#         # AND COMPUTE NEW RESIDUUM
#         # -----------------------------------------------------        
#         # assign stress values from y to sigma_tmp
#         sigma_tmp[ind_dsig_undef] = y[ind_dsig_undef]
# #         print '\t sigma_tmp: ', sigma_tmp
# 
# #         eps_test = np.linalg.solve(C_ep_tmp_rot, sigma_tmp)
# #         print '\teps_test: '
# #         print eps_test
#         
#         # assign strain values from y to eps_tmp
#         eps_tmp[ind_deps_undef] = y[ind_deps_undef]
# #         print '\t eps_tmp: ', eps_tmp
#         
#         # increase counter
#         k+=1
# 
#     
# #     # derive number of inner newton iterations
# #     no_inner_it = statev_tmp[umat_aux_vars['nstatv']-2]
#         
#     # return epsilon, sigma, C_ep_tmp, state variables, number of outer newton iterations
#     # and array with number of inner newton iterations per outer iteration
#     return (sigma_tmp, eps_tmp, statev_tmp, k, no_inner_it, t_all)