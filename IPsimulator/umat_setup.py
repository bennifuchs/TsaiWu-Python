'''
Created on Apr 11, 2014
defines function for setup of all variables needed by UMAT
@author: benni
'''
import numpy as np
from error_functions import MaterialInputError, MaterialParameterError


def set_variables(mat_file_name):
    # reads material parameter file
    # and initializes auxiliary variables needed by UMAT
    # returns dictionary with all neccessary data
    
    # read and load material data
    mat_data = np.array(read_material_data(mat_file_name), order='F')
        
    # create and initialize auxiliary variables    
    # -----------------------------------------
    
    n_yf = mat_data[1]              # number of yield functions
    n_int_vars = mat_data[2]        # number of internal variables    
    n_int_vars_dam = mat_data[3]    # number of internal damage variables
    n_omega = mat_data[4]           # number of damage variables
    
    nshr = 3                                                # number of shear stress components (3D)
    ndi = 3                                                 # number of direct stress components (3D)
    ntens = nshr + ndi                                      # number of total stress components

    nstatv = 1+6+6+ \
                n_yf+n_int_vars+n_yf+ \
                n_int_vars_dam+6+n_omega+ \
                36+ \
                1+1+1                                       # number of state variables ( n_act_yf, eps_pl, eps_el, 
                                                            #                             delta_lambda, alpha, yield_flags,
                                                            #                             alpha_dam, delta_eps_pl, omega,
                                                            #                             Cep,
                                                            #                             n_yf_comb_tested, n_inner_newton_it,
                                                            #                             status_var)

    props = mat_data[1:]             # material parameters array
    nprops = props.shape[0]          # number of material parameters
    statev = np.zeros((nstatv,), order='F')     # define and initialize state variables array with zeros    
    
    # auxiliary input/output variables which are not important for the current material model
    sse = 1.
    spd = 1.
    scd = 1. 
    time = np.array([1., 1.], order='F')
    dtime = 1.
    temp = 1.
    dtemp = 1.
    predef = 1.
    dpred = 1.
    kinc = 1
    kstep = 1
    kspt = 1
    layer = 1
    npt = 1
    noel = 1
    dfgrd0 = np.identity(3)
    dfgrd1 = np.identity(3)
    celent = 108.
    drot = np.identity(3)
    coords = np.array([0., 0., 0.], order='F')
    
    # set material name according to the material number
    if mat_data[0] == 200:
        cmname = 'mod_leon'
    elif mat_data[0] == 300:
        cmname = 'mod_leon_damage'
    elif mat_data[0] == 400:
        cmname = 'mod_leon_damage_rock'
    elif mat_data[0] == 500:
        cmname = 'mod_leon_damage_plast_ansitr_rock'
    elif mat_data[0] == 600:
        cmname = 'mod_leon_damage_el-pl_ansitr_rock'
    elif mat_data[0] == 700:
        cmname = 'tsai_wu_transv_iso'
    
    # create dictionary for all auxiliary variables needed by umat
    aux_vars = {
                     'nshr' : nshr,
                     'ndi' : ndi,
                     'ntens' : ntens,
                     'nstatv' : nstatv,
                     'nprops' : nprops,
                     'props' : props,
                     'sse' : sse,
                     'spd' : spd,
                     'scd' : scd,
                     'time' : time,
                     'dtime' : dtime,
                     'temp' : temp,
                     'dtemp' : dtemp,
                     'predef' : predef,
                     'dpred' : dpred,
                     'kinc' : kinc,
                     'kstep' : kstep,
                     'kspt' : kspt,
                     'layer' : layer,
                     'npt' : npt,
                     'noel' : noel,
                     'dfgrd0' : dfgrd0,
                     'dfgrd1' : dfgrd1,
                     'celent' : celent,
                     'drot' : drot,
                     'coords' : coords,
                     'cmname' : cmname
                     }
    
    return aux_vars

def read_material_data(mat_file_name):
    # function that reads the material input und the analysis type from an input file
    # and returns the required material parameters depending on the analysis type
    # input:
    # mat_file_name ...    string with name of the material file
    # output:
    # parameters ...        dictionary with number aof the analysis type and the required material parameters
    
    # open input file
    fobj = open(mat_file_name, "r")
    
    # create empty dictionary
    data = {}
    
    # read lines of input file
    for line in fobj:
        
        # check for comment line
        # and empty line
        # if yes, skip line
        if (line[0] == "#" or line.isspace()):
            continue
        
        # split line at whitespaces and store line segments in temporary list
        tmp = line.split()
        
        # check for incomplete line, parameter-value-pair
        try:
            # add parameter-value-pair to dictionary
            data.update({tmp[0] : float(tmp[1])})
        except IndexError:
            print "Error: missing parameter/value pair!! \nchange inputfile!!!"
            fobj.close()
            raise IndexError
    
    # close input file    
    fobj.close()
    
    # check for missing analysis type declaration
    try:        

        # isotropic modified leon with damage in rock notation , i.e. m0 and e instead of ftu and fbu
        if data["analysis_type"] == 400:
            #check if required parameters for linear elastic ideal plastic analysis are specified
            try:

                #return required material parameters as tuple
                return (data["analysis_type"],
                         data["n_yf"], data["n_int_vars"], data["n_int_vars_dam"], data["n_omega"], 
                         data["E"], data["nu"],
                         data["fcy"], data["fcu"], data["m0"], data["e"],
                         data["mg0"], 
                         data["Ah"], data["Bh"], data["Ch"], data["Dh"],
                         data["GfI"], data["As"],
                         data["Ad"], data["Bd"])
            
            except KeyError:
                print "Error: Some material parameters for isotropic modified leon model with damage, formulated for rock, are missing!! \
                            \nrequired parameters: n_yf, n_int_vars, n_int_vars_dam, n_omega, E, nu, fcy, fcu, m0, e, mg0, Ah, Bh, Ch, Dh, Gf1, As, Ad and Bd \
                            \nadd parameters to input file!!"
                raise MaterialParameterError

        # plastic anisotropic, elastic isotropic modified leon with damage in rock notation , i.e. m0 and e instead of ftu and fbu
        elif data["analysis_type"] == 500:
            #check if required parameters for linear elastic ideal plastic analysis are specified
            try:

                #return required material parameters as tuple
                return (data["analysis_type"],
                         data["n_yf"], data["n_int_vars"], data["n_int_vars_dam"], data["n_omega"], 
                         data["E"], data["nu"],
                         data["fcy"], data["fcu_mean"], data["m0"], data["e"],
                         data["mg0"], 
                         data["Ah"], data["Bh"], data["Ch"], data["Dh"],
                         data["GfI"], data["As"],
                         data["Ad"], data["Bd"],
                         data["Aa1"], data["ba1"], data["ba2"], data["ba3"])
            
            except KeyError:
                print "Error: Some material parameters for plastic anisotropic modified leon model with damage, formulated for rock, are missing!! \
                            \nrequired parameters: n_yf, n_int_vars, n_int_vars_dam, n_omega, E, nu, fcy, fcu/fcu_mean, m0, e, mg0, Ah, Bh, Ch, Dh, Gf1, As, Ad and Bd \
                            \nand plastic anisotropic parameters Aa1, ba1, ba2, ba3 \
                            \nadd parameters to input file!!"
                raise MaterialParameterError
        
        # plastic and elastic anisotropic modified leon with damage in rock notation , i.e. m0 and e instead of ftu and fbu
        elif data["analysis_type"] == 600:
            #check if required parameters for linear elastic ideal plastic analysis are specified
            try:

                #return required material parameters as tuple
                return (data["analysis_type"],
                         data["n_yf"], data["n_int_vars"], data["n_int_vars_dam"], data["n_omega"], 
                         data["E1"], data["nu1"],
                         data["fcy"], data["fcu_mean"], data["m0"], data["e"],
                         data["mg0"], 
                         data["Ah"], data["Bh"], data["Ch"], data["Dh"],
                         data["GfI"], data["As"],
                         data["Ad"], data["Bd"],
                         data["Aa1"], data["ba1"], data["ba2"], data["ba3"],
                         data["E2"], data["G2"], data["nu2"])
            
            except KeyError:
                print "Error: Some material parameters for elastic-plastic anisotropic modified leon model with damage, formulated for rock, are missing!! \
                            \nrequired parameters: n_yf, n_int_vars, n_int_vars_dam, n_omega, E1, nu1, fcy, fcu/fcu_mean, m0, e, mg0, Ah, Bh, Ch, Dh, Gf1, As, Ad and Bd \
                            \nand plastic anisotropic parameters Aa1, ba1, ba2, ba3 \
                            \nand elastic anisotropic parameters E2, G2, nu2 \
                            \nadd parameters to input file!!"
                raise MaterialParameterError
        
        # transversely isotropic linear elastic, ideal plastic tsai-wu material model  
        elif data["analysis_type"] == 700:
            
            try:
                
                return(data["analysis_type"],
                       data["n_yf"], data["n_int_vars"], data["n_int_vars_dam"], data["n_omega"], 
                       data["E1"], data["nu1"], data["E2"], data["G2"], data["nu2"],
                       data["F2"], data["F3"],
                       data["F22"], data["F33"], data["F44"], data["F12"], data["F23"])

            except KeyError:
                print "Error: some parameters for the transversely isotropic linear elastic, ideal plastic tsai-wu material model are missing!! \
                            \nrequired parameters: n_yf, n_int_vars, n_int_vars_dam, n_omega, E1, nu1, E2, G2, nu2 \
                            \nand plastic anisotropic parameters F2, F3, F22, F33, F44, F12, F23 \
                            \nadd parameters to input file!!"
                raise MaterialParameterError
            
        
        #unsupported material number
        else:
            print "Error: The specified analysis_type number is currently not supported!! \
                    \nchange the analysis_type number in the input file"
            raise MaterialInputError
    
    except KeyError:
        print "Error: Analysis type is missing!! \
                \nadd analysis type to input file!!!"
        raise KeyError


def read_loadhistory(loadhistory_file_name):
    pass