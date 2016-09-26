'''
Created on Jun 15, 2016

@author: benni
'''

import numpy as np
import numpy.ma as ma
import scipy as scp
import scipy.optimize as scop
import matplotlib.pyplot as plt

from IPsimulator import tsaiwu as ml



class yfTsaiWu():
    
    def __init__(self, sigma, stress_csys, 
                 Fi, Fi_csys, 
                 Fij, Fij_csys, 
                 beta, internal_csys='local'):        
        
        
        # define class members and initialize them
        self.sigma = self.stress_init_(sigma)
           
        # define and initialize remaining class members    
        self.Fi = np.zeros(6)
        self.Fij = np.zeros((6,6))        
        self.csys = internal_csys
        self.beta = beta        # always defined as rotation TOWARDS the internal csys

        self.__stress_csys_inp = stress_csys
        self.__Fi_csys_inp = Fi_csys
        self.__Fij_csys_inp = Fij_csys
            
        
        # in case of internal representation in local csys
        # transform all input variables into this csys
        if internal_csys == 'local':
            
            if stress_csys == 'local':
                self.sigma[:] = sigma
            elif stress_csys == 'global':
                self.sigma[:] = self.stressTransform(sigma, beta)
                
            if Fi_csys == 'local':
                self.Fi[:] = Fi
            elif Fi_csys == 'global':
                self.Fi[:] = self.FiTransform(Fi, beta)
                
            if Fij_csys == 'local':
                self.Fij[:,:] = Fij
            elif Fij_csys == 'global':
                self.Fij[:,:] = self.FijTransform(Fij, beta)
                
        # in case of internal representation in global csys
        # transform all input variables into this csys
        elif internal_csys == 'global':
            
            if stress_csys == 'local':
                self.sigma[:] = self.stressTransform(sigma, beta)
            elif stress_csys == 'global':
                self.sigma[:] = sigma
            
            if Fi_csys == 'local':
                self.Fi[:] = self.FiTransform(Fi, beta)
            elif Fi_csys == 'global':
                self.Fi[:] = Fi
                
            if Fij_csys == 'local':
                self.Fij[:,:] = self.FijTransform(Fij, beta)
            elif Fij_csys == 'global':
                self.Fij[:,:] = Fij
        
        # compute yf value and store as class member    
        self.yf = self.comp_yf(self.sigma, self.Fi, self.Fij)
                

    def stress_init_(self,sigma):
        # initialize internal stress variable as 1D-numpy array with 6 entries
        # this function is intended for adjustments in derived classes
        
        return np.zeros(6)

            
    def stressLocal(self):
        # returns sigma in local coordinates
        
        if self.csys == 'local':
            return self.sigma
        elif self.csys == 'global':
            return self.stressTransform(self.sigma, -self.beta)

    
    def stressGlobal(self):
        # returns sigma in global coordinates
        
        if self.csys == 'local':
            return self.stressTransform(self.sigma, -self.beta)
        elif self.csys == 'global':
            return self.sigma

        
    def FiLocal(self):
        # returns Fi in local coordinates

        if self.csys == 'local':
            return self.Fi
        elif self.csys == 'global':
            return self.FiTransform(self.Fi, -self.beta)

        
    def FiGlobal(self):
        # returns Fi in global coordinates
        
        if self.csys == 'local':
            return self.FiTransform(self.Fi, -self.beta)
        elif self.csys == 'global':
            return self.Fi

        
    def FijLocal(self):
        # returns Fij in local coordinates
        
        if self.csys == 'local':
            return self.Fij
        elif self.csys == 'global':
            return self.FijTransform(self.Fij, -self.beta)

        
    def FijGlobal(self):    
        # returns Fij in global coordinates
        
        if self.csys == 'local':
            return self.FijTransform(self.Fij, -self.beta)
        elif self.csys == 'global':
            return self.Fij
        
    
    def stressTransform(self, sig, beta):
        # perform coord transformation of stress tensor in vectorial notation
        
        # call fortran stress transformation routine
        sig_rot = ml.stress_routines.rotate_sigma(sig, 0., beta, 0.)
        
        # return transformed stress tensor in vectorial notation
        return sig_rot
        
    
    def FiTransform(self, Fi, beta):
        # perform coord transformation of structural vector Fi
        
        # transform structural vector with same relationship as strain vector from fortran
        Fi_rot = ml.stress_routines.rotate_epsilon(Fi, 0., beta, 0.)
        
        # return transformed structural vector
        return Fi_rot
       

    def FijTransform(self, Fij, beta):
        # perform coord transformation of structural tensor Fij
        
        # transform structural tensor with same relationship as inverse material matrix from fortran
        Fij_rot = ml.stress_routines.rotate_inv_material_matrix(Fij, 0., beta, 0.)
        
        # return transformed structural tensor
        return Fij_rot
        
    
    def structTensTransform(self, Fi, Fij,
                            beta):
        
        # transform structural vector with some relationship as strain vector
        Fi_rot = ml.stress_routines.rotate_epsilon(Fi, 0., beta, 0.)
        
        # transform structural tensor with same relationship as inverse material matrix
        Fij_rot = ml.stress_routines.rotate_inv_material_matrix(Fij, 0., beta, 0.)
        
        return (Fi_rot, Fij_rot)
    
    
    def comp_yf(self, sig, Fi, Fij):
                
        # computation of yf via vector-matrix and matrix-matrix operations
        yf = np.dot(Fi, sig) + np.dot(sig,np.dot(Fij,sig)) - 1
        
        return yf


    def update(self, **kwargs):
        # updates yield function object
        # by recalling the __init__ routine
        # kwargs must be some of the arguments (passed in as keyword arguments) required by the __init__ routine
        
        # check if kwargs are empty
        if len(kwargs) == 0:
            return
        
        # define default keyword args dictionary
        default_kwargs = {'sigma': self.sigma, 'stress_csys': self.__stress_csys_inp, 
                          'Fi': self.Fi, 'Fi_csys': self.__Fi_csys_inp, 
                          'Fij': self.Fij, 'Fij_csys': self.__Fij_csys_inp,
                          'beta': self.beta, 'internal_csys': self.csys}
        
        # update default keyword dictionary with new passed in kwargs
        default_kwargs.update(kwargs)
        
        # call __init__ routine with updated default keyword dictionary
        self.__init__(**default_kwargs)


class yfTsaiWuLists(yfTsaiWu):
    # same functionality as yfTsaiWu class but intended for stress input as list of numpy arrays
    # hence, the computation of yf values for multiple stress states is possible
    # this is useful for contour plots via meshgrid arrays and for optimization
    # format of sigma passed into __init__ method: sigma = [array, array, array, array, array, array,]
    #                                                         ^      ^      ^      ^      ^      ^
    #                                                         |      |      |      |      |      |
    #                                                       sig11  sig22  sig33  sig12  sig13  sig23
    
            
    def stress_init_(self, sigma):
        # initialize internal stress variable as list of same objects as the ones in the sigma input list and set to zero
        # those objects are intended to be numpy arrays
        
        return [component*0. for component in sigma]
    
    
    def stressTransform(self, sig, beta):
        # perform coordinate transformation for stress as 6 items list of numpy arrays
        # with the i-th array containing multiple values of the i-th stress component
        
        # number of elements per numpy array
        nel = sig[0].size
        
        # create and initialize sig list for output of same size and type as input sig
        sig_rot_list = self.stress_init_(sig)
        
        # create and initialize temporary stress variables
        # as 1D numpy arrays of 
        sig_tmp = np.zeros(6)
        sig_rot_tmp = np.zeros(6)
        
        if (type(beta) == float or type(beta) == np.float64):

            # iterate through all elements of flattened stress component arrays in input list
            for i in range(nel):
                
                # assign i-th entry of each flattened stress component array to the
                # right position in the temporary array
                sig_tmp[:] = [sig[0].flat[i], sig[1].flat[i], sig[2].flat[i], 
                              sig[3].flat[i], sig[4].flat[i], sig[5].flat[i]]
                
                # perform coordinate transformation of temporary array
                sig_rot_tmp[:] = ml.stress_routines.rotate_sigma(sig_tmp, 0., beta, 0.)
                
                # assign transformed stress components to the right positions in the output list
                for k in range(6):                
                    sig_rot_list[k].flat[i] = sig_rot_tmp[k]
        
        elif type(beta) == np.ndarray:

            # iterate through all elements of flattened stress component arrays in input list
            for i in range(nel):
                
                # assign i-th entry of each flattened stress component array to the
                # right position in the temporary array
                sig_tmp[:] = [sig[0].flat[i], sig[1].flat[i], sig[2].flat[i], 
                              sig[3].flat[i], sig[4].flat[i], sig[5].flat[i]]
                
                # perform coordinate transformation of temporary array
                sig_rot_tmp[:] = ml.stress_routines.rotate_sigma(sig_tmp, 0., beta.flat[i], 0.)
                
                # assign transformed stress components to the right positions in the output list
                for k in range(6):                
                    sig_rot_list[k].flat[i] = sig_rot_tmp[k]

#         elif type(beta) == list:
# 
#             # iterate through all elements of flattened stress component arrays in input list
#             for i in range(nel):
#                 
#                 # assign i-th entry of each flattened stress component array to the
#                 # right position in the temporary array
#                 sig_tmp[:] = [sig[0].flat[i], sig[1].flat[i], sig[2].flat[i], 
#                               sig[3].flat[i], sig[4].flat[i], sig[5].flat[i]]
#                 
#                 # perform coordinate transformation of temporary array
#                 sig_rot_tmp[:] = ml.stress_routines.rotate_sigma(sig_tmp, 0., beta[i], 0.)
#                 
#                 # assign transformed stress components to the right positions in the output list
#                 for k in range(6):                
#                     sig_rot_list[k].flat[i] = sig_rot_tmp[k]
            
        
        # return output list with transformed stress components in the list entry arrays                        
        return sig_rot_list
    
    
    def comp_yf(self, sig, Fi, Fij):
        
        s11, s22, s33, s12, s13, s23 = sig
                         
        # computation of yf via multiplication and summation of single stress components and structural tensor components
        yf = Fi[0]*s11 + Fi[1]*s22 + Fi[2]*s33 + Fi[3]*s12 + Fi[4]*s13 + Fi[5]*s23 +\
            Fij[0,0]*s11**2 + Fij[1,1]*s22**2 + Fij[2,2]*s33**2 +\
            Fij[3,3]*s12**2 + Fij[4,4]*s13**2 + Fij[5,5]*s23**2 +\
            2*(Fij[0,1]*s11*s22 + Fij[0,2]*s11*s33 + Fij[0,3]*s11*s12 + Fij[0,4]*s11*s13 + Fij[0,5]*s11*s23 +\
               Fij[1,2]*s22*s33 + Fij[1,3]*s22*s12 + Fij[1,4]*s22*s13 + Fij[1,5]*s22*s23 + \
               Fij[2,3]*s33*s12 + Fij[2,4]*s33*s13 + Fij[2,5]*s33*s23 +\
               Fij[3,4]*s12*s13 + Fij[3,4]*s12*s23 +\
               Fij[4,5]*s13*s23) - 1
         
        return yf

    
    def plot_yf_contour(self,index1, index2, scale1=1., scale2=1.,
                        style='contour', scaled=False, color=None, **kwargs):
        # plots yf values as contour plot over grid of stress components of index1 and index2
        # only works if the tow stress components are rank 2 arrays with the same ranges
        # should be tested at the beginning of the function (STILL TO DO)
        
        
#         # check if usage of non default csys is specified
#         if 'csys' in kwargs:
#             
#             if kwargs['csys'] == self.csys:
#                 pass
                
            
        # stress index labels used for axis labeling in plots  
        ind_labels = ['11', '22', '33', '12', '13', '23']
       
        # find min and max values in stress arrays
        range_x = [self.sigma[index1].min()*scale1, self.sigma[index1].max()*scale1]
        range_y = [self.sigma[index2].min()*scale2, self.sigma[index2].max()*scale2]
#         range_min = np.amin([self.sigma[index1].min()*scale1, self.sigma[index2].min()*scale2])
#         range_max = np.amax([self.sigma[index1].max()*scale1, self.sigma[index2].max()*scale2])     
        
        # axis labels for plots
        xlabel='$\sigma_{' + ind_labels[index1] + '}$'
        ylabel='$\sigma_{' + ind_labels[index2] + '}$'
        
        # create plot objects
#         plt.title(r"Compr. meridian of yf, $\beta=%3.0f^{\circ}$" % self.beta_deg, family='serif', size=12)
#         plt.title(r"Compr. meridian of yf, $\beta=%3.0f^{\circ}$" % self.beta_deg)
        
        if style == 'contour':
            # default plot contour line of yf=0
            tmp = plt.contour(self.sigma[index1]*scale1, self.sigma[index2]*scale2, self.yf, 
                              levels=[0], antialiased=True, colors=color)
        elif style == 'contourf':
            # plot filled contours
            tmp = plt.contourf(self.sigma[index1]*scale1, self.sigma[index2]*scale2, self.yf, 
                               levels=np.linspace(-50.,0.,20), antialiased=True)
            # create colorbar
            cb = plt.colorbar(mappable=tmp)
            cb.set_label('$f$')
#             cb.ax.set_yticklabels(np.array(cb.ax.get_yticks())*(-50.), size=10, family='cmr10')
        
        # set axes labels
#         plt.xlabel(xlabel, family='serif', size=12)
#         plt.ylabel(ylabel, family='serif', size=12)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        
        # plot hydrostatic axis and zero axes
        plt.plot(np.array(range_x), np.array(range_y), color='black')
        plt.axhline(0., color='black', linestyle='--')      # zero x-axis
        plt.axvline(0., color='black', linestyle='--')      # zero y-axis
#         plt.plot(np.array([range_min, range_max])*np.sqrt(2), np.array([range_min, range_max]), color='black')
#         plt.axhline(0., color='black', linestyle='--')      # zero x-axis
#         plt.axvline(0., color='black', linestyle='--')      # zero y-axis
        
        # scale axes if desired
        if scaled == True:
            plt.axis('scaled')
        
        # set x- and y-limits
        plt.xlim(range_x)
        plt.ylim(range_y)
#         plt.xlim(range_min, range_max )
#         plt.ylim(range_min, range_max)
        plt.ticklabel_format(style='sci', scilimits=(0,0), axis='both')
        
        # set tick label size and font
#         a = plt.gca()
#         a.set_xticklabels(a.get_xticks(), size=10, family='serif')
#         a.set_yticklabels(a.get_yticks(), size=10, family='serif')
        
        return tmp
        
            
    
if __name__ == '__main__':
    
#     # loading angle in rad
#     beta = np.pi/180.*60.
#     
#     # test structural vector
#     strucVec = np.array((10., 10., 20., 0., 0., 0.))
# 
#     # test structural tensor
#     strucTens2 = np.zeros((6,6))
#     strucTens2[0:3,0:3] = 1.
#     strucTens2[2,2] = 1.
#     strucTens2[3,3] = 1.
#     strucTens2[4,4] = 1.
#     strucTens2[5,5] = 1.
#     
#     # generate object from TsaiWu class
#     yf1 = yfTsaiWu([0., 0., -10., 0., 0., 0.], 'global', 
#               strucVec, 'local', 
#               strucTens2, 'local', 
#               beta)
#      
#     print 'internal representation csys of yf1: '
#     print yf1.csys, '\n'
#      
#     print 'internal representation of sigma: '
#     print yf1.sigma, '\n'
#      
#     print 'internal representation of Fi: '
#     print yf1.Fi, '\n'
#      
#     print 'internal representation of Fij: '
#     print yf1.Fij, '\n'

    
    n = 100
    range_max = 1.
    range_min = -10.
    range1 = np.linspace(range_max, range_min, n)
    
    sig0_mesh, sig33_mesh = np.meshgrid(range1, range1)
    
    # loading angle in rad
    beta = np.pi/180.*90.
    
    # structural vector
    F2 = 3.
    F3 = 10.
    strucVec = np.array([F2, F2, F3, 0., 0., 0.])      # transverse isotropy
    
    # structural tensor
    F12 = -5.
    F23 = -6.
    F22 = -F12-F23
    F33 = -2.*F23
    F44 = 0.
    strucTens = np.array([[F22, F12, F23, 0., 0., 0.],
                          [F12, F22, F23, 0., 0., 0.],
                          [F23, F23, F33, 0., 0., 0.],
                          [0., 0., 0., F44, 0., 0.],
                          [0., 0., 0., 0., F44, 0.],
                          [0., 0., 0., 0., 0., 2.*(F22-F12)]])  # transverse isotropy
    

    yfList = yfTsaiWuLists([sig0_mesh, sig0_mesh, sig33_mesh, 
                            np.zeros_like(sig0_mesh), np.zeros_like(sig0_mesh), np.zeros_like(sig0_mesh)], 'global', 
                           strucVec, 'local',
                           strucTens, 'local',
                           beta, internal_csys='global')
    
    TestContour = yfList.plot_yf_contour(0, 2, scale1= np.sqrt(2.), scaled=True)
    
    # test to access coordinate list defining the zero contour of the contour plot object
    # -----------------------------------------------------------------------------------
    print TestContour.__doc__
#     print TestContour.collections[0].get_paths()[0].vertices
    print TestContour.allsegs[0][0]
    # -----------------------------------------------------------------------------------
    
    # structural vector
    F2 = 3.
    F3 = 10.
    strucVec = np.array([F2, F2, F3, 0., 0., 0.])      # transverse isotropy
    # structural tensor
    F12 = 5.
    F23 = 6.
    F22 = -F12-F23
    F33 = -2.*F23
    F44 = 0.
    strucTens = np.array([[F22, F12, F23, 0., 0., 0.],
                          [F12, F22, F23, 0., 0., 0.],
                          [F23, F23, F33, 0., 0., 0.],
                          [0., 0., 0., F44, 0., 0.],
                          [0., 0., 0., 0., F44, 0.],
                          [0., 0., 0., 0., 0., 2.*(F22-F12)]])  # transverse isotropy
    yfList.update(Fi=strucVec,Fij=strucTens)
    yfList.plot_yf_contour(0, 2, scale1= np.sqrt(2.), scaled=True, color='r')

#     strucTens[0,1] = -7.
#     strucTens[1,0] = -7.
#     yfList.update(Fij=strucTens)
#     yfList.plot_yf_contour(0, 2, scale1= np.sqrt(2.), scaled=True, color='g')
#     
#     strucTens[0,1] = -8.
#     strucTens[1,0] = -8.
#     yfList.update(Fij=strucTens)
#     yfList.plot_yf_contour(0, 2, scale1= np.sqrt(2.), scaled=True, color='orange')
    
    plt.show()        
    
        