'''
Created on Sep 30, 2014

@author: benni
'''
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from ScrolledFrameModule import VerticalScrolledFrame

import numpy as np
from Tkinter import *

        

class FunctionExplorer(Frame):      # class inherits from frame superclass
    
    def __init__(self, labels, minRang, maxRang, alpha_meas, fcu_meas, parent=None):
        
        Frame.__init__(self, parent)        # initialization of parent class Frame
        self.pack(expand=YES, fill=BOTH)    # packing of self
        
        
        # DATA PREPARATION
        #--------------------------
        n = len(labels)                         # number of parameters        
        self.pars = []                          # empty list for initial function parameters
        for i in range(n):                      # create double var for each function par and assign minimum value to it
            self.pars.append(DoubleVar())
            self.pars[i].set(minRang[i]+(maxRang[i]-minRang[i])/2.)
        
          
        self.alpha = np.linspace(0., np.pi/2., 100)         # create alpha array with linearly spaced values between 0 and pi/2
        self.fcu = self.comp_fcu(self.alpha)                # comp fcu values with initial parameter values
        self.alpha_meas = alpha_meas
        self.fcu_meas = fcu_meas
        self.norm_err = self.comp_error()
        
        
        # PLOTTING OF DATA
        #--------------------------
        f = Figure(figsize=(5,3), dpi=100)                  # create matplotlib figure
        a = f.add_subplot(111)                              # add subplot to it
        a.grid()                                            # set grid
        a.plot(alpha_meas, fcu_meas)                        # plot measured values to subplot
        self.line, = a.plot(self.alpha, self.fcu)           # plot function values with intial params to subplot and store function line
        a.axis([0.,np.pi/2.,min(fcu_meas), max(fcu_meas)])  # set x-axis between 0 and pi/2 and y-axis between 0 and 60
        

        # CREATION OF VISIBLE GUI ELEMENTS
        #-----------------------------------
        self.canvas = FigureCanvasTkAgg(f, master=self)         # create canvas object from matplotlib to tkinter class
        self.canvas.show()                                      # show canvas
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=YES)   # pack canvas

        toolbar = NavigationToolbar2TkAgg( self.canvas, self)
        toolbar.update()
        self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=YES)        

        scrolledframe = VerticalScrolledFrame(self)
        scrolledframe.pack(side=TOP, fill=BOTH, expand=YES)
        
        sliderframe = LabelFrame(scrolledframe.interior, text='Params', borderwidth=1)            # create frame for sliders with visible frame borders
        sliderframe.pack(side=TOP, expand=YES)    # pack frame
#         sliderframe.grid(row=0, column=0)   # grid frame
        self.scales = []
        i = 0
        for name in labels:                                     # create scales for each function parameter in frame
            sc = Scale(sliderframe, label=name,
                  variable=self.pars[i],
                  command=self.update_plot,
                  from_=minRang[i], to=maxRang[i],
                  length=400, tickinterval=(maxRang[i]-minRang[i]),
                  resolution=(maxRang[i]-minRang[i])/3000.,
                  orient='horizontal', width=6)                # connect sliders to initial variables values and connect command to plot update function
            sc.grid(row=i, column=0)                           # position sliders in grid
            self.scales.append(sc)
            i+=1
        
        self.error = Label(sliderframe, text='Error: '+str(self.norm_err))
        self.error.grid(row=i)
        
        minmaxFrame = LabelFrame(scrolledframe.interior, text='Ranges', borderwidth=1)
        minmaxFrame.pack(side=TOP, expand=YES)
#         minmaxFrame.grid(row=0, column=1)
        self.minRangChar = []
        self.maxRangChar = []
        for i in range(n):                          # create string var for each function par and assign minimum range value and max range value to it
            self.minRangChar.append(StringVar())
            self.maxRangChar.append(StringVar())
            self.minRangChar[i].set(minRang[i])
            self.maxRangChar[i].set(maxRang[i])
        
        for i in range(n):
            Label(minmaxFrame, text='min '+labels[i]).grid(row=i, column=0, sticky=E, padx=2)
            tmp = Entry(minmaxFrame, textvariable=self.minRangChar[i], width=12)
            tmp.grid(row=i, column=1)
            tmp.bind('<KeyPress-Return>', self.changeLimits)
            Label(minmaxFrame, text='max '+labels[i]).grid(row=i, column=2, sticky=E, padx=2)
            tmp = Entry(minmaxFrame, textvariable=self.maxRangChar[i], width=12)
            tmp.grid(row=i, column=3)
            tmp.bind('<KeyPress-Return>', self.changeLimits)
        
        
    def comp_fcu(self, alpha):       
        vals = []                   # empty list for parameter values
        for par in self.pars:       # extract floats from pars list which contains DoubleVar objects
            vals.append(par.get())

#         return vals[0]*( 1. + vals[1]*(1.-3.*np.cos(alpha)**2) + \
#                              vals[2]*vals[1]**2*(1.-3.*np.cos(alpha)**2)**2 + \
#                              vals[3]*vals[1]**3*(1.-3.*np.cos(alpha)**2)**3 + \
#                              vals[4]*vals[1]**4*(1.-3.*np.cos(alpha)**2)**4)           # compute fcu values and return

        n = len(vals)    
        f_mean = vals[0]
        A1 = vals[1]    
        a = A1*(1.-3.*np.cos(alpha)**2)
        b = A1*(1.-3.*np.cos(alpha)**2)
        for i in range(2,n):
            b = b + vals[i]*a**i
        
        return f_mean * (1.+b)
        
    
    def comp_error(self):
        
        err = self.fcu_meas - self.comp_fcu(self.alpha_meas)
        return np.linalg.norm(err)**2 
        
    
    def update_plot(self, val):
        self.fcu = self.comp_fcu(self.alpha)          # call fcu computation
        self.line.set_ydata(self.fcu)       # update function line
        self.canvas.draw()                  # redraw canvas object
        self.norm_err = self.comp_error()
        self.error.config(text='Error: '+str(self.norm_err))
        
    def changeLimits(self, event):
        minRangNew = []
        for val in self.minRangChar:
            minRangNew.append(float(val.get()))
        maxRangNew = []
        for val in self.maxRangChar:
            maxRangNew.append(float(val.get()))
        
        i=0
        for sc in self.scales:
            sc.config(from_=minRangNew[i], to=maxRangNew[i],
                      tickinterval=maxRangNew[i]-minRangNew[i],
                      resolution=(maxRangNew[i]-minRangNew[i])/3000.)
            i+=1

class FunctionExplorerDP(FunctionExplorer):
    # function explorer class for anisotropic uniaxial compressive strength 
    # based on anisotropic drucker prager criterion with anisotropic cohesion
    # class inherits from FunctionExplorer base class, only the comp_fcu method was redefined
    
    def comp_fcu(self, alpha):       
        vals = []                   # empty list for parameter values
        for par in self.pars:       # extract floats from pars list which contains DoubleVar objects
            vals.append(par.get())

#         return vals[0]*( 1. + vals[1]*(1.-3.*np.cos(alpha)**2) + \
#                              vals[2]*vals[1]**2*(1.-3.*np.cos(alpha)**2)**2 + \
#                              vals[3]*vals[1]**3*(1.-3.*np.cos(alpha)**2)**3 + \
#                              vals[4]*vals[1]**4*(1.-3.*np.cos(alpha)**2)**4)           # compute fcu values and return

        n = len(vals)    
        A_phi = vals[0]
        B_phi_mean = vals[1]
        A1 = vals[2]    
        a = A1*(1.-3.*np.cos(alpha)**2)
        b = A1*(1.-3.*np.cos(alpha)**2)
        for i in range(3,n):
            b = b + vals[i]*a**i
        
        return B_phi_mean*(1.+b)/(np.sqrt(2./3.)-A_phi/3.)


class FunctionExplorerDP_friction(FunctionExplorer):
    # function explorer class for anisotropic uniaxial compressive strength 
    # based on anisotropic drucker prager criterion with anisotropic friction
    # class inherits from FunctionExplorer base class, only the comp_fcu method was redefined
    
    def comp_fcu(self, alpha):       
        vals = []                   # empty list for parameter values
        for par in self.pars:       # extract floats from pars list which contains DoubleVar objects
            vals.append(par.get())

#         return vals[0]*( 1. + vals[1]*(1.-3.*np.cos(alpha)**2) + \
#                              vals[2]*vals[1]**2*(1.-3.*np.cos(alpha)**2)**2 + \
#                              vals[3]*vals[1]**3*(1.-3.*np.cos(alpha)**2)**3 + \
#                              vals[4]*vals[1]**4*(1.-3.*np.cos(alpha)**2)**4)           # compute fcu values and return

        n = len(vals)    
        B_phi = vals[0]
        A_phi_mean = vals[1]
        A1 = vals[2]    
        a = A1*(1.-3.*np.cos(alpha)**2)
        b = A1*(1.-3.*np.cos(alpha)**2)
        for i in range(3,n):
            b = b + vals[i]*a**(i-1)
        
        return B_phi/(np.sqrt(2./3.)-A_phi_mean*(1.+b)/3.)



class FunctionExplorerGao(FunctionExplorer):
    # function explorer class for anisotropic uniaxial compressive strength 
    # based on anisotropic modified Leon criterion with anisotropic fcu based on Gaos approach
    # class inherits from FunctionExplorer base class, only the comp_fcu method was redefined
    
    def comp_fcu(self, alpha):       
        vals = []                   # empty list for parameter values
        for par in self.pars:       # extract floats from pars list which contains DoubleVar objects
            vals.append(par.get())

        n = len(vals)
        f_mean = vals[0]
        
        #A = 3./2.*np.sin(alpha)**2 - 1.            # gaos approach based on positive sign convention for compressive stresses
        #A = (3./2.*np.sin(alpha)**2 - 1.)*(-1.)     # adjustment for negative sign convention for compressive stresses and assuming that d>0
        A = -0.5*(1.-3.*np.cos(alpha)**2)     # adjustment for negative sign convention for compressive stresses and assuming that d>0
        
        g = 0.
        for i in range(1,n):
            #g = g + vals[i]*(1.+A)**i            # gaos approach based on positive sign convention for compr. stresses 
                                                  # (addition of 1 to A assures that g is always zero for horizontal foliation planes, ie alpha=0)
            g = g + vals[i]*(-1.+A)**i            # adjustment for negative sign convention of compr. stresses
                                                  # assures that g is zero for alpha=0)
            
        return f_mean*np.exp(g)


if __name__ == '__main__':
    
    root = Tk()

    # pietruczcack
#     alpha_meas = np.array([0., 0.25924, 0.523353, 0.784487, 1.0467, 1.21965, 1.5707963267949])   # values of inclination angle alpha
#     fcu_meas = np.array([46., 30.78, 24.7746, 22.2834, 20.2456, 29.8398, 34.842])           # values of fcu depending on alpha
#     fu = FunctionExplorer(['f_mean', 'A1', 'b1', 'b2', 'b3'], 
#                           [20.0, 0.0, 0., 0., 0.], [50.0, 0.2, 100., 1000., 4000.], 
#                           alpha_meas, fcu_meas, root)
    
    #cho 2012 gneiss
#     alpha_meas = np.pi/180. * np.array([0., 15., 30., 45., 60., 75., 90.])   # values of inclination angle alpha
#     fcu_meas = np.array([183.8, 202.4, 158.8, 111.4, 110.4, 98.6, 223.2])           # values of fcu depending on alpha
# #     fu = FunctionExplorer(['f_mean', 'A1', 'b1', 'b2', 'b3'], 
# #                           [90.0, -0.1, -1000., -100000., -800000.], [230.0, 0.1, 1000., 100000., 800000.], 
# #                           alpha_meas, fcu_meas, root)
#     fu = FunctionExplorer(['f_mean', 'A1', 'b1', 'b2', 'b3', 'b4'], 
#                           [90.0, -0.4, -20., -50., -200., -500.], [230.0, 0.4, 20., 50., 200., 500.], 
#                           alpha_meas, fcu_meas, root)


    #cho 2012 gneiss with scalar anisotropic variable according to Gao et al
#     alpha_meas = np.pi/180. * np.array([0., 15., 30., 45., 60., 75., 90.])   # values of inclination angle alpha
#     fcu_meas = np.array([183.8, 202.4, 158.8, 111.4, 110.4, 98.6, 223.2])           # values of fcu depending on alpha
#     fu = FunctionExplorerGao(['f_mean', 'A1', 'b1', 'b2', 'b3', 'b4'], 
#                           [90.0, -1., -2, -2., -0.5, -0.1], [230.0, 1., 2., 2., 0.5, 0.1], 
#                           alpha_meas, fcu_meas, root)


    #cho 2012 shale
#     alpha_meas = np.pi/180. * np.array([0., 15., 30., 45., 60., 75., 90.])   # values of inclination angle alpha
#     fcu_meas = np.array([89.2, 91.0, 83.4, 64.8, 53.5, 96.5, 126.2])           # values of fcu depending on alpha
#     fu = FunctionExplorer(['f_mean', 'A1', 'b1', 'b2', 'b3'], 
#                           [50., -0.1, -10000., -100000., -800000.], [130., 0.1, 10000., 100000., 800000.], 
#                           alpha_meas, fcu_meas, root)


    #cho 2012 shale with scalar anisotropic variable according to Gao et al
#     alpha_meas = np.pi/180. * np.array([0., 15., 30., 45., 60., 75., 90.])   # values of inclination angle alpha
#     fcu_meas = np.array([89.2, 91.0, 83.4, 64.8, 53.5, 96.5, 126.2])           # values of fcu depending on alpha
#     fu = FunctionExplorerGao(['f_mean', 'A1', 'b1', 'b2', 'b3'], 
#                           [30., -1., -2, -2., -0.5, -0.1], [120., 1., 2., 2., 0.5, 0.1], 
#                           alpha_meas, fcu_meas, root)


    #cho 2012 schist
#     alpha_meas = np.pi/180. * np.array([0., 15., 30., 45., 60., 75., 90.])   # values of inclination angle alpha
#     fcu_meas = np.array([85., 119.6, 48.1, 8.9, 27.7, 27.5, 124.7])           # values of fcu depending on alpha
#     fu = FunctionExplorer(['f_mean', 'A1', 'b1', 'b2', 'b3'], 
#                           [0., -0.1, -50000., -500000., -800000.], [130., 0.1, 50000., 500000., 800000.], 
#                           alpha_meas, fcu_meas, root)


#     #cho 2012 schist with scalar anisotropic variable according to Gao et al
#     alpha_meas = np.pi/180. * np.array([0., 15., 30., 45., 60., 75., 90.])   # values of inclination angle alpha
#     fcu_meas = np.array([85., 119.6, 48.1, 8.9, 27.7, 27.5, 124.7])           # values of fcu depending on alpha
#     fu = FunctionExplorerGao(['f_mean', 'A1', 'b1', 'b2', 'b3'], 
#                           [30., -1., -2, -2., -0.5, -0.1], [120., 1., 2., 2., 0.5, 0.1], 
#                           alpha_meas, fcu_meas, root)


#     #graz 2014 innsbruck quartzphyllite
#     alpha_meas = np.pi/180. * np.array([0., 60., 90.])      # values of inclination angle alpha
#     fcu_meas = np.array([112.311097, 24.5162162, 16.99811333])           # values of fcu depending on alpha
#     fu = FunctionExplorer(['f_mean', 'A1', 'b1', 'b2'], 
#                           [0., -0.1, -50000., -10000.], [130., 0.1, 50000., 10000.], 
#                           alpha_meas, fcu_meas, root)
    
#     #graz 2014 innsbruck quartzphyllite DP test
#     alpha_meas = np.pi/180. * np.array([0., 60., 90.])      # values of inclination angle alpha
#     fcu_meas = np.array([112.311097, 24.5162162, 16.99811333])           # values of fcu depending on alpha
#     fu = FunctionExplorerDP(['A_phi', 'B_phi_mean', 'A1', 'b1', 'b2'], 
#                           [0., 0., -0.1, -50000., -10000.], [2.5, 130., 0.1, 50000., 10000.], 
#                           alpha_meas, fcu_meas, root)

    #graz 2014 innsbruck quartzphyllite DP test
    alpha_meas = np.pi/180. * np.array([0., 60., 90.])      # values of inclination angle alpha
    fcu_meas = np.array([112.311097, 24.5162162, 16.99811333])           # values of fcu depending on alpha
    fu = FunctionExplorerDP_friction(['B_c', 'A_phi_mean', 'A1', 'b1', 'b2'], 
                          [0., 0., -2., -10., -10.], [100., 2.5, 2., 10., 10.], 
                          alpha_meas, fcu_meas, root)

#     #graz 2014 innsbruck quartzphyllite with scalar anisotropic variable according to Gao et al
#     alpha_meas = np.pi/180. * np.array([0., 60., 90.])      # values of inclination angle alpha
#     fcu_meas = np.array([112.311097, 24.5162162, 16.99811333])           # values of fcu depending on alpha
#     fu = FunctionExplorerGao(['f_mean', 'c1', 'c2', 'c3'], 
#                           [0., -10., -5., -5.], [130., 10., 5., 5.], 
#                           alpha_meas, fcu_meas, root)
    
    root.mainloop()