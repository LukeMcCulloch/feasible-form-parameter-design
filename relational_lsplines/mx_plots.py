# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 18:45:16 2015

@author: lukemcculloch

souces:
    http://www.randalolson.com/2014/06/28/
        how-to-make-beautiful-data-visualizations-in-python-with-matplotlib/
"""
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.axislines import SubplotZero
import matplotlib.patches as patches
import matplotlib.path    as path

from polyMaker import polymake
#from polyMaker_traditional_stern import poly_prep, polymake

import curve as spline
#
from ADILS import interval_bounds
from ADILS import Lagrangian
from ADILS import lagrangian_bounds
from ADILS import IntervalLagrangeSpline
from FormParameter import FormParameterDict
#

def default_input( message, defaultVal ):
    """http://stackoverflow.com/
    questions/5403138/how-to-set-a-default-string-for-raw-input
    """
    if defaultVal:
        return raw_input( "%s [%s]:" % (message,defaultVal) ) or defaultVal
    else:
        return raw_input( "%s " % (message) )
        
def check_collinear(vertices):
    v1 = list(vertices[0])
    v2 = list(vertices[1])
    v3 = list(vertices[2])
    v1.insert(0,1.)
    v2.insert(0,1.)
    v3.insert(0,1.)
    test = np.asarray([v1,v2,v3])
    if np.linalg.det(test) ==0:
        return True
    else:
        return False

class Canvas(object):
    def __init__(self):
        self.fig            = plt.figure(1)
        self.ax             = self.fig.add_subplot(111)
        self.fig.add_subplot(self.ax)
        return
       
class Plotter(object):
    
    def __init__(self, title = '', num_plots = 1,
                 x_axis_name = 'y',
                 y_axis_name = 'z'):
        """ self.legend_loc set in function set_range
        """
        self.fig            = None
        self.verbose        = False
        self.range_set      = False
        self.tableau_set    = False
        self.ticks_set      = False
        self.title          = title
        self.num_plots      = num_plots
        self.x_axis_name    = x_axis_name
        self.y_axis_name    = y_axis_name
        plt.rc('text', usetex=True)
        plt.rc('font', family='sans-serif')
        if self.num_plots == 1:        
            self.fig            = plt.figure(1)
            self.ax             = self.fig.add_subplot(111)
            #self.ax = SubplotZero(self.fig, 111)
            self.fig.add_subplot(self.ax)
        elif self.num_plots > 1:
            if self.verbose: print 'warning, using autoplotting with n>1 plots'
            #self.ax := (ax1,ax2,..axn)
            self.fig, self.ax = plt.subplots(1, 2)
        self.tableau_set    = self.set_tableau()
        #self.range_set      = self.set_range()
        #self.ticks_set      = self.set_ticks()
        return
    

    def __call__(self,  composite_curve, ADLsplines, mytitle = 'default title'):
        """After init
        call on the class to ploi a list of ADIALSplines
        """
        self.display(composite_curve, ADLsplines, mytitle)
        return
        
        
    def display(self, composite_curve, ADLsplines, mytitle):

        xmin,xmax,ymin,ymax = composite_curve.extreme_C0()
        self.set_range(xi = xmin, xe = xmax, dx = xmax-xmin,
                       yi = ymin, ye = ymax, dy = ymax-ymin)
        self.set_ticks()
        for adialsp in ADLsplines:
            self.plot_this(plotfunc = adialsp.curve.plot_subplot_curve )
            for el, hue in zip(adialsp.Lagrangian.equality, self.tableau20):
                self.plot_this(
                    plotfunc = adialsp.Lagrangian.equality[el].object.plot,
                    object_ = adialsp.Lagrangian.equality[el],
                    Lspline_ = adialsp,
                    const_num = el)#, 
                    #color_='green')
        self.legend_loc[1] -=1.
        self.plot_composite(curve = composite_curve)
        self.save_pdf(filename=mytitle, ftype = '.png')
        plt.show()
        return
        
    def set_tableau(self):
        """ Setup the "Tableau 20" colors as RGB.    
        """
        self.tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
                          (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
                          (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
                          (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
                          (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)] 
        for i in range(len(self.tableau20)):    
            r, g, b = self.tableau20[i]    
            self.tableau20[i] = (r / 255., g / 255., b / 255.) 
        return True
        
        
    def set_range(self,
                  xi,xe, dx, 
                  yi,ye, dy, 
                  font_size = 14):
        """Limit the range of the plot to only where the data is.
        Avoid unnecessary whitespace.
        """  
        self.xi = int(xi)
        self.xe = int(xe)
        self.dx = int(dx)
        self.yi = int(yi)
        self.ye = int(ye)
        self.dy = int(dy)
        self.min_x = min(self.xi,self.xe)
        self.min_y = min(self.yi,self.ye)
        self.max_x = max(self.xi,self.xe)
        self.max_y = max(self.yi,self.ye)
        #dx = xe-xi
        #dy = ye-yi
        #excessx = .5
        #dummy = xi -excessx*dx
        #other_x_loc = -5.
        self.max_x  += 1.5
        other_x_loc = self.max_x + 1.
        self.legend_loc = [other_x_loc, self.max_y, other_x_loc]
        
        self.font_size = font_size
        #plt.ylim(yi, ye)    #quasi global variables???
        #plt.xlim(xi, xe) 
        #plt.xlim(-12., xe) 
        
        self.range_set = True
        if self.verbose:
            return True
        else:
            return
            
    
    
    def set_ticks(self):
        # remove framelines
        self.ax.spines["top"].set_visible(False)    
        self.ax.spines["bottom"].set_visible(True)    
        self.ax.spines["right"].set_visible(False)    
        self.ax.spines["left"].set_visible(True)  
        
        # ticks on the bottom and left of the plot  only
        self.ax.get_xaxis().tick_bottom()    
        self.ax.get_yaxis().tick_left() 
        max_arrow = max(self.max_x,self.max_y) +1.
        self.ax.arrow(0., 0., max_arrow, 0., head_width=0.05, head_length=0.1, fc='k', ec='k')
        self.ax.arrow(0., 0., 0., max_arrow, head_width=0.05, head_length=0.1, fc='k', ec='k')
        
        return True
        
    def plot_this(self,  
                  plotfunc, 
                  object_ = None,
                  Lspline_=None, 
                  color_ = 'grey',
                  const_num = None,
                  docurvature=True):
        if object_ != None:
            if docurvature == False:
                axi = plotfunc(object_, 
                               Lspline = Lspline_,
                               canvas = self.ax,
                               legendloc = self.legend_loc,
                               cn = const_num,
                               docurvature = False)#, color = color_)
                self.legend_loc[1] -= 1.
            else:
                
                axi = plotfunc(object_, 
                               Lspline = Lspline_,
                               canvas = self.ax,
                               legendloc = self.legend_loc,
                               cn = const_num)#, color = color_)
                self.legend_loc[1] -= 1.
            
        else:
            axi = plotfunc( Lspline = Lspline_,
                           canvas = self.ax)#, color = color_)
        self.fig.add_subplot(axi)
        return
        
    def plot(self, 
             Lspline = None, 
             excessx = .75,
             excessy = .75):
        plt.axis('equal')
        if Lspline != None:
            xi,xe,yi,ye = Lspline.curve.extreme_C0()
            dx = xe-xi
            dy = ye-yi
            #plt.xlim(xi-.5*excessx*dx,xe+excessx*dx)
            #plt.ylim(yi-excessy*dy,ye+excessy*dy)
            name = str(self.x_axis_name)
            plt.annotate(name, xy = (self.max_x,-1.), xycoords='data')
            name = str(self.y_axis_name)
            plt.annotate(name, xy = (-1.,self.max_y), xycoords='data')
        plt.show()
        return
        
    def plot_composite(self,
                       curve,
                       excessx = .75,
                       excessy = .75):
        #plt.axis('equal')
        xi,xe,yi,ye = curve.extreme_C0()
        dx = xe-xi
        dy = ye-yi
        plt.xlim(xi-excessx*dx,xe+excessx*dx)
        plt.ylim(yi-excessy*dy,ye+excessy*dy)
        name = str('x')#self.xe)
        plt.annotate(name, xy = (self.max_x,-1.), xycoords='data')
                        
        name = str('y')#self.ye)
        plt.annotate(name, xy = (-1.,self.max_y), xycoords='data')
        
        plt.show()
        return
    
    def save_pdf(self, filename = None, ftype = '.pdf'):
        """ save pdf file.
        No file extension needed.
        """
        if filename == None:
            filename = default_input('please enter a name for the picture', 'curve')
        plt.savefig(filename+ftype, bbox_inches = 'tight')
        plt.close()
        return
    

class SAC(object):
    k       = 4
    nump    = 30
    xb      = 0.
    yb      = 0.
    def __init__(self, Amax, #v1,v2,v3, 
                 l1,l2,l3, 
                 Xe, Xm, Xr, Xlcb, 
                 Afp, Aap, Cpe, Cpr):
        #self.V      = vol
        self.Amax   = Amax
        #self.Ve     = v1
        #self.Vm     = v2
        #self.Vr     = v3
        self.Le     = l1
        self.Lm     = l2
        self.Lr     = l3
        self.L      = l1+l2+l3
        self.Xe     = Xe
        self.Xm     = Xm
        self.Xr     = Xr
        self.Xlcb   = Xlcb
        self.Afp    = Afp
        self.Aap    = Aap
        self.Cpe    = Cpe
        self.Cpr    = Cpr
        self.go()
        self.aggregate_SAC_curves_Small()
        return
        
    def plot(self):
        """Plot the composite SAC curve
        with all Form Parameters used in
        its construction
        """
        self.nice_plot = Plotter(x_axis_name = 'x',
                                 y_axis_name = 'z')
        ADLsplines = [self.entrance,
                        self.midbody,
                        self.run]
        self.nice_plot(self.SAC,
                       ADLsplines,
                       mytitle = 'SAC')
        
        return
        
    def go(self):
        self.entrance   = self.compute_entrance_SAC()
        self.midbody    = self.compute_midbody_SAC()
        self.run        = self.compute_run_SAC()
        return
        
    def aggregate_SAC_curves_regular(self):
        """put fwd, mid, and run
        into one SAC curve
        """
        n1 = self.run.curve.n
        n2 = self.midbody.curve.n
        n3 = self.entrance.curve.n
        self.midbody.translate(self.Le+SAC.xb)
        self.entrance.translate(self.Lm + self.Le+SAC.xb)
        nv = n1 + n2 + n3
        
        Cv = np.zeros((nv,2),float)
        Cv[0:n1,:] = self.run.curve.vertices
        Cv[n1:n1+n2,:] = self.midbody.curve.vertices
        Cv[n1+n2:n1+n2+n3,:] = self.entrance.curve.vertices
        
        self.SAC = spline.Bspline(Cv, SAC.k, nump = SAC.nump*3)
        return

    def aggregate_SAC_curves_Small(self):
        """put fwd, mid, and run
        into one SAC curve
        """
        n1 = self.run.curve.n
        n2 = self.midbody.curve.n
        n3 = self.entrance.curve.n
        self.midbody.translate(self.Le+SAC.xb)
        self.entrance.translate(self.Lm + self.Le+SAC.xb)
        nv = n1 + n2-2 + n3
        
        Cv = np.zeros((nv,2),float)
        Cv[0:n1,:] = self.run.curve.vertices
        Cv[n1:n1+n2-2,:] = self.midbody.curve.vertices[1:-1]
        Cv[n1+n2-2:n1+n2+n3,:] = self.entrance.curve.vertices
        
        self.SAC = spline.Bspline(Cv, SAC.k, nump = SAC.nump*3)
        return
    
    #def compute_entrance_SAC(self):
    def compute_run_SAC(self):
        k       = SAC.k
        nump    = SAC.nump
        nCV     = 7
        xb      = SAC.xb
        xe      = xb + self.Le
        yb      = SAC.yb
        ye      = self.Amax
        yfp     = self.Afp
        Xc      = self.Xe
        ab = 0.
        ae = 0.
        ini_v = InitializeControlVertices(xb,yb,xe,ye,
                                          alphae=ae, 
                                          alphab=ab,
                                          Cab_given=0.,
                                          Cae_given=0.,
                                          nCV = 7, 
                                          slope = 'up')
        curve = spline.Bspline(ini_v.vertices, k, nump)  
        
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        #FPD.add_AreaConstraint(kind='equality', value = curve_area)
        FPD.add_AngleConstraint(kind = 'equality', value = 0., location = 0.)
        FPD.add_AngleConstraint(kind = 'equality', value = 0., location = 1.)
        FPD.add_CurvatureConstraint(kind = 'equality', value = 0., location = 0.)
        FPD.add_CurvatureConstraint(kind = 'equality', value = 0., location = 1.)
        FPD.add_XcConstraint(kind = 'equality', value = 8.5)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        
        Lspline.compute_lagrangian()
        Lspline.display(mytitle = 'SAC_run_ini', 
                        x_axis_name = 'x',
                        y_axis_name = 'z')
        Lspline.optimize()
        Lspline.display(mytitle = 'SAC_run', 
                        x_axis_name = 'x',
                        y_axis_name = 'z')
        return Lspline
        
        
    def compute_midbody_SAC(self):
        k       = SAC.k
        nump    = SAC.nump
        n       = 7
        xb = SAC.xb
        xe = xb + self.Lm
        yb = self.Amax
        ye = self.Amax
        Xc = self.Xm
        ab = 0.
        ae = 0.
        vertices = linear_vertices([xb,yb],
                                   [xe,ye],
                                   num=n)
        curve = spline.Bspline(vertices, k, nump)  
        
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        #FPD.add_AreaConstraint(kind='equality', value = curve_area)
        FPD.add_AngleConstraint(kind = 'equality', value = 0., location = 0.)
        FPD.add_AngleConstraint(kind = 'equality', value = 0., location = 1.)
        FPD.add_CurvatureConstraint(kind = 'equality', value = 0., location = 0.)
        FPD.add_CurvatureConstraint(kind = 'equality', value = 0., location = 1.)
        #FPD.add_XcConstraint(kind = 'equality', value = 6.)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        
        Lspline.compute_lagrangian()
        Lspline.display(mytitle = 'SAC_mid_ini', 
                        x_axis_name = 'x',
                        y_axis_name = 'z')
        Lspline.optimize()
        Lspline.display(mytitle = 'SAC_mid', 
                        x_axis_name = 'x',
                        y_axis_name = 'z')
        return  Lspline
        
    
    #def compute_run_SAC(self):
    def compute_entrance_SAC(self):
        k       = SAC.k
        nump    = SAC.nump
        nCV     = 7
        xb  = SAC.xb
        xe  = xb + self.Lr
        yb  = self.Amax
        ye  = 0.
        yap = self.Aap
        Xc  = self.Xr
        ab = 0.
        ae = 0.
        ini_v = InitializeControlVertices(xb,yb,xe,ye,
                                          alphae=ae, 
                                          alphab=ab,
                                          Cab_given=0.,
                                          Cae_given=0.,
                                          nCV = 7, 
                                          slope = 'down')
        curve = spline.Bspline(ini_v.vertices, k, nump)  
        
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        #FPD.add_AreaConstraint(kind='equality', value = curve_area)
        FPD.add_AngleConstraint(kind = 'equality', value = 0., location = 0.)
        FPD.add_AngleConstraint(kind = 'equality', value = 0., location = 1.)
        FPD.add_CurvatureConstraint(kind = 'equality', value = 0., location = 0.)
        FPD.add_CurvatureConstraint(kind = 'equality', value = 0., location = 1.)
        FPD.add_XcConstraint(kind = 'equality', value = 3.5)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        
        Lspline.compute_lagrangian()
        Lspline.display(mytitle = 'SAC_fwd_ini', 
                        x_axis_name = 'x',
                        y_axis_name = 'z')
        Lspline.optimize()
        Lspline.display(mytitle = 'SAC_fwd', 
                        x_axis_name = 'x',
                        y_axis_name = 'z')
        return  Lspline
        
        
def initial_curve(start=(0.,12.), end=(12.,0.), num=7, k=4,nump=50):
    vertices = linear_vertices(start, end, num)
    k=4
    nump=50
    curve = spline.Bspline(vertices, k, nump)
    return curve        
    
def linear_vertices(start, end, num):
    start = np.asarray(start)
    end = np.asarray(end)
    D = end - start
    dim = len(D)
    vertices = []
    for i in range(dim):
        vi = np.linspace(start[i],end[i],num)
        vertices.append(vi)
    vertices = np.asarray(vertices)
    vertices = np.transpose(vertices)
    return vertices        


DIM = 3        
def read_maxsurf(filename = 'simple_surf.txt'):
        maxsurf_file  = open(filename)
        lines =  maxsurf_file.readlines()
        maxsurf_file.close()
        data = []
        surfaces = {}
        #surface_name = 'surf'
        for line in lines:
            data.append(line.split())
        for el in data:
            assert(isinstance(el[0],str))
            try:
                surface_name = el[0]
                cn = int(el[1])
            except:
                surface_name = el[0]+el[1]
                el.pop(1)
                cn = int(el[1])
            
            if surface_name in surfaces:
                if cn in surfaces[surface_name]:
                    surfaces[surface_name][cn].append([ float(a) for a in el[3:3+DIM]])
                else:
                    surfaces[surface_name][cn]=[]
                    surfaces[surface_name][cn].append([ float(a) for a in el[3:3+DIM]])
            else:
                surfaces[surface_name] = {}
                if cn in surfaces[surface_name]:
                    surfaces[surface_name][cn].append([ float(a) for a in el[3:3+DIM]])
                else:
                    surfaces[surface_name][cn]=[]
                    surfaces[surface_name][cn].append([ float(a) for a in el[3:3+DIM]])
                
        for s in surfaces:
            for key in surfaces[s]:
                surfaces[s][key] = np.asarray(surfaces[s][key])
            
            
        return surfaces

def make_tnet(surfaces, k=4, nump=30):
    for s in surfaces:
        nl = len(surfaces[s].keys())
        nt = len(surfaces[s][0])
        tnet = []
        
        tarray = np.zeros((nt,nl,DIM),float)
        for i in range(nt):
            for j in range(nl):
                #print i,j
                tarray[i,j] = surfaces[s][j][i]
        if k<nl:
            for i in range(nt):
                tnet.append(spline.Bspline(tarray[i], k, nump))
            surfaces[s]['tcurve_net'] = tnet #TLM backspace error
        surfaces[s]['tnet'] = tarray
    return surfaces
 
        

def make_maxsurf_long_curves(surfaces, k=4, nump = 30):
    for s in surfaces:
        lcurve_net = []
        for key in surfaces[s]:
            try:
                lcurve_net.append(spline.Bspline(surfaces[s][key],k, nump))
            except:
                print 'could not make curves for {}'.format(s)
        surfaces[s]['lcurve_net'] = lcurve_net
    return surfaces
    

def plot_mx(tcurve_net):
    """plots a single max surface tcurve_net
    """
    #surfaces['Stern']['tcurve_net'][3].plot()
    tcurve_net[0].plot3DmultiList(tcurve_net,[])
    return
    
def plot_surfaces(surfaces):
    t_list = []
    l_list = []
    for s in surfaces:
        if 'tcurve_net' in surfaces[s]:
            for curve in surfaces[s]['tcurve_net']:
                t_list.append(curve)
        if 'lcurve_net' in surfaces[s]:
            for curve in surfaces[s]['lcurve_net']:
                l_list.append(curve)
    curve.plot3DmultiList(t_list,l_list, unlimited=True, limx=80., limy=80., limz = 80.)
    return

def swap_indices(curve,i,j):
    if isinstance(curve, spline.Bspline):
        new_vertices = copy.deepcopy(curve.vertices)
        new_vertices[:,i] = curve.vertices[:,j]
        new_vertices[:,j] = curve.vertices[:,i]
    elif isinstance(curve, np.ndarray):
        new_vertices = copy.deepcopy(curve)
        new_vertices[:,i] = curve[:,j]
        new_vertices[:,j] = curve[:,i]
    return new_vertices
    
def reverse_curve(curve):
    if isinstance(curve, spline.Bspline):
        new_vertices = []
        for vertex in reversed(curve.vertices):
            new_vertices.append(vertex)
        new_vertices = np.asarray(new_vertices)
        new_curve = spline.Bspline(new_vertices, curve.k, curve.nump)
    return new_curve

def flip_sign(vertices,i):
    vertices[:,i] = -vertices[:,i]
    return vertices
    

def define_hull():
    compute_panels      = False
    normalize           = True
    compute_volumes     = True
    compute_volKernels  = False
    interpolated_x      = False
    interpolated_y      = False
    
    # number of hull panels and points (small numbers to test program)
    nhl = 32  # lengthwise 
    nhg = 16   # girthwise below WL  ex: for Fn 0.3 wavelength hull of length 1, I need # panels
    nhf = 1   # girthwise above WL
    
    hFraction=.8
    xFraction=.65
    Compute_Curvatures=False
    L=100
    return
    
class surface(object):
    def __init__(self, surface):
        self.surface = surface
        self.tnet = surface['tnet']
        self.lnet = surface['lcurve_net']
        self.interpolated_x = False
        self.interpolated_y = False
        self.nhl = 32  # lengthwise 
        self.nhg = 16   # girthwise below WL  ex: for Fn 0.3 wavelength hull of length 1, I need # panels
        self.nhf = 1   # girthwise above WL
        self.hFraction=.8
        self.xFraction=.65
        self.Compute_Curvatures=False
        return
        
    def old_surf(self):
        """We must condition the geometry:
        move the hull bow to the zero point
        scale the hull by its max length
        this is how the nonlinear panel method
        is formated
        """
        #Points networks:
        tnet    = [] #Transverse Curve Net:
        lnet    = [] #Longitudinal Curve Net
        tnet = self.tnet
        L=0.
        lmin = 0.
        B=0.
        T=0.
        #normalize the ship by submerged max length
        for el in tnet:
            for a in el:
                L=max(L,abs(a[0]))
                lmin = min(lmin,abs(a[0]))
                B=max(B,abs(a[1]))
                T=max(T,abs(a[2]))
        L=abs(L)
        lmin = abs(lmin)
        L = L - lmin
        B=abs(B)
        T=abs(T)
        # Move the ship to have the foremost submerged part be at x=0
        # and extend in the potitive x direction to the stern.
        # maybe counterintuitive geometrically...
        # but sets up a flow in the +x direction!
        for el in tnet:
            for a in el:
                a[0]=a[0]+L
        # rid ourselves of that pesky 'about zero' flip flop error
        # which actually causes us to underestimate the
        # length of the ship
        """
        for el in tnet:
            for a in el:
                L=max(L,abs(a[0]))
        for el in tnet:
            for a in el:
                a[0]=a[0]/L
                a[1]=a[1]/L
                a[2]=a[2]/L
        B=B/L
        T=T/L
        L=L/L
        #"""
        
        
        ## --- --- --- Optimize Transverse Curves Here --- --- --- ##
        ##
        
        ## call routines...
        
        
        
        ## --- --- --- Make Transverse curve Net --- --- --- ##
        ##
        tCurveNet = []
        
        if self.interpolated_x == True:
            for el in tnet:
                tCurveNet.append(spline.interpolatedBspline(el, k, nump))
                
        elif self.interpolated_x == False:
            for el in tnet:
                tCurveNet.append(spline.Bspline(el, k, nump))
        
        
        ## --- --- --- Optimize Longi Curves Here --- --- --- ##
        ##
        
        ## call routines...
        
            
        
        ## --- --- --- Make Longitudinal curve Net --- --- --- ##
        ##
        #"""
        tray = np.asarray(tnet)
        for i in range(len(tray[0])):
            print i
            lnet.append(tray[:,i])
        
        
        lCurveNet = []
        
        if self.interpolated_y == True:
            for ii,vertices in enumerate(lnet):
                tempCurve = spline.interpolatedBspline(vertices, k, nump)
        
                lCurveNet.append(tempCurve)
        elif self.interpolated_y == False:
            for ii,vertices in enumerate(lnet):
                tempCurve = spline.Bspline(vertices, k, nump)
                lCurveNet.append(tempCurve)
                
        surf1 = makeSurf(tnet,4,4,nump)
        #if compute_panels ==True:
        Lfs = 4.*L   #free surface discretization length
        nfsl = 4*self.nhl
        dx = Lfs/float(nfsl)
        deltax = self.xFraction*dx
        panelShift  = deltax
        
        print '------------------------------------------------------------------------------------------------'
        #nhl = int(raw_input('hull  panels  lengthwise --> '))
        #nhg = int(raw_input('hull  panels  girthwise --> '))
        #nhf = int(raw_input('# panels girthwise ABL --> '))
        
        print 'Set new free surface panel height fraction:'
        print 'to position of free surface panels above the calm water level, z=0.'
        print 'H. Raven suggests: about 0.8 of panel length zfspanels = 0.8 * L/ nhl \n TLM suggests .65 with his code'
        
        print 'Current setting : {}'.format(self.hFraction)
        hFraction = .65 #float(raw_input('New Setting --> '))
        FSH         = hFraction*L / float(self.nhl)
        
        print '\n\n'
        
        
        x = Lfs/float(nfsl)
        print 'Set the position of the collocation points with a x point shift fraction \n.'
        print 'Collocation points on the free surface must be in front of the parent panel to ensure waves stay behind the ship. \n TLM suggests .8 with his code'
        print 'current setting : {}'.format(self.xFraction)
        self.xFraction = 0.8 #float(raw_input('New Setting --> '))
        deltax = self.xFraction*dx
        panelShift  = deltax
        
        
        points, panels,Bfs,yfs,npoints,npanels,nfspanels,ntransompanels,ntransomt,ntransoml,tnpt =\
        polymake(surf1, FSH, panelShift, self.nhl, self.nhg, self.nhf, self.Compute_Curvatures, L,B,T)

        #points, panels = polymake(surf1, FSH, panelShift, self.nhl, self.nhg, self.nhf, self.Compute_Curvatures, L,B,T)
        
    
    
        ##*******************************************************************************
        #if Compute_Curvatures==True:
        
        return surf1
        
def makeSurf(allVertices,kT,kL,nump):

            
    ## --- --- --- Make Transverse curve Net --- --- --- ##
    ##
    tCurveNet = []
    for el in allVertices:
        tCurveNet.append(spline.Bspline(el, kT, nump))


    ## --- --- --- Make Longitudinal curve Net --- --- --- ##
    ##
    lnet    = []
    tray = np.asarray(allVertices)
    for i in range(len(tray[0])):
        lnet.append(tray[:,i])
        
    lCurveNet = []
    for ii,vertices in enumerate(lnet):
        tempCurve = spline.Bspline(vertices, kL, nump)

        lCurveNet.append(tempCurve)

    ## --- --- --- Make the Surface --- --- --- ##
    ##
    surf = spline.BsplineSurface(tCurveNet,lCurveNet)
    
    return surf
    
    
    
    
class gfunc(object):
    """Match a function
    between param_strt and param_end
    
        storage:
                stored in the Lagrangian 
                as data in a dict:
                sp[name] -> function(AD_vertices)
        usage:
                for el in self.Lagrangian.sp:
                    f += self.Lagrangian.sp[el](vertices)
                    
    """
    def __init__(self, *args):
        self.curve          = args[0] #curve to be designed
        self.func           = args[1] #target function
        self.param_strt     = args[2]
        self.param_end      = args[3]
        self.trange         = args[4] #target function parmeterization
        self.weight         = args[5]
        self.nump           = len(self.trange)
        #assert( len(self.ptarget ) == self.nump )
        self.k              = self.curve.k
        self.p              =self.curve.p
        self.span           = np.zeros((self.nump),int)
        self.prange         = np.linspace(self.param_strt, self.param_end, self.nump)
        self.gf_basis       = np.zeros((self.nump,self.curve.n),float)
        self.target         = []
        self.analysis_type = 'LS'
        for i, s, targ in zip(range(self.nump), self.prange, self.trange):
            span = self.curve.FindSpan(s)
            self.span[i] = span
            self.curve.BasisFuns(span,s,
                                 self.gf_basis[i,span-self.p:span+1] )
            self.target.append(self.func(targ))
        return
        
    def __call__(self, *args):
        """
        func = b.CurvePoint
        range = np.linspace(0.,.4139,20)
        """
        verts   = args[0]
        xpts    = verts[0]
        ypts    = verts[1]
        #func    = args[0]
        #prange  = args[1]
        arc_e = 0.
        #Xcenter = self.cp[0]
        #Ycenter = self.cp[1]
        
            

        basiso = self.gf_basis[0]
        #targo = self.target[0]
        qxo = np.dot(xpts,basiso)
        qyo = np.dot(ypts,basiso)
      
        for i,basis in enumerate(self.gf_basis):
            #print i, basis
            #targ = self.target[i]
            qx = np.dot(xpts,basis)
            qy = np.dot(ypts,basis)
            #arc_e +=  (qx-targ[0])**2 + (qy-targ[1])**2 
            
            dx = qx.value - qxo.value
            dy = qy.value - qyo.value
            dr = (dx*dx + dy*dy)**.5
            targ = self.target[i]
            test = self.weight*dr*( (qx-targ[0])**2 + (qy-targ[1])**2 )
            
            if self.analysis_type == 'min' and test.value > 0:
                arc_e +=  self.weight*dr*( (qx-targ[0])**2 + (qy-targ[1])**2 )#**.5)
                #arc_e +=  (np.sqrt((qx-Xcenter)**2 + (qy-Ycenter)**2) - radius)**2.
            elif self.analysis_type =='LS':
                arc_e +=  self.weight*dr*( (qx-targ[0])**2 + (qy-targ[1])**2 )#**.5)
                
        
        self.result = arc_e
        if self.analysis_type == 'LS':  
            return self.weight*arc_e 
        elif self.analysis_type == 'min':
            #if arc_e <0.:
            #    return arc_e*0. 
            #else:
            return 0.*arc_e  #self.weight* 
                
            

class circle(object):
    def __init__(self, origin, radius):
        self.origin = origin
        self.radius = radius
        return
    def circle_pt(self, s):
        return self.origin[0] + self.radius*np.cos(s), self.origin[1] + self.radius*np.sin(s)
    def plot(self, prange):
        c = np.asarray(self.circle_pt(prange)).T
        plt.plot(c[:,0],c[:,1])
        plt.show()
        return

class ellipse(object):
    def __init__(self, origin, ra,rb):
        self.origin = origin
        self.ra = ra
        self.rb = rb
        return
    def circle_pt(self, s):
        return self.origin[0] + self.ra*np.cos(s),self.origin[1] + self.rb*np.sin(s)
    def plot(self, prange):
        c = np.asarray(self.circle_pt(prange)).T
        plt.plot(c[:,0],c[:,1])
        plt.show()
        return
    
    
class BBowSAC(object):
    def __init__(self, 
                 Xb, Xe, height, 
                 Ccapr = .05,
                 Np=6):
        self.start  = [Xb,0.]
        self.end    = [Xe,height]
        self.height = height
        self.length = abs(Xe-Xb)
        self.Ccapr  = Ccapr #Coeff : cap to total length ratio
        self.h      = self.Ccapr*self.length
        radius      = self.height*.5
        self.Lflat  = self.length*(1.-Ccapr) #length of the barrel-like 
                                            #portion of bbow
        self.rLflat = self.length*(Ccapr) #inverse that.
        #volume of the bulbous bow:
        self.area   = self.height*self.length*(1.-Ccapr) + \
                        (.333333)*np.pi*(self.h**2)*\
                        (3.*radius - self.h )
        self.k      = 4
        self.nump   = 30
        self.Np     = Np
        self.Lspline = self.make_SAC()
        self.SAC = self.Lspline.curve
        return
    
    def __call__(self, u):
        return self.SAC(u)
    
    def make_SAC(self):
        verts = linear_vertices(self.start,
                                self.end,self.Np)
        curve = spline.Bspline(vertices = verts,
                               k=self.k,
                               nump=self.nump)
        #
        #******************************************************************
        # 
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve)
        #
        #******************************************************************
        # 
        FPD.add_AreaConstraint(kind = 'equality',
                               value = self.area)
        FPD.add_yFixity(index=-1,
                        value=self.height)
        FPD.add_yFixity(index=-2,
                        value=self.height)
        FPD.add_yFixity(index=-3,
                        value=self.height)
        FPD.add_yFixity(index=3,
                        value=self.height)
        FPD.add_xFixity(index=-2,
                        value=self.rLflat)
        FPD.add_xFixity(index=-3,
                        value=self.rLflat)
        #
        FPD.add_E1(kind='LS', weight = .1)
        FPD.add_E2(kind='LS', weight = .1)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        #
        #******************************************************************
        # 
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data,
                                                     small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        #
        #******************************************************************
        # Solve
        Lspline.optimize()
        Lspline.curve.compute_curve()
        return Lspline
    
    def plot(self, prange):
        c = np.asarray(self.circle_pt(prange)).T
        plt.plot(c[:,0],c[:,1])
        plt.show()
        return
    

def translate(*args):
    curve = args[0]
    if len(args)==3:
        translate = [args[1],args[2]]
    else:
        translate = args[1]
    v1 = copy.deepcopy(curve.vertices)
    v1 = v1 + translate
    curve.vertices = v1
    return curve
    
class frame(object):
    def __init__(self):
        return
    @classmethod
    def add1D(self, array, which_index=2):
        size,dim = np.shape(array) 
        new_array = np.zeros((size,dim+1),float)
        if which_index == 0:
            new_array[:,1:] = array[:]
        if which_index == 1:
            new_array[:,0] = array[:,0]
            new_array[:,2] = array[:,1]
        if which_index == 2:
            new_array[:,0] = array[:,0]
            new_array[:,1] = array[:,1]
        return new_array

if __name__ == "__main__":
    # TLM B-Spline Curve class
    import curve             as     spline
    from   initialValues     import InitializeControlPoints, InitializeControlVertices
    import copy
    from   ADILS             import IntervalLagrangeSpline, Lagrangian
    from   FormParameter     import FormParameterDict
    from   initialValues     import interval_bounds, lagrangian_bounds


    A_mid = 10.080798436245335
    BulbLength = 4.341062137105493
    self = BBowSAC(Xb=0.,Xe = BulbLength,
                               height=A_mid,
                               Ccapr=.1)
    self.SAC.plotcurve_detailed()
    
    option  = 'other'
    #test    = None
    #test    = 'container'
    #test    = 'destroyer'
    #test    = 'destroyer_bow_pic'
    test    = 'destroyer_bow_seq_area'
    test    = None
    do_all = False
        
    
        
    if do_all:
        test = 'destroyer'
    if test == 'destroyer':
        filename = 'destroyer.dat'
        #"""
        #-----Curve 2
        nump = 30
        k = 4
        
        surfaces = read_maxsurf(filename)
        s= 'HullSide'
        self = surfaces
        surfaces = make_tnet(surfaces)
        surfaces = make_maxsurf_long_curves(surfaces)
        #plot_surfaces(surfaces)
        
        surfs = []
        for s in surfaces:
            try:
                si = spline.BsplineSurface(surfaces[s]['tcurve_net'],surfaces[s]['lcurve_net'])
                surfs.append(si)
                surfaces[s]['sf'] = si
            except:
                pass
         
        ## Make 1 Hull from 2 surfaces...
        #
        #---------------------------------------------------------------------
        #
        """
        s = 'Bow'
        n = 10
        k_=4
        nump_=30
        match_bow = []
        itcv = []
        count_list = []
        #for i in range(len(s1.tcurveNet)): #19
        for i in range(len(surfaces['Bow']['tcurve_net'])): #19
            #xv  = s1.tcurveNet[i].vertices[:,0]
            #bv  = s1.tcurveNet[i].vertices[:,1:3] 
            
            xv = surfaces['Bow']['tcurve_net'][i].vertices[:,0]
            bv = surfaces['Bow']['tcurve_net'][i].vertices[:,1:3] 
            
            
            b   = spline.Bspline(vertices = bv, k=k_, nump=nump_)
            curve = initial_curve(bv[0],bv[-1],n)
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve) 
            for s in np.linspace(0.01,.99,10):
                value_ = b.CurvePoint(s)
                FPD.add_xPointConstraint(kind='LS', location=s, value=value_[0], weight=10000.)#*1./(s))
                FPD.add_yPointConstraint(kind='LS', location=s, value=value_[1], weight=10000.)#*1./(s))
            
            FPD.add_E1(kind='LS', weight = .01)
            FPD.add_E2(kind='LS', weight = .01)
            FPD.add_E3(kind='LS', weight = .01)
            FPD.add_ArcLengthApprox(kind='LS', weight = 1.0)
            nump = 50
            mp = .4
            #gfc = gfunc(curve, b.CurvePoint, 0., mp, np.linspace(0., mp, nump), 1000000. )
            L = Lagrangian(FPD)
            #L.sp['crv'] = gfc
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            #Lspline.curve.plot()
            Lspline.optimize(stop = 50)
            itcv.append(Lspline.curve)
            count_list.append([Lspline.nitr,Lspline.curve])
        
        bbow3d=[]
        for i in range(len(itcv)):
            nv = frame.add1D(itcv[i].vertices, 0)
            #ax = np.average(s1.tcurveNet[i].vertices[:,0])
            #nv[:,0] = ax
            #for i in range(n):
            #curve = initial_curve(s1.tcurveNet[i].vertices[0,0:2],s1.tcurveNet[i].vertices[-1,0:2],n) 
            curve = initial_curve(surfaces['Bow']['tcurve_net'][i].vertices[0,0:2],
                                  surfaces['Bow']['tcurve_net'][i].vertices[-1,0:2], n)
            nv[:,0] = curve.vertices[:,0]
            bbow3d.append(spline.Bspline(vertices = nv, k=k_, nump=nump_))
            
        newhull = surfaces['Mainhull']['tcurve_net']
        for el in bbow3d:
            newhull.append(el)
        #
        #---------------------------------------------------------------------
        #
        tarray = np.zeros((len(newhull),newhull[0].n,DIM),float)
        surfaces['newull'] = {}
        surfaces['newull']['tcurve_net'] = newhull
        surfaces['newull']['lcurve_net'] = []
        surfaces['newull']['tnet'] = []
        for i, curve in enumerate(surfaces['newull']['tcurve_net']):
            surfaces['newull']['tnet'].append(curve.vertices)
            #tarray[i,:] = curve.vertices
        #surfaces['newull']['tnet'] = tarray
        dh = surface(surfaces['Mainhull'])
        self = dh
        db = surface(surfaces['Bow'])
        dt = surface(surfaces['Transom'])
        s1 = db.old_surf()
        
        dnew = surface(surfaces['newull'])
        dt = dnew.old_surf()
        
        
        ##
        ## Make 1 Hull from 2 surfaces...
        ##
        curve.plot3DmultiList(surfaces['newull']['tcurve_net'],[])
        dt.plotSurface()
        curve.plot3DmultiList(dt.tcurveNet,[])
        #"""
        #"""
        k_=4
        nump_ = 30
        db = surface(surfaces['Bow'])
        s1 = db.old_surf()
        f=s1.tcurveNet[16]
        ##
        ## Bow Transverse with Bulb
        ##
        a = spline.Bspline(vertices = s1.tcurveNet[16].vertices, k=4, nump=30)
        
        bv  = 30*s1.tcurveNet[16].vertices[:,1:3] + [0.,2.0934466019417477]
        b   = spline.Bspline(vertices = bv, k=4, nump=30)
        #b.plotcurve_detailed(curvature = 'yes', scale = 0.001 )
          
        n = 8
        curve = initial_curve(bv[0],bv[-1],n) 
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        for s in np.linspace(0.01,1.,10):
            value_ = b.CurvePoint(s)
            FPD.add_xPointConstraint(kind='LS', location=s, value=value_[0], weight=10000.*1./s)
            FPD.add_yPointConstraint(kind='LS', location=s, value=value_[1], weight=10000.*1./s)
        FPD.add_E1(kind='LS', weight = .1)
        FPD.add_E2(kind='LS', weight = .1)
        FPD.add_E3(kind='LS', weight = .1)
        FPD.add_ArcLengthApprox(kind='LS', weight = 100.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
        curve1 = Lspline.curve
        #---------------------------------------------------
        # DONE - Now Plot:
        #
        bp = copy.deepcopy(b)
        bp.vertices = [5.,0] + bp.vertices #move a curve to the right for pic
        #curve1.plotcurve_detailed()
        #bp.plotcurve_detailed()
        ###
        curve2 = copy.deepcopy(curve1)
        Lspline.translate(0.,-Lspline.curve.vertices[0,1])
        curve21 = copy.deepcopy(Lspline.curve)
        #"""
        
        #"""
        #--------------------------------------------------------
        # usegfunc to match just the bulb
        #
        b   = spline.Bspline(vertices = bv, k=k_, nump=nump_)
        curve = initial_curve(bv[0],bv[-1],n)
        #"""
        #tane = b.compute_tangent(1.)
        #curve2.knot_insertion(.58)
        interval_data, small = interval_bounds(curve2)
        FPD = FormParameterDict(curve2) 
        #FPD.add_xVertexConstraint(kind= 'equality', value = b.vertices[-4,0] , index=7)
        #FPD.add_yVertexConstraint(kind= 'LS', value = 3.0 , index=4)
        #FPD.add_AngleConstraint(kind='equality', location = 1., value = -tane)
        for s in np.linspace(0.01,1.,10):
            value_ = b.CurvePoint(s)
            FPD.add_xPointConstraint(kind='LS', location=s, value=value_[0], weight=10000.)#*1./(s))
            FPD.add_yPointConstraint(kind='LS', location=s, value=value_[1], weight=10000.)#*1./(s))
        FPD.add_E1(kind='LS', weight = .01)
        FPD.add_E2(kind='LS', weight = .01)
        FPD.add_E3(kind='LS', weight = .01)
        FPD.add_ArcLengthApprox(kind='LS', weight = .01)
        nump = 50
        mp = .4
        gfc = gfunc(curve2, b.CurvePoint, 0., mp, np.linspace(0., mp, nump), 1000000. )
        L = Lagrangian(FPD)
        L.sp['crv'] = gfc
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve2, L, data = interval_data)
        #Lspline.curve.plot()
        Lspline.optimize(stop = 20) #50
        curve3 = copy.deepcopy(Lspline.curve)
        match_bulb = copy.deepcopy(Lspline.curve) 
        #curve3.plotcurve_detailed()
        #"""
        #-------------------------------------------------------
        # use gfunc to match a circle at the bulb
        #
        #"""
        bbow = []
        radius = 1.33
        rb_ = 1.0
        #for ra_ in [.5,.8,1.,.8,.4]:
        for ra_ in [1.]:
            print ra_
            curve3 = copy.deepcopy(curve21)
            #curve3.knot_insertion(.15)
            #curve3.knot_insertion(.31)
            #curve3.knot_insertion(.45)
            interval_data, small = interval_bounds(curve3)
            FPD = FormParameterDict(curve3) 
            for s in np.linspace(.4,1.,10):
                value_ = b.CurvePoint(s)
                FPD.add_xPointConstraint(kind='LS', location=s, value=value_[0], weight=10000.)#*1./(s))
                FPD.add_yPointConstraint(kind='LS', location=s, value=value_[1], weight=10000.)#*1./(s))
            FPD.add_E1(kind='LS', weight = .1)
            FPD.add_E2(kind='LS', weight = .1)
            FPD.add_E3(kind='LS', weight = .1)
            FPD.add_ArcLengthApprox(kind='LS', weight = .1)
            L = Lagrangian(FPD)
            nump = 50
            mp = .4
            c1 = circle(origin=[0.,1.],radius=1.)
            #c1 = ellipse(origin=[0.,radius],ra=ra_, rb=rb_)
            circ_p = np.linspace(-np.pi/2., np.pi/2., nump)
            gfc3 = gfunc(curve3, c1.circle_pt, 0.,mp,circ_p, 10000.)
            L.sp['crv'] = gfc3
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve3, L, data = interval_data)
            Lspline.optimize(stop = 50)
            bbow.append(copy.deepcopy(Lspline.curve))
        bbowi = []
        for curve, x in zip(bbow,[0.,-.1,-.2,-.3,-.4]):
            nv = frame.add1D(curve.vertices, 0)
            nv[:,0] = x
            bbowi.append(spline.Bspline(vertices=nv, k=4, nump=30 ))
            
        curve4 = copy.deepcopy(Lspline.curve)
        curve4.plot()
        
        c1.plot(np.linspace(-np.pi/2.,np.pi/2.,30))
        #"""
    #-----------------------------------------------------------------
    if test == 'container':
        filename = 'container.dat'

        nump = 30
        k = 4
        
        surfaces = read_maxsurf(filename)
        self = surfaces
        surfaces = make_tnet(surfaces)
        surfaces = make_maxsurf_long_curves(surfaces)
        plot_surfaces(surfaces)
        
        surfs = []
        for s in surfaces:
            try:
                si = spline.BsplineSurface(surfaces[s]['tcurve_net'],surfaces[s]['lcurve_net'])
                surfs.append(si)
                surfaces[s]['sf'] = si
            except:
                pass
        
        
    if test == 'destroyer_bow_pic':
        filename = 'destroyer.dat'
        #"""
        #-----Curve 2
        nump = 30
        k = 4
        
        surfaces = read_maxsurf(filename)
        s= 'HullSide'
        self = surfaces
        surfaces = make_tnet(surfaces)
        surfaces = make_maxsurf_long_curves(surfaces)
        #plot_surfaces(surfaces)
        
        surfs = []
        for s in surfaces:
            try:
                si = spline.BsplineSurface(surfaces[s]['tcurve_net'],surfaces[s]['lcurve_net'])
                surfs.append(si)
                surfaces[s]['sf'] = si
            except:
                pass
        #"""
                
        k_=4
        nump_ = 30
        db = surface(surfaces['Bow'])
        s1 = db.old_surf()
        f=s1.tcurveNet[16]
        ##
        ## Bow Transverse with Bulb
        ##
        a = spline.Bspline(vertices = s1.tcurveNet[16].vertices, k=4, nump=30)
        
        bv  = 30*s1.tcurveNet[16].vertices[:,1:3] + [0.,2.0934466019417477]
        b   = spline.Bspline(vertices = bv, k=4, nump=30)
        b   = translate(b, 0., -bv[0,1])
        bv = copy.deepcopy(b.vertices)
        #b.plotcurve_detailed(curvature = 'yes', scale = 0.001 )
          
        n = 8
        curve = initial_curve(bv[0],bv[-1],n) 
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        for s in np.linspace(0.01,1.,10):
            value_ = b.CurvePoint(s)
            FPD.add_xPointConstraint(kind='LS', location=s, 
                                     value=value_[0], weight=10000.*1./s)
            FPD.add_yPointConstraint(kind='LS', location=s, 
                                     value=value_[1], weight=10000.*1./s)
        FPD.add_E1(kind='LS', weight = .1)
        FPD.add_E2(kind='LS', weight = .1)
        FPD.add_E3(kind='LS', weight = .1)
        FPD.add_ArcLengthApprox(kind='LS', weight = 100.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
        curve1 = Lspline.curve
        #---------------------------------------------------
        # DONE - Now Plot:
        #
        bp = copy.deepcopy(b)
        bp.vertices = [5.,0] + bp.vertices #move a curve to the right for pic
        #curve1.plotcurve_detailed()
        #bp.plotcurve_detailed()
        ###
        curve2 = copy.deepcopy(curve1)
        Lspline.translate(0.,-Lspline.curve.vertices[0,1])
        curve21 = copy.deepcopy(Lspline.curve)
        #"""
        
        #"""
        #--------------------------------------------------------
        # usegfunc to match just the bulb
        #
        #b   = spline.Bspline(vertices = bv, k=k_, nump=nump_)
        #curve = initial_curve(bv[0],bv[-1],n)
        #"""
        #tane = b.compute_tangent(1.)
        #curve2.knot_insertion(.58)
        interval_data, small = interval_bounds(curve2)
        FPD = FormParameterDict(curve2) 
        #FPD.add_xVertexConstraint(kind= 'equality', value = b.vertices[-4,0] , index=7)
        #FPD.add_yVertexConstraint(kind= 'LS', value = 3.0 , index=4)
        #FPD.add_AngleConstraint(kind='equality', location = 1., value = -tane)
        for s in np.linspace(0.01,1.,10):
            value_ = b.CurvePoint(s)
            FPD.add_xPointConstraint(kind='LS', location=s, 
                                     value=value_[0], weight=10000.)#*1./(s))
            FPD.add_yPointConstraint(kind='LS', location=s, 
                                     value=value_[1], weight=10000.)#*1./(s))
        FPD.add_E1(kind='LS', weight = .01)
        FPD.add_E2(kind='LS', weight = .01)
        FPD.add_E3(kind='LS', weight = .01)
        FPD.add_ArcLengthApprox(kind='LS', weight = .01)
        nump = 50
        mp = .4
        gfc = gfunc(curve2, b.CurvePoint, 0., mp, np.linspace(0., mp, nump), 1000000. )
        L = Lagrangian(FPD)
        L.sp['crv'] = gfc
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve2, L, data = interval_data)
        #Lspline.curve.plot()
        Lspline.optimize(stop = 20) #50
        curve3 = copy.deepcopy(Lspline.curve)
        match_bulb = copy.deepcopy(Lspline.curve) 
        #curve3.plotcurve_detailed()
        #"""
        #-------------------------------------------------------
        # use gfunc to match a circle at the bulb
        #
        #"""
        b.compute_area_to_y()
        carea = b.area_to_y
        bbow =[]
        radius = 2.25
        rb_ = 1.8
        #for ra_ in [.5,.8,1.,.8,.4]:
        for ra_ in [2.25]:
            print ra_
            curve3 = copy.deepcopy(curve21)
            #curve3.knot_insertion(.15)
            #curve3.knot_insertion(.31)
            #curve3.knot_insertion(.45)
            interval_data, small = interval_bounds(curve3)
            FPD = FormParameterDict(curve3) 
            for s in np.linspace(.4,1.,10):
                value_ = b.CurvePoint(s)
                FPD.add_xPointConstraint(kind='LS', location=s, value=value_[0], weight=10000.*1./(s))
                FPD.add_yPointConstraint(kind='LS', location=s, value=value_[1], weight=10000.*1./(s))
            #FPD.add_xAttractor(kind='LS', value= 100.,index =4, weight=100.)
            
            #FPD.add_xAttractor(kind='LS', value= 100.,index =5, weight=-100.)
            #FPD.add_xAttractor(kind='equality', value= .35,index =6, weight=-100.)
            FPD.add_xVertexConstraint(kind='min', value=.25, index =3, weight=10000.)
            FPD.add_xVertexConstraint(kind='min', value=.25, index =4, weight=10000.)
            FPD.add_xVertexConstraint(kind='min', value=.25, index =5, weight=10000.)
            FPD.add_xVertexConstraint(kind='min', value=.25, index =6, weight=10000.)
            #FPD.add_yVertexConstraint(kind='LS', value=17., index =6, weight=10000.)
            #FPD.add_xVertexConstraint(kind='LS', value=-25., index =6, weight=10000.)
            FPD.add_Area_to_y_Constraint(kind = 'min', value = carea-1.5)
            #FPD.add_Area_to_y_Constraint(kind = 'equality', value = carea, weight=1.)
            FPD.add_Area_to_y_Constraint(kind = 'max', value = carea+1.5)
            FPD.add_E1(kind='LS', weight = .1)
            FPD.add_E2(kind='LS', weight = .1)
            FPD.add_E3(kind='LS', weight = .1)
            FPD.add_ArcLengthApprox(kind='LS', weight = .1)
            L = Lagrangian(FPD)
            nump = 50
            mp = .4
            gfc = gfunc(curve3, b.CurvePoint, 0.,1., np.linspace(0., mp, nump), 1. )
            nump = 50
            mp = .4
            c1 = circle(origin=[0.,2.25],radius=2.25)
            e1 = ellipse(origin=[0.,rb_],ra=ra_, rb=rb_)
            circ_p = np.linspace(-np.pi/2., np.pi/4., nump)
            gfc3 = gfunc(curve3, e1.circle_pt, 0.,mp,circ_p, 10000.)
            #L.sp['crv']     = gfc
            L.sp['circ']    = gfc3
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve3, L, data = interval_data)
            Lspline.optimize(stop = 50)
            bbow.append(copy.deepcopy(Lspline.curve))
        bbowi = []
        for curve, x in zip(bbow,[0.,-.1,-.2,-.3,-.4]):
            nv = frame.add1D(curve.vertices, 0)
            nv[:,0] = x
            bbowi.append(spline.Bspline(vertices=nv, k=4, nump=30 ))
            
        curve4 = copy.deepcopy(Lspline.curve)
        #curve4.plot()
        
        pspace = np.linspace(-np.pi/2.,np.pi/2.,30)
        c1.plot(pspace)
        e1.plot(pspace)
        #"""        
        
        AL = initial_curve([5.,0.],[5.,b.vertices[-1,1]],7)
        CL = initial_curve([0.,0.],[0.,b.vertices[-1,1]],7)
        BL = initial_curve([-5.,0.],[-5.,b.vertices[-1,1]],7)
        
        ep = ellipse(origin=[4.,rb_],ra=ra_, rb=rb_)
        b.compute_area_to_y()
        curve4.compute_area_to_y()
        b=translate(b, [-5.,0.])
        c = translate(curve4, [5.,0])
        ep.plot(pspace)
        
        b.plotcurve_detailed(color_='black')
        c.plotcurve_detailed(color_='black')
        curve2.plotcurve_detailed(color_='black')
        AL.plot(color_='black')
        BL.plot(color_='black')
        CL.plot(color_='black')
    ##
    ##-------------------------------------------------------------------------
    ##-------------------------------------------------------------------------
    ##-------------------------------------------------------------------------
    ##-------------------------------------------------------------------------
    ##
    if test == 'destroyer_bow_seq_area':
        filename = 'destroyer.dat'
        #"""
        #-----Curve 2
        nump = 30
        k = 4
        
        surfaces = read_maxsurf(filename)
        s= 'HullSide'
        self = surfaces
        surfaces = make_tnet(surfaces)
        surfaces = make_maxsurf_long_curves(surfaces)
        #plot_surfaces(surfaces)
        
        surfs = []
        for s in surfaces:
            try:
                si = spline.BsplineSurface(surfaces[s]['tcurve_net'],surfaces[s]['lcurve_net'])
                surfs.append(si)
                surfaces[s]['sf'] = si
            except:
                pass
        #"""
                
        k_=4
        nump_ = 30
        db = surface(surfaces['Bow'])
        s1 = db.old_surf()
        f=s1.tcurveNet[16]
        ##
        ## Bow Transverse with Bulb
        ##
        a = spline.Bspline(vertices = s1.tcurveNet[16].vertices, k=4, nump=30)
        
        bv  = 30*s1.tcurveNet[16].vertices[:,1:3] + [0.,2.0934466019417477]
        b   = spline.Bspline(vertices = bv, k=4, nump=30)
        b.compute_area_to_y()
        carea = b.area_to_y
        b   = translate(b, 0., -bv[0,1])
        bv = copy.deepcopy(b.vertices)
        #b.plotcurve_detailed(curvature = 'yes', scale = 0.001 )
          
        n = 8
        curve = initial_curve(bv[0],bv[-1],n) 
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        for s in np.linspace(0.01,1.,10):
            value_ = b.CurvePoint(s)
            FPD.add_xPointConstraint(kind='LS', location=s, value=value_[0], weight=10000.*1./s)
            FPD.add_yPointConstraint(kind='LS', location=s, value=value_[1], weight=10000.*1./s)
        FPD.add_Area_to_y_Constraint(kind = 'equality', value = carea)
        FPD.add_E1(kind='LS', weight = .1)
        FPD.add_E2(kind='LS', weight = .1)
        FPD.add_E3(kind='LS', weight = .1)
        FPD.add_ArcLengthApprox(kind='LS', weight = 100.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
        curve1 = Lspline.curve
        #---------------------------------------------------
        # DONE - Now Plot:
        #
        bp = copy.deepcopy(b)
        bp.vertices = [5.,0] + bp.vertices #move a curve to the right for pic
        #curve1.plotcurve_detailed()
        #bp.plotcurve_detailed()
        ###
        curve2 = copy.deepcopy(curve1)
        Lspline.translate(0.,-Lspline.curve.vertices[0,1])
        curve21 = copy.deepcopy(Lspline.curve)
        #"""
        
        #"""
        #--------------------------------------------------------
        # usegfunc to match just the bulb
        #
        #b   = spline.Bspline(vertices = bv, k=k_, nump=nump_)
        #curve = initial_curve(bv[0],bv[-1],n)
        #"""
        #tane = b.compute_tangent(1.)
        #curve2.knot_insertion(.58)
        interval_data, small = interval_bounds(curve2)
        FPD = FormParameterDict(curve2) 
        #FPD.add_xVertexConstraint(kind= 'equality', value = b.vertices[-4,0] , index=7)
        #FPD.add_yVertexConstraint(kind= 'LS', value = 3.0 , index=4)
        #FPD.add_AngleConstraint(kind='equality', location = 1., value = -tane)
        for s in np.linspace(0.01,1.,10):
            value_ = b.CurvePoint(s)
            FPD.add_xPointConstraint(kind='LS', location=s, 
                                     value=value_[0], weight=10000.)#*1./(s))
            FPD.add_yPointConstraint(kind='LS', location=s, 
                                     value=value_[1], weight=10000.)#*1./(s))
        FPD.add_Area_to_y_Constraint(kind = 'equality', value = carea)
        FPD.add_E1(kind='LS', weight = .01)
        FPD.add_E2(kind='LS', weight = .01)
        FPD.add_E3(kind='LS', weight = .01)
        FPD.add_ArcLengthApprox(kind='LS', weight = .01)
        nump = 50
        mp = .4
        gfc = gfunc(curve2, b.CurvePoint, 0., mp, np.linspace(0., mp, nump), 1000000. )
        L = Lagrangian(FPD)
        L.sp['crv'] = gfc
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve2, L, data = interval_data)
        #Lspline.curve.plot()
        Lspline.optimize(stop = 20) #50
        curve3 = copy.deepcopy(Lspline.curve)
        match_bulb = copy.deepcopy(Lspline.curve) 
        #curve3.plotcurve_detailed()
        #"""
        #-------------------------------------------------------
        # use gfunc to match a circle at the bulb
        #
        #"""
        bbow =[]
        radius = 2.25
        rb_ = 1.8
        #for ra_ in [.5,.8,1.,.8,.4]:
        for ra_ in [2.25]:
            print ra_
            curve3 = copy.deepcopy(curve21)
            #curve3.knot_insertion(.15)
            #curve3.knot_insertion(.31)
            #curve3.knot_insertion(.45)
            interval_data, small = interval_bounds(curve3)
            FPD = FormParameterDict(curve3) 
            for s in np.linspace(.4,1.,10):
                value_ = b.CurvePoint(s)
                FPD.add_xPointConstraint(kind='LS', location=s, value=value_[0], weight=10000.*1./(s))
                FPD.add_yPointConstraint(kind='LS', location=s, value=value_[1], weight=10000.*1./(s))
            FPD.add_Area_to_y_Constraint(kind = 'equality', value = carea)
            FPD.add_E1(kind='LS', weight = .1)
            FPD.add_E2(kind='LS', weight = .1)
            FPD.add_E3(kind='LS', weight = .1)
            FPD.add_ArcLengthApprox(kind='LS', weight = .1)
            L = Lagrangian(FPD)
            nump = 50
            mp = .4
            gfc = gfunc(curve3, b.CurvePoint, 0.,1., np.linspace(0., mp, nump), 1. )
            nump = 50
            mp = .4
            c1 = circle(origin=[0.,2.25],radius=2.25)
            e1 = ellipse(origin=[0.,rb_],ra=ra_, rb=rb_)
            circ_p = np.linspace(-np.pi/2., np.pi/4., nump)
            gfc3 = gfunc(curve3, e1.circle_pt, 0.,mp,circ_p, 10000.)
            #L.sp['crv']     = gfc
            L.sp['circ']    = gfc3
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve3, L, data = interval_data)
            Lspline.optimize(stop = 50)
            bbow.append(copy.deepcopy(Lspline.curve))
        bbowi = []
        for curve, x in zip(bbow,[0.,-.1,-.2,-.3,-.4]):
            nv = frame.add1D(curve.vertices, 0)
            nv[:,0] = x
            bbowi.append(spline.Bspline(vertices=nv, k=4, nump=30 ))
            
        curve4 = copy.deepcopy(Lspline.curve)
        #curve4.plot()
        
        pspace = np.linspace(-np.pi/2.,np.pi/2.,30)
        c1.plot(pspace)
        e1.plot(pspace)
        #"""        
        
        ep = ellipse(origin=[4.,rb_],ra=ra_, rb=rb_)
        b=translate(b, [-5.,0.])
        c = translate(curve4, [5.,0])
        ep.plot(pspace)
        bdry = initial_curve([5.,0.],[5.,b.vertices[-1,1]],7)
        
        b.plotcurve_detailed()
        curve2.plotcurve_detailed()
        c.plotcurve_detailed()
        bdry.plot()
        
        
        
        
        
        
        
        bbow =[]
        radius = 2.25
        rb_ = 1.8
        #for ra_ in [.5,.8,1.,.8,.4]:
        ra_ = 2.25
        print ra_
        curve5 = copy.deepcopy(curve21)
        #curve3.knot_insertion(.15)
        #curve3.knot_insertion(.31)
        #curve3.knot_insertion(.45)
        interval_data, small = interval_bounds(curve5)
        FPD = FormParameterDict(curve5) 
        for s in np.linspace(0.1,.99,10):
            value_ = b.CurvePoint(s)
            FPD.add_xPointConstraint(kind='LS', location=s, value=value_[0], weight=10000.*1./(s))
            FPD.add_yPointConstraint(kind='LS', location=s, value=value_[1], weight=10000.*1./(s))
        #FPD.add_xVertexConstraint(kind='min', value= 0.1,index =4, weight=1.)
        #FPD.add_xVertexConstraint(kind='min', value= 0.1,index =3, weight=1.)
        #FPD.add_xVertexConstraint(kind='min', value= 0.1,index =5, weight=1.)
        #FPD.add_xVertexConstraint(kind='max', value= 2.1,index =5, weight=1.)
        
        #FPD.add_xAttractor(kind='LS', value= 100.,index =4, weight=100.)
        #FPD.add_xAttractor(kind='LS', value= 100.,index =5, weight=-100.)
        #FPD.add_xAttractor(kind='equality', value= .35,index =6, weight=-100.)
        #FPD.add_xVertexConstraint(kind='equality', value=.25, index =5, weight=10000.)
        #FPD.add_yVertexConstraint(kind='LS', value=17., index =6, weight=10000.)
        #FPD.add_xVertexConstraint(kind='LS', value=-25., index =6, weight=10000.)
        FPD.add_Area_to_y_Constraint(kind = 'min', value = carea-1.5)
        FPD.add_Area_to_y_Constraint(kind = 'equality', value = carea, weight=1.)
        FPD.add_Area_to_y_Constraint(kind = 'max', value = carea+1.5)
        FPD.add_E1(kind='LS', weight = .1)
        FPD.add_E2(kind='LS', weight = 10.1)
        FPD.add_E3(kind='LS', weight = 10.1)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.1)
        L = Lagrangian(FPD)
        nump = 50
        mp = .4
        gfc = gfunc(curve5, b.CurvePoint, 0.,1., np.linspace(0., mp, nump), 1. )
        nump = 50
        mp = .2
        c1 = circle(origin=[0.,2.25],radius=2.25)
        e1 = ellipse(origin=[0.,rb_],ra=ra_, rb=rb_)
        circ_p = np.linspace(-np.pi/2., np.pi/4., nump)
        gfc3 = gfunc(curve5, e1.circle_pt, 0.,mp,circ_p, 10000.)
        #L.sp['crv']     = gfc
        L.sp['circ']    = gfc3
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve5, L, data = interval_data)
        Lspline.optimize(stop = 50)
        bbow.append(copy.deepcopy(Lspline.curve))
        bbowi = []
        for curve, x in zip(bbow,[0.,-.1,-.2,-.3,-.4]):
            nv = frame.add1D(curve.vertices, 0)
            nv[:,0] = x
            bbowi.append(spline.Bspline(vertices=nv, k=4, nump=30 ))
            
        curve6 = copy.deepcopy(Lspline.curve)
        
        CL = initial_curve([0.,0.],[0.,b.vertices[-1,1]],7)
        BL = initial_curve([-5.,0.],[-5.,b.vertices[-1,1]],7)
        
        curve6.plot()
        b.plot()
        BL.plot()
        CL.plot()
        
        
        
        
        
        
        
        
        ## minimal
        
        radius = 2.25
        rb_ = 1.8
        #for ra_ in [.5,.8,1.,.8,.4]:
        ra_ = 2.25
        #print ra_
        curve5 = copy.deepcopy(curve21)
        
        interval_data, small = interval_bounds(curve5)
        FPD = FormParameterDict(curve5) 
        #for s in np.linspace(0.1,.99,10):
        #    value_ = b.CurvePoint(s)
        #    FPD.add_xPointConstraint(kind='LS', location=s, value=value_[0], weight=10000.*1./(s))
        #    FPD.add_yPointConstraint(kind='LS', location=s, value=value_[1], weight=10000.*1./(s))


        FPD.add_Area_to_y_Constraint(kind = 'min', value = carea-1.5)
        #FPD.add_Area_to_y_Constraint(kind = 'equality', value = carea, weight=1.)
        FPD.add_Area_to_y_Constraint(kind = 'max', value = carea+1.5)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .1)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.1)
        L = Lagrangian(FPD)
        nump = 50
        mp = .4
        gfc = gfunc( copy.deepcopy(curve21), 
                    b.CurvePoint, 0.,1., np.linspace(0., mp, nump), 1. )
        nump = 50
        mp = .2
        c1 = circle(origin=[0.,2.25],radius=2.25)
        e1 = ellipse(origin=[0.,rb_],ra=ra_, rb=rb_)
        circ_p = np.linspace(-np.pi/2., np.pi/4., nump)
        gfc3 = gfunc(curve5, e1.circle_pt, 0.,mp,circ_p, 1000000.)
        L.sp['crv']     = gfc
        #L.sp['circ']    = gfc3
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve5, L, data = interval_data)
        Lspline.optimize(stop = 50)
        
        
        
        k_=4
        nump_ = 30
        db = surface(surfaces['Bow'])
        s1 = db.old_surf()
        f=s1.tcurveNet[16]
        ##
        ## Bow Transverse with Bulb
        ##
        a = spline.Bspline(vertices = s1.tcurveNet[16].vertices, k=4, nump=30)
        
        bv  = 30*s1.tcurveNet[16].vertices[:,1:3] + [0.,2.0934466019417477]
        b   = spline.Bspline(vertices = bv, k=4, nump=30)
        b.compute_area_to_y()
        carea = b.area_to_y
        b   = translate(b, 0., -bv[0,1])
        bv = copy.deepcopy(b.vertices)
        #b.plotcurve_detailed(curvature = 'yes', scale = 0.001 )
          
        n = 8
        curve = initial_curve(bv[0],bv[-1],n) 
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        #for s in np.linspace(0.01,1.,10):
        #    value_ = b.CurvePoint(s)
        #    FPD.add_xPointConstraint(kind='LS', location=s, value=value_[0], weight=10000.*1./s)
        #    FPD.add_yPointConstraint(kind='LS', location=s, value=value_[1], weight=10000.*1./s)
        FPD.add_Area_to_y_Constraint(kind = 'equality', value = carea)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .1)
        FPD.add_E3(kind='LS', weight = .1)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        
        c1 = circle(origin=[0.,2.25],radius=2.25)
        circ_p = np.linspace(-np.pi/2., np.pi/4., nump, endpoint=False)
        gfc3 = gfunc(curve5, c1.circle_pt, 0.,mp,circ_p, 100.)
        L.sp['circ']     = gfc3
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()