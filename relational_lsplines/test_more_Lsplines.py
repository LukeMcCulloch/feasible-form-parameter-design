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
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.axislines import SubplotZero
import matplotlib.patches as patches
from matplotlib.patches import Ellipse
import matplotlib.path    as path
mpl.rcParams['legend.handlelength'] = 0

import curve as spline

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
        self.fig.tight_layout()
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
                  font_size = 23):
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
        #self.ax.arrow(0., 0., max_arrow, 0., head_width=0.05, head_length=0.1, fc='k', ec='k')
        #self.ax.arrow(0., 0., 0., max_arrow, head_width=0.05, head_length=0.1, fc='k', ec='k')
        
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
                               docurvature = False, 
                               font_size = self.font_size)#, color = color_)
                self.legend_loc[1] -= 1.
            else:
                
                axi = plotfunc(object_, 
                               Lspline = Lspline_,
                               canvas = self.ax,
                               legendloc = self.legend_loc,
                               cn = const_num, 
                               font_size = self.font_size)#, color = color_)
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
            #name = str(self.x_axis_name)
            #plt.annotate(name, xy = (self.max_x,-1.), xycoords='data')
            #name = str(self.y_axis_name)
            #plt.annotate(name, xy = (-1.,self.max_y), xycoords='data')
        plt.legend(markerscale=0.)
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
        #plt.xlim(xi-excessx*dx,xe+excessx*dx)
        #plt.ylim(yi-excessy*dy,ye+excessy*dy)
        #plt.ylim(-.5,20.)        
        #name = str('x')#self.xe)
        #plt.annotate(name, xy = (self.max_x,-1.), xycoords='data')
                        
        #name = str('y')#self.ye)
        #plt.annotate(name, xy = (-1.,self.max_y), xycoords='data')
        
        plt.show()
        return
    
    def save_pdf(self, filename = None, ftype = '.pdf', closeit=True):
        """ save pdf file.
        No file extension needed.
        """
        if filename == None:
            filename = default_input('please enter a name for the picture', 'curve')
        plt.savefig(filename+ftype, bbox_inches = 'tight')
        if closeit:
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
        
if __name__ == "__main__":
    np.set_printoptions(threshold= 200)
    # TLM B-Spline Curve class
    import curve             as     spline
    from   initialValues     import InitializeControlPoints, InitializeControlVertices
    import copy
    from   ADILS             import IntervalLagrangeSpline, Lagrangian
    from   FormParameter     import FormParameterDict
    from   initialValues     import interval_bounds, lagrangian_bounds
    option  = 'ini_test'
    test    = 'basic'
        
    if option =='ini_test':
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 12.
        xe = 12.
        ye = 0.
        alphab = 0.#-5.
        alphae = 0.#-15.
        Cab_given = 0.
        Cae_given = 0.
        xc = 4.
        yc = 4.
        curve_area = 72.#84.#
        slope = 'down'
        
        ae = alphae
        ab = alphab
        
        numcv = 8 #7
        ini_v = InitializeControlVertices(alphae=ae, alphab=ab,
                                          Cab_given=0.,Cae_given=0.,
                                                              nCV = numcv)
        curve = spline.Bspline(ini_v.vertices, k, nump)  
        wi = 3.5 #2.5#
        w  = .5 #.65
        ep = 1.e-10
        sp = 1.e-4
        x = ini_v.vertices[:,0]
        y = ini_v.vertices[:,1]
        interval_data, small = interval_bounds(curve)
        
    if test == 'basic': 
        FPD = FormParameterDict(curve) 
        #**********************************************************************
        FPD.add_AreaConstraint(kind='equality', 
                               value = curve_area)#value = AF.fuzzyNumber(60.,70.,80.) )#
        FPD.add_AngleConstraint(kind='equality', 
                                location = 0., value = alphab)#AF.fuzzyNumber(-5.,-2.,0.))#
        FPD.add_AngleConstraint(kind='equality', 
                                location = 1., value = alphae)
        FPD.add_CurvatureConstraint(kind='equality', 
                                    location = 0., value = Cab_given)
        FPD.add_CurvatureConstraint(kind='equality', 
                                    location = 1., value = Cae_given)
        #**********************************************************************
        # inequalities are 'sticky' - should this be the case???
        #
        # in the absence of these maximums, the best curve
        # has x[3] = 4.80015552e+00
#        FPD.add_xVertexConstraint(kind = 'max',
#                                  value = 4.2,
#                                  index = 3)
        #
        # 5 is 'sticky' - until Feb 2018 change
        FPD.add_xVertexConstraint(kind = 'max',
                                  value = 5.0,
                                  index = 3)
        #
        # 6 we likley never make contact - not sticky
        #
#        FPD.add_xVertexConstraint(kind = 'max',
#                                  value = 6.0,
#                                  index = 3)
        #
#        FPD.add_xVertexConstraint(kind = 'min',
#                                  value = 4.9,
#                                  index = 3)
        #
#        FPD.add_xVertexConstraint(kind = 'min',
#                                  value = 4.9,
#                                  index = 3)
        #
        FPD.add_xVertexConstraint(kind = 'min',
                                  value = 4.6,
                                  index = 3)
        #
        #**********************************************************************
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        #FPD.add_ArcLength(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        #**********************************************************************
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        vertices = [Lspline.curve.xpts,Lspline.curve.ypts]
        Lspline.compute_lagrangian(vertices, Lspline.Lagrangian)
        Lspline.display(mytitle = 'initial_basic_curve')
        
        Lspline.curve.verbose = True
        self = Lspline
        vertices = [self.curve.xpts, self.curve.ypts]
        
        
        vertices = Lspline.optimize(vertices, stop = 100, Lagrangian = self.Lagrangian) 
        #Lspline.display(mytitle = 'final basic curve')
        
        Lspline.curve.plotcurve_detailed()
        
        print Lspline.curve.vertices