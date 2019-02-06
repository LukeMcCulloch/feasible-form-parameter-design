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
    #option  = 'ini_test'
    #test    = 'basic'
    option  = 'other'
    #test    = None
    #test    = 'fwd1'
    #test    = 'aft_section'
    #test    = 's_section'
    #test    = 'mid_section'
    #test    = 'bow_section2'
    #test    = 'SAC_1'
    #test    = 'SAC_2'  # curved sac Xcg variation
    #test    = 'SAC_3' #possible used byt SAC curves and SAC Xcg variation
    #test    = 'SAC_4'
    #test    = 'cusp1'
    #test    = 'exp1'
    #test    = 'exp2'
    #test    = 'exp3'
    #test    = 'exp4' # - like Birk's bow curve
    #test    = 'exp4_jan_25_16'
    #test    = 'exp5' #aft run
    #test    = 'exp6'#aft run, area to right edge..an idea that may not be necessary
    #test    = 'container'
    test    = 'dist_contraint'
    #test    = 'exp7' #fancy run 1  the old (fixed in space) run plot
    #test    = 'exp8'
    do_all = False
        
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
        FPD.add_AreaConstraint(kind='equality', value = curve_area)#value = AF.fuzzyNumber(60.,70.,80.) )#
        FPD.add_AngleConstraint(kind='equality', location = 0., value = alphab)#AF.fuzzyNumber(-5.,-2.,0.))#
        FPD.add_AngleConstraint(kind='equality', location = 1., value = alphae)
        FPD.add_CurvatureConstraint(kind='equality', location = 0., value = Cab_given)
        FPD.add_CurvatureConstraint(kind='equality', location = 1., value = Cae_given)
        #"""
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        #FPD.add_ArcLength(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        vertices = [Lspline.curve.xpts,Lspline.curve.ypts]
        Lspline.compute_lagrangian(vertices, Lspline.Lagrangian)
        Lspline.display(mytitle = 'initial_basic_curve')
        
        Lspline.curve.verbose = True
        self = Lspline
        vertices = [self.curve.xpts, self.curve.ypts]
        
        
        vertices = Lspline.optimize(vertices, stop = 35, Lagrangian = self.Lagrangian) 
        Lspline.display(mytitle = 'final basic curve')
        
    if test == 'fwd1':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 12.
        xe = 12.
        ye = 0.
        alphab = 75.#-5.
        alphae = 0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 58.#84.#
        slope = 'down'
        
        ae = alphae
        ab = alphab
        
        numcv = 8 #7
        ini_v = InitializeControlVertices(alphae=ae, alphab=ab,
                                          Cab_given=Cab,Cae_given=Cae,
                                           nCV = numcv, slope = 'down')
        #self = ini_v
        #ini_v.buildvertices()
        self = ini_v
        curve = spline.Bspline(ini_v.vertices, k, nump)  
        wi = 3.5 #2.5#
        w  = .5 #.65
        ep = 1.e-10
        sp = 1.e-4
        x = ini_v.vertices[:,0]
        y = ini_v.vertices[:,1]
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = curve_area)#value = AF.fuzzyNumber(60.,70.,80.) )#
        FPD.add_AngleConstraint(kind='equality', location = 0., value = -alphab)#AF.fuzzyNumber(-5.,-2.,0.))#
        FPD.add_AngleConstraint(kind='equality', location = 1., value = alphae)
        FPD.add_AngleConstraint(kind='equality', location = .5, value = -20.)
        FPD.add_xPointConstraint(kind='equality', location=.6, value = 7.)
        FPD.add_yPointConstraint(kind='equality', location=.6, value = 5.)
        #FPD.add_CurvatureConstraint(kind='equality', location = .5, value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 0., value = Cab_given)
        #FPD.add_CurvatureConstraint(kind='equality', location = 1., value = Cae_given)
        FPD.add_YcConstraint(kind='equality', value = 3.2)
        #"""
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        #FPD.add_ArcLength(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        
#        for el in Lspline.Lagrangian.equality:
#            Lspline.Lagrangian.equality[el].
        vertices = [Lspline.curve.xpts,Lspline.curve.ypts]
        Lspline.compute_lagrangian(vertices, Lspline.Lagrangian)
        Lspline.display(mytitle = 'initial_bow_curve')
        Lspline.optimize()
        Lspline.display(mytitle = 'final_bow_curve')
        #Lspline.display(mytitle = 'test')
        
    if do_all:
        test = 'aft_section'
    if test == 'aft_section':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 12.
        xe = 12.
        ye = 0.
        alphab = 75.#-5.
        alphae = 89.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 60.#58.#84.#
        slope = 'down'
        
        ae = alphae
        ab = alphab
        
        numcv = 8 #7
        ini_v = InitializeControlVertices(alphae=ae, alphab=ab,
                                          Cab_given=Cab,Cae_given=Cae,
                                           nCV = numcv, slope = 'down')
        #self = ini_v
        #ini_v.buildvertices()
        self = ini_v
        curve = spline.Bspline(ini_v.vertices, k, nump)  
        x = ini_v.vertices[:,0]
        y = ini_v.vertices[:,1]
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = curve_area)
        
        FPD.add_AngleConstraint(kind='equality', location = 0., value = -alphab)
        FPD.add_verticalAngleConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_YcConstraint(kind='equality', value = 3.0 )
        #"""
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        #vertices = [Lspline.curve.xpts,Lspline.curve.ypts]

        Lspline.compute_lagrangian()
        Lspline.display(mytitle = 'initial_aft_section')
        Lspline.optimize()
        Lspline.display(mytitle = 'aft_section')
        Lspline.rotate(180.)
    if test == 'fwd3':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 0.
        xe = 12.
        ye = 12.
        alphab = -75.#-5.
        alphae = 0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 58.#84.#
        slope = 'up'
        
        ae = alphae
        ab = alphab
        
        numcv = 8 #8
        ini_v = InitializeControlVertices(xb,yb,xe,ye,alphae=ae, alphab=ab,
                                          Cab_given=Cab,Cae_given=Cae,
                                           nCV = 8, slope = 'up')
#        ini_v.vertices = ini_v.vertices[::-1]
#        xv = copy.deepcopy(ini_v.vertices[:,1])
#        yv = copy.deepcopy(ini_v.vertices[:,0])
#        ini_v.vertices[:,0] = xv
#        ini_v.vertices[:,1] = yv
        
        curve = spline.Bspline(ini_v.vertices, k, nump)  

#
#        x = ini_v.vertices[:,0]
#        y = ini_v.vertices[:,1]
#        interval_data, small = interval_bounds(curve)
#        
#        FPD = FormParameterDict(curve) 
#        FPD.add_AreaConstraint(kind='equality', value = curve_area)#value = AF.fuzzyNumber(60.,70.,80.) )#
#        FPD.add_AngleConstraint(kind='equality', location = 0., value = -alphab)#AF.fuzzyNumber(-5.,-2.,0.))#
#        FPD.add_AngleConstraint(kind='equality', location = 1., value = alphae)
#        FPD.add_AngleConstraint(kind='equality', location = .5, value = -20.)
#        FPD.add_xPointConstraint(kind='equality', location=.6, value = 7.)
#        FPD.add_yPointConstraint(kind='equality', location=.6, value = 5.)
#        #FPD.add_CurvatureConstraint(kind='equality', location = .5, value = 0.)
#        #FPD.add_CurvatureConstraint(kind='equality', location = 0., value = Cab_given)
#        #FPD.add_CurvatureConstraint(kind='equality', location = 1., value = Cae_given)
#        FPD.add_YcConstraint(kind='equality', value = 3.2)
#        #"""
#        FPD.add_E1(kind='LS', weight = 1.)
#        FPD.add_E2(kind='LS', weight = 1.)
#        FPD.add_E3(kind='LS', weight = 1.)
#        #FPD.add_ArcLength(kind='LS', weight = 1.)
#        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
#        L = Lagrangian(FPD)
#        
#        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
#        
#        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
#        
#
#        vertices = [Lspline.curve.xpts,Lspline.curve.ypts]
#        Lspline.compute_lagrangian(vertices, Lspline.Lagrangian)
#        Lspline.display(mytitle = 'initial_bow_curve3')
#        Lspline.optimize()
#        Lspline.display(mytitle = 'final_bow_curve3')
#        #Lspline.display(mytitle = 'test')
    if do_all:
        test = 's_section'
    if test == 's_section':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 0.
        xe = 12.
        ye = 12.
        alphab = 15.#-5.
        alphae = 15.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 72.#84.#
        slope = 'up'
        
        ae = alphae
        ab = alphab
        
        numcv = 8 #10
        ini_v = InitializeControlVertices(xb,yb,xe,ye,alphae=ae, alphab=ab,
                                          Cab_given=Cab,Cae_given=Cae,
                                           nCV = numcv, slope = 'up')

        
        curve = spline.Bspline(ini_v.vertices, k, nump)  
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = curve_area)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = alphab)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = alphae)
        FPD.add_CurvatureConstraint(kind='equality', location = 0.0, value = Cab)
        FPD.add_CurvatureConstraint(kind='equality', location = 0.1, value = Cab)
        FPD.add_CurvatureConstraint(kind='equality', location = 1.0, value = Cae)        
        FPD.add_CurvatureConstraint(kind='equality', location = 0.9, value = Cae)
        #FPD.add_xPointConstraint(kind='equality', location=.75, value = 9.)
        #FPD.add_xPointConstraint(kind='equality', location=.25, value = 3.)
        #FPD.add_YcConstraint(kind='equality', value = 3.0 )
        #"""
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        
        Lspline.compute_lagrangian()
        Lspline.display(mytitle = 'initial_s_section')
        Lspline.optimize()
        Lspline.display(mytitle = 's_section')

    if do_all:
        test = 'mid_section'
    if test == 'mid_section':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 0.
        xe = 12.
        ye = 12.
        alphab = 5.#-5.
        alphae = 75.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 35.
        slope = 'up'
        
        ae = alphae
        ab = alphab
        
        numcv = 8 #10
        ini_v = InitializeControlVertices(xb,yb,xe,ye,alphae=ae, alphab=ab,
                                          Cab_given=Cab,Cae_given=Cae,
                                           nCV = numcv, slope = 'up')

        
        curve = spline.Bspline(ini_v.vertices, k, nump)  
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        #FPD.add_AreaConstraint(kind='equality', value = 25.66)
        FPD.add_Area_to_y_Constraint(kind='equality', value = 118.34)
        
        FPD.add_AngleConstraint(kind='equality', location = 0., value = alphab)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = alphae)
        
        FPD.add_CurvatureConstraint(kind='equality', location = 0.0, value = Cab)
        FPD.add_CurvatureConstraint(kind='equality', location = 0.15, value = Cab)
        
        FPD.add_CurvatureConstraint(kind='equality', location = 1.0, value = Cae)        
        FPD.add_CurvatureConstraint(kind='equality', location = 0.85, value = Cae)
        
        #FPD.add_xPointConstraint(kind='equality', location=.2, value = 1.)
        #FPD.add_yPointConstraint(kind='equality', location=.1, value = .33333333333333333)
        
        #FPD.add_xPointConstraint(kind='equality', location=.35, value = 3.)
        #FPD.add_yPointConstraint(kind='equality', location=.35, value = 1.)
        
        #FPD.add_xPointConstraint(kind='equality', location=.75, value = 9.)
        #FPD.add_yPointConstraint(kind='equality', location=.75, value = 9.)
        
        #FPD.add_YcConstraint(kind='equality', value = 3.0 )
        #"""
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)#, fritzjohn_conditions=False)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)#, fritz_john_conditions = False)
        
        Lspline.compute_lagrangian()
        Lspline.display(mytitle = 'initial_mid_section')
        Lspline.optimize()
        #Lspline.curve.plot()
        Lspline.display(mytitle = 'mid_section')
        
        #for xc in []
        
    if do_all:
        test = 'bow_section2'
    if test == 'bow_section2':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 0.
        xe = 12.
        ye = 12.
        alphab = 5.#-5.
        alphae = 75.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 35.
        slope = 'up'
        
        ae = alphae
        ab = alphab
        
        numcv = 8 #10
        ini_v = InitializeControlVertices(xb,yb,xe,ye,alphae=ae, alphab=ab,
                                          Cab_given=Cab,Cae_given=Cae,
                                           nCV = numcv, slope = 'up')

        
        curve = spline.Bspline(ini_v.vertices, k, nump)  
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        #FPD.add_AreaConstraint(kind='equality', value = 25.66)
        FPD.add_Area_to_y_Constraint(kind='equality', value = 98.34)
        
        FPD.add_AngleConstraint(kind='equality', location = 0., value = alphab)
        FPD.add_verticalAngleConstraint(kind='equality', location = 1., value = 0.)
        
        #FPD.add_CurvatureConstraint(kind='equality', location = 0.0, value = Cab)
        
        ###FPD.add_CurvatureConstraint(kind='equality', location = 0.15, value = Cab)

        FPD.add_CurvatureConstraint(kind='equality', location = 0.5, value = -5.)
        
        #FPD.add_CurvatureConstraint(kind='equality', location = 1.0, value = Cae)        
        
        ###FPD.add_CurvatureConstraint(kind='equality', location = 0.85, value = Cae)
        
        #FPD.add_xPointConstraint(kind='equality', location=.2, value = 1.)
        #FPD.add_yPointConstraint(kind='equality', location=.1, value = .33333333333333333)
        
        #FPD.add_xPointConstraint(kind='equality', location=.35, value = 3.)
        #FPD.add_yPointConstraint(kind='equality', location=.35, value = 1.)
        
        #FPD.add_xPointConstraint(kind='equality', location=.75, value = 9.)
        #FPD.add_yPointConstraint(kind='equality', location=.75, value = 9.)
        
        #FPD.add_YcConstraint(kind='equality', value = 3.0 )
        #"""
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        
        Lspline.compute_lagrangian()
        #Lspline.display(mytitle = 'initial_bow_section2')
        Lspline.optimize()
        #Lspline.curve.plot()
        #Lspline.display(mytitle = 'bow_section2')

    if do_all:
        test = 'aft_section2'
    if test == 'aft_section2':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 0.
        xe = 12.
        ye = 12.
        alphab = 75.#-5.
        alphae = 89.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 60.#58.#84.#
        slope = 'down'
        
        ae = alphae
        ab = alphab
        
        numcv = 8
        ini_v = InitializeControlVertices(xb,yb,xe,ye,alphae=ae, alphab=ab,
                                          Cab_given=Cab,Cae_given=Cae,
                                           nCV = numcv, slope = 'up')
        #self = ini_v
        #ini_v.buildvertices()
        self = ini_v
        curve = spline.Bspline(ini_v.vertices, k, nump)  
        x = ini_v.vertices[:,0]
        y = ini_v.vertices[:,1]
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        FPD.add_Area_to_y_Constraint(kind='equality', value = curve_area)
        
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = alphab)
        
        
        #FPD.add_Yc_to_yConstraint(kind='equality', value = 7.0 )
        #"""
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        #vertices = [Lspline.curve.xpts,Lspline.curve.ypts]

        Lspline.compute_lagrangian()
        Lspline.display(mytitle = 'initial_aft_section2')
        Lspline.optimize()
        Lspline.display(mytitle = 'aft_section2')
        
    if do_all:
        test = 'SAC_1'
    if test == 'SAC_1':
        print test
        Amax_   = 12.
        l1_     = 12.
        l2_     = 12.
        l3_     = 12.
        Xe_     = 9.
        Xm_     = 6.
        Xr_     = 3.
        Xlcb_   = 18.
        Afp_    = 1.
        Aap_    = 1.
        Cpe_    = 0.
        Cpr_    = 0.
        sac = SAC(Amax=Amax_, #v1,v2,v3, 
                 l1=l1_,l2=l2_,l3=l3_, 
                 Xe=Xe_, Xm=Xm_, 
                 Xr=Xr_, Xlcb=Xlcb_, 
                 Afp=Afp_, Aap=Aap_, 
                 Cpe=Cpe_, Cpr=Cpr_)
        sac.plot()
        
        store1 = copy.deepcopy(sac.SAC)
        curve = copy.deepcopy(sac.SAC)
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)
        

        FPD.add_xPointConstraint(kind='equality', location = 0.333, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.333, value = 12.)
        FPD.add_xPointConstraint(kind='equality', location = 0.666, value = 24.)
        FPD.add_yPointConstraint(kind='equality', location = 0.666, value = 12.)
        
        FPD.add_xPointConstraint(kind='equality', location = 0.5, value = 18.)
        FPD.add_yPointConstraint(kind='equality', location = 0.5, value = 12.)
        
        FPD.add_CurvatureConstraint(kind='equality', location = .333, value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .666, value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = .333, value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = .666, value = 0.)
        
        FPD.add_CurvatureConstraint(kind='equality', location = .5, value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = .5, value = 0.)
        
        FPD.add_XcConstraint(kind = 'equality', value = 18.)
        
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        
        #FPD.add_verticalAngleConstraint(kind='equality', location = 1., value = 0.)
        
        #FPD.add_Yc_to_yConstraint(kind='equality', value = 7.0 )
        #"""
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        #vertices = [Lspline.curve.xpts,Lspline.curve.ypts]

        Lspline.optimize()
        Lspline.display(mytitle = 'flat_integrated_SAC')
        
        store2 = copy.deepcopy(sac.SAC)

        curve = copy.deepcopy(sac.SAC)
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)
        

        FPD.add_XcConstraint(kind = 'equality', value = 16.)
        #"""
        #FPD.add_xPointConstraint(kind='equality', location = 0.333, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.333, value = 12.)
        #FPD.add_xPointConstraint(kind='equality', location = 0.666, value = 24.)
        FPD.add_yPointConstraint(kind='equality', location = 0.666, value = 12.)
        
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 0., value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        #"""
        #FPD.add_verticalAngleConstraint(kind='equality', location = 1., value = 0.)
        
        #FPD.add_Yc_to_yConstraint(kind='equality', value = 7.0 )
        #"""
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)

        Lspline.optimize()
        Lspline.display(mytitle = 'curved_integrated_SAC_xcg_16test_simple')
    
    if do_all:
        test = 'SAC_2'
    if test == 'SAC_2':
        print test
        Amax_   = 12.
        l1_     = 12.
        l2_     = 12.
        l3_     = 12.
        Xe_     = 9.
        Xm_     = 6.
        Xr_     = 3.
        Xlcb_   = 18.
        Afp_    = 1.
        Aap_    = 1.
        Cpe_    = 0.
        Cpr_    = 0.
        
        #        ini_v = InitializeControlVertices(xb=0.,yb=0.,
        #                                          xe=0.,
        #                                          ye=0.,
        #                                          nCV = 10, 
        #                                          slope = 'up')
        curve = initial_curve((0.,0.),(36.,0.), num=10, k=4,nump=30)
        v = copy.deepcopy(curve.vertices)
        #v[4:6,1] = 12.
        v[3,1]=12.
        v[4,1]=12.
        v[2,1]=5.
        v[5,1]=5.
        curve.vertices = v #updates the curve automagically
        
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=18.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 0.)
        #FPD.add_yPointConstraint(kind='equality', location = 0.333, value = 12.)
        #FPD.add_yPointConstraint(kind='equality', location = 0.666, value = 12.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 0., value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_E1(kind='LS', weight = .1)
        FPD.add_E2(kind='LS', weight = .5)
        #FPD.add_E3(kind='LS', weight = .1)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)

        Lspline.optimize()
        #Lspline.curve.plot()
        Lspline.display(mytitle = 'curved_SAC_xcg_18test_simple')
        
        #Lspline.conv = 1.
        #Lspline.Lagrangian.equality[1].pass_value = 16.
        #Lspline.Lagrangian.equality[1].value = 16.
        #Lspline.Lagrangian.equality[1].computed_value = 0.
        save1 = copy.deepcopy(Lspline.curve)
        
        curve = Lspline.curve
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=16.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = -0.)
        #FPD.add_yPointConstraint(kind='equality', location = 0.333, value = 12.)
        #FPD.add_yPointConstraint(kind='equality', location = 0.666, value = 12.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 0., value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_E1(kind='LS', weight = .1)
        FPD.add_E2(kind='LS', weight = .5)
        #FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)

        Lspline.optimize()
        #Lspline.curve.plot()
        Lspline.display(mytitle = 'curved_SAC_xcg_16test_simple')
        
        save2 = copy.deepcopy(Lspline.curve)
        curve = Lspline.curve
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=20.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = -0.)
        #FPD.add_yPointConstraint(kind='equality', location = 0.333, value = 12.)
        #FPD.add_yPointConstraint(kind='equality', location = 0.666, value = 12.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 0., value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_E1(kind='LS', weight = .1)
        FPD.add_E2(kind='LS', weight = .5)
        #FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
        Lspline.display(mytitle = 'curved_SAC_xcg_20test_simple')
        
        
    if test == 'SAC_3':
        print test

        curve = initial_curve((0.,0.),(36.,0.), num=12, k=4,nump=30)
        v = copy.deepcopy(curve.vertices)
        v[4:7,1] = 15.
        
        curve.vertices = v #updates the curve automagically
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=18.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = -0.)
        FPD.add_yPointConstraint(kind='equality', location = 0.333, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.666, value = 12.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 0., value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_E1(kind='LS', weight = .1)
        FPD.add_E2(kind='LS', weight = .5)
        #FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
        Lspline.display(mytitle = 'curved_SAC_xcg_18test')
        save1 = copy.deepcopy(Lspline.curve)
        
        
        curve = Lspline.curve
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=18.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = -0.)
        FPD.add_yPointConstraint(kind='equality', location = 0.333, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.666, value = 12.)
        FPD.add_AngleConstraint(kind='equality', location = .333, value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = .666, value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 0., value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
        Lspline.display(mytitle = 'flat_SAC_xcg_18test')
        save2 = copy.deepcopy(Lspline.curve)
        
        curve = Lspline.curve
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=18.)
        #FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        #FPD.add_AngleConstraint(kind='equality', location = 1., value = -0.)
        FPD.add_yPointConstraint(kind='equality', location = 0.333, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.666, value = 12.)
        FPD.add_AngleConstraint(kind='equality', location = .333, value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = .666, value = 0.)
        FPD.add_yPointConstraint(kind='equality', location = .5, value = 12.)
        #FPD.add_AngleConstraint(kind='equality', location = .5, value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = .5, value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = .333, value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = .666, value = 0.)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
        Lspline.display(mytitle = 'flat_SAC_xcg_18test2')
        save3 = copy.deepcopy(Lspline.curve)
        
        
        curve = Lspline.curve
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=18.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = -0.)
        FPD.add_yPointConstraint(kind='equality', location = 0.333, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.666, value = 12.)
        
        FPD.add_xPointConstraint(kind='equality', location = 0.333, value = 12.)
        FPD.add_xPointConstraint(kind='equality', location = 0.666, value = 24.)
        
        FPD.add_AngleConstraint(kind='equality', location = .333, value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = .666, value = 0.)
        FPD.add_yPointConstraint(kind='equality', location = .5, value = 12.)
        #FPD.add_AngleConstraint(kind='equality', location = .5, value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = .333, value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = .666, value = 0.)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
        Lspline.display(mytitle = 'flat_SAC_xcg_18test3')
        save4 = copy.deepcopy(Lspline.curve)
        
        curve = Lspline.curve
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=18.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = -0.)
        FPD.add_yPointConstraint(kind='equality', location = 0.3333333333, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.6666666666, value = 12.)
        FPD.add_xPointConstraint(kind='equality', location = 0.3333333333, value = 12.)
        FPD.add_xPointConstraint(kind='equality', location = 0.6666666666, value = 24.)
        FPD.add_AngleConstraint(kind='equality', location = .3333333333, value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = .6666666666, value = 0.)
        FPD.add_yPointConstraint(kind='equality', location = .5, value = 12.)
        #FPD.add_AngleConstraint(kind='equality', location = .5, value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .3333333333, value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .6666666666, value = 0.)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
        Lspline.display(mytitle = 'flat_SAC_xcg_18test4')
        save5 = copy.deepcopy(Lspline.curve)

    if test == 'SAC_4':
        print test

        curve = initial_curve((0.,0.),(36.,0.), num=12, k=4,nump=30)
        v = copy.deepcopy(curve.vertices)
        v[4:7,1] = 15.
        
        curve.vertices = v #updates the curve automagically
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=18.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_yPointConstraint(kind='equality', location = 0.3333333333, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.6666666666, value = 12.)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
        #Lspline.display(mytitle = 'curved_SAC_xcg_18test_4_0')
        save0 = copy.deepcopy(Lspline.curve)
        
        
        curve = Lspline.curve
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=18.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 0.)        
        FPD.add_yPointConstraint(kind='equality', location = 0.3333333333, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.5, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.6666666666, value = 12.)
        FPD.add_xPointConstraint(kind='equality', location = 0.3333333333, value = 12.)
        FPD.add_xPointConstraint(kind='equality', location = 0.6666666666, value = 24.)
        FPD.add_AngleConstraint(kind='equality', location = .3333333333, value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = .6666666666, value = 0.)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'flat_SAC_xcg_18test_4_1')
        save1 = copy.deepcopy(Lspline.curve) 
        
        
        curve = Lspline.curve
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=18.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_yPointConstraint(kind='equality', location = 0.3333333333, value = 12.)
        #FPD.add_yPointConstraint(kind='equality', location = 0.5, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.6666666666, value = 12.)
        FPD.add_xPointConstraint(kind='equality', location = 0.3333333333, value = 12.)
        FPD.add_xPointConstraint(kind='equality', location = 0.6666666666, value = 24.)
        FPD.add_AngleConstraint(kind='equality', location = .3333333333, value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = .6666666666, value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .3333333333, value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .6666666666, value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .0, value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'flat_SAC_xcg_18test_4_2')
        save2 = copy.deepcopy(Lspline.curve) 
        
        curve = Lspline.curve
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=16.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_yPointConstraint(kind='equality', location = 0.3333333333, value = 12.)
        #FPD.add_yPointConstraint(kind='equality', location = 0.5, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.6666666666, value = 12.)
        FPD.add_xPointConstraint(kind='equality', location = 0.3333333333, value = 12.)
        FPD.add_xPointConstraint(kind='equality', location = 0.6666666666, value = 24.)
        FPD.add_AngleConstraint(kind='equality', location = .3333333333, value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = .6666666666, value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .3333333333, value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .6666666666, value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .0, value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'flat_SAC_xcg_16test_4_3')
        save3 = copy.deepcopy(Lspline.curve)
        
        curve = copy.deepcopy(save2)
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=20.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_yPointConstraint(kind='equality', location = 0.3333333333, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.6666666666, value = 12.)
        FPD.add_xPointConstraint(kind='equality', location = 0.3333333333, value = 12.)
        FPD.add_xPointConstraint(kind='equality', location = 0.6666666666, value = 24.)
        FPD.add_AngleConstraint(kind='equality', location = .3333333333, value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = .6666666666, value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .3333333333, value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .6666666666, value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .0, value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'flat_SAC_xcg_20test_4_4')
        save4 = copy.deepcopy(Lspline.curve)
        
    if do_all:
        test = 'cusp1'
    if test == 'cusp1':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 12.
        xe = 12.
        ye = 0.
        alphab = 0.#-5.
        alphae = 0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 72.#58.#84.#
        slope = 'down'
        
        ae = alphae
        ab = alphab
        
        curve = initial_curve([xb,yb],[xe,ye],11) 

        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = curve_area)
        
        FPD.add_xVertexConstraint(kind = 'equality', value=6.5, index =4)
        FPD.add_yVertexConstraint(kind = 'equality', value=7.5, index =4)
        
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 0.)
        
        
        #FPD.add_Yc_to_yConstraint(kind='equality', value = 7.0 )
        #"""
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 0.)
        #FPD.add_ArcLength(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        FPD.add_xAttractor(kind='LS', weight = 1.05, index = 2)
        FPD.add_xAttractor(kind='LS', weight = 1.05, index = 3)
        FPD.add_yAttractor(kind='LS', weight = 1.05, index = 2)
        FPD.add_yAttractor(kind='LS', weight = 1.05, index = 3)

        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        
        #curve.data['area']={'desired': 65.,'weight':.01} #300
        """
        curve.data['force'] = {'vertices':[2,3],'weights':[(1.05,1.05),(1.05,1.05)]} #0.4,0.05,
        #"""
        #curve.data['force'] = {'vertices':[1,2,3],'weights':[(.5,.05,.5),(.5,.05,.5)]} #0.4,0.05,
        curve.data['knuckle']=[{'coincident':[5,6],'freevertex':4}]
        curve.data['barrier'] = False
        #curve.data['knuckle']=[{'coincident':[3,5],'freevertex':4},
        #                       {'coincident':[7],'freevertex':6}] #should make point 2 and 3 coincident
        
        
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        #vertices = [Lspline.curve.xpts,Lspline.curve.ypts]

        Lspline.compute_lagrangian()
        Lspline.display(mytitle = 'initial_cusp_section1')
        Lspline.optimize()
        Lspline.display(mytitle = 'cusp_section1')
        
    if do_all:
        test = 'exp1'
    if test == 'exp1':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 0.
        xe = 12.
        ye = 12.
        alphab = 0.#-5.
        alphae = 0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 72.#58.#84.#
        slope = 'down'
        
        ae = alphae
        ab = alphab
        
        curve = initial_curve([xb,yb],[xe,ye],11) 

        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        #FPD.add_AreaConstraint(kind='equality', value = curve_area)
        FPD.add_Area_to_y_Constraint(kind='equality', value = curve_area)
        
        #FPD.add_xVertexConstraint(kind = 'equality', value=6.5, index =4)
        #FPD.add_yVertexConstraint(kind = 'equality', value=7.5, index =4)
        
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 0.)
        
        
        #FPD.add_Yc_to_yConstraint(kind='equality', value = 7.0 )
        #"""
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 0.)
        #FPD.add_ArcLength(kind='LS', weight = 1.)
        #"""
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        FPD.add_xAttractor(kind='LS', weight = 1.05, index = 2)
        FPD.add_xAttractor(kind='LS', weight = 1.05, index = 3)
        FPD.add_yAttractor(kind='LS', weight = 1.05, index = 2)
        FPD.add_yAttractor(kind='LS', weight = 1.05, index = 3)
        #"""
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        
        #curve.data['area']={'desired': 65.,'weight':.01} #300
        """
        curve.data['force'] = {'vertices':[2,3],'weights':[(1.05,1.05),(1.05,1.05)]} #0.4,0.05,
        #"""
        #curve.data['force'] = {'vertices':[1,2,3],'weights':[(.5,.05,.5),(.5,.05,.5)]} #0.4,0.05,
        curve.data['knuckle']=[{'coincident':[5,6],'freevertex':4}]
        curve.data['barrier'] = False
        #curve.data['knuckle']=[{'coincident':[3,5],'freevertex':4},
        #                       {'coincident':[7],'freevertex':6}] #should make point 2 and 3 coincident
        
        
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        #vertices = [Lspline.curve.xpts,Lspline.curve.ypts]

        Lspline.compute_lagrangian()
        Lspline.display(mytitle = 'exp1_ini')
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'exp1')
    if do_all:
        test = 'exp2'
    if test == 'exp2':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 0.
        xe = 12.
        ye = 12.
        alphab = 0.#-5.
        alphae = 0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 72.#58.#84.#
        slope = 'down'
        
        ae = alphae
        ab = alphab
        
        curve = initial_curve([xb,yb],[xe,ye],11) 

        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 

        
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 0.)
        
        #"""
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        #FPD.add_ArcLength(kind='LS', weight = 1.)
        #"""
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        FPD.add_xAttractor(kind='LS', weight = 1.05, index = 2)
        FPD.add_xAttractor(kind='LS', weight = 1.05, index = 3)
        FPD.add_yAttractor(kind='LS', weight = 1.05, index = 2)
        FPD.add_yAttractor(kind='LS', weight = 1.05, index = 3)
        #"""
        FPD.add_AreaConstraint(kind='LS',value = 72.)
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
               
        curve.data['knuckle']=[{'coincident':[5,6],'freevertex':4}]
        curve.data['barrier'] = False

        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)

        Lspline.compute_lagrangian()
        Lspline.display(mytitle = 'exp2_ini')
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'exp2')
        
        
    if do_all:
        test = 'exp3'
    if test == 'exp3':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 0.
        xe = 15.
        ye = 12.
        alphab = 0.#-5.
        alphae = 0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 90.#58.#84.#
        slope = 'down'
        
        ae = alphae
        ab = alphab
        
        curve = initial_curve([xb,yb],[xe,ye],11) 

        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = curve_area)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 20.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 30.)
        FPD.add_CurvatureConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp30')
        curve = Lspline.curve
        
        
        
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = curve_area)
        FPD.add_xVertexConstraint(kind = 'equality', value=6., index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=5., index =5)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 20.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 30.)
        FPD.add_CurvatureConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp31')
        curve = Lspline.curve
        
        
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 85)
        FPD.add_xVertexConstraint(kind = 'equality', value=7., index =6)
        FPD.add_yVertexConstraint(kind = 'equality', value=8., index =6)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 20.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 30.)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'exp32')
        curve = Lspline.curve
        
        
        
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 90)
        #FPD.add_Area_to_y_Constraint(kind='equality', value = curve_area)
        
        FPD.add_xVertexConstraint(kind = 'equality', value=4.5, index =4)
        FPD.add_yVertexConstraint(kind = 'equality', value=3., index =4)
        FPD.add_xVertexConstraint(kind = 'equality', value=4., index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=7., index =5)
        FPD.add_xVertexConstraint(kind = 'equality', value=7., index =6)
        FPD.add_yVertexConstraint(kind = 'equality', value=7.5, index =6)
        FPD.add_xVertexConstraint(kind = 'equality', value=6.5, index =7)
        FPD.add_yVertexConstraint(kind = 'equality', value=4.5, index =7)
        
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 20.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 30.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 0., value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        
        #FPD.add_Yc_to_yConstraint(kind='equality', value = 7.0 )
        #"""
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        """
        #FPD.add_xAttractor(kind='LS', weight = 1.05, index = 3)
        #FPD.add_yAttractor(kind='LS', weight = 1.05, index = 3)
        
        FPD.add_xAttractor(kind='LS', weight = -1.05, index = 4 )
        FPD.add_yAttractor(kind='LS', weight = -1.05, index = 4)
        
        FPD.add_xAttractor(kind='LS', weight = -1.05, index = 5 )
        FPD.add_yAttractor(kind='LS', weight = -1.05, index = 5)
        #"""
        #FPD.add_AreaConstraint(kind='LS',value = 72.)
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        
        #curve.data['area']={'desired': 65.,'weight':.01} #300
        """
        curve.data['force'] = {'vertices':[2,3],'weights':[(1.05,1.05),(1.05,1.05)]} #0.4,0.05,
        #"""
        #curve.data['force'] = {'vertices':[1,2,3],'weights':[(.5,.05,.5),(.5,.05,.5)]} #0.4,0.05,
        """        
        curve.data['knuckle']=[{'coincident':[5,6],'freevertex':4}]
        curve.data['barrier'] = False
        #"""
        #curve.data['knuckle']=[{'coincident':[3,5],'freevertex':4},
        #                       {'coincident':[7],'freevertex':6}] #should make point 2 and 3 coincident
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'exp33')
        curve = Lspline.curve
    
    if do_all:
        test = 'exp4_jan_25_16'
    if test == 'exp4_jan_25_16':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 0.
        xe = 15.
        ye = 12.
        alphab = 0.#-5.
        alphae = 0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 80.#58.#84.#
        slope = 'down'
        
        ae = alphae
        ab = alphab
        
        curve = initial_curve([xb,yb],[xe,ye],8) 
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        #FPD.add_AreaConstraint(kind='equality', value = 96)
        FPD.add_xVertexConstraint(kind = 'equality', value=4., index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=3., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=8., index =3)
        FPD.add_xVertexConstraint(kind = 'equality', value=9., index =4)
        FPD.add_yVertexConstraint(kind = 'equality', value=9.5, index =4)
        FPD.add_xVertexConstraint(kind = 'equality', value=6.5, index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=4.5, index =5)
        
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 20.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 60.)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp4')
        curve = Lspline.curve
        
        
        #"""
        FPD = FormParameterDict(curve) 
        FPD.add_Area_to_y_Constraint(kind='equality', value = 84.)
        FPD.add_xVertexConstraint(kind = 'LS', value=4., index =2, weight = 10.)
        FPD.add_yVertexConstraint(kind = 'LS', value=3., index =2, weight = 10.)
        FPD.add_xVertexConstraint(kind = 'LS', value=3.5, index =3, weight = 1000000.)
        FPD.add_yVertexConstraint(kind = 'LS', value=8., index =3, weight = 1000000.)
        FPD.add_xVertexConstraint(kind = 'LS', value=9., index =4, weight = 1000000.)
        FPD.add_yVertexConstraint(kind = 'LS', value=10., index =4, weight = 1000000.)
        FPD.add_xVertexConstraint(kind = 'LS', value=5.5, index =5, weight = 1000000.)
        FPD.add_yVertexConstraint(kind = 'LS', value=4.5, index =5, weight = 1000000.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 20.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 60.)
        FPD.add_E1(kind='LS', weight = .2)
        FPD.add_E2(kind='LS', weight = .2)
        FPD.add_E3(kind='LS', weight = 2.)
        FPD.add_ArcLengthApprox(kind='LS', weight = .2)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp42')
        curve = Lspline.curve
        #"""
        curve.compute_area_to_y()
        curve.compute_moments_to_y()
        curve.computeXc_to_y()
        curve.computeYc_to_y()
        curve.compute_area()
        curve.compute_moments()
        curve.computeXc()
        curve.computeYc()
        store1 = copy.deepcopy(curve)
        
        ##
        ## Best fwd curve in green imediately below
        ##
        
        #"""
        #""" plot it
        curve.plotcurve_detailed()
        #curve.plot_area_to_x_shiftaxis(shift_axis = 12.)
        Xc_y = curve.Xc_to_y
        Yc_y = curve.Yc_to_y
        Xc = curve.Xc
        Yc = curve.Yc
        plt.annotate('$Area_y$', xy = (Xc_y,Yc_y),xycoords = 'data')
        plt.plot(Xc_y,Yc_y, marker = 'o', alpha=.4)
        plt.plot([Xc_y,Xc_y],[Yc_y+.2,Yc_y-.2], color = 'black')
        plt.plot([Xc_y+.2,Xc_y-.2],[Yc_y,Yc_y], color = 'black')
        
        plt.annotate(r'$\alpha_B$', xy=(0.,0.), xycoords='data',
                xytext=(0, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate(r'$\alpha_E$', xy=(15.,12.), xycoords='data',
                xytext=(0, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
                
        plt.annotate('$Area$', xy = (Xc,Yc),xycoords = 'data',
                     xytext=(20, -20), textcoords='offset points')
        plt.annotate('$Xc$,$Yc$', xy = (Xc,Yc),xycoords = 'data',
                     xytext=(-50, -20), textcoords='offset points')
        plt.plot(Xc,Yc, marker = 'o', alpha=.4)
        plt.plot([Xc,Xc],[Yc+.2,Yc-.2], color = 'black')
        plt.plot([Xc+.2,Xc-.2],[Yc,Yc], color = 'black')
        #"""
        #
        # Now shift the bottom area   Xc - to 9.5
        #
        #"""
        FPD = FormParameterDict(curve) 
        FPD.add_Area_to_y_Constraint(kind='equality', value = 84.)
        FPD.add_XcConstraint(kind='equality', value = 9.5)
        FPD.add_xVertexConstraint(kind = 'LS', value=4., index =2, weight = 10.)
        FPD.add_yVertexConstraint(kind = 'LS', value=3., index =2, weight = 10.)
        FPD.add_xVertexConstraint(kind = 'LS', value=3.5, index =3, weight = 1000000.)
        FPD.add_yVertexConstraint(kind = 'LS', value=8., index =3, weight = 1000000.)
        FPD.add_xVertexConstraint(kind = 'LS', value=9., index =4, weight = 1000000.)
        FPD.add_yVertexConstraint(kind = 'LS', value=10., index =4, weight = 1000000.)
        FPD.add_xVertexConstraint(kind = 'LS', value=5.5, index =5, weight = 1000000.)
        FPD.add_yVertexConstraint(kind = 'LS', value=4.5, index =5, weight = 1000000.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 20.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 60.)
        FPD.add_E1(kind='LS', weight = .2)
        FPD.add_E2(kind='LS', weight = .2)
        FPD.add_E3(kind='LS', weight = 2.)
        FPD.add_ArcLengthApprox(kind='LS', weight = .2)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp42')
        curve = Lspline.curve
        curve.compute_area_to_y()
        curve.compute_moments_to_y()
        curve.computeXc_to_y()
        curve.computeYc_to_y()
        curve.compute_area()
        curve.compute_moments()
        curve.computeXc()
        curve.computeYc()
        store2 = copy.deepcopy(curve)
        #""" plot it
        curve.plotcurve_detailed()
        #curve.plot_area_to_x_shiftaxis(shift_axis = 12.)
        Xc_y = curve.Xc_to_y
        Yc_y = curve.Yc_to_y
        Xc = curve.Xc
        Yc = curve.Yc
        plt.annotate('$Area_y$', xy = (Xc_y,Yc_y),xycoords = 'data')
        plt.plot(Xc_y,Yc_y, marker = 'o', alpha=.4)
        plt.plot([Xc_y,Xc_y],[Yc_y+.2,Yc_y-.2], color = 'black')
        plt.plot([Xc_y+.2,Xc_y-.2],[Yc_y,Yc_y], color = 'black')
        
        plt.annotate(r'$\alpha_B$', xy=(0.,0.), xycoords='data',
                xytext=(0, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate(r'$\alpha_E$', xy=(15.,12.), xycoords='data',
                xytext=(0, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
                
        plt.annotate('$Area$', xy = (Xc,Yc),xycoords = 'data',
                     xytext=(20, -20), textcoords='offset points')
        plt.annotate('$Xc$,$Yc$', xy = (Xc,Yc),xycoords = 'data',
                     xytext=(-50, -20), textcoords='offset points')
        plt.plot(Xc,Yc, marker = 'o', alpha=.4)
        plt.plot([Xc,Xc],[Yc+.2,Yc-.2], color = 'black')
        plt.plot([Xc+.2,Xc-.2],[Yc,Yc], color = 'black')
        #"""
        
        #
        # Now shift the bottom area   Xc - to 9.5   Yc = 4.0
        #
        #"""
        FPD = FormParameterDict(curve) 
        FPD.add_Area_to_y_Constraint(kind='equality', value = 84.)
        FPD.add_XcConstraint(kind='equality', value = 9.5)
        FPD.add_YcConstraint(kind='equality', value = 4.0)
        FPD.add_xVertexConstraint(kind = 'LS', value=4., index =2, weight = 10.)
        FPD.add_yVertexConstraint(kind = 'LS', value=3., index =2, weight = 10.)
        FPD.add_xVertexConstraint(kind = 'LS', value=3.5, index =3, weight = 1000000.)
        FPD.add_yVertexConstraint(kind = 'LS', value=8., index =3, weight = 1000000.)
        FPD.add_xVertexConstraint(kind = 'LS', value=9., index =4, weight = 1000000.)
        FPD.add_yVertexConstraint(kind = 'LS', value=10., index =4, weight = 1000000.)
        FPD.add_xVertexConstraint(kind = 'LS', value=5.5, index =5, weight = 1000000.)
        FPD.add_yVertexConstraint(kind = 'LS', value=4.5, index =5, weight = 1000000.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 20.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 60.)
        FPD.add_E1(kind='LS', weight = .2)
        FPD.add_E2(kind='LS', weight = .2)
        FPD.add_E3(kind='LS', weight = 2.)
        FPD.add_ArcLengthApprox(kind='LS', weight = .2)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp42')
        curve = Lspline.curve
        curve.compute_area_to_y()
        curve.compute_moments_to_y()
        curve.computeXc_to_y()
        curve.computeYc_to_y()
        curve.compute_area()
        curve.compute_moments()
        curve.computeXc()
        curve.computeYc()
        store3 = copy.deepcopy(curve)
        
        #"""
        curve.plotcurve_detailed()
        #curve.plot_area_to_x_shiftaxis(shift_axis = 12.)
        Xc_y = curve.Xc_to_y
        Yc_y = curve.Yc_to_y
        Xc = curve.Xc
        Yc = curve.Yc
        plt.annotate('$Area_y$', xy = (Xc_y,Yc_y),xycoords = 'data')
        plt.plot(Xc_y,Yc_y, marker = 'o', alpha=.4)
        plt.plot([Xc_y,Xc_y],[Yc_y+.2,Yc_y-.2], color = 'black')
        plt.plot([Xc_y+.2,Xc_y-.2],[Yc_y,Yc_y], color = 'black')
        
        plt.annotate(r'$\alpha_B$', xy=(0.,0.), xycoords='data',
                xytext=(0, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate(r'$\alpha_E$', xy=(15.,12.), xycoords='data',
                xytext=(0, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
                
        plt.annotate('$Area$', xy = (Xc,Yc),xycoords = 'data',
                     xytext=(20, -20), textcoords='offset points')
        plt.annotate('$Xc$,$Yc$', xy = (Xc,Yc),xycoords = 'data',
                     xytext=(-50, -20), textcoords='offset points')
        plt.plot(Xc,Yc, marker = 'o', alpha=.4)
        plt.plot([Xc,Xc],[Yc+.2,Yc-.2], color = 'black')
        plt.plot([Xc+.2,Xc-.2],[Yc,Yc], color = 'black')
        #"""
        tlm = []
        xcrange = [9.6,9.4,9.0,8.79]
        for Xc in xcrange:
            FPD = FormParameterDict(curve) 
            FPD.add_Area_to_y_Constraint(kind='equality', value = 84.)
            FPD.add_XcConstraint(kind='equality', value = Xc)
            #FPD.add_YcConstraint(kind='equality', value = 4.0)
            FPD.add_xVertexConstraint(kind = 'LS', value=4., index =2, weight = 10.)
            FPD.add_yVertexConstraint(kind = 'LS', value=3., index =2, weight = 10.)
            FPD.add_xVertexConstraint(kind = 'LS', value=3.5, index =3, weight = 1000000.)
            FPD.add_yVertexConstraint(kind = 'LS', value=8., index =3, weight = 1000000.)
            FPD.add_xVertexConstraint(kind = 'LS', value=9., index =4, weight = 1000000.)
            FPD.add_yVertexConstraint(kind = 'LS', value=10., index =4, weight = 1000000.)
            FPD.add_xVertexConstraint(kind = 'LS', value=5.5, index =5, weight = 1000000.)
            FPD.add_yVertexConstraint(kind = 'LS', value=4.5, index =5, weight = 1000000.)
            FPD.add_AngleConstraint(kind='equality', location = 0., value = 20.)
            FPD.add_AngleConstraint(kind='equality', location = 1., value = 60.)
            FPD.add_E1(kind='LS', weight = .2)
            FPD.add_E2(kind='LS', weight = .2)
            FPD.add_E3(kind='LS', weight = 2.)
            FPD.add_ArcLengthApprox(kind='LS', weight = .2)
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            Lspline.optimize(stop = 50)
            curve.compute_area_to_y()
            curve.compute_moments_to_y()
            curve.computeXc_to_y()
            curve.computeYc_to_y()
            curve.compute_area()
            curve.compute_moments()
            curve.computeXc()
            curve.computeYc()
            tlm.append(copy.deepcopy(Lspline.curve))
        #"""
        
        for cv in tlm:
            plt.plot(cv.r[:,0],cv.r[:,1], color='black')
            Xc_y = cv.Xc_to_y
            Yc_y = cv.Yc_to_y
            Xc = cv.Xc
            Yc = cv.Yc
            #plt.annotate('$Area$', xy = (Xc,Yc),xycoords = 'data')
            
            plt.plot(Xc,Yc, marker = 'o', alpha=.4, color = 'black')
            plt.plot([Xc,Xc],[Yc+.2,Yc-.2], color = 'black')
            plt.plot([Xc+.2,Xc-.2],[Yc,Yc], color = 'black')
        plt.annotate(r'$\alpha_B$', xy=(0.,0.), xycoords='data',
                xytext=(0, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate(r'$\alpha_E$', xy=(15.,12.), xycoords='data',
                xytext=(0, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
                
        plt.annotate(' $Area$', xy = (6.,1.5),xycoords = 'data',
                     xytext=(20, -20), textcoords='offset points')
        plt.annotate('Controlled Variation in $Xc$', xy = (6.,2.),xycoords = 'data')
        plt.plot([9.6,9.6],[3.,3.7], color = 'black')
        plt.plot([8.79,8.79],[3.,3.7], color = 'black')
        plt.annotate('', 
                     xy=(9.6,3.), 
                     xycoords='data',
                     xytext=(8.79,3.), textcoords='data',
                     arrowprops=dict(arrowstyle="<->",
                                      connectionstyle="bar",
                                      ec="k",
                                      shrinkA=5, shrinkB=5,
                                      ))   

    if do_all:
        test = 'exp4'
    if test == 'exp4':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 0.
        xe = 15.
        ye = 12.
        alphab = 0.#-5.
        alphae = 0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 80.#58.#84.#
        slope = 'down'
        
        ae = alphae
        ab = alphab
        
        curve = initial_curve([xb,yb],[xe,ye],8) 
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        #FPD.add_AreaConstraint(kind='equality', value = 96)
        FPD.add_xVertexConstraint(kind = 'equality', value=4., index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=3., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=8., index =3)
        FPD.add_xVertexConstraint(kind = 'equality', value=9., index =4)
        FPD.add_yVertexConstraint(kind = 'equality', value=9.5, index =4)
        FPD.add_xVertexConstraint(kind = 'equality', value=6.5, index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=4.5, index =5)
        
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 20.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 60.)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp4')
        curve = Lspline.curve
        
        #"""
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 96)
        FPD.add_xVertexConstraint(kind = 'equality', value=4., index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=3., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=8., index =3)
        FPD.add_xVertexConstraint(kind = 'equality', value=9., index =4)
        FPD.add_yVertexConstraint(kind = 'equality', value=9.5, index =4)
        FPD.add_xVertexConstraint(kind = 'equality', value=6.5, index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=4.5, index =5)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 20.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 60.)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp41')
        curve = Lspline.curve
        #"""
        #"""
        FPD = FormParameterDict(curve) 
        FPD.add_Area_to_y_Constraint(kind='equality', value = 84.)
        FPD.add_xVertexConstraint(kind = 'equality', value=4., index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=3., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=8., index =3)
        FPD.add_xVertexConstraint(kind = 'equality', value=9., index =4)
        FPD.add_yVertexConstraint(kind = 'equality', value=9.5, index =4)
        FPD.add_xVertexConstraint(kind = 'equality', value=6.5, index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=4.5, index =5)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 20.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 60.)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp42')
        curve = Lspline.curve
        #"""
        
        ##
        ## Best fwd curve in green imediately below
        ##
        """Now MOVE this curve over 1 or 2 units right in x
        and increase the y - area accordingly
        """
        #"""
        ### Best fwd curve to date!
        verts = curve.vertices
        verts = verts + [2.,2.]
        curve.vertices = verts
        FPD = FormParameterDict(curve) 
        #FPD.add_AreaConstraint(kind='equality', value = 126.)
        FPD.add_Area_to_y_Constraint(kind='equality', value = 108.)
        FPD.add_xVertexConstraint(kind = 'equality', value=6., index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=5., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=5.5, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        FPD.add_xVertexConstraint(kind = 'equality', value=11., index =4)
        FPD.add_yVertexConstraint(kind = 'equality', value=11.5, index =4)
        FPD.add_xVertexConstraint(kind = 'equality', value=8.5, index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=6.5, index =5)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 20.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = 60.)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp431')
        curve = Lspline.curve
        #"""
        """
        curve.plotcurve_detailed()
        curve.plot_area_to_x_shiftaxis(shift_axis = 14.)
        curve.compute_area_to_y()
        curve.compute_moments_to_y()
        curve.computeXc_to_y()
        curve.computeYc_to_y()
        Xc_y = curve.Xc_to_y
        Yc_y = curve.Yc_to_y
        plt.annotate('$Area$', xy = (Xc_y,Yc_y),xycoords = 'data')
        plt.annotate('$tangent_b$', xy=(2.,2.), xycoords='data',
                xytext=(0, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate('$tangent_b$', xy=(17.,14.), xycoords='data',
                xytext=(0, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
                
        #"""
    if do_all:
        test = 'exp5'
    if test == 'exp5':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        print '\n\n\n\n CURVE 1'
        k=4
        nump=30
        xb = 0.
        yb = 12.
        xe = 15.
        ye = 0.
        alphab = 0.#-5.
        alphae = 0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 80.#58.#84.#
        slope = 'down'
        
        ae = alphae
        ab = alphab
        
        curve = initial_curve([xb,yb],[xe,ye],12) 
        #-------------------------------------------------------------------
        #"""        
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        
        #FPD.add_xPointConstraint(kind = 'equality', value= .0, location=.3)
        #FPD.add_yPointConstraint(kind = 'equality', value=10., location=.3)
        
        FPD.add_xVertexConstraint(kind = 'equality', value=8.0, index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=8.,  index =5)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =6)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =6)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =7)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =7)

        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .88, value = 0.)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        
        
        Lspline.curve.pts_M_pts()
        Lspline.curve.compute_arclength()
        print 'initial Initial run curve arc length = ',Lspline.curve.AL
        save_ini_al = copy.deepcopy(Lspline.curve.AL)
        
        Lspline.compute_lagrangian()
        save_ini_curve = copy.deepcopy(Lspline)
        keys = Lspline.Lagrangian.obj.keys()
        #Lspline.Lagrangian.obj
        ini_ans = {}
        for key in keys:
            print Lspline.Lagrangian.obj[key].type, key, ' : ',Lspline.Lagrangian.obj[key].computed_value
            ini_ans[Lspline.Lagrangian.obj[key].type] = Lspline.Lagrangian.obj[key].computed_value
        
        Lspline.curve.compute_arclength()
        print 'initial run curve arc length = ',Lspline.curve.AL
        
        
        
        
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp50')
        curve1 = copy.deepcopy(Lspline.curve)
        #curve1 = Lspline.curve
        save1 = copy.deepcopy(Lspline.curve)
        
        
            
        print 'Completed Initial run curve Curve Solve'
        keys = Lspline.Lagrangian.obj.keys()
        #Lspline.Lagrangian.obj
        fini_ans = {}
        for key in keys:
            print Lspline.Lagrangian.obj[key].type, key, ' : ',Lspline.Lagrangian.obj[key].computed_value
            fini_ans[Lspline.Lagrangian.obj[key].type] = Lspline.Lagrangian.obj[key].computed_value
        
        anskeys = fini_ans.keys()
        for key in anskeys:
           print key, 100.*( fini_ans[key].value   - ini_ans[key].value  ) / ini_ans[key].value 
        for key in anskeys:
           print key, 100.*( ini_ans[key].value  - fini_ans[key].value  ) / ini_ans[key].value 
        Lspline.curve.compute_arclength()
        print 'final Initial run curve arc length = ',Lspline.curve.AL
        
        
        
        """
        #
        #
        #
        #
        #-------------------------------------------------------------------
        #
        #
        #
        #
        #
        #"""
        print '\n\n\n\n CURVE 2'
        curve = copy.deepcopy(curve1)        
        curve.knot_insertion(.58)
        curve.knot_insertion(.6)
        #curve.knot_insertion(.62)
        
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        FPD.add_xVertexConstraint(kind = 'equality', value=8.0, index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=8.,  index =5)
        
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =6)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =6)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =7)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =7)
        
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =8)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =8)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =9)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =9)
    
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .88, value = 0.)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        
        
        
        
        Lspline.curve.pts_M_pts()
        Lspline.compute_lagrangian()
        keys = Lspline.Lagrangian.obj.keys()
        ini_ans2 = {}
        for key in keys:
            print Lspline.Lagrangian.obj[key].type, key, ' : ',Lspline.Lagrangian.obj[key].computed_value
            ini_ans2[Lspline.Lagrangian.obj[key].type] = Lspline.Lagrangian.obj[key].computed_value
        Lspline.curve.pts_M_pts()
        Lspline.curve.compute_arclength()
        print 'ini run curve 2 arc length = ',Lspline.curve.AL
        
        
        
        
        
        
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp51')
        curve2 = copy.deepcopy(Lspline.curve)
        save2 = copy.deepcopy(Lspline.curve)
        #curve2 = Lspline.curve
        
        
        
        print 'Completed FLAT SAC Curve Solve 1 (total curves solve 2)'
        print 'with Tangent flat endpoints'
        keys = Lspline.Lagrangian.obj.keys()
        fini_ans2 = {}
        for key in keys:
            print Lspline.Lagrangian.obj[key].type, key, ' : ',Lspline.Lagrangian.obj[key].computed_value
            fini_ans2[Lspline.Lagrangian.obj[key].type] = Lspline.Lagrangian.obj[key].computed_value
        
        anskeys = fini_ans2.keys()
        for key in anskeys:
           print key, 100.*( fini_ans2[key].value  - ini_ans[key].value  ) / ini_ans[key].value 
        for key in anskeys:
           print key, 100.*( ini_ans[key].value  - fini_ans2[key].value  ) / ini_ans[key].value 
        Lspline.curve.compute_arclength()
        print 'final run curve 2 arc length = ',Lspline.curve.AL
        
        """
        #
        #
        #
        #
        #-------------------------------------------------------------------
        #
        #
        #
        #
        #
        #"""
        print '\n\n\n\n CURVE 3'
        curve = copy.deepcopy(curve2)      
        curve.knot_insertion(.62)
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        #FPD.add_AreaConstraint(kind='equality', value = 60.)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        
        
        FPD.add_xVertexConstraint(kind = 'equality', value=8.0, index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=8.,  index =5)
        
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =6)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =6)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =7)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =7)
         
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =8)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =8)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =9)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =9)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =10)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =10)
        
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .88, value = 0.)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        
        
        Lspline.curve.pts_M_pts()
        Lspline.compute_lagrangian()
        keys = Lspline.Lagrangian.obj.keys()
        ini_ans3 = {}
        for key in keys:
            print Lspline.Lagrangian.obj[key].type, key, ' : ',Lspline.Lagrangian.obj[key].computed_value
            ini_ans3[Lspline.Lagrangian.obj[key].type] = Lspline.Lagrangian.obj[key].computed_value
        Lspline.curve.pts_M_pts()
        Lspline.curve.compute_arclength()
        print 'Initial Run curve 3 arc length = ',Lspline.curve.AL
        
        
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp52')
        curve3 = Lspline.curve
        save3 = copy.deepcopy(Lspline.curve) 
            
        print 'Completed Run curve 3 Solve'
        fini_ans3 = {}
        for key in keys:
            print Lspline.Lagrangian.obj[key].type, key, ' : ',Lspline.Lagrangian.obj[key].computed_value
            fini_ans3[Lspline.Lagrangian.obj[key].type] = Lspline.Lagrangian.obj[key].computed_value
        
        anskeys = fini_ans3.keys()
        for key in anskeys:
           print key, 100.*( fini_ans3[key].value   - ini_ans[key].value  ) / ini_ans[key].value 
        for key in anskeys:
           print key, 100.*( ini_ans[key].value  - fini_ans3[key].value  ) / ini_ans[key].value 
        Lspline.curve.compute_arclength()
        print 'final Run curve 3  arc length = ',Lspline.curve.AL
        
        print 'check percents: '
        print 'compute_percent_decrease','(3.696e3  & 2.721e3)'
        
        
        
        """
        #
        #
        #
        #
        #-------------------------------------------------------------------
        #
        #
        #
        #
        #
        #"""
        print '\n\n\n\n CURVE 4'
        curve = copy.deepcopy(curve3)      
        curve.knot_insertion(.59)
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        #FPD.add_AreaConstraint(kind='equality', value = 60.)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        
        
        FPD.add_xVertexConstraint(kind = 'equality', value=8.0, index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=8.,  index =5)
        
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =6)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =6)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =7)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =7)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =8)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =8)
         
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =9)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =9)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =10)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =10)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =11)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =11)
        
        #FPD.add_AngleConstraint(kind='equality', location = 0., value = -20.)
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .88, value = 0.)
        
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        
        
            
        Lspline.curve.pts_M_pts()
        Lspline.compute_lagrangian()
        keys = Lspline.Lagrangian.obj.keys()
        ini_ans4 = {}
        for key in keys:
            print Lspline.Lagrangian.obj[key].type, key, ' : ',Lspline.Lagrangian.obj[key].computed_value
            ini_ans4[Lspline.Lagrangian.obj[key].type] = Lspline.Lagrangian.obj[key].computed_value
        Lspline.curve.pts_M_pts()
        Lspline.curve.compute_arclength()
        print 'initial Run curve 4 arc length = ',Lspline.curve.AL
        
        
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp53')
        curve4 = Lspline.curve
            
        print 'Completed Run curve 4 Solve'
        fini_ans4 = {}
        for key in keys:
            print Lspline.Lagrangian.obj[key].type, key, ' : ',Lspline.Lagrangian.obj[key].computed_value
            fini_ans4[Lspline.Lagrangian.obj[key].type] = Lspline.Lagrangian.obj[key].computed_value
        
        anskeys = fini_ans4.keys()
        for key in anskeys:
           print key, 100.*( fini_ans4[key].value   - ini_ans[key].value  ) / ini_ans[key].value 
        for key in anskeys:
           print key, 100.*( ini_ans[key].value  - fini_ans4[key].value  ) / ini_ans[key].value 
        Lspline.curve.compute_arclength()
        print 'final Run curve 4 arc length = ',Lspline.curve.AL
        
        
        """
        #
        #
        #
        #
        #-------------------------------------------------------------------
        #
        #
        #
        #
        #
        #"""
        if False:
            curve = copy.deepcopy(curve4)      
            #curve.knot_insertion(.59)
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve) 
            FPD.add_AreaConstraint(kind='equality', value = 60.)
            #FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =1)
            #FPD.add_yVertexConstraint(kind = 'equality', value=10., index =1)
            FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
            FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
            FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
            FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
            
            
            FPD.add_xVertexConstraint(kind = 'equality', value=8.0, index =5)
            FPD.add_yVertexConstraint(kind = 'equality', value=8.,  index =5)
            
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =6)
            FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =6)
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =7)
            FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =7)
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =8)
            FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =8)
             
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =9)
            FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =9)
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =10)
            FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =10)
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =11)
            FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =11)
            
            #FPD.add_AngleConstraint(kind='equality', location = 0., value = -20.)
            FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
            FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
            FPD.add_CurvatureConstraint(kind='equality', location = .88, value = 0.)
            
            FPD.add_E1(kind='LS', weight = 1.)
            FPD.add_E2(kind='LS', weight = 1.)
            FPD.add_E3(kind='LS', weight = 1.)
            FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            Lspline.optimize(stop = 50)
            #Lspline.display(mytitle = 'exp54')
            curve5 = Lspline.curve
        """
        #
        #
        #
        #
        #-------------------------------------------------------------------
        #
        #
        #
        #
        #
        #"""
        if False:
            curve = copy.deepcopy(curve5)      
            curve.knot_insertion(.1)
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve) 
            FPD.add_AreaConstraint(kind='LS', value = 60.)
            FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
            FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
            FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
            FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
            FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =4)
            FPD.add_yVertexConstraint(kind = 'equality', value=10., index =4)
            
            
            FPD.add_xVertexConstraint(kind = 'equality', value=8.0, index =6)
            FPD.add_yVertexConstraint(kind = 'equality', value=8.,  index =6)
            
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =7)
            FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =7)
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =8)
            FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =8)
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =9)
            FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =9)
             
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =10)
            FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =10)
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =11)
            FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =11)
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =12)
            FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =12)
            
            FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
            FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
            FPD.add_CurvatureConstraint(kind='equality', location = .88, value = 0.)
            
            FPD.add_E1(kind='LS', weight = 1.)
            FPD.add_E2(kind='LS', weight = 1.)
            FPD.add_E3(kind='LS', weight = 1.)
            FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            Lspline.optimize(stop = 50)
            #Lspline.display(mytitle = 'exp55')
            curve6 = Lspline.curve
        """
        #
        #
        #
        #
        #-------------------------------------------------------------------
        #
        #
        #
        #
        #
        #"""
        if False:
            curve = copy.deepcopy(curve6)   
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve) 
            FPD.add_AreaConstraint(kind='equality', value = 60.)
            
            FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
            FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
            FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
            FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
            FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =4)
            FPD.add_yVertexConstraint(kind = 'equality', value=10., index =4)
            
            #FPD.add_xPointConstraint(kind = 'equality', value=.0,  location=.3)
            #FPD.add_yPointConstraint(kind = 'equality', value=10., location=.3)
            
            FPD.add_yVertexConstraint(kind = 'equality', value=10., index =5)
            
            
            FPD.add_xVertexConstraint(kind = 'equality', value=8.0, index =6)
            FPD.add_yVertexConstraint(kind = 'equality', value=8.,  index =6)
            
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =7)
            FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =7)
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =8)
            FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =8)
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =9)
            FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =9)
             
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =10)
            FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =10)
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =11)
            FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =11)
            FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =12)
            FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =12)
            
            FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
            FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
            FPD.add_CurvatureConstraint(kind='equality', location = .88, value = 0.)
            
            FPD.add_E1(kind='LS', weight = 1.)
            FPD.add_E2(kind='LS', weight = 1.)
            FPD.add_E3(kind='LS', weight = 1.)
            FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            Lspline.optimize(stop = 50)
            #Lspline.display(mytitle = 'exp56')
            curve7 = Lspline.curve
        """
        #
        #
        #
        #
        #-------------------------------------------------------------------
        #
        #
        #
        #
        #
        #"""
        if False:
            plt.close()
            
            shift_axis = 12.
            el = Ellipse((2, -1), 0.5, 0.5)
            color = 'green'
            curve.plotcurve_detailed()
            plt.fill_between(self.r[:,0],shift_axis,self.r[:,1],
                                facecolor = color, alpha=.1,
                                label='area')
            plt.annotate('$Area$', xy = (10.,6.),xycoords = 'data')
            
            plt.annotate('$tangent_b$', xy=(0.,12.), xycoords='data',
                    xytext=(-30, 20), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->")
                    )
                    
            plt.annotate('$knuckle_0$', xy=(0.,10.), xycoords='data',
                    xytext=(20, 30), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->")
                    )
            plt.annotate('$tangent_e$', xy=(15.,0.), xycoords='data',
                    xytext=(-50, 30), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->")
                    )
            plt.annotate('$knuckle_1$', xy=(3.5,5.), xycoords='data',
                    xytext=(20, 0), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->")
                    )
            plt.annotate('$knuckle_2$', xy=(3.5,3.), xycoords='data',
                    xytext=(20, 0), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->")
                    )
            plt.close()
    
    if do_all:
        test = 'container'
    if test == 'container':
        pass
    
    if do_all:
        test = 'exp6'
    if test == 'exp6':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 12.
        xe = 12.
        ye = 0.
        alphab = 0.
        alphae = 0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 72.#58.#84.#
        slope = 'up'
        
        ae = alphae
        ab = alphab
        
        curve = initial_curve([xb,yb],[xe,ye],12) 
        #-------------------------------------------------------------------
        #"""        
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 

        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .0, value = 0.)
        
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp50')
        curve1 = Lspline.curve
        #"""
        #-------------------------------------------------------------------
        #"""
        curve = copy.deepcopy(curve1)   
        
    if do_all:
        test = 'dist_contraint'
    if test == 'dist_contraint':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 12.
        xe = 12.
        ye = 0.
        alphab = 0.
        alphae = 0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 72.
        slope = 'up'
        
        ae = alphae
        ab = alphab
        
        curve = initial_curve([xb,yb],[xe,ye],7) 
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 

        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 0., value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_AreaConstraint(kind = 'equality', value = 72.)
        """
        FPD.add_xVertexConstraint(kind = 'equality', 
                                  location=None, 
                                  value = 4., 
                                  weight = 1.0, 
                                  index = 4 )
        #"""
        #"""
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                       value = 1., weight = 1.0, 
                                       index = 3, index2 = 4, 
                                       seperation = 1. )
        FPD.add_relative_yVertexConstraint(kind = 'equality', location=None, 
                                       value = 7., weight = 1.0, 
                                       index = 3, index2 = 4, 
                                       seperation = 7. )
        #"""
        #"""
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                       value = -1.2, weight = 1.0, 
                                       index = 2, index2 = 3, 
                                       seperation = -1.2 )
        FPD.add_relative_yVertexConstraint(kind = 'equality', location=None, 
                                       value = 1., weight = 1.0, 
                                       index = 2, index2 = 3, 
                                       seperation = 1. )
        #"""
        FPD.add_E1(kind='LS', weight = .01)
        FPD.add_E2(kind='LS', weight = .01)
        FPD.add_E3(kind='LS', weight = .001)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 25)
        
        Lspline.curve.plotcurve_detailed()
        Lspline.curve.plot_area_to_x()
        Lspline.curve.compute_moments()
        Lspline.curve.computeCentroid()
        plt.annotate('$Area$', 
                     xy = (float(Lspline.curve.Xc.value),
                           float(Lspline.curve.Yc.value)),
                     xycoords = 'data')
        plt.annotate('$tangent_b$', xy=(0.,12.), xycoords='data',
                xytext=(-30, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate('$tangent_e$', xy=(12.,0.), xycoords='data',
                xytext=(-50, 30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.plot(float(Lspline.curve.Xc.value),
                 float(Lspline.curve.Yc.value),marker = 'o',
                 color = 'black', alpha=.4)
                 
        """
        plt.annotate('', 
                     xy=(Lspline.curve.vertices[3,0],
                         Lspline.curve.vertices[2,1]), 
                     xycoords='data',
                     xytext=(Lspline.curve.vertices[3,0],
                             Lspline.curve.vertices[3,1]), textcoords='data',
                     arrowprops=dict(arrowstyle="<->",
                                      connectionstyle="bar",
                                      ec="k",
                                      shrinkA=5, shrinkB=5,
                                      ))
        plt.annotate('', 
                     xy=(Lspline.curve.vertices[4,0],
                         Lspline.curve.vertices[4,1]), 
                     xycoords='data',
                     xytext=(Lspline.curve.vertices[3,0],
                             Lspline.curve.vertices[4,1]), 
                     textcoords='data',
                     arrowprops=dict(arrowstyle="<->",
                                      connectionstyle="bar",
                                      ec="k",
                                      shrinkA=5, shrinkB=5,
                                      ))
        #"""
        plt.annotate('', 
                     xy=Lspline.curve.vertices[2], 
                     xycoords='data',
                     xytext=Lspline.curve.vertices[3], textcoords='data',
                     arrowprops=dict(arrowstyle="<->",
                                      connectionstyle="bar",
                                      ec="k",
                                      shrinkA=5, shrinkB=5,
                                      ))
        plt.annotate('', 
                     xy=Lspline.curve.vertices[3], 
                     xycoords='data',
                     xytext=(Lspline.curve.vertices[4]), 
                     textcoords='data',
                     arrowprops=dict(arrowstyle="<->",
                                      connectionstyle="bar",
                                      ec="k",
                                      shrinkA=5, shrinkB=5,
                                      ))
                                      
    
    if do_all:
        test = 'exp7'
    if test == 'exp7':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 12.
        xe = 15.
        ye = 0.
        alphab = 0.#-5.
        alphae = 0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 80.#58.#84.#
        slope = 'down'
        
        ae = alphae
        ab = alphab
        
        curve = initial_curve([xb,yb],[xe,ye],12) 
        #-------------------------------------------------------------------
        #"""        
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        
        #FPD.add_xPointConstraint(kind = 'equality', value= .0, location=.3)
        #FPD.add_yPointConstraint(kind = 'equality', value=10., location=.3)
        
        FPD.add_xVertexConstraint(kind = 'equality', value=8.0, index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=8.,  index =5)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =6)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =6)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =7)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =7)

        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .88, value = 0.)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'exp50')
        curve1 = Lspline.curve
        #"""
        #-------------------------------------------------------------------
        #"""
        curve = copy.deepcopy(curve1)        
        curve.knot_insertion(.58)
        curve.knot_insertion(.6)
        #curve.knot_insertion(.62)
        
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        FPD.add_xVertexConstraint(kind = 'equality', value=8.0, index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=8.,  index =5)
        
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =6)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =6)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =7)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =7)
        
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =8)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =8)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =9)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =9)
    
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .88, value = 0.)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'exp51')
        curve2 = Lspline.curve
        #"""
        #-------------------------------------------------------------------
        #"""
        curve = copy.deepcopy(curve2)      
        curve.knot_insertion(.62)
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        #FPD.add_AreaConstraint(kind='equality', value = 60.)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        
        
        FPD.add_xVertexConstraint(kind = 'equality', value=8.0, index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=8.,  index =5)
        
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =6)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =6)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =7)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =7)
         
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =8)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =8)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =9)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =9)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =10)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =10)
        
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .88, value = 0.)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'exp52')
        curve3 = Lspline.curve
        #"""
        #-------------------------------------------------------------------
        #"""
        curve = copy.deepcopy(curve3)      
        curve.knot_insertion(.59)
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        #FPD.add_AreaConstraint(kind='equality', value = 60.)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        
        
        FPD.add_xVertexConstraint(kind = 'equality', value=8.0, index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=8.,  index =5)
        
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =6)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =6)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =7)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =7)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =8)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =8)
         
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =9)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =9)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =10)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =10)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =11)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =11)
        
        #FPD.add_AngleConstraint(kind='equality', location = 0., value = -20.)
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .88, value = 0.)
        
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'exp53')
        curve4 = Lspline.curve
        #"""
        #-------------------------------------------------------------------
        #"""
        curve = copy.deepcopy(curve4)      
        #curve.knot_insertion(.59)
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 60.)
        #FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =1)
        #FPD.add_yVertexConstraint(kind = 'equality', value=10., index =1)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        
        
        FPD.add_xVertexConstraint(kind = 'equality', value=8.0, index =5)
        FPD.add_yVertexConstraint(kind = 'equality', value=8.,  index =5)
        
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =6)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =6)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =7)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =7)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =8)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =8)
         
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =9)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =9)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =10)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =10)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =11)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =11)
        
        #FPD.add_AngleConstraint(kind='equality', location = 0., value = -20.)
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .88, value = 0.)
        
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'exp54')
        curve5 = Lspline.curve
        #"""
        #-------------------------------------------------------------------
        #"""
        curve = copy.deepcopy(curve5)      
        curve.knot_insertion(.1)
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='LS', value = 60.)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =4)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =4)
        
        
        FPD.add_xVertexConstraint(kind = 'equality', value=8.0, index =6)
        FPD.add_yVertexConstraint(kind = 'equality', value=8.,  index =6)
        
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =7)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =7)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =8)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =8)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =9)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =9)
         
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =10)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =10)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =11)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =11)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =12)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =12)
        
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .88, value = 0.)
        
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'exp55')
        curve6 = Lspline.curve
        #"""
        #-------------------------------------------------------------------
        #"""
        curve = copy.deepcopy(curve6)   
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 60.)
        """
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =4)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =4)
        #"""
        #"""
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 0, index2 = 2, 
                                           seperation = 0. )
        FPD.add_relative_yVertexConstraint(kind = 'equality', location=None, 
                                           value = 2., weight = 1.0, 
                                           index = 0, index2 = 2, 
                                           seperation = 2. )
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 0, index2 = 3, 
                                           seperation = 0. )
        FPD.add_relative_yVertexConstraint(kind = 'equality', location=None, 
                                           value = 2., weight = 1.0, 
                                           index = 0, index2 = 3, 
                                           seperation = 2. )
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 0, index2 = 4, 
                                           seperation = 0. )
        FPD.add_relative_yVertexConstraint(kind = 'equality', location=None, 
                                           value = 2., weight = 1.0, 
                                           index = 0, index2 = 4, 
                                           seperation = 2. )
        #"""
        
        #FPD.add_xPointConstraint(kind = 'equality', value=.0,  location=.3)
        #FPD.add_yPointConstraint(kind = 'equality', value=10., location=.3)
        
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =5)
        
        
        FPD.add_xVertexConstraint(kind = 'equality', value=8.0, index =6)
        FPD.add_yVertexConstraint(kind = 'equality', value=8.,  index =6)
        
        #"""
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =7)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =7)
        #"""
        """
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =8)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =8)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =9)
        FPD.add_yVertexConstraint(kind = 'equality', value=5.,  index =9)
        #"""
        #"""
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 7, index2 = 8, 
                                           seperation = 0. )
        FPD.add_relative_yVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 7, index2 = 8, 
                                           seperation = 0. )
        #"""
        #"""
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 7, index2 = 9, 
                                           seperation = 0. )
        FPD.add_relative_yVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 7, index2 = 9, 
                                           seperation = 0. )
        #"""
        """
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =10)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =10)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =11)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =11)
        FPD.add_xVertexConstraint(kind = 'equality', value=3.5, index =12)
        FPD.add_yVertexConstraint(kind = 'equality', value=3.0, index =12)
        #"""
        #"""
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 7, index2 = 10, 
                                           seperation = 0. )
        FPD.add_relative_yVertexConstraint(kind = 'equality', location=None, 
                                           value = 2., weight = 1.0, 
                                           index = 7, index2 = 10, 
                                           seperation = 2. )
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 7, index2 = 11, 
                                           seperation = 0. )
        FPD.add_relative_yVertexConstraint(kind = 'equality', location=None, 
                                           value = 2., weight = 1.0, 
                                           index = 7, index2 = 11, 
                                           seperation = 2. )
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 7, index2 = 12, 
                                           seperation = 0. )
        FPD.add_relative_yVertexConstraint(kind = 'equality', location=None, 
                                           value = 2., weight = 1.0, 
                                           index = 7, index2 = 12, 
                                           seperation = 2. )
        #"""
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        FPD.add_CurvatureConstraint(kind='equality', location = .88, value = 0.)
        
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #Lspline.display(mytitle = 'exp56')
        curve7 = Lspline.curve
        #"""
        plt.close()
        
        shift_axis = 12.
        el = Ellipse((2, -1), 0.5, 0.5)
        color = 'green'
        plt.close()
        curve7.plotCurvature(scale = 3., color_='grey')
        curve7.plotcurve_detailed()
        plt.fill_between(curve7.r[:,0],shift_axis,curve7.r[:,1],
                            facecolor = color, alpha=.1,
                            label='area')
        plt.annotate('$Area$', xy = (10.,6.),xycoords = 'data')
        
        plt.annotate('$tangent_b$', xy=(0.,12.), xycoords='data',
                xytext=(-30, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        
        plt.annotate('$curvature_0$', xy=curve7.CurvePoint(.88), xycoords='data',
                xytext=(-30, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
                
        plt.annotate('$knuckle_0$', xy=(0.,10.), xycoords='data',
                xytext=(20, 30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate('$tangent_e$', xy=(15.,0.), xycoords='data',
                xytext=(-50, 30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate('$knuckle_1$', xy=(3.5,5.), xycoords='data',
                xytext=(20, 0), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate('$knuckle_2$', xy=(3.5,3.), xycoords='data',
                xytext=(20, 0), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate('', 
                     xy=Lspline.curve.vertices[10], 
                     xycoords='data',
                     xytext=Lspline.curve.vertices[7], textcoords='data',
                     arrowprops=dict(arrowstyle="<->",
                                      connectionstyle="bar",
                                      ec="k",
                                      shrinkA=5, shrinkB=5,
                                      ))
        plt.annotate('$Vertical$ $Offset$', xy=(2.8,4.), xycoords='data',
                xytext=(-70, 40), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        #plt.close()
                
                
    if do_all:
        test = 'exp8'
    if test == 'exp8':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 12.
        xe = 12.
        ye = 0.
        alphab = 0.
        alphae = 0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 82.
        slope = 'up'
        
        x_offset1 = 3.
        x_offset2 = -2.
        y_offset = 2.
        
        ae = alphae
        ab = alphab
        
        curve = initial_curve([xb,yb],[xe,ye],10) 
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 

        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        FPD.add_CurvatureConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        #FPD.add_AreaConstraint(kind = 'equality', value = curve_area)
        FPD.add_AreaToAnyXAxisConstraint(kind = 'equality', 
                                         value = curve_area, 
                                         x_axis_loc = 12.)

        #"""
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                       value = x_offset1, weight = 1.0, 
                                       index = 3, index2 = 4, 
                                       seperation = x_offset1 )
        FPD.add_relative_yVertexConstraint(kind = 'equality', location=None, 
                                       value = 0., weight = 1.0, 
                                       index = 3, index2 = 4, 
                                       seperation = 0. )
        #"""
        #"""
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                       value = 0., weight = 1.0, 
                                       index = 4, index2 = 5, 
                                       seperation = 0. )
        FPD.add_relative_yVertexConstraint(kind = 'equality', location=None, 
                                       value = y_offset, weight = 1.0, 
                                       index = 4, index2 = 5, 
                                       seperation = y_offset )
        #"""
        #"""
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                       value = x_offset2, weight = 1.0, 
                                       index = 5, index2 = 6, 
                                       seperation = x_offset2 )
        FPD.add_relative_yVertexConstraint(kind = 'equality', location=None, 
                                       value = 0., weight = 1.0, 
                                       index = 5, index2 = 6, 
                                       seperation = 0. )
        #"""
        FPD.add_E1(kind='LS', weight = .01)
        FPD.add_E2(kind='LS', weight = .01)
        FPD.add_E3(kind='LS', weight = .001)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 25)
        
        
        curve1 = copy.deepcopy(Lspline.curve)
        
        curve1.plotcurve_detailed()
        curve1.plot_area_to_x_shiftaxis(shift_axis=12.)
        curve1.plotCurvature()
        plt.annotate('$Area$', xy = (9.,6.),xycoords = 'data')
        
        plt.annotate('$tangent_b$', xy=(0.,12.), xycoords='data',
                xytext=(-30, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        
        plt.annotate('$curvature_b$', xy=curve1.CurvePoint(.0), xycoords='data',
                xytext=(30, 30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
                
        plt.annotate('$tangent_e$', xy=(12.,0.), xycoords='data',
                xytext=(-30, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
                
        plt.annotate('$curvature_e$', xy=curve1.CurvePoint(1.0), xycoords='data',
                xytext=(30, 30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate('', 
                     xy=Lspline.curve.vertices[4], 
                     xycoords='data',
                     xytext=Lspline.curve.vertices[3], textcoords='data',
                     arrowprops=dict(arrowstyle="<->",
                                      connectionstyle="bar",
                                      ec="k",
                                      shrinkA=5, shrinkB=5,
                                      ))
        plt.annotate('', 
                     xy=Lspline.curve.vertices[5], 
                     xycoords='data',
                     xytext=Lspline.curve.vertices[4], textcoords='data',
                     arrowprops=dict(arrowstyle="<->",
                                      connectionstyle="bar",
                                      ec="k",
                                      shrinkA=5, shrinkB=5,
                                      ))
        plt.annotate('', 
                     xy=Lspline.curve.vertices[6], 
                     xycoords='data',
                     xytext=Lspline.curve.vertices[5], textcoords='data',
                     arrowprops=dict(arrowstyle="<->",
                                      connectionstyle="bar",
                                      ec="k",
                                      shrinkA=5, shrinkB=5,
                                      ))   
        plt.annotate('$relation$ $constraint$', xy=(4.0,7.88), xycoords='data',
                xytext=(-70, 50), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate('$relation$ $constraint$', xy=(1.88,6.0), xycoords='data',
                xytext=(-70, 50), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate('$relation$ $constraint$', xy=(3.5,4.4), xycoords='data',
                xytext=(-70, -50), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )