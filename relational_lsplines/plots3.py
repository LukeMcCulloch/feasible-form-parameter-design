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
    nump_=nump
    curve = spline.Bspline(vertices, k, nump_)
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
    test    = None
    test    = 'good_run'  # this one makes the run plot with free stern tube
    #test    = 'exp7'
    #test    = 'exp8'
    do_all = False
    
    if do_all:
        test = 'good_run'
    if test == 'good_run':
        print 'test = {}'.format(test)
        print 'option = {}'.format(option)
        cl = []
        k=4
        nump=100
        xb = 0.
        yb = 12.
        xe = 15.
        ye = 0.
        alphab = 0.
        alphae = 0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 80.
        slope = 'down'
        x_offset1 = 3.
        x_offset2 = -2.
        y_offset = 2.
        ae = alphae
        ab = alphab
        
        curve = initial_curve([xb,yb],[xe,ye],16,nump=100) 
        interval_data, small = interval_bounds(curve)
        
        FPD = FormParameterDict(curve) 
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =4)
        FPD.add_yPointConstraint( kind = 'equality', location=.15, value=10.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        FPD.add_E1(kind='LS', weight = .1)
        FPD.add_E2(kind='LS', weight = .1)
        FPD.add_E3(kind='LS', weight = .1)
        FPD.add_ArcLengthApprox(kind='LS', weight = .1)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'good_run_1')
        cl.append(copy.deepcopy(Lspline.curve))
        #
        #----------------------------------------------------
        #
        FPD = FormParameterDict(Lspline.curve) 
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        #FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =1)
        #FPD.add_yVertexConstraint(kind = 'equality', value=10., index =1)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        FPD.add_yPointConstraint( kind = 'equality', location=.15, value=10.)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =4)
        FPD.add_AngleConstraint(kind='equality', location = .150001, value =  0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        FPD.add_E1(kind='LS', weight = .1)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .1)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'good_run_2')
        cl.append(copy.deepcopy(Lspline.curve))
        #
        #----------------------------------------------------
        #
        FPD = FormParameterDict(Lspline.curve) 
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        #FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =1)
        #FPD.add_yVertexConstraint(kind = 'equality', value=10., index =1)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        FPD.add_yPointConstraint( kind = 'equality', location=.15, value=10.)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =4)
        FPD.add_AngleConstraint(kind='equality', location = .150001, value =  0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
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
                                           
        FPD.add_E1(kind='LS', weight = .1)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .1)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'good_run_3')
        cl.append(copy.deepcopy(Lspline.curve))
        #
        #----------------------------------------------------
        #
        FPD = FormParameterDict(Lspline.curve)
        FPD.add_AreaConstraint(kind='equality', value = 60.)
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        #FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =1)
        #FPD.add_yVertexConstraint(kind = 'equality', value=10., index =1)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        FPD.add_yPointConstraint( kind = 'equality', location=.15, value=10.)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =4)
        FPD.add_AngleConstraint(kind='equality', location = .150001, value =  0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        #
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 7, index2 = 5, 
                                           seperation = 0. )
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
                                           
        FPD.add_E1(kind='LS', weight = .1)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .1)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'good_run_4')
        cl.append(copy.deepcopy(Lspline.curve))
        #
        #----------------------------------------------------
        #
        FPD = FormParameterDict(Lspline.curve)
        FPD.add_AreaToAnyXAxisConstraint(kind='equality',
                                         value = 120.,
                                         x_axis_loc = 12.)
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        FPD.add_yPointConstraint( kind = 'equality', location=.15, value=10.)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =4)
        FPD.add_AngleConstraint(kind='equality', location = .150001, value =  0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = .95, value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = .955, value = 0.)
        #
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 7, index2 = 5, 
                                           seperation = -4. )
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 7, index2 = 6, 
                                           seperation = -4. )
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
                                           
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'good_run_6')
        cl.append(copy.deepcopy(Lspline.curve))
        #
        #----------------------------------------------------
        #
        Lspline.curve.knot_insertion(.95)
        interval_data, small = interval_bounds(Lspline.curve)
        FPD = FormParameterDict(Lspline.curve)
        FPD.add_AreaToAnyXAxisConstraint(kind='equality',
                                         value = 120.,
                                         x_axis_loc = 12.)
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        FPD.add_yPointConstraint( kind = 'equality', location=.15, value=10.)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =4)
        FPD.add_AngleConstraint(kind='equality', location = .150001, value =  0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = .35, value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = .95, value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = .955, value = 0.)
        #
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 7, index2 = 4, 
                                           seperation = -4. )
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 7, index2 = 5, 
                                           seperation = -4. )
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 7, index2 = 6, 
                                           seperation = -4. )
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
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'good_run_7')
        cl.append(copy.deepcopy(Lspline.curve))
        #
        #----------------------------------------------------
        #
        interval_data, small = interval_bounds(Lspline.curve)
        FPD = FormParameterDict(Lspline.curve)
        FPD.add_AreaToAnyXAxisConstraint(kind='equality',
                                         value = 120.,
                                         x_axis_loc = 12.)
        FPD.add_verticalAngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =2)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =2)
        FPD.add_xVertexConstraint(kind = 'equality', value=0.0, index =3)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =3)
        FPD.add_yPointConstraint( kind = 'equality', location=.15, value=10.)
        FPD.add_yVertexConstraint(kind = 'equality', value=10., index =4)
        FPD.add_AngleConstraint(kind='equality', location = .150001, value =  0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value =  -5.)
        #FPD.add_CurvatureConstraint(kind='equality', location = .35, value = 0.)
        
        ## LB says take out curvature (curved is more realistic):
        #FPD.add_CurvatureConstraint(kind='equality', location = .95, value = 0.)

        #FPD.add_CurvatureConstraint(kind='equality', location = .955, value = 0.)
        #
        #FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
        #                                   value = 0., weight = 1.0, 
        #                                   index = 7, index2 = 4, 
        #                                   seperation = -4. )
        
        #"""
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 7, index2 = 5, 
                                           seperation = -4. )
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = 0., weight = 1.0, 
                                           index = 7, index2 = 6, 
                                           seperation = -4. )
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
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        Lspline.display(mytitle = 'good_run_9')
        cl.append(copy.deepcopy(Lspline.curve))
        
        shift_axis = 12.
        color = 'green'
        Lspline.curve.plotcurve_detailed()
        #Lspline.curve.plotCurvature()
        #Lspline.curve.plotCurvature_spines()
        #Lspline.curve.plotCurvature_nospines(alpha=0.)
        Lspline.curve.plotCurvature_expensive(scale = 1., 
                                              nump=2,
                                              flip=True)
        
        plt.fill_between(Lspline.curve.r[:,0],shift_axis,Lspline.curve.r[:,1],
                            facecolor = color, alpha=.1,
                            label='area')
        plt.annotate('$k_0$', xy=(0.,10.), xycoords='data',
                xytext=(20, 30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate('$k_1$', xy=(Lspline.curve.vertices[7]), 
                     xycoords='data',
                     xytext=(20, 0), textcoords='offset points',
                     arrowprops=dict(arrowstyle="->")
                     )
        plt.annotate('$k_2$', xy=(Lspline.curve.vertices[11]), 
                     xycoords='data',
                     xytext=(20, 0), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->")
                    )
                    
        plt.annotate('', 
                     xy=Lspline.curve.vertices[5], 
                     xycoords='data',
                     xytext=Lspline.curve.vertices[6], 
                             textcoords='data',
                     arrowprops=dict(arrowstyle="<->",
                                      connectionstyle="bar",
                                      ec="k",
                                      shrinkA=5, shrinkB=5,
                                      ))
                    
        plt.annotate('', 
                     xy=Lspline.curve.vertices[7], 
                     xycoords='data',
                     xytext=Lspline.curve.vertices[6], 
                             textcoords='data',
                     arrowprops=dict(arrowstyle="<->",
                                      connectionstyle="bar",
                                      ec="k",
                                      shrinkA=5, shrinkB=5,
                                      ))
        plt.annotate('', 
                     xy=Lspline.curve.vertices[11], 
                     xycoords='data',
                     xytext=Lspline.curve.vertices[7], textcoords='data',
                     arrowprops=dict(arrowstyle="<->",
                                      connectionstyle="bar",
                                      ec="k",
                                      shrinkA=5, shrinkB=5,
                                      ))

        plt.annotate('$r_{1}$', xy=(4.0,6.88), xycoords='data',
                xytext=(-70, 50), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate('$r_{2}$', xy=(1.88,5.0), xycoords='data',
                xytext=(-70, 50), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        plt.annotate('$r_{3}$', xy=(3.5,3.4), xycoords='data',
                xytext=(-70, -50), textcoords='offset points',
                arrowprops=dict(arrowstyle="->")
                )
        
    if do_all:
        test = 'exp7' #real run, not free
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
        test = 'exp8' # not a real run, but is "free"
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