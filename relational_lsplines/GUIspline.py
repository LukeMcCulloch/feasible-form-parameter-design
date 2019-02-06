# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 19:56:19 2015

@author: lukemcculloch
"""

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.artist import Artist
from matplotlib.mlab import dist_point_to_segment

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.widgets import Button

from menu import Menu, MenuItem, ItemProperties
#from quat1 import CubeAxesInteractive

import curve as spline
from initialValues     import InitializeControlPoints, InitializeControlVertices
#from FormParameter     import FormParameterDict, generalized_aattractor
from   ADILS             import IntervalLagrangeSpline, Lagrangian
from   initialValues     import interval_bounds, lagrangian_bounds
from   FormParameter     import FormParameterDict
from plots import Plotter

from HierarchicalBspline import BsplineHierarchy, LevelSpline

class Index:
    ind = 0
    def next(self, event):
        self.ind += 1
        i = self.ind % len(freqs)
        ydata = np.sin(2*np.pi*freqs[i]*t)
        l.set_ydata(ydata)
        plt.draw()

    def prev(self, event):
        self.ind -= 1
        i = self.ind % len(freqs)
        ydata = np.sin(2*np.pi*freqs[i]*t)
        l.set_ydata(ydata)
        plt.draw()
        
class use_plot(object):
    def __init__(self):
        return
    def add_contraints(interactor):
        interactor.Lspline.display(mytitle = 'Test 1')
        return

class FormParameterInteractor:
    def __init__(self, Lspline):
        self.Lspline = Lspline
        self.Lspline.Lagrangian = None
        self.FPD = FormParameterDict(self.Lspline.curve)
        self.verbose = True
        #self.fig = figure()
        self.fig, self.ax = plt.subplots()
        self.fig.subplots_adjust(left=0.3)
        prop = ItemProperties(labelcolor='black', bgcolor='white',
                           fontsize=15, alpha=0.2)
        hoverprop = ItemProperties(labelcolor='white', bgcolor='blue',
                                fontsize=15, alpha=0.2)
        menuitems = []
        for label in ('tangent', 
                      'curvature', 
                      'area', 
                      'Xc', 
                      'Yc',
                      'clear',
                      'done'):
            def onselect(item):
                print('you selected %s' % item.labelstr)
                getattr(self, item.labelstr)()
                return
            item = MenuItem(self.fig, label, props=prop, 
                            hoverprops=hoverprop,
                            on_select=onselect)
            menuitems.append(item)
        menu = Menu(self.fig, menuitems)
        plt.show(block = False)
        return
        
    def __call__(self):
        return self.Lspline
        
    def tangent(self):
        if self.verbose: print 'tangent FP'
        loc = -1.
        while not (0.<=loc<=1.):
            loc = float(raw_input('enter location for form parameter [0,1] \n'))
        val = float(raw_input('enter desired tangent angle in degrees \n'))        
        self.FPD.add_AngleConstraint(kind='equality', 
                                location = loc, 
                                value = val)
        return
    def curvature(self):
        if self.verbose: print 'curvature FP'
        loc = -1.
        while not (0.<=loc<=1.):
            loc = float(raw_input('enter location for form parameter [0,1] \n'))
        val = float(raw_input('enter desired curvature \n(units are inverse length) \n')) 
        self.FPD.add_CurvatureConstraint(kind='equality', 
                                    location = loc, 
                                    value = val)
        return
    def area(self):
        if self.verbose: print 'area FP'
        val = float(raw_input('enter desired area \n')) 
        self.FPD.add_AreaConstraint(kind='equality', value = val)
        return
    def Xc(self):
        if self.verbose: print 'Xc FP'
        val = float(raw_input('enter desire x centroid \n')) 
        self.FPD.add_XcConstraint(kind='equality', value = val , weight = 2000.)
        return
    def Yc(self):
        if self.verbose: print 'Yc FP'
        val = float(raw_input('enter desire y centroid \n'))
        self.FPD.add_YcConstraint(kind='equality', value = val)
        return
    def clear(self):
        if self.verbose: print 'clear FP'
        self.FPD = FormParameterDict(self.Lspline.curve)
        return
    def done(self):
        if self.verbose: print 'Done'
        interval_data, small = interval_bounds(self.Lspline.curve)
        self.FPD.add_E1(kind='LS', weight = 1.)
        self.FPD.add_E2(kind='LS', weight = .5)
        self.FPD.add_E3(kind='LS', weight = .5)
        #self.FPD.add_ArcLength(kind='LS', weight = 1.)
        self.FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(self.FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        
        self.Lspline = IntervalLagrangeSpline(self.Lspline.curve, L, data = interval_data)
        #self.Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        plt.close(3)
        return
    def set_initial_values(self):
        return



        
        

class BsplineInteractor(object):
    """
    An polygon editor.

    Key-bindings

      't' toggle vertex markers on and off.  When vertex markers are on,
          you can move them, delete them

      'd' delete the vertex under point

      'i' insert a vertex at point.  You must be within epsilon of the
          line connecting two existing vertices
    """
    
    #showverts  = True
    #epsilon = 5 
    def __init__(self, Lspline, 
                 mouse_update_listener=None, 
                 optimization_listener=None):
        self.Lspline    = Lspline
        self.curve      = Lspline.curve
        self.mouse_update_listener = [mouse_update_listener]
        self.optimization_listener = [optimization_listener]
        self.showverts  = True
        self.epsilon    = 5 
        self.verbose    = True
        self.fig        = plt.figure()#figsize=(4, 4)
        #self.fig2       = plt.figure()
        self.ax         = self.fig.add_subplot(111)
        self.fig.add_axes(self.ax)
        self.controlpoly    = Polygon(list(zip( self.curve.vertices[:,0], 
                                                self.curve.vertices[:,1])), 
                                                animated=True, 
                                                fill = False, 
                                                closed = False)
        self.curve_r        = Polygon(list(zip( self.curve.r[:,0], 
                                                self.curve.r[:,1])), 
                                                animated=True, 
                                                fill = False, 
                                                closed = False)
        self.ax.add_patch(self.controlpoly)
        self.ax.add_patch(self.curve_r)
                                                 
        self.line_poly = Line2D(self.curve.vertices[:,0], 
                           self.curve.vertices[:,1], 
                           color='blue',
                           marker='o', 
                           markerfacecolor='r',
                           alpha = .5,
                           animated=True)
                           
        self.line_r = Line2D(self.curve.r[:,0], 
                            self.curve.r[:,1], 
                            color='black',
                            alpha = .75,
                            markerfacecolor='r', 
                            animated=True)
        
        self.ax.add_line(self.line_poly)
        self.ax.add_line(self.line_r)
        
        
        self.cid = self.controlpoly.add_callback(self.poly_changed)
        self.did = self.curve_r.add_callback(self.poly_changed)
        self._ind = None # the active vert

        canvas = self.controlpoly.figure.canvas
        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        
        self.setconstraints = plt.axes([0.50, 0.05, 0.1, 0.075])
        self.sconstraints = Button(self.setconstraints, 'SET')
        self.sconstraints.on_clicked(self.set_contraints)
        
        self.axconstraints = plt.axes([0.61, 0.05, 0.1, 0.075])
        self.bconstraints = Button(self.axconstraints, 'Form Parm')
        self.bconstraints.on_clicked(self.add_constraints)
        
        self.axprev = plt.axes([0.72, 0.05, 0.1, 0.075])
        self.bprev = Button(self.axprev, 'Optimize')
        self.bprev.on_clicked(self.do_optimization)
        
        self.canvas = canvas
        self.ax.set_title('Lspline Plot')
        xmin,xmax,ymin,ymax = self.Lspline.extreme_C0()
        self.ax.set_xlim((xmin-2.,xmax+2.))
        self.ax.set_ylim((ymin-2.,ymax+2.))
        
        return

    def add_constraints(self, event):
        if self.verbose == True: print 'adding constraints'
        self.FPmaker = FormParameterInteractor(self.Lspline)
        return
        
    def set_contraints(self, event):
        if self.verbose == True: print 'setting constraints'
        self.Lspline = self.FPmaker()
        return
        
    def do_optimization(self, event):
        if self.verbose == True: print 'doing optimization'
        self.Lspline.optimize()
        self.compute_curve_plot()
        if self.optimization_listener[0] is not None:
            for func in self.optimization_listener:
                func(self.Lspline.curve.vertices)
        return

    def plot_contraints(self, event):
        if self.verbose == True: print 'plotting constraints'
        self.Lspline.display(mytitle = 'Test 1')
        return    
        
    def compute_curve_plot(self):
        self.controlpoly.xy = self.curve.vertices
        self.line_poly.set_data(zip(*self.curve.vertices)) 
        self.canvas.restore_region(self.background)
        #self.curve.vertices = self.controlpoly.xy
        self.curve.allCurveArray()
        self.curve_r = Polygon( self.curve.r, 
                                animated=True, 
                                fill = False, 
                                closed = False)
        self.curve_r.yx = self.curve.r
        self.line_r.set_data(self.curve.r[:,0], 
                            self.curve.r[:,1])
        self.ax.draw_artist(self.controlpoly)
        self.ax.draw_artist(self.line_poly)
        self.ax.draw_artist(self.curve_r)
        self.ax.draw_artist(self.line_r)
        self.canvas.blit(self.ax.bbox)
        return
    
    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        
        self.ax.draw_artist(self.controlpoly)
        self.ax.draw_artist(self.line_poly)
        
        self.ax.draw_artist(self.curve_r)
        self.ax.draw_artist(self.line_r)
        
        #self.canvas.blit(self.ax.bbox)
        return
        
    def poly_changed(self, controlpoly, curve_r):
        vis     = self.line_poly.get_visible()
        Artist.update_from(self.line_poly, controlpoly)
        self.line_poly.set_visible(vis) 
        vis     = self.line_r.get_visible()
        Artist.update_from(self.line_r, curve_r)
        self.line_r.set_visible(vis) 
        return
    
    def get_ind_under_point(self, event):
        """get the index of the vertex under point if within epsilon tolerance
        """        
        xy = np.asarray(self.controlpoly.xy)
        xyt = self.controlpoly.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        print xt, yt
        print event.x, event.y
        d = np.sqrt((xt-event.x)**2 + (yt-event.y)**2)
        indseq = np.nonzero(np.equal(d, np.amin(d)))[0]
        ind = indseq[0]

        if d[ind]>=self.epsilon:
            ind = None
        return ind
        
    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if self.verbose:
            print 'button_press_callback GUISPLINE'
        print not self.showverts
        print event.inaxes==None
        print event.button != 1
        if not self.showverts: return
        if event.inaxes==None: return
        if event.button != 1: return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        'whenever a mouse button is reileased'
        if not self.showverts: return
        if event.button != 1: return
        self._ind = None
    
    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes: return
        if event.key=='t':
            self.showverts = not self.showverts
            self.line_poly.set_visible(self.showverts)
            if not self.showverts: self._ind = None
        elif event.key=='d':
            ind = self.get_ind_under_point(event)
            if ind is not None:
                self.controlpoly.xy = [tup for i,tup in enumerate(self.controlpoly.xy) if i!=ind]
                self.line_poly.set_data(zip(*self.controlpoly.xy))
        elif event.key=='i':
            xys = self.controlpoly.get_transform().transform(self.controlpoly.xy)
            p = event.x, event.y # display coords
            for i in range(len(xys)):
                s0 = xys[i-1]
                s1 = xys[i]
                d = dist_point_to_segment(p, s0, s1)
                if d<=self.epsilon:
                    self.controlpoly.xy = np.array(
                        list(self.controlpoly.xy[:i]) +
                        [(event.xdata, event.ydata)] +
                        list(self.controlpoly.xy[i:]))
                    self.line_poly.set_data(zip(*self.controlpoly.xy))
                    break
        self.canvas.draw() 
        
    def motion_notify_callback(self, event):
        'on mouse movement'
        if self.verbose:print 'motion_notify_callback'
        if not self.showverts: return
        if self._ind is None: return
        if event.inaxes is None: return
        if event.button != 1: return
        x,y = event.xdata, event.ydata
        self.curve.vertices[self._ind,0] = x
        self.curve.vertices[self._ind,1] = y
        self.compute_curve_plot()
        if self.mouse_update_listener[0] is not None:
            for func in self.mouse_update_listener:
                func(x,y,self._ind)
        return
        
        
class HBsplineInteractor(object):
    def __init__(self,SplineHierarchy,  
                                 mouse_update_listener=None, 
                                 optimization_listener=None):
        #super(BsplineInteractor, self).__init__( 
        #         Lspline, 
        #         mouse_update_listener=None, 
        #         optimization_listener=None)
        self.mouse_update_listener = [mouse_update_listener]
        self.optimization_listener = [optimization_listener]
        self.fig        = plt.figure()
        self.ax         = self.fig.add_subplot(111)
        self.fig.add_axes(self.ax)
        self.SplineHierarchy = SplineHierarchy.SplineHierarchy
        self.nlevels = len(self.SplineHierarchy)
        self.verbose = True
        self.epsilon    = 5
        
        
        # plot the vertices
        # do the update through the offsets
        
        self.controlpolys    = [ Polygon(curve.vertices[:], 
                                    animated=True, 
                                    fill = False, 
                                    closed = False) for curve in self.SplineHierarchy]
        self.curve_rs         = [ Polygon(curve.r[:], 
                                    animated=True, 
                                    fill = False, 
                                    closed = False) for curve in self.SplineHierarchy]
        self.controlpoly = self.controlpolys[0]
        self.line_poly = self.curve_rs[0]
        for patch in self.controlpolys:
            self.ax.add_patch(patch)
        for patch in self.curve_rs:
            self.ax.add_patch(patch)
        
        # plot the vertices
        # do the update through the offsets
        # it will work...!  Say it in German for luck.
                                                 
        self.line_polys =  [ Line2D(curve.vertices[:,0], 
                                   curve.vertices[:,1], 
                                   color='blue',
                                   marker='o', 
                                   markerfacecolor='r',
                                   alpha = .5,
                                   animated=True) for curve in self.SplineHierarchy]
                           
        self.line_rs = [ Line2D(curve.r[:,0], 
                                curve.r[:,1], 
                                color='black',
                                alpha = .75,
                                markerfacecolor='r', 
                                animated=True) for curve in self.SplineHierarchy]
        
        for curve in self.line_polys:
            self.ax.add_line(curve)
        for curve in self.line_rs:
            self.ax.add_line(curve)
            
        #self.cid = self.controlpoly.add_callback(self.poly_changed)
        #self.did = self.curve_r.add_callback(self.poly_changed)
        self._ind = None
        canvas = self.controlpolys[0].figure.canvas
        canvas.mpl_connect('draw_event', 
                           self.draw_callback)
        canvas.mpl_connect('button_press_event', 
                           self.button_press_callback)
        #canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_release_event', 
                           self.button_release_callback)
        canvas.mpl_connect('motion_notify_event', 
                           self.motion_notify_callback)
        self.canvas = canvas
        self.ax.set_title('LHspline Plot')
        self.ax.set_xlim((-2,15))
        self.ax.set_ylim((-2,15))
        
    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        
        for cp, lp in zip(self.controlpolys,self.line_polys):
            self.ax.draw_artist(cp)
            self.ax.draw_artist(lp)
        for  cr, lr in zip(self.curve_rs, self.line_rs):
            self.ax.draw_artist(cr)
            self.ax.draw_artist(lr)
        
        self.canvas.blit(self.ax.bbox)
        return
        
    def compute_curve_plot(self):
        print 'compute curve plot!'
        for i in range(self.nlevels):
            #polygons:
            self.controlpolys[i].set_xy(
                self.SplineHierarchy[i].vertices)            

            self.curve_rs[i].set_xy(
                self.SplineHierarchy[i].r)
            #Line2D
            self.line_polys[i].set_data(
                zip(*self.SplineHierarchy[i].vertices))
            self.line_rs[i].set_data(
                zip(*self.SplineHierarchy[i].r))
            
        #self.canvas.blit(self.ax.bbox)
        self.canvas.draw()
        #self.ax.figure.show()
        #self.draw_callback()
        #self.ax.sh
        return
            
    def get_ind_under_point(self, event):
        """get the index of the vertex under point if within epsilon tolerance
        """    
        ind_list    = []
        dlist       = []
        which_poly  = []
        for i, controlpoly in enumerate (self.controlpolys):
            xy = np.asarray(controlpoly.xy)
            xyt = controlpoly.get_transform().transform(xy)
            xt, yt = xyt[:, 0], xyt[:, 1]
            #print xt, yt
            #print event.x, event.y
            d = np.sqrt((xt-event.x)**2 + (yt-event.y)**2)
            indseq = np.nonzero(np.equal(d, np.amin(d)))[0] #index of min distance.
            #ind = indseq[0]
            ind_list.append(indseq[0])  #best c.p. index for this curve
            dlist.append(d)             #distance for this c.p. index
            which_poly.append(i)        #index of this curve
        
        ind             = None
        the_poly        = None
        self._poly      = None
        self._ind       = None
        #if d[ind]>=self.epsilon:
        #    ind = None
        #return ind 
        for test, d, ply in zip(ind_list,dlist,which_poly):
            if d[test]>=self.epsilon:
                pass
            else:
                ind = test
                the_poly = ply
        
        for i, curve in enumerate(self.SplineHierarchy):
            if ind in curve.offset_reference.keys():
                the_poly = i
                return ind, the_poly 
        
    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if self.verbose:
            print 'button_press_callback GUISPLINE'
        #print not self.showverts
        print event.inaxes==None
        print event.button != 1
        #if not self.showverts: return
        if event.inaxes==None: return
        if event.button != 1: return
        print 'event = ', event
        print self.get_ind_under_point(event)
        self._ind, self._poly = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        'whenever a mouse button is reileased'
        #if not self.showverts: return
        if event.button != 1: return
        self._ind = None
        
    def motion_notify_callback(self, event):
        'on mouse movement'
        if self.verbose:
            print 'motion_notify_callback'
            print event.xdata, event.ydata
            print  '_ind = ', self._ind,'_poly = ', self._poly
        #if not self.showverts: return
        if self._ind is None: return
        if event.inaxes is None: return
        if event.button != 1: return
        x,y = event.xdata, event.ydata
        
        hmap = self.SplineHierarchy[self._poly].offset_reference[self._ind]
        print 'hmap = ',hmap
        self.SplineHierarchy[self._poly].vertices[hmap[0],0] = x
        self.SplineHierarchy[self._poly].vertices[hmap[0],1] = y
        self.SplineHierarchy[self._poly].allCurveArray()
        #self.SplineHierarchy[self._poly].compute_curve()
        #if self._ind in self.SplineHierarchy[self._poly].keys():
        for obsvr in self.SplineHierarchy[self._poly].hierarchical_obsvrs:
                obsvr.propagate_update()
        if self.mouse_update_listener[0] is not None:
            for func in self.mouse_update_listener:
                func(x,y,self._ind)
        self.compute_curve_plot()
        return
        
    
        
    def poly_changed(self, controlpoly, curve_r):
        vis     = self.line_poly.get_visible()
        Artist.update_from(self.line_poly, controlpoly)
        self.line_poly.set_visible(vis) 
        vis     = self.line_r.get_visible()
        Artist.update_from(self.line_r, curve_r)
        self.line_r.set_visible(vis) 
        return    
    
    
        
def setup_optimization(Lspline, 
                       optimization_listener=None,
                       mouse_update_listener=None):
    ILspline = BsplineInteractor(Lspline, 
                                 mouse_update_listener,
                                 optimization_listener)
    
    setconstraints = plt.axes([0.50, 0.05, 0.1, 0.075])
    sconstraints = Button(setconstraints, 'SET')
    sconstraints.on_clicked(ILspline.set_contraints)
    
    axconstraints = plt.axes([0.61, 0.05, 0.1, 0.075])
    bconstraints = Button(axconstraints, 'Form Parm')
    bconstraints.on_clicked(ILspline.add_constraints)
    
    axprev = plt.axes([0.72, 0.05, 0.1, 0.075])
    bprev = Button(axprev, 'Optimize')
    bprev.on_clicked(ILspline.do_optimization)
    
    axnext = plt.axes([0.83, 0.05, 0.1, 0.075])
    bnext = Button(axnext, 'plot')
    bnext.on_clicked(ILspline.plot_contraints)
    #return ILspline

if __name__ == '__main__':
    from   initialValues     import InitializeControlPoints, InitializeControlVertices
    import copy
    
    option  = 'ini_test'
    test    = 'basic'
    interactive = False
    Lspline = None
    interactive = False
    hsp = True
    
    do_widget_example = True#False#
    
    if do_widget_example:
        freqs = np.arange(2, 20, 3)
        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.2)
        t = np.arange(0.0, 1.0, 0.001)
        s = np.sin(2*np.pi*freqs[0]*t)
        #l, = plt.plot(t, s, lw=2)
    
    if option =='ini_test':
        k=4
        nump=30
        xb = 0.
        yb = 12.
        xe = 12.
        ye = 0.
        
        alphae = 0.
        alphab = 0.
        Cab_given = 0.
        Cae_given = 0
        xc = 4.
        yc = 4.
        curve_area = 72.
        slope = 'down'
        
        ini_v = InitializeControlVertices(alphae = 0., alphab = 0.,Cab_given = 0.,Cae_given =0.)
            
        curve = spline.Bspline(ini_v.vertices, k, nump)  
        ae = alphae
        ab = alphab
        
        ini_v = InitializeControlVertices(alphae=ae, alphab=ab,
                                          Cab_given=0.,Cae_given=0.,
                                                              nCV = 7)
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
        Lspline.curve.verbose = True
        self = Lspline
        
        
        
    if interactive == True:
        assert(Lspline !=None)
        #self = setup_optimization(Lspline)
        self = BsplineInteractor(Lspline)
        
        
        
        setconstraints = plt.axes([0.50, 0.05, 0.1, 0.075])
        sconstraints = Button(setconstraints, 'SET')
        sconstraints.on_clicked(self.set_contraints)
        
        axconstraints = plt.axes([0.61, 0.05, 0.1, 0.075])
        bconstraints = Button(axconstraints, 'Form Parm')
        bconstraints.on_clicked(self.add_constraints)
        
        axprev = plt.axes([0.72, 0.05, 0.1, 0.075])
        bprev = Button(axprev, 'Optimize')
        bprev.on_clicked(self.do_optimization)
        
        axnext = plt.axes([0.83, 0.05, 0.1, 0.075])
        bnext = Button(axnext, 'plot')
        bnext.on_clicked(self.plot_contraints)
        
    
    if hsp:
        k=4
        nump=300
        start = [0.,12.]
        end = [12.,0.]
        num = 10
        vv = spline.curvalinear_vertices(start,end,num)
        
        vb = spline.Bspline(vv,k,nump)
        h0spline = LevelSpline(vv,k,nump)
        
        h1spline = h0spline.bounds_refinement(bounds=[.75,.8],
                                              knots=[.77,.775,.78])
                                              
        bsh = BsplineHierarchy([h0spline, h1spline])
        
        
        self = HBsplineInteractor(SplineHierarchy = bsh)