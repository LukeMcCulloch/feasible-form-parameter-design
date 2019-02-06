#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 15:45:57 2018

@author: luke
"""


#
import numpy as np
import copy
import matplotlib.pyplot as plt
import cPickle as pickle
#
from extended_interval_arithmetic import ia
#import sqKanren as lp
#
import curve             as     spline
from   ADILS             import IntervalLagrangeSpline, \
                                Lagrangian, \
                                LazyFPDConstraintInterface
#
from   FormParameter     import FormParameterDict
from   initialValues     import interval_bounds, lagrangian_bounds
#
import thbsurface as thbspline
from myspline import rbspline, BoundsList
import wavelet_matrices as wm
#import THButilities as THBU #only needed by thb curve and surface themselves
import THBderivedUtilities as THBDU #helper to use THB curve and surface
#
import utility_optimization as uopt
#




class MultiSurf(object):
    def __init__(self, **surfaces):
        for key, value in surfaces.iteritems():
            setattr(self, key, value)
        self.hull = self.surfaces[1]
        self.bulb = self.surfaces[2]
            
    
    def match_edge_curve(self, curve_to_change, curve_to_match,dof=1):
        """
        curve_to_change = s1.tcurvelist[0]
        curve_to_match = bulb.tcurvelist[-1]
        
        consider winslow/elliptic/safe laplacian/smoothing
        """
        n = curve_to_match.n
        curve_to_change.vertices[:n] = curve_to_match.vertices
        
        top = curve_to_change.vertices[-1]
        stop = curve_to_change.vertices[n]
        diff = top-stop
        step = diff/(curve_to_change.n-n-1)
        
        start = curve_to_match.vertices[-1]
        total_dist = np.linalg.norm(top - start)
        #print 'total_dist = ',total_dist
        for i in range(n,curve_to_change.n-1):
            local_dis = np.linalg.norm(curve_to_change.vertices[i] - start)
            #print 'local_dis = ',local_dis
            weight = 1.-(local_dis/total_dist)**.5
            #assert(weight>0.),'wtf i = {}'.format(i)
            curve_to_change.vertices[i] = curve_to_change.vertices[i]-weight*step
        return
    
    
    def match_surface_edges(self, s1,s2):
        hull = s1#self.surfaces[s1]
        bulb = s2#self.surfaces[s2]
        #curve = hull.surface.child.tcurvelist[0]
        #curve.vertices[:7] = bulb.surface.tcurvelist[-1].vertices
        self.match_edge_curve(curve_to_change = hull.tcurvelist[0],
                              curve_to_match = bulb.tcurvelist[0])
        hull.set_vertices_from_tcurves(override_v=True)
        hull.compute_surface()
        return
    
    
    
    def return_outline(self, surf):
        """bespoke programming at its finest.
        --that is not a compliment.
        --this is an utter annoyance.
        #
        #
        #
        surf = bulb
        #
        #
        #
        #
        #
        """
        #
        #
        outline= list(surf.lcurvelist[0].vertices)
        for vert in surf.tcurvelist[-1].vertices[1:-1]:
            outline.append(vert)
        top = list(surf.lcurvelist[-1].vertices)
        top.reverse()
        outline+= top
        #
        #
        return outline
    
            
    def plot(self, 
             fancy=False,
             start_hullcurve=0,
             canvas = None):
        keys = self.surfaces.keys()
        
        
        surf = self.surfaces[keys[0]]
        
        if canvas is None:
            ax = surf.surface.plotSurface(fancy=fancy,
                                          show_THBvertices=True,
                                          start_hullcurve=start_hullcurve)
#            ax = surf.surface.plotSurface(fancy=False,
#                                            start_hullcurve=0,
#                                            show_hullcurves=False,
#                                            show_THBcurves=True,
#                                            show_THBvertices=True)
            itthis = keys[1:]
        else:
            ax = canvas
            itthis = keys
            
        for key in itthis:
            surf = self.surfaces[key]
            ax = surf.surface.plotSurface(fancy=fancy,
                                          start_hullcurve=start_hullcurve,
                                          show_THBvertices=True,
                                          canvas = ax)
#            ax = surf.surface.plotSurface(fancy=False,
#                                            start_hullcurve=0,
#                                            show_hullcurves=False,
#                                            show_THBcurves=True,
#                                            show_THBvertices=True,
#                                            canvas = ax)
        return ax
    
    def intersect_top(self, box1, box2):
        """
        Non-general
        intersect the top right pt of  b1 
        with the far left curve of b2
        
        box1 = bulb
        box2 = hull
        
        
        dev
        ---------
        
        b1 = box2.surface
        b2 = box1.surface
        """
        
        pt = np.asarray([box1.xb,box1.ye,0.])
        
        b2 = box2.surface
        b1 = box1.surface
        #pt = b1.tcurvelist[-1](1.)[0]
        
        result = b2.tcurvelist[0].FindPoint(pt[1],1)
        uloc = result[0]
        return uloc
            

class box_surf(object):
    def __init__(self, ku=4, kv=4, nu=4,nv=4, 
                 xb=0.,yb=0.,xe=1.,ye=1., dim=2):
        self.ku = ku
        self.kv = kv
        self.nump = 30
        self.dim = dim
        self.nu = nu # = tcurve.n
        self.nv = nv # = lcurve.n
        self.xb = xb
        self.yb = yb
        self.xe = xe
        self.ye = ye
        self.z = 0.
        
    def make_box(self):
        tverts  = []
        tcurves = []
        lverts  = []
        lcurves = []
        tverts.append( uopt.linear_vertices([self.xb,self.yb,self.z],
                                           [self.xb,self.ye,self.z],
                                           self.nu))
        tcurves.append(
                    rbspline(tverts[0], k=self.ku, nump=self.nump)
                    )
        lverts.append( uopt.linear_vertices([self.xb,self.yb,self.z],
                                           [self.xe,self.yb,self.z],
                                           self.nv))
        lcurves.append(
                    rbspline(lverts[0], k=self.kv, nump=self.nump)
                    )
        for i in range(1,self.nv):
            x = lverts[0][i,0]
            tverts.append(uopt.linear_vertices([x,self.yb,self.z],
                                               [x,self.ye,self.z],
                                               self.nv))
            tcurves.append(
                    rbspline(tverts[-1], k=self.ku, nump=self.nump)
                        )
            
        for j in range(1,self.nu):
            y = tverts[0][j,1]
            lverts.append(uopt.linear_vertices([self.xb,y,self.z],
                                               [self.xe,y,self.z],
                                               self.nv))
            lcurves.append(
                    rbspline(lverts[-1], k=self.ku, nump=self.nump)
                    )
            
        self.tverts = tverts
        self.tcurves = tcurves
        self.lverts = lverts
        self.lcurves = lcurves
        
        self.surface = THBDU.make_surface(self.lcurves,
                                          refine=False)
        
        return
    


if __name__ == """__main__""":
    box1 = box_surf(nu=7,nv=7)
    box1.make_box()
    s0 = box1.surface
    
    ax = box1.surface.plotSurface(fancy=False,
                            start_hullcurve=0)
    
    box2 = box_surf(nu=7,nv=5,
                    xb=.0,xe=-.5,
                    yb=0.,ye=.55)
    box2.make_box()
    bulb = box2.surface
    surfaces = MultiSurf( surfaces={1:box1,2:box2})
#    
    ax = box2.surface.plotSurface(fancy=False,
                            start_hullcurve=0,
                            canvas = ax)
    
    
    box1.surface.dyadic_loft_refinement()
    s1 = box1.surface.child
    
    surfaces.plot()
    
    ##
    ##*******************************8
    ## practice:  get bounds to cover the fwd edge
    #
    # s0.mastercurve.get_bounds_given_index
    #
    
    loc = surfaces.intersect_top(box2,box1) #0.5
    basis_enclosed = box1.surface.child.get_basis_which_are_enclosed(
            thbspline.Box(ia(0.,loc),ia(0.,.125)) )
    
    #basis_enclosed    box1.surface.child.get_basis_which_are_active(
    #            thbspline.Box(ia(0.,loc),ia(0.,.125)) )
    
    box1_indices_of_interest = [el[0] for el in basis_enclosed]
    #basis_for_ji = [b0.mastercurve_u.basis()]
    #bdji = thbspline.BoundsList(box1_indices_of_interest)
    
    #box1_indices_of_interest = [0,1,2,3,4,5]
    #box1_indices_of_interest = [0,1]
    ranges = []
    boxes = []
    for index in box1_indices_of_interest:
        thisb = box1.surface.child.get_bounds_given_indices(index,0)
        ranges.append(thisb )
        boxes.append(thbspline.Box(*thisb))
    
    bp = thbspline.Box( ia(.0,1.),    
                        ia(.0,.25) ) 
    b0 = thbspline.Box( ia(.0,.25),    
                        ia(.0,.25) ) 
    b1 = thbspline.Box( ia(.75,1.),    
                        ia(.75,1.) ) 
    
    
    
    BL = thbspline.BoxList(*boxes )
    #BL = thbspline.BoxList(bp)
    #BL = thbspline.BoxList(b0,b1)
    box1.surface.child.bounds = BL
    surfaces.plot()
    
    
    self = surfaces
    """
    before:
        
        s1.tcurvelist[0].vertices
        Out[437]: 
        array([[0.        , 0.        , 0.        ],
               [0.        , 0.08333333, 0.        ],
               [0.        , 0.20833333, 0.        ],
               [0.        , 0.32291667, 0.        ],
               [0.        , 0.41666667, 0.        ],
               [0.        , 0.5       , 0.        ],
               [0.        , 0.58333333, 0.        ],
               [0.        , 0.67708333, 0.        ],
               [0.        , 0.79166667, 0.        ],
               [0.        , 0.91666667, 0.        ],
               [0.        , 1.        , 0.        ]])
        
        self.bulb.surface.tcurvelist[-1].vertices
        Out[438]: 
        array([[0.        , 0.        , 0.        ],
               [0.        , 0.09166667, 0.        ],
               [0.        , 0.18333333, 0.        ],
               [0.        , 0.275     , 0.        ],
               [0.        , 0.36666667, 0.        ],
               [0.        , 0.45833333, 0.        ],
               [0.        , 0.55      , 0.        ]])
    
        s1.tcurvelist[0].vertices[:self.bulb.surface.nu] - \
             self.bulb.surface.tcurvelist[-1].vertices
             
            
    
    #"""
    
    #"""
    #self.match_surface_edges(s1=1,s2=2)
    self.match_surface_edges(surfaces.surfaces[1].surface.child,
                             surfaces.surfaces[2].surface)
    s1.bounds = thbspline.BoxList(thbspline.Box(ia(0.,1.),ia(0.,1.)))
    s0.compute_surface()
    self.plot()
    #"""
    
    
    
    # NEXT LEVEL:
    #s2 = box1.surface.child.dyadic_loft_refinement()
    #s2.bounds = thbspline.BoxList(thbspline.Box(ia(0.,0.),ia(0.,0.)))
    
    
    
    
    """
    #self.match_surface_edges(s1=1,s2=2)
    self.match_surface_edges(s2,
                             bulb)
    s2.bounds = thbspline.BoxList(thbspline.Box(ia(0.,1.),ia(0.,1.)))
    s0.compute_surface()
    self.plot()
    #"""
    
    oo = self.return_outline(bulb)
    print len(oo)