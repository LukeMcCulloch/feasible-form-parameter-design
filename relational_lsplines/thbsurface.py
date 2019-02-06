#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 21:54:18 2017

@author: luke

Notes:  refined surface curves can't know about individual
lower and upper curve vertices

coarse and fine are properties of the surfaces proper.

refined curves are different in number than coarse curves
thanks to bi-directional refinement


Perhaps the errors in evaluation are due to the loft?
Maybe a straight ahead refinement would yield better vertices?


TODO:
    
    USE the MRTHBLspline
    multi-resolution THBspline Lagrangian solver
    to fine tune the individual longitudinals
    of the THB surface.
    
    (thb-refine the hullcurves by themselves 
    if 'in-between' vertices are needed)
    
    Then just update the vertices and bounds
    wherever they have changed 
        - This goes for vertices directly 
            from those solved curves and
            their corresponding knot intervals!
    howto:
    
        -Use active_bounds_given_cv_index
        (say a (v-dir)lcurve 'lives' 
        at (u-dir param val = u)
        
        here is the trick:
            
        lcurve(i).children 
        lives at
        loc = tcruve[any].active_bounds_given_cv_index(i)
        
        a box = BoundsBox([
                lcurve.children.bounds,
                tcurve[any].active_bounds_given_cv_index(i)
        ])
        (or is it just BoundsBox([tcurve,lcurve]) ?)
         
         
         
         stuff you dont need:
             ga = tcurve.greville_abscissa(i)
             stuff used internally:
             tcruve[any].get_range_of_basis(i)
         
         
        active_bounds_given_cv_index
        )

So far we have ignored the effect on the 'in between vertices'
produced during dyadic refinemet...
    
how to make THB surface from refined solved Longitudinal net:
    

    -The 'in between curves' from dyadic refinement in the u dir
        *must reflect the actual position of the refined-original-longitudinals
        *if they are there at all!
    
    -any trick to the results in the u dir
        as a result of refinement in the v dir?
        
    Here is the real trick:
        -Each refined position generated from an old longitudinal
        in fact is a vertex on a new transverse curve at the
        same old level as the transverses before.
        -This is true even when the longitudinals themselves 
        have dyadically more control vertices
        and have been deformed by the solver
        -So, generate new transverse curves at the old level
        but running through all these new longitudinal control points
        at this level...
        -And guess what we do?  REFINE THEM AS NORMAL!
        
    What does this amount too?
        -It amounts to the same lofting as ever, as long 
        as we loft on the deformed positions!
   
    
--The only question is what should the bounds be on the resulting surface?
    -meaning:  how much of the refined parametric domain, away from 
    the deformations, should be swept into the ACTIVE domain
    at that level?
    
    -I thought perhaps the answer is that it has to be all of it, if the span 
    is beteen two curves deformed at 'this' level?
    
    -But actually, 
    SINCE dyadic refinement procedes by:
        1.) dyadically refine the longitudinals
        2.) cast transverses through the extra longitudinal vertices
        3.) dyadically refine the new transverses
    
    THEN:
        1.)  If we merely enforce that all transverses 
            in step 2 pass through the refined, deformed
            longitudinal vertices (well, deformed or not,
            just pass these transverses through those curves!)
                *Then the geometry is 'already there' in the 
                transverse direction at the intermediate level,
                *Just refine what's there in the transverse
                direction to complete the pass up to 
                the dyadically refined surface level!
    
    -The only question now, is what about the ends of the bounds?
    
        *we know that the refined bounds must encompass the entire 
        bounds of activity of any and every refined vertex.  
        (vertex-area of corresponding knot non-zero-ness to the uninitiated)
        *To be conservative, run the 'domain right up to the bounds 
        of the inactive vertices'
        -now, how precisely is that defined, when by definition, 
        they are inactive?
        -This is encompassed by the prior work on 
        'getting the bounds of an index at a coarse level when projected to 
        the next finer level.'
        
        
        You need to double check that function as you have 
        'thought your way into double use of it - 
        -using it in thought experiments where it supplies relevant 
        bounds for a coarse index on the coarse level, 
        a coarse index on the refined level, and 
        possibly even a refined index on the refined level??'
        
        anyway, double check that one!
        
        
        
*******************************************************************************
    ##
    ## Dyadic Loft Refinement Design of THB-spline surfaces
    #       from a set of refined longitudinals THB-spline curves
    #       which may or may not span the refined surface space
    #
    # 1.) make the THB surface as normal, but just the first level
    #
    # 2.) In this already established surface, loop over the longitudinals
    #       solving them/improving them with something like that seen in this file.
    #
    # 3.)  Dyadic lofting
    #
    #       -truly, the only thing that is importatnt is
    #       to match the original longitudinal vertex net.
    #
    #       -which is now multilevel-
    #       
    #       -so during dyadic refinement:
    #       
    #       IF    the longidutinal is just 1 level refined
    #       THEN  The current procedure 
    #           (using the refined longis at one level to make next)
    #           should work
    #           (because the longitudinals span the loft refinement
    #            initial space)
    #
    #       ELIf  longitudinals are refined at multiple levels,
    #               (the spline curves space is way short of the full space)
    #               (and yet it IS enough data to span the first level)
    #               (and thus it IS enough data to prolongate to any level)
    #               (we just want to make that prolonged data come from)
    #               (the extra info on the finer levels of the Lcurves)
    #       THEN   *Build the finest level first:
    #               (prolong all data to the finest curves and reconstruct a surf)
    #               via
    #                   1.)  project Lcurve to finest
    #                   2.) build a vertex net equiv to the desired
    #                       surface net at that level
    #                       (using only the lcurve data)
    #                       
    #                       -build out a surface by sending 
    #                       transverse curves across all vertices
    #                       of this fine level
    #                       
    #                       -refine these transverses to the desired 
    #                       number of times so there is a level
    #                       to match each Lcurve level
    #                       
    #                       -pass new longitudinals through the gaps
    #                       where the refined Tcurves at the finest level
    #                       have made room for new curves
    #                       
    #                       -Finest Net Complete
    #                       
    #              *Restrict the fine net to reconstruct
    #                 each coarser net
    #                       
    #                 -(Or use lower level lcurve details)
    #                  (to incorporate some differences down there??)
    #                       
    #
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from itertools import islice, cycle
import copy
from extended_interval_arithmetic import ia
import extended_interval_arithmetic as iafuncs
from iaBox import BoxList, Box

from myspline import rbspline, BoundsList
import wavelet_matrices as wm
import curve as spline

import  utility_optimization as uopt


from Quaternion import Quaternion
from frames     import Frame,DeepVector


np.set_printoptions(threshold= 200)
#
import cPickle as pickle

#
import THButilities
#
from itertools import product  
#
import os #check if file exists

#def get_3D_instantiation_from_1part_shape
def instantiation3D_from_1part_shape(array,dim=3):
    array_shp = np.shape(array)
    return np.zeros((array_shp[0],array_shp[1],dim))




class level_properties(object):
    def __init__(self, nlevel,nu,nv):
        self.nlevel = nlevel
        self.nu = nu
        self.nv = nv
        
class levels(object):
    def __init__(self, *levels):
        self.levels = []
        for level in levels:
            self.levels.append(level)
            
class level_basis_data(object):
    """Class to hold 
    THBspline method compute_apb(self, u,v)
    data in nice and descriptive form
    
    This is form gathering active basis
    functions at all levels for a given (u,v)
    
    -which is still not yet specific enough
    """
    def __init__(self, surf, data):
        self.surf = surf
        self.level = surf.level
        self.uactive = data[0]
        self.vactive = data[1]
        self.uinactive = data[2]
        self.vinactive = data[3]

    def __repr__(self):
        return 'level{}basis_data()'.format(self.level)

    def __str__(self):
        return 'level{}basis_data()'.format(self.level)
    
    def print_(self):
        print 'uactive:'
        print self.uactive
        print 'vactive:'
        print self.vactive
        print 'uinactive:'
        print self.uinactive
        print 'vinactive:'
        print self.vinactive
        return

class Index(object):
    """
    This object just maps from index (level) to surface at that level
    """
    def __init__(self, surface):
        self.index = {}
        start_surf = surface.get_coarsest_level()
        self.index = self.index_hierarchical_surface(start_surf, 
                                                     level=0,
                                                     index = self.index)
        return
    
    def __call__(self, level):
        return self.index[level]
    
    def __getitem__(self, index):
        return self.index[index]

    def __repr__(self):
        return 'Index({})'.format(self.index)

    def __str__(self):
        return 'Index({})'.format(self.index)
    
    def index_hierarchical_surface(self, surface, level, index):
        if surface.child is None:
            index[level] = surface
            return index
        else:
            index[level] = surface
            return self.index_hierarchical_surface(surface.child, level+1)
    
    def addLevel(self, surface):
        # this was never used:
        #top_surf = surface.get_finest_level()
        #
        #these were 'always' commented out:
        #bot_surf = surface.get_coarsest_level()
        #assert(surface.parent == bot_surf)
        #assert(top_surf == self.index[-1]) #dict keys, not list!!
        next = max(self.index.keys())+1
        """
        assert(top_surf == self.index[next-1])
        #"""
        self.index[next] = surface
        return next
    

def print_(mssg):
    print mssg

def noprint(mssg):
    return


class THBsurface(spline.BsplineSurface):
    """this class builds it's own hierarchy
    where each class instance is a THBsurface 
    that lives within the hierarchy.
    
    The hierarchy is effectively a linked list of surfaces
    
    DONE: add facility for gathering the minimal
    active set of basis from each level, given some (u,v)
    -right now this just happens implicitly 
    during surface evaluation
    
    
    
    surface vertex net (self.vertices):
        
        |----------------------------|
        | <---longitudinal rows----> |
        |           and              |
        | ^   longtitudinal curves   |
        | |                          |
        | |                          |
        | |                          |
        | |                          |
        | transverse     transverse  |
        | | columns  and   curves    |
        | |                          |
        | |                          |
        | |                          |
        | v                          |
        |----------------------------|
        
        THB self hierarchy:
            each "level" of the THBspline surface hierarchy
            is both a surface in the spline sense
            and a node in the 'subdivision tower'
            such that a surface is 'self-organizing'
            i.e. - there is no need for a hierarchy tree
            
    parameters
    ----------
    ulevel = the dyadic 'level' 
            corresponding to nu, 
            the number of vertices in the u direction
    vlevel = the dyadic 'level'
            corresponding to nv,
            the number of vertices in the v direction
    """
    """Nov 2017 desired change (DONE):
        have ulevel and vlevel
        separately (DONE)
        (i.e. u might be n=7, v might be n=11)
        
        this was 'easy enough'
        the tough part would be 
        -intermediate numbers of control points-
        -that is, curves with curve.n
        somewhere in between two dyadic levels
        -also, if we decide to insert vertices to 
        match the next highest level,
        this might as well be done at initial Lspline time
        since otherwise it would require an optimization
        to re-parameterize the curve to the correct knot vector
        (or to dig up some algo to do it another way...)
        
        **The real solution is to change to 
            second generation B-spline 
            wavelet projection operators**
        
    """
    """
    immediate efficiency speedups at creation time:
        -get the projection matrices for all used levels from pickle
        -every time you create a curve or surface pass it's projection
        matrix - the level is always known!
        TODO-currently this has not been done 
        for the THBspline curves only.
        
    things to cleanup:
        -self.level was originally going to be the index into a dict of surface levels
        it makes more sense to return the tuple
        of the 
        (dyadic level in the u directio,
        dyadic level in the v direction)
        -because surface relationships between levels
        are handled more naturally by the tree (linked list)
    """
    def __init__(self, 
                 hullcurvelist = None, #hull transverse curves, Lagrangian designed
                 tcurvelist=None, #surface definition transverses
                 lcurvelist=None, #surface definition longitudinals
                 vertices = None,
                 ku = 4,
                 kv = 4,
                 uknots = None,
                 vknots = None,
                 #numpu = 120,
                 #numpv = 120,
                 numpu = 30,
                 numpv = 30,
                 level = 0,
                 ulevel = 0,
                 vlevel = 0,
                 #
                 index=None,
                 indexu = None,
                 indexv = None,
                 #
                 bounds = BoxList(  Box( ia(0.,1.),ia(0.,1.) )  ),
                 #
                 parent=None,
                 child=None,
                 #
                 trm=None,
                 lrm=None,
                 verbose=False,
                 current_rot = None,
                 tol = 1.e-8,
                 xax = (1.,0.,0.), 
                 yax = (0.,1.,0.), 
                 zax = (0.,0.,1.)  ):
        
        self._print = {0:noprint,1:print_}
        if verbose:
            self.print_switch = 1
        else:
            self.print_switch = 0
        self._verbose = verbose
        self.verbose = verbose
        if self._verbose: print 'TODO add setter to verbose, right now, if you update verbose it will not take'
            
        self.wgen = wm.FindP
        self.trm = None
        self.lrm = None
        self.Pget = THButilities.ProjectionGetter()
        self.ku = ku
        self.kv = kv
        self.uknots = uknots
        self.vknots = vknots
        self.numpu = numpu
        self.numpv = numpv
        self.u = np.linspace(0.,1.,self.numpu)
        self.v = np.linspace(0.,1.,self.numpv)
        self.bounds = bounds
        self.tol = tol
        
        #
        #******************************** curves and vertices
        #
        self.tcurvelist = tcurvelist
        self.lcurvelist = lcurvelist 
        self.vertices_copy = vertices
        self.vertices = vertices
        if vertices is not None:
            self.n       = len(vertices)
        #
        self.surface     = np.zeros((self.numpu,
                                     self.numpv,
                                     3),float) 
        self.flats_outline_exists = False
        #
        #********************************
        #
        ##
        ##**************************************************
        ## Coordinates System
        Origin = np.asarray([0.,0.,0.])
        x, y, z = np.identity(3)
        self.ccs = Frame( np.asarray([x,y,z]) ) #curve coordinate system
        self.origin = DeepVector(v = Origin,
                                B=self.ccs)
        
        ## axes of rotation
        self._ax_yaw        = yax           # self._ax_LR = yax #(0, -1, 0) #Left Right
        self._ax_alt_yaw_   = (0., 0., 1.)     # self._ax_LR_alt = (0, 0, 1)
        self._ax_pitch      = xax           # self._ax_UD = xax #(1,0,0)
        self._ax_roll       = zax
        
        ## initial conditions:
        self.rots = Quaternion.from_v_theta(z, 0.) #null rotation vector
        # define the current rotation axis and angle
        if current_rot is None:
            self._current_rot = Quaternion.from_v_theta((0., 0., 1.), 0.)
        else:
            self._current_rot = current_rot
        self._current_trans_x = 0.
        self._current_trans_y = 0.
        self._current_trans_z = 0.
        ##
        ##**************************************************
        ##
        
        #self.tlevels = [] # real transverse curves
        #self.tcurvelist = [] # surface forming u curves
        #self.lcurvelist = [] # surface forming v curves
        
        self.level_vertices = {}
        self.real_boxes = {} #
        self.knot_boxes = {}
        
        
        self.parent = parent
        self.child  = child
        
        if index is None:
            self.index = Index(self)
            self.level = level
        else:
            self.index = index
            level = self.index.addLevel(self)
            self.level = level
            
            
        #        if indexu is None:
        #            self.indexu = Index(self)
        #            self.levelu = levelu
        #        else:
        #            self.indexu = indexu
        #            levelu = self.indexu.addLevel(self)
        #            self.levelu = levelu
        #            
        #            
        #        if index is None:
        #            self.index = Index(self)
        #            self.level = level
        #        else:
        #            self.index = index
        #            level = self.index.addLevel(self)
        #            self.level = level
            
        
        if vertices is not None:
            #self.vertices = vertices
            self.level_vertices[level] = vertices
        elif lcurvelist is not None:
            nu = len(lcurvelist)
            nv = lcurvelist[0].n
            self.nu = nu
            self.nv = nv
            vertices = np.zeros((nu,nv,3),float)
            self.vertices = vertices
            for j in range(len(lcurvelist)):
                self.vertices[j] = lcurvelist[j].vertices[:]
            
            
        if tcurvelist is None and \
        lcurvelist is None:
            assert(self.vertices is None),'WARNING: THBsurface initialized with vertices but no curves.'
            self.initialize_grid(nx=self.ku,
                                 ny=self.kv)
        else:
            tcurvelist,lcurvelist = self.check_fix_initial_curves(
                                                            tcurvelist,ulevel,
                                                            lcurvelist,vlevel)
            
            self.tcurvelist = tcurvelist
            self.lcurvelist = lcurvelist 
            
            self.uknots = tcurvelist[0].t
            self.vknots = lcurvelist[0].t
            self.ubasis = tcurvelist[0].basis
            self.vbasis = lcurvelist[0].basis
            self.mastercurve_u = tcurvelist[0]
            self.mastercurve_v = lcurvelist[0]
        
        if hullcurvelist is None:
            self.hullcurvelist = copy.deepcopy(self.tcurvelist)
        else:
            self.hullcurvelist = hullcurvelist
        #""" TLM testing Nov 25 2017
        self.set_level(level, ulevel,vlevel)
        #"""
        #""" TLM testing Nov 25 2017
        self.set_rm(tcurve_rm=trm,
                    lcurve_rm=lrm)
        self.expand_initial_curve_list() #correct base tnet to include
        #all lcurve vertices that went into the self.vertex net!
        
        self.set_parent()
        self.nu = self.tcurvelist[0].n
        self.nv = self.lcurvelist[0].n
        
        
        #"""
        #self.set_child()
    ### End of __init__
    ##*****************************************  
    ##
    

    def compute_surface(self):
        """assumes surfNet is already in place (but shall change?)
        
        TODO: make update on change
        
        """
        self.re_sample(self.numpu,
                       self.numpv)
        
        self.surface    = np.zeros((self.numpu,
                                    self.numpv,
                                    3),float)   
        P = self.project_all_to_finest()
        for w1,u in enumerate(self.u):
            for w2,v in enumerate(self.v):
                self.surface[w1,w2,:] += self.eval_smart(u,v,P)[0,0,:]
        return
    
    
    def __call__(self, u,v):
        """THB evaluation for surfaces
        """
        return self.eval_smart(u,v)
    
    
    def make_lines_curves(self, num_curves = 2, k=4, nump=30):
        """
        -descretize one boundary in its native direction, say u.
        -evaluate the surface at each point on the descretization
        -use the y values from this as the points of constant y
        to search for in v 
        as you progress in the z direction
        
        c-s:
            x: transverse
            y: vertical
            z: longitudinal
        """
        p_vector = np.linspace(0.,1.,num_curves,endpoint=True)
        self.waterlines = []
        for p in p_vector:
            start = self(0.,p)[0][0]
            wlc = self.make_waterlines(yloc=start[1],
                                       nump=nump)
            self.waterlines.append( 
                    spline.Bspline(vertices=wlc, 
                                   k=k, nump=nump) )
            
        self.elevations = []
        for p in p_vector:
            start = self(p,0.)[0][0]
            wlc = self.make_elevations(xloc=start[0],
                                       nump=nump)
            self.elevations.append( 
                    spline.Bspline(vertices=wlc, 
                                   k=k, nump=nump) )
        
        self.sections = []
        p_vector = np.linspace(0.,1.,60,endpoint=True)
        for p in p_vector:
            start = self(0.,p)[0][0]
            wlc = self.make_sections(zloc=start[2],
                                     nump=nump)
            self.sections.append( 
                    spline.Bspline(vertices=wlc, 
                                   k=k, nump=nump) )
        """
        np.asarray(polycurve)
        self.sections.append( 
                    spline.Bspline(vertices=np.asarray(polycurve), 
                                   k=k, nump=nump) )
        """
        return
    
    def plot_lines(self, canvas=None):
        if canvas is None:
            ax = None
        else:
            ax = canvas
        try:
            for curve in self.waterlines:
                ax = curve.plot3D_hull_system(canvas=ax,
                                              color='.4')
            for curve in self.elevations:
                ax = curve.plot3D_hull_system(canvas=ax,
                                              color='.4')
            for curve in self.sections:
                ax = curve.plot3D_hull_system(canvas=ax,
                                              color='red')
            self.canvas = ax
        except:
            self.make_lines_curves()
            self.plot_lines()
        return ax
    
    def make_elevations(self, xloc, nump=30):
        """the trick here is that
        -the issue of x indeterminacy
        -is solved by the issue of y indeterminacy
        
        KEEP X constant
        
        """
        p_vector = np.linspace(0.,1.,nump,endpoint=True)
        polycurve = []
        for p in p_vector:
            ans_v = self.findPtloc(pt=xloc,
                                   s=p,
                                   uv_search_Index=0,
                                   r3_index=0)
            polycurve.append(ans_v[1])
        return np.asarray(polycurve)
    
    
    def make_waterlines(self, yloc, nump=30):
        """the trick here is that
        -the issue of x indeterminacy
        -is solved by the issue of y indeterminacy
        
        KEEP Y constant
        
        """
        p_vector = np.linspace(0.,1.,nump,endpoint=True)
        polycurve = []
        for p in p_vector:
            ans_v = self.findPtloc(pt=yloc,
                                   s=p,
                                   uv_search_Index=1,
                                   r3_index=1)
            polycurve.append(ans_v[1])
        return np.asarray(polycurve)
    
    
    def make_sections(self, zloc, nump=30):
        """the trick here is that
        -the issue of x indeterminacy
        -is solved by the issue of y indeterminacy
        
        KEEP Z constant
        
        """
        p_vector = np.linspace(0.,1.,nump,endpoint=True)
        polycurve = []
        for p in p_vector:
            ans_v = self.findPtloc(pt=zloc,
                                   s=p,
                                   uv_search_Index=1,
                                   r3_index=2)
            polycurve.append(ans_v[1])
        return np.asarray(polycurve)
    


    def findPtloc(self, pt, s, uv_search_Index, r3_index):
        """Function to find the v parameterization on a surface, given:
                -pt location in real space which we want to locate the parameterization of
                -s location in parameter space

            inputs:
                surface
                pt x or y or z location
                u or v location
                uv_search_Index = [0,1]   maps => [u,v]
                    0 to search for u
                    1 to search for v
                r3_index = [0,1,2] maps => [x,y,z]
                    0 to search using x
                    1 to search using y
                    2 to search using z

            outputs:
                mid     = u or v location
                q       = [x,y,z] surf evaluation
                found   = True or False
                


            An example:

                say you want to find the v location
                where:
                x=.45
                u=.4

                then:
                pt=.45
                s=.4
                uv_search_Index=1
                r3_index=0
                
        dev
        ----------
        pt=.45
        s=.4
        uv_search_Index=1
        r3_index=0
        """

        span = [[self.uknots[0],self.uknots[-1]],[self.vknots[0],self.vknots[-1]]]

        # algorithm: bianary search
        # range in either u or v:
        high = span[uv_search_Index][-1]
        low  = span[uv_search_Index][0]
        mid  = (low + high)/2. #the initial guess

        search = [s]
        search.insert(uv_search_Index,mid)
        
        q=self.eval(search[0],search[1])[0,0]

        count=0
        while((abs(q[r3_index]-pt))>.00000001 and count<101):
            if pt<q[r3_index]:
                high = mid
            else:
                low = mid
            mid = (low+high)/2.
            search = [s]
            search.insert(uv_search_Index,mid)
            q=self.eval(search[0],search[1])[0,0]
            count += 1
        if count > 100:
            print 'error, point not found'
            found = False
            #return mid, q, found
        else:
            found = True
            #return mid, q, found    

        #mid = param loc
        #q = point in R2 or R3
        #found = True or False
        # count = number of iterations required
        return mid, q, found, count
    
    
    def set_bounds(self, bounds):
        """
        
        
        DEV
        -------
            BoxList()
            from iaBox import BoxList
        """
        self.bounds = bounds
    
    def set_level(self, level, ulevel=None, vlevel=None):
        """desired change:
            have ulevel and vlevel
            separately
            
            this is 'easy enough'
            the tough part will be 
            -intermediate numbers of control points-
            -that is, curves with curve.n
            somewhere in between two dyadic levels
            -also, if we decide to insert vertices to 
            match the next highest level,
            this might as well be done at initial Lspline time
            since otherwise it would require an optimization
            to re-parameterize the curve to the correct knot vector
            (or to dig up some algo to do it another way...)
        """
        if ulevel is None:
            ulevel = level
        if vlevel is None:
            vlevel = level
        self.level = level
        self.ulevel = ulevel
        self.vlevel = vlevel
        self.mastercurve_u.level = ulevel 
        self.mastercurve_v.level = vlevel
        #self.set_hullcurvelist(ulevel)
        self.set_tcurve_level(ulevel)
        self.set_lcurve_level(vlevel)
        return
    
    def set_tcurve_level(self,level):
        for curve in self.tcurvelist:
            curve.level = level
        return
    def set_lcurve_level(self, level):
        for curve in self.lcurvelist:
            curve.level = level
        return
    def set_hullcurvelist(self, level):
        for curve in self.hullcurvelist:
            curve.level = level
        return
    
    def set_a_list_by_an_attr(self, listattr, fieldattr, val):
        selflist = self.__getattribute__(listattr)
        for el in selflist:
            el.__setattr__(fieldattr, val)
        return
    
    def set_rm(self,
               tcurve_rm=None,
               lcurve_rm=None):
        """Set the basis prjection matrices
        these project the u or v directions up one level
        to a higher (finer) level of refinement.
        """
        self._print[self._verbose](
                'func: set_rm')
        print 'setting rm:'
        print 'self.trm = self.wgen({}-1,{}+1) '.format(self.ku,self.ulevel)
        
        if self.trm is None and tcurve_rm is None:
            print 'self.trm is indeed none, calc based on the above'
            try:
                self.trm = self.Pget(self.ulevel)
                print 'FAST Ulevel = ',self.ulevel
            except:
                self.trm = self.wgen(self.ku-1,self.ulevel+1)
            print 'self.trm done, shape = {}'.format(np.shape(self.trm))
        else:
            self.trm = tcurve_rm
            print 'self.trm already found, shape = {}'.format(np.shape(self.trm))
            
        if self.lrm is None and lcurve_rm is None:
            try:
                self.lrm = self.Pget(self.vlevel)
                print 'FAST Vlevel = ',self.vlevel
            except:
                self.lrm = self.wgen(self.kv-1,self.vlevel+1)
        else:
            self.lrm = lcurve_rm
            
        #self.set_a_list_by_an_attr('hullcurvelist','rm',self.trm)
        if tcurve_rm is not None:
            self.set_a_list_by_an_attr('tcurvelist','rm',tcurve_rm)
            self.mastercurve_u.rm = tcurve_rm
        else:
            self.set_a_list_by_an_attr('tcurvelist','rm',self.trm)
            self.mastercurve_u.rm = self.trm
            
        if lcurve_rm is not None:
            self.set_a_list_by_an_attr('lcurvelist','rm',lcurve_rm)
            self.mastercurve_v.rm = lcurve_rm
        else:
            self.set_a_list_by_an_attr('lcurvelist','rm',self.lrm)
            self.mastercurve_v.rm = self.lrm
        
        return
    
    
    def set_top(self):
        """warning:
            only use this if you know what you are doing!
        """
        self.parent = None
        self.mastercurve_u.parent = None
        self.mastercurve_v.parent = None
        self.set_a_list_by_an_attr('hullcurvelist','parent',None)
        self.set_a_list_by_an_attr('tcurvelist','parent',None)
        self.set_a_list_by_an_attr('lcurvelist','parent',None)
        return
    
    def set_parent(self):
        """What if mastercurve_s need to parent and child
        their way around??
        """
        if self.parent is not None:
            self.mastercurve_u.parent = self.parent.mastercurve_u
            self.mastercurve_v.parent = self.parent.mastercurve_v
            self.set_a_list_by_an_attr('hullcurvelist',
                                       'parent',
                                       self.parent.mastercurve_u)
            self.set_a_list_by_an_attr('tcurvelist',
                                       'parent',
                                       self.parent.mastercurve_u)
            self.set_a_list_by_an_attr('lcurvelist',
                                       'parent',
                                       self.parent.mastercurve_v)
            #self.parent.set_child()
            #self.parent.set_parent()
        return
    
    def set_child(self):
        """TODO:
            consider always having the hullcurvelist
            curves parent and child themselves
            since they stay the same number...
            or do they??
        """
        if self.child is not None:
            self.mastercurve_u.children = self.child.mastercurve_u
            self.mastercurve_v.children = self.child.mastercurve_v
            self.set_a_list_by_an_attr('hullcurvelist',
                                       'children',
                                       self.child.mastercurve_u)
            self.set_a_list_by_an_attr('tcurvelist',
                                       'children',
                                       self.child.mastercurve_u)
            self.set_a_list_by_an_attr('lcurvelist',
                                       'children',
                                       self.child.mastercurve_v)
            #self.child.set_child()
            #self.child.set_parent()
        return
    
    def get_tcurvenet(self):
        return self.tcurvelist
    def get_lcurvenet(self):
        return self.lcurvelist
    
    def check_fix_initial_curves(self, 
                                 tcurvenet,ulevel,
                                 lcurvenet,vlevel):
        assert(ulevel is not None),'ERROR, init ulevel is None'
        assert(vlevel is not None),'ERROR, init vlevel is None'
        
        #convert the t curve net to THB form
        if not isinstance(tcurvenet[0],rbspline):
            nnetT = []
            for curve in tcurvenet:
                vtx=curve.vertices
                k = curve.k
                nump = curve.nump
                knots = curve.t
                nnetT.append(rbspline(vtx,k,nump,
                                      t=knots,level=ulevel) )
            tcurvenet = nnetT
        # now convert the lcurve net to THB form
        if not isinstance(lcurvenet[0],rbspline):
            nnetL = []
            for curve in lcurvenet:
                vtx=curve.vertices
                k = curve.k
                nump = curve.nump
                knots = curve.t
                nnetL.append(rbspline(vtx,k,nump,
                                      t=knots,level=vlevel) )
            lcurvenet = nnetL
        
        
        return tcurvenet, lcurvenet
    
    def re_sample(self, numpu,numpv):
        self.numpu = numpu
        self.numpv = numpv
        self.u = np.linspace(0.,1.,self.numpu)
        self.v = np.linspace(0.,1.,self.numpv)
        return
    
    def expand_initial_curve_list(self, 
                                  override = False,
                                  override_u=False,
                                  override_v=False):
        """
        for some reason, i've setup the initialization
        of the tnet of the base space to match the hullcurves
        and -not- the number of vertices in the lcurve net!
        
        This is not right because in actuality 
        the tcurves are divorced from
        the underlying hull-curves.
        
        In fact always and everywhere, there should be 
        a tcurve for every ith lcurve vertex at a given level
        
        this fixes the oversight for a level
        
        -should only ever need be used on the base space!
        -assuming lcurves are the more dense base
        and that dyadic lofting always uses the lcurves to 
        'start off the process to get the next level'
        """
        if override:
            override_u=True
            override_v=True
            
        nu,nv,dim = np.shape(self.vertices) 
        ku = self.ku
        kv = self.kv
        numpu = self.numpu
        numpv = self.numpv
        ulevel = self.mastercurve_u.level
        vlevel = self.mastercurve_v.level
        
        trm = self.mastercurve_u.rm
        lrm = self.mastercurve_v.rm
        
        
        if len(self.tcurvelist)<self.mastercurve_v.n or override_u:
            tlist = []
            for j in range(nv):
                curve = rbspline(self.vertices[:,j],
                                 ku,
                                 numpu,
                                 level = ulevel,
                                 rm = trm)
                tlist.append(curve)
            self.tcurvelist = tlist
        
        if len(self.lcurvelist)<self.mastercurve_u.n or override_v:
            llist = []
            for i in range(nu):
                curve = rbspline(self.vertices[i,:],
                                 kv,
                                 numpv,
                                 level = vlevel,
                                 rm = lrm)
                llist.append(curve)
            self.lcurvelist = llist
        return
    
    def set_tcurves_from_vertices(self):
        """
        let's say there is a process which moves the vertices of
        the tcurvelist at the highest level
        this then sets that level's vertices accordingly.
        """
        for i in range(self.nv):
            self.tcurvelist[i].vertices[:] = self.vertices[:,i]
        return
    
    def set_vertices_from_tcurves(self, override_v=False):
        """
        let's say there is a process which moves the vertices of
        the tcurvelist at the highest(finest?) (how about any) level
        this then sets that level's vertices accordingly.
        """
        for i in range(self.nv):
            self.vertices[:,i] = self.tcurvelist[i].vertices[:]
        if override_v:
            self.set_THB_curves_from_vertex_net(override_v=True)
        return
    
    def set_vertices_from_lcurves(self, override_u=False):
        """
        let's say there is a process which moves the vertices of
        the tcurvelist at the highest(finest?) (how about any) level
        this then sets that level's vertices accordingly.
        """
        for i in range(self.nu):
            self.vertices[i,:] = self.lcurvelist[i].vertices[:]
        if override_u:
            self.set_THB_curves_from_vertex_net(override_u=True)
        return
    
    
    
    def set_THB_curves_from_vertex_net(self,
                                      override_u=False,
                                      override_v=False):
        """This is specifically meant to set the 
        tcurvelist and lcurvelist curves (nothing to do with hull curves!)
        once the THBsurface vertex net has been modified at _THIS_ level
        """
        if override_v:
            for i in range(self.nu):
                self.lcurvelist[i].vertices = self.vertices[i,:]
                self.lcurvelist[i].compute_curve()
                self.lcurvelist[i].vertices_copy = copy.deepcopy(self.lcurvelist[i].vertices)
        if override_u:
            for j in range(self.nv):
                self.tcurvelist[j].vertices = self.vertices[:,j]
                self.tcurvelist[j].compute_curve()
                self.tcurvelist[j].vertices_copy = copy.deepcopy(self.tcurvelist[j].vertices)
        return
    
    def initialize_grid(self, 
                        nx=4, ny=4,
                        vertices=None,
                        nump=30, box=None):
        """Not really used when
        hull curves are availible
        
        see dyadic_loft
        """
        self._print[self._verbose]('func: initialize_grid')
        ku = self.ku
        kv = self.kv
        if box is None:
            box = Box(ia(0.,12.),
                      ia(0.,12.), 
                      ia(0.,12.))
        self.real_boxes[0] = box
            
        knot_box = Box(ia(0.,1.),
                       ia(0.,1.))
        self.knot_boxes[0] = knot_box
        
        
        
        x = np.linspace(box[0][0],box[0][1],nx,endpoint=True)
        y = np.linspace(box[1][0],box[1][1],ny,endpoint=True)
        #z = np.linspace(box[2][0],box[2][1],nz,endpoint=True)
        
        #mesh = np.outer(x,y)
        
        mx,my = np.meshgrid(x,y)
        mz = np.zeros_like(mx)
        s = np.linspace(0.,1.,nx,endpoint=True)
        sy = np.linspace(0.,1.,nx,endpoint=True)
        mz[:] = np.sin(2.*s*np.pi)
        mz[0] = mz[0]*np.cos(sy[0]*2.*np.pi)
        mz[1] = mz[1]*np.cos(sy[1]*2.*np.pi)
        mz[2] = mz[2]*np.cos(sy[2]*2.*np.pi)
        mz[3] = mz[3]*np.cos(sy[3]*2.*np.pi)
        self.mx = mx
        self.my = my
        self.mz = mz
        
        #
        #******************************** start new loft procedure
        #
        self.hullcurvelist = []
        self.tcurvelist = []
        self.lcurvelist = []
        #store real boat refined transverse curves:
        for i in range(ny):
            pts = np.asarray([mx[i],my[i],mz[i]]).T
            lvl = len(pts)-ku 
            #problem with this: it only works for slen(pts) = 5 
            # when k is 4.
            self.hullcurvelist.append(rbspline(pts, 
                                               ku,
                                               nump, 
                                               level=lvl))
        
        #loft them longitudinally. 
#        for i in range(ny):
#            pts = np.asarray([mx[:,i],my[:,i],mz[:,i]]).T
#            lvl = len(pts)-kv
#            self.lcurvelist.append(
#                    rbspline(pts, 
#                             kv,
#                             nump,
#                             interpolate=True, 
#                             level=lvl))
            
        #TODO: use Lspline and make these nice
        self.lcurvelist = uopt.make_THB_lcurvenet(self.hullcurvelist)
            
        # now set vertices.... for fast surface eval later.
        nv_pts = self.lcurvelist[0].n
        nu_pts = len(self.lcurvelist)
        vertices = np.zeros((nu_pts,nv_pts,3),float) #possibly backwards
        
        for row_i, curve in enumerate(self.lcurvelist):
            vertices[row_i,:] =  curve.vertices
            
        #now finalize the transverse net.
        #
        #this transverse curve set
        # must use the vertices of the new 
        #longitudinal set directly
        l2_ucurves = []
        for j_col in range(nv_pts):
            new_curve = rbspline(vertices[:,j_col],
                                 ku,
                                 nump,
                                 level=0)
            #new_curve.bounds = umc.bounds
            #new_curve.parent = umc
            l2_ucurves.append(new_curve)
        #store real transverses:
        self.tcurvelist = l2_ucurves  
        
        #if vertices is not None:
        self.n       = len(vertices)
        self.vertices_copy = vertices
        self.vertices = vertices
        
        ulvl = 0#len(vertices[:,0])-ku
        self.mastercurve_u = rbspline(vertices[:,0], 
                                      ku, 
                                      nump,
                                      level=ulvl)
        vlvl = 0#len(vertices[0,:])-kv
        self.mastercurve_v = rbspline(vertices[0,:], 
                                      ku, 
                                      nump,
                                      level=vlvl)
        self.uknots = self.mastercurve_u.t
        self.vknots = self.mastercurve_v.t
        self.ubasis = self.mastercurve_u.basis
        self.vbasis = self.mastercurve_v.basis
        assert(ulvl == 0),'ERROR, default 4x4 unet is not at level 0!?'
        assert(vlvl == 0),'ERROR, default 4x4 vnet is not at level 0!?'
        self.set_level(level = 0,
                       ulevel = self.mastercurve_u.level,
                       vlevel = self.mastercurve_v.level)
        return
    
    
    
    def get_coarsest_surface(self):
        return self.get_coarsest_level()
    def get_coarsest_level(self):
        """return the coarsest surface level in the hierarchy
        """
        if self.parent is None:
            return self
        else:
            return self.parent.get_coarsest_level()
        
    
    def get_finest_surface(self):
        return self.get_finest_level()
    def get_finest_level(self):
        """return the finest surface level in the hierarchy
        """
        if self.child is None:
            return self
        else:
            return self.child.get_finest_level()
        
        
    def refine_finest(self,
                    nL_refine_depth=1,
                    ubounds=None,
                    vbounds=None):
        surface = self.get_finest_level()
        surface.refine_grid(nL_refine_depth,
                            ubounds,
                            vbounds)
        return
    
        
    def dyadic_loft_refinement(self,
                          bounds = BoxList(  Box( ia(0.,1.), 
                                                 ia(0.,1.) ) ),
                           trm = None,
                           lrm = None): 
        """Dyadic refinment of the THB surface
        without refining the underlying 
        HULL curve net
        
        trm : opportunity to drop in pre-computed projection
        matrices for the next transverse (u) level
        
        lrm : opportunity to drop in pre-computed projection
        matrices for the next longidudinal (v) level
        
        """
        boundslist = [] #bounds for the new level
        use_basic_bounds = True
        
        self._print[self._verbose]('func: dyadic_refinement')
        level = self.level
        ulevel = self.ulevel
        vlevel = self.vlevel
        numpu = self.numpu
        numpv = self.numpv
        
        hullcurves, ucurves, vcurves = self.get_level_all_curves(level)
        vertices = self.vertices
        
        
        """
        Critically for THB-Lspline adaptive Longtitudinal
        optimization in Form Parameter Design
        Followed by consruction of a THB-spline 
        surface from the resulting muiltilevel
        longitudinals,
        dyadically refine only if the next level
        does not already exist!
        
        
        Box = thbspline.Box
        BoxList = thbspline.BoxList
        #"""
        #loft longitudinals:- create intermediate longitudinals
        self.fix_influence_matrix()
        lcurvenet_ini = []
        for li, curve in enumerate(vcurves):
            if curve.children is  None:
                new_curve = curve.dyadic_refinement()
            else:
                use_basic_bounds = False
                #thisga = self.longitudinal_curve_to_ga[li]
                new_curve = curve.children #not this has been a misnomer for a while.  should be curve.child
                width_t = curve.get_range_of_basis(li)
                width_l = new_curve.bounds
                for bds in width_l:
                    boundslist.append(
                            Box(width_t, bds) 
                            )
            lcurvenet_ini.append(new_curve)
        
        if not use_basic_bounds:
            bounds = BoxList(*boundslist)
            
        #prepare intermediate transverse vertices
        """
        Critically for THB-Lspline adaptive Longtitudinal
        optimization in Form Parameter Design
        Followed by consruction of a THB-spline 
        surface from the resulting muiltilevel
        longitudinals,
        vertices here comes from the intermediate 
        longitudinals.
        """
        itverts = np.zeros((len(lcurvenet_ini),
                            lcurvenet_ini[0].n,3),float)
        for i,cv in enumerate(lcurvenet_ini):
            itverts[i,:] = cv.vertices #uses refined vertices, if any!
            
        #make new refined transverses:
        #transverse_curvenet = [rbspline() for verts in itverts]
        N = lcurvenet_ini[0].n
        k = ucurves[0].k
        nump = ucurves[0].nump
        tcurvenet = []
        for i in range(N):
            vtx = itverts[:,i]
            new_curve = rbspline(vtx, k, nump).dyadic_refinement() #rm = self.trm).dyadic_refinement()
            new_curve.level = ulevel + 1
            tcurvenet.append(new_curve)
        
        """Ah, here is the problem - 
        How do we prepare the final curves
        so that the next level of dyadic refinement
        would be able to incorporate the 
        longitudinal optimization data
        from yet another level in the refinement
        hierarchy?
        -Those curves are.....now split into 2 in some cases
        and kept as they were in others?
        -For this to go, there must be some way to 
        tag a curve as belonging to its parent 
        even after transformation
        equivalent to that it has just gone
        through above.
        """
        #**************************************************
        #prep final longitudianl vertices (control vertex net)
        #
        ilverts = np.zeros((tcurvenet[0].n,
                            len(tcurvenet),3),float)
        for i, cv in enumerate(tcurvenet):
            ilverts[:,i] = cv.vertices
        
        #**************************************************
        # and final longitudinal curves
        #
        N = tcurvenet[0].n
        lcurvenet = []
        for i in range(N):
            vty = ilverts[i,:]
            new_curve = rbspline(vty, k, nump)
            new_curve.level = vlevel + 1
            lcurvenet.append(new_curve)#,rm = self.lrm))
        
        
        self.child = THBsurface(hullcurves,
                                tcurvenet,
                                lcurvenet,
                                ilverts,
                                ku = self.ku,
                                kv = self.kv,
                                level = self.level+1,
                                ulevel = ulevel+1,
                                vlevel = vlevel+1,
                                trm = trm,
                                lrm = lrm,
                                index = self.index,
                                bounds=bounds,
                                parent=self,
                                verbose = self._verbose)
#        
        self.child.mastercurve_u.parent = self.mastercurve_u
        self.mastercurve_u.children = self.child.mastercurve_u
#        
        self.child.mastercurve_v.parent = self.mastercurve_v
        self.mastercurve_v.children = self.child.mastercurve_v
        
        
        
        
#        #""" TLM testing Nov 25 2017
#        self.child.set_rm()
#       
        """
        self.set_parent()
        self.set_child() #try taking this out Nov 27 2017
        self.child.set_parent()
        self.child.set_child()
#        #"""
#        self.child.mastercurve_u.level = ulevel+1
#        self.child.mastercurve_v.level = vlevel+1
#        
#        self.child.assign_bounds_to_master_curves()
        
        return self.child
    
    
    
    def master_curve_refinement(self, 
                               level,
                               ubounds=BoundsList([ia(0.0, 1.0)]),
                                vbounds=BoundsList([ia(0.0, 1.0)]) ):
        
        umaster_curve, vmaster_curve = self.get_level_base_curves(level)
        
        
        
        unew_curve = umaster_curve.dyadic_refinement()
        unew_curve.bounds = ubounds
        unew_curve.parent = umaster_curve
        umaster_curve.children = unew_curve
        #self.nHLevels['u'] +=1
        
        
        vnew_curve = vmaster_curve.dyadic_refinement()
        vnew_curve.bounds = vbounds
        vnew_curve.parent = vmaster_curve
        vmaster_curve.children = vnew_curve
        #self.nHLevels['v'] +=1
        
        self.set_level_master_curves(unew_curve, vnew_curve, 
                                     ubounds,vbounds,
                                     level+1)
        return
    
    
    def loft_through_transverse(self, tcurve_set):#, unew_curve):
        """
        use this  --or something like it
        within dyadic loft refinement to
        loft:
            
            longitudinals 
            through freshly dyadically refined transverse curves

        
        vertices[:,j] = the transvere pts in column j
        
        lvertices[i,:] = the longitudinal pts in row i
        """
        if self._verbose: print 'thbsurface function: loft_through_transverse'
        nu_pts = tcurve_set[0].n
        nv_pts = len(tcurve_set)
        vertices = np.zeros((nu_pts,nv_pts,3),float) 
        
        #vertices[:,j] = the transvere pts in column j:
        loft_set = []
        for col_j, curve in enumerate(tcurve_set):
            vertices[:,col_j] = curve.vertices 

        
        #lvertices[i,:] = the longitudinal pts in row i:
        k = curve.k
        numpv = curve.nump
        for row_i in range(nu_pts):
            loft_set.append(spline.interpolatedBspline(
                                vertices[row_i,:],
                                k,numpv))
            
        self.save_loft_set = loft_set
        self.save_vertices = vertices
        self.save_nu_pts = nu_pts
        nn = []
        for curve,row_i in zip(loft_set, range(nu_pts)):
            ipts = vertices[row_i,:]
            #print 'using verts = ',
            #print vertices
            nn.append(uopt.interpolate_vertices(
                                            initial_vertices = ipts,
                                            helper_curve = curve,
                                            return_Lspline = False,
                                            make_thb=True)
                                                )
        
        loft_set = nn
        return loft_set
    
    
    
    #    def refine_grid_u(self, 
    #                    nL_refine_depth=None, 
    #                    bounds=1):
    #        self.refine_grid('u',nL_refine_depth,bounds)
    #        return
    #    
    #    def refine_grid_v(self, 
    #                    nL_refine_depth=None, 
    #                    bounds=1):
    #        self.refine_grid('v',nL_refine_depth,bounds)
    #        return
    def active(self, span,p):
        return range(span-p,span+1)
    
    def get_surface_by_level(self, level):
        return self.index(level)
    
    
        
    def return_transverse_vertices(self, level, index):
        surface = self.index(level)
        return surface.vertices[index,:]
        #return self.level_vertices[level][index,:]
    
        
    def return_longitudinal_vertices(self, level, index):
        surface = self.index(level)
        return surface.vertices[:,index]
        #return self.level_vertices[level][:,index]
        
    
    def get_level_all_curves(self, level):
        surface = self.index(level)
        return surface.hullcurvelist, surface.tcurvelist, surface.lcurvelist
        #return self.ulevels[level], self.vlevels[level]
    
    def get_level_base_curves(self, level):
        return self.get_level_master_curves(level)
    
    def get_level_master_curves(self, level):
        surface = self.index(level)
        return surface.mastercurve_u, surface.mastercurve_v
        #return self.HL['u'][level], self.HL['v'][level]
    
    def set_level_master_curves(self, 
                                ucurve, 
                                vcurve,
                                ubounds, 
                                vbounds,
                                level):
        ucurve.bounds = ubounds
        vcurve.bounds = vbounds
        #self.HL['u'][level] = ucurve
        #self.HL['v'][level] = vcurve
        self.mastercurve_u = ucurve
        self.mastercurve_v = vcurve
        return
        
    
    def get_level_vertices(self, level):
        """Adding a level of indirection here
        because I am not sure how I want these
        stored and accessed.
        """
        #get the vertices
        vertices = self.level_vertices[level]
        return vertices
    
    def surface_point(self, u,v):
        """
        Done?: find out why this 
        becomes wrong at higher levels of refinement
        
        answer:  it is not wrong
        
        >>> gg.surface_point(.5,.5,0)
        array([[ 6.,  0.]])
        
        >>> gg.HL['v'][0]
        <class curve.Bspline>
        >>> gg.HL['v'][0].eval(.5)
        
        array([[ 0.,  6.]])
        >>> gg.HL['u'][0].eval(.5)
        array([[ 6.,  0.]])
        
        """
        ucurve,vcurve = self.get_level_base_curves(self.level)
        p = ucurve.p
        q = vcurve.p
        uspan = ucurve.FindSpan(u)
        vspan = vcurve.FindSpan(v)
        uBasisFuns = ucurve.TLMBasisFuns(u)
        vBasisFuns = vcurve.TLMBasisFuns(v)
        
        
        P = self.vertices #self.get_level_vertices(level)
        
        uind = uspan-p
        #ushape = np.shape(uBasisFuns)
        #vshape = np.shape(vBasisFuns)
        #Sarray = np.zeros((ushape[0],ushape[1],3))
        S = 0.
        ckS = 0.
        for l in range(q+1):
            temp = 0.
            vind = vspan-q+l
            #Sarray[l] = np.dot( uBasisFuns.T, P)
            for k in range(0,p+1):
                temp += uBasisFuns[uind+k]*P[uind+k][vind]
            S += vBasisFuns[vspan-q+l]*temp
                
        ckS = ckS + np.dot(uBasisFuns.T, np.dot(vBasisFuns.T,P))
        assert(np.linalg.norm(ckS[0][0:] - S) <.0001),'ERROR: ckS != S::  {} ??= {}'.format(ckS[0], S)
        return S
    
    
    #def assign_bounds_to_master_curves(self, u,v,this_surf=None):
    def assign_bounds_to_master_curves(self, this_surf=None):
        if this_surf is None:
            sf = self.get_finest_level()
            this_surf = sf
            self.assign_bounds_to_master_curves(this_surf)
        else:
            #bounds = this_surf.bounds.get_box_containing(u,v)
            ub = BoundsList([el[0] for el in this_surf.bounds])
            vb = BoundsList([el[1] for el in this_surf.bounds])
            this_surf.mastercurve_u.bounds = ub#BoundsList([bounds[0]])
            this_surf.mastercurve_v.bounds = vb#BoundsList([bounds[1]])
            if this_surf.parent is None:
                return
            else:
                this_surf = this_surf.parent
                self.assign_bounds_to_master_curves(this_surf)
        return
                
    
    
    #"""----------------------------------------------------------
    #"""----------------------------------------------------------
    #"""----------------------------------------------------------
    #curve = ucurve
    def recursive_to_bottom(self, 
                            ucurve, 
                            vcurve, 
                            u, v, 
                            ubasis_levels=None, 
                            vbasis_levels=None,
                            uinactive_bl=None,
                            vinactive_bl=None):
        """u=.7
        """
        if ubasis_levels is None:
            ubasis_levels = []
            vbasis_levels = []
            uinactive_bl = []
            vinactive_bl = []
        if ucurve.children is None:
            pkg = self.bottom_level(ucurve, u, vcurve, v)
            ubasis, vbasis, uinactive_Basis, vinactive_Basis = pkg
            
            ubasis_levels.append(ubasis)
            vbasis_levels.append(vbasis)
            
            uinactive_bl.append(uinactive_Basis)
            vinactive_bl.append(vinactive_Basis)
            
            return ubasis_levels, vbasis_levels, uinactive_bl, vinactive_bl, ucurve, vcurve
        else:
            return self.recursive_to_bottom(ucurve.children,
                                            vcurve.children, 
                                            u,v,
                                            ubasis_levels,
                                            vbasis_levels,
                                            uinactive_bl,
                                            vinactive_bl)
    #"""
    
    #"""
    def bottom_level(self, 
                     ucurve, u,
                     vcurve, v):
        #assert(self.children is None)
        uoBasis = ucurve.TLMBasisFuns(u)
        uactive_parents, apb = ucurve.get_one_level(u)
        
        
        voBasis = vcurve.TLMBasisFuns(v)
        vactive_parents, vapb = vcurve.get_one_level(v)
        
        
        if len(uactive_parents)==0:
            vactive_parents=[]
        if len(vactive_parents)==0:
            uactive_parents=[]
#        
        vpset = set(vactive_parents)
        vapbset = set(vapb)
        upset = set(uactive_parents)
        uapbset = set(apb)
        
        if len(vapbset.intersection(vpset))==0:
            uactive_parents = []
            vactive_parents = []
        if len(uapbset.intersection(upset))==0:
            vactive_parents = []
            uactive_parents = []
#        
        uactiveBasis = self.clip_basis_outside_refinement(uoBasis, 
                                                uactive_parents)
        vactiveBasis = self.clip_basis_outside_refinement(voBasis, 
                                                vactive_parents)
#        usave = uactiveBasis.copy()
#        vsave = vactiveBasis.copy()
#        uactiveBasis = uactiveBasis*sum(vsave)
#        vactiveBasis = vactiveBasis*sum(usave)
        
        uinactive_Basis = uoBasis - uactiveBasis
        vinactive_Basis = voBasis - vactiveBasis
        
#        if (uactiveBasis == 0.).all():
#            vactiveBasis[:] = 0.
#            vinactive_Basis = voBasis - vactiveBasis
#            #assert((voBasis == 0.).all())
#            #voBasis = self.clip_basis_outside_refinement(voBasis, 
#            #                                    [])
#        elif (vactiveBasis == 0.).all():
#            uactiveBasis[:] = 0.
#            uinactive_Basis = uoBasis - uactiveBasis
#            #assert((uoBasis == 0.).all())
#            #uoBasis = self.clip_basis_outside_refinement(uoBasis, 
#            #                                    [])
        return uactiveBasis, vactiveBasis, uinactive_Basis, vinactive_Basis
    
    
    def get_one_level(self, u, v, curve, vcurve):
        """
        TODO: uses bounds[0]
        
        """
        active_parent_bounds_index = curve.bounds.contains(u)
        len_apbi = len(active_parent_bounds_index)
        assert(len(active_parent_bounds_index)<2),'error: parent bounds active in two spans at once'
         
        vactive_parent_bounds_index = vcurve.bounds.contains(u)
        vlen_apbi = len(vactive_parent_bounds_index)
        assert(len(vactive_parent_bounds_index)<2),'error: parent bounds active in two spans at once'
        
        
        if len_apbi==1 and vlen_apbi==1:
            parent_bounds = curve.bounds[active_parent_bounds_index[0]]
        else:
            parent_bounds = ia(0.,0.)#curve.bounds[0]#dummy_bounds#
        
        active_parents = curve.enclosed_basis_in_range(
                                        parent_bounds
                                        ) #basis within Omega at this level
        apb = curve.bases(curve.active(u)) #active basis at u at this level
        
        
        if vlen_apbi==1 and len_apbi==1:
            vparent_bounds = vcurve.bounds[vactive_parent_bounds_index[0]]
        else:
            vparent_bounds = ia(0.,0.)#vcurve.bounds[0]#dummy_bounds#
        
        vactive_parents = vcurve.enclosed_basis_in_range(
                                        vparent_bounds
                                        ) #basis within Omega at this level
        vapb = vcurve.bases(vcurve.active(u)) #active basis at u at this level
        
        return active_parents, apb, vactive_parents, vapb
    #"""----------------------------------------------------------
    #"""----------------------------------------------------------
    #"""----------------------------------------------------------
    def get_all_active_indices(self):
        """
        -returns all active indices in u and v directions
            per level. (period, as in not dependent on u,v))
        -returns all inactive indices in the u and v directions
            per level. (period, as in not dependent on u,v)
            
        >> Jn 26 2018 
        Conjecture:
            perhaps you need to set bounds at a higher level
        and then plot the bounded pieces of a level
        as seperate entities??
        """
        #bounds = self.bounds.get_box_containing(u,v)
        tot_uai = []
        tot_vai = []
        tot_uii = []
        tot_vii = []
        for bounds in self.bounds:
            ucurve = self.mastercurve_u
            vcurve = self.mastercurve_v
            ucurve.bounds = BoundsList([bounds[0]])
            vcurve.bounds = BoundsList([bounds[1]])
            
            uai = vai = set() 
            uii = vii = set()
            utotal_indices = range(ucurve.n)
            vtotal_indices = range(vcurve.n)
            
            #u = bounds[0].midpoint()
            #uactive, unon_zero = ucurve.get_one_level(u)
            uactive = ucurve.get_all_one_level()
            uti = set(utotal_indices)
            uai = set(uactive)
            uii = uti - uai
        
            #v = bounds[1].midpoint()
            #vactive, vnon_zero = vcurve.get_one_level(v)
            vactive = vcurve.get_all_one_level()
            vti = set(vtotal_indices)
            vai = set(vactive)
            vii = vti - vai
                
            tot_uai += list(uai)
            tot_vai += list(vai)
            tot_uii += list(uii)
            tot_vii += list(vii)
        return list(set(tot_uai)), list(set(tot_vai)), list(set(tot_uii)), list(set(tot_vii))
        #return list((tot_uai)), list((tot_vai)), list((tot_uii)), list((tot_vii))
    
    
    
    def get_active_inactice_indices(self, u,v):
        """**Important Function
        
        -returns all _active indices in u and v directions
            per level at a given (u,v)
        -returns all _inactive indices in the u and v directions
            per level at a given (u,v)
            
        -different that get_all_active_indeices 
        in that it gets only those active at a particular u,v.
        
        note the use of this code:
            bounds = self.bounds.get_box_containing(u,v)
            which means that this function acts on one
            bounds box.
        """
        bounds = self.bounds.get_box_containing(u,v)
        ucurve = self.mastercurve_u
        vcurve = self.mastercurve_v
        ucurve.bounds = BoundsList([bounds[0]])
        vcurve.bounds = BoundsList([bounds[1]])
        
        uai = vai = set() 
        uii = vii = set()
        utotal_indices = range(ucurve.n)
        vtotal_indices = range(vcurve.n)
        
        if bounds.contains(u,v):
            uactive, unon_zero = ucurve.get_one_level(u)
            uti = set(utotal_indices)
            uai = set(uactive) #this is different to influential_indices function below
            uii = uti - uai
        
        
            vactive, vnon_zero = vcurve.get_one_level(v)
            vti = set(vtotal_indices)
            vai = set(vactive)
            vii = vti - vai
        else:
            uii = set(utotal_indices)
            vii = set(vtotal_indices)
            
        return list(uai), list(vai), list(uii), list(vii)
    
    
    
    
    def get_influential_indices(self, u,v):
        """
        -returns all _possibly_ nonzero indices in u and v directions
            per level at a given (u,v)
            
            (regardless of their actual use)
            
        note the use of this code:
            bounds = self.bounds.get_box_containing(u,v)
            which means that this function acts on one
            bounds box.
        """
        bounds = self.bounds.get_box_containing(u,v)
        ucurve = self.mastercurve_u
        vcurve = self.mastercurve_v
        ucurve.bounds = BoundsList([bounds[0]])
        vcurve.bounds = BoundsList([bounds[1]])
        
        uai = vai = set() 
        uii = vii = set()
        utotal_indices = range(ucurve.n)
        vtotal_indices = range(vcurve.n)
        
        if bounds.contains(u,v):
            uactive, unon_zero = ucurve.get_one_level(u)
            uti = set(utotal_indices)
            uai = set(unon_zero) #this is different to active_inactive function above
            uii = uti - uai
        
        
            vactive, vnon_zero = vcurve.get_one_level(v)
            vti = set(vtotal_indices)
            vai = set(vnon_zero)
            vii = vti - vai
        else:
            uii = set(utotal_indices)
            vii = set(vtotal_indices)
        
        return list(uai), list(vai), list(uii), list(vii)
    
    
    def get_basis_which_are_enclosed(self, box):
        """Irrespective of the actual bounds in place
        """
        urange = self.mastercurve_u.enclosed(box[0])
        vrange = self.mastercurve_v.enclosed(box[1])
        return list(product(urange,vrange))
    
    def get_basis_which_are_active(self, box):
        """Irrespective of the actual bounds in place
        """
        urange = self.mastercurve_v.bases(
                    self.mastercurve_u.active(box[0]) )
        vrange = self.mastercurve_v.bases(
                    self.mastercurve_v.active(box[1]) )
        return list(product(urange,vrange))
    
    
    def get_bounds_given_indices(self, i,j):
        """
        Irrespective of the actual bounds in place
        
        surface counterpart to the following identical thbspline curve
        funtions:
            -get_bounds_given_index
            -get_bounds_given_index
            -active_bounds_given_cv_index
            
            
        see also  
        ----------
        get_box_containing(u,v) : return bounds box given u,v
        
        """
        #urange = self.mastercurve_u.get_range_of_basis(i)
        #vrange = self.mastercurve_v.get_range_of_basis(j)
        
        # use a bit of abstraction (in name only):
        urange = self.mastercurve_u.get_bounds_given_index(i)
        vrange = self.mastercurve_v.get_bounds_given_index(j)
        return urange,vrange
    
    
    
    def compute_apb(self, u,v):
        """
        return a dictionary of level_basis_data objects
        1 for each level
        
        each object contais lists of these properties 
        for a particular level:
            uactive : the u basis functions active and non zero 
            vactive : the v basis function active and non zero
            uinactive : just what you think
            vinactive : just what you think
        (and note that the particular level also happens
        to be the dictionary key)
        """
        dapb = {}
        surf = self.get_finest_level()
        while surf is not None:
            apb = surf.get_influential_indices(u,v)
            dapb[surf.level] = level_basis_data(surf, apb)
            surf = surf.parent
        return dapb
    
    
    
    def compute_active_by_level(self, u,v):
        """
        return a dictionary of level_basis_data objects
        1 for each level
        
        each object contais lists of these properties 
        for a particular level:
            uactive : the u basis functions active (not necessarily non-zero?) 
            vactive : the v basis function active (not necessarily non-zero?)
            uinactive : just what you think
            vinactive : just what you think
        (and note that the particular level also happens
        to be the dictionary key)
        """
        dapb = {}
        surf = self.get_finest_level()
        while surf is not None:
            active = surf.get_active_inactice_indices(u,v)
            dapb[surf.level] = level_basis_data(surf, active)
            surf = surf.parent
        return dapb
    
    
    
    #"""----------------------------------------------------------
    #"""----------------------------------------------------------
    #"""----------------------------------------------------------
    def compute_total_active_by_level(self):
        """
        Uses:
            (1) find the active control vertices for plotting
        
        Function:
            return a dictionary of level_basis_data objects
            1 for each level
            
            each object contais lists of these properties 
            for a particular level:
                uactive : the u basis functions active 
                vactive : the v basis function active
                uactive : the u basis functions inactive 
                vactive : the v basis function inactive
            (and note that the particular level also happens
            to be the dictionary key)
        """
        dapb = {} #Dictionary Active Parent Basis
        # a name which is derived from
        # the way I originally created 
        # truncated basis functions
        surf = self.get_finest_level()
        while surf is not None:
            active = surf.get_all_active_indices()
            dapb[surf.level] = level_basis_data(surf, active)
            surf = surf.parent
        return dapb
    #"""----------------------------------------------------------
    #"""----------------------------------------------------------
    #"""----------------------------------------------------------
    def compute_true_nonzero_basis(self, u,v):
        """
        return a level_basis_object
        categorizing the real basis functions
        on each level
        which are really non zero 
        at this (u,v)
        
        This one collects only the active basis 
        actually used in computation 
        for a given u,v
        
        access: 
            dapb[3].uactive     = [3,4,5]
            dapb[3].vactive     = [3,4,5]
            dapb[3].uinactive   = [0,1,2,6]
            dapb[3].vinactive   = [0,1,2,6]
        
        """
        this = self.get_coarsest_level()
        ufunc = this.mastercurve_u
        vfunc = this.mastercurve_v
        ub = ufunc.rthb_basis(u)
        vb = vfunc.rthb_basis(v)
        #now gather them up:
        dapb = {}
        count = 0
        uai = []
        vai = []
        uii = []
        vii = []
        this = self.get_finest_level()
        while this is not None:
            #uai = [el for el in np.nonzero(ub[count])[0]]
            #vai = [el for el in np.nonzero(vb[count])[0]]
            #uai = [el for el in np.nonzero(ub[count] != 0.)[0]]
            uai = list( np.nonzero(ub[count] != 0.)[0])
            vai = list( np.nonzero(vb[count] != 0.)[0])
            uii = list( np.nonzero(ub[count] == 0.)[0])
            vii = list( np.nonzero(vb[count] == 0.)[0])
            dapb[this.level] = level_basis_data(this, (uai,
                                                        vai,
                                                        uii,
                                                        vii))
            count+=1
            this = this.parent
        return dapb
        
    #"""----------------------------------------------------------
    #"""----------------------------------------------------------
    #"""----------------------------------------------------------
    def get_current_level(self, u,v):
        """
            input   : (u,v) pair denoting a location on
                        the B-spline surface parameterization
            output  : top surface level containing this (u,v) pair
            
        Note:
            see 
            'Adaptive CAD model (re-)construction 
                with THB-splines p.7'
            for use in adaptive CAD fitting
            
        
        Usage:
            approximatly interpolate a set of points.
            compute the L2 error.
            given whatever representation levels, l, 
            that we are using...
            if [(x-s(p))**2]**.5 > threshold, 
            then add the area 
            around p
            (i.e. the interval of the
            interior of the knot span 
            where we have p)
            to level (l+1)
        """
        
        level = self.level
        self_contains_bounds = self.bounds.contains(u,v)
        
        if self_contains_bounds:
            if self.child is not None:
                child_contains_bounds = self.child.bounds.contains(u,v)
                if child_contains_bounds:
                    return self.child.get_current_level(u,v)
                else:
                    return level
            else:
                return level
        else:
            if self.parent is not None:
                return self.parent.get_current_level(u,v)
            else:
                return False
    #"""----------------------------------------------------------
    #"""----------------------------------------------------------
    #"""----------------------------------------------------------
    
    def THBasis(self,u,v,ucurve,vcurve):
        """
        Not used
        """
        dummy_bounds = ia(0.,0.)
        pkg = self.recursive_to_bottom(ucurve,vcurve,u,v)
        ubasis_levels, vbasis_levels, uinactive_basis_levels, vinactive_basis_levels, ulcurve, vlcurve = pkg
        #vbasis_levels, vlcurve = vcurve.recursive_to_bottom(v)
        ulevel_mats = []
        vlevel_mats = []
        #TODO in THBsurface: make ulcurve.parent make sense.
        #note that it is only needed for basis, -whew-
        while ulcurve.parent is not None: #TODO: make this work for surface curves!
        
            ulcurve = ulcurve.parent
            vlcurve = vlcurve.parent
            
            #uiBasis =  ubasis_levels[-1]
            #viBasis =  vbasis_levels[-1]
            #uiBasis = uinactive_basis_levels[-1]
            #viBasis = vinactive_basis_levels[-1]
            
            pkg = self.get_one_level(u,v,
                                     ulcurve, vlcurve) #!uses bounds to get apb

            uactive_parents, uapb, vactive_parents, vapb = pkg
            #uactive_parents, uapb = ulcurve.get_one_level(u) #!uses bounds to get apb
            uoBasis = ulcurve.TLMBasisFuns(u)
            uchild = ulcurve.children
            
            
            #vactive_parents, vapb = vlcurve.get_one_level(v) #!uses bounds to get apb
            voBasis = vlcurve.TLMBasisFuns(v)
            vchild = vlcurve.children
            
            if len(uapb)==0:
                vapb = []
                vactive_parents=[]
            elif len(vapb)==0:
                uapb = []
                uactive_parents=[]
                
            if len(uactive_parents)==0:
                vactive_parents=[]
                vapb = []
                
            if len(vactive_parents)==0:
                uactive_parents=[]
                uapb = []
                
                
            vpset = set(vactive_parents)
            vapbset = set(vapb)
            upset = set(uactive_parents)
            uapbset = set(uapb)
            
            if len(vapbset.intersection(vpset))==0:
                uactive_parents = []
                vactive_parents = []
                #uapb = []
            if len(uapbset.intersection(upset))==0:
                vactive_parents = []
                uactive_parents = []
                #vapb = []

            """end change
            """
            
            uactive_child_bounds_index = uchild.bounds.contains(u)
            ulen_acbi = len(uactive_child_bounds_index)
            assert(len(uactive_child_bounds_index)<2),'error: child bounds active in two spans at once'
            
            
            
            vactive_child_bounds_index = vchild.bounds.contains(v)
            vlen_acbi = len(vactive_child_bounds_index)
            assert(len(vactive_child_bounds_index)<2),'error: child bounds active in two spans at once'
            
            
            
            if ulen_acbi==1 and vlen_acbi==1:
                uchild_bounds = uchild.bounds[uactive_child_bounds_index[0]]
                uinner = uchild.enclosed_basis_in_range(uchild_bounds)
                uouter = uchild.parent.enclosed_basis_in_range(uchild_bounds) #to be eliminated
                
            else:
                uchild_bounds = dummy_bounds
                uinner = []
                uouter = [] #none to be eliminated
                # if none are encl
                #uapb=[]
                
            
            if vlen_acbi==1 and ulen_acbi==1:
                vchild_bounds = vchild.bounds[vactive_child_bounds_index[0]]
                vinner = vchild.enclosed_basis_in_range(vchild_bounds)
                vouter = vchild.parent.enclosed_basis_in_range(vchild_bounds) #to be eliminated
                
            else:
                vchild_bounds = dummy_bounds
                vinner = []
                vouter = []
                #vapb = []
            
        
            if ulcurve.rm is None:
                ulcurve.rm = ulcurve.wgen(ulcurve.p,ulcurve.level+1)#
        
            if vlcurve.rm is None:
                vlcurve.rm = vlcurve.wgen(vlcurve.p,vlcurve.level+1)#
            
            
            
            
            ulevel_mats.append(ulcurve.rm)
            vlevel_mats.append(vlcurve.rm)
            
            
            
#            cl = len(ubasis_levels)
#            matrix_reducer = np.matmul(ulevel_mats[-1].T,uiBasis)# + np.matmul(vlevel_mats[-1].T,viBasis)
#            #matrix_reducer = np.matmul(vlevel_mats[-1].T,matrix_reducer)
#            for i in range(cl-1):
#                cBasis = ubasis_levels[i]
#                local_reducer = np.matmul(ulevel_mats[i].T,cBasis)
#                for j in range(i+1,cl):
#                    local_reducer = np.matmul(ulevel_mats[j].T,local_reducer)
#                    #local_reducer = np.matmul(vlevel_mats[j].T,local_reducer)
#                matrix_reducer = matrix_reducer + local_reducer
#            uoBasis[uapb] = uoBasis[uapb] - matrix_reducer[uapb]
            
            
            iBasis = uinactive_basis_levels[-1]
            uoBasis = np.matmul(ulevel_mats[-1].T,iBasis)
            
            uactiveBasis = ulcurve.clip_basis_outside_refinement(uoBasis, 
                                                    uactive_parents)
            
            
            
            
#            cl = len(vbasis_levels)
#            matrix_reducer = np.matmul(vlevel_mats[-1].T,viBasis) #+ np.matmul(ulevel_mats[-1].T,uiBasis)
#            #matrix_reducer = np.matmul(ulevel_mats[-1].T,matrix_reducer)
#            for i in range(cl-1):
#                cBasis = vbasis_levels[i]
#                local_reducer = np.matmul(vlevel_mats[i].T,cBasis)
#                for j in range(i+1,cl):
#                    #transBasis = basis_levels[j]
#                    local_reducer = np.matmul(vlevel_mats[j].T,local_reducer)
#                    #local_reducer = np.matmul(ulevel_mats[j].T,local_reducer)
#                matrix_reducer = matrix_reducer + local_reducer
#            voBasis[vapb] = voBasis[vapb] - matrix_reducer[vapb]
            
            
            iBasis = vinactive_basis_levels[-1]
            voBasis = np.matmul(ulevel_mats[-1].T,iBasis)
            
            vactiveBasis = vlcurve.clip_basis_outside_refinement(voBasis, 
                                                    vactive_parents)
            
            ##TLM additional fix?
            #may need further tweaking because it needs to look at all levels
#            if len(vapb)==0 or len(vactive_parents)==0 or (voBasis == 0.).all():
#                #pass
#                #uoBasis[uouter] = 0.
#                uoBasis[:] = 0.
#                #voBasis[vouter] = 0.
#            else:
#                #voBasis[vouter] = 0.
#                pass
                #uoBasis[uouter] = 0.
            vactiveBasis[vouter] = 0.
            ##
            ##
            ####
            ##
            ##
            ##TLM additional fix?
            #may need further tweaking because it needs to look at all levels
#            if len(uapb)==0 or len(uactive_parents)==0 or (uoBasis == 0.).all():
#                #pass
#            #    uoBasis[uouter] = 0.
#                voBasis[:] = 0.
#            else:
#                #uoBasis[uouter] = 0.
#                pass
                #voBasis[uouter] = 0.
                #voBasis[vouter] = 0.
            #ubasis_levels.append(uoBasis)
            uactiveBasis[uouter] = 0.
            
            
#            if (uoBasis == 0.).all():
#                voBasis = self.clip_basis_outside_refinement(voBasis, 
#                                                    [])
#            elif (voBasis == 0.).all():
#                uoBasis = self.clip_basis_outside_refinement(uoBasis, 
#                                                    [])


            """-------------------------------------------------------------"""
            #if (uoBasis == 0.).all():
            #    vactiveBasis[:] = 0.
                #assert((voBasis == 0.).all()),'uo: u = {}, v={}'.format(u,v)
                #voBasis = self.clip_basis_outside_refinement(voBasis, 
                #                                    [])
            #elif (voBasis == 0.).all():
            #    uactiveBasis[:] = 0.
                #assert((uoBasis == 0.).all()),'vo: u = {}, v={}'.format(u,v)
                #uoBasis = self.clip_basis_outside_refinement(uoBasis, 
                #                                    [])
            
#            usave = uactiveBasis.copy()
#            vsave = vactiveBasis.copy()
#            uactiveBasis = uactiveBasis*sum(vsave)
#            vactiveBasis = vactiveBasis*sum(usave)
            
            uinactive_Basis = uoBasis - uactiveBasis
            vinactive_Basis = voBasis - vactiveBasis
            
            vbasis_levels.append(vactiveBasis)
            ubasis_levels.append(uactiveBasis)
            uinactive_basis_levels.append(uinactive_Basis)
            vinactive_basis_levels.append(vinactive_Basis)
            
            
        return ubasis_levels, vbasis_levels
    
    def clip_basis_outside_refinement(self, basis, interval):
        save = np.zeros_like(basis)
        save[interval] = copy.copy(basis[interval])
        basis = basis*0.
        basis[interval] = save[interval]
        return basis
    
    def set_row_column_to_basis_range(self):
        """
           |-+-+-+-+..
           | | | | |
        u_i|.+-+-+-+.     <- this row at u_i is a longitudinal 
           |..-+-+-+         at the ith TCURVE.greville_abscissa(i)
           | +
           | |
           __+_+_+_+
           
                ^
                |
                v_j   This column at v_j is a transvserse
                        at the jth LCURVE.greville_abscissa(j)
        """
        self.transverse_curve_to_ga = []
        self.longitudinal_curve_to_ga = []
        for j in range(self.nv):
            self.transverse_curve_to_ga.append(
                    self.mastercurve_v.greville_abscissa(j))
        for i in range(self.nu):
            self.longitudinal_curve_to_ga.append(
                    self.mastercurve_u.greville_abscissa(i))
        return
    
    def fix_influence_matrix(self):
        """self.Xch : characteristic matrix
        """
        self.set_row_column_to_basis_range()
        self.Xch = np.zeros((self.nu,self.nv),bool) 
        self.Xch_inactive = np.zeros((self.nu,self.nv),bool) 
        for j in range(self.nv):
            v = self.transverse_curve_to_ga[j]
            for i in range(self.nu):
                u = self.longitudinal_curve_to_ga[i]
                #
                bounds = self.bounds.get_box_containing(u,v)
                if bounds.isempty:
                    self.Xch[i,j] = False
                    self.Xch_inactive[i,j] = True
                else:
                    self.Xch[i,j] = True
                    self.Xch_inactive[i,j] = False
        return
    
    def project_matrix_to_finest(self,matrix):
        """THB projection of a surface matrix
        independent of u,v
        to finest level
        
        Problem here:
        Can no longer just find a curve's worth of active inactive
        
        you must operate on the matrix of indices at once.
        """
        this = self.get_coarsest_level()
        this.fix_influence_matrix() # don't actually use this one
        
        P = matrix.copy()
        
        #P[this.Xch_inactive] = 0. #avoid truncation below? -no.
        while this is not None:
            
            if this.parent is None:
                pass
            else:
                
                this.fix_influence_matrix() 
                
    
                intermediate_dv = np.zeros_like(this.vertices)
                final_dv        = np.zeros_like(this.vertices)
                
                
                # PROJECTION: parent vertices to the child vertices: 
                #**************************************************************
                for i in range(this.parent.nv): #loop over the u curves
                    intermediate_dv[:,i] = np.matmul(this.parent.trm,
                                                       P[:,i]) 
                for j in range(this.nu): #loop over the v curves
                    #print 'j = ',j
                    final_dv[j,:] = np.matmul(this.parent.lrm,
                                        intermediate_dv[j,:this.parent.nv])
                #**************************************************************
                
                lP = final_dv.vertices.copy()  
                # Truncation:
                #**************************************************************
                #
                # (inactive this level):
                lP[this.Xch_inactive] = 0. 
                #
                # (parents, excluded on this level's active interval):
                final_dv[this.Xch] = 0.
                #**************************************************************
                
                P = lP + final_dv
                
                self.assign_bounds_to_master_curves()
                
            this = this.child
        return P
    
    def project_all_to_finest(self):
        """THB projection of a surface matrix
        independent of u,v
        to finest level
        
        Problem here:
        Can no longer just find a curve's worth of active inactive
        
        you must operate on the matrix of indices at once.
        """
        this = self.get_coarsest_level()
        this.fix_influence_matrix() # don't actually use this one
        
        P = this.vertices.copy()
        
        #P[this.Xch_inactive] = 0. #avoid truncation below? -no.
        while this is not None:
            
            if this.parent is None:
                pass
            else:
                
                this.fix_influence_matrix() 
                
                lP = this.vertices.copy()  
                
                # Truncation 
                # (inactive this level):
                lP[this.Xch_inactive] = 0.
                #
                this.assign_bounds_to_master_curves()
    
                intermediate_dv = np.zeros_like(this.vertices)
                final_dv        = np.zeros_like(this.vertices)
                
                
                # PROJECTION: parent vertices to the child vertices: 
                #**************************************************************
                for i in range(this.parent.nv): #loop over the u curves
                    intermediate_dv[:,i] = np.matmul(this.parent.trm,
                                                       P[:,i]) 
                for j in range(this.nu): #loop over the v curves
                    #print 'j = ',j
                    final_dv[j,:] = np.matmul(this.parent.lrm,
                                        intermediate_dv[j,:this.parent.nv])
                #**************************************************************
                
                # Truncation 
                # (parents, excluded on this level's active interval):
                final_dv[this.Xch] = 0.
                
                lP += final_dv#.copy()
                P = lP#.copy() #update P to this level
                
                this.assign_bounds_to_master_curves()
                
            this = this.child
        return P
    
    def eval_smart(self, u,v,P=None):
        """what would be smart is to 
        project only when something changes!
        """
        if P is None:
            P = self.project_all_to_finest()
        this    = self.get_finest_level()
        ucurve,vcurve = this.get_level_master_curves(this.level)
        ubasis  = ucurve.TLMBasisFuns(u)
        vbasis  = vcurve.TLMBasisFuns(v) 
        S = np.dot(ubasis.T, np.dot(vbasis.T,P))
        return S
    
    #2nd version of eval
    def eval(self, u,v):#, ans = 0., level=0):
        """Basic idea of THB surfaces:
            -project all active knots (control vertex, basis function pairs)
            to the highest level of detail
            -evaluate the surface at this level of detail only
            using standard B-spline surface evaluation.
            
        with curves, anything goes.
        
        for the (a?) curve version see 
        def rthbeval(self, u) in 'myspline.py'
        
         unrefined 
         space
         over here
         (this side of Bu_{i+1})
        .           .       .
        |           |       |           refined space
        .           .       .           (this side of Bu{i+1})
        |           |       |
        Bu_{i}      i1      Bu_{i+1}  i2
        
        if i1 is in between two knots points,
        and one is refined and the other is not
        
        i.e. 
        
        
        
        then we are ~not enclosed~ on the refined basis
        take the non-refined answer 
        for
            S = A U R
            we would be in A, not R
            we would be in the unrefined space
        """
        #udepth = self.nHLevels['u'], u,v, level):
        """
        
        >>> gg.surface_point(.5,.5,0)
        array([[ 6.,  0.]])
        
        >>> gg.HL['v'][0]
        <class curve.Bspline>
        >>> gg.HL['v'][0].eval(.5)
        
        array([[ 0.,  6.]])
        >>> gg.HL['u'][0].eval(.5)
        array([[ 6.,  0.]])
        
        """
        #bounds = self.bounds.get_box_containing(u,v)
        self.assign_bounds_to_master_curves()
        
        #dv = {}
        #rep_dv = {}
        #finelevel = self.get_finest_level()
        
        this = self.get_coarsest_level()
        ucurve,vcurve = this.get_level_master_curves(this.level)
        count = 0
        P = 0.
        lP = 0.
        #tol = self.tol
        tol = 1.e-2
        #P = this.vertices.copy() 
        while this is not None:
            
            
            if this.parent is None:
                P = this.vertices.copy()  #vertices from the coarse level
                pass#lP = this.vertices.copy()  #?? TLM Nov added - not needed
            else:
                lP = this.vertices.copy()  #what the heck??
            
            
                #print 'count = ',count
                #v[count] = this.vertices.copy() 
                ulcurve,vlcurve = this.get_level_master_curves(this.level)
                
                pkg = this.get_active_inactice_indices(u,v)
                uactive     = pkg[0]
                uinactive   = pkg[2]
                vactive     = pkg[1]
                vinactive   = pkg[3]
                #self.assign_bounds_to_master_curves()
                
                #"""#****************************
#                for i in range(this.nv):  #loop over the u curves
#                    lP[uinactive,i] = 0.
#                for j in range(this.nu): #loop over the v curves
#                    lP[j,vinactive] = 0.


                for i in vinactive:  #loop over the u curves
                    lP[uinactive,i] = 0.
                    #print 'level ',this.level
                    #print 'lp : ',lP
                    self.lP = lP
                for j in uinactive: #loop over the v curves
                    lP[j,vinactive] = 0.
                    
                    
                self.assign_bounds_to_master_curves()
    
                    
                #"""
                ##
                ##********************************
                ##********************************
                ##    project the basis up one level:
                ##   (to the next finer level)
                intermediate_dv = np.zeros_like(this.vertices)
                final_dv        = np.zeros_like(this.vertices)
                
                ulcurve,vlcurve    = this.get_level_master_curves(this.level)
                #uactive, uinactive = ulcurve.children.get_active_inactice_indices(u)
                #vactive, vinactive = vlcurve.children.get_active_inactice_indices(v)
                
                #silly attempt:
#                pkg = this.parent.get_active_inactice_indices(u,v)
#                iuactive     = pkg[0]
#                iuinactive   = pkg[2]
#                ivactive     = pkg[1]
#                ivinactive   = pkg[3]
                #do not use these to cancel stuff down here
                
                #step up one level:
                pkg = this.get_active_inactice_indices(u,v)
                uactive     = pkg[0]
                uinactive   = pkg[2]
                vactive     = pkg[1]
                vinactive   = pkg[3]
                #self.assign_bounds_to_master_curves()
                
                
                self.store_P = P
                """
                # Q. why parent(?):-
                # A. because we are projecting the parent vertices
                #   to the child vertices 
                #   i.e. (coarse to fine)
                #
                #loop over the original u curve level vertices
                #"""
                for i in range(this.parent.nv): #loop over the u curves
                    intermediate_dv[:,i] = np.matmul(this.parent.trm,
                                                       P[:,i]) 
                #for i in range(this.parent.nv): #loop over the u curves
                #    intermediate_dv[uactive,i] = 0.
                #for j in range(this.parent.nu): #loop over the v curves
                #    final_dv[j,ivactive] = 0.
                    
                #intermediate_dv[uactive,i] = 0. 
                #now you have to extend vertices in to additional x directed extent
                #so that the y direction has somehting to operate on
                
                self.store_lP = lP
                self.store_P = P
                self.store_intermediate_dv = intermediate_dv
                self.store_final_dv = final_dv
                self.store_nu = this.parent.nu
                self.store_nv = this.parent.nv
                self.store_lrm = this.parent.lrm
                self.store_trm = this.parent.trm
                
                #loop over the fine level v curve vertices
                """switched 
                
                
                    final_dv[j,:] = np.matmul(this.parent.lrm,
                                        intermediate_dv[j,:this.parent.nu])
                
                to
                
                    final_dv[j,:] = np.matmul(this.parent.lrm,
                                        intermediate_dv[j,:this.parent.nv])
                
                remember, u direction is 'vertical' == 'transverse'
                v direction is 'longitudinal' == 'horizontal'
                
                mind your (u,v)'s!
                
                
                """
                for j in range(this.nu): #loop over the v curves
                    #print 'j = ',j
                    final_dv[j,:] = np.matmul(this.parent.lrm,
                                        intermediate_dv[j,:this.parent.nv])
                    #
                    #final_dv[j,:] = np.matmul(this.parent.lrm,
                    #                    intermediate_dv[j,:this.parent.nu])
                #
                #final_dv[uactive,vactive] = 0. #think - what if u is all good
                #but v has some inactives??
                #ergo this immediately above cannot be right
                #"""#****************************supposedly good
#                for i in range(this.nv): #loop over the u curves
#                    final_dv[uactive,i] = 0.
#                for j in range(this.nu): #loop over the v curves
#                    final_dv[j,vactive] = 0.
                """#****************************
                for i in vactive: #loop over the u curves
                    final_dv[uactive,i] = 0.
                for j in uactive: #loop over the v curves
                    final_dv[j,vactive] = 0.
                #"""#****************************
                
                ##ulcurve,vlcurve    = this.get_level_master_curves(this.level)
                #assert(len(ulcurve.bounds.bounds)==1)
                ##assert(len(ulcurve.bounds.bounds)==1),'len = {}'.format(len(ulcurve.bounds.bounds))
                #uinner = ulcurve.enclosed_basis_in_range(ulcurve.bounds)
                #uouter = ulcurve.enclosed_basis_in_range(ulcurve.bounds[0]) 
                #vinner = vlcurve.enclosed_basis_in_range(vlcurve.bounds)
                #assert(len(vlcurve.bounds.bounds)==1)
                ##assert(len(ulcurve.bounds.bounds)==1),'len = {}'.format(len(ulcurve.bounds.bounds))
                #vouter = vlcurve.enclosed_basis_in_range(vlcurve.bounds[0]) 
                
                
                #"""#Succeeds:
                #self.assign_bounds_to_master_curves()#reset bounds to everything
                for j in range(this.nu):  #   uouter: #
                    #final_dv[j,vouter] = 0.
                    final_dv[j,vactive] = 0.
                for i in range(this.nv): #   vouter:  #    #loop over the u curves
                    #final_dv[uouter,i] = 0.
                    final_dv[uactive,i] = 0.
                ###
                lP += final_dv.copy()
                
                P = lP.copy() #update P to this level, after projection and truncation
                
                
                self.assign_bounds_to_master_curves()
                
            ##
            ##********************************##********************************
            ## NExt level
            ##
            #ckS = ckSx+ckSy
            #assert(np.linalg.norm(ckS - S) <.0001),'ERROR: ckS != S::  {} ??= {}'.format(ckS[0], S)
            this = this.child
            count +=1
        


        
        ckS = 0.
        S   = 0.
        this    = self.get_finest_level()
        ucurve,vcurve = this.get_level_master_curves(this.level)
        ubasis  = ucurve.TLMBasisFuns(u)
        vbasis  = vcurve.TLMBasisFuns(v) 
            
        p = ucurve.p
        q = vcurve.p
        
        uspan = ucurve.FindSpan(u)
        vspan = vcurve.FindSpan(v)
        uind  = uspan-p

            
        """  eval the P&T way:
        for l in range(q+1):
            temp = 0.
            vind = vspan-q+l
            for k in range(p+1):  
                temp += ubasis[uind+k]*P[uind+k][vind] 
            S += vbasis[vspan-q+l]*temp
        #"""
            
        #"""just dot it:
        S = np.dot(ubasis.T, np.dot(vbasis.T,P))
        #"""
        
        """The old dydic loft soft check:
        ckS = np.dot(ubasis.T, np.dot(vbasis.T,P))
        assert(np.linalg.norm(ckS - S) <tol),'ERROR: ckS != S::  {} ??= {}'.format(ckS[0], S)
        assert(abs(np.linalg.norm(S-self.child.child.surface_point(u,v))) < tol),\
                'error! eval({},{}) = {} != surface point = {}'.format(u,v,S,self.surface_point(u,v))
        #"""
        return S 
    
    
    #1st version of eval:
    def eval_old(self, u,v):#, ans = 0., level=0):
        """
        THE ISSUE
        
        Is that -active- vs -inactive-
        basis functions on level L
        have to be aware of -active- -inactive-
        functions on u AND v
        at level L+1
        
        """
        """
        Solution: cast them all to the highest level
        applying tuncation along the way
        evaluate the resultant there.
        
        
        for the curve version see 
        def rthbeval(self, u):
            but note that it does not project to
            highest
            curves can safely evaluate on every level.
        
        
        if i1 is in between two knots points,
        and one is refined and the other is not
        then we are ~not enclosed~ on the refined basis
        take the non-refined answer 
        for
            S = A U R
            we would be in A, not R
            we would be in the unrefined space
        """
        #udepth = self.nHLevels['u'], u,v, level):
        """
        
        >>> gg.surface_point(.5,.5,0)
        array([[ 6.,  0.]])
        
        >>> gg.HL['v'][0]
        <class curve.Bspline>
        >>> gg.HL['v'][0].eval(.5)
        
        array([[ 0.,  6.]])
        >>> gg.HL['u'][0].eval(.5)
        array([[ 6.,  0.]])
        
        """
        dv = {}
        #rep_dv = {}
        #finelevel = self.get_finest_level()
        
        this = self.get_coarsest_level()
        count = this.level
        #ucurve,vcurve = this.get_level_master_curves(this.level)

        while this.child is not None:
            count = this.level
            #print 'count = ',count
            dv[count] = this.vertices.copy() 
            #print count
            #""" Could exclude this by moving to full dot products.
            #ulcurve,vlcurve = this.get_level_master_curves(this.level) #deleted TLM July 24 2017

            #uactive, uinactive = ulcurve.get_active_inactice_indices(u)
            #vactive, vinactive = vlcurve.get_active_inactice_indices(v)
            
            pkg = this.get_active_inactice_indices(u,v)
            uactive     = pkg[0]
            uinactive   = pkg[2]
            vactive     = pkg[1]
            vinactive   = pkg[3]
            self.assign_bounds_to_master_curves()
            
            #dv[count][uinactive,:] = 0.
            #dv[count][:,vinactive] = 0.
            #            for i in uinactive:
            #                #for j in vinactive:
            #                dv[count][i,vinactive][:]  = 0.
            #            for j in vinactive:
            #                dv[count][uinactive,j][:] = 0.
            #for i in range(this.nu): #loop over the v curves
            #    dv[count][uinactive,i] = 0.
            #for j in range(this.nv):
            #    dv[count][:,vinactive] = 0.
            #""" 
            #Do this at the next level
            
            #dv[count][uinactive,vinactive] = 0. #think - what if u is all good
            #but v has some inactives??
            #ergo this immediately above cannot be right
            
            #Fails:  yes indeed, it fails!
            """#****************************
            for i in range(this.nv):  #loop over the u curves
                dv[count][uinactive,i] = 0.
            for j in range(this.nu): #loop over the v curves
                dv[count][j,vinactive] = 0.
            #"""#****************************
            
            #SUCCEEDS!:
            #"""#****************************
            for i in vinactive:#loop over the u curves where v is active
                dv[count][uinactive,i] = 0.
            for j in uinactive:#loop over the v curves where u is active
                dv[count][j,vinactive] = 0.
            #"""#****************************
            """
            for i in range(this.nv): 
                dv[count][i,vinactive] = 0.
            for j in range(this.nu):
                dv[count][uinactive,j] = 0.
            #"""
            loc_surf = this
            
            while loc_surf.child is not None:
                #up_one = this.child
                
                
                
                intermediate_dv = np.zeros_like(loc_surf.child.vertices)
                final_dv        = np.zeros_like(loc_surf.child.vertices)
                
                #ulcurve,vlcurve    = loc_surf.get_level_master_curves(loc_surf.level) #deleted TLM july 24 2017
                #uactive, uinactive = ulcurve.children.get_active_inactice_indices(u)
                #vactive, vinactive = vlcurve.children.get_active_inactice_indices(v)
                
                #jgpkg = loc_surf.get_active_inactice_indices(u,v) #fix bounds
                
                pkg = loc_surf.child.get_active_inactice_indices(u,v)
                uactive     = pkg[0]
                vactive     = pkg[1]
                uinactive   = pkg[2]
                vinactive   = pkg[3]
                
                ulcurve,vlcurve    = loc_surf.get_level_master_curves(loc_surf.child.level) #deleted TLM july 24 2017
                
                assert(len(ulcurve.bounds.bounds)==1)
                #uinner = ulcurve.enclosed_basis_in_range(ulcurve.bounds)
                uouter = ulcurve.parent.enclosed_basis_in_range(ulcurve.bounds[0]) 
                #vinner = vlcurve.enclosed_basis_in_range(vlcurve.bounds)
                assert(len(vlcurve.bounds.bounds)==1)
                vouter = vlcurve.parent.enclosed_basis_in_range(vlcurve.bounds[0]) 
                
                
                """ Not needed
                #self.assign_bounds_to_master_curves()#reset bounds to everything
                for j in range(loc_surf.nu):  #   uouter: #
                    dv[count][j,vouter] = 0.
                for i in range(loc_surf.nv): #   vouter:  #    #loop over the u curves
                    dv[count][uouter,i] = 0.
                #"""
                    
                for i in range(loc_surf.nv): #loop over the u curves
                    intermediate_dv[:,i] = np.matmul(loc_surf.trm,
                                                       dv[count][:,i]) 
                #for i in vactive:  #loop over the u curves where v is active - which do not all necessarily exist yet
                #    intermediate_dv[uactive,i] = 0.
                
                
                #intermediate_dv[uactive,i] = 0. 
                #now you have to extend vertices in to additional x directed extent
                #so that the y direction has somehting to operate on
                
                for j in range(loc_surf.child.nu): #loop over the v curves
                    #print 'j = ',j
                    final_dv[j,:] = np.matmul(loc_surf.lrm,
                                                intermediate_dv[j,:loc_surf.nv])
                #                    final_dv[j,:] = np.matmul(loc_surf.lrm,
                #                                                intermediate_dv[j,:loc_surf.nu])
                    
                    
                
                assert(len(ulcurve.bounds.bounds)==1)
                #uinner = ulcurve.enclosed_basis_in_range(ulcurve.bounds)
                uouter = ulcurve.enclosed_basis_in_range(ulcurve.bounds[0]) 
                #vinner = vlcurve.enclosed_basis_in_range(vlcurve.bounds)
                assert(len(vlcurve.bounds.bounds)==1)
                vouter = vlcurve.enclosed_basis_in_range(vlcurve.bounds[0]) 
                
                
                #"""#Succeeds:
                #self.assign_bounds_to_master_curves()#reset bounds to everything
                for j in range(loc_surf.child.nu):  #   uouter: #
                    final_dv[j,vouter] = 0.
                for i in range(loc_surf.child.nv): #   vouter:  #    #loop over the u curves
                    final_dv[uouter,i] = 0.
                #"""#
                
                """#fails:
                for j in uouter:  
                    final_dv[j,vouter] = 0.
                for i in vouter:   #loop over the u curves
                    final_dv[uouter,i] = 0.
                #"""
                
                    #final_dv[j,vactive] = 0. 
                #for i in range(loc_surf.nu):
                #  final_dv[uactive,i] = 0. 
                   
                #final_dv[uactive,:] = 0. 
                #final_dv[:,vactive] = 0. 
                #for i in uactive:
                #    #final_dv[i,vactive][:] = 0.
                #    for j in vactive:
                #        final_dv[i,j][:] = 0.
                #for j in vactive:
                #    final_dv[uactive,j][:] = 0.
                """
                # of Lcurve vertices == # of Tcurves
                discard active parts of refined coefficients
                on this particular level
                #"""
                ### let me run this in terminal!
                #
                #
                #final_dv[uactive,vactive] = 0. #think - what if u is all good
                #but v has some inactives??
                #ergo this immediately above cannot be right
                
                #FAILS?:  nope, this actually Succeeds!?? -- and is apparently not needed.
                """#****************************supposedly good
                for i in range(loc_surf.child.nv): #loop over the u curves
                    final_dv[uactive,i] = 0.
                for j in range(loc_surf.child.nu): #loop over the v curves
                    final_dv[j,vactive] = 0.
                #"""#****************************
                
                #Succeeds!: -- and is apparently not needed.
                """
                for i in vactive:  #loop over the u curves where v is active
                    final_dv[uactive,i] = 0.
                for j in uactive:  #loop over the v curves where u is active
                    final_dv[j,vactive] = 0.
                #"""
                ###
                """
                for i in range(loc_surf.child.nv): 
                    final_dv[i,vactive] = 0.
                for j in range(loc_surf.child.nu):
                    final_dv[uactive,j] = 0.
                #"""
                dv[count] = final_dv
                loc_surf = loc_surf.child
                
                
            #ckS = ckSx+ckSy
            #assert(np.linalg.norm(ckS - S) <.0001),'ERROR: ckS != S::  {} ??= {}'.format(ckS[0], S)
            this = this.child
            count +=1
        
        #this = self.get_finest_level()
        assert(this == self.get_finest_level()),'Error this is not the finest level of rep'
        dv[count] = this.vertices.copy() 
        ulcurve,vlcurve = this.get_level_master_curves(this.level)
        #uactive, uinactive = ulcurve.get_active_inactice_indices(u)
        #vactive, vinactive = vlcurve.get_active_inactice_indices(v)
            
        pkg = this.get_active_inactice_indices(u,v)
        uactive     = pkg[0]
        uinactive   = pkg[2]
        vactive     = pkg[1]
        vinactive   = pkg[3]
        self.assign_bounds_to_master_curves()
        
        #dv[count][uinactive,:] = 0.
        #dv[count][:,vinactive] = 0.
        #        for i in uinactive:
        #            #for j in vinactive:
        #            dv[count][i,vinactive][:] = 0.
        #        for j in vinactive:
        #            dv[count][uinactive,j][:] = 0.
        """
            Discard inactive from finest level:
        #"""
        #dv[count][uinactive,vinactive] = 0. #no good - 
        #
        #        #"""#****************************supposedly good
        #        for i in range(this.nv): 
        #            dv[count][uinactive,i] = 0.
        #        for j in range(this.nu):
        #            dv[count][j,vinactive] = 0.
        #        #****************************
        
        
        #"""#succeeds, 
        for i in vinactive:
            dv[count][uinactive,i] = 0.
        for j in uinactive:
            dv[count][j,vinactive] = 0.
        #"""
        #****************************
        """#Fails!
        for i in range(this.nv):
            dv[count][uinactive,i] = 0.
        for j in range(this.nu):
            dv[count][j,vinactive] = 0.
        #"""
        
        
        #****************************
        #"""
        """
        for i in range(this.nv): 
            dv[count][i,vinactive] = 0.
        for j in range(this.nu):
            dv[count][uinactive,j] = 0.
        #"""
        
        ckS = 0.
        S   = 0.
        this    = self.get_finest_level()
        ucurve,vcurve = this.get_level_master_curves(this.level)
        ubasis  = ucurve.TLMBasisFuns(u)
        vbasis  = vcurve.TLMBasisFuns(v) 
            
        p = ucurve.p
        q = vcurve.p
        
        uspan = ucurve.FindSpan(u)
        vspan = vcurve.FindSpan(v)
        uind  = uspan-p
        
        P=0.
        for key in dv:
            P += dv[key]
            
        
        """  eval the P&T way:
        for l in range(q+1): #vcurve.n-1):#
            temp = 0.
            vind = vspan-q+l
            #Sarray[l] = np.dot( uBasisFuns.T, P)
            for k in range(p+1):  #ucurve.n-1):#
                #print 'P[{},{}] = {}'.format(vind,uind+k,P[vind,uind+k])
                temp += ubasis[uind+k]*P[uind+k][vind] #[uind+k][vind] #P[vind][uind+k] #
            S += vbasis[vspan-q+l]*temp
        #"""
        
        #ckSx = ckSx + np.dot(uBasisFuns_lel.T, np.dot(uBasisFuns_lel.T,P))[0][0]#[0,0,0]
        #ckSy = ckSy + np.dot(vBasisFuns_lel.T, np.dot(vBasisFuns_lel.T,P))[0][0]#[0,0,1]
        
        #ckS = np.dot(vbasis.T, np.dot(ubasis.T,P))
        
        #""" just dot it:
        S = np.dot(ubasis.T, np.dot(vbasis.T,P))
        #"""
        
        """ The old dydic loft soft check:
        ckS = np.dot(ubasis.T, np.dot(vbasis.T,P))
        assert(np.linalg.norm(ckS - S) <.0001),'ERROR: ckS != S::  {} ??= {}'.format(ckS[0], S)
        
        assert(abs(np.linalg.norm(S-self.surface_point(u,v))) < .01),\
                'error! eval({},{}) = {} != surface point = {}'.format(u,v,S,self.surface_point(u,v))
        #"""
        return S
    
    
    #1st version eval-like to get vertices as seen by a point:
    def compute_subdivision_vertices(self, u,v):#, ans = 0., level=0):
        """
        """
        dv = {}
        
        this = self.get_coarsest_level()
        
        while this.child is not None:
            count = this.level
            
            #print 'count = ',count
            dv[count] = this.vertices.copy() 
            
            pkg = this.get_active_inactice_indices(u,v)
            uactive     = pkg[0]
            uinactive   = pkg[2]
            vactive     = pkg[1]
            vinactive   = pkg[3]
            
            self.assign_bounds_to_master_curves()
            
            #dv[count][uinactive,:] = 0.
            #dv[count][:,vinactive] = 0.
            #            for i in uinactive:
            #                #for j in vinactive:
            #                dv[count][i,vinactive][:]  = 0.
            #            for j in vinactive:
            #                dv[count][uinactive,j][:] = 0.
            #for i in range(this.nu): #loop over the v curves
            #    dv[count][uinactive,i] = 0.
            #for j in range(this.nv):
            #    dv[count][:,vinactive] = 0.
            #""" 
            #Do this at the next level
            
            #dv[count][uinactive,vinactive] = 0. #think - what if u is all good
            #but v has some inactives??
            #ergo this immediately above cannot be right
            
            #Fails: 
            """#****************************
            for i in range(this.nv):  #loop over the u curves
                dv[count][uinactive,i] = 0.
            for j in range(this.nu): #loop over the v curves
                dv[count][j,vinactive] = 0.
            #"""#****************************
            
            #SUCCEEDS!:
            #"""#****************************
            for i in vinactive:#loop over the u curves where v is active
                dv[count][uinactive,i] = 0.
            for j in uinactive:#loop over the v curves where u is active
                dv[count][j,vinactive] = 0.
            #"""#****************************
            """
            for i in range(this.nv): 
                dv[count][i,vinactive] = 0.
            for j in range(this.nu):
                dv[count][uinactive,j] = 0.
            #"""
            loc_surf = this
            
            while loc_surf.child is not None:
                #up_one = this.child
                
                
                
                intermediate_dv = np.zeros_like(loc_surf.child.vertices)
                final_dv        = np.zeros_like(loc_surf.child.vertices)
                
                #ulcurve,vlcurve    = loc_surf.get_level_master_curves(loc_surf.level) #deleted TLM july 24 2017
                #uactive, uinactive = ulcurve.children.get_active_inactice_indices(u)
                #vactive, vinactive = vlcurve.children.get_active_inactice_indices(v)
                
                #jgpkg = loc_surf.get_active_inactice_indices(u,v) #fix bounds
                
                pkg = loc_surf.child.get_active_inactice_indices(u,v)
                uactive     = pkg[0]
                vactive     = pkg[1]
                uinactive   = pkg[2]
                vinactive   = pkg[3]
                
                ulcurve,vlcurve    = loc_surf.get_level_master_curves(loc_surf.child.level) #deleted TLM july 24 2017
                
                assert(len(ulcurve.bounds.bounds)==1)
                #uinner = ulcurve.enclosed_basis_in_range(ulcurve.bounds)
                uouter = ulcurve.parent.enclosed_basis_in_range(ulcurve.bounds[0]) 
                #vinner = vlcurve.enclosed_basis_in_range(vlcurve.bounds)
                assert(len(vlcurve.bounds.bounds)==1)
                vouter = vlcurve.parent.enclosed_basis_in_range(vlcurve.bounds[0]) 
                
                
                """ Not needed
                #self.assign_bounds_to_master_curves()#reset bounds to everything
                for j in range(loc_surf.nu):  #   uouter: #
                    dv[count][j,vouter] = 0.
                for i in range(loc_surf.nv): #   vouter:  #    #loop over the u curves
                    dv[count][uouter,i] = 0.
                #"""
                    
                for i in range(loc_surf.nv): #loop over the u curves
                    intermediate_dv[:,i] = np.matmul(loc_surf.trm,
                                                       dv[count][:,i]) 
                #for i in vactive:  #loop over the u curves where v is active - which do not all necessarily exist yet
                #    intermediate_dv[uactive,i] = 0.
                
                
                #intermediate_dv[uactive,i] = 0. 
                #now you have to extend vertices in to additional x directed extent
                #so that the y direction has somehting to operate on
                
                for j in range(loc_surf.child.nu): #loop over the v curves
                    #print 'j = ',j
                    final_dv[j,:] = np.matmul(loc_surf.lrm,
                                                intermediate_dv[j,:loc_surf.nu])
                    
                    
                
                assert(len(ulcurve.bounds.bounds)==1)
                #uinner = ulcurve.enclosed_basis_in_range(ulcurve.bounds)
                uouter = ulcurve.enclosed_basis_in_range(ulcurve.bounds[0]) 
                #vinner = vlcurve.enclosed_basis_in_range(vlcurve.bounds)
                assert(len(vlcurve.bounds.bounds)==1)
                vouter = vlcurve.enclosed_basis_in_range(vlcurve.bounds[0]) 
                
                
                #"""#Succeeds:
                #self.assign_bounds_to_master_curves()#reset bounds to everything
                for j in range(loc_surf.child.nu):  #   uouter: #
                    final_dv[j,vouter] = 0.
                for i in range(loc_surf.child.nv): #   vouter:  #    #loop over the u curves
                    final_dv[uouter,i] = 0.
                #"""#
                
                """#fails:
                for j in uouter:  
                    final_dv[j,vouter] = 0.
                for i in vouter:   #loop over the u curves
                    final_dv[uouter,i] = 0.
                #"""
                
                    #final_dv[j,vactive] = 0. 
                #for i in range(loc_surf.nu):
                #  final_dv[uactive,i] = 0. 
                   
                #final_dv[uactive,:] = 0. 
                #final_dv[:,vactive] = 0. 
                #for i in uactive:
                #    #final_dv[i,vactive][:] = 0.
                #    for j in vactive:
                #        final_dv[i,j][:] = 0.
                #for j in vactive:
                #    final_dv[uactive,j][:] = 0.
                """
                # of Lcurve vertices == # of Tcurves
                discard active parts of refined coefficients
                on this particular level
                #"""
                #
                #
                #final_dv[uactive,vactive] = 0. #think - what if u is all good
                #but v has some inactives??
                #ergo this immediately above cannot be right
                
                #FAILS?:  nope, this actually Succeeds!?? -- and is apparently not needed.
                """#****************************supposedly good
                for i in range(loc_surf.child.nv): #loop over the u curves
                    final_dv[uactive,i] = 0.
                for j in range(loc_surf.child.nu): #loop over the v curves
                    final_dv[j,vactive] = 0.
                #"""#****************************
                
                #Succeeds!: -- and is apparently not needed.
                """
                for i in vactive:  #loop over the u curves where v is active
                    final_dv[uactive,i] = 0.
                for j in uactive:  #loop over the v curves where u is active
                    final_dv[j,vactive] = 0.
                #"""
                ###
                """
                for i in range(loc_surf.child.nv): 
                    final_dv[i,vactive] = 0.
                for j in range(loc_surf.child.nu):
                    final_dv[uactive,j] = 0.
                #"""
                dv[count] = final_dv
                loc_surf = loc_surf.child
                
                
            #ckS = ckSx+ckSy
            #assert(np.linalg.norm(ckS - S) <.0001),'ERROR: ckS != S::  {} ??= {}'.format(ckS[0], S)
            this = this.child
            count +=1
        
        #this = self.get_finest_level()
        assert(this == self.get_finest_level()),'Error this is not the finest level of rep'
        dv[count] = this.vertices.copy() 
        ulcurve,vlcurve = this.get_level_master_curves(this.level)
        #uactive, uinactive = ulcurve.get_active_inactice_indices(u)
        #vactive, vinactive = vlcurve.get_active_inactice_indices(v)
            
        pkg = this.get_active_inactice_indices(u,v)
        uactive     = pkg[0]
        uinactive   = pkg[2]
        vactive     = pkg[1]
        vinactive   = pkg[3]
        self.assign_bounds_to_master_curves()
        
        #dv[count][uinactive,:] = 0.
        #dv[count][:,vinactive] = 0.
        #        for i in uinactive:
        #            #for j in vinactive:
        #            dv[count][i,vinactive][:] = 0.
        #        for j in vinactive:
        #            dv[count][uinactive,j][:] = 0.
        """
            Discard inactive from finest level:
        #"""
        #dv[count][uinactive,vinactive] = 0. #no good - 
        #
        #        #"""#****************************supposedly good
        #        for i in range(this.nv): 
        #            dv[count][uinactive,i] = 0.
        #        for j in range(this.nu):
        #            dv[count][j,vinactive] = 0.
        #        #****************************
        
        
        #"""#succeeds, 
        for i in vinactive:
            dv[count][uinactive,i] = 0.
        for j in uinactive:
            dv[count][j,vinactive] = 0.
        #"""
        #****************************
        """#Fails!
        for i in range(this.nv):
            dv[count][uinactive,i] = 0.
        for j in range(this.nu):
            dv[count][j,vinactive] = 0.
        #"""
        
        
        #****************************
        #"""
        """
        for i in range(this.nv): 
            dv[count][i,vinactive] = 0.
        for j in range(this.nu):
            dv[count][uinactive,j] = 0.
        #"""
        
        #        ckS = 0.
        #        S   = 0.
        #        this    = self.get_finest_level()
        #        ucurve,vcurve = this.get_level_master_curves(this.level)
        #        ubasis  = ucurve.TLMBasisFuns(u)
        #        vbasis  = vcurve.TLMBasisFuns(v) 
        #            
        #        p = ucurve.p
        #        q = vcurve.p
        #        
        #        uspan = ucurve.FindSpan(u)
        #        vspan = vcurve.FindSpan(v)
        #        uind  = uspan-p
        #        
        #        P=0.
        #        for key in dv:
        #            P += dv[key]
        #            
        #            
        #        for l in range(q+1): #vcurve.n-1):#
        #            temp = 0.
        #            vind = vspan-q+l
        #            #Sarray[l] = np.dot( uBasisFuns.T, P)
        #            for k in range(p+1):  #ucurve.n-1):#
        #                #print 'P[{},{}] = {}'.format(vind,uind+k,P[vind,uind+k])
        #                temp += ubasis[uind+k]*P[uind+k][vind] #[uind+k][vind] #P[vind][uind+k] #
        #            S += vbasis[vspan-q+l]*temp
        #        #ckSx = ckSx + np.dot(uBasisFuns_lel.T, np.dot(uBasisFuns_lel.T,P))[0][0]#[0,0,0]
        #        #ckSy = ckSy + np.dot(vBasisFuns_lel.T, np.dot(vBasisFuns_lel.T,P))[0][0]#[0,0,1]
        #        
        #        ckS = np.dot(ubasis.T, np.dot(vbasis.T,P))
        #        #ckS = np.dot(vbasis.T, np.dot(ubasis.T,P))
        #        
        #        assert(np.linalg.norm(ckS - S) <.0001),'ERROR: ckS != S::  {} ??= {}'.format(ckS[0], S)
        #        #S[0] = ckSx
        #        #S[1] = ckSy
        #        assert(abs(np.linalg.norm(S-self.surface_point(u,v))) < .01),\
        #                'error! eval({},{}) = {} != surface point = {}'.format(u,v,S,self.surface_point(u,v))
        return dv #S#  ckS# 
    
    def make_thb_equivalent_curves(self):
        """for each level in the THB surface
        make rbspline (thb-spline) curves
        that correspond to each u row and v column 
        of vertices
        """
        surf = self.get_coarsest_level()
        nu = surf.nu
        nv = surf.nv
        surf.equiv_ucurves = []
        surf.equiv_vcurves = []
        for i in range(nu):
            surf.equiv_ucurves.append()
        return
    
    
    def plotpts(self, up=None, vp=None):
        if up is None and vp is None:
            up = self.tcurvelist[0]
            vp = self.lcurvelist[0]
        for el in up:
            el.plotcurve()
        for el in vp:
            el.plotcurve()
        return
    
    def plot_all_levels(self, color_ci=0):
        """TODO:
            redo with surfaces
        """
        ci = color_ci
        for el in self.tcurvelist:
            #umc = self.HL['u'][el]
            #for uel in self.tcurvelist[el]:
                #assert(uel.bounds == umc.bounds),'umc bounds != uel bounds!'
            el.plotcurve_onelevel(color_ci=ci)
        for el in self.lcurvelist:
            #vmc = self.HL['v'][el]
            #for vel in self.lcurvelist[el]:
                #assert(vel.bounds == vmc.bounds),'vmc bounds != vel bounds!'
            el.plotcurve_onelevel(color_ci=ci)
            
        if self.child is not None:
            ci = (ci+1)%3
            self.child.plot_all_levels(color_ci=ci)
        return
    
    def plot_THBcurves(self, level=None):
        """TODO:  make interior real curve
        aware of their THB basis!
        
        mapping problem:  dyadic refinement in u
        gives new curves in v that don't belong to 
        previous curves in v
        
        give them a u location based on the knot loc
        in the u space, the knot that they -represent-?
        """
        #if level is None:
        #    el = self.level
        #umc = self.HL['u'][el]
        for el in self.tcurvelist:
            #for uel in self.tcurvelist[el]:
            #assert(uel.bounds == umc.bounds),'umc bounds != uel bounds!'
            el.plotcurve()
        #el = level
        #vmc = self.HL['v'][el]
        for el in self.lcurvelist:
            #for vel in self.lcurvelist[el]:
            #assert(vel.bounds == vmc.bounds),'vmc bounds != vel bounds!'
            el.plotcurve()
        return

    def get_active_vertices_per_level(self, u,v,surf):
        """
        TODO: uses bounds[0]
        
        """
        #urange = surf.mastercurve_u.bases(
        #        surf.mastercurve_u.active(u))
        #vrange = surf.mastercurve_v.bases(
        #        surf.mastercurve_v.active(v))
        active_parent_bounds_index = surf.mastercurve_u.bounds.contains(u)
        len_apbi = len(active_parent_bounds_index) 
        if len_apbi==1:
            parent_bounds = surf.mastercurve_u.bounds[active_parent_bounds_index[0]]
        else:
            parent_bounds = surf.mastercurve_u.bounds[0]       
        urange = surf.mastercurve_u.enclosed_basis_in_range(
                parent_bounds)
#                surf.mastercurve_u.bounds[0])
        
        
        
        active_parent_bounds_index = surf.mastercurve_v.bounds.contains(v)
        len_apbi = len(active_parent_bounds_index)
        if len_apbi==1:
            parent_bounds = surf.mastercurve_v.bounds[active_parent_bounds_index[0]]
        else:
            parent_bounds = surf.mastercurve_v.bounds[0]
        vrange = surf.mastercurve_v.enclosed_basis_in_range(
                parent_bounds)
#                surf.mastercurve_v.bounds[0])
        return urange,vrange
    
    
    def get_active_vertices(self, u,v,indices=None, surf=None):
        if surf is None:
            indices={}
            surf = self.get_finest_level()
            indices[surf.level] = surf.get_active_vertices_per_level(u,v,surf)
            return self.get_active_vertices(u,v,indices, surf.parent)
        else:
            indices[surf.level] = surf.get_active_vertices_per_level(u,v,surf)
            if surf.parent is not None:
                return self.get_active_vertices(u,v,indices, surf.parent)
            else:
                return indices
            
            
            
    def get_surfaces_by_level(self, level=None):
        if level is None:
            surf = self.get_coarsest_level()
            return surf
        else:
            return self.level
    def get_surfacelist(self):
        surf_list = []
        surf = self.get_coarsest_level()
        surf_list.append(surf)
        while surf.child is not None:
            surf = surf.child
            surf_list.append(surf)
        return surf_list
    
    def establish_bounds(self, level=None):
        surf = self.get_coarsest_level()
        surf.mastercurve_u.bounds = BoundsList([surf.bounds[0][0]])
        surf.mastercurve_v.bounds = BoundsList([surf.bounds[0][1]])
        while surf.child is not None:
            surf = surf.child
            surf.mastercurve_u.bounds = BoundsList([surf.bounds[0][0]])
            surf.mastercurve_v.bounds = BoundsList([surf.bounds[0][1]])
        return
    
    
    
    #    def update_tcurves(self,Rs):
    #        nnet = []
    #        for curve in self.tcurvelist:
    #            vertices = np.dot(curve.vertices,Rs)
    #            curve.vertices = vertices
    #            nnet.append(curve)
    #        tcurvelist = nnet
    #        return tcurvelist
    #    
    #    def update_lcurves(self,Rs):
    #        nnet = []
    #        for curve in self.lcurvelist:
    #            vertices = np.dot(curve.vertices[:],Rs)
    #            curve.vertices = vertices
    #            nnet.append(curve)
    #        lcurvelist = nnet
    #        return lcurvelist
    
    
    
    def translate_tcurves(self,Ts):
        nnet = []
        for curve in self.tcurvelist:
            vertices = curve.vertices.copy()
            vertices[:,0] = curve.vertices[:,0] + Ts[0]
            vertices[:,1] = curve.vertices[:,1] + Ts[1]
            vertices[:,2] = curve.vertices[:,2] + Ts[2]
            curve.vertices = vertices
            nnet.append(curve)
        tcurveNet = nnet
        return tcurveNet
    
    def translate_lcurves(self,Ts):
        nnet = []
        for curve in self.lcurvelist:
            vertices = curve.vertices.copy()
            vertices[:,0] = curve.vertices[:,0] + Ts[0]
            vertices[:,1] = curve.vertices[:,1] + Ts[1]
            vertices[:,2] = curve.vertices[:,2] + Ts[2]
            curve.vertices = vertices
            nnet.append(curve)
        lcurveNet = nnet
        return lcurveNet
    
    def rotate_tcurves(self,Rs):
        nnet = []
        for curve in self.tcurvelist:
            vertices = np.dot(curve.vertices,Rs)
            curve.vertices = vertices
            nnet.append(curve)
        tcurveNet = nnet
        return tcurveNet
    
    def rotate_lcurves(self,Rs):
        nnet = []
        for curve in self.lcurvelist:
            vertices = np.dot(curve.vertices[:],Rs)
            curve.vertices = vertices
            nnet.append(curve)
        lcurveNet = nnet
        return lcurveNet
    
    def rotate_vertices(self,Rs):
        return np.dot(self.vertices[:,:],Rs[:])
    def translate_vertices(self,Ts):
        vertices = self.vertices.copy()
        vertices[:,:,0] = self.vertices[:,:,0] + Ts[0]
        vertices[:,:,1] = self.vertices[:,:,1] + Ts[1]
        vertices[:,:,2] = self.vertices[:,:,2] + Ts[2]
        return vertices
        
    
    def update_loft(self,Rs=None,Ts=None):
        """TODO:
            make this more efficient
            by similar methods 
            as was done/started for B-spline curves
            i.e. don't instantiate all over again!
        """
        if Rs is not None:
            tcurveNet = self.rotate_tcurves(Rs)
            lcurveNet = self.rotate_lcurves(Rs)
            verts = self.rotate_vertices(Rs)
        if Ts is not None:
            tcurveNet = self.translate_tcurves(Ts)
            lcurveNet = self.translate_lcurves(Ts)
            verts = self.translate_vertices(Ts)
        if Ts is None and Rs is None:
            tcurveNet = self.tcurveNet
            
        return THBsurface(tcurvelist = tcurveNet,
                          ulevel = self.ulevel,
                          lcurvelist = lcurveNet,
                          vlevel = self.vlevel,
                          vertices = verts,
                          trm = self.trm,
                          lrm = self.lrm)
        
    
    def plotVertices(self, level, active_dict, ax,
                     annotate=False,
                     colorize_levels=False):
        """
        inputs
        ----------
        ax : canvas
        
        
        dev
        ----------
        
        self = SD.THBhull
        
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        import thbsurface as thbspline
        
        import numpy as np
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        
        
        
        self.establish_bounds()
        #active_dict = self.get_active_vertices(ui,vi)
        #
        #*******************************************************
        # get the active basis all at once, to avoid 
        # massive duplication in plotting points
        ads = self.compute_total_active_by_level()
        #
        active_dict = {}
        for level in ads:
            this = ads[level]
            active_dict[level] = (this.uactive,
                                   this.vactive)
            
        
            
        fig = plt.figure(figsize=(9.5,5.0))
        rect = fig.add_subplot(1, 1, 1).get_position()
        ax = Axes3D(fig, rect)
        biggest = max(self.get_extremes())
        ax.set_xlim3d(-biggest,biggest)
        ax.set_ylim3d(-biggest,biggest)
        ax.set_zlim3d(-biggest,biggest)
        
        level = 0
        color = False
        """
        #for key in active_dict:
        key = level
        surf = self.index[key]
        active_u, active_v = active_dict[key]
        vertices = surf.vertices
        
        
        if colorize_levels:
            colors = {0:'red',
                  1:'green',
                  2:'blue',
                  3:'white'}
        else:
            colors = {0:'.65',
                      1:'.55',
                      2:'.28',
                      3:'.1'}
            
        cc = colors[key%3]
        for i in active_u:
            #connects all active vertices in the v direction
            #issue - possible to connect across disjoint portions of the space.
            ax.plot(vertices[i,active_v,2], 
                    vertices[i,active_v,0],
                    vertices[i,active_v,1],
                    marker = "o",  
                    alpha = .5,
                    color=cc)
        
        for i in active_v:
            #this just connects all active vertices in the u direction
            ax.plot(vertices[active_u,i,2], 
                    vertices[active_u,i,0],
                    vertices[active_u,i,1],
                    marker = "o",  
                    alpha = .5,
                    color=cc)
        return ax
    
    def plotonecurve(self, curve, level, 
                     canvas,
                     greville_abscissa,
                     isu=False,isv=False, 
                     nump=30, dim=3,
                     colorize_levels=False,
                     color_override=None,
                     plot_vertices=False):
        """
        greville_abscissa : the parametric location in the perpindicular 
        parametric direction at which this curve 'lives'
        
        
        greville_abscissa = ga
        
        level = surf.level
        
        TODO:  make consistent plotting of curves
        between u and v directions.
        Right now it is overly exhuberent in plotting 
        vertices though (seems) correct(????) 
        in plotting curves?
        """
        if colorize_levels:
            colors = {0:'purple',
                  1:'green',
                  2:'blue',
                  3:'white'}
        else:
            colors = {0:'.65',
                      1:'.55',
                      2:'.28',
                      3:'.1'}
        cc = colors[level%3]
        if color_override is not None:
            cc = color_override
        
        if isu:
            Direction = 'u'
        elif isv:
            Direction = 'v'
            
        surf = self.get_surface_by_level(level)
        
        ilist = surf.get_active_disjoint_bounds(direction=Direction,
                                                ploc = greville_abscissa)
        
        if len(ilist)>0:
            for intval in ilist:
                a = intval[0]
                b = intval[1]
                s = np.linspace(a,b,nump,endpoint=True)
                cvec = np.zeros((nump,curve.dim))
                for i,u in enumerate(s):
                    #this check must be done outside 
                    # to ensure a consistent plotting interval 
                    # once we get to this point
                    #ck = self.child.bounds.get_box_containing()
                    #if ck.isempty:
                    cvec[i] = curve.CurvePoint(u) #assumes that this is a THB curve
                canvas.plot(cvec[:,2], 
                            cvec[:,0],
                            cvec[:,1],
                            linestyle = "-",
                            color = cc)
                if plot_vertices:
                    
                    # plot all active:
                    vi = curve.bases(curve.active_basis(intval))
                    # enclosed only:
                    #vi = curve.enclosed(intval)
                    
                    cpts = curve.vertices[vi]
                    canvas.plot(cpts[:,2], 
                                cpts[:,0],
                                cpts[:,1],
                                marker = "o",  
                                alpha = .5,
                                color=cc)
                    #except:
                    #    pass
        return canvas
    
    def get_active_disjoint_bounds(self, direction,ploc,nump=30):
        """returns the disjoint set
        of active intervals 
        in ONE DIRECTION (u or v)
        on ONE LEVEL of the surface
        for one parametric location in the other direction
        
        
        assumptions
        ----------
        -nested spaces!
        
        
        inputs:
        ----------
        direction : this is the direction 
                    in which a curve or ray runs.
                    --------------------------
                    type: string, either u, or v
                    you can also pass an int, 0 for u or 1 for v.
                    
                    (this direction becomes dt)
                    
             ploc : real between 0 and 1, inclusive
                    This is the location in_the_other_direction
                    at which you wish to have the disjoint active intervals
                    for this level
                    
            ploc = ga #(greville_abscissa)
        
        
        returns:
        ----------
        directed bounds list:
        whatever boxes the ray along ploc intersects, take the portions
        of those boxes that are active on this level only,
        and return the direction, dt for those boxes
        i.e.
        return the intervals of the active boxes 
        that the ray passes through on this level.
            
        
        dev
        ----------
        Box( ia(u,u) , ia(v,v) )
        BoxList([Box1, Box2,.... BoxN])
        
        direction = 'x'
        ploc = .44
        nump = 30
        
        BoxList = thbspline.BoxList
        box = thbspline.Box
        
        """
        if direction == 'u' or direction == 0:
            dt = 0
        else:
            assert(direction ==  'v' or direction == 1),'ERROR: bad direction for THB parametric bounding'
            dt = 1
            
        surf    = self
            
        ray     = np.linspace(0.,1.,nump)
        #give me the boxes at this level 
        # which can possibly be relevant:
        this_bdl = []
        val = 1 - dt
        for r in ray:
            # run along direction dt
            # keep the other direction fixed at ploc
            ckr = [0.,0.]
            ckr[val] = ploc
            ckr[dt] = r
            #i checked: the code below is backwards:
            #            ckr[dt] = ploc
            #            ckr[val] = r
            bx = surf.bounds.get_box_containing(*ckr)
            if not bx.isempty:
                if this_bdl:
                    pre_box = this_bdl[-1]
                    ckand = pre_box & bx
                    if not ckand.isempty:
                        this_bdl[-1] = pre_box | bx
                    else:
                        this_bdl.append(bx)
                    #new = [ s1 & s2 for s2 in this_bdl ]
                else:
                    this_bdl.append(bx)
                    
        #take just the [dt] direction for the bounds:
        this_directedbounds = [ bds[dt] for bds in this_bdl]
        
        child = surf.child
        if child is not None:
            child_directedbounds = [ bds[dt] for bds in child.bounds if bds[val].contains(ploc)]
            
            if child_directedbounds and this_directedbounds:
                this_directedbounds = iafuncs.interval_list_diff(this_directedbounds,
                                                                 child_directedbounds)
                
        
        return this_directedbounds
    
    
    
    
    

    def plotSurface(self, 
                    innerBox=None, 
                    OuterBox=None, 
                    limx=None,limy=None,limz=None,
                    minx=0.,miny=0.,minz=0.,
                    Rstride=100,
                    Cstride=100,
                    Shade=False,
                    colorize_vertices = False,
                    colorize_curves = False,
                    show_hullcurves=True,
                    show_THBcurves = True,
                    show_THBvertices = True,
                    fancy=True,
                    start_hullcurve=1,
                    simplefancy=False,
                    alpha=.75,
                    canvas=None):
        """
            Plot a surface
            
            Note:  self.compute_surface() uses eval to 
                populate the self.surface array
        
        
        Developer testing:
        --------------------
            self = SD.THBhull
            
            #---------------------
            import utility_optimization as uopt
            
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
            import thbsurface as thbspline
                    
            import numpy as np
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
            
            limx = None
            
            Shade = True
            Cstride = 100
            Rstride = 100
            colorize_vertices = False
            colorize_curves = False
            
            show_THBcurves = True
            show_hullcurves = True
            
            fancy = True
            
            nump=30
            minx=miny=minz=-100.
            limx=limy=limz=100.
            #-------------------------
            
            
            
            fig = plt.figure(figsize=(9.5,5.0))
            rect = fig.add_subplot(1, 1, 1).get_position()
            ax = Axes3D(fig, rect)
            
            
            bounds = thbspline.BoxList(  
                                thbspline.Box( ia(.0,.8125),
                                               ia(.0,.34375) ) )
            
            bounds = thbspline.BoxList(  
                                thbspline.Box( ia(.0,.8125),
                                               ia(.0,1.) ) )
            
            bounds = thbspline.BoxList(
                            thbspline.Box( ia(.0,.8125),   
                                               ia(.0,.34375) ),
                            thbspline.Box( ia(0.,.4375),      
                                           ia(.34375,1.) )
                            )
                            
            DO the second set of bounds activate?
            let's find out:
                
            a=SD.THBhull.compute_true_nonzero_basis(0.01,.99)
            
            L2 bounds:
                
            bounds = thbspline.BoxList(  
                                thbspline.Box( ia(.0,.8125),
                                               ia(.0,1.) ) )
            
            a[2].uactive
            Out[35]: [0, 1, 2, 3]
            
            a[2].vactive
            Out[36]: [31, 32, 33, 34]
            
            
            L2 bounds:
                
            SD.THBhull.child.child.bounds
            Out[48]: BoxList((Box((ia(0.0, 0.8125), ia(0.0, 0.34375))), 
                       Box((ia(0.0, 0.4375), ia(0.34375, 1.0)))))

            
            Good result:
                a[2].uactive
                Out[50]: [0, 1, 2, 3]
                
                a[2].vactive
                Out[51]: [31, 32, 33, 34]
            
            plot thickens:
                SD.THBhull.child.child.mastercurve_u.bounds
                Out[53]: BoundsList([ia(0.0, 0.8125)])
                
                SD.THBhull.child.child.mastercurve_v.bounds
                Out[54]: BoundsList([ia(0.0, 1)])
            
        """
        #
        #
        def aaa():
            from matplotlib.colors import LogNorm
            fig = plt.figure(figsize=(9.5,5.0))
            rect = fig.add_subplot(1, 1, 1).get_position()
            ax = Axes3D(fig, rect)
            #
            biggest = max(self.get_extremes())
            ax.set_xlim3d(-biggest*.0,1.*biggest)
            ax.set_ylim3d(-biggest*.5,.5*biggest)
            ax.set_zlim3d(-biggest*.5,.5*biggest)
            return ax
        #ax = aaa()
        #
        #
        """ #The wireframe algorithm is not prepared
        #to handle complicated spline surfaces
         --do not use  --->     ax.plot_wireframe(self.surface[:,:,0], 
                                                  self.surface[:,:,1], 
                                                  self.surface[:,:,2], 
                                                  rstride=3, 
                                                  cstride=3, 
                                                  color = 'b',
                                                  alpha=.8)
        #"""
        #""" #traditional
        NiceCmap='viridis'
        NiceCmap='binary'
        
        if canvas is None:
            ax = aaa()
        else:
            ax = canvas
        
        if not fancy:
            ax.plot_surface(self.surface[:,:,2],  
                            self.surface[:,:,0], 
                            self.surface[:,:,1], 
                            rstride=Rstride, 
                            cstride=Cstride, 
                            color='0.75',
                            alpha=alpha,
                            shade=Shade)
            
        elif simplefancy:
            ax.plot_surface(self.surface[:,:,2],  
                            self.surface[:,:,0], 
                            self.surface[:,:,1], 
                            rstride=1, 
                            cstride=1, 
                            color='0.75',
                            alpha=alpha,
                            shade=False)
                    
            if show_hullcurves:
                for curve in self.hullcurvelist[start_hullcurve:]:  
                    ax.plot(curve.r[:,2], 
                            curve.r[:,0], 
                            curve.r[:,1], 
                            label='parametric curve', 
                            color = 'r')
            return ax
        else:
            ax.plot_surface(self.surface[:,:,2],  
                            self.surface[:,:,0], 
                            self.surface[:,:,1], 
                            rstride=1, 
                            cstride=1, 
                            color='0.75',
                            alpha=alpha,
                            shade=False)
            
            
        
        

        #for vi in self.v:
        #    for ui in self.u:
        #    self.assign_bounds_to_master_curves(ui,vi)
        """
        s0.mastercurve_u.enclosed_basis_in_range(
                                                ia(0.2,.8)
                                                )
        #"""
        self.establish_bounds()
        #active_dict = self.get_active_vertices(ui,vi)
        #
        #*******************************************************
        # get the active basis all at once, to avoid 
        # massive duplication in plotting points
        ads = self.compute_total_active_by_level()
        #
        active_dict = {}
        for level in ads:
            this = ads[level]
            active_dict[level] = (this.uactive,
                                   this.vactive)
            
            
            
        if show_THBcurves:
                
                
            for level in active_dict:
                surf = self.get_surface_by_level(level)
                for ci, curve in enumerate(surf.lcurvelist):  
                    ga = surf.mastercurve_u.greville_abscissa(ci)
                    ax = surf.plotonecurve(curve,
                                           level,
                                           ax,
                                           ga,
                                           isv=True,
                                           colorize_levels=colorize_curves,
                                           plot_vertices = show_THBvertices)
                for ci, curve in enumerate(surf.tcurvelist):  
                    ga = surf.mastercurve_v.greville_abscissa(ci)
                    ax = surf.plotonecurve(curve,
                                           level,
                                           ax,
                                           ga,
                                           isu=True,
                                           colorize_levels=colorize_curves,
                                           plot_vertices = show_THBvertices)
                    
        if show_hullcurves:
            #start_hullcurve says
            # how far aft to start this list 
            # (probably just input the info)
            for curve in self.hullcurvelist[start_hullcurve:]:  
                ax.plot(curve.r[:,2], 
                        curve.r[:,0], 
                        curve.r[:,1], 
                        label='parametric curve', 
                        color = 'r')
            """
            if self.flats_outline_exists:
                ax.plot(self.flats_outline.r[:,2], 
                        self.flats_outline.r[:,0], 
                        self.flats_outline.r[:,1], 
                        label='parametric curve', 
                        color = '.1')
                #ax = self.flats_outline.plot3D_hull_system(
                #                            canvas=ax,color='black')
            else:
                try:
                    self.issue_flats_outline()
                    ax.plot(self.flats_outline.r[:,2], 
                            self.flats_outline.r[:,0], 
                            self.flats_outline.r[:,1], 
                            label='parametric curve', 
                            color = '.1')
                    #ax = self.flats_outline.plot3D_hull_system(
                    #                            canvas=ax,color='black')
                except:
                    print 'Warning, exception with curvature outline'
            #"""
            
        plt.show()
        return ax
    
    
    
    def issue_flats_outline(self):
        vertices = []
        for curve in self.hullcurvelist[:3]:
            vertices.append(curve.vertices[2])
        for curve in self.hullcurvelist[3:5]:
            vertices.append(curve.vertices[3])
        for curve in self.hullcurvelist[5:]:
            vertices.append(curve.vertices[2])
        vertices = np.asarray(vertices)
        self.flats_outline = spline.Bspline(vertices, 
                                            k=self.ku,
                                            nump=self.numpu)
        self.flats_outline_exists = True
        return
    
    def compatibility(self,surf=None):
        if surf is None:
            surf = self
        surf.tcurveNet = surf.tcurvelist
        surf.lcurveNet = surf.lcurvelist
        surf.surfNet = surf.vertices
        if surf.child is not None:
            surf.compatibility(surf.child)
        return
    
    def get_extremes(self):
        minx = self.tcurvelist[0].vertices[0,0]
        maxx = minx
        miny = self.tcurvelist[0].vertices[0,1]
        maxy = miny
        minz = self.tcurvelist[0].vertices[0,2]
        maxz = minz
        
        xtuple1 = [el for el in self.mastercurve_u.extremes1d(0)]
        xtuple2 = [el for el in self.mastercurve_v.extremes1d(0)]
        minx = min( [minx] + xtuple1 + xtuple2)
        maxx = max( [maxx] + xtuple1 + xtuple2)
        
        ytuple1 = [el for el in self.mastercurve_u.extremes1d(1)]
        ytuple2 = [el for el in self.mastercurve_v.extremes1d(1)]
        miny = min( [miny] + ytuple1 + ytuple2)
        maxy = max( [maxy] + ytuple1 + ytuple2)
        
        ztuple1 = [el for el in self.mastercurve_u.extremes1d(2)]
        ztuple2 = [el for el in self.mastercurve_v.extremes1d(2)]
        minz = min( [minz] + ztuple1 + ztuple2)
        maxz = max( [maxz] + ztuple1 + ztuple2)
        return minx,maxx,miny,maxy,minz,maxz
    
    
    def IssueRhinoCurveVertices(self, 
                                verticeslist, 
                                the_filename, 
                                usage='a'):
        
        with open(the_filename, usage) as f:
            for el in verticeslist:
                for s in el:
                    st = '{}, {}, {}\n'.format(s[0],s[1],s[2])
                    f.write(st)
                f.write('new curve   \n')
        return
    
    
    def IssueRhinoCurveKnots(self, 
                             knotslist, 
                             the_filename,
                             usage='a'):
        
        with open(the_filename, usage) as f:
            for el in knotslist:
                
                st = ''
                for s in el:
                    st = st+' '+str(s)
                f.write(st)
                f.write('\n new knot vector \n')
                
        return
    
    def IssueRhinoTcurves(self, 
                          prefix='defaultsurf_',
                          the_filename='transverse'):
        tlist = []
        tknots = []
        for curve in self.tcurvelist:
            tlist.append([ list(pt)  for pt in curve.vertices ] )
            tknots.append(list(curve.t[1:-1])) #Rhino style knots.
            
        the_filename = the_filename.split('.')[0]
            
        self.IssueRhinoCurveVertices(verticeslist = tlist,
                    the_filename = prefix+the_filename+'_curves.txt',
                                     usage = 'w')
        self.IssueRhinoCurveKnots(knotslist = tknots,
                    the_filename = prefix+the_filename+'_knots.txt',
                                  usage = 'w')
        return
    
    
    def IssueRhinoLcurves(self, 
                          prefix='defaultsurf_',
                          the_filename='longitudinal'):
        llist = []
        lknots = []
        for curve in self.lcurvelist:
            llist.append([ list(pt)  for pt in curve.vertices ] )
            lknots.append(list(curve.t[1:-1]))
            
        the_filename = the_filename.split('.')[0]
            
        self.IssueRhinoCurveVertices(verticeslist = llist,
                    the_filename = prefix+the_filename+'_curves.txt',
                                     usage = 'w')
        self.IssueRhinoCurveKnots(knotslist = lknots,
                    the_filename = prefix+the_filename+'_knots.txt',
                                  usage = 'w')
        return
    
    
    def IssueRhinoHcurves(self, 
                          prefix='defaultsurf_',
                          the_filename='hullfp'):
        hlist = []
        hknots = []
        for curve in self.hullcurvelist:
            hlist.append([ list(pt)  for pt in curve.vertices ] )
            hknots.append(list(curve.t[1:-1]))
            
        the_filename = the_filename.split('.')[0]
            
        self.IssueRhinoCurveVertices(verticeslist = hlist,
                    the_filename = prefix+the_filename+'_curves.txt',
                                     usage = 'w')
        self.IssueRhinoCurveKnots(knotslist = hknots,
                    the_filename = prefix+the_filename+'_knots.txt',
                                  usage = 'w')
        return
    
    
    def IssueRhinoSurfaceVertices(self, 
                                  verticeslist, 
                                  the_filename,
                                  usage='a'):
        
        if os.path.exists(the_filename):
            usage = 'a' # append if already exists
        else:
            usage = 'w' # make a new file if not
            
        with open(the_filename, usage) as f:
            f.write('\n new vertices   \n')
            for el in verticeslist:
                #f.write('\n new uRow   \n')
                for s in el:
                    st = '{} {} {}\n'.format(s[0],s[1],s[2])
                    f.write(st)
        return
    
    def IssueRhinoSurfaceKnots(self, 
                               uknots,
                               vknots, 
                               the_filename,
                               usage='a'):
        
        if os.path.exists(the_filename):
            usage = 'a' # append if already exists
        else:
            usage = 'w' # make a new file if not
        
        with open(the_filename, usage) as f:
            f.write('\n new uknot vector \n')
            st = ''
            for s in uknots[1:-1]:
                st = st+' '+str(s)
            f.write(st)
            f.write('\n new vknot vector \n')
            st = ''
            for s in vknots[1:-1]:
                st = st+' '+str(s)
            f.write(st)
        return
    
    
    def IssueRhinoSurface(self, the_filename='surface_data.txt'):
        #if os.path.exists(the_filename):
            
        f = open(the_filename, 'w')
        f.write('\n new surface \n')
        #
        #**********************************************************************
        #
        vertices = self.project_all_to_finest()
        surf = self.get_finest_surface()
        #
        #**********************************************************************
        #
        cpts = []
        for i in range(surf.nu):
            cpts.append( list(vertices[i,:]) )
        #
        #**********************************************************************
        #
        f.write('\n point count \n')
        pc = np.shape(vertices)
        p_c = str(pc[0]) +' '+ str(pc[1])
        f.write(p_c)
        #
        #**********************************************************************
        #
        f.close()
        #
        #**********************************************************************
        #
        self.IssueRhinoSurfaceVertices(verticeslist = cpts,
                                       the_filename = the_filename,
                                       usage = 'a')
        #
        #**********************************************************************
        #
        self.IssueRhinoSurfaceKnots(uknots = surf.uknots,
                                    vknots = surf.vknots,
                                    the_filename = the_filename,
                                    usage = 'a')
        #
        #**********************************************************************
        #
        f = open(the_filename, 'a')
        f.write('\n degree u \n')
        f.write(str(surf.ku-1))
        f.write('\n degree v \n')
        f.write(str(surf.kv-1))
        f.write('\n END of Surface \n')
        f.close()
        #
        #**********************************************************************
        # Curves!
        self.IssueRhinoTcurves(prefix = 'complex_Tcurve',
                               the_filename = the_filename)
        self.IssueRhinoLcurves(prefix = 'complex_Lcurve',
                               the_filename = the_filename)
        self.IssueRhinoHcurves(prefix = 'complex_Hullcurve',
                               the_filename = the_filename)
        return 
    

def check_uv(u,v):
    print '\n\n----------------------------'
    print 'Check_uv ',u,v
    print '\nsurface_point:'
    print s0.surface_point(u,v)
    print 'eval:'
    print s0.eval(u,v)
    #array([[[ 8.4,  3.6,  0. ]]])
    print '\nu:' 
    print s0.mastercurve_u.rthbeval(u)
    #print s1.mastercurve_u.rthbeval(u)
    #print s2.mastercurve_u.rthbeval(u)
    print '\nv:'
    print s0.mastercurve_v.rthbeval(v)
    #print s1.mastercurve_v.rthbeval(v)
    #print s2.mastercurve_v.rthbeval(v)
    #array([[ 3.6,  0. ,  0. ]])
    #>>>array([[ 0. ,  8.4,  0. ]])
    return

def check_bounds(self):
    if self.child is None:
        print self.bounds
        print self.mastercurve_u.bounds
        print self.mastercurve_v.bounds
        return 
    else:
        print self.bounds
        print self.mastercurve_u.bounds
        print self.mastercurve_v.bounds
        return check_bounds(self.child)
    
def check_parents(s1,s0):
    assert(s1.parent==s0),'Fail: s1.parent != s0'
    assert(s1.mastercurve_u.parent==s0.mastercurve_u),'Fail: s1.mastercurve_u.parent != s0.mastercurve_u'
    assert(s1.mastercurve_v.parent==s0.mastercurve_v),'Fail: s1.mastercurve_v.parent != s0.mastercurve_v'

    
    assert(s1.mastercurve_u==s0.mastercurve_u.children),'Fail: s1.mastercurve_u != s0.mastercurve_u.children'
    assert(s1.mastercurve_v==s0.mastercurve_v.children),'Fail: s1.mastercurve_v != s0.mastercurve_v.children'
    print 'check complete'
    return

def print_array(a):
    for el in a:
        print el
    return



def test_2scale_basis():
    h2 = s2.mastercurve_u.TLMBasisFuns(.6)
    h1 = s1.mastercurve_u.TLMBasisFuns(.6)
    h1_from_h2 = np.dot(s1.trm.T,h2)
    h2_from_h1 = np.dot(s1.trm,h1) #
    print h1_from_h2
    print h2_from_h1
    return

def test_2scale_basisderivatives():
    h2 = s2.mastercurve_u.TLMDBasisFuns(.6)[0]
    h1 = s1.mastercurve_u.TLMDBasisFuns(.6)[0]
    h1_from_h2 = np.dot(s1.trm.T,h2)
    h2_from_h1 = np.dot(s1.trm,h1) #
    print h1_from_h2
    print h2_from_h1
    return
    


if __name__ == "__main__":
    checkall = False
    s0 = THBsurface(verbose=True)
    #s0 = THBsurface(verbose=True,ku=4,kv=7)
    s0 = s0.dyadic_loft_refinement()
    s0 = s0.dyadic_loft_refinement()
    s0 = s0.dyadic_loft_refinement()
    s0.set_top()
    #s0.numpu = 60
    #s0.numpv = 60
    print 'right now we assume at least one level of refinement'
    #s1 = s0.dyadic_loft_refinement()
    s1 = s0.dyadic_loft_refinement(bounds = BoxList(  
            Box( ia(.125,.875), ia(.125,1.) ) )    )
    
    
    ub1 = s1.uknots[s1.ku]
    ub2 = s1.uknots[s1.n-1]
    #s1 = s0.dyadic_loft_refinement(bounds = BoxList(  
    #        Box( ia(.4,.6), ia(.4,.6) ) )    )
    #s1 = s0.dyadic_loft_refinement()
    #s2 = s1.dyadic_loft_refinement(bounds = BoxList(  
    #        Box( ia(.21,.79), ia(.21,.79) ) )    )
    #s1.parent = s0
    #s0.child = s1
    #s1.set_parent()
    #s0.set_child()
    #s0.set_top()
    #s1.parent = s0
    #s0.child = s1
    #TODO: make this get the right master curves, etc.
    #assert(s0.child  == s1),'failed assertion: s0.child  == s1 ?'
    #assert(s1.parent == s0),'failed assertion: s1.parent == s0'
    
    """
    OLD STUFF
    
    s2 = s1.dyadic_loft_refinement(bounds = BoxList(  
            Box( ia(.37,.63), ia(.37,.63) ) )    )
    
    """
    
    
    """
    NEW STUFF
    
    
    s2 = s1.dyadic_loft_refinement(bounds = BoxList(  
            Box( ia(.2,.4), ia(.2,.4) ),
            Box( ia(.6,.75), ia(.5,.75) ))    
            )
    
    
    
    """
    s2 = s1.dyadic_loft_refinement(bounds = BoxList(  
            Box( ia(.0,.3), ia(.0,.3) ),
            Box( ia(.4,.6), ia(.4,.6) ),
            Box( ia(.7,.9), ia(.8,.9) ))    
            )
    
    #s2 = s1.dyadic_loft_refinement(bounds = BoxList(  
    #        Box( ia(.4,.6), ia(.4,.6) ) )    )
    #s2 = s1.dyadic_loft_refinement(bounds = BoxList(  
    #        Box( ia(.40625,.6875), ia(.40625,.6875) ) )    )
    #s2 = s1.dyadic_loft_refinement(bounds = BoxList(  
    #        Box( ia(.3,.7), ia(.3,.7) ) )    )
    #assert(s1.child  == s2),'failed assertion: s1.child  == s2 ?'
    #assert(s0.child.child  == s2),'failed assertion: s0.child.child  == s2 ?'
    #s3 = s2.dyadic_loft_refinement(bounds = BoxList(  
    #        Box( ia(.4,.6), ia(.4,.6) ) )    )
    
    #"""
    self = s0
    #self.assign_bounds_to_master_curves()
    
    check_uv(.3,.7)
    check_uv(.5,.5)
    check_uv(.7,.3)
    #check_uv(.1,.1)
    #check_uv(.2,.2)
    check_uv(.4,.4)
    
    check_uv(.5,1.)
    #"""
    #check_parents(s1,s0)
    #check_parents(s2,s1)
    
    
    ui = .5
    vi = .5
    dd=self.get_active_vertices(ui,vi)
    print dd
    
    s2.vertices
    #s0.compute_surface()
    #s0.plotSurface()
    self = s0
    
    
    
    #"""
    print 'minimalistically modify the THBsurface at one point!'
    u=.5
    v=.5
    this    = self.get_finest_level()
    ucurve,vcurve = this.get_level_master_curves(this.level)
    ubasis  = ucurve.TLMBasisFuns(u)
    vbasis  = vcurve.TLMBasisFuns(v) 
    
    avb = self.compute_active_by_level(.5,.5)
    #av = self.compute_apb(.5,.5) #do not use this one
    av = self.compute_true_nonzero_basis(.5,.5) 
    print 'level 3:'
    av[3].print_()
    print 'level 4:'
    av[4].print_()
    print 'level 5:'
    av[5].print_()
    dv=  self.compute_subdivision_vertices(.5,.5) #useful for checking eval
    #"""
    #"""
    print 'NOTE:'
    print 'compute_apb picks up apb that are'
    print 'seemingly active down in the lower levels'
    print 'but not really active! -so dont use this kind of thing'
    P = dv[3]
    ckS_0 = np.dot(ubasis.T, np.dot(vbasis.T,P))
    P = dv[4]
    ckS_1 = np.dot(ubasis.T, np.dot(vbasis.T,P))
    P = dv[5]
    ckS_2 = np.dot(ubasis.T, np.dot(vbasis.T,P))
    """
    for level in av:
        this_apo = av[level]
        surf = this_apo.surf
        #        if surf == s0 or surf == s1:
        #            pass
        #        else:
        ulevel = this_apo.uactive
        vlevel = this_apo.vactive
        #for el in av[5][0]: #loop over u curves
        for el in ulevel: #loop over u curves
            #print s2.vertices[el,av[5][1]][:,2] #plot v group
            #s2.vertices[el,av[5][1][0]:av[5][1][-1]+1][:,2] = 5. # v group
            surf.vertices[el,vlevel[0]:vlevel[-1]+1][:,2] = 5. # v group
            
        #for el in av[5][1]: #loop over v curves
        for el in vlevel: #loop over u curves
            #print s2.vertices[av[5][0],el][:,2] #plot u group
           # s2.vertices[av[5][0][0]:av[5][0][-1]+1,el][:,2] = 5. # u group
            surf.vertices[ulevel[0]:ulevel[-1]+1,el][:,2] = 5. # u group
    #"""
    s0.compute_surface()
    s0.plotSurface()
    print 'THBsurface modification complete!'
    #"""
    """
    s2.vertices[av[5][0],:][:,2] = 20.
    s2.vertices[av[5][1],:][:,2] = 20.
    
    s0.vertices[:,:][:,2] = 0.
    s0.vertices[:,:][:,2] = 0.
    s1.vertices[:,:][:,2] = 0.
    s1.vertices[:,:][:,2] = 0.
    s2.vertices[:,:][:,2] = 0.
    s2.vertices[:,:][:,2] = 0.
    s0.compute_surface()
    s0.plotSurface()
    #"""
    u=.5
    v=1.
    
    u=.5
    v=.25
    check_uv(.25,.5)
    u=.25
    v=.5
    
    #u=.4
    #v=.4
    
    u=.3
    v=.7
    
    
    """check the difference between B-spline and THBspline
    #suface evaluation for a given THBsurface called hull:
    #"""
    
    """
    s1.bounds = bounds = BoxList(  
            Box( ia(0.,1.), ia(0.,1.) ) )  
    s2.bounds = BoxList(  
            Box( ia(0.,1.), ia(.0,1.) ) )  
            
            
            
    s1.bounds = bounds = BoxList(  
            Box( ia(0.,1.), ia(0.,1.) ) )  
    s2.bounds = BoxList(  
            Box( ia(0.25,.75), ia(.25,.75) ) )  
    
    
    s1.bounds = bounds = BoxList(  
            Box( ia(.125,.875), ia(.125,1.) ) )
    s2.bounds = BoxList(  
            Box( ia(.125,.3), ia(.125,.3) ),
            Box( ia(.4,.6), ia(.4,.6) ),
            Box( ia(.7,.875), ia(.8,.875) ))  
            
    s1.bounds = bounds = BoxList(  
            Box( ia(.125,.975), ia(.025,1.) ) )
    s2.bounds = BoxList(  
            Box( ia(.125,.3), ia(.125,.3) ),
            Box( ia(.4,.6), ia(.4,.6) ),
            Box( ia(.7,.875), ia(.8,.875) ))  
    
    
    u=0.9655172413793103 
    
    v=0.1724137931034483
    
    u = 0.7931034482758621 
    v = 0.9655172413793103
    s0.eval_smart(u,v) - s0.surface_point(u,v)
    """
    if checkall:
        #uopt.THBspline_to_Bspline_ckdiff_evaluator(s0)
        uopt.THBspline_evaluator_diffs(s0)
    
    #    s2.bounds = BoxList(  
    #            Box( ia(.0,.3), ia(.0,.3) ),
    #            Box( ia(.4,.6), ia(.4,.6) ),
    #            Box( ia(.7,.9), ia(.8,.9) ))  
    
    s1.bounds = BoxList(
            Box( ia(0.125, 0.9375), ia(0.0625, .9375) )
            )
    
    s2.bounds = BoxList(  
            Box( ia(.125,.21875), ia(.125,.21875) ),
            Box( ia(.21875,.34375), ia(.21875,.34375) ),
            Box( ia(.34375,.46875), ia(.34375,.46875) ),
            Box( ia(.46875,.59375), ia(.46875,.59375) ),
            Box( ia(0.59375, 0.6875), ia(0.59375, 0.6875) ),
            Box( ia(0.6875, 0.78125), ia(0.6875, 0.78125) ),
            Box( ia(0.78125, 0.9375), ia(0.78125, 0.9375) )
            )
    """
    s2.bounds.contains(.5,.5)
    Out[18]: True
    
    s0.eval(.5,.5)
    Out[19]: array([[[6.       , 6.       , 1.5239735]]])
    
    s0.eval_smart(.5,.5)
    Out[20]: array([[[6., 6., 5.]]])
    """
    
    if checkall:
        uopt.THBspline_evaluator_diffs(s0)
    
    s0.compute_surface()
    s0.plotSurface()
    
    
    #
    #*****************************************
    # simple as 'pi'
    s1.bounds = BoxList(
            Box( ia(.3,.7), ia(.3,.7) ),
            )
    
    s2.bounds = BoxList(  
            Box( ia(.3,.7), ia(.3,.7) ),
            )
     
    
    if checkall:
        uopt.THBspline_evaluator_diffs(s0)
    
    """
    s0.numpv = 30
    s0.numpu = 30
    s0.numpv = 60
    s0.numpu = 60
    s0.compute_surface()
    s0.plotSurface()
    
    #"""
    
    """
        
    s1.bounds = BoxList(
            Box( ia(0.125, 0.9375), ia(0.0625, .9375) )
            )
    
    s2.bounds = BoxList(  
            Box( ia(.125,.21875), ia(.125,.21875) ),
            Box( ia(.21875,.34375), ia(.21875,.34375) ),
            Box( ia(.34375,.46875), ia(.34375,.46875) ),
            Box( ia(0.59375, 0.6875), ia(0.59375, 0.6875) ),
            Box( ia(0.6875, 0.78125), ia(0.6875, 0.78125) ),
            Box( ia(0.78125, 0.9375), ia(0.78125, 0.9375) )
            )
        
        
    s2.bounds.contains(.5,.5)
    Out[18]: False
    
    s0.eval(.5,.5)
    Out[23]: array([[[ 6.00000000e+00,  6.00000000e+00, -1.30104261e-16]]])
    
    s0.eval_smart(.5,.5)
    Out[24]: array([[[6.        , 6.        , 0.13964385]]])
    
    s0.compute_surface()
    s0.plotSurface()
    
    
    s0.numpv = 300
    s0.numpu = 300
    s0.compute_surface()
    s0.plotSurface()
    s0.plotSurface(colorize_curves=True)
    #"""