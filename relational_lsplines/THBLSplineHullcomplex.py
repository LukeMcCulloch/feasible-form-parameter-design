#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 18:19:16 2018

@author: luke

MultiResolution
Truncated
Hierarchical
Basis
Lagrange
Spline


Testing 3 part bare hull (THB redux)

This one does a single surface around the bulb and  fwd bare hull

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
from myspline import rbspline, BoundsList
import wavelet_matrices as wm
#
import utility_optimization as uopt
#
import mx_plots as spfuncs # circle, gfunc etc
initial_curve = uopt.initial_curve
linear_vertices = uopt.linear_vertices
#
k= 4
nump=30
#
import thbsurface as thbspline
#
def make_curve(section):
    interval_data, small = interval_bounds(section)
    FPD = FormParameterDict(section)
    return



def get_tcurves(fname=None):
    """unpickle a set of transverse curve vertices
    (standard B-spline curves)
    """
    if fname is None: fname='transverse_curves'
    with open(fname,'r') as f:
        tcurvenet = pickle.load(f)
    dl = []
    for vertices in tcurvenet:
        dl.append(rbspline(vertices,k,nump) )
    return dl

def get_lcurves(fname=None, postfix=None):
    """unpickle a set of longitudinal curve vertices
    (standard B-spline curves)
    """
    if fname is None: fname='longitudinal_curves'
    if postfix is not None:
        fname='longitudinal_curves'+postfix
    with open(fname,'r') as f:
        lcurvenet = pickle.load(f)
    dl = []
    for vertices in lcurvenet:
        dl.append(rbspline(vertices,k,nump) )
    return dl

def make_vertex_matrix(tcurvenet):
    """
    tvertnet[0] == tcurvenet[0].vertices
    """
    tvertnet = []
    for c in tcurvenet:
        tvertnet.append(c.vertices)
    return np.asarray(tvertnet)






def nthb_vs_nvi(ti,vi):
    """
    ti = thbcurve with n vertices
    vi = full set of transverse control vertices
        that a single longitudinal curve should interpolate
    """
    return ti.n<np.shape(vi)[0]



def evaluate_gradient(Lspline):
    """
    mark gradient indices 
    that are in need or refinement
    """
    
    return 


def issue_longitudinals(tcurvenet,
                        k=4,nump=30,
                       refinemax=1,
                       maxlevel=4,
                       normalize=False):
    """THB Optimize longitudinal curves
    
    #
    self = thbcurve
    #
    
    """
    #
    tvertnet = make_vertex_matrix(tcurvenet)
    lcurvenet = []
    curve = tcurvenet[0]
    N = curve.n
    nnew = N-2 #just optimize the interior THB curves?
    #
    #
    i=1
    #
    #
    for i in range(    3 ,4):    #   1,nnew+1):   #    
        #**************************************
        #linear interpolation
        #
        #vanilla_spline = spline.interpolatedBspline(
        #                tvertnet[:,i], k, nump)
        #**************************************
        # simple 4 c.v. starting linear curve
        lvts = linear_vertices(tvertnet[:,i][0],
                               tvertnet[:,i][-1],
                               4)
        vtinterp = tvertnet[:,i] 
        thbcurve = rbspline(lvts,k,nump)
        #
        print 'new curve',i
        lcurvenet.append(uopt.THBsolver(
                                thbcurve,
                                vtinterp,
                                refinemax=refinemax,
                                maxlevel=maxlevel,
                                normalize=normalize))
        
    return lcurvenet




def make_surface(lcurvenet,local_tcurvenet=None,
                 ul=2,vl=3,
                 refine=True):
    """
    lcurvenet = self.lcurvenet
    local_tcurvenet=None
    ul=2  #7 lcurves, tcurve.n=7
    vl=0  #4 tcurves, lcurve.n=4
    
    surfNet = []
    for el in lcurvenet_fwd:
        print np.shape(el.vertices)
        surfNet.append(el.vertices)
    surfNet = np.asarray(surfNet)
    #
    
    
    surfNet = []
    for el in lcurvenet_mid:
        print np.shape(el.vertices)
        surfNet.append(el.vertices)
    surfNet = np.asarray(surfNet)
    #
    art_tnet = []
    for i in range(lcurvenet_mid[0].n):
        art_tnet.append(rbspline(surfNet[:,i],
                                 k=4,
                                 nump=30))
        #
    s0 = thbspline.THBsurface(vertices = surfNet,
                              lcurvelist = lcurvenet_mid,
                              tcurvelist = art_tnet,
                              hullcurvelist=None,
                              ulevel = 0,
                              vlevel = 0)
    from iaBox import BoxList, Box
    """
    surfNet = []
    for el in lcurvenet:
        surfNet.append(el.vertices)
    surfNet = np.asarray(surfNet)
    #
    art_tnet = []
    for i in range(lcurvenet[0].n):
        art_tnet.append(rbspline(surfNet[:,i],
                                 k=4,
                                 nump=30))
    #
    s0 = thbspline.THBsurface(vertices = surfNet,
                              lcurvelist = lcurvenet,
                              tcurvelist = art_tnet,
                              hullcurvelist=local_tcurvenet,
                              ulevel = ul,
                              vlevel = vl)
    
    if refine:
        s1 = s0.dyadic_loft_refinement()
    s0.compute_surface()
    return s0
#
#******************************************************************************
#
# thb solver setup for complex hull form
#
#
#******************************************************************************
#
def setup_transverse_0th_curve(curve):
    return 

#
#******************************************************************************
#
# COMPLEX HULL MAKING BELOW
#
#
#******************************************************************************
#
class Bulb(object):
    def __init__(self):
        self.hullbulb = None
        self.get_curves()
    def get_curves(self,tname=None,lname=None):
        self.tcurvenet = self.get_tcurves()
        self.lcurvenet = self.get_lcurves()
        return
    
    def get_tcurves(self, fname=None):
        """unpickle a set of transverse curve vertices
        (standard B-spline curves)
        """
        if fname is None: fname='transverse_bulb_curves'
        with open(fname,'r') as f:
            tcurvenet = pickle.load(f)
        dl = []
        for vertices in tcurvenet:
            dl.append(rbspline(vertices,k,nump) )
        return dl
    
    def get_lcurves(self, fname=None, postfix=None):
        """unpickle a set of longitudinal curve vertices
        (standard B-spline curves)
        """
        if fname is None: fname='longitudinal_bulb_curves'
        if postfix is not None:
            fname='longitudinal_curves'+postfix
        with open(fname,'r') as f:
            lcurvenet = pickle.load(f)
        dl = []
        for vertices in lcurvenet:
            dl.append(rbspline(vertices,k,nump) )
        return dl
    
    
    def make_surface(self, ul=2,vl=1):
        """see opt_simple_hull.py
        function: def make_THB_bulb
        """
        self.hullbulb = make_surface(self.lcurvenet,
                                self.tcurvenet,
                                ul=ul,vl=vl,
                                refine=False)
        return
    
    
    
    def return_tcurve(self,index):
        return self.tcurvenet[index]
    
    def return_outline_vertices(self,num2use=3):
        """
        bub tcurve[0] is aft
        tcurve[0].vertices[0] is keel
        """
        #if isinstance(use, int):
        #    use = [use]
        #else:
        #    use.sort()
        #------
        # if we are using 2 curves, use the fwd most 2,
        # which translates into the 'last' 2 tcurves in the bbow net
        #"""
        use = [-el-1 for el in range(num2use)]
        #----
        fmc = self.return_tcurve(use.pop(0))
        fmv = fmc.vertices
        #----
        vnettop = []
        vnetbot = []
        #
        use.reverse() 
        for el in use:  #loop from bow ([-1]) to aft [0]
            vnetbot.append(self.tcurvenet[el].vertices[0])
        use.reverse()
        for el in use:
            vnettop.append(self.tcurvenet[el].vertices[-1])
            
        outline = []
        for v in vnetbot:
            outline.append(v)
        for v in fmv:
            outline.append(v)
        for v in vnettop:
            outline.append(v)
        outline = np.asarray(outline)
        return outline
        #"""
        ## different idea below, instead of outlining all,
        # just split the fwd tcurve around the bbow.
        #return copy.deepcopy(self.tcurvenet[0].vertices)
    
    def mask_verts(self, num2use=3):
        mask=[]
        mask = np.array(mask, dtype=bool)
        return mask
    
    def return_fwd_boundary(self):
        return
    
    
    def translate(self, vector):
        """
        
        vector = Vector_hull_bulb
        """
        self.hullbulb = self.hullbulb.translate_surface(
                vector[0],vector[1],vector[2])
        self.hullbulb.compute_surface()
        return
        
    
    def plot(self,
             fancyish   = False,
             color      = False,
             simplefancy= False,
             canvas     = None): 
        if canvas is None:
            ax = self.hullbulb.plotSurface(
                            start_hullcurve=0,
                            fancy=fancyish,
                            colorize_curves=color,
                            colorize_vertices=color,
                            simplefancy=simplefancy)
        else:
            ax = self.hullbulb.plotSurface(
                            start_hullcurve=0,
                            fancy=fancyish,
                            colorize_curves=color,
                            colorize_vertices=color,
                            simplefancy=simplefancy,
                            canvas = canvas)
        return ax
    
    def plotbulbcurves(self,
                       plottcurves = True,
                       plotlcurves = False,
                       canvas = None):
        
        if canvas is None:
            cc = self.tcurvenet[0]
            ax = cc.plotcurvehierarchy()#lots of mess lying around... like this
        else:
            ax = canvas
            
        if plottcurves:
            for curve in self.tcurvenet:
                    ax.plot(curve.r[:,2], 
                        curve.r[:,0], 
                        curve.r[:,1], 
                        label='parametric curve', 
                        color = 'r')
        if plotlcurves:
            for curve in self.lcurvenet:
                    ax.plot(curve.r[:,2], 
                        curve.r[:,0], 
                        curve.r[:,1], 
                        label='parametric curve', 
                        color = 'r')
        return ax
        
        

class BareHull(object):
    def __init__(self):
        self.lcurvenet_fwd = None
        self.lcurvenet_mid = None
        self.lcurvenet_aft = None
        self.tcurvenet_fwd = None
        self.tcurvenet_mid = None
        self.tcurvenet_aft = None
        #
        self.hullfwd = None
        self.hullmid = None
        self.hullaft = None
        
    def make(self):
        """make a 3 part bare hull from
        a set of longitudinals conforming to 
        the 3 aft, 2 mid, 3 fwd spec.
        """
        tcurvenet = get_tcurves()
        self.hullcurves = tcurvenet
        c = tcurvenet[0]
        #
        original_lcurves = get_lcurves()
        self.original_lcurves=original_lcurves
        #"""
        # LAUNCH THE MULTIGRID SOLVER:
        """
        dev
        ----------
        cpk=original_lcurves[0]
        dwl=original_lcurves[-1]
        min_fwd_work=False
        
        k=4
        nump=30
        maxlevel=4
        normalize=False
        
        tcurve_flat_ini_index=3
        tcurve_flat_fin_index=4
        Longindex=2
        
        from utility_optimization import thbmatch2D
        """
        lc = uopt.issue_split_longitudinals(tcurvenet,
                                            cpk=original_lcurves[0],
                                            dwl=original_lcurves[-1],
                                            min_fwd_work=False)
        #
        #
        #fwd:
        lcurvenet_fwd = lc[0]#[curve for curve in lc[0]]
        tcurvenet_fwd = tcurvenet[:4]
        #mid:
        lcurvenet_mid = lc[1]
        #lcurvenet_mid = [curve.dyadic_refinement() for curve in lc[1]]
        #lcurvenet_mid = [curve.dyadic_refinement() for curve in lcurvenet_mid]
        tcurvenet_mid = None
        #aft:
        lcurvenet_aft = lc[2]#[curve for curve in lc[2]]
        tcurvenet_aft = tcurvenet[4:]
        #
        sfwd = make_surface(lcurvenet_fwd,
                            tcurvenet_fwd,
                            ul=2,vl=2)
        """
        sfwd.lcurvelist[-1]
        Out[38]: <class curve.Bspline>
        
        sfwd.lcurvelist[-1](locs[0])
        Out[39]: array([[12.35339247, 17.33641345, 17.70896871]])
        
        lcurvenet_fwd[-1](locs[0])
        Out[40]: array([[12.35339247, 17.33641345, 17.70896871]])
        
        lcurvenet_fwd[-1].CurvePoint(locs[0])
        Out[41]: array([12.3705051 , 17.33641345, 17.47426893])
        
        tcurvenet_fwd[1].vertices[-1]
        Out[42]: array([12.35339247, 17.33641345, 17.70896871])
        
        its working.
        """
        smid = make_surface(lcurvenet_mid,
                            tcurvenet_mid,
                            ul=2,vl=0)
        saft = make_surface(lcurvenet_aft,
                            tcurvenet_aft,
                            ul=2,vl=2)
        #
        self.lcurvenet_fwd = lcurvenet_fwd
        self.lcurvenet_mid = lcurvenet_mid
        self.lcurvenet_aft = lcurvenet_aft
        self.tcurvenet_fwd = tcurvenet_fwd
        self.tcurvenet_mid = tcurvenet_mid
        self.tcurvenet_aft = tcurvenet_aft
        #
        self.hullfwd = sfwd
        self.hullmid = smid
        self.hullaft = saft
        return
    
    
    
    def plot(self,
             fancyish       = False,
             color          = False,
             simplefancy    = False,
             canvas         = None,
             plothullcurves = True): 
        """
             fancyish       = True
             color          = False
             simplefancy    = True
             canvas         = None
             plothullcurves = True
        """
        
        if canvas is None:
            ax = self.hullfwd.plotSurface(
                            start_hullcurve=0,
                            fancy=fancyish,
                            colorize_curves=color,
                            colorize_vertices=color,
                            simplefancy=simplefancy)
        else:
            ax = self.hullfwd.plotSurface(
                            start_hullcurve=0,
                            fancy=fancyish,
                            colorize_curves=color,
                            colorize_vertices=color,
                            simplefancy=simplefancy,
                            canvas = canvas)
            
        
        ax = self.hullmid.plotSurface(
                        start_hullcurve=0,
                        fancy=fancyish,
                        colorize_curves=color,
                        colorize_vertices=color,
                        simplefancy=simplefancy,
                        canvas = ax)
        
        ax = self.hullaft.plotSurface(
                        start_hullcurve=0,
                        fancy=fancyish,
                        colorize_curves=color,
                        colorize_vertices=color,
                        simplefancy=simplefancy,
                        canvas = ax)
        
        
        if plothullcurves:
            try:
                for curve in self.hullcurves:
                        ax.plot(curve.r[:,2], 
                            curve.r[:,0], 
                            curve.r[:,1], 
                            label='parametric curve', 
                            color = 'r')
            except:
                mssg = 'Warning, hull curves are not attached to '+\
                            'this hull object, and you tried to plot them!'
                print mssg
        return ax
        

class ComplexHull(object):
    def __init__(self):
        self.k = 4
        self.nump=30
        self.barehull       = BareHull()
        self.bulb           = Bulb()
        
        self.hullcurves2    = get_tcurves()
        #
        #**********************************************************************
        # hi res:  solve the bulb longitudinals on the 35x35 level
        self.make_hires = True #  True #  
        self.highres_num = 19 #35
        self.lowres_num = 11 #19
        #---------------------------
        self.highres_num = 11
        self.lowres_num = 7
        #---------------------------
        self.maxlevel_thbsolver = 3
        #
        # how about the starting curves for blending:
        self.makefresh_blends = True #start from scratch, 11 vertices
        # else: use the existing longitudinal
        #
        #**********************************************************************
        # 
        self.num2use    = 4 #number of  transverse bulb curves 
                            #to use for the bulb
                            # use e.g. 2 with 19x19
        #
        #**********************************************************************
        # number of transition curves added to the bulb solve:
        self.num_extra_boundary=0 #3
        #
        #**********************************************************************
        # equality constrained bulb
        self.barehullcurveblend = 'equality'
        self.bulbblend = 'equality'
        #
        #**********************************************************************
        # 
        
        #we could.....
        #self.ulevel_of_blend = 4
        #self.vlevel_of_blend = 4
        #self.nu_blend = 19
        #self.nv_blend = 19
        #dictate the level of the blend by fiat?
        #(not doing this way)
        
        pass
    
    def make_bare_hull_3part(self):
        self.barehull.make()
        return 
    
    
    def make_bulb(self,ul=2,vl=1):
        self.bulb.make_surface(ul=ul,vl=vl)
        return
    
    def make_system0(self, bulblevels = None):
        if bulblevels is None:
            bulblevels = [0,2]
            
        self.make_bare_hull_3part() #make the bare hull
        self.make_bulb(ul=bulblevels[0],
                       vl=bulblevels[1])
        
        self.add_bulb()
        return
    
    
    def make_system(self, bulblevels = None):
        """
        
        
        dev
        ----------
        self = ch
        """
        if bulblevels is None:
            bulblevels = [0,2]
            
        self.make_bare_hull_3part() #make the bare hull
        
        self.make_bulb(ul=bulblevels[0],
                       vl=bulblevels[1])
        """
        
        self.plot_initial_situation()
        
        """
        
        self.add_bulb2()
        return
    
    
    def add_bulb2(self):
        """
        notes:
        ----------
        *the bare hull tcurves are n=7
        but you've already made a thb surface with 1 level
        refined longitudinals, 
        gauranteeing that the finest tcurves in that 
        surface are n=11
        
        *so dont refine them if the intercept is high.
        
        """
        #"""
        self.move_bulb_to_bare_hull()
        self.split_fwd_tcurve_at_bulb()
        self.refine_to_fair_bulb_to_bare_hull()
        self.make_longitudinals_bulb2()
        #"""
        return
    
    
    def move_bulb_to_bare_hull(self, offset='Not Used At The Moement!'):
        """
        Assumes bbow has already been placed flush
        with the bare hull keel.
        
        dev
        ----------
        
        cc = rbspline( fwd_hull_curve.vertices[:,:2], k=4, nump=30)
        
        tp = fwd_hull_curve.CurvePoint(.5)
        fwd_hull_curve.FindPoint(tp,0)[0]
        
        Status: Done, except perhaps add an offset to smooth things??
        """
        bhull = self.barehull.hullfwd.get_finest_surface()
        #tcurvelist = self.barehull.original_lcurves
        bbow = self.bulb.hullbulb
        
        #
        #*******************************************************
        # intercept bare hull fwd tcurve with bbow top aft.
        bbow_top_aft_pt = bbow.tcurvelist[0].vertices[-1]
        yb = bbow_top_aft_pt[1]
        
        fwd_hull_curve = bhull.tcurvelist[0]
        
        p_intercept = fwd_hull_curve.FindPoint(x=yb,
                                               index=1,
                                               maxiter=200,
                                               tol=1.e-10)
        self.p_intercept = p_intercept
        #
        #*******************************************************
        #
        
        barehull_real_intercept = fwd_hull_curve.CurvePoint(p_intercept[0])
        
        Vector_hull_bulb = barehull_real_intercept - bbow_top_aft_pt
        self.bulb.translate(vector = Vector_hull_bulb)
        #self.bulb.hullbulb
        
        self.pfactor = self.estimate_bulbfractionalheight(p_intercept[0])
        self.barehull_real_intercept = barehull_real_intercept
        #bhull = bhull.parent
        return
    
    def estimate_bulbfractionalheight(self, p_intercept):
        """
        Ensure there are enough vertices to fair in the bulb
        
        * refine_to_fair_bulb_to_bare_hull
            -will use this directly to refine pfactor number of times
            where pfactor is the integer we return here.
            
        status: DONE
        """
        if p_intercept > .75:
            return 0  #the (assumed 7 vertex) bulb is nearly as tall as the hull. only refine once
                        #to capture the 7 vertices + a little extra
        elif p_intercept > .3:
            return 1
        else:
            return 2  #the bulb less than 1/3 height of the bulb
                      #we need more vertices 'down low'
    
    def refine_to_fair_bulb_to_bare_hull(self):
        """
        Ensure there are enough vertices to fair in the bulb
        
        status: DONE
        """
        #bhull = self.barehull.hullfwd
        bhull = self.barehull.hullfwd.get_finest_surface()
        bbow = self.bulb.hullbulb
        
        
        #while bhull.ulevel<(self.pfactor*bbow.ulevel +1) or \
        #                bhull.vlevel<(self.pfactor*bbow.vlevel +1):
        
        tcurvelist = copy.deepcopy(bhull.tcurvelist)
        for i in range(self.pfactor):
            #bhull = bhull.dyadic_loft_refinement()
            for curve in tcurvelist:
                curve = curve.dyadic_refinement()
        self.ini_tcurvelist = []
        for el in tcurvelist:
            cv = el.get_finest_curve()
            cv.parent = None
            self.ini_tcurvelist.append(cv)
        #self.barehull.hullfwd = bhull
        return
    
    
    
    
    
    def split_fwd_tcurve_at_bulb(self):
        """
        Find the location of the intersection
        betwee the fwd bare hull tcurve (refined)
        and the bulbous bow
        
        
        This function assumes more or less, 
        that the incoming tcurves barehull are length 11
        and the incoming tcurfve on the bulb are length 7
        
        It will form an n=11 tcurve at the fwd 
        end that .... does NOT QUITE match with the bulb...
        due to knot vector issues.  That could be fixed.
        in the finest level..
        with more work...
        :(
        
        status:  Done?
        
        
        
        
        dev
        ----------
        
        Q>  Does a curve really keep the same parametric location
        after refinement? 
        
        tc = rbspline(fwd_hull_curve.vertices, k=4,nump=3)
        tc.FindPoint(x=yb,
                    index=1,
                    maxiter=200,
                    tol=1.e-10)
        >>> (0.8501319609422353, array([0.        , 9.93214032, 1.47979534]), True)
        yes, a curve really keeps the same parametric location
        after refinement.
        
        
        
        ax = fwd_hull_curve.plotcurvehierarchy()

        ax = nuc.plotcurvehierarchy(canvas = ax)
        
        Status: eh...
        
        There was a bunching problem:
            -due to the need for a cusp as the upper edge interface
            for the bulb and hull
            - the tcurve[0] must have a cusp condition at 
            whatever level it is defined at.
            -alternative is to break the surface there.
        
        """
        bhull = self.barehull.hullfwd.get_finest_surface()
        bbow = self.bulb.hullbulb
        #
        #*******************************************************
        # intercept bare hull (refined) fwd tcurve with bbow top aft.
        bbow_top_aft_pt = bbow.tcurvelist[0].vertices[-1]
        yb = bbow_top_aft_pt[1]
        
        fwd_hull_curve = bhull.tcurvelist[0]
        fwd_hull_curve.bounds = BoundsList([ia(0.,1.)])
        
        p_intercept = fwd_hull_curve.FindPoint(x=yb,
                                               index=1,
                                               maxiter=200,
                                               tol=1.e-10)
        #
        #*******************************************************
        # Split the fwd bare hull tcurve
        lower_curve, upper_curve = fwd_hull_curve.CurveSplit(p_intercept[0])
        #
        #*******************************************************
        # Copy the aft most bulbous bow tcurve 
        new_bottom_curve = copy.deepcopy(bbow.tcurvelist[0])
        self.len_bulb_boundary = new_bottom_curve.n
        #
        #*******************************************************
        # 4 vertex curve emulating the upper_curve
        vtxup = linear_vertices(upper_curve.vertices[0],
                               upper_curve.vertices[-1],
                               bhull.tcurvelist[0].n-self.len_bulb_boundary 
                               )
        #nuc = rbspline(vtxup, k = self.k, nump = self.nump)
        #-----------------------------------------------------
        #nuc = nuc.dyadic_refinement()
        #-----------------------------------------------------
        #nuc.vertices[1] = nuc.vertices[0]
        #vtxup[1] = vtxup[0]
        #vtxup[2] = vtxup[0]
        #
        #*******************************************************
        # join curves
        """IT WILL CAUSE A vertex eating CUSP
        if you try to enforce strict curve compliance 
        without a real boundary.
        """
        jl = []
        for vtx in new_bottom_curve.vertices:#[:-1]:
            jl.append(vtx)
        for vtx in vtxup: #nuc.vertices[1:]:
            jl.append(vtx)
        jl = np.asarray(jl)
        jcv = rbspline(jl, k=self.k, nump=self.nump)
        
        if False:
            ax = jcv.plotcurvehierarchy()
            for tc in self.bulb.hullbulb.tcurvelist:
                ax = tc.plotcurvehierarchy(canvas=ax)
            
            """
            ax = jcv.plotcurvehierarchy(canvas = ax)
            
            
            #
            bhull = self.barehull.hullfwd.get_finest_surface()
            bbow = self.bulb.hullbulb
            jcv = bhull.tcurvelist[0] 
            #
            
            #"""
        
        while jcv.n < self.barehull.hullfwd.tcurvelist[0].n:
            jcv = jcv.dyadic_refinement()
            
        bhull.tcurvelist[0] = jcv
        return
    
    
    
    
    
    
    
    def make_longitudinals_bulb2(self):
        """
        
        Status: InProgress
        """
        #bhull = self.barehull.hullfwd.get_finest_surface()
        #bbow = self.bulb.hullbulb
        
        tcurvelist = self.ini_tcurvelist
        #
        self.solver_records = {}
        #self.tnum = bhull.tcurvelist[0].n #number of lcurves to solve
        self.tnum = tcurvelist[0].n #number of lcurves to solve (all!)
        #
        #**********************************************************************
        #  check that we have enough vertices everywhere:
        #for i,el in enumerate(bhull.tcurvelist):
        for i,el in enumerate(tcurvelist):
            print el.n, self.tnum
            assert(el.n == self.tnum),'ERROR: bhull.tcurvelist[{}].n = {}, expected = {}'.format(i,el.n ,self.tnum)
        
        
        #
        #**********************************************************************
        # make sure the transverse target curves have vertices 
        # for every longitudinal.
        IntThis = copy.deepcopy(self.barehull.hullcurves[:4])
        newc = []
        for cv in IntThis:
            while cv.n < self.tnum :
                cv = cv.dyadic_refinement()
            newc.append(cv)
        IntThis = newc
        self.fwdtargetcurves = IntThis
        
        
        actually_interp = self.fwdtargetcurves[1:-1]
        for i in range(self.tnum):
            #
            #--------------------------------------------------
            # Here is the new fwd tcurve!:
            #
            #    NEW FWD CURVE HERE:
            #--------------------------------------------------
            #vini = self.barehull.hullfwd.tcurvelist[0].vertices[i]
            #vini = bhull.tcurvelist[0].vertices[i]
            vini = tcurvelist[0].vertices[i]
            #
            #
            #--------------------------------------------------
            vfini = self.fwdtargetcurves[-1].vertices[i]
            #--------------------------------------------------
            
            linverts = linear_vertices(vini,vfini,4)
            curve = rbspline(linverts,k=4,nump=30)
            curve = curve.dyadic_refinement() #5
            curve = curve.dyadic_refinement() #7
            curve = curve.dyadic_refinement() #11 #
            curve.parent = None
            if not self.makefresh_blends:
                print 'ERROR: we do not have longidutinal that'
                print 'meet the new hull fwd tcurve[0] !'
            #    curve = bhull.lcurvelist[i]
            #--------------------------------------------------
            #
            self.solver_records[i]    = records(curve)
            #    
            #--------------------------------------------------
            # bulb vertices to match with bulb style
            #
            #  INTERPOLATION of bulb boundary!
            #
            vtinterp = []
            for j in range(len(actually_interp)):
                vtinterp.append(actually_interp[j].vertices[i])
            #
            self.solver_records[i].vtinterp = vtinterp
            #    
            #--------------------------------------------------
            # 
            if i > self.len_bulb_boundary-1:
                self.solver_records[i].nosetangent = True
            elif i > self.len_bulb_boundary-1:
                self.solver_records[i].blendtan = True
        
        ##
        ##****************************************************************
        ## SOLVE
        self.complex_fwd_Lnet = []
        self.complex_fwd_lcurves = []
        for i in range(self.tnum):
            print '------------------------------------------'
            print '\n\nBegin Solving New Longitudinal({}):\n\n'.format(i)
            print '------------------------------------------'
            curve = self.solver_records[i].curve
            Lspline = self.setup_fwd_bulb_interpolation_solver2(
                                        thbcurve    = curve,
                                        bareverts   = self.solver_records[i].vtinterp,
                                        nosetan     = self.solver_records[i].nosetangent,
                                        blendtan     = self.solver_records[i].blendtan )
            
            Lspline = uopt.THBsolver(thbcurve=curve,
                                     Lspline=Lspline,
                                     maxlevel=self.maxlevel_thbsolver,
                                     normalize=False)
            # done
            self.complex_fwd_Lnet.append(Lspline)
            self.complex_fwd_lcurves.append(Lspline.curve)
        
        ##
        ##****************************************************************
        ## set the hull curves
        """
        cv = self.bulb.hullbulb.tcurvelist[0]
        ax = cv.plotcurvehierarchy()
        for cv in self.bulb.hullbulb.tcurvelist:
            ax = cv.plotcurvehierarchy(canvas = ax)
        for cv in self.complex_fwd_lcurves:
            ax = cv.plotcurvehierarchy(canvas = ax)
            
        ax = self.barehull.hullfwd.tcurvelist[0].plotcurvehierarchy(canvas = ax)
        
        """
        #level map:
        self.lmap = {4:0,
                     5:1,
                     7:2,
                     11:3,
                     19:4,
                     35:5}
            
#        uln = self.lmap[len(self.complex_fwd_lcurves)]
#        vln = self.lmap[self.complex_fwd_lcurves[0].n]
#        
#        
#        fs = make_surface(lcurvenet = self.complex_fwd_lcurves,
#                          local_tcurvenet = self.hullcurves2[:4],
#                          ul=uln,
#                          vl=vln,
#                          refine=False)
        
        
        test_idea = []
        for curve in self.complex_fwd_lcurves:
            test_idea.append(curve.get_finest_curve())
            
        uln = self.lmap[len(test_idea)]
        vln = self.lmap[test_idea[0].n]
        
        fs = make_surface(lcurvenet = test_idea,
                          local_tcurvenet = self.hullcurves2[:4],
                          ul=uln,
                          vl=vln,
                          refine=True)
        
        fs.compute_surface()
        
        
        self.complexhull = self.barehull
        self.complexhull.newsurf = fs
        
        
        #"""
        #self.complexhull.hullfwd = fs
        
        bhull = self.barehull.hullfwd.get_finest_surface()
        if bhull.ulevel ==  fs.ulevel and bhull.vlevel == fs.vlevel:
            self.complexhull.hullfwd = bhull.parent
            bhull.parent.child = fs
            fs.parent = bhull.parent
            bhull.parent = None
            self.complexhull.hullfwd = fs
        elif bhull.ulevel ==  fs.ulevel+1 and bhull.vlevel == fs.vlevel+1:
            self.complexhull.hullfwd = bhull
            bhull.child = fs
            fs.parent = bhull
            self.complexhull.hullfwd = fs
        else:
            self.complexhull.hullfwd = fs
        
        self.complexhull.hullfwd.compute_surface()
        self.plotcomplex_hull(fancyish=True,
                              simplefancy=True)
        #alt change to fs.set_THB_curves_from_vertex_net(override_u=True, override_v = True)
        
        #fs = self.complexhull.hullfwd.get_finest_surface()
        #fs.bounds[0][0].sup = .75 #.8125 #.75 #.6875 #.46875 #.4375 #.40625#.375#.34375#3125
        #fs.bounds[0][1].sup = .4
        fs.bounds[0][0].sup = 1.
        fs.bounds[0][1].sup = 1.
        fs.compute_surface()
        self.complexhull.hullmid.numpu = 20
        self.complexhull.hullmid.numpv = 5
        self.complexhull.hullmid.compute_surface()
        self.complexhull.hullfwd.numpu = 30
        self.complexhull.hullfwd.numpv = 30
        self.complexhull.hullfwd.compute_surface()
        
        #self.complexhull.hullaft.numpu = 30
        #self.complexhull.hullaft.numpv = 30
        #self.complexhull.hullaft.compute_surface()
        
        fs.tcurvelist[0].bounds.bounds[0].sup = 1.0
        fs.tcurvelist[0].compute_curve()
        self.plotcomplex_hull(fancyish=True,color=True,simplefancy=True)
        self.plotcomplex_hull()
        
        #
        #"""
        ax = fs.tcurvelist[0].plotcurvehierarchy()
        for cv in fs.tcurvelist:
            cv = rbspline(cv.vertices, k=cv.k, nump=cv.nump)
            cv.compute_curve()
            ax = cv.plotcurvehierarchy(canvas = ax)
        #"""
        return
    
    
    def add_bulb(self):
        """setup the vertex system
        and make the longitudinals which define the 
        bbow+lower bare hull portion of the complex hull.
        
        
        ideas for how to do this:
        
        0.) dyadiaclly refine the bare hull
            so that it has the DOFs to incorporate the bulb
            
        ------------ option 1-------------------------------------
        *.1) find projection of bulb aft tcurve
            onto bare hull surface (at the level used to swallow it)
            
        *.2) find parametric locations of bulb boundary
            projection to bare hull surface
            
        *.3) create a (prjective) boundary curve here which may prove useful...
            (for lofting a blending surface)
            
        ------------ option 2-------------------------------------
        *.1)  find enough dofs on the bare hull (e.g. (7x7)) 
            to interpolate the bulb transverse curves
            
        ----------------------------------------------------------
        
        1.) activate the parametric locations on the fine bare hull
            surface so that they encompass the (7x7) bulb DOFs.
        
        
        dev
        ----------
        
        SD.bbow.hullsurf.surfNet[:,0] == SD.bbow.tcurvenet[0].vertices
        
        
        SD.bbow.hullsurf.surfNet[0,:] == SD.bbow.lcurvenet[0].vertices
        
        
                
        fs is self.barehull.hullfwd.child.child.child
        Out[301]: True
        
        fs.level
        Out[302]: 3
        
        fs.ulevel
        Out[303]: 5
        
        fs.vlevel
        Out[304]: 5
        
        """
        #
        k = self.barehull.hullfwd.tcurvelist[0].k
        nump = self.barehull.hullfwd.tcurvelist[0].nump
        #
        numbcs = self.num2use #the fwdmost numbcs (number) of  transverse bulb curves to use for the bulb
        #5 would be all of them
        #
        #**********************************************************************
        # save the base fwd bare hull
        self.original_barehullfwd = copy.deepcopy(self.barehull.hullfwd) #slooooooow
        #
        #**********************************************************************
        # get the finest, fine (enough) level bare hull fwd -> complex hull to be
        #
        fs = self.setup_addbulb()  #dyadic refinement of the bare hull fwd surface
        #
        #
        #**********************************************************************
        # exterior tcurve boundary vertices of the bulb, 
        # usde the fwdmost 3 transverse bulb
        # curves in the initial experiment
        #
        # now just take the aft of the bulb!:
        outline = self.bulb.return_outline_vertices(num2use = numbcs) 
        
        #TODO:
        """
          use outline to design the initial longitudinals.
        That is, make outline the real fwd transverse
        then interpolate it when attaching the bbow.
        """
        
        #
        # add a connective vertex or two or threee back on the tcurve
#        index=0
#        while fs.tcurvelist[0].vertices[index,1] < outline[-1,1]:
#            index+=1
#        outline = list(outline)
#        #
#        #**********************************************************************
#        # transition curves added to the bulb solve:
#        #self.num_extra_boundary=3 moved to init
#        #
#        #**********************************************************************
#        #
#        for i in range(self.num_extra_boundary):
#            outline.append(fs.tcurvelist[0].vertices[index+i])
#        outline = np.asarray(outline)
        
        #use all of the curves now.
        self.tnum = self.barehull.hullfwd.tcurvelist[0].n #
        #self.tnum = len(outline)
        #
        #**********************************************************************
        # boundary curves and vertices for clarity:
        #(non-essential)
        self.fwd_boundary_vertices      = outline
        self.fwd_boundary_curve         = rbspline(self.fwd_boundary_vertices,
                                                   k=k,
                                                   nump=nump)
        self.aft_boundary_vertices      = fs.tcurvelist[-1].vertices
        self.aft_boundary_curve         = rbspline(self.aft_boundary_vertices,
                                                   k=k,
                                                   nump=nump)
        self.aft_boundary_coarse_curve  = self.barehull.tcurvenet_fwd[-1]
        #
        #**********************************************************************
        # setup the complex hull's fwd most outline transverse 
        #  and
        # initialize the complex lcurves with 
        # the right boundaries (i.e. the new outline!)
        self.set_fwdboundary_complexhull(outline,fs)
        #
        #**********************************************************************
        # doit
        self.makeblendedsurface(outline, fs)
        #
        return
    
    def makeblendedsurface(self, outline, fs):
        """make the bulb interpolating longitudinals
        of the complex hull.
        
        
        parameters
        ---------
        tnum : the nubmer of longitudinal curves that
        will go into the making of this bulb
        (the bulb will then live as a part of fs's level, by the way
        (complex hullw/bulb to be on fs's level, as in fact, fs itself.))
        
        todo
        ---------
        * add an ignore 1 bare hull tcurve aft of the bow
        -because that was designed without real intention of
        being a bbow one day!
        -alt mod hull maker to stay flat with cpkeel all the way
        to the final rise.  Actually this would be better anyway.
        """
        print 'setup solver constraints for complex bulb integration Lspline solving'
        # 
        #maybe the key number here:
        tnum = self.tnum #number of lsplines making the bulb a part of the hull on this level
        #
        #*******************************
        # how many interior curves of the fwd bare hull to use?
        use_both_bare_hull_curves = True
        self.use_both_bare_hull_curves = use_both_bare_hull_curves
        #
        if use_both_bare_hull_curves:
            numberof_bare_hull_tcurves_to_interpolate = \
                            len(self.barehull.tcurvenet_fwd[1:])
        else:
            numberof_bare_hull_tcurves_to_interpolate = \
                            len(self.barehull.tcurvenet_fwd[2:])
        #
        #*******************************
        #
        self.nmberbhtcti = numberof_bare_hull_tcurves_to_interpolate
        #
        num2use = self.num2use #number of bulb tcurves to use
        """
        #stuff in here never gets used
        bhtnum = len(self.barehull.hullfwd.tcurvelist)#bare hullfwd: num of tcurves
        #
        len_one_bulb_tcurve = self.bulb.tcurvenet[0].n
        #----------------
        k = fs.tcurvelist[0].k
        nump = fs.tcurvelist[0].nump
        #----------------
        nt = fs.tcurvelist[0].n
        #"""
        lvl = fs.tcurvelist[0].level
        #
        #****************************************
        #  * initialize list of tcurve from bare hull -> to be interpolated
        # why do you want this one?
        # Answer: because you have not yet fixed it to the right place?
        """
        c1 = rbspline(fs.tcurvelist[0].vertices, k,nump)
        c1.level = lvl
        c1.bounds = fs.tcurvelist[0].bounds
        #"""
        #
        #****************************************
        #  * chief tcurves from the bare hull => interpolate these!
        #FPcurvenet = [c1] #see above -why ?double count?
        FPcurvenet = [] 
        #TODO: need inital solver to be more aware
        # of what we are doing in order to use this curve!:
        if self.use_both_bare_hull_curves:
            for el in self.barehull.tcurvenet_fwd[1:-1]:  #fixing bug here: :-1 so as to leave off the last!
                el.level = 2 #level =2 because its a 7 vertex curve.
                lel = el.level
                while lel<lvl:
                    el = el.dyadic_refinement()
                    lel = el.level
                FPcurvenet.append(el)
        else:
            for el in self.barehull.tcurvenet_fwd[2:-1]: #fixing bug here: :-1 so as to leave off the last!
                el.level = 2 #level =2 because its a 7 vertex curve.
                lel = el.level
                while lel<lvl:
                    el = el.dyadic_refinement()
                    lel = el.level
                FPcurvenet.append(el)
        #tvertnet = uopt.make_vertex_matrix(FPcurvenet) #bare hull vertices
        self.FPcurvenet = FPcurvenet
            
            
        #
        #****************************************
        # Begin solver Lagrange setup
        #only solving the 11 curves on the bottom
        
        
        #(1) here are the traditional hull curves fwd, 
        #    updated to have 35 vertices:
        #for i in range(tnum):  #tnum=len(outline)
        #    curve = fs.lcurvelist[i]
        """ just use them directly, below """
        
        #(2) next we need to properly attach the 
        #    bulb transverses to longi solver curves
        #
        # working in the natural way a human would get these
        # since they are running aft to front,
        # due to construction, etc...
        #
        blist = range(num2use)
        blist = [-el-1 for el in blist[1:]]
        iblist = range(num2use)
        iblist = [-el-1 for el in iblist[1:]]
        iblist.reverse()
        # you do not need -1 as that has been taken care of
        # [-2,-3,...,-num2use-1-1] so we sweep in 
        # all bulb tcurves that we want to use
        
        #************************************
        # make a map from this reverse indexed stuff back to fwd indices
        # starting from the bow.  This will be how we find our fixity 
        # indices later
        invertmap = range(num2use)
        invertmap.pop(0)
        ivm = {}
        for m,im in zip(blist,invertmap):
            ivm[m] = im
        self.ivm = ivm
        #************************************
        
        #
        # do not grab the exterior curve; i.e. bulb tcurvelist[-1]
        # in the initial conditions 
        # i.e. [-3,-2,-1] for num2use 
        # would double count the fwd end of the bulb tcurve:: the boundary.
        
        mapverts = {}
        for i in range(tnum):
            if i < num2use:
                ai = i
                mapverts[i] = iblist[:ai]
                mapverts[i].reverse()
            elif i >= tnum - self.num_extra_boundary - num2use:
                if i < tnum - self.num_extra_boundary:
                    print 'i = ',i
                    ai = tnum-self.num_extra_boundary-i-1
                    mapverts[i] = iblist[:ai]
                    mapverts[i].reverse()
                else:
                    mapverts[i] = []
            else:
                mapverts[i] = blist
        
        
        
        
        
        
        # loop longitudinals to solve
        self.solver_records = {}
        for i in range(tnum): #each i is to be a solve
            #
            curve = fs.lcurvelist[i] #solve this curve
            
            self.solver_records[i] = records(curve)
            self.solver_records[i].bulbmap = mapverts[i]
            #
            #****************************************
            vtinterp = [] #transverse vertices which the longitudinal 
            # curves must interpolate
            #****************************************
            indices = [] #for checking integers - not really used 
            #****************************************
            FPLongFix = [] #longitudinal vertices
            #  which denote the present longitudinals 
            # -which already interpolate the tcurve vertices
            # but are not global to the complex hull
            FPLongboundary = [] #maybe we just want to fix the boundary
            #****************************************
            
            # tangent condition for the vertices just above the bulb:
            flatnose = False #(vertices at interface to rest of surf need to match tangency - would be nice to "blend tangency")
            if i >= tnum - self.num_extra_boundary:
                flatnose = True
            #loop bulb-curves to interpolate
            #*****************************************************
            #case structure may need to consume this:
            for el in mapverts[i]: #[::-1]:
                indices.append(el)
                """
                implementing the degeneracy conditions
                results in some complex logic to get the
                vertices (to be interpolated) in the right place
                """
                #*****************************************************
                # case structure needed again
                if i < num2use: #[0,1,2] #(for a 3 t-segment bulb)
                    #---------------------------------------------------------
                    # we are at the bottom boundary of the bulb:
                    # bottom lsplines catch the 0th vertex of the elth bulb tcurve
                    # these lsplines are degenerate
                    '''lets try fixing the bulb'''
                    vt = None#self.bulb.tcurvenet[el].vertices[0] # num2use, el , 1  #do nothing!, we are going to fix it!
                    #
                    Lv = self.bulb.lcurvenet[0].vertices[el]
                    iL = el
                    #Lb = self.bulb.tcurvenet[el].vertices[0] #these are boundary points, emphasize them!
                    Lb = self.bulb.lcurvenet[0].vertices[el]
                    #---------------------------------------------------------
                elif i >= tnum - self.num_extra_boundary - num2use: #[8,9,10]
                    #---------------------------------------------------------
                    if i < tnum - self.num_extra_boundary:
                        # we are at the top boundary of the bulb:
                        #top lsplines catch the '6th' vertex of the el'th bulb tcurve
                        #again degenerate
                        '''lets try fixing the bulb'''
                        vt = None#self.bulb.tcurvenet[el].vertices[self.bulb.tcurvenet[0].n-1] #do nothing!, we are going to fix it!
                        #
                        Lv = self.bulb.lcurvenet[self.bulb.tcurvenet[0].n-1].vertices[el]
                        iL = el
                        #Lb = self.bulb.tcurvenet[el].vertices[self.bulb.tcurvenet[0].n-1]  #these are boundary points, emphasize them!
                        Lb = self.bulb.lcurvenet[self.bulb.tcurvenet[0].n-1].vertices[el]
                        #---------------------------------------------------------
                    else:
                        #we are 'above' the bulb
                        #regular interior bulb interpolation vertices:
                        '''lets try fixing the bulb'''
                        vt = None
                        Lv = []
                        iL = el
                        Lb = [] #no special (degenerate) boundary
                        flatnose = True
                        assert(False),'we never get here'
                    #---------------------------------------------------------
                else:
                    # we are on the bulb interior
                    #---------------------------------------------------------
                    '''lets try fixing the bulb'''
                    vt = self.bulb.tcurvenet[el].vertices[i-num2use+1]  #do nothing!, we are going to fix it!
                    Lv = self.bulb.lcurvenet[i-num2use+1].vertices[el]
                    iL = el
                    Lb = [] #no special (degenerate) boundary
                    #---------------------------------------------------------
                # end cases
                #*****************************************************
                if vt is not None:
                    vtinterp.append(vt) #vertices the longi's must interpolate
                FPLongFix.append({iL:Lv}) #vertices that the longis can use to do so.
                FPLongboundary.append({iL:Lb})
            #*****************************************************
            #loop bare hull curves to interpolate
            """ #use barevts instead!  --to differentiate these
            for j in range(len(FPcurvenet)):#number of what?
                #
                vtinterp.append(FPcurvenet[j].vertices[i])
                self.solver_records[i].barehullfwd_tcurveindices.append(j)
            #"""
            #
            #vtinterp = np.asarray(vtinterp)
            self.solver_records[i].vtinterp = None#vtinterp #vertices to interpolate
            self.solver_records[i].lfix = None #FPLongFix #fix interior bulb lcurve vertices
            self.solver_records[i].lboundary = None #FPLongboundary #fix bulb boundary
            if flatnose:
                self.solver_records[i].nosetangent = True
        
        ##
        ##**********************************************
        ## Generate the curves via solver
        self.complex_fwd_Lnet = []
        self.complex_fwd_lcurves = []
        for i in range(tnum): #each i is to be a solve
            vini        = self.fwd_boundary_curve.vertices[i]
            vfini       = self.aft_boundary_curve.vertices[i]
            #--------------------------------------------------
            vtinterp    = self.solver_records[i].vtinterp
            lfix        = self.solver_records[i].lfix
            lbounds     = self.solver_records[i].lboundary
            nosetangent = self.solver_records[i].nosetangent
            #--------------------------------------------------
            bvts = []
            for curve in self.FPcurvenet:
                bvts.append(curve.vertices[i])
            bvts = np.asarray(bvts)
            #--------------------------------------------------
            vl          = [] #list(self.solver_records[i].vtinterp)
            #
            vl.append(vfini)
            vl.insert(0,vini)
            #
            # initialize curve
            vlarray = np.asarray(vl)
            if self.makefresh_blends:
                linverts = linear_vertices(vl[0],vl[-1],4)
                curve = rbspline(linverts,k=4,nump=30)
                curve = curve.dyadic_refinement() #5
                curve = curve.dyadic_refinement() #7
                curve = curve.dyadic_refinement() #11 #so there is enough to catch all the degenerate bounds
            else:
                curve = self.solver_records[i].curve #curve might have level >3.
            curve.parent = None
            # initialize solver and Lagrangian
            Lspline = self.setup_fwd_bulb_interpolation_solver(thbcurve = curve,
                                                               vtinterp = vtinterp,
                                                               lfix = lfix,
                                                               lbounds = lbounds,
                                                               bareverts = bvts,
                                                               nosetan = nosetangent,
                                                               this_kind = self.bulbblend,
                                                               fixedkind = self.barehullcurveblend)
            # use THB solver process
            Lspline = uopt.THBsolver(thbcurve=curve,
                                     vtinterp=vlarray,
                                     Lspline=Lspline,
                                     maxlevel=self.maxlevel_thbsolver,
                                     normalize=False)
            # done
            self.complex_fwd_Lnet.append(Lspline)
            self.complex_fwd_lcurves.append(Lspline.curve)
            
        ##
        ##****************************************************************
        ## set the hull curves
        fine_curves = []
        #----------------------------------------
        if self.make_hires: #ad hoc: hi res is 35 vertex, i.e. dyadic level 5 [0:4,1:5,2:7,3:11,4:19,5:35]
            for curve in self.complex_fwd_lcurves:
                this = curve.get_finest_curve()
                while this.n < self.highres_num:
                    this = this.dyadic_refinement()
                fine_curves.append(this)
        else: #ad hoc: low res is "19" i.e. dyadic level 4
            for curve in self.complex_fwd_lcurves:
                this = curve.get_finest_curve()
                while this.n < self.lowres_num:
                    this = this.dyadic_refinement()
                fine_curves.append(this)
        #-----------------------------------------
        
        self.complexhull = self.barehull
        fs = self.complexhull.hullfwd.get_finest_surface()
        #
        for i in range(len(self.complex_fwd_lcurves)):
            while fine_curves[i].level < fs.lcurvelist[i].level:
                fine_curves[i] = fine_curves[i].dyadic_refinement()
        #
        for i in range(len(self.complex_fwd_lcurves)):
            fs.lcurvelist[i].vertices = fine_curves[i].vertices
            fs.lcurvelist[i].compute_curve()
            fs.vertices[i,:] = fine_curves[i].vertices
        
        for i in range(fs.nv):
            fs.tcurvelist[i].vertices[:] = fs.vertices[:,i]
            fs.tcurvelist[i].compute_curve()
        
        fs.compute_surface()
        self.complexhull.hullfwd.compute_surface()
        self.complexhull.plot(fancyish=True,color=False)
        self.plotcomplex_hull(fancyish=True,
                              simplefancy=True)
        #alt change to fs.set_THB_curves_from_vertex_net(override_u=True, override_v = True)
        
        #fs = self.complexhull.hullfwd.get_finest_surface()
        fs.bounds[0][0].sup = 1. #.8125 #.75 #.6875 #.46875 #.4375 #.40625#.375#.34375#3125
        fs.bounds[0][1].sup = 1.
        fs.compute_surface()
        self.complexhull.hullmid.numpu = 20
        self.complexhull.hullmid.numpv = 5
        self.complexhull.hullmid.compute_surface()
        self.complexhull.hullfwd.numpu = 30
        self.complexhull.hullfwd.numpv = 30
        self.complexhull.hullfwd.compute_surface()
        
        #self.complexhull.hullaft.numpu = 30
        #self.complexhull.hullaft.numpv = 30
        #self.complexhull.hullaft.compute_surface()
        
        #
        fs.tcurvelist[0].bounds.bounds[0].sup = 1.0
        fs.tcurvelist[0].compute_curve()
        self.plotcomplex_hull(fancyish=True,color=True,simplefancy=True)
        self.plotcomplex_hull()
        
        #
        #"""
        ax = fs.tcurvelist[0].plotcurvehierarchy()
        for cv in fs.tcurvelist:
            cv = rbspline(cv.vertices, k=cv.k, nump=cv.nump)
            cv.compute_curve()
            ax = cv.plotcurvehierarchy(canvas = ax)
        #"""
        return
    
    def makeblendedsurface_NEWER(self, outline, fs):
        """make the bulb interpolating longitudinals
        of the complex hull.
        
        
        parameters
        ---------
        tnum : the nubmer of longitudinal curves that
        will go into the making of this bulb
        (the bulb will then live as a part of fs's level, by the way
        (complex hullw/bulb to be on fs's level, as in fact, fs itself.))
        
        todo
        ---------
        * add an ignore 1 bare hull tcurve aft of the bow
        -because that was designed without real intention of
        being a bbow one day!
        -alt mod hull maker to stay flat with cpkeel all the way
        to the final rise.  Actually this would be better anyway.
        """
        print 'setup solver constraints for complex bulb integration Lspline solving'
        # 
        #maybe the key number here:
        tnum = self.tnum #number of lsplines making the bulb a part of the hull on this level
        #
        #*******************************
        # how many interior curves of the fwd bare hull to use?
        use_both_bare_hull_curves = True
        self.use_both_bare_hull_curves = use_both_bare_hull_curves
        #
        if use_both_bare_hull_curves:
            numberof_bare_hull_tcurves_to_interpolate = \
                            len(self.barehull.tcurvenet_fwd[1:])
        else:
            numberof_bare_hull_tcurves_to_interpolate = \
                            len(self.barehull.tcurvenet_fwd[2:])
        #
        #*******************************
        #
        self.nmberbhtcti = numberof_bare_hull_tcurves_to_interpolate
        #
        num2use = self.num2use #number of bulb tcurves to use
        """
        #stuff in here never gets used
        bhtnum = len(self.barehull.hullfwd.tcurvelist)#bare hullfwd: num of tcurves
        #
        len_one_bulb_tcurve = self.bulb.tcurvenet[0].n
        #----------------
        k = fs.tcurvelist[0].k
        nump = fs.tcurvelist[0].nump
        #----------------
        nt = fs.tcurvelist[0].n
        #"""
        lvl = fs.tcurvelist[0].level
        #
        #****************************************
        #  * initialize list of tcurve from bare hull -> to be interpolated
        # why do you want this one?
        # Answer: because you have not yet fixed it to the right place?
        """
        c1 = rbspline(fs.tcurvelist[0].vertices, k,nump)
        c1.level = lvl
        c1.bounds = fs.tcurvelist[0].bounds
        #"""
        #
        #****************************************
        #  * chief tcurves from the bare hull => interpolate these!
        #FPcurvenet = [c1] #see above -why ?double count?
        FPcurvenet = [] 
        #TODO: need inital solver to be more aware
        # of what we are doing in order to use this curve!:
        if self.use_both_bare_hull_curves:
            for el in self.barehull.tcurvenet_fwd[1:-1]:  #fixing bug here: :-1 so as to leave off the last!
                el.level = 2 #level =2 because its a 7 vertex curve.
                lel = el.level
                while lel<lvl:
                    el = el.dyadic_refinement()
                    lel = el.level
                FPcurvenet.append(el)
        else:
            for el in self.barehull.tcurvenet_fwd[2:-1]: #fixing bug here: :-1 so as to leave off the last!
                el.level = 2 #level =2 because its a 7 vertex curve.
                lel = el.level
                while lel<lvl:
                    el = el.dyadic_refinement()
                    lel = el.level
                FPcurvenet.append(el)
        #tvertnet = uopt.make_vertex_matrix(FPcurvenet) #bare hull vertices
        self.FPcurvenet = FPcurvenet
            
            
        #
        #****************************************
        # Begin solver Lagrange setup
        #only solving the 11 curves on the bottom
        
        
        #(1) here are the traditional hull curves fwd, 
        #    updated to have 35 vertices:
        #for i in range(tnum):  #tnum=len(outline)
        #    curve = fs.lcurvelist[i]
        """ just use them directly, below """
        
        #(2) next we need to properly attach the 
        #    bulb transverses to longi solver curves
        #
        # working in the natural way a human would get these
        # since they are running aft to front,
        # due to construction, etc...
        #
        blist = range(num2use)
        blist = [-el-1 for el in blist[1:]]
        iblist = range(num2use)
        iblist = [-el-1 for el in iblist[1:]]
        iblist.reverse()
        # you do not need -1 as that has been taken care of
        # [-2,-3,...,-num2use-1-1] so we sweep in 
        # all bulb tcurves that we want to use
        
        #************************************
        # make a map from this reverse indexed stuff back to fwd indices
        # starting from the bow.  This will be how we find our fixity 
        # indices later
        invertmap = range(num2use)
        invertmap.pop(0)
        ivm = {}
        for m,im in zip(blist,invertmap):
            ivm[m] = im
        self.ivm = ivm
        #************************************
        
        #
        # do not grab the exterior curve; i.e. bulb tcurvelist[-1]
        # in the initial conditions 
        # i.e. [-3,-2,-1] for num2use 
        # would double count the fwd end of the bulb tcurve:: the boundary.
        
        mapverts = {}
        for i in range(tnum):
            if i < num2use:
                ai = i
                mapverts[i] = iblist[:ai]
                mapverts[i].reverse()
            elif i >= tnum - self.num_extra_boundary - num2use:
                if i < tnum - self.num_extra_boundary:
                    print 'i = ',i
                    ai = tnum-self.num_extra_boundary-i-1
                    mapverts[i] = iblist[:ai]
                    mapverts[i].reverse()
                else:
                    mapverts[i] = []
            else:
                mapverts[i] = blist
        
        
        
        
        
        
        # loop longitudinals to solve
        self.solver_records = {}
        for i in range(tnum): #each i is to be a solve
            #
            curve = fs.lcurvelist[i] #solve this curve
            
            self.solver_records[i] = records(curve)
            self.solver_records[i].bulbmap = mapverts[i]
            #
            #****************************************
            vtinterp = [] #transverse vertices which the longitudinal 
            # curves must interpolate
            #****************************************
            indices = [] #for checking integers - not really used 
            #****************************************
            FPLongFix = [] #longitudinal vertices
            #  which denote the present longitudinals 
            # -which already interpolate the tcurve vertices
            # but are not global to the complex hull
            FPLongboundary = [] #maybe we just want to fix the boundary
            #****************************************
            
            # tangent condition for the vertices just above the bulb:
            flatnose = False #(vertices at interface to rest of surf need to match tangency - would be nice to "blend tangency")
            if i >= tnum - self.num_extra_boundary:
                flatnose = True
            #loop bulb-curves to interpolate
            #*****************************************************
            #case structure may need to consume this:
            # DELETED STUFF
            #*****************************************************
            #loop bare hull curves to interpolate
            """ #use barevts instead!  --to differentiate these
            for j in range(len(FPcurvenet)):#number of what?
                #
                vtinterp.append(FPcurvenet[j].vertices[i])
                self.solver_records[i].barehullfwd_tcurveindices.append(j)
            #"""
            #
            #vtinterp = np.asarray(vtinterp)
            self.solver_records[i].vtinterp = None#vtinterp #vertices to interpolate
            self.solver_records[i].lfix = None #FPLongFix #fix interior bulb lcurve vertices
            self.solver_records[i].lboundary = None #FPLongboundary #fix bulb boundary
            if flatnose:
                self.solver_records[i].nosetangent = True
        
        ##
        ##**********************************************
        ## Generate the curves via solver
        self.complex_fwd_Lnet = []
        self.complex_fwd_lcurves = []
        for i in range(tnum): #each i is to be a solve
            vini        = self.fwd_boundary_curve.vertices[i]
            vfini       = self.aft_boundary_curve.vertices[i]
            #--------------------------------------------------
            vtinterp    = self.solver_records[i].vtinterp
            lfix        = self.solver_records[i].lfix
            lbounds     = self.solver_records[i].lboundary
            nosetangent = self.solver_records[i].nosetangent
            #--------------------------------------------------
            bvts = []
            for curve in self.FPcurvenet:
                bvts.append(curve.vertices[i])
            bvts = np.asarray(bvts)
            #--------------------------------------------------
            vl          = [] #list(self.solver_records[i].vtinterp)
            #
            vl.append(vfini)
            vl.insert(0,vini)
            #
            # initialize curve
            vlarray = np.asarray(vl)
            if self.makefresh_blends:
                linverts = linear_vertices(vl[0],vl[-1],4)
                curve = rbspline(linverts,k=4,nump=30)
                curve = curve.dyadic_refinement() #5
                curve = curve.dyadic_refinement() #7
                curve = curve.dyadic_refinement() #11 #so there is enough to catch all the degenerate bounds
            else:
                curve = self.solver_records[i].curve #curve might have level >3.
            curve.parent = None
            # initialize solver and Lagrangian
            Lspline = self.setup_fwd_bulb_interpolation_solver(thbcurve = curve,
                                                               vtinterp = vtinterp,
                                                               lfix = lfix,
                                                               lbounds = lbounds,
                                                               bareverts = bvts,
                                                               nosetan = nosetangent,
                                                               this_kind = self.bulbblend,
                                                               fixedkind = self.barehullcurveblend)
            # use THB solver process
            Lspline = uopt.THBsolver(thbcurve=curve,
                                     vtinterp=vlarray,
                                     Lspline=Lspline,
                                     maxlevel=self.maxlevel_thbsolver,
                                     normalize=False)
            # done
            self.complex_fwd_Lnet.append(Lspline)
            self.complex_fwd_lcurves.append(Lspline.curve)
            
        ##
        ##****************************************************************
        ## set the hull curves
        fine_curves = []
        #----------------------------------------
        if self.make_hires: #ad hoc: hi res is 35 vertex, i.e. dyadic level 5 [0:4,1:5,2:7,3:11,4:19,5:35]
            for curve in self.complex_fwd_lcurves:
                this = curve.get_finest_curve()
                while this.n < self.highres_num:
                    this = this.dyadic_refinement()
                fine_curves.append(this)
        else: #ad hoc: low res is "19" i.e. dyadic level 4
            for curve in self.complex_fwd_lcurves:
                this = curve.get_finest_curve()
                while this.n < self.lowres_num:
                    this = this.dyadic_refinement()
                fine_curves.append(this)
        #-----------------------------------------
        
        self.complexhull = self.barehull
        fs = self.complexhull.hullfwd.get_finest_surface()
        #
        for i in range(len(self.complex_fwd_lcurves)):
            while fine_curves[i].level < fs.lcurvelist[i].level:
                fine_curves[i] = fine_curves[i].dyadic_refinement()
        #
        for i in range(len(self.complex_fwd_lcurves)):
            fs.lcurvelist[i].vertices = fine_curves[i].vertices
            fs.lcurvelist[i].compute_curve()
            fs.vertices[i,:] = fine_curves[i].vertices
        
        for i in range(fs.nv):
            fs.tcurvelist[i].vertices[:] = fs.vertices[:,i]
            fs.tcurvelist[i].compute_curve()
        
        fs.compute_surface()
        self.complexhull.hullfwd.compute_surface()
        self.complexhull.plot(fancyish=True,color=False)
        self.plotcomplex_hull(fancyish=True,
                              simplefancy=True)
        #alt change to fs.set_THB_curves_from_vertex_net(override_u=True, override_v = True)
        
        #fs = self.complexhull.hullfwd.get_finest_surface()
        fs.bounds[0][0].sup = .375 #.8125 #.75 #.6875 #.46875 #.4375 #.40625#.375#.34375#3125
        fs.bounds[0][1].sup = .1875
        fs.compute_surface()
        self.complexhull.hullmid.numpu = 20
        self.complexhull.hullmid.numpv = 5
        self.complexhull.hullmid.compute_surface()
        self.complexhull.hullfwd.numpu = 30
        self.complexhull.hullfwd.numpv = 30
        self.complexhull.hullfwd.compute_surface()
        
        #self.complexhull.hullaft.numpu = 30
        #self.complexhull.hullaft.numpv = 30
        #self.complexhull.hullaft.compute_surface()
        
        self.plotcomplex_hull(fancyish=True,color=True,simplefancy=True)
        #
        return
    def check_solver(self, fs=None):
        if fs is None:
            fs = self.complexhull.hullfwd.get_finest_surface()
            
        if self.hullcurves2 is None:
            self.hullcurves2 = get_tcurves()
            
        
        # solved curves:
        ax = self.aft_boundary_curve.plotcurvehierarchy()
        ax = self.fwd_boundary_curve.plotcurvehierarchy(canvas=ax)
        for i in range(self.tnum):
            #fs.lcurvelist[i].compute_curve()
            curve = rbspline(fs.lcurvelist[i].vertices,k=4,nump=30)
            #ax = fs.lcurvelist[i].plotcurvehierarchy(canvas = ax)
            ax = curve.plotcurvehierarchy(canvas = ax)
            
        for tcurve in self.bulb.tcurvenet:
            ax = tcurve.plotcurvehierarchy(canvas=ax)
        
        
        ax = self.aft_boundary_curve.plotcurvehierarchy()
        ax = self.fwd_boundary_curve.plotcurvehierarchy(canvas=ax)
        for lcurve in fs.lcurvelist:
            #lcurve.compute_curve()
            curve = rbspline(lcurve.vertices,k=4,nump=30)
            ax = curve.plotcurvehierarchy(canvas=ax)
        for tcurve in self.hullcurves2:
            ax = tcurve.plotcurvehierarchy(canvas=ax)
        
        
        for tcurve in self.bulb.tcurvenet:
            ax = tcurve.plotcurvehierarchy(canvas=ax)
        
        return
    
    def setup_fwd_bulb_interpolation_solver(self, 
                                           thbcurve, 
                                           vtinterp, #bulb vertices
                                           lfix, #no longer used?
                                           lbounds, #boundary bulb outline
                                           bareverts, #bare hull vertices
                                           nosetan,
                                           this_kind = 'equality', #bulb vertices
                                           fixedkind = 'equality'): #bare hull vertices
        # 'equality' 'LS'
        """
        """
        print '   complexhull.setup_fwd_bulb_interpolation_solver'
        #*********************************************
        interval_data, small = interval_bounds(thbcurve)
        FPD = FormParameterDict(thbcurve)
        
        LSkind = 'LS'
        
        #big = 1.
        big  = 100000.
        mig =  50000.
        #mig =  10000.
        smallbig = 1000.
        small = .1
        reg = 1.
        scalar = 1.
        #*********************************************
        FPD.add_E1(kind='LS', weight = 1.)  #reg) #
        FPD.add_E2(kind='LS', weight = .5  )  # reg) #
        FPD.add_E3(kind='LS', weight = .05  )  #  reg) #
        #FPD.add_ArcLengthApprox(kind='LS',weight = big)
        FPD.add_ArcLengthApprox(kind='LS',weight = mig)
        #FPD.add_ArcLengthApprox(kind='LS',weight = smallbig)
        #*********************************************
        #ntinterp = len(vtinterp)
        nbvts = len(bareverts)
        #ukbar = np.linspace(0.,1.,ntinterp+nbvts+2)[1:-1]
        ntot = thbcurve.n
        ukbar = np.linspace(0.,1.,ntot)[-1-nbvts:-1]
        #nth = 1
        #totn = ntinterp
        xb = thbcurve.vertices[0,0]
        yb = thbcurve.vertices[0,1]
        zb = thbcurve.vertices[0,2]
        #*********************************************
        # C2 continuity to bulb.
        FPD.add_xFixity(index = 1,
                        value = xb,
                        track_prolongation=False)
        FPD.add_xFixity(index = 2,
                        value = xb,
                        track_prolongation=False)
        FPD.add_yFixity(index = 1,
                        value = yb,
                        track_prolongation=False)
        FPD.add_yFixity(index = 2,
                        value = yb,
                        track_prolongation=False)
        
        """
        try:
            for pt,loc in zip(vtinterp,ukbar[:ntinterp]):
                #*********************************************
                FPD.add_xPointConstraint(kind= this_kind, 
                                         value=pt[0], 
                                         location=loc, 
                                         weight = big ) 
                FPD.add_yPointConstraint(kind= this_kind, 
                                         value=pt[1], 
                                         location=loc, 
                                         weight = big  )
                FPD.add_zPointConstraint(kind= this_kind, 
                                         value=pt[2], 
                                         location=loc, 
                                         weight = big  )
        except:
            print 'no bulb vertices to interpolate'
        #
        #"""
        #for pt,loc in zip(bareverts,ukbar[ntinterp:]):
        for pt,loc in zip(bareverts,ukbar):
            #*********************************************
            FPD.add_xPointConstraint(kind= fixedkind, 
                                     value=pt[0], 
                                     location=loc, 
                                     weight = big ) 
            FPD.add_yPointConstraint(kind= fixedkind, 
                                     value=pt[1], 
                                     location=loc, 
                                     weight = big  )
            FPD.add_zPointConstraint(kind= fixedkind, 
                                     value=pt[2], 
                                     location=loc, 
                                     weight = big  )
        
            #*********************************************
        #*********************************************
        #*********************************************
        """
        for key in self.solver_records:
            lbounds = self.solver_records[key].lboundary
            for record in lbounds:
                for el in record:
                    rc = record[el]
                    #if rc:
                    index = self.ivm[el]
                    vtx = rc
                    try:
                        print 'index = ', el,self.ivm[el]
                        print 'vtx = ',record[el]
                        print key, vtx[0],vtx[1],vtx[2]
                    except:
                        print 'skip key = ',key
            
        """
        # degenerate boundary fixities:
        """
        for record in lbounds:
            for el in record:
                rc = record[el]
                index = self.ivm[el]
                vtx = rc
                try:
                    #*********************************************
                    # add bulb fixing...
                    FPD.add_xFixity(index = index,
                                    value = vtx[0],
                                    track_prolongation=False)
                    FPD.add_yFixity(index = index,
                                    value = vtx[1],
                                    track_prolongation=False)
                    FPD.add_zFixity(index = index,
                                    value = vtx[2],
                                    track_prolongation=False)
                except:
                    print 'no boundary for ',record
        #"""
        ## fix the entire bulb?
        """
        for record in lfix:
            for el in record:
                index = self.ivm[el]
                vtx = record[el]
                #*********************************************
                # add bulb fixing...
                FPD.add_xFixity(index = index,
                                value = vtx[0],
                                track_prolongation=False)
                FPD.add_yFixity(index = index,
                                value = vtx[1],
                                track_prolongation=False)
                FPD.add_zFixity(index = index,
                                value = vtx[2],
                                track_prolongation=False)
        #"""
                #*********************************************
        
        if nosetan:
                FPD.add_yFixity(index = 1,
                                value = yb,
                                track_prolongation=False)
                FPD.add_zFixity(index = 1,
                                value = zb,
                                track_prolongation=False)
                FPD.add_yFixity(index = 2,
                                value = yb,
                                track_prolongation=False)
                FPD.add_zFixity(index = 2,
                                value = zb,
                                track_prolongation=False)
            
        #*********************************************
        #*********************************************
            
            
            
        #*********************************************
        #*********************************************
        #*********************************************
#        # add fore and aft fixiy and or derivative constraints
        #(not really needed because we are not looking at these curves this
        # far away from the bulbous bow)
        ## TODO: bring this back:
        """
        FPD.add_xFixity(index = curve.n-2,
                        value = curve.vertices[-1,0],
                        track_prolongation=False)
        FPD.add_yFixity(index = curve.n-2,
                        value = curve.vertices[-1,1],
                        track_prolongation=False)
#        #------------------------------------------
#        FPD.add_xFixity(index = curve.n-3,
#                        value = curve.vertices[-1,0],
#                        track_prolongation=False)
#        FPD.add_yFixity(index = curve.n-3,
#                        value = curve.vertices[-1,1],
#                        track_prolongation=False)
        #"""
        #*********************************************
        #*********************************************
        #*********************************************
#        FPD.add_AngleConstraint(
#                                      kind      = this_kind,
#                                      location  = 0.,
#                                      value     = 0.,
#                                      weight    = scalar)
#        FPD.add_CurvatureConstraint(
#                                      kind  = this_kind,
#                                      location    = 0.,
#                                      value       = 0.,
#                                      weight      = scalar)
        #*********************************************
        # add bulb fixing... old guess way
#        FPD.add_xFixity(index = 1,
#                        value = curve.vertices[-1,0],
#                        track_prolongation=False)
#        FPD.add_yFixity(index = 1,
#                        value = curve.vertices[-1,1],
#                        track_prolongation=False)
#        #------------------------------------------
#        FPD.add_xFixity(index = 2,
#                        value = curve.vertices[-1,0],
#                        track_prolongation=False)
#        FPD.add_yFixity(index = 2,
#                        value = curve.vertices[-1,1],
#                        track_prolongation=False)
        #*********************************************
        #
        #
        #*********************************************
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                            interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(thbcurve, L, 
                                         data = interval_data)
        #
        #*********************************************
        return Lspline
    
    
    
    
    def setup_fwd_bulb_interpolation_solver2(self, 
                                           thbcurve, 
                                           bareverts, #interior bare hull vertices
                                           nosetan = False,
                                           blendtan = False,
                                           this_kind = 'equality'): #bare hull vertices
        # 'equality' 'LS'
        """
        
        dev
        ----------
        
        thbcurve = self.solver_records[i].curve
        
        
        for i in range(self.tnum):
            print self.solver_records[i].nosetangent
            
            
        for i in range(self.tnum):
            print self.solver_records[i].blendtan
            
            
        #bareverts == vtinterp
        for i in range(self.tnum):
            print self.solver_records[i].vtinterp
            
            
        nbvts = len(bareverts)
        ntot = thbcurve.n
        ukbar = np.linspace(0.,1.,ntot)[-1-nbvts-1:-1-1]
        """
        print '   complexhull.setup_fwd_bulb_interpolation_solver'
        #*********************************************
        interval_data, small = interval_bounds(thbcurve)
        FPD = FormParameterDict(thbcurve)
        
        
        #big = 1.
        big  = 100000.
        mig =  50000.
        #mig =  10000.
        smallbig = 1000.
        small = .1
        reg = 1.
        scalar = 1.
        #*********************************************
        FPD.add_E1(kind='LS', weight = 1.)  #reg) #
        FPD.add_E2(kind='LS', weight = .5  )  # reg) #
        FPD.add_E3(kind='LS', weight = .05  )  #  reg) #
        FPD.add_ArcLengthApprox(kind='LS',weight = mig)
        #*********************************************
        nbvts = len(bareverts)
        ntot = thbcurve.n
        ukbar = np.linspace(0.,1.,ntot)[-1-nbvts-1:-1-1]
        #*********************************************
        # fair ini
        xb = thbcurve.vertices[0,0]
        yb = thbcurve.vertices[0,1]
        zb = thbcurve.vertices[0,2]
        # fair fini
        xe = thbcurve.vertices[-1,0]
        ye = thbcurve.vertices[-1,1]
        ze = thbcurve.vertices[-1,2]
        if nosetan:
            #*********************************************
            # nose tangent/curvature
            FPD.add_yFixity(index = 1,
                            value = yb,
                            track_prolongation=False)
            FPD.add_yFixity(index = 2,
                            value = yb,
                            track_prolongation=False)
            #--------------------------
            FPD.add_zFixity(index = 1,
                            value =zb,
                            track_prolongation=False)
            FPD.add_zFixity(index = 2,
                            value = zb,
                            track_prolongation=False)
        elif blendtan:
            pass
        else:
            #*********************************************
            # C2 continuity to bulb.
            FPD.add_xFixity(index = 1,
                            value = xb,
                            track_prolongation=False)
            FPD.add_xFixity(index = 2,
                            value = xb,
                            track_prolongation=False)
            #--------------------------
            FPD.add_yFixity(index = 1,
                            value = yb,
                            track_prolongation=False)
            FPD.add_yFixity(index = 2,
                            value = yb,
                            track_prolongation=False)
        #*********************************************
        # C2 continuity to midship.
        FPD.add_xFixity(index = -2,
                        value = xe,
                        track_prolongation=False)
        FPD.add_xFixity(index = -2,
                        value = xe,
                        track_prolongation=False)
        FPD.add_yFixity(index = -3,
                        value = ye,
                        track_prolongation=False)
        FPD.add_yFixity(index = -3,
                        value = ye,
                        track_prolongation=False)
        
        
        #for pt,loc in zip(bareverts,ukbar[ntinterp:]):
        for pt,loc in zip(bareverts,ukbar):
            #*********************************************
            FPD.add_xPointConstraint(kind= this_kind, 
                                     value=pt[0], 
                                     location=loc, 
                                     weight = big ) 
            FPD.add_yPointConstraint(kind= this_kind, 
                                     value=pt[1], 
                                     location=loc, 
                                     weight = big  )
            FPD.add_zPointConstraint(kind= this_kind, 
                                     value=pt[2], 
                                     location=loc, 
                                     weight = big  )
        #
        #*********************************************
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                            interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(thbcurve, L, 
                                         data = interval_data)
        #
        #*********************************************
        return Lspline
    
    def bulbschema(self):
        """defunct idea that needed
        variying numbers of tcurve vertices on the bulb
        and yet not 
        """
        bs = {'-1':7,
              '-2':9,
              '-3':11,
              '-4':13,
              '-5':15}
        mapb = {'-1':[0,1,2,3,4,5,6],
                '-2':[0,2,3,4,5,6,7,8,10],
                '-3':[0,1,2,3,4,5,6,7,8,9,10,11,12],
                '-4':13,
                '5':15}
        schema = [mapb,bs]
        return schema
    
    
    
    def setup_addbulb(self):
        """set bounds
        and make sure there is enough 
        resolution to add the bulb
        
        return
        ---------
        hullfwd : the finest fwd bare hull surface
        """
        hullfwd = self.barehull.hullfwd.get_finest_level()
        nu = len(self.bulb.tcurvenet)
        nv = len(self.bulb.lcurvenet)
        #this arbitrary guess
        #could be changed to depend on
        # 1.) the physical height of the bulb
        #       -which depends on the rules
        # 2.)  directly on the bulb-rules themselves
        
        #
        # could change to
        # def estimate_bulbfractionalheight
        #
        # from ComplexHull.py
        #
        while hullfwd.nu <2*nu or hullfwd.nv < 2*nv:
            print 'refining fwd surface for bulb'
            hullfwd = hullfwd.dyadic_loft_refinement()
        #
        #*****************************************************
        # hi resolution solve - 35x35... or whatever it is defined as
        # at init time.
        if self.make_hires: 
            hullfwd = hullfwd.dyadic_loft_refinement() 
        #
        #*****************************************************
        #
        # need extra pts in v to handle bulb vertices
        # need extra pts in u to handle complex interface
        # imagine an 11 tcurve fwd surface:  there would be no
        # leftover curves to define the upper portion!
        #
        # idea:  cut the surface at the degenerate boundary?
        #   (and an associated longi line "somewhere" in the degenerate set...)
            
        # better idea: just deform the bare hull to the bulb
        # meet C2 and keep both surfaces independent.
        #
        # bulb nu == number of curves in the v direction
        # bulb nv == number of curves in the u direction
        active_v_bounds = hullfwd.mastercurve_u.active_bounds_given_cv_index(nu)
        active_u_bounds = hullfwd.mastercurve_v.active_bounds_given_cv_index(nv)
        
        active_v_bounds.inf = 0.
        active_u_bounds.inf = 0.
        box1 = thbspline.Box( active_u_bounds, active_v_bounds )
        hullfwd.bounds = thbspline.BoxList(box1)
        return hullfwd
    
    
    
    def optimize_add_bulb(self, fs, tcurvenet):
        fs.tcurvelist[0].vertices #becomes the exterior boundary of the bulb
        # portion of the hull
        return
    
    def set_fwdboundary_complexhull(self, outline,fs):
        """Initial tcurve outline of the bbow complex hull
        and send the 
        """
        print '\n\n set_fwdboundary_complexhull\n\n'
        bhull = self.barehull.hullfwd.get_finest_surface()
        bbow = self.bulb.hullbulb
        
        
        nvts = len(outline)#number of vertices 
                            #in the original tcurve[0]
        
        
        #
        #*******************************************************
        #
        #
        #*******************************************************
        # intercept bare hull (refined) fwd tcurve with bbow top aft.
        bbow_top_aft_pt = bbow.tcurvelist[0].vertices[-1]
        yb = bbow_top_aft_pt[1]
        
        fwd_hull_curve = bhull.tcurvelist[0]
        fwd_hull_curve.bounds = BoundsList([ia(0.,1.)])
        
        p_intercept = fwd_hull_curve.FindPoint(x=yb,
                                               index=1,
                                               maxiter=200,
                                               tol=1.e-10)
        #
        #*******************************************************
        # Split the fwd bare hull tcurve
        lower_curve, upper_curve = fwd_hull_curve.CurveSplit(p_intercept[0])
        #
        #*******************************************************
        # 4 vertex curve emulating the upper_curve
        # now
        # expecting the degenerate 
        # interpolation of the entire bulb.
        vtxup = linear_vertices(upper_curve.vertices[0],
                               upper_curve.vertices[-1],
                               bhull.tcurvelist[0].n-nvts 
                               )
        #
        #*******************************************************
        #
        #
        #*******************************************************
        #
        
        #
        #*********************************************
        # bbow outline (bottom of fwd hull tcurve[0])
        fs.tcurvelist[0].vertices[:nvts] = outline[:]
        fs.vertices[:nvts,0]  = outline[:]
        for el in range(nvts):
            fs.lcurvelist[el].vertices[0] = outline[el]
        
        #
        #*********************************************
        # bbow highline
        fs.tcurvelist[0].vertices[nvts:] = vtxup[:]
        fs.tcurvelist[0].bounds.bounds[0].sup = 1.0
        fs.tcurvelist[0].compute_curve()
        
        fs.set_vertices_from_tcurves()
        fs.set_THB_curves_from_vertex_net(override_u=False,
                                                        override_v=True)
            
        # set hull vertices and tcurves based on lcurves:
#        for i in range(len(fs.lcurvelist)):
#            fs.lcurvelist[i].compute_curve()
#            fs.vertices[i,:] = fs.lcurvelist[i].vertices
#        
#        for i in range(fs.nv):
#            fs.tcurvelist[i].vertices[:] = fs.vertices[:,i]
#            fs.tcurvelist[i].compute_curve()
        fs.compute_surface()
        return 
    
    def set_fwdboundary_complexhull_old(self, outline,fs):
        """Initial tcurve outline of the bbow complex hull
        
        Does not know about the upper portion of the hull
        """
        nvts = len(outline)#number of vertices 
                            #in the original tcurve[0]
        
        #
        #*********************************************
        # bbow outline (bottom of fwd hull tcurve[0])
        fs.tcurvelist[0].vertices[:nvts] = outline[:]
        fs.vertices[:nvts,0]  = outline[:]
        for el in range(nvts):
            fs.lcurvelist[el].vertices[0] = outline[el]
            
        for i in range(len(fs.lcurvelist)):
            fs.lcurvelist[i].compute_curve()
            fs.vertices[i,:] = fs.lcurvelist[i].vertices
        
        for i in range(fs.nv):
            fs.tcurvelist[i].vertices[:] = fs.vertices[:,i]
            fs.tcurvelist[i].compute_curve()
        return 
    
    
    def plot_initial_situation(self,
                                 fancyish = True,
                                 color=False,
                                 simplefancy=True):
        """
        dev
        ----------
        fancyish = False
        color=False
        simplefancy=False
        
        
        fancyish = True
        color=True
        
        
        fancyish = False
        color=True
        
        
        
        fancyish = True
        simplefancy=True
        
        """
        ax = self.barehull.plot(fancyish,color,
                            simplefancy=simplefancy)
        ax = self.bulb.plot(fancyish,color,canvas=ax,
                            simplefancy=simplefancy)
        #for el in self.bulb.tcurvenet:
        #    ax = el.plotcurvehierarchy(canvas = ax)
        return
    
    def plotcomplex_hull(self,
                         fancyish = False,
                         color=False,
                         simplefancy=False,
                         bulbtcurves=True):
        """
        fancyish=True
        color = True
        bubtcurves=True
        simplefancy=True
        """
        ax = self.complexhull.plot(fancyish,color,
                                   simplefancy=simplefancy,
                                   plothullcurves=False)
        for curve in self.complexhull.hullcurves:
                ax.plot(curve.r[:,2], 
                        curve.r[:,0], 
                        curve.r[:,1], label='parametric curve', color = 'r')
        ax = self.bulb.plotbulbcurves(plottcurves = True,
                                      plotlcurves = False,
                                      canvas      = ax)
        ax = self.bulb.plot(fancyish,color,canvas=ax,
                            simplefancy=simplefancy)
        
        #boundary curve (not necesarily part of the hull)
        #ax = self.aft_boundary_curve.plotcurvehierarchy() #this one is part of the hull
        
        
        #ax = self.fwd_boundary_curve.plotcurvehierarchy(canvas=ax) 
        ax = self.barehull.hullfwd.tcurvelist[0].plotcurvehierarchy(canvas=ax) 
        return

class records(object):
    def __init__(self, curve):
        self.curve = curve
        self.indices = []
        self.bulbmap = None
        self.barehullfwd_tcurveindices = []
        self.vtinterp = None
        self.lfix = None
        self.nosetangent = False
        self.blendtan = False
        self.Lspline = None
        return

if __name__ == """__main__""":
    ch = ComplexHull()
    
    self = ch
    ch.make_system0()
    
    #method2:
    #ch.make_system()
    
    
    
    fs = self.complexhull.hullfwd.get_finest_surface()
    
    #method2:
    #fs = self.complexhull.newsurf
    
    
    self.plot_initial_situation(fancyish=True,
                              simplefancy=True)
    
    
    
    
    """
    self.make_bulb()
    self.bulb.hullbulb.compute_surface()
    """
    