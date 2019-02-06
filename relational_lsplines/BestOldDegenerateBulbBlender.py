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




def setup_interpolation_solver(thbcurve,
                               vtinterp,
                               this_kind='equality'):
    #                           this_kind='equality'):
    #                           this_kind='LS'):
    """setup THBcurve Lspline vertex interpolation solver
    
    
    inputs
    ----------
        thbcurve    = starting base level THB curve (n=4)
        vinterp     = the transverse vertices to be interpolated
                                        (including the endpoints!)
                                        
    dev
    ----------
    Lspline = setup_interpolation_solver(thbcurve,
                                         vtinterp)
    this_kind = 'equality'
    """
    big = 100000.
    interval_data, small = interval_bounds(thbcurve)
    FPD = FormParameterDict(thbcurve)
    
    FPD.add_E1(kind='LS', weight = 1.)
    FPD.add_E2(kind='LS', weight = .5)
    FPD.add_E3(kind='LS', weight = .05)
    FPD.add_ArcLength(kind='LS',
                          weight = big)
    ntinterp = len(vtinterp)
    ukbar = np.linspace(0.,1.,ntinterp)[1:-1]
    nth = 1
    
    start   = vtinterp[3]
    end     = vtinterp[6]
    dx      = np.linalg.norm(end-start)
    for pt,loc in zip(vtinterp[1:-1],ukbar):
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
        #*********************************************
        # X Fixing
#        if nth in [4,5,6,7]:
#            FPD.add_xFixity(index=nth,
#                            value=pt[0])
#        else:
#            FPD.add_xPointConstraint(kind= this_kind, 
#                                 value=pt[0], 
#                                 location=loc, 
#                                 weight = big ) 
#        #*********************************************
#        # Y Fixing
#        if nth in [3,4,5,6]:
#            FPD.add_xFixity(index=nth,
#                            value=pt[1])
#        else:
#            FPD.add_yPointConstraint(kind= this_kind, 
#                                 value=pt[1], 
#                                 location=loc, 
#                                 weight = big  )
#        #*********************************************
#        # Z Fixing
#        if nth in [3,4,5,6,7]:
#            if nth == 4:
#                #assert(pt[2]==start[2]),'You got START wrong!'
#                FPD.add_zFixity(index=nth,
#                                value=start[2]-.1*dx)
#            elif nth == 5:
#                FPD.add_zFixity(index=nth,
#                                value=start[2]-.05*dx)
#                
#            elif nth == 6:
#                FPD.add_zFixity(index=nth,
#                                value=end[2]+.05*dx)
#            elif nth == 7:
#                #assert(pt[2]==end[2]),'You got END wrong!'
#                FPD.add_zFixity(index=nth,
#                                value=end[2]+.1*dx)
#            #else:
#                #assert(1==0),'Logic Failure'
#        else:
#            FPD.add_zPointConstraint(kind= this_kind, 
#                                 value=pt[2], 
#                                 location=loc, 
#                                 weight = big  )
        #*********************************************
#        if nth in [3,4,5,6]:
#            FPD.add_CurvatureXoverZConstraint(value=0.,
#                                        location=loc,
#                                        weight=big)
        #*********************************************
#        if nth == 1:
#            FPD.add_zFixity(index=nth,
#                            value=0.)
#        else:
#            FPD.add_zPointConstraint(kind= this_kind, 
#                                 value=pt[2], 
#                                 location=loc, 
#                                 weight = big  )
            
        #*********************************************
        if nth ==3:
            FPD.add_CurvatureXoverZConstraint(kind=this_kind,
                                        location=loc,
                                        value = 0.,
                                        weight =big)
            FPD.add_AngleXoverZConstraint(kind=this_kind,
                                            location = loc,
                                            value = loc,
                                            weight=big)
#        if nth ==4:
#            FPD.add_CurvatureXoverZConstraint(kind=this_kind,
#                                        location=loc+.01,
#                                        value = 0.,
#                                        weight =big)
#        if nth ==5:
#            FPD.add_CurvatureXoverZConstraint(kind=this_kind,
#                                        location=loc-.01,
#                                        value = 0.,
#                                        weight =big)
#            FPD.add_CurvatureXoverZConstraint(kind=this_kind,
#                                        location=loc,
#                                        value = 0.,
#                                        weight =big)
        if nth ==4:
#            FPD.add_CurvatureXoverZConstraint(kind=this_kind,
#                                        location=loc,
#                                        value = 0.,
#                                        weight =big)
            #-----------------------------
#            test = FPD.curve.greville_abscissa(7)
#            FPD.add_CurvatureXoverZConstraint(kind=this_kind,
#                                        location=loc-.1,
#                                        value = 0.,
#                                        weight =big)
#            FPD.add_CurvatureXoverZConstraint(kind=this_kind,
#                                        location=loc,
#                                        value = 0.,
#                                        weight =big)
            FPD.add_CurvatureXoverZConstraint(kind=this_kind,
                                        location=loc,
                                        value = 0.,
                                        weight =big)
#            FPD.add_CurvatureXoverZConstraint(kind=this_kind,
#                                        location=test,
#                                        value = 0.,
#                                        weight =big)
            #---- --- --- --
#            FPD.add_AngleXoverZConstraint(kind=this_kind,
#                                            location = loc,
#                                            value = 0.,
#                                            weight=big)
#            FPD.add_AngleXoverZConstraint(kind=this_kind,
#                                            location = loc-.01,
#                                            value = 0.,
#                                            weight=big)
            FPD.add_AngleXoverZConstraint(kind=this_kind,
                                            location = loc,
                                            value = 0.,
                                            weight=big)
        nth +=1
    FPD.add_verticalXoverZAngleConstraint(kind=this_kind,
                                            location = 0.,
                                            value = 0.,
                                            weight=big)

#    FPD.add_AngleXoverZConstraint(kind      = this_kind,
#                                 location   = 0.,
#                                 value      = 90.,
#                                 weight     = 10*big)
    L = Lagrangian(FPD)
    interval_data, small = lagrangian_bounds(L, 
                        interval_data, small, 1.e4)
    Lspline = IntervalLagrangeSpline(thbcurve, L, 
                                     data = interval_data)
    return Lspline



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

def THBsolver(thbcurve,
              vtinterp,
              Lspline=None):
    """interntion : demonstrate an adaptive solver
    for THB-spline optimization of a longitudinal 
    hull curve!
    
    inputs
    ----------
        thbcurve    = starting base level THB curve (n=4)
        vinterp     = the transverse vertices to be interpolated
                                        (including the endpoints!)
    process
    ----------
        1.  solve one iteration on the base level
        2.  check the gradient
        
    notes
    ----------
        11 control vertices and 6 interpolation vertices
        makes for 51 elements of the gradient
        11*3 + 6*3 = 51
    """
    print 'Local THBsolver'
    #because we are using equality constraints
    # we must begin with a soluble number of vertices
    # to get any information out of our solver!
    #**************************************************
    # refine until the ini curve is feasible 
    # given the number of vertices to 
    # be interpolated
    #
    #**************************************************
    #
    while nthb_vs_nvi(thbcurve,vtinterp):
        thbcurve = thbcurve.dyadic_refinement()
    #**************************************************
    thbcurve.parent = None
    #
    print '   setup'
    if Lspline is None:
        Lspline = setup_interpolation_solver(thbcurve,
                                             vtinterp)
    #******************************************
    print 'switch equality to LS'
    Lspline.switch_all_of_one_Lagrangian_kind(
            ini_kind='equality',
            fini_kind='LS')
    #******************************************
    print '   optimize'
    Lspline.optimize_thb(stop=1,
                         ceiling = 1.e-3,
                         relax_grad_bounds=False)
    #******************************************
    print 'restore LS to equality'
    Lspline.restore_Lagrangian()
    #    
    #******************************************
    # manually call refine:
    # auto-revert back to 
    # original intentions for constraints
    Lspline.refine() #need to check if there is anything to refine yet
    
    this = Lspline.get_finest_Lspline()
    #    
    #******************************************
    print 'switch equality to LS'
    this.switch_all_of_one_Lagrangian_kind(
            ini_kind='equality',
            fini_kind='LS')
    #******************************************
    cmax = 1
    count =0
    while this is not None and count<cmax:
        print '   optimize'
        this.tol = 1.e-4
        #******************************************
        #******************************************
        # April 4
        #******************************************
        #do some smoothing:
        print 'switch equality to LS'
        this.switch_all_of_one_Lagrangian_kind(
                ini_kind='equality',
                fini_kind='LS')
        #******************************************
        print '   optimize'
        this.optimize_thb(stop=7,
                          ceiling = 1.e-3,
                          relax_grad_bounds=False)
        #******************************************
        print 'restore LS to equality'
        this.restore_Lagrangian()
        #******************************************
        print '   optimize'
        this.optimize_thb(stop=5,
                          ceiling = 1.e-3,
                          relax_grad_bounds=False)
        #******************************************
        # End April 4
        #**************************************************
        # manually call refine:
        if count<1:
            this.refine()
        this = Lspline.get_finest_Lspline()
        #******************************************
        count +=1
    #
    #**************************************************
    #print '   optimize'
    #this.parent.optimize_thb(stop=20,
    #                     relax_grad_bounds=True)
    #print '   check for refinement'
    #Lspline = set_refinement(Lspline,
    #                         vtinterp)
    #print '   optimize'
    #Lspline.optimize_thb(stop=5,
    #                     relax_grad_bounds=True)
    print 'Done with this longitudinal'
    return this.get_coarsest_Lspline()



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
        #lcurvenet.append(THBsolver(thbcurve,vtinterp))
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
        self.hull = None
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
    
    
    def make_surface(self):
        """see opt_simple_hull.py
        function: def make_THB_bulb
        """
        self.hullbulb = make_surface(self.lcurvenet,
                                self.tcurvenet,
                                ul=2,vl=1,
                                refine=False)
        return
    
    
    
    def return_tcurve(self,index):
        return self.tcurvenet[index]
    
    def return_outline_vertices(self,num2use=3):
        #if isinstance(use, int):
        #    use = [use]
        #else:
        #    use.sort()
        #------
        # if we are using 2 curves, use the fwd most 2,
        # which translates into the 'last' 2 tcurves in the bbow net
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
    
    def mask_verts(self, num2use=3):
        mask=[]
        mask = np.array(mask, dtype=bool)
        return mask
    
    def return_fwd_boundary(self):
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
        lc = uopt.issue_split_longitudinals(tcurvenet,
                                            cpk=original_lcurves[0],
                                            dwl=original_lcurves[-1])
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
                            curve.r[:,1], label='parametric curve', color = 'r')
            except:
                mssg = 'Warning, hull curves are not attached to '+\
                            'this hull object, and you tried to plot them!'
                print mssg
        return ax
        

class ComplexHull(object):
    def __init__(self):
        self.barehull       = BareHull()
        self.bulb           = Bulb()
        
        self.hullcurves2    = get_tcurves()
        #
        #**********************************************************************
        # hi res:  solve the bulb longitudinals on the 35x35 level
        self.make_hires = False#True
        #
        #**********************************************************************
        # 
        self.num2use    =  5 #number of  transverse bulb curves 
                            #to use for the bulb
        #
        #**********************************************************************
        # number of transition curves added to the bulb solve:
        self.num_extra_boundary=-1 #0 #3  #all of them, always
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
    def make_bulb(self):
        self.bulb.make_surface()
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
        #self.original_barehullfwd = copy.deepcopy(self.barehull.hullfwd) #slooooooow
        #
        #**********************************************************************
        # get the finest, fine (enough) level bare hull fwd -> complex hull to be
        #
        fs = self.setup_addbulb()  #dyadic refinement of the bare hull fwd surface
        self.fsini = fs
        #
        #
        #**********************************************************************
        # exterior tcurve boundary vertices of the bulb, 
        # use the fwdmost 3 transverse bulb
        # curves in this initial experiment
        outline = self.bulb.return_outline_vertices(num2use = numbcs) 
        #
        # add a connective vertex or two or threee back on the tcurve
        index=0
        while fs.tcurvelist[0].vertices[index,1] < outline[-1,1]:
            index+=1
        outline = list(outline)
        #
        #**********************************************************************
        # transition curves added to the bulb solve:
        #self.num_extra_boundary=3 moved to init
        #
        #**********************************************************************
        #
        #if self.num_extra_boundary == -1:
        #    self.num_extra_boundary = fs.tcurvelist[0].n
        #for i in range(self.num_extra_boundary):
        #    outline.append(fs.tcurvelist[0].vertices[index+i])
        for i in range(index, self.num_extra_boundary):
            outline.append(fs.tcurvelist[0].vertices[i])
        outline = np.asarray(outline)
        self.tnum = len(outline)
        #self.tnum = fs.tcurvelist[0].n
        self.outline = outline
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
        for i in range(tnum+self.num_extra_boundary):
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
        
        self.solver_records = {}
        
        
        
        
        
        
        
        # loop longitudinals to solve
        self.solver_records = {}
        #for i in range(tnum): #each i is to be a solve
        for i in range(fs.lcurvelist[0].n):
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
                    #bottom lsplines catch the 0th vertex of the elth bulb tcurve
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
            vtinterp = np.asarray(vtinterp)
            self.solver_records[i].vtinterp = vtinterp #vertices to interpolate
            self.solver_records[i].lfix = FPLongFix #fix interior bulb lcurve vertices
            self.solver_records[i].lboundary = FPLongboundary #fix bulb boundary
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
            for curve in FPcurvenet:
                bvts.append(curve.vertices[i])
            bvts = np.asarray(bvts)
            #--------------------------------------------------
            vl          = list(self.solver_records[i].vtinterp)
            #
            vl.append(vfini)
            vl.insert(0,vini)
            #
            # initialize curve
            vlarray = np.asarray(vl)
            linverts = linear_vertices(vl[0],vl[-1],4)
            curve = rbspline(linverts,k=4,nump=30)
            curve = curve.dyadic_refinement() #5
            curve = curve.dyadic_refinement() #7
            curve = curve.dyadic_refinement() #11 #so there is enough to catch all the degenerate bounds
            curve.parent = None
            # initialize solver and Lagrangian
            Lspline = self.setup_fwd_bulb_interpolation_solver(thbcurve = curve,
                                                               vtinterp = vtinterp,
                                                               lfix = lfix,
                                                               lbounds = lbounds,
                                                               bareverts = bvts,
                                                               nosetan = nosetangent)
            # use THB solver process
            Lspline = uopt.THBsolver(thbcurve=curve,
                                     vtinterp=vlarray,
                                     Lspline=Lspline,
                                     normalize=False)
            # done
            self.complex_fwd_Lnet.append(Lspline)
            self.complex_fwd_lcurves.append(Lspline.curve)
            
        ##
        ##****************************************************************
        ##
        fine_curves = []
        #----------------------------------------
        if self.make_hires: #ad hoc: hi res is 35 vertex, i.e. dyadic level 5 [0:4,1:5,2:7,3:11,4:19,5:35]
            for curve in self.complex_fwd_lcurves:
                this = curve.get_finest_curve()
                while this.n < 35:
                    this = this.dyadic_refinement()
                fine_curves.append(this)
        else: #ad hoc: low res is "19" i.e. dyadic level 4
            for curve in self.complex_fwd_lcurves:
                this = curve.get_finest_curve()
                while this.n < 19:
                    this = this.dyadic_refinement()
                fine_curves.append(this)
        #-----------------------------------------
        
        self.complexhull = self.barehull
        fs = self.complexhull.hullfwd.get_finest_surface()
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
        #
        #
        fs.bounds[0][0].sup = .5 #.8125 #.75 #.6875 #.46875 #.4375 #.40625#.375#.34375#3125
        fs.bounds[0][1].sup = .15
        #
        #"""
        self.endvtx_index = 2*self.num2use+self.bulb.hullbulb.nu -2
        bounds = fs.mastercurve_u.get_bounds_given_index(self.endvtx_index)
        
        box1 = thbspline.Box( ia(.0,.8125),    
                             ia(.0,.34375) ) 
        #"""
        #
        fs.compute_surface()
        self.complexhull.hullmid.numpu = 20
        self.complexhull.hullmid.numpv = 5
        self.complexhull.hullmid.compute_surface()
        self.complexhull.hullfwd.numpu = 100
        self.complexhull.hullfwd.numpv = 100
        self.complexhull.hullfwd.compute_surface()
        
        self.plotcomplex_hull(fancyish=True,color=True,simplefancy=True)
        #
        #
        """
        
        self.barehull.hullcurves = get_tcurves()
        
        
        
        c = self.complex_fwd_lcurves[0]
        ax = c.plotcurvehierarchy()
        for cv in self.complex_fwd_lcurves:
            ax = cv.plotcurvehierarchy(canvas=ax)
            
        for el in self.barehull.hullcurves:
            ax = el.plotcurvehierarchy(canvas = ax)
        
        
        for el in self.barehull.hullfwd.tcurvelist:
            ax = el.plotcurvehierarchy(canvas = ax)
        for el in self.barehull.hullmid.tcurvelist:
            ax = el.plotcurvehierarchy(canvas = ax)
        for el in self.barehull.hullaft.tcurvelist:
            ax = el.plotcurvehierarchy(canvas = ax)
            
        for el in self.bulb.tcurvenet:
            ax = el.plotcurvehierarchy(canvas = ax)
            
            
            
        for i in range(16):
            record = self.solver_records[i]
            print i,len(record.lfix)
            
        i=10
        print self.solver_records[i]
        
        
        c = self.barehull.lcurvenet_fwd[0]
        cf = c.get_finest_curve()
        cf.n #35
        
        bh = self.barehull.hullfwd.get_finest_surface()
        bh.nu #35
        bh.nv # 35
        
        
        fs = self.complexhull.hullfwd.get_finest_surface()
        #>>> fs.bounds
        #>>> Out[4]: BoxList((Box((ia(0.0, 0.25), ia(0.0, 0.1875))),))
        fs.bounds[0][0].sup =.5
        fs.bounds[0][1].sup = .5
        #
        
        fs.bounds[0][0].sup =.25
        fs.bounds[0][1].sup = .25
        #
        
        self.complexhull.hullfwd.numpu = 30
        self.complexhull.hullfwd.numpv = 30
        self.complexhull.hullmid.numpu = 15
        self.complexhull.hullmid.numpv = 15
        self.complexhull.hullmid.compute_surface()
        #fs.numpu = 30
        #fs.numpv = 30
        self.complexhull.hullfwd.numpu = 20
        self.complexhull.hullfwd.numpv = 20
        #
        
        self.complexhull.hullfwd.compute_surface()
        self.plotcomplex_hull(fancyish=True,color=True,simplefancy=True)
        """
        return
    
    
    def check_solver(self, fs=None):
        if fs is None:
            fs = self.complexhull.hullfwd.get_finest_surface()
            
        if self.hullcurves2 is None:
            self.hullcurves2 = get_tcurves()
        #aftbdry = fs.tcurvelist[-1].vertices
        #lck = []
        #        for i in range(self.tnum):
        #            vl=list(self.solver_records[i].vtinterp)
        #            vini = self.fwd_boundary_curve.vertices[i]
        #            vfini = self.aft_boundary_curve.vertices[i]
        #            vl.append(vfini)
        #            vl.insert(0,vini)
        #            lck.append(rbspline(np.asarray(vl),k=4,nump=30))
        
        # setup curves:
        #ax = self.aft_boundary_curve.plotcurvehierarchy()
        #ax = self.fwd_boundary_curve.plotcurvehierarchy(canvas=ax)
        #for i in range(self.tnum):
        #    ax = lck[i].plotcurvehierarchy(canvas=ax)
            
        
        # solved curves:
        ax = self.aft_boundary_curve.plotcurvehierarchy()
        ax = self.fwd_boundary_curve.plotcurvehierarchy(canvas=ax)
        for i in range(self.tnum):
            #fs.lcurvelist[i].compute_curve()
            curve = rbspline(fs.lcurvelist[i].vertices,k=4,nump=30)
            #ax = fs.lcurvelist[i].plotcurvehierarchy(canvas = ax)
            ax = curve.plotcurvehierarchy(canvas = ax)
        
        
        ax = self.aft_boundary_curve.plotcurvehierarchy()
        ax = self.fwd_boundary_curve.plotcurvehierarchy(canvas=ax)
        for lcurve in fs.lcurvelist:
            #lcurve.compute_curve()
            curve = rbspline(lcurve.vertices,k=4,nump=30)
            ax = curve.plotcurvehierarchy(canvas=ax)
        for tcurve in self.hullcurves2:
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
        ntinterp = len(vtinterp)
        nbvts = len(bareverts)
        ukbar = np.linspace(0.,1.,ntinterp+nbvts+2)[1:-1]
        nth = 1
        totn = ntinterp
        xb = thbcurve.vertices[0,0]
        yb = thbcurve.vertices[0,1]
        zb = thbcurve.vertices[0,2]
        #*********************************************
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
        for pt,loc in zip(bareverts,ukbar[ntinterp:]):
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
        #"""
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
        
        while hullfwd.nu <2*nu or hullfwd.nv < 2*nv:
            print 'refining fwd surface for bulb'
            nlevel = hullfwd.dyadic_loft_refinement()
            hullfwd = nlevel
        #
        #*****************************************************
        # hi resolution solve - 35x35
        if self.make_hires: 
            hullfwd = hullfwd.dyadic_loft_refinement() #35x35 surface
        # or don't and instead try with 19x19 surface.
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
        lower_curve, upper_curve = fwd_hull_curve.CurveSplit(
                                                p_intercept[0])
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
        
        
        if self.num_extra_boundary == -1:
            self.num_extra_boundary = bhull.tcurvelist[0].n - \
                                                    len(outline)
        fs.compute_surface()
        return 
    
    def set_fwdboundary_complexhull_old(self, outline,fs):
        """Initial tcurve outline of the bbow complex hull
        and send the 
        """
        nvts = len(outline)
        fs.tcurvelist[0].vertices[:nvts] = outline[:]
        fs.vertices[:nvts,0]  = outline[:]
        for el in range(nvts):
            fs.lcurvelist[el].vertices[0] = outline[el]
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
        
        #boundary curve (not necesarily part of the hull)
        #ax = self.aft_boundary_curve.plotcurvehierarchy() #this one is part of the hull
        ax = self.fwd_boundary_curve.plotcurvehierarchy(canvas=ax) 
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
        return

if __name__ == """__main__""":
    ch = ComplexHull()
    ch.make_bare_hull_3part() #make the bare hull
    ch.make_bulb()
    ch.plot_initial_situation(fancyish=True,
                              simplefancy=True)
    
    self = ch
    
    self.add_bulb()
    fs = self.complexhull.hullfwd.get_finest_surface()
    
    oldexperiments = False
    if oldexperiments:
        k=4
        nump=100
        start = [0.,12.]
        end = [12.,0.]
        num = 4
        vv = spline.curvalinear_vertices(start,end,num)
    
        c0c = rbspline(vv,k,nump)
    
        
        c0c = c0c.dyadic_refinement()
        c0c = c0c.dyadic_refinement()
        c0c = c0c.dyadic_refinement()
        
        c0c.parent = None #TLM crucial to make it less deep a tree!
        c1 = c0c.dyadic_refinement()
        c2 = c1.dyadic_refinement()
        
        
        #
        #**************************************************************************
        #
        #Xi = get_test_vertices()
        #T0,T1,T2 = get_test_thb_vertices()
        #
        #**************************************************************************
        # get transverse curves to interpolate:
        tcurvenet = get_tcurves()
        c = tcurvenet[0]
        #
        original_lcurves = get_lcurves()
        #"""
        # LAUNCH THE MULTIGRID SOLVER:
        lc = uopt.issue_split_longitudinals(tcurvenet,
                                            cpk=original_lcurves[0],
                                            dwl=original_lcurves[-1])
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
        c.plot3DmultiList(tcurvenet,lcurvenet_fwd+lcurvenet_mid+lcurvenet_aft) 
        #
        ll = tcurvenet+lcurvenet_fwd+lcurvenet_mid+lcurvenet_aft
        ax = c.plotcurvehierarchy()
        for el in ll:
            ax = el.plotcurvehierarchy(canvas=ax)
        curve = lc[0][0]
        ##
        ##*************************************************************************
        ##
        
        
        
        sfwd = make_surface(lcurvenet_fwd,
                            tcurvenet_fwd,
                            ul=2,vl=2)
        smid = make_surface(lcurvenet_mid,
                            tcurvenet_mid,
                            ul=2,vl=0)
        saft = make_surface(lcurvenet_aft,
                            tcurvenet_aft,
                            ul=2,vl=2)
        
        fancyish = False
        #fancyish = True
        ax = sfwd.plotSurface(
                        start_hullcurve=0,
                        fancy=fancyish,
                        colorize_curves=True,
                        colorize_vertices=True)
        
        ax = smid.plotSurface(
                        start_hullcurve=0,
                        fancy=fancyish,
                        colorize_curves=True,
                        colorize_vertices=True,
                        canvas = ax)
        
        ax = saft.plotSurface(
                        start_hullcurve=0,
                        fancy=fancyish,
                        colorize_curves=True,
                        colorize_vertices=True,
                        canvas = ax)