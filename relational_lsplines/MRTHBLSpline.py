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
    for i in range(      1,nnew+1):   #  3 ,4):    #   
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


def tests():
    """quick test of cleanup of rbspline 
    (thbspline curves)
    
    dev
    --------
    
    the list, ll, computed in the sub function
    'verify_uneffected'
    is checked with total summations
    
    checking locally with a few 
    vertices instead, we have:
        
    ll[5][3]
    Out[152]: array([-0.08166504, 10.31408691])
    
    ll[0][3]
    Out[153]: array([-0.08166504, 10.31408691])
    
    ll[5][9]
    Out[154]: array([0.27685547, 6.08837891])
    
    ll[0][9]
    Out[155]: array([0.27685547, 6.08837891])
    
    
    
    what about after editing the fine levels?
    
    
    c2.restrict(use_active=False) - c1.vertices #good
    
    c2.restrict() - c1.vertices #good
     
    c1.prolong() - c2.vertices #good
    
    c2.vertices[12:20] += .3
    
    c2.vertices[12:20]
    Out[159]: 
    array([[1.10712891, 4.71552734],
           [1.33833008, 4.2206543 ],
           [1.59736328, 3.75654297],
           [1.88459473, 3.32282715],
           [2.20039062, 2.91914063],
           [2.54511719, 2.54511719],
           [2.91914062, 2.20039063],
           [3.32282715, 1.88459473]])
    
    c1.prolong() - c2.vertices
    >>>  array of zeros (using 1 level accurate prolongation)
    
    
    c2.restrict() - c1.vertices
    >>> array of zeros
    
    c2.restrict(use_active=False) 
    >>>
    array([[ 1.44250091e-03,  1.20014425e+01],
           [-4.51083090e-02,  1.14132250e+01],
           [-8.81422756e-02,  1.03024827e+01],
           [-8.22544401e-02,  8.74782368e+00],
           [ 5.42588730e-02,  7.35894637e+00],
           [ 2.45653380e-01,  6.05424713e+00],
           [ 6.02625875e-01,  4.93856337e+00],
           [ 1.43486036e+00,  4.31571974e+00],
           [ 1.80838219e+00,  3.24588219e+00],
           [ 2.57622381e+00,  2.57622381e+00],
           [ 3.27786487e+00,  1.84036487e+00],
           [ 3.80845660e+00,  9.27597222e-01],
           [ 4.97235451e+00,  6.36417010e-01],
           [ 6.04771222e+00,  2.39118470e-01],
           [ 7.36124756e+00,  5.65600581e-02],
           [ 8.74660327e+00, -8.34748525e-02],
           [ 1.03034145e+01, -8.72105037e-02],
           [ 1.14126620e+01, -4.56713618e-02],
           [ 1.20016772e+01,  1.67718702e-03]])
    
    c2.vertices[12:20] += .3
    vector = c0c.project_vertices_all_active()
    v2 = c1.prolong()
    vector - v2
    >>> vector of all zeros
    
    more succinct
    
    v = c0c.prolong()
    vv = c1.prolong(vector=v)
    vv - c0c.project_vertices_all_active()
    
    now go the other way
    
    vr1 = c2.restrict()
    vr0 = c1.restrict(vector = vr1)
    vr0 - c0c.vertices
    
    
    make operators to speed the process
    
    c2.prolong_vertices_to_this_level() - c2.vertices
    """
    tol = 1.e-14
    adj1 =  ia(0.125, 0.375)
    adj2 =  ia(0.1875, 0.25)
    adj3 =  ia(0.0, 0.25)
    
    c1.bounds = BoundsList([ia(0.01, 0.75)])
    
    c2.bounds = BoundsList([ia(0.03125, 0.09375),
                            ia(0.0625, 0.125),
                            ia(0.09375, 0.15625),
                            ])

    c2.bounds = BoundsList([ia(.1,.28),
                            ia(.3,.5),
                            ia(.51,.72)])
    
    self = c0c
    for u in np.linspace(0.,1.,30):
        print self(u)-self.eval_new(u) # old THB - new THB
        test = self.CurvePoint(u)-self.eval_new(u) 
        test2 = self.CurvePoint(u)-self.eval_new2(u) 
        print test #Bspline - new THB
        print test2 #Bspline - new THB
        print self.CurvePoint(u)-self(u) #Bspline - old THB
        assert(np.linalg.norm(test)<tol),'ERROR: curves dont match. diff= {}'.format(test)
        assert(np.linalg.norm(test2)<tol),'ERROR: curves dont match. diff= {}'.format(test2)
        print '--- --- --- ---'
        print 'curve evaluation is identical across methods'
        print '--- --- --- ---'
        
    c0c.plotcurve()
    c0c.plot_thbasis()
    
    def verify_uneffected():
        ll=[]
        for i in [0.,.2,.4,.6,.8,.99]:
            #ll.append(c0c.project_vertices_sensible(i))
            ll.append(c0c.project_vertices_all_active())
        ck1 = ll[0]
        for el in ll[1:]:
            print np.linalg.norm(ck1-el)
        print 'The difference in vertices is (now!!!)'
        print 'zero, no matter where they are '
        print 'computed'
        return
    verify_uneffected()
    
    
    def verify_thb_fine_edit():
        """here we can no longer check against the old B-spline
        But we may still compare THB evaluation methods
        against each other.
        
        AND
        
        Verify that is does not matter where one computes the 
        projected vertives from.
        """
        c2.vertices[12:20] += .3
        self = c0c
        for u in np.linspace(0.,1.,30):
            test =  self(u)-self.eval_new(u) # old THB - new THB
            #
            # this one is not going to be the same:
            test2 = self.children.children.CurvePoint(u)-self.eval_new2(u) 
            #unless c2.bounds = ia(0.,1.)
            print test #Bspline - new THB
            print test2
            assert(np.linalg.norm(test)<tol),'ERROR: curves dont match. diff= {}'.format(test)
            #assert(np.linalg.norm(test2)<tol),'ERROR: curves dont match. diff= {}'.format(test2)
            print '--- --- --- ---'
            
        print 'Now check constancy of vertices'
        ll=[]
        for i in [0.,.2,.3,.4,.5,.6,.8,.99]:
            ll.append(c0c.project_vertices_sensible(i))
        ck1 = ll[0]
        for el in ll[1:]:
            print np.linalg.norm(ck1-el)
        return
    verify_thb_fine_edit()
    return

if __name__ == """__main__""":
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
    ###
    ###
    ###
    #"""
    # LAUNCH THE MULTIGRID SOLVER:
    lc = issue_longitudinals(tcurvenet)
    #lc = uopt.issue_longitudinals(tcurvenet)
    #
    # failed stepwise longitudinal maker experiment
    #lc = uopt.issue_longitudinals_stepwise(
    #                                tcurvenet)
    #
    """
    nlist = [] 
    for Lspline in lc:
        this = Lspline.get_finest_Lspline()
        #******************************************
        print 'restore LS to equality'
        this.restore_Lagrangian()
        #    
        #******************************************
        this.optimize_thb(stop=5,
                          ceiling = 1.e-10,
                          relax_grad_bounds=False,
                          do_thb=True)
        #**************************************************
        nlist.append(Lspline)
        
    lcn = [Lspline.curve for Lspline in nlist]
    lcn.insert(0,original_lcurves[0])
    lcn.append(original_lcurves[-1])
    c.plot3DmultiList(tcurvenet,lcn)
    #"""
    #
    #
    #lcurvenet2 = [Lspline.curve.children for Lspline in lc]
    #lcurvenet2.insert(0,original_lcurves[0])
    #lcurvenet2.append(original_lcurves[-1])
    
    
    #lcurvenet3 = [Lspline.curve.get_coarsest_curve() for Lspline in lc]
    #lcurvenet3.insert(0,original_lcurves[0])
    #lcurvenet3.append(original_lcurves[-1])
    #
    
    #
    #
    lcurvenet = [Lspline.curve for Lspline in lc]
    #
    Lspline = lc[0]
    #
    # SD.hull.DWL
    #
    #"""
    """
    #*******************************************
    # use some example curves from ShipDesigner
    #
    hull1 = rbspline(SD.hull.CProfile.vertices,k,nump)
    lcurvenet.insert(0,hull1)
    #
    #
    hull2 = rbspline(SD.hull.DWL.vertices,k,nump)
    lcurvenet.append(hull2)
    #"""
    
    #*******************************************
    # make a full suite of longitudinal curves
    # by attaching the DWL and CProfile
    #"""
    lcurvenet.insert(0,original_lcurves[0])
    lcurvenet.append(original_lcurves[-1])
    #
    # PLOT EVERYTHING:
    c.plot3DmultiList(tcurvenet,lcurvenet)
    #c.plot3DmultiList(tcurvenet,lcurvenet2)
    #
    #c.plot3DmultiList(tcurvenet,lcurvenet3)
    #"""
    ##
    ##*************************************************************************
    ## Plan to make a THB surface from these longitudinals
    #
    #
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
    #       If      longitudinals are refined at multiple levels,
    #               (the spline curves space is way short of the full space)
    #               (and yet it IS enough data to span the first level)
    #               (and thus it IS enough data to prolongate to any level)
    #               (we just want to make that prolonged data come from)
    #               (the extra info on the finer levels of the Lcurves)
    #       THEN    Build the finest level first:
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
    #                   3.) Restrict the fine net to reconstruct
    #                       each coarser net
    #                       
    #                       -(Or use lower level lcurve details)
    #                       (to incorporate some differences down there??)
    #                       
    #
    
    ##
    ##*************************************************************************
    ##
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
                              hullcurvelist=tcurvenet,
                              ulevel = 2,
                              vlevel = 3)
    
    s1 = s0.dyadic_loft_refinement()
    s0.compute_surface()
    s0.plotSurface(start_hullcurve=0,
                   fancy=True,
                   colorize_curves=True,
                   colorize_vertices=True)
    #"""