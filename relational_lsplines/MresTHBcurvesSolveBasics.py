#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 05:44:39 2018

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


from MRTHBLSpline import get_tcurves, get_lcurves


def init_Lspline(thbcurve, FPD):
    """setup THBcurve Lspline solver
    parameters
    ----------
    thbcurve : a myspline.rbspline THB curve
    FRP : form parameter design constraint object
    
    dev
    ----------
    thbcurve = curve
    """
    interval_data, small = interval_bounds(thbcurve)

    L = Lagrangian(FPD)
    interval_data, small = lagrangian_bounds(L, 
                        interval_data, small, 1.e4)
    Lspline = IntervalLagrangeSpline(thbcurve, L, 
                                     data = interval_data)
    return Lspline

def THBsolver(Lspline):
    #**************************************************
    #
    print 'switch equality to LS'
    Lspline.switch_all_of_one_Lagrangian_kind(
                                            ini_kind='equality',
                                            fini_kind='LS')
    #**************************************************
    #
    print '   optimize'
    Lspline.optimize_thb(stop=5,
                         ceiling = 1.e-18,
                         relax_grad_bounds=False)
    print Lspline.refine_by_these_grad_indices
    #Lspline = Lspline.get_finest_Lspline()
    #**************************************************
    # RESTORE to EQUALITY CONSTRAINTS 
    print 'restore LS to equality'
    Lspline.restore_Lagrangian()
    print 'again = ',Lspline.refine_by_these_grad_indices
    #**************************************************
    # refine:
    Lspline.refine() 
    #Lspline.refine_by_force(
    #        boundslst = [ia(0.,.2),ia(.8,1.)]) 
    #**************************************************
    #
    this = Lspline.get_finest_Lspline()
    print 'this.curve.bounds = ',this.curve.bounds
    #**************************************************
    #******************************************
    #do some smoothing:
    print 'switch equality to LS'
    this.switch_all_of_one_Lagrangian_kind(
                                    ini_kind='equality',
                                    fini_kind='LS')
    #******************************************
    print '   optimize'
    this.optimize_thb(stop=8,
                      ceiling = 1.e-10,
                      relax_grad_bounds=False)
    #this = Lspline.get_finest_Lspline()
    #******************************************
    #**************************************************
    cmax = 3
    count =0
    while this is not None and count<cmax:
        print '   optimize'
        #******************************************
        print 'restore LS to equality'
        this.restore_Lagrangian()
        this.tol = 1.e-7
        #******************************************
        print '   optimize equality constraints'
        this.optimize_thb(stop=20,
                          ceiling = 1.e-10,
                          relax_grad_bounds=False)
        #this = Lspline.get_finest_Lspline()
        #******************************************
        ckit = uopt.ckLspline(this)
        if ckit and this.conv<this.tol:
            setconv = True
            break
        #******************************************
        # refine:
        if count<cmax:
            this.refine()
        this = Lspline.get_finest_Lspline()
        this.verbose = True
        #******************************************
        #do some smoothing:
        print 'switch equality to LS'
        this.switch_all_of_one_Lagrangian_kind(
                ini_kind='equality',
                fini_kind='LS')
        #******************************************
        print '   optimize'
        this.optimize_thb(stop=5,
                          ceiling = 1.e-3,
                          relax_grad_bounds=False)
        #this = Lspline.get_finest_Lspline()
        #******************************************
        #print 'restore LS to equality'
        #this.restore_Lagrangian()
        setconv = False
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
    if setconv:
        print 'convergence is true'
    else:
        print 'convergence is false'
    return this.get_coarsest_Lspline()
    
    
def make1():
    """
    testing
    ----------
    dev
    ----------
    stop=5
    ceiling = 1.e-12,
    relax_grad_bounds=False
    
    self = Lspline
    Lagrangian = self.Lagrangian
    L = self.Lagrangian
    vertices = self.curve.ADvertices
    
    
    vtx = Lspline.curve.project_vertices_all_active()
    
    print Lspline.child.Lagrangian.equality[3].type
    fp = Lspline.child.Lagrangian.equality[3]
    
    Lspline.child.curve.compute_area(vertices=vtx)
    """
    k=4
    nump=30
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
    
    curve = initial_curve([xb,yb],[xe,ye],4,nump=200) 
    curve = uopt.make_thb(curve)
    #curve.parent = None
    #curve = curve.dyadic_refinement()
    #curve = curve.dyadic_refinement()
    #curve = curve.dyadic_refinement()
    curve.parent = None
    
    FPD = FormParameterDict(curve) 
    
    
    #
    #**************************************************************************
    # derivative
#    FPD.add_verticalAngleConstraint(kind='equality', 
#                                    location = 0., value = 0.)
    #FPD.add_AngleConstraint(kind='equality', 
    #                        location = .5, value = 45.)
#    FPD.add_AngleConstraint(kind='equality', 
#                            location = 1., value =  0.)
    # curvature
#    FPD.add_CurvatureConstraint(kind='equality',
#                                location=.0,
#                                value = 0.,
#                                weight=1.)
    #
    #**************************************************************************
    # integral
#    FPD.add_AreaConstraint(kind='equality',
#                           value=43.,
#                           weight=1.)
#    FPD.add_XcConstraint(kind       = 'equality',
#                           value    = 3.58,
#                           weight   = 1.)
    #
    #**************************************************************************
    # positional
#    FPD.add_xPointConstraint(kind       = 'equality',
#                             location   = .5,
#                             value      = 3.5)
#    FPD.add_yPointConstraint(kind       = 'equality',
#                             location   = .5,
#                             value      = 5.)
#    FPD.add_yPointConstraint( kind      = 'equality', 
#                             location   = .15, 
#                             value      = 10.)
    #FPD.add_xFixity(value=0. , index=1, track_prolongation=False)
    #FPD.add_xFixity(value=0. , index=2, track_prolongation=False)
    #
    FPD.add_xFixity(value=1. , index=2, track_prolongation=False)
    FPD.add_yFixity(value=1. , index=2, track_prolongation=True) #last vertices might be trickier
    #
    #**************************************************************************
    # fairness
    FPD.add_E1(kind='LS', weight = .5)
    FPD.add_E2(kind='LS', weight = .5)
    FPD.add_E3(kind='LS', weight = .5)
    FPD.add_ArcLengthApprox(kind='LS', weight = .5)
    #
    #**************************************************************************
    # setup
    #
    Lspline = init_Lspline(curve, FPD)
    #
    #Lspline.optimize_thb(stop = 50,
    #                     ceiling = 1.e-12,
    #                     relax_grad_bounds=False)
    Lspline.verbose = True
    #Lspline = THBsolver(Lspline=Lspline)
    #
    #**************************************************************************
    # solve
    Lspline = uopt.THBsolver(Lspline=Lspline,
                             refinemax=4, 
                             maxlevel=3,
                             normalize=True)
    return Lspline


    
def make2():
    """
    testing
    ----------
    dev
    ----------
    stop=5
    ceiling = 1.e-12,
    relax_grad_bounds=False
    
    self = Lspline
    Lagrangian = self.Lagrangian
    L = self.Lagrangian
    vertices = self.curve.ADvertices
    
    
    vtx = Lspline.curve.project_vertices_all_active()
    
    print Lspline.child.Lagrangian.equality[3].type
    fp = Lspline.child.Lagrangian.equality[3]
    
    Lspline.child.curve.compute_area(vertices=vtx)
    """
    this_kind = 'equality'
    big = 100000.  #1. #
    mig = big/2.  #.5 #
    k=4
    nump=30
    xb = 0.
    yb = 17.33641345
    zb = 0.
    
    xe = 14.15239481
    ye = 17.33641345
    ze = 49.19157975
    
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
    #
    curve = initial_curve([xb,yb,zb],
                          [xe,ye,ze],
                          11,nump=30) 
    curve = uopt.make_thb(curve)
    #curve.parent = None
    #curve = curve.dyadic_refinement()
    #curve = curve.dyadic_refinement()
    #curve = curve.dyadic_refinement()
    curve.parent = None
    #
    FPD = FormParameterDict(curve) 
    #
    vtinterp = [np.asarray([12.35339247, 17.33641345, 17.70896871]),
                np.asarray([14.04681792, 17.33641345, 36.89368482])]
    ntinterp = len(vtinterp)
    ukbar = np.linspace(0.,1.,ntinterp+2)[1:-1]
    #
    count = 0
    for pt,loc in zip(vtinterp,ukbar):
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
        count += 1
    #
    FPD.add_yFixity(index = 1,
                    value = yb,
                    track_prolongation=False)
    #FPD.add_yFixity(index = 2,
    #                value = yb,
    #                track_prolongation=False)
    #--------------------------
    FPD.add_zFixity(index = 1,
                    value =zb,
                    track_prolongation=False)
    #
    # C2 continuity to midship.
    FPD.add_xFixity(index = -2,
                    value = xe,
                    track_prolongation=False)
    FPD.add_yFixity(index = -2,
                    value = ye,
                    track_prolongation=False)
    #--------------------------------------------------
    #    FPD.add_zFixity(index = -2,
    #                    value = ze,
    #                    track_prolongation=False)
    #--------------------------------------------------
    #FPD.add_xFixity(index = -3,
    #                value = xe,
    #                track_prolongation=False)
    #FPD.add_yFixity(index = -3,
    #                value = ye,
    #                track_prolongation=False)
    #    FPD.add_CurvatureXoverZConstraint(  kind        = this_kind,
    #                                        location    = 1.,
    #                                        value       = 0.,
    #                                        weight      = mig)
    #    FPD.add_CurvatureConstraint(        kind        = this_kind,
    #                                        location    = 1.,
    #                                        value       = 0.,
    #                                        weight      = mig)
    #
    #**************************************************************************
    # fairness
    FPD.add_E1(kind='LS', weight = 1.)
    FPD.add_E2(kind='LS', weight = .5)
    FPD.add_E3(kind='LS', weight = .05)
    FPD.add_ArcLengthApprox(kind='LS', weight =  mig) #1.5)  # 
    #
    #**************************************************************************
    # setup
    #
    Lspline = init_Lspline(curve, FPD)
    #
    #Lspline.optimize_thb(stop = 50,
    #                     ceiling = 1.e-12,
    #                     relax_grad_bounds=False)
    Lspline.verbose = True
    #Lspline = THBsolver(Lspline=Lspline)
    #
    #**************************************************************************
    # solve
    #Lspline.optimize()
    Lspline = uopt.THBsolver(Lspline=Lspline,
                             refinemax=4, 
                             maxlevel=3,
                             normalize=False)
    
    print ukbar
    print vtinterp[0]
    print vtinterp[1]
    print 'ukbar 0'
    print Lspline.curve(ukbar[0])
    print 'ukbar 1'
    print Lspline.curve(ukbar[1])
    print ''
    print 'diff 0'    
    print Lspline.curve(ukbar[0])[0] - vtinterp[0]
    print ''
    print 'diff 1'
    print Lspline.curve(ukbar[1])[0] - vtinterp[1]
    return Lspline

if __name__ == """__main__""":
    #Lspline = make1()
    Lspline = make2() #fwd upper hull test - can it hit the tcurve vertices?
    Lspline.curve.plotcurvehierarchy()
    self = Lspline.get_finest_Lspline()
    vtx = self.curve.project_vertices_all_active()
    print 'area = ',self.curve.area
    print 'Xc = ',self.curve.Xc
    print 'tan(0.) = ',self.curve.compute_tangent(0.,vtx)
    print 'tan(1.) = ',self.curve.compute_tangent(1.,vtx)
    print 'curvature(0.) = ',self.curve.compute_curvature(0.,vtx)
    print 'curvature(1.) = ',self.curve.compute_curvature(1.,vtx[:,0::2])
    print 'curvature(1.) = ',self.curve.compute_curvature(1.,vtx[:,1:])