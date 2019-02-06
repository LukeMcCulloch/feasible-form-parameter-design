#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 13:54:01 2017

@author: luke
"""
import sqKanren as lp #logic programming, etc.
np = lp.np #numpy
ia = lp.ia #interval arithmetic
import copy
#
#Optimization Machinery:
from initialValues import InitializeControlVertices
from initialValues import interval_bounds, \
                          lagrangian_bounds
from plots         import Plotter
from adials_gui3d  import frame, \
                          DrawCurveInteractive, \
                          hull_frame
from ADILS         import IntervalLagrangeSpline, \
                          Lagrangian
from FormParameter import FormParameterDict
#
#Basic Bspline Machinery
import curve as spline
from myspline import rbspline, BoundsList
import thbsurface as thbspline

verbose = True


    
##
##*****************************************************************************
## Check an Lspline solve:
def ckLspline(Lspline):
    if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
        ckit =False
    else:
        ckit = True 
    return ckit

##
##*****************************************************************************
##
def make_lcurvenet(tcurvenet):#,kind='InterpolatedBspline'):
    """issues a full complement 
    to the tcurvenet.
    -all ofthe duals of the tcurvenet
    
    Ensure that the knot vectors are -nice-
    """
    if verbose:print 'utility function: make_lcurvenet'
    #    kind = {'Bspline':spline.Bspline,
    #            'InterpolatedBspline':spline.interpolatedBspline,
    #            'THBspline':thbspline.rbspline}
    #    assert(isinstance(kind,(spline.Bspline,
    #                      spline.interpolatedBspline,
    #                      thbspline.rbspline))),'Kind is not availible'
    k=tcurvenet[0].k
    nump = tcurvenet[0].nump
    N = tcurvenet[0].n
    
    tvertnet = np.asarray(
            [el.vertices for el in tcurvenet])
    lcurvenet = []
    for i in range(N):
        lcurvenet.append(
            spline.interpolatedBspline(
            tvertnet[:,i], k, nump))
        
    return lcurvenet


##
##*****************************************************************************
##
##
##*****************************************************************************
##
def make_thb(curve):
    """
    parameters
    -----------------
        Bspline curve
    returns:
    -----------------
        THBspline curve
    """
    vts     = curve.vertices
    k       = curve.k
    nump    = curve.nump
    ncurve = rbspline(vts, k, nump)
    if hasattr(curve, 'ivertices'):
        ncurve.ivertices = curve.ivertices
    return ncurve

def make_thb_3D(curve,
                rdict,
                parent,
                k=4,nump=30):
    index = rdict['index']
    vert = rdict['vertices']
    bds = curve.bounds
    nvts = np.zeros((curve.n,3),float)
    for i in [0,1,2]:
        if i == index:
            nvts[:,i] = vert
        else:
            nvts[:,i] = curve.vertices[:,rdict[i]]
    #
    ncurve = rbspline(nvts, k, nump)
    ncurve.bounds = bds
    ncurve.parent = parent
    if parent is not None:
        parent.children = ncurve
    return ncurve
##
##*****************************************************************************
##
##
##*****************************************************************************
##
    
def make_THB_lcurvenet(tcurvenet):
    """issues a full complement 
    to the tcurvenet.
    -all ofthe duals of the tcurvenet
    """
    if verbose:print 'utility function: make_THB_lcurvenet'
    k=tcurvenet[0].k
    nump = tcurvenet[0].nump
    N = tcurvenet[0].n
    
    tvertnet = np.asarray(
            [el.vertices for el in tcurvenet])
    
    lcurvenet = []
    for i in range(N):
        lcurvenet.append(
            spline.interpolatedBspline(
            tvertnet[:,i], k, nump))
        
    nnet = []
    for curve in lcurvenet:
        vtx = curve.vertices
        nnet.append(
            thbspline.rbspline(
                    vtx, curve.k, curve.nump))
        
    lcurvenet = nnet
    
    return lcurvenet

##
##*****************************************************************************
##
##
##
##
##
##   Adaptive Longitudinal Solver
##
##
##   setup (optional)
##
## curve.level  = 0, 1, 2,  3,  4,  5
## curve.n      = 4, 5, 7, 11, 19, 35
##
##*****************************************************************************
def setup_bare_hull_solver(thbcurve,
                           vtinterp,
                           d1=3,d2=4):
    """THB-spline Bare hull in 3 parts
    """
    forsec = [0,1,2,3]  # h1
    midsec = [3,4]      # h2
    aftsec = [4,5,6,7]  # h3
    h1 = h2 = h3 = None
    return

def setup_interpolation_solver(thbcurve,
                               vtinterp,
                               this_kind='equality'):
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
    print '   uopt.setup_interpolation_solver'
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
    FPD.add_E1(kind='LS', weight = 1.)  #reg) #
    FPD.add_E2(kind='LS', weight = .5  )  # reg) #
    FPD.add_E3(kind='LS', weight = .05  )  #  reg) #
    #FPD.add_ArcLengthApprox(kind='LS',weight = big)
    FPD.add_ArcLengthApprox(kind='LS',weight = mig)
    #FPD.add_ArcLengthApprox(kind='LS',weight = smallbig)
    ntinterp = len(vtinterp)
    ukbar = np.linspace(0.,1.,ntinterp)[1:-1]
    nth = 1
    totn = ntinterp
    zb = thbcurve.vertices[0,2]
    for pt,loc in zip(vtinterp[1:-1],ukbar):
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
        #*********************************************
        if nth == 3:
            start_flat = pt
            print ' got an integer for  fwd end of flat portion nth = ',nth
            FPD.add_AngleXoverZConstraint(kind      = LSkind,
                                          location  = loc,
                                          value     = 0.,
                                          weight    = big)
            FPD.add_CurvatureXoverZConstraint(kind  = LSkind,
                                          location  = loc,
                                          value     = 0.,
                                          weight    = big)
            #-----------------------------------------------
            #        if nth == 4:
            #            print ' got an integer for  fwd end of flat portion nth = ',nth
            #            FPD.add_AngleXoverZConstraint(kind='LS',
            #                                          location = loc,
            #                                          value = loc,
            #                                          weight=big)
            #            FPD.add_CurvatureXoverZConstraint(kind='LS',
            #                                        location=loc,
            #                                        value = 0.,
            #                                        weight =big)
##        #*********************************************
##        #*********************************************
        if nth == totn-3-1:
            end_flat = pt
            print ' got an integer for aft end of flat portion nth = ',nth
            FPD.add_AngleXoverZConstraint(kind      = LSkind,
                                          location  = loc,
                                          value     = 0.,
                                          weight    = big)
            FPD.add_CurvatureXoverZConstraint(kind  = LSkind,
                                          location  = loc,
                                          value     = 0.,
                                          weight    = big)
            #-----------------------------------------------
#        if (nth == totn-3-2) and (totn-3-2 > 3+1):
#            print ' got an integer for aft end of flat portion nth = ',nth
#            FPD.add_AngleXoverZConstraint(kind = this_kind,
#                                          location = loc,
#                                          value = 0.,
#                                          weight=big)
#            FPD.add_CurvatureXoverZConstraint(kind=this_kind,
#                                        location=loc,
#                                        value = 0.,
#                                        weight =big)
        #*********************************************
        #*********************************************
        nth +=1
    
    #*********************************************
    #tan at 0.:
    FPD.add_zFixity(index=1,
                    value=zb,
                    track_prolongation=False)
#    FPD.add_zFixity(index=2,
#                    value=zb,
#                    track_prolongation=False)
#-------------------------------------
#    for i in [5]:
#        print 'fwd flat = ',vtinterp[3,2]
#        FPD.add_xFixity(index=i,
#                    value=vtinterp[3,0],
#                    track_prolongation=False)
#        FPD.add_yFixity(index=i,
#                    value=vtinterp[3,1],
#                    track_prolongation=False)
#        FPD.add_zFixity(index=i,
#                    value=vtinterp[3,2],
#                    track_prolongation=False)
#-------------------------------------
#    for i in [5]:
#        FPD.add_xFixity(index=i,
#                    value=vtinterp[4,0])
#        FPD.add_yFixity(index=i,
#                    value=vtinterp[4,1])
#        FPD.add_zFixity(index=i,
#                    value=vtinterp[4,2])
#-------------------------------------
    for i,diff in zip([5,6],[.5,0.]):
        print 'fwd flat = ',vtinterp[3,2]
        FPD.add_xFixity(index=i,
                    value=vtinterp[3,0],
                    track_prolongation=False)
        FPD.add_yFixity(index=i,
                    value=vtinterp[3,1],
                    track_prolongation=False)
        FPD.add_zFixity(index=i,
                    value=vtinterp[3,2]-diff,
                    track_prolongation=False)
#-------------------------------------
#    for i in [6]:
#        print 'fwd flat = ',vtinterp[3,2]
#        FPD.add_xFixity(index=i,
#                    value=vtinterp[4,0],
#                    track_prolongation=False)
#        FPD.add_yFixity(index=i,
#                    value=vtinterp[4,1],
#                    track_prolongation=False)
#        FPD.add_zFixity(index=i,
#                    value=vtinterp[4,2],
#                    track_prolongation=False)
#-------------------------------------
# obviously, you need to make prolongated 
# fixities stay next to each other!
    for i,diff in zip([6,7],[0.,1.]):
        print 'aft flat = ',vtinterp[4,2]
        FPD.add_xFixity(index=i,
                    value=vtinterp[4,0],
                    track_prolongation=False)
        FPD.add_yFixity(index=i,
                    value=vtinterp[4,1],
                    track_prolongation=False)
        FPD.add_zFixity(index=i,
                    value=vtinterp[4,2]+diff,
                    track_prolongation=False)    
#        
        
        ####################
#    for i,diff in zip([4,5,6],[0.,.1,.2]):
#        FPD.add_xFixity(index=i,
#                    value=vtinterp[3,0])
#        FPD.add_yFixity(index=i,
#                    value=vtinterp[3,1])
#        FPD.add_zFixity(index=i,
#                    value=vtinterp[3,2]-diff)
#    for i,diff in zip([7,8,9],[0.,.1,.2]):
#        FPD.add_xFixity(index=i,
#                    value=vtinterp[4,0])
#        FPD.add_yFixity(index=i,
#                    value=vtinterp[4,1])
#        FPD.add_zFixity(index=i,
#                    value=vtinterp[4,2]-diff)
        
    #*********************************************
    #*********************************************
###
#    FPD.add_CurvatureXoverZConstraint(kind      = LSkind,
#                                      location  = 0.5,
#                                      value     = 0.,
#                                      weight    = big)
#    FPD.add_AngleXoverZConstraint(    kind      = LSkind,#this_kind,
#                                      location  = 0.5,
#                                      value     = 0.,
#                                      weight    = big)
##
    
    #*********************************************
    #*********************************************
    #*********************************************
    #*********************************************
#    FPD.add_verticalXoverZAngleConstraint(kind      = this_kind,
#                                          location  = 0.,
#                                          value     = 0.,
#                                          weight    = big)
    #*********************************************
    #*********************************************
    #*********************************************
    #*********************************************
#    FPD.add_CurvatureXoverZConstraint(kind = LSkind,
#                                      location  = 0.,
#                                      value     = 0.,
#                                      weight    = big)
##        
    #*********************************************
    #*********************************************
    zmid = .5*(end_flat[2] - start_flat[2])
    print 'midpoint = ',zmid
#    FPD.add_xPointConstraint(kind       = LSkind, 
#                             value      = start_flat[0], 
#                             location   = 0.5, 
#                             weight     = big ) 
#    FPD.add_yPointConstraint(kind       = LSkind, 
#                             value      = start_flat[1], 
#                             location   = 0.5, 
#                             weight     = big  )
#    FPD.add_zPointConstraint(kind       = LSkind, 
#                             value      = zmid, 
#                             location   = loc, 
#                             weight     = big  )
        
##  #*********************************************
    #*********************************************
    L = Lagrangian(FPD)
    interval_data, small = lagrangian_bounds(L, 
                        interval_data, small, 1.e4)
    Lspline = IntervalLagrangeSpline(thbcurve, L, 
                                     data = interval_data)
    return Lspline














##
##***************************************************************
##   2D Longitudinals-fore and aft specific constraints

def setup_intrp_solve_fwd_DWL(FPD,
                             original_curve,
                             live_vertices,
                             this_kind='equality',
                             scalar=1.):
    """2D DWL
    longitudinal curve generation
    for interior longitudinals
    
    
    -interpolate vertices (done in main setup)
    -forward and aft flats.
    """
    fullcurve = original_curve
    #
    area = fullcurve.area
    Xc = fullcurve.Xc
    #
    ti = fullcurve.compute_tangent(0.,vertices=live_vertices)
    te = fullcurve.compute_tangent(1.,vertices=live_vertices)
    ci = fullcurve.compute_curvature(0.,vertices=live_vertices)
    ce = fullcurve.compute_curvature(1.,vertices=live_vertices)
    #-----------------------------------------------
    # beggining curvature, method 1
#    FPD.add_verticalAngleConstraint(
#                                  kind      = this_kind,
#                                  location  = 0.,
#                                  value     = 0.,
#                                  weight    = scalar)
#    FPD.add_CurvatureConstraint(kind  = this_kind,
#                                  location  = 0.,
#                                  value     = 0.,
#                                  weight    = scalar)
#    #-----------------------------------------------
#    # ending curvature, method 1
#    FPD.add_AngleConstraint(
#                                  kind      = this_kind,
#                                  location  = 1.,
#                                  value     = 0.,
#                                  weight    = scalar)
#    FPD.add_CurvatureConstraint(
#                                  kind  = this_kind,
#                                  location    = 1.,
#                                  value       = 0.,
#                                  weight      = scalar)
    #-----------------------------------------------
    # beggining curvature, method 2 (also works)
    FPD.add_xFixity(index=1,
                value = live_vertices[0,0],
                track_prolongation=False)
    FPD.add_xFixity(index=2,
                value = live_vertices[0,0],
                track_prolongation=False)
    #-----------------------------------------------
    # end curvature, method 2  (also works)
#    FPD.add_xFixity(index=-2,
#                value = live_vertices[-1,0],
#                track_prolongation=False)
#    FPD.add_xFixity(index=-3,
#                value = live_vertices[-1,0],
#                track_prolongation=False)
    FPD.add_yFixity(index=-2,
                value = live_vertices[-1,1],
                track_prolongation=False)
    FPD.add_yFixity(index=-3,
                value = live_vertices[-1,1],
                track_prolongation=False)
#    FPD.add_yFixity(index=-4,
#                value = live_vertices[-1,1],
#                track_prolongation=False)
    #-----------------------------------------------
    # need to work on 
    # the original DWL maker
    # to ensure equality 
    # can be maintained while
    # gettting a flat transition to midship aft
    # oh wait, these are 2D.
    FPD.add_AreaConstraint(
                                  kind      = this_kind,
                                  value     = area,
                                  weight    = scalar)
    FPD.add_XcConstraint(
                                  kind      = this_kind,
                                  value     = Xc,
                                  weight    = scalar)
    return FPD


def setup_intrp_solve_aft_DWL(FPD,
                             original_curve,
                             live_vertices,
                             this_kind='equality',
                             scalar=1.):
    """2D DWL
    longitudinal curve generation
    for interior longitudinals
    
    
    -interpolate vertices (done in main setup)
    -forward and aft flats.
    """
    fullcurve = original_curve
    #
    area = fullcurve.area
    Xc = fullcurve.Xc
    #
    ti = fullcurve.compute_tangent(0.,vertices=live_vertices)
    te = fullcurve.compute_tangent(1.,vertices=live_vertices)
    ci = fullcurve.compute_curvature(0.,vertices=live_vertices)
    ce = fullcurve.compute_curvature(1.,vertices=live_vertices)
    #
    FPD.add_AngleConstraint(
                                  kind      = this_kind,
                                  location  = 0.,
                                  value     = 0.,
                                  weight    = scalar)
    FPD.add_CurvatureConstraint(kind  = this_kind,
                                  location  = 0.,
                                  value     = 0.,
                                  weight    = scalar)
    #-----------------------------------------------
    FPD.add_AngleConstraint(
                                  kind      = this_kind,
                                  location  = 1.,
                                  value     = 0.,
                                  weight    = scalar)
    FPD.add_CurvatureConstraint(
                                  kind  = this_kind,
                                  location    = 1.,
                                  value       = 0.,
                                  weight      = scalar)
    #-----------------------------------------------
#    FPD.add_yFixity(index=1,
#                value = live_vertices[0,1],
#                track_prolongation=False)
#    FPD.add_yFixity(index=2,
#                value = live_vertices[0,1],
#                track_prolongation=False)
#    #-----------------------------------------------
#    FPD.add_yFixity(index=-1,
#                value = live_vertices[-1,1],
#                track_prolongation=False)
#    FPD.add_yFixity(index=-2,
#                value = live_vertices[-1,1],
#                track_prolongation=False)
    #-----------------------------------------------
    FPD.add_AreaConstraint(
                                  kind      = this_kind,
                                  value     = area,
                                  weight    = scalar)
    return FPD


def setup_intrp_solve_fwd_CPK(FPD,
                              original_curve,
                              live_vertices,
                              this_kind='equality',
                              scalar=1.):
    """2D DWL
    longitudinal curve generation
    for interior longitudinals
    
    
    -interpolate vertices (done in main setup)
    -forward and aft flats.
    """
    print 'setup_intrp_solve_fwd_CPK'
    fullcurve = original_curve
    #
    area = fullcurve.area
    Xc = fullcurve.Xc
    #
    ti = fullcurve.compute_tangent(0.,vertices=live_vertices)
    te = fullcurve.compute_tangent(1.,vertices=live_vertices)
    ci = fullcurve.compute_curvature(0.,vertices=live_vertices)
    ce = fullcurve.compute_curvature(1.,vertices=live_vertices)
    #
    print 'area ',area
    print 'Xc ',Xc
    
    # this is not the same as the original CPK!  
    #it does not 'go all the way up' to the DWL so 
    # the old stem requirement is allready fulfilled by the transverse
    #-----------------------------------------------
#    FPD.add_AngleConstraint(
#                                  kind      = this_kind,
#                                  location  = 0.,
#                                  value     = ti,
#                                  weight    = scalar)
#    FPD.add_CurvatureConstraint(
#                                  kind      = this_kind,
#                                  location  = 0.,
#                                  value     = ci,
#                                  weight    = scalar)
    #-----------------------------------------------
    #no need for curveture 0 at s=.04 due to this not bing the full cLP anymore
    #-----------------------------------------------
    FPD.add_AngleConstraint(
                                  kind      = this_kind,
                                  location  = 1.,
                                  value     = 0.,
                                  weight    = scalar)
    FPD.add_CurvatureConstraint(
                                  kind  = this_kind,
                                  location    = 1.,
                                  value       = 0.,
                                  weight      = scalar)
    #-----------------------------------------------
    FPD.add_yFixity(index=-2,
                value = FPD.curve.vertices[-1,1],
                track_prolongation=False)
    FPD.add_yFixity(index=-3,
                value = FPD.curve.vertices[-1,1],
                track_prolongation=False)
    FPD.add_yFixity(index=-4,
                value = FPD.curve.vertices[-1,1],
                track_prolongation=False)
    FPD.add_AreaConstraint(
                                  kind      = this_kind,
                                  value     = area,
                                  weight    = scalar)
#    FPD.add_XcConstraint(
#                                  kind      = this_kind,
#                                  value     = Xc,
#                                  weight    = scalar)
    return FPD


def setup_intrp_solve_aft_CPK(FPD,
                              original_curve,
                              live_vertices,
                              this_kind='equality',
                              scalar=1.):
    """2D DWL
    longitudinal curve generation
    for interior longitudinals
    
    
    -interpolate vertices (done in main setup)
    -forward and aft flats.
    """
    fullcurve = original_curve
    #
    area = fullcurve.area
    Xc = fullcurve.Xc
    #
    ti = fullcurve.compute_tangent(0.,vertices=live_vertices)
    te = fullcurve.compute_tangent(1.,vertices=live_vertices)
    ci = fullcurve.compute_curvature(0.,vertices=live_vertices)
    ce = fullcurve.compute_curvature(1.,vertices=live_vertices)
    #
#    FPD.add_AngleConstraint(
#                                  kind      = this_kind,
#                                  location  = 0.,
#                                  value     = 0.,
#                                  weight    = scalar)
#    FPD.add_CurvatureConstraint(kind  = this_kind,
#                                  location  = 0.,
#                                  value     = 0.,
#                                  weight    = scalar)
    #-----------------------------------------------
    #FPD.add_AngleConstraint( ABSOLUTELY NOT!!!location  = 1.)
    
#    FPD.add_CurvatureConstraint(
#                                  kind  = this_kind,
#                                  location    = 1.,
#                                  value       = 0.,
#                                  weight      = scalar)
    #-----------------------------------------------
    
    FPD.add_yFixity(index=1,
                value = FPD.curve.vertices[0,1],
                track_prolongation=False)
    FPD.add_yFixity(index=2,
                value = FPD.curve.vertices[0,1],
                track_prolongation=False)
    FPD.add_AreaConstraint(
                                  kind      = this_kind,
                                  value     = area,
                                  weight    = scalar)
    return FPD

##
##***************************************************************
##***************************************************************
##   3D Longitudinals-fore and aft specific constraints


def setup_intrp_solve_fwd(FPD,
                          this_kind='LS',
                          scalar=1.):
    """3D longitudinals
    FWD section longitudinal curve generation
    for interior longitudinals
    
    
    -interpolate vertices (done in main setup)
    -forward and aft flats.
    
    this_kind='equality',
    this_kind='LS',
    """
    print 'setup_intrp_solve_fwd'
    vertices = FPD.curve.vertices
    #-----------------------------------------------
#    FPD.add_verticalXoverZAngleConstraint(
#                                  kind      = this_kind,
#                                  location  = 0.,
#                                  value     = 0.,
#                                  weight    = scalar)
#    FPD.add_CurvatureXoverZConstraint(kind  = this_kind,
#                                  location  = 0.,
#                                  value     = 0.,
#                                  weight    = scalar)
    #-----------------------------------------------
#    FPD.add_AngleXoverZConstraint(kind      = this_kind,
#                                  location  = 1.,
#                                  value     = 0.,
#                                  weight    = scalar)
#    FPD.add_CurvatureXoverZConstraint(kind  = this_kind,
#                                  location    = 1.,
#                                  value       = 0.,
#                                  weight      = scalar)
    #-----------------------------------------------
    # replace fwd part of fwd continuity condtions
    FPD.add_zFixity(index=1,
                value = vertices[0,2],
                track_prolongation=False)
    FPD.add_zFixity(index=2,
                value = vertices[0,2],
                track_prolongation=False)
    #-----------------------------------------------
    # replace rear X (breatdth) part of fwd continuity condtions
    FPD.add_xFixity(index=-2,
                value = vertices[-1,0],
                track_prolongation=False)
#    FPD.add_xFixity(index=-3,
#                value = vertices[-1,0],
#                track_prolongation=False)
    #-----------------------------------------------
    # replace rear Y (depth) part of fwd continuity condtions
    FPD.add_yFixity(index=-2,
                value = vertices[-1,1],
                track_prolongation=False)
#    FPD.add_yFixity(index=-3,
#                value = vertices[-1,1],
#                track_prolongation=False)
    #-----------------------------------------------
    #idea: push the final 1 or two z vertices fwd a bit to make 
    # things less discontinuous
#    FPD.add_zFixity(index=-1,
#                value = vertices[-1,2],
#                track_prolongation=False)
#    FPD.add_zFixity(index=-2,
#                value = vertices[-1,2],
#                track_prolongation=False)
    #-----------------------------------------------
    return FPD

def setup_intrp_solve_mid(FPD,
                          this_kind='equality',
                          scalar=1.):
    """
    FWD section longitudinal curve generation
    for interior longitudinals
    
    
    -interpolate vertices (done in main setup)
    -forward and aft flats.
    """
    print 'setup_intrp_solve_mid'
    #-----------------------------------------------
    return FPD

def setup_intrp_solve_aft(FPD,
                          this_kind='LS',
                          scalar=1.):
    """3D longitudinals
    FWD section longitudinal curve generation
    for interior longitudinals
    
    
    -interpolate vertices (done in main setup)
    -forward and aft flats.
    
    this_kind='equality',
    this_kind='LS',
    """
    print 'setup_intrp_solve_aft'
    vertices = FPD.curve.vertices
    #-----------------------------------------------
#    FPD.add_AngleXoverZConstraint(
#                                  kind      = this_kind,
#                                  location  = 0.,
#                                  value     = 0.,
#                                  weight    = scalar)
#    FPD.add_CurvatureXoverZConstraint(
#                                  kind      = this_kind,
#                                  location  = 0.,
#                                  value     = 0.,
#                                  weight    = scalar)
#    #-----------------------------------------------
#    FPD.add_AngleXoverZConstraint(kind      = this_kind,
#                                  location  = 1.,
#                                  value     = 0.,
#                                  weight    = scalar)
#    FPD.add_CurvatureXoverZConstraint(kind  = this_kind,
#                                  location    = 1.,
#                                  value       = 0.,
#                                  weight      = scalar)
    #-----------------------------------------------
    FPD.add_xFixity(index=1,
                value = vertices[0,0],
                track_prolongation=False)
#    FPD.add_xFixity(index=2,
#                value = vertices[0,0],
#                track_prolongation=False)
    #-----------------------------------------------
    FPD.add_yFixity(index=1,
                value = vertices[0,1],
                track_prolongation=False)
#    FPD.add_yFixity(index=2,
#                value = vertices[0,1],
#                track_prolongation=False)
    #-----------------------------------------------
#    FPD.add_xFixity(index=-2,
#                value = vertices[-1,0],
#                track_prolongation=False)
#    FPD.add_xFixity(index=-3,
#                value = vertices[-1,0],
#                track_prolongation=False)
    #-----------------------------------------------
    return FPD



def setup_interpolation_solver_parsi(thbcurve,
                               vtinterp,
                               fwd=False,
                               mid=False,
                               aft=False,
                               this_kind='equality'):
    """setup 3D THBcurve Lspline vertex interpolation solver
    interior longitudinals
    
    Lspline solver 3 part
    
    
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
    print '   uopt.setup_interpolation_solver'
    interval_data, small = interval_bounds(thbcurve)
    FPD = FormParameterDict(thbcurve)
    
    LSkind = 'LS'
    
    #big = 1.
    big  = 100000. #   1.  #
    mig =  50000.
    #mig =  10000.
    smallbig = 1000.
    small = .1
    reg = 1.
    FPD.add_E1(kind='LS', weight = 1.)  #reg) #
    FPD.add_E2(kind='LS', weight = .5  )  # reg) #
    FPD.add_E3(kind='LS', weight = .05  )  #  reg) #
    FPD.add_ArcLengthApprox(kind='LS',weight = big)
    ntinterp = len(vtinterp)
    #
    #**************************************************
    #
    """briefly misguided here:
        Do not treat the vertices to be interpolated
        as if they had the endpoints amond them.
        They do not.
        
    assert(vtinterp do not contain the end points!)
    for pt,loc in zip(vtinterp[1:-1],ukbar):
    """
    ukbar = np.linspace(0.,1.,ntinterp+2)[1:-1]
    targets = vtinterp
    nth = 1
    totn = ntinterp
    zb = thbcurve.vertices[0,2]
    for pt,loc in zip(targets,ukbar):
        #*********************************************
        #print '\nfixing a vertex\n'
        #print 'vertex = ',pt
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
    if fwd:
        FPD = setup_intrp_solve_fwd(FPD)
    elif aft:
        FPD = setup_intrp_solve_aft(FPD)

    #*********************************************
    L = Lagrangian(FPD)
    interval_data, small = lagrangian_bounds(L, 
                        interval_data, small, 1.e4)
    Lspline = IntervalLagrangeSpline(thbcurve, L, 
                                     data = interval_data)
    return Lspline
##
##*****************************************************************************
##
##
##
##
##
##   Adaptive Longitudinal Solver
##
##
##   solver
##
## curve.level  = 0, 1, 2,  3,  4,  5
## curve.n      = 4, 5, 7, 11, 19, 35
##
##*****************************************************************************
def THBsolver(thbcurve=None,
              vtinterp=None,
              Lspline=None,
              refinemax=1,
              maxlevel=4,
              normalize=True):
    """Lspline solver (adaptive) for THB splines
    
    intention : demonstrate an adaptive solver
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
#    
#    print '   uopt.THBsolver'
#    print '\n inputs:'
#    print '   refinemax = ',refinemax, ' no longer used'
#    print '   maxlevel  = ',maxlevel
#    print '   normalize = ',normalize
    num_refined = 0
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
    if vtinterp is not None and Lspline is None:
        if thbcurve.n<11:
            while nthb_vs_nvi(thbcurve,vtinterp) and thbcurve.n<11:
                print 'pre-refinement'
                thbcurve = thbcurve.dyadic_refinement()
        #else: we chance that it will refine in time!
        #**************************************************
        thbcurve.parent = None
        #
        print '   setup'
        Lspline = setup_interpolation_solver(thbcurve,
                                             vtinterp)
        #******************************************
    #print 'switch equality to LS'
    Lspline.switch_all_of_one_Lagrangian_kind(
                                        ini_kind='equality',
                                        fini_kind='LS')
    #******************************************
    #
    print '   initial-optimize'
    Lspline.tol = 1.e-10
    Lspline.optimize_thb(stop=5,
                         ceiling = 1.e-18,
                         relax_grad_bounds=False)
    print Lspline.refine_by_these_grad_indices
    #******************************************
    #print 'restore LS to equality'
    Lspline.restore_Lagrangian()
    print 'again = ',Lspline.refine_by_these_grad_indices
    #------------------
    #normalize
    if normalize: Lspline.normalize()
    #    
    #******************************************
    # manually call refine:
    #Lspline.refine_by_force(
    #        boundslst = [ia(0.,.2),ia(.8,1.)]) 
    if Lspline.curve.level<maxlevel:
        Lspline.refine()
        this = Lspline.get_finest_Lspline()
        if this is not Lspline:
            num_refined +=1
    else:
        this = Lspline
#    #******************************************
    #print 'switch equality to LS'
    this.switch_all_of_one_Lagrangian_kind(
                                    ini_kind='equality',
                                    fini_kind='LS')
    #******************************************
    #
    print '   secondary smooth'
    this.tol = 1.e-5
    this.optimize_thb(stop=5,
                         ceiling = 1.e-10,
                         relax_grad_bounds=False)
    #******************************************
    cmax = maxlevel+1-this.curve.level #max num times to iterate
    print 'maxlevel = ',maxlevel
    print 'this.curve.level = ',this.curve.level
    print 'cmax while iter num = ',cmax
    count =0
    while this is not None and count<cmax:
        print '   optimize'
        #******************************************
        # retore equality constraints
        print 'restore LS to equality'
        this.restore_Lagrangian()
        #------------------
        #normalize
        if normalize: this.normalize()
        #******************************************
        print '   optimize eq constraints'
        this.verbose = True
        """(1) save and 
        restoring the vertices allows 
        THB hierarchical opti to work (1 of 2)"""
        save_vertices = copy.copy(this.curve.vertices)
        #*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        this.tol = 1.e-7
        this.verbose=True
        this.optimize_thb(stop=15,
                          ceiling = 1.e-10,
                          relax_grad_bounds=False)
        this.tol = 1.e-7
        #**************************************************
        ckit = ckLspline(this)
        if ckit and this.conv<this.tol:
            print 'converged: ',this.conv
            setconv = True
            break
        elif not ckit:
            print 'restore vertices'
            """(2) makes THB hierarchical opti work! (2 of 2)"""
            this.curve.vertices = save_vertices
            #*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        #**************************************************
        # manually call refine:
        #******************************************
        #print 'restore LS to equality'
        #this.restore_Lagrangian() #not necessary to restore as we have already restored
        print 'again = ',this.refine_by_these_grad_indices
        if this.curve.level<maxlevel:
            #if num_refined < refinemax:
            print 'Not converged, refining again...'
            this.refine()
            newlevel = this.get_finest_Lspline()
            if newlevel is not this:
                num_refined +=1
            this = newlevel
#        #******************************************
#            this = Lspline.get_finest_Lspline()
#            this.tol = 1.e-10
#            this.verbose = True
        #******************************************
        #do some smoothing:
        print 'switch equality to LS'
        this.switch_all_of_one_Lagrangian_kind(
                ini_kind='equality',
                fini_kind='LS')
        #******************************************
        print '   optimize'
        this.optimize_thb(stop=7,
                          ceiling = 1.e-7,
                          relax_grad_bounds=False)
        setconv = False
#        #******************************************
        count +=1

    
    #******************************************
    # retore equality constraints
    print 'restore LS to equality'
    this.restore_Lagrangian()
    #******************************************
    #print '   optimize'
    #Lspline.optimize_thb(stop=5,
    #                     relax_grad_bounds=True)
    print 'Done with this longitudinal'
    if setconv:
        print 'convergence is true'
    else:
        print 'convergence is false'
    print 'refined {} times'.format(num_refined)
    print 'while count = {}  at end'.format(count)
    return this.get_coarsest_Lspline()


##
##*****************************************************************************
##
##
##
##
##
##   Adaptive Longitudinal Solver
##
##
##   Driving Function
##
## curve.level  = 0, 1, 2,  3,  4,  5
## curve.n      = 4, 5, 7, 11, 19, 35
##
##*****************************************************************************
def issue_longitudinals(tcurvenet,
                        k=4,nump=30,
                       refinemax=1,
                       maxlevel=4,
                       normalize = False,
                       refine1st = False):
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
    for i in range(1,nnew+1):
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
        # speed things up by going ahead and refining here?
        if refine1st:
            thbcurve = thbcurve.dyadic_refinement()#5
            thbcurve = thbcurve.dyadic_refinement()#7
        #
        print 'new curve ', i
        lcurvenet.append(THBsolver(thbcurve,
                                   vtinterp,
                                   refinemax=refinemax,
                                   maxlevel=maxlevel,
                                   normalize=normalize))
    return lcurvenet



def thbmatch2D(fullcurve,
               matchpts,
               locations = None,
               n=4,
               #d0=2,d1=0,rep=1,
               k=4,nump=30,
               refinemax=1,
               maxlevel=4,
               normalize=False,
               isdwl=False,
               iscpk=False,
               yindex = 0,
               Longindex=2,
               fwd=False,
               aft=False):
    """
    
    dev
    ----------
    
    
    import FormParameter as fmp
    FormParameterDict = fmp.FormParameterDict
    
    
    n            = 4
    d0           = 2
    d1=0
    rep=1
    k=4
    nump=30
    refinemax=1
    maxlevel=4
    normalize=False
    isdwl=False
    iscpk=False
    yindex = 0
    Longindex=2
    fwd      = False
    aft      = False
    
    
    Longindex=2
    n=7
    k=4
    nump=30
    d0=2
    d1=0
    rep=1
    #
    yindex = 0, # x in 3d acts in the y direction here in 2D
    Longindex=2 # y in 3D acts in the longitudinal (x) direction here in 2D
    
    
    fullcurve    = fwd_cpk_a
    matchpts     = tvertnet[1:3,i]
    
    yindex = 1
    Longindex=2
    
    """
    #d0=Longindex
    #d1=yindex
    #rep=1
    this_kind = 'equality'
    big = 1.e5
    #
    #n = fullcurve.n
    cv0 = fullcurve.vertices[:,Longindex]
    cv1 = fullcurve.vertices[:,yindex]
    # not something that is THB
    # just a curve from a splitting:
#    #---------------------------
    cv = np.zeros((fullcurve.n,2),float)
    cv[:,0] = cv0
    cv[:,1] = cv1
    #**************************
    # convert matchpoints to 2D
    npts = np.zeros((len(matchpts),2),float)
    print 'npts = '
    print npts
    print 'matchpts = '
    print matchpts
    print 'Longindex = ',Longindex
    npts[:,0] = matchpts[:,Longindex]
    npts[:,1] = matchpts[:,yindex]
    matchpts = npts
    print 'matchpts = '
    print matchpts
    #
    yindex = 1
    Longindex=0
    #*********************************************
    #
    lvts = linear_vertices(cv[0],
                           cv[-1],
                           n)
    #
    #*********************************************
    curve = rbspline(lvts,k,nump)
    #
    curve = curve.dyadic_refinement()
    curve = curve.dyadic_refinement()
    curve.parent = None
    #curve.plotcurvehierarchy()
    curve.basis_matrix()
    curve.compute_area()
    curve.compute_moments()
    curve.computeXc()
    #
    # ---- fullcurve - original curve - FPDs
    #
    fullcurve.basis_matrix()
    fullcurve.compute_area(vertices = cv)
    fullcurve.compute_moments(vertices = cv)
    fullcurve.computeXc()
    #
    #
    #
    #**************************************************************************
    # 
    interval_data, small = interval_bounds(curve)
    FPD = FormParameterDict(curve)
    #
    #**************************************************************************
    # normalize
    curve.fairness()
    curve.pts_M_pts()
    curve.compute_arclength()
    curve.compute_arclengthApprox()
    E1_0 = curve.E1
    E2_0 = curve.E2
    E3_0 = curve.E3
    S_0  = curve.ALapprox
    #
    #**************************************************************************
    # Form Parameters
    wt = 1.
    #
    #**************************************************************************
    # positional parameters
    if locations is None:
        for pt in matchpts:
            #*********************************************
            loc = curve.FindPoint(pt[Longindex],Longindex)[0]
            #*********************************************
            if fwd and isdwl:
                print '*********************************************'
                print 'fwd DWL'
                print '\n constraints:'
                print 'location',loc
                print 'x,y =  ',pt[Longindex],pt[yindex]
                print 'kind = ', this_kind
                print '*********************************************'
            #*********************************************
            #
            FPD.add_xPointConstraint(kind= this_kind, 
                                     value=pt[Longindex], 
                                     location=loc, 
                                     weight = big ) 
            FPD.add_yPointConstraint(kind= this_kind, 
                                     value=pt[yindex], 
                                     location=loc, 
                                     weight = big  )
            #*********************************************
    else:
        for pt, loc in zip(matchpts,locations):
            #*********************************************
            if fwd and isdwl:
                print '*********************************************'
                print 'fwd DWL'
                print '\n constraints:'
                print 'location',loc
                print 'x,y =  ',pt[Longindex],pt[yindex]
                print 'kind = ', this_kind
                print '*********************************************'
            #*********************************************
            #
            print 'loc = ',loc
            FPD.add_xPointConstraint(kind= this_kind, 
                                     value=pt[Longindex], 
                                     location=loc, 
                                     weight = big ) 
            FPD.add_yPointConstraint(kind= this_kind, 
                                     value=pt[yindex], 
                                     location=loc, 
                                     weight = big  )
            #*********************************************
    if fwd and isdwl:
        print 'fwd DWL'
        FPD = setup_intrp_solve_fwd_DWL(FPD,
                                        original_curve=fullcurve,
                                        live_vertices = cv)
    elif fwd and iscpk:
        print 'fwd cpkeel profile'
        FPD = setup_intrp_solve_fwd_CPK(FPD,
                                        original_curve=fullcurve,
                                        live_vertices = cv)
    elif aft and isdwl:
        print 'aft DWL'
        FPD = setup_intrp_solve_aft_DWL(FPD,
                                        original_curve=fullcurve,
                                        live_vertices = cv)
    elif aft and iscpk:
        print 'aft cpkeel profile'
        FPD = setup_intrp_solve_aft_CPK(FPD,
                                        original_curve=fullcurve,
                                        live_vertices = cv)
        
    #print 'constraints:'
    #for el in FPD.FormParam:
    #    fp = FPD.FormParam[el]
    #    print 'type = ',fp.type
    #    print 'pass val = ',fp.pass_value
    #
    #**************************************************************************
    # Fairness
#    FPD.add_E1(kind='LS', weight = 1./E1_0)
#    FPD.add_E2(kind='LS', weight = 1./E2_0)
#    FPD.add_E3(kind='LS', weight = 1./E3_0)
#    FPD.add_ArcLengthApprox(kind='LS', weight = 1./S_0)
    FPD.add_E1(kind='LS', weight = .5)
    FPD.add_E2(kind='LS', weight = .1)
    FPD.add_E3(kind='LS', weight = .01)
    FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
    #
    #**************************************************************************
    #
    interval_data, small = interval_bounds(curve)
    L = Lagrangian(FPD)
    interval_data, small = lagrangian_bounds(L, 
                        interval_data, small, 1.e4)
    Lspline = IntervalLagrangeSpline(curve, L, 
                                     data = interval_data)
    #
    #**************************************************************************
    #
    Lspline = THBsolver(curve,
                        Lspline = Lspline,
                        refinemax=refinemax,
                        maxlevel=maxlevel,
                        normalize=normalize)
    return Lspline


def lifting_2DDWL(verts, elevation,k,nump):
        """
        Parameters
        --------------------
            verts       : 2D control vertices
            elevation   : 1 dimension will have a constant value; this is it.
            k           : curve.k
            nump        : curve.nump
        
        Returns
        --------------------
            lifted curve (3D)
        
        
        notes
        --------------------
            *automatically lifts only one way: using xy_to_zyx_elevate_noflp
            *Knot removal for reparameterization
            is currently more reliable in 2D
            *We use it to reparameterize now
            *ergo we want to abstract out the CProfile lifting
        """
        vnew  = frame.xy_to_zxy_elevate(verts,
                                        elevation=elevation)
        return spline.Bspline(vnew, k, nump)


def lifting_2DCProfile(verts, elevation,k,nump):
        """
        Parameters
        --------------------
            verts       : 2D control vertices
            elevation   : 1 dimension will have a constant value; this is it.
            k           : curve.k
            nump        : curve.nump
        
        Returns
        --------------------
            lifted curve (3D)
        
        
        notes
        --------------------
            *automatically lifts only one way: using xy_to_zyx_elevate_noflp
            *Knot removal for reparameterization
            is currently more reliable in 2D
            *We use it to reparameterize now
            *ergo we want to abstract out the CProfile lifting
        """
        vnew  = frame.xy_to_zyx_elevate_noflp(verts,
                                            elevation=elevation)
        return spline.Bspline(vnew, k, nump)

def issue_split_longitudinals(tcurvenet,
                              dwl,cpk,
                              k=4,nump=30,
                              maxlevel=4,
                              normalize=False,
                              min_fwd_work=True,
                              tcurve_flat_ini_index=3,
                              tcurve_flat_fin_index=4,
                              Longindex=2):
    """THB Optimize longitudinal curves
    
    #
    tcurvenet - network of tranverse hull curves
    DWL - hull DWL          curve (3D)
    CPK - hull cLProfile    curve (3D)
    # why 3d?  You are splitting and solving in 3D
    
    
    dev
    ----------
    from utility_optimization import *
    
    cpk=original_lcurves[0]
    dwl=original_lcurves[-1]
    
    maxlevel=4
    
    tcurve_flat_ini_index=3
    tcurve_flat_fin_index=4
    
    
    #cpk
    Longindex   = 2
    yindex      = 1
    
    #dwl
    Longindex   = 2
    yindex      = 0
    
    long_ini = tcurvenet[tcurve_flat_ini_index].vertices[0,Longindex]
    long_fin = tcurvenet[tcurve_flat_fin_index].vertices[0,Longindex]
    
    """
    #
    tvertnet = make_vertex_matrix(tcurvenet)
    #
    lcurvenet = [] #all longitudinal curves
    lcurvenet_fwd = [] 
    lcurvenet_mid = [] 
    lcurvenet_aft = [] 
    #
    curve = tcurvenet[0]
    N = curve.n
    nnew = N-2 #just optimize the interior THB curves?
    #
    #**************************************************************************
    # first split the boundary curves into three parts
    #-------------------------------------------------------------
    # longi location of the splitting (real space)
    long_ini = tcurvenet[tcurve_flat_ini_index].vertices[0,Longindex]
    long_fin = tcurvenet[tcurve_flat_fin_index].vertices[0,Longindex]
    # dwl top
    long_ini_dwl_top = tcurvenet[tcurve_flat_ini_index].vertices[-1,Longindex]
    long_fin_dwl_top = tcurvenet[tcurve_flat_fin_index].vertices[-1,Longindex]
    #-------------------------------------------------------------
    # longi location of the splitting (parametric space)
    dwl_split_ini_s = dwl.FindPoint(long_ini_dwl_top,
                                    Longindex)
    dwl_split_fin_s = dwl.FindPoint(long_fin_dwl_top,
                                    Longindex)
    #
    cpk_split_ini_s = cpk.FindPoint(long_ini,
                                    Longindex)
    cpk_split_fin_s = cpk.FindPoint(long_fin,
                                    Longindex)
    #-------------------------------------------------------------
    # start split curves
    fwd_dwl_a, fwd_dwl_b = dwl.CurveSplit(dwl_split_ini_s[0])
    fwd_cpk_a, fwd_cpk_b = cpk.CurveSplit(cpk_split_ini_s[0])
    #-------------------------------------------------------------
    # end split curves
    aft_dwl_a, aft_dwl_b =  dwl.CurveSplit(dwl_split_fin_s[0])
    aft_cpk_a, aft_cpk_b =  cpk.CurveSplit(cpk_split_fin_s[0])
    #-------------------------------------------------------------
    # flat curves
    dwlflatvert = linear_vertices(fwd_dwl_b.vertices[0],
                                  aft_dwl_b.vertices[0],
                                  k)
    cpkflatvert = linear_vertices(fwd_cpk_b.vertices[0],
                                  aft_cpk_b.vertices[0],
                                  k)
    #
    dwl_flat = rbspline(dwlflatvert,k,nump) #n=4
    cpk_flat = rbspline(cpkflatvert,k,nump) #n=4
    assert(dwl_flat.n == 4),'dwl_flat.n not equal to 4'
    assert(cpk_flat.n == 4),'cpk_flat.n not equal to 4'
    #
    #**************************************************************************
    # Now issue the three part longitudinals:  
    #
    #**************************************************************************
    # FWD
    # fwd cLProfile
    if min_fwd_work:
        fwdlvl = 3
    else:
        fwdlvl = maxlevel
    i=0
    print 'fwd curve {}'.format(i)
    cpkLspline = thbmatch2D(fwd_cpk_a,
                            tvertnet[1:3,i],
                            maxlevel    = fwdlvl,
                            yindex      = 1,
                            Longindex   = 2,
                            iscpk       = True,
                            fwd         = True)
    rstrv = fwd_cpk_a.vertices[0,0]
    cpkLspline.restore = {'index':0,
                          'vertices':rstrv, # centerline x
                          2:0,
                          1:1}
    Lspline = cpkLspline.get_coarsest_Lspline()
    Lspline.restore = cpkLspline.restore
    #append the fwd cLProfile Lspline curve
    lcurvenet_fwd.append(
            fix_Lspline_23D(Lspline) )
    #  
    #--------------------------
    # fwd DWL
    #
    i=6
    print 'fwd curve {}'.format(i)
    loc1x = fwd_dwl_a.FindPoint(
                            tvertnet[1,i][0],0)
    loc1y = fwd_dwl_a.FindPoint(
                            tvertnet[1,i][1],1)
    loc2x = fwd_dwl_a.FindPoint(
                            tvertnet[2,i][0],0)
    loc2y = fwd_dwl_a.FindPoint(
                            tvertnet[2,i][1],1)
    locs = [loc1x[0],loc2x[0]]
    dwlLspline = thbmatch2D(fwd_dwl_a,
                            tvertnet[1:3,i],
                            locations = locs,
                            maxlevel=fwdlvl,
                            isdwl       = True,
                            fwd         = True)
    rstrv = fwd_dwl_a.vertices[0,1]
    dwlLspline.restore = {'index':1,
                          'vertices':rstrv, #dwl height y
                          0:1,
                          2:0}  
    Lspline_ = dwlLspline.get_coarsest_Lspline()
    Lspline_.restore = dwlLspline.restore
    dwlLspline = Lspline_
    #
    #******************************************
    # 3D solving interior longitudinals
    for i in range(1,nnew+1):
        print 'fwd curve {}'.format(i)
        #**************************************
        lvts = linear_vertices(tvertnet[:4,i][0],
                               tvertnet[:4,i][-1],4)
        thbcurve = rbspline(lvts,k,nump)
        thbcurve = thbcurve.dyadic_refinement()
        thbcurve = thbcurve.dyadic_refinement()
        thbcurve.parent = None
        Lspline = setup_interpolation_solver_parsi(thbcurve,
                                                   tvertnet[1:3,i],
                                                   fwd=True)
        Lspline = THBsolver(thbcurve,
                               tvertnet[1:3,i],
                               Lspline  = Lspline,
                               maxlevel = fwdlvl,
                               normalize = normalize)
        lcurvenet_fwd.append(Lspline.curve)
    #******************************************
    #append the fwd DWL Lspline curve
    lcurvenet_fwd.append(
            fix_Lspline_23D(dwlLspline) )
    #
    #**************************************************************************
    # MID
    #---------------------
    i=0
    lvts = linear_vertices(tvertnet[3:5,i][0],
                           tvertnet[3:5,i][-1],4)
    thbcurve = rbspline(lvts,k,nump)
    lcurvenet_mid.append(thbcurve)
    #---------------------
    # i: 1-5
    for i in range(1,nnew+1):
        lvts = linear_vertices(tvertnet[3:5,i][0],
                               tvertnet[3:5,i][-1],4)
        thbcurve = rbspline(lvts,k,nump)
        lcurvenet_mid.append(thbcurve)
    #-----------------------
    i=6
    lvts = linear_vertices(tvertnet[3:5,i][0],
                           tvertnet[3:5,i][-1],4)
    thbcurve = rbspline(lvts,k,nump)
    lcurvenet_mid.append(thbcurve)
    #
    #**************************************************************************
    # AFT
    # aft cLProfile
    i=0
    print 'aft curve {}'.format(i)
    cpkLspline = thbmatch2D(aft_cpk_b,
                            tvertnet[5:7,i],
                            maxlevel=maxlevel,
                            yindex      = 1,
                            Longindex   = 2,
                            iscpk       = True,
                            aft         = True)
    rstrv = aft_cpk_b.vertices[0,0] #x vertex is removed from cLP2D
    cpkLspline.restore = {'index':0,
                          'vertices':rstrv, #vector of 0. centerline x
                          2:0,
                          1:1} 
    Lspline = cpkLspline.get_coarsest_Lspline()
    Lspline.restore = cpkLspline.restore
    #append the aft cLProfile Lspline.curve
    lcurvenet_aft.append(
                    fix_Lspline_23D(Lspline) )
    # 
    #--------------------------
    # aft DWL
    i=6
    print 'aft curve {}'.format(i)
    dwlLspline = thbmatch2D(aft_dwl_b,
                            tvertnet[5:7,i],
                            maxlevel    = maxlevel,
                            isdwl       = True,
                            aft         = True)
    """
    This made a THB curve.
    It is correct, but the coarse curve is not the thing to look at!
    
    e.g. use the THB evaluator:
        
    dwlLspline.curve(0.5610096883028746)
    Out[23]: array([[17.70896871, 12.35339247]])
    
    tvertnet[1:3,i]
    Out[24]: 
    array([[12.35339247, 17.33641345, 17.70896871],
           [14.04681792, 17.33641345, 36.89368482]])
    """
    rstrv = aft_dwl_b.vertices[0,1] #y vertex is removed from DWL2D
    dwlLspline.restore = {'index':1,
                          'vertices':rstrv, #dwl height y
                          0:1,
                          2:0} 
    Lspline_ = dwlLspline.get_coarsest_Lspline()
    Lspline_.restore = dwlLspline.restore
    dwlLspline = Lspline_
    #
    #******************************************
    # Now back to 3D solving
    for i in range(1,nnew+1):
        print 'aft curve {}'.format(i)
        #**************************************
        #linear interpolation
        #
        #vanilla_spline = spline.interpolatedBspline(
        #                tvertnet[:,i], k, nump)
        #**************************************
        # simple 4 c.v. starting linear curve
        lvts = linear_vertices(tvertnet[4:,i][0],
                               tvertnet[4:,i][-1],
                               4)
        thbcurve = rbspline(lvts,k,nump)
        thbcurve = thbcurve.dyadic_refinement()
        thbcurve = thbcurve.dyadic_refinement()
        thbcurve.parent = None
        #
        print 'new aft curve ', i
        Lspline = setup_interpolation_solver_parsi(thbcurve,
                                                   tvertnet[5:7,i],
                                                   aft=True)
        Lspline = THBsolver(thbcurve,
                               tvertnet[5:7,i],
                               Lspline = Lspline,
                               maxlevel=maxlevel,
                               normalize=normalize)
        lcurvenet_aft.append(Lspline.curve)
    #******************************************
    #append the aft DWL Lspline.curve
    lcurvenet_aft.append(
            fix_Lspline_23D(dwlLspline) )
    
    return [lcurvenet_fwd,lcurvenet_mid,lcurvenet_aft]




def fix_Lspline_23D(Lspline):
    """
    this : Lspline (needed for the restore {dictionary})
    basecurve : coarsest curve
    """
    this = Lspline.get_coarsest_Lspline()
    restore = this.restore
    #this = Lspline
    #assert(this.parent is None),'error called fix_lspline23d on non-base Lspline'
    lcurve = this.curve
    basecurve = make_thb_3D(lcurve,
                            rdict = restore,
                            parent=None)
    
    while lcurve.children is not None:
        lcurve = lcurve.children
        nlcurve = make_thb_3D(lcurve,
                         rdict  = restore,
                         parent = basecurve)
        basecurve = nlcurve
    return basecurve.get_coarsest_curve()
        
    
##
##*****************************************************************************
##
##
##
##
##
##  End of THB adaptive algorithms
##
##
## 
##
## 
## 
##
##*****************************************************************************
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
##
##*****************************************************************************
## Really find a location of interest on a curve
##
##*****************************************************************************
##
#def find_xy_max(curve):
    

##
##*****************************************************************************
## Curve-Curve matching (possibly to be replaced by knot vector solver)
##
##*****************************************************************************
##
def match_curve(initial_curve,
                LazyFPD = None,
                return_Lspline=False,
                make_thb = False,
                ncpts=None,
                special_func = None,
                nmatch=50,
                wt=1.):
    """Reparameterize a curve to standard form via optimization
    
        input: 
        --------
            initial_curve : curve with bad parameterization 
                            that we'd like to match closely
            
            LazyFPD:
        aka exactmatch_params : dictionary containing
                                name:[parameter,value] for CurvePoint constraints
            
            special_func: dictionary of closures which
                            wrap 'special functions' - special objectives
                            which, in the Lspline optimization loop,
                            only require being passes the vertices
                            and added to the total fairness functional, f.
            
        output: 
        --------
            curve or Lspline object matching the 
            shape of the input but now with a standard knot vector
        
        this is a speciallization of the way
        to match any curve of equal dimension.
        
        Idea:
        --------
            Knot Deformation Solver:
                meta solver which smoothly transitions from 
                one knot vector to another whilst holding the
                curve the same shape.
                
        Example Issue which this match and its LazyFPD aims to solve:
        --------
            In the bare hull generation program, we split the CPkeel.
            This gives a curve of suitable geometry to be the bottom
            lcurve in the longitudinal curve net, 
            but  this curve
            not suitable in terms of knots - in terms of parameterization.
            -Knot removal is failing for this example :(
            -So we optimize.
            -The program there knows what constraints the 'old' curve must meet
            -the function here is intended to be generic (we'll see how close we come)
            -so we've used the LazyFPDConstraintInterface class to 
            turn the constraint setters into closures which await FPD.
            FPD can only be constructed when we have the requisite starting curve.
            -because FPD will, in general, precompute things that are basis dependent.
            
        dev
        -------
            import curve as spline
            from myspline import rbspline
            from ADILS import interval_bounds
            from FormParameter import FormParameterDict as FPD
        
            initial_curve = copy.deepcopy(self.CPK_fwd)
            LazyFPD = CI
            return_Lspline=True
            make_thb = True
            ncpts=self.example_transverse.n
            nmatch=None
            wt=10
    """
    k = initial_curve.k
    nump = initial_curve.nump
    
    if verbose:print 'utility function: match_curve'
    
    if ncpts is None: 
        ncpts = initial_curve.n
    elif ncpts !=  initial_curve.n:
        assert(ncpts > initial_curve.n),'initial curve has to many points, cannot match (using requested number of points) without reduced fidelity'
        diff = ncpts - initial_curve.n
        for i in range(diff):
            initial_curve.knot_insertion(.5)
    
    if make_thb:
        curve = rbspline(initial_curve.vertices,
                         k=k,
                         nump=nump)
    else:
        curve = spline.Bspline(initial_curve.vertices,
                               k=k,
                               nump=nump)
    #
    #***********************************
    # 
    interval_data, small = interval_bounds(curve)
    #
    #***********************************
    # 
    FPD = FormParameterDict(curve)
    if LazyFPD is not None:
        #FPD = FPDin 
        #a preset FPD would have a ------------bad knot vector--------------.
        #FPD.curve = #setting it one level up in hull_from_simple_designspace... 
        #       ---------------to messy--------------
        # again functional continuation would be nice:
        # here set the form parameters in the hull program
        # whilst delaying choice of curve and delaying precomputations
        # which depend on knot vectors.  Oh we have a class for that
        # ADILS.LazyFPDConstraintInterface
        # should have done it this way all along
        # as it seperates constraints from curves
        # and facilitates fast restarts with new information 
        # (when curve knots change, or other such)
        # meanwhile there is still a need for restart  
        # keeping basis data the same (more efficient restart than this)
        for constraint in LazyFPD.lazyconstraints:
            FPD = constraint(FPD)
        
        
    #
    #***********************************
    # 
    FPD.add_E1(kind='LS', weight = .5)
    FPD.add_E2(kind='LS', weight = .5)
    FPD.add_E3(kind='LS', weight = .5)
    FPD.add_ArcLengthApprox(kind='LS', 
                            weight = .5)
    #
    #***********************************
    # 
    if nmatch is not None:
        print 'warning: doing approximate match'
        for pt in np.linspace(.1,.9,nmatch,endpoint=True):
            loc = initial_curve.CurvePoint(pt)
            FPD.add_xPointConstraint(kind= 'LS', 
                                     value=loc[0], 
                                     location=pt ,
                                     weight=wt)
            FPD.add_yPointConstraint(kind= 'LS', 
                                     value=loc[1], 
                                     location=pt ,
                                     weight=wt)
            if initial_curve.dim == 3:
                FPD.add_zPointConstraint(kind= 'LS', 
                                         value=loc[2], 
                                         location=pt ,
                                         weight=wt)
    #
    #***********************************
    # 
    L = Lagrangian(FPD)
    #
    #***********************************
    # reify the special function 
    if special_func is not None:
        for key in special_func:
            L.sp[key] = special_func[key](curve)
    #
    #***********************************
    # 
    interval_data, small = lagrangian_bounds(L, 
                        interval_data, small, 1.e4)
    Lspline = IntervalLagrangeSpline(curve, L, 
                                     data = interval_data)
    #
    #***********************************
    # 
    Lspline.optimize(stop=30)
    
    if return_Lspline:
        return Lspline
    else:
        return Lspline.curve
    
    

##
##*****************************************************************************
##  Curve-Knot-Vector Operations
##
##*****************************************************************************
##
# helper function
def match_curve_knots(curve_to_fix, good_curve, tol=.1):
    """transform the knot vetor of one curve
    to match another knot vector 
    while attempting to leave the initial curve unchanged
    (don't try this without adult supervision, folks.)
    
    taken from hull_from_simple_designspace (just so the names match...)
    """
    retcv = match_knots(initial_curve = copy.copy(curve_to_fix),
                        curve_w_good_properties=good_curve,
                        TOL=.1)
    return retcv
#
# 
def reparameterize(initial_curve, new_knot_vector,
                   CI, steps = 5):
    """Reparameterize a curve 
    via knot (hopefully homotopic) deformation transform 
    transformation from one knot vector to another
    
    ToDo:
        an outer loop which changes the knot vector 
        and reconverges the curve
    
        
    Parameters
    ----------
        CI : LazyFPDConstraintInterface containing the
            Form Parameters which essentially describe the curve
            with bad parameterization
            that we intend to fix 
            with good parameterization
    
    dev:
    ----------
        from ADILS import LazyFPDConstraintInterface
        
        initial_curve = copy.deepcopy(self.CPK_cdr_b)
        
        
        
        self.SAC.plot3DmultiList([self.CPK_cdr_b],
                                 [self.CProfile])
        
        
        self.CProfile.plot3DmultiList([matched_curve],
                                      [self.CProfile])
        
        self.CPK_cdr_b.plotcurve()
        self.CPK_fwd.plotcurve()
        self.CPK_aft.plotcurve()
        self.CProfile2D.plotcurve()
        
        self.SAC.plot3DmultiList(
            [self.CPK_cdr_b,self.CPK_fwd,self.CPK_aft],
            [self.CProfile])
        
        import curve as spline
        
        initial_curve = copy.copy(self.CPK_cdr_b)  ####COPY
        
        new_knot_vector = self.CProfile.t
        CI = CI
        steps = 5
        
        
    """
    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return idx #array[idx]

    badknots = initial_curve.t
    gknots = new_knot_vector
    
    lenbad = len(badknots)
    lengood = len(gknots)
    ith = initial_curve.k
    n = initial_curve.n
    k = initial_curve.k
    candidate_knots = copy.deepcopy(gknots[k:k+n+1])
    extant = badknots[k:n]
    for knot in extant:
        item = find_nearest(candidate_knots, knot )
        a=list(candidate_knots)
        a.pop(item)
        candidate_knots = np.asarray(a)
        
    
    ith = 0
    while ( lenbad < lengood) and ith<steps:
        knot = candidate_knots[ith]
        if knot in badknots:
            pass
        else:
            initial_curve.knot_insertion(knot)
            badknots = initial_curve.t
        lenbad = len(badknots)
        lengood = len(gknots)
        ith +=1
    
    
    knot_deformation_generator = spline.knot_deformation(initial_curve.t,
                                            new_knot_vector,
                                            steps)
    nc = initial_curve
    for step in range(1,steps):
        ithknots = knot_deformation_generator(step)
        #
        #**********************************************************************
        # solve for the reparameterization, optimization method:
        #
        # note: FPD.curve <=> initial_curve
        #
        """#solver style
        Lspline = match_curve(initial_curve,
                                LazyFPD = CI,
                                return_Lspline=True,
                                make_thb = True,
                                ncpts=self.CProfile.n,
                                nmatch=200,
                                wt=10,)
        self.Lsplines.long_keel = Lspline
        matched_curve = Lspline.curve
        """
        nc = match_knots(initial_curve = nc, 
                         goodknots = ithknots,
                         n = initial_curve.n,
                         k = initial_curve.k,
                         TOL=.01)
        
    return nc
    
##
##*****************************************************************************
## Curve-Knot-Vector Operations
##
##*****************************************************************************
##
def match_knots(initial_curve,
                curve_w_good_properties=None,
                goodknots=None,
                n=None,
                k=None,
                TOL=.001):
    """transform the knot vetor of one curve
    to match another knot vector 
    while attempting to leave the initial curve unchanged
    (don't try this without adult supervision, folks.)
    
    Parameters
    ----------
        initial_curve : 
                    curve with 'wrong' parameterization
        curve_w_good_properties : 
                    curve with the right parameterization
        
            
    Returns
    ----------
        curve which matches the initial curve,
        but has the knot vector good_knots.
        
    Notes
    ----------
        It is left to the user of this function
        to make sure that what it attempts is possible.
        
        
    dev
    ----------
        initial_curve = copy.copy(self.CPK_cdr_b) 
        curve_w_good_properties = self.CProfile
        
        
        initial_curve = nc, 
        goodknots = ithknots
        n = initial_curve.n
        k = initial_curve.k
    """
    if curve_w_good_properties is None:
        gknots = goodknots
        #n = n
        #k = k
    else:      
        gknots = curve_w_good_properties.t
        n = curve_w_good_properties.n
        k = curve_w_good_properties.k
        
    def find_knot_to_delete(index, counter=0, maxit=1):
        if initial_curve.t[index] not in gknots:
            return index
        elif (counter<maxit): #remember knot counts go up after insertion!
            c = counter+1
            if index <len(initial_curve.t)-2:
                return find_knot_to_delete(index+1,
                                           counter=c,
                                           maxit=maxit)
            else:
                return find_knot_to_delete(initial_curve.k,
                                           counter=c,
                                           maxit=maxit)
        else:
            print 'could not find a knot to delete'
            print 'ERROR: warning, multiplicities not accounted for!'
            return None
    
    delete_guess = initial_curve.k 
                    # we are going to search
                    # for knots to remove 
                    # because this is more lazy than
                    # making knot_insertion return __where__ the new
                    # knots were placed in the first place
                    # because
                    # programmer time is becoming
                    # the prime mover  ;-) ... so I thought  :(
    miter = initial_curve.n+initial_curve.k
    for i,gknot in enumerate(gknots[k:n]):
        if gknot in initial_curve.t:
            pass
        else:
            print 'substituting new {} '.format(gknot)
            dindex = find_knot_to_delete(delete_guess,counter=0,
                                maxit=miter)
            if dindex is not None and initial_curve.t[dindex] != 0. \
                                    and initial_curve.t[dindex] != 1.:
                
                print ' in place of old {}'.format(initial_curve.t[dindex])
                if gknot < initial_curve.t[dindex]:
                    initial_curve.knot_insertion(gknot)
                    initial_curve.knot_removal(dindex+1,
                                               tol=TOL)
                elif gknot > initial_curve.t[dindex]:
                    initial_curve.knot_insertion(gknot)
                    initial_curve.knot_removal(dindex,
                                               tol=TOL)
                
    return initial_curve



##
##*****************************************************************************
## Point Interpolation
##
##*****************************************************************************
##      
def interpolate_vertices(initial_vertices,
                         helper_curve,
                         dk = 4,
                         dnump = 30,
                         return_Lspline=False,
                         make_thb = False,
                         nose_constraint=False,
                         Bmax=None,
                         hullblockextremes=None,
                         exact_fit = True,
                         dropout = None,
                         ignore = None,
                         fixities = None):
    """
    Lspline optimization to pass a new curve through 
        a set of initial vertices
        (quite often this function returns 'lcurvenet' curves)
        

    Parameters
    ----------
    
        initial_vertices: -vertices we want to interpolate
        
        helper_curve:  -curve similar to the one we want
                        in that it interpolates the initial vertices
                        but is not as nice - due to parameterization
                        and possibly smoothness
        dnump:      curve.nump number of points to compute on the curve
                    (standard data from curve.Bsplines class)
        
        return_Lspline:     -return the curve solver entire
        make_thb:           -return the curve as a THB style curve
        nose_constraint:    -flat end conditions for the bow of a bare hull
        Bmax:               -max breadth 
                            (to be replaced by hullblockextremes)
        hullblockextremes:  -comprehensive block extremes of the bare hull
                            non-symmetric 'half'
                            
        exact_fit: was playing around with approximate fitting (again)
        
        dropout: realized the issue (duh) was more likely to be
                    in the use of 11 tranverse curves being fit by
                    longitudinals with 11 control vertices total!
                    tune in (to the problem) and drop out! (a control vertex)
                    
    Notes
    ----------
        longitudinals are defined fore to aft
        
    Returns
    ----------
        curve   (reparameterized with a 'standard' knot vector)
    
    this is a speciallization of the way
    to match any curve of equal dimension.
    """
    #if verbose:
    print '\n utility function: interpolate_vertices \n'
    # get the vertices to interpolate:
    ukbar = helper_curve.ukbar  #ukbar are the parameters 
                                #at which the interpolation happens
    #print 'the ukbar are'
    #print ukbar
    #   
    if make_thb:
        curve = rbspline(helper_curve.vertices,
                         k=dk,
                         nump=dnump)
    else:
        curve = spline.Bspline(helper_curve.vertices,
                                 k=dk,
                                 nump=dnump)
    
    #try this on LS all the way here:
    if exact_fit:
        this_kind = 'equality'
    else:
        this_kind = 'LS'
        
    if dropout is None:
        dropout = [1]
    if ignore is None:
        ignore = []
    
    inicurve = copy.deepcopy(curve)
    def initial_form(curve):
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve)
        #
        #**********************************************************************
        # Fairness
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .05)
        #FPD.add_ArcLengthApprox(kind='LS', 
        #                        weight = .5)
        #
        #**********************************************************************
        # nose constraint (actually smooths things - depinding on availibel DOFS):
        FPD.add_ArcLength(kind='LS',
                          weight = 1.5)
        #
        #**********************************************************************
        # interoplation constraints
        #
        #FPD.add_zFixity(index=1,
        #                value=curve.vertices[0,2])
        # initial_vertices : vertices to interpolate (CurvePoint constraints - not vertex constraints ;-)
        # ukbar : 
        #
        start = 1
        count = 1
        for pt,loc in zip(initial_vertices[start:-1],ukbar[start:-1]):
            if count in dropout:
                #
                FPD.add_xPointConstraint(kind= 'LS', 
                                         value=pt[0], 
                                         location=loc, 
                                         weight = 50.) 
                FPD.add_yPointConstraint(kind= 'LS', 
                                         value=pt[1], 
                                         location=loc, 
                                         weight = 50.)
                FPD.add_zPointConstraint(kind= 'LS', 
                                         value=pt[2], 
                                         location=loc, 
                                         weight = 50. )
            elif count not in ignore:
                #                continue #thou shal not pass!!!! -because it  (pass) will exit the loop, hehe
                #            else:
                #
                FPD.add_xPointConstraint(kind= this_kind, 
                                         value=pt[0], 
                                         location=loc, 
                                         weight = 100. ) #these weights do nothing unless kind = 'LS' (least squares)
                FPD.add_yPointConstraint(kind= this_kind, 
                                         value=pt[1], 
                                         location=loc, 
                                         weight = 100.  )
                FPD.add_zPointConstraint(kind= this_kind, 
                                         value=pt[2], 
                                         location=loc, 
                                         weight = 100.  )
            count +=1
        
        #
        #**************************************************************************
        # 
        
        if hullblockextremes is not None:
            for i,pt in enumerate(curve.vertices[1:-1]):
                this = i+1
                #--------------------------------------x
#                diff = abs(hullblockextremes.be - hullblockextremes.bb)
#                FPD.add_xVertexConstraint(kind='min',
#                                          value=hullblockextremes.bb-.1*diff,
#                                          index=this)
#                FPD.add_xVertexConstraint(kind='max',
#                                          value=hullblockextremes.be*1.3,
#                                          index=this)
                #--------------------------------------y
#                diff = abs(hullblockextremes.de-hullblockextremes.db)
#                FPD.add_yVertexConstraint(kind='min',
#                                          value=hullblockextremes.db-.15*diff,
#                                          index=this)
#                FPD.add_yVertexConstraint(kind='max',
#                                          value=hullblockextremes.de*1.15,
#                                          index=this)
                #--------------------------------------z
#                diff = abs(hullblockextremes.le - hullblockextremes.lb)
#                FPD.add_zVertexConstraint(kind='min',
#                                          value=hullblockextremes.lb-.15*diff,
#                                          index=this)
#                FPD.add_zVertexConstraint(kind='max',
#                                          value=hullblockextremes.le*1.15,
#                                          index=this)
          
        #
        #**************************************************************************
        # Fixities
        #
        if fixities is not None:
            for ix,valx in fixities[0]:
                FPD.add_xFixity(index=ix,
                                value=valx)
        #
        #**************************************************************************
        #
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                            interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, 
                                         data = interval_data)
        #Lspline.curve.pts_M_pts()
        #Lspline.curve.compute_arclength()
        #
        #**************************************************************************
        #
        print 'lcurve initial_form solve:'
        Lspline.optimize(stop=15)
        return Lspline
    
    
    
    def simplified_form(Lspline):
        norms = package_norms(Lspline)
        E1_0 = norms[0]
        E2_0 = norms[1]
        E3_0 = norms[2]
        S_0  = norms[3]
        curve = Lspline.curve
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve)
        #
        #**************************************************************************
        # interoplation constraints
        #
        this_kind = 'LS'
        #
        start = 1
        count = 1
        for pt,loc in zip(initial_vertices[start:-1],ukbar[start:-1]):
            if count in dropout:
                #
                FPD.add_xPointConstraint(kind= 'LS', 
                                         value=pt[0], 
                                         location=loc, 
                                         weight = 1.5/E1_0 ) #
                FPD.add_yPointConstraint(kind= 'LS', 
                                         value=pt[1], 
                                         location=loc, 
                                         weight = 1.5/E1_0 )
                FPD.add_zPointConstraint(kind= 'LS', 
                                         value=pt[2], 
                                         location=loc, 
                                         weight = 1.5/E1_0  )
            elif count not in ignore:
                #                continue #thou shal not pass!!!! - because it (pass) will exit the loop, hehe
                #            else:
                FPD.add_xPointConstraint(kind= this_kind,
                                         value=pt[0], 
                                         location=loc,
                                         weight=10./E1_0)
                FPD.add_yPointConstraint(kind= this_kind, 
                                         value=pt[1], 
                                         location=loc,
                                         weight=10./E1_0 )
                FPD.add_zPointConstraint(kind= this_kind, 
                                         value=pt[2], 
                                         location=loc,
                                         weight=10./E1_0 )
            count+=1
                
        #
        #**********************************************************************
        # block extrema?
        #
        # be careful here:
        #  the extrema of the hull are not necessarily the max extents
        #  that the lcurve vertices should be allowed to excurt to
        #  -maybe better to change this to max extents for 
        #   __actual_points_on_the_curves_at_even_intervals__
        #  -for now we just give it a little lee-way
        
        
        if hullblockextremes is not None:
            for i,pt in enumerate(curve.vertices[1:-1]):
                this = i+1
                #--------------------------------------x
#                diff = abs(hullblockextremes.be - hullblockextremes.bb)
#                FPD.add_xVertexConstraint(kind='min',
#                                          value=hullblockextremes.bb-.5*diff,
#                                          index=this)
#                FPD.add_xVertexConstraint(kind='max',
#                                          value=hullblockextremes.be*1.5,
#                                          index=this)
                #--------------------------------------y
#                diff = abs(hullblockextremes.de-hullblockextremes.db)
#                FPD.add_yVertexConstraint(kind='min',
#                                          value=hullblockextremes.db-.5*diff,
#                                          index=this)
#                FPD.add_yVertexConstraint(kind='max',
#                                          value=hullblockextremes.de*1.5,
#                                          index=this)
#                #--------------------------------------z
#                diff = abs(hullblockextremes.le - hullblockextremes.lb)
#                FPD.add_zVertexConstraint(kind='min',
#                                          value=hullblockextremes.lb-.5*diff,
#                                          index=this)
#                FPD.add_zVertexConstraint(kind='max',
#                                          value=hullblockextremes.le*1.5,
#                                          index=this)
                #--------------------------------------done
        #
        #**************************************************************************
        # Fairness
        FPD.add_E1(kind='LS', weight = .5/E1_0)
        FPD.add_E2(kind='LS', weight = .5/E2_0)
        FPD.add_E3(kind='LS', weight = .01/E3_0)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1./S_0)
        #
        #
        #**************************************************************************
        #
        interval_data, small = interval_bounds(curve)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                            interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, 
                                         data = interval_data)
        #
        #**************************************************************************
        # Optimize
        print 'lcurve fixing solve:'
        Lspline.optimize(stop=15)
        return Lspline
    
    
    
    def nose_constraint_form(Lspline):
        """ leading edge of the bare hull shape 
        expressed in the longitudinal curve list
        
        Actually, how about solving the nose constraint
        in the THB refined space only!
        
        """
        message = 'ERROR in utility_optimization: send Bmax to: \n function <interpolate_vertices>'
        assert(Bmax is not None),message
        #
        norms = package_norms(Lspline)
        E1_0 = norms[0]
        E2_0 = norms[1]
        E3_0 = norms[2]
        S_0  = norms[3]
        curve = Lspline.curve
        #
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve)
        #
        #**************************************************************************
        # make the nose:
        #
        #        FPD.add_xVertexConstraint(kind='equality',
        #                                 index = 1,
        #                                 value = 0.)
        #
        #        FPD.add_relative_xVertexConstraint(
        #                                 kind='equality',
        #                                 location=None,
        #                                 index = 0,
        #                                 index2 = 1,
        #                                 value = 0.,
        #                                 weight=1.)
        FPD.add_xFixity(index=1,
                        value=curve.vertices[0,0])
        #
        FPD.add_yVertexConstraint(kind='min',
                                 index = 1,
                                 value = 0.)
        #        FPD.add_yVertexConstraint(kind='LS',
        #                                 index = 1,
        #                                 value = Bmax)
        #        FPD.add_yVertexConstraint(kind='LS',
        #                                 index = 2,
        #                                 value = Bmax)
        #
        #**********************************************************************
        # interoplation constraints
        for pt,loc in zip(initial_vertices[1:-1],ukbar[1:-1]):
            #print 'pt = ', pt
            FPD.add_xPointConstraint(kind= this_kind, 
                                     value=pt[0], 
                                     location=loc )
            FPD.add_yPointConstraint(kind= this_kind, 
                                     value=pt[1], 
                                     location=loc )
            FPD.add_zPointConstraint(kind= this_kind, 
                                     value=pt[2], 
                                     location=loc )   
        #
        #**********************************************************************
        # block extrema?
        #
        # be careful here:
        #  the extrema of the hull are not necessarily the max extents
        #  that the lcurve vertices should be allowed to excurt to
        #  -maybe better to change this to max extents for 
        #   __actual_points_on_the_curves_at_even_intervals__
        #  -for now we just give it a little lee-way
        if hullblockextremes is not None:
            for i,pt in enumerate(curve.vertices[1:-1]):
                this = i+1
                #--------------------------------------x
                FPD.add_xVertexConstraint(kind='min',
                                          weight=.5,
                                          value=hullblockextremes.bb,
                                          index=this)
                FPD.add_xVertexConstraint(kind='max',
                                          weight=.5,
                                          value=hullblockextremes.be*1.5,
                                          index=this)
                #--------------------------------------y
                FPD.add_yVertexConstraint(kind='min',
                                          weight=.5,
                                          value=hullblockextremes.db,
                                          index=this)
                FPD.add_yVertexConstraint(kind='max',
                                          weight=.5,
                                          value=hullblockextremes.de*1.5,
                                          index=this)
                #--------------------------------------z
                FPD.add_zVertexConstraint(kind='min',
                                          weight=.5,
                                          value=hullblockextremes.lb,
                                          index=this)
                FPD.add_zVertexConstraint(kind='max',
                                          weight=.5,
                                          value=hullblockextremes.le*1.5,
                                          index=this)
                #--------------------------------------done
        
        #
        #**************************************************************************
        # Fairness
        FPD.add_E1(kind='LS', weight = .5/E1_0)
        FPD.add_E2(kind='LS', weight = .5/E2_0)
        FPD.add_E3(kind='LS', weight = .5/E3_0)
        FPD.add_ArcLengthApprox(kind='LS', weight = 5./S_0)
        #
        #**************************************************************************
        #
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                            interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, 
                                         data = interval_data)
        #
        #**************************************************************************
        # Optimize
        print 'lcurve stage 2 solve (nose constraint):'
        Lspline.optimize(stop=25)
        #
        #**************************************************************************
        # Done
        return Lspline
    
    
    def getckit(Lspline):
        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
            ckit =False
        else:
            ckit = True 
        return ckit
        
    
    #
    #**********************************************************************
    # initialize
    ckit = False
    #
    #**********************************************************************
    # solve Lcurve stage 1
    Lspline = initial_form(curve)
    #
    # save
    Lsave1 = copy.deepcopy(Lspline)
    Lsave_ini = copy.deepcopy(Lspline)
    #
    #**********************************************************************
    # check Lcurve stage 1
    ckit = getckit(Lspline)
    if ckit:
        pass
    else:
        print 'switching to linear interpolation curves'
        #Lspline.curve = inicurve 
        Lspline = simplified_form(Lspline)  
    #
    #**********************************************************************
    # Lcurve nose constraint  (LEAVE for THB!)
    #    if False:
    #        if nose_constraint:
    #            Lspline = nose_constraint_form(Lspline)
    #            if Lspline.conv > Lspline.tol: #TODO: fix need for this due to stupid delta_z = NAN bail out:
    #                Lspline.error_code = 'max_iter'
    #                
    #                
    #            if Lspline.error_code is 'NAN':
    #                if ckit:
    #                    Lspline = Lsave1
    #                    Lspline = simplified_form(Lspline)
    #                else: #initial Lspline is bad
    #                    Lspline = initial_form(curve) #give up
    #            elif Lspline.error_code is 'max_iter':
    #                if Lspline.conv>.01:
    #                    if ckit:
    #                        Lspline = Lsave1
    #                        Lspline = simplified_form(Lspline)
    #                    else: #initial Lspline is bad
    #                        Lspline = initial_form(curve) #give up
    #                else:
    #                    pass #take Lspline from solution 2!
                
    #
    #**********************************************************************
    #  no nose constraints at this level of detail!
    #    if ckit:
    #        pass #go with these
    #        #        Lspline = simplified_form(Lspline)
    #        #        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
    #        #            ckit =False
    #        #            #Lspline = Lsave1
    #        #        else:
    #        #            ckit = True 
    #        #            #Lspline is simplified form Lspline
    #    else:
    #        Lspline = simplified_form(Lspline)  
    #        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
    #            print 'switching to linear interpolation curves'
    #            ckit =False
    #            Lspline.curve = inicurve #use linear interpolation
    #            Lspline.Lagrangian.obj = {}
    #        else:
    #            ckit = True 
    #            #Lspline is simplified form Lspline

    #---
    #TODO: only put the nose on THB fine level(s)!
            
            
    
    print 'Lspline interpolated curve done'
            
    if return_Lspline:
        return Lspline
    else:
        return Lspline.curve


##
##
##  Extended fixities
  
def interpolate_vertices_fixany(initial_vertices,
                         helper_curve,
                         dk = 4,
                         dnump = 30,
                         return_Lspline=False,
                         make_thb = False,
                         nose_constraint=False,
                         Bmax=None,
                         hullblockextremes=None,
                         exact_fit = True,
                         dropout = None,
                         ignore = None,
                         fixities = None):
    """
    Lspline optimization to pass a new curve through 
        a set of initial vertices
        (quite often this function returns 'lcurvenet' curves)
        

    Parameters
    ----------
    
        initial_vertices: -vertices we want to interpolate
        
        helper_curve:  -curve similar to the one we want
                        in that it interpolates the initial vertices
                        but is not as nice - due to parameterization
                        and possibly smoothness
        dnump:      curve.nump number of points to compute on the curve
                    (standard data from curve.Bsplines class)
        
        return_Lspline:     -return the curve solver entire
        make_thb:           -return the curve as a THB style curve
        nose_constraint:    -flat end conditions for the bow of a bare hull
        Bmax:               -max breadth 
                            (to be replaced by hullblockextremes)
        hullblockextremes:  -comprehensive block extremes of the bare hull
                            non-symmetric 'half'
                            
        exact_fit: was playing around with approximate fitting (again)
        
        dropout: realized the issue (duh) was more likely to be
                    in the use of 11 tranverse curves being fit by
                    longitudinals with 11 control vertices total!
                    tune in (to the problem) and drop out! (a control vertex)
                    
    Notes
    ----------
        longitudinals are defined fore to aft
        
    Returns
    ----------
        curve   (reparameterized with a 'standard' knot vector)
    
    this is a speciallization of the way
    to match any curve of equal dimension.
    """
    #if verbose:
    print '\n utility function: interpolate_vertices \n'
    # get the vertices to interpolate:
    ukbar = helper_curve.ukbar  #ukbar are the parameters 
                                #at which the interpolation happens
    #print 'the ukbar are'
    #print ukbar
    #   
    if make_thb:
        curve = rbspline(helper_curve.vertices,
                         k=dk,
                         nump=dnump)
    else:
        curve = spline.Bspline(helper_curve.vertices,
                                 k=dk,
                                 nump=dnump)
    
    #try this on LS all the way here:
    if exact_fit:
        this_kind = 'equality'
    else:
        this_kind = 'LS'
        
    if dropout is None:
        dropout = [1]
    if ignore is None:
        ignore = []
    
    inicurve = copy.deepcopy(curve)
    def initial_form(curve):
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve)
        #
        #**********************************************************************
        # Fairness
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .05)
        #FPD.add_ArcLengthApprox(kind='LS', 
        #                        weight = .5)
        #
        #**********************************************************************
        # nose constraint (actually smooths things - depinding on availibel DOFS):
        FPD.add_ArcLength(kind='LS',
                          weight = 1.5)
        #
        #**********************************************************************
        # interoplation constraints
        #
        #FPD.add_zFixity(index=1,
        #                value=curve.vertices[0,2])
        # initial_vertices : vertices to interpolate (CurvePoint constraints - not vertex constraints ;-)
        # ukbar : 
        #
        start = 1
        count = 1
        for pt,loc in zip(initial_vertices[start:-1],ukbar[start:-1]):
            if count in dropout:
                #
                FPD.add_xPointConstraint(kind= 'LS', 
                                         value=pt[0], 
                                         location=loc, 
                                         weight = 50.) 
                FPD.add_yPointConstraint(kind= 'LS', 
                                         value=pt[1], 
                                         location=loc, 
                                         weight = 50.)
                FPD.add_zPointConstraint(kind= 'LS', 
                                         value=pt[2], 
                                         location=loc, 
                                         weight = 50. )
            elif count not in ignore:
                #                continue #thou shal not pass!!!! -because it  (pass) will exit the loop, hehe
                #            else:
                #
                FPD.add_xPointConstraint(kind= this_kind, 
                                         value=pt[0], 
                                         location=loc, 
                                         weight = 100. ) #these weights do nothing unless kind = 'LS' (least squares)
                FPD.add_yPointConstraint(kind= this_kind, 
                                         value=pt[1], 
                                         location=loc, 
                                         weight = 100.  )
                FPD.add_zPointConstraint(kind= this_kind, 
                                         value=pt[2], 
                                         location=loc, 
                                         weight = 100.  )
            count +=1
        
        #
        #**************************************************************************
        # 
        
        if hullblockextremes is not None:
            for i,pt in enumerate(curve.vertices[1:-1]):
                this = i+1
                #--------------------------------------x
#                diff = abs(hullblockextremes.be - hullblockextremes.bb)
#                FPD.add_xVertexConstraint(kind='min',
#                                          value=hullblockextremes.bb-.1*diff,
#                                          index=this)
#                FPD.add_xVertexConstraint(kind='max',
#                                          value=hullblockextremes.be*1.3,
#                                          index=this)
                #--------------------------------------y
#                diff = abs(hullblockextremes.de-hullblockextremes.db)
#                FPD.add_yVertexConstraint(kind='min',
#                                          value=hullblockextremes.db-.15*diff,
#                                          index=this)
#                FPD.add_yVertexConstraint(kind='max',
#                                          value=hullblockextremes.de*1.15,
#                                          index=this)
                #--------------------------------------z
#                diff = abs(hullblockextremes.le - hullblockextremes.lb)
#                FPD.add_zVertexConstraint(kind='min',
#                                          value=hullblockextremes.lb-.15*diff,
#                                          index=this)
#                FPD.add_zVertexConstraint(kind='max',
#                                          value=hullblockextremes.le*1.15,
#                                          index=this)
          
        #
        #**************************************************************************
        # Fixities
        #
        if fixities is not None:
            for ix,valx in fixities[0]:
                FPD.add_xFixity(index=ix,
                                value=valx)
            for ix,valx in fixities[1]:
                FPD.add_yFixity(index=ix,
                                value=valx)
            for ix,valx in fixities[2]:
                FPD.add_zFixity(index=ix,
                                value=valx)
        #
        #**************************************************************************
        #
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                            interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, 
                                         data = interval_data)
        #Lspline.curve.pts_M_pts()
        #Lspline.curve.compute_arclength()
        #
        #**************************************************************************
        #
        print 'lcurve initial_form solve:'
        Lspline.optimize(stop=15)
        return Lspline
    
    
    
    def simplified_form(Lspline):
        norms = package_norms(Lspline)
        E1_0 = norms[0]
        E2_0 = norms[1]
        E3_0 = norms[2]
        S_0  = norms[3]
        curve = Lspline.curve
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve)
        #
        #**************************************************************************
        # interoplation constraints
        #
        this_kind = 'LS'
        #
        start = 1
        count = 1
        for pt,loc in zip(initial_vertices[start:-1],ukbar[start:-1]):
            if count in dropout:
                #
                FPD.add_xPointConstraint(kind= 'LS', 
                                         value=pt[0], 
                                         location=loc, 
                                         weight = 1.5/E1_0 ) #
                FPD.add_yPointConstraint(kind= 'LS', 
                                         value=pt[1], 
                                         location=loc, 
                                         weight = 1.5/E1_0 )
                FPD.add_zPointConstraint(kind= 'LS', 
                                         value=pt[2], 
                                         location=loc, 
                                         weight = 1.5/E1_0  )
            elif count not in ignore:
                #                continue #thou shal not pass!!!! - because it (pass) will exit the loop, hehe
                #            else:
                FPD.add_xPointConstraint(kind= this_kind,
                                         value=pt[0], 
                                         location=loc,
                                         weight=10./E1_0)
                FPD.add_yPointConstraint(kind= this_kind, 
                                         value=pt[1], 
                                         location=loc,
                                         weight=10./E1_0 )
                FPD.add_zPointConstraint(kind= this_kind, 
                                         value=pt[2], 
                                         location=loc,
                                         weight=10./E1_0 )
            count+=1
                
        #
        #**********************************************************************
        # block extrema?
        #
        # be careful here:
        #  the extrema of the hull are not necessarily the max extents
        #  that the lcurve vertices should be allowed to excurt to
        #  -maybe better to change this to max extents for 
        #   __actual_points_on_the_curves_at_even_intervals__
        #  -for now we just give it a little lee-way
        
        
        if hullblockextremes is not None:
            for i,pt in enumerate(curve.vertices[1:-1]):
                this = i+1
                #--------------------------------------x
#                diff = abs(hullblockextremes.be - hullblockextremes.bb)
#                FPD.add_xVertexConstraint(kind='min',
#                                          value=hullblockextremes.bb-.5*diff,
#                                          index=this)
#                FPD.add_xVertexConstraint(kind='max',
#                                          value=hullblockextremes.be*1.5,
#                                          index=this)
                #--------------------------------------y
#                diff = abs(hullblockextremes.de-hullblockextremes.db)
#                FPD.add_yVertexConstraint(kind='min',
#                                          value=hullblockextremes.db-.5*diff,
#                                          index=this)
#                FPD.add_yVertexConstraint(kind='max',
#                                          value=hullblockextremes.de*1.5,
#                                          index=this)
#                #--------------------------------------z
#                diff = abs(hullblockextremes.le - hullblockextremes.lb)
#                FPD.add_zVertexConstraint(kind='min',
#                                          value=hullblockextremes.lb-.5*diff,
#                                          index=this)
#                FPD.add_zVertexConstraint(kind='max',
#                                          value=hullblockextremes.le*1.5,
#                                          index=this)
                #--------------------------------------done
          
        #
        #**************************************************************************
        # Fixities
        #
        if fixities is not None:
            for ix,valx in fixities[0]:
                FPD.add_xFixity(index=ix,
                                value=valx)
            for ix,valx in fixities[1]:
                FPD.add_yFixity(index=ix,
                                value=valx)
            for ix,valx in fixities[2]:
                FPD.add_zFixity(index=ix,
                                value=valx)
        #
        #**************************************************************************
        # Fairness
        FPD.add_E1(kind='LS', weight = .5/E1_0)
        FPD.add_E2(kind='LS', weight = .5/E2_0)
        FPD.add_E3(kind='LS', weight = .01/E3_0)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1./S_0)
        #
        #
        #**************************************************************************
        #
        interval_data, small = interval_bounds(curve)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                            interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, 
                                         data = interval_data)
        #
        #**************************************************************************
        # Optimize
        print 'lcurve fixing solve:'
        Lspline.optimize(stop=15)
        return Lspline
    
    
    
    def nose_constraint_form(Lspline):
        """ leading edge of the bare hull shape 
        expressed in the longitudinal curve list
        
        Actually, how about solving the nose constraint
        in the THB refined space only!
        
        """
        message = 'ERROR in utility_optimization: send Bmax to: \n function <interpolate_vertices>'
        assert(Bmax is not None),message
        #
        norms = package_norms(Lspline)
        E1_0 = norms[0]
        E2_0 = norms[1]
        E3_0 = norms[2]
        S_0  = norms[3]
        curve = Lspline.curve
        #
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve)
        #
        #**************************************************************************
        # make the nose:
        #
        #        FPD.add_xVertexConstraint(kind='equality',
        #                                 index = 1,
        #                                 value = 0.)
        #
        #        FPD.add_relative_xVertexConstraint(
        #                                 kind='equality',
        #                                 location=None,
        #                                 index = 0,
        #                                 index2 = 1,
        #                                 value = 0.,
        #                                 weight=1.)
        FPD.add_xFixity(index=1,
                        value=curve.vertices[0,0])
        #
        FPD.add_yVertexConstraint(kind='min',
                                 index = 1,
                                 value = 0.)
        #        FPD.add_yVertexConstraint(kind='LS',
        #                                 index = 1,
        #                                 value = Bmax)
        #        FPD.add_yVertexConstraint(kind='LS',
        #                                 index = 2,
        #                                 value = Bmax)
        #
        #**********************************************************************
        # interoplation constraints
        for pt,loc in zip(initial_vertices[1:-1],ukbar[1:-1]):
            #print 'pt = ', pt
            FPD.add_xPointConstraint(kind= this_kind, 
                                     value=pt[0], 
                                     location=loc )
            FPD.add_yPointConstraint(kind= this_kind, 
                                     value=pt[1], 
                                     location=loc )
            FPD.add_zPointConstraint(kind= this_kind, 
                                     value=pt[2], 
                                     location=loc )   
        #
        #**********************************************************************
        # block extrema?
        #
        # be careful here:
        #  the extrema of the hull are not necessarily the max extents
        #  that the lcurve vertices should be allowed to excurt to
        #  -maybe better to change this to max extents for 
        #   __actual_points_on_the_curves_at_even_intervals__
        #  -for now we just give it a little lee-way
        if hullblockextremes is not None:
            for i,pt in enumerate(curve.vertices[1:-1]):
                this = i+1
                #--------------------------------------x
                FPD.add_xVertexConstraint(kind='min',
                                          weight=.5,
                                          value=hullblockextremes.bb,
                                          index=this)
                FPD.add_xVertexConstraint(kind='max',
                                          weight=.5,
                                          value=hullblockextremes.be*1.5,
                                          index=this)
                #--------------------------------------y
                FPD.add_yVertexConstraint(kind='min',
                                          weight=.5,
                                          value=hullblockextremes.db,
                                          index=this)
                FPD.add_yVertexConstraint(kind='max',
                                          weight=.5,
                                          value=hullblockextremes.de*1.5,
                                          index=this)
                #--------------------------------------z
                FPD.add_zVertexConstraint(kind='min',
                                          weight=.5,
                                          value=hullblockextremes.lb,
                                          index=this)
                FPD.add_zVertexConstraint(kind='max',
                                          weight=.5,
                                          value=hullblockextremes.le*1.5,
                                          index=this)
                #--------------------------------------done
        
        #
        #**************************************************************************
        # Fairness
        FPD.add_E1(kind='LS', weight = .5/E1_0)
        FPD.add_E2(kind='LS', weight = .5/E2_0)
        FPD.add_E3(kind='LS', weight = .5/E3_0)
        FPD.add_ArcLengthApprox(kind='LS', weight = 5./S_0)
        #
        #**************************************************************************
        #
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                            interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, 
                                         data = interval_data)
        #
        #**************************************************************************
        # Optimize
        print 'lcurve stage 2 solve (nose constraint):'
        Lspline.optimize(stop=25)
        #
        #**************************************************************************
        # Done
        return Lspline
    
    
    def getckit(Lspline):
        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
            ckit =False
        else:
            ckit = True 
        return ckit
        
    
    #
    #**********************************************************************
    # initialize
    ckit = False
    #
    #**********************************************************************
    # solve Lcurve stage 1
    Lspline = initial_form(curve)
    #
    # save
    Lsave1 = copy.deepcopy(Lspline)
    Lsave_ini = copy.deepcopy(Lspline)
    #
    #**********************************************************************
    # check Lcurve stage 1
    ckit = getckit(Lspline)
    if ckit:
        pass
    else:
        print 'switching to linear interpolation curves'
        #Lspline.curve = inicurve 
        Lspline = simplified_form(Lspline)  
    #
    #**********************************************************************
    # Lcurve nose constraint  (LEAVE for THB!)
    #    if False:
    #        if nose_constraint:
    #            Lspline = nose_constraint_form(Lspline)
    #            if Lspline.conv > Lspline.tol: #TODO: fix need for this due to stupid delta_z = NAN bail out:
    #                Lspline.error_code = 'max_iter'
    #                
    #                
    #            if Lspline.error_code is 'NAN':
    #                if ckit:
    #                    Lspline = Lsave1
    #                    Lspline = simplified_form(Lspline)
    #                else: #initial Lspline is bad
    #                    Lspline = initial_form(curve) #give up
    #            elif Lspline.error_code is 'max_iter':
    #                if Lspline.conv>.01:
    #                    if ckit:
    #                        Lspline = Lsave1
    #                        Lspline = simplified_form(Lspline)
    #                    else: #initial Lspline is bad
    #                        Lspline = initial_form(curve) #give up
    #                else:
    #                    pass #take Lspline from solution 2!
                
    #
    #**********************************************************************
    #  no nose constraints at this level of detail!
    #    if ckit:
    #        pass #go with these
    #        #        Lspline = simplified_form(Lspline)
    #        #        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
    #        #            ckit =False
    #        #            #Lspline = Lsave1
    #        #        else:
    #        #            ckit = True 
    #        #            #Lspline is simplified form Lspline
    #    else:
    #        Lspline = simplified_form(Lspline)  
    #        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
    #            print 'switching to linear interpolation curves'
    #            ckit =False
    #            Lspline.curve = inicurve #use linear interpolation
    #            Lspline.Lagrangian.obj = {}
    #        else:
    #            ckit = True 
    #            #Lspline is simplified form Lspline

    #---
    #TODO: only put the nose on THB fine level(s)!
            
            
    
    print 'Lspline interpolated curve done'
            
    if return_Lspline:
        return Lspline
    else:
        return Lspline.curve



        

##
##*****************************************************************************
##
##   THB Optimization of the Nose Constraint
##
##
##*****************************************************************************
def interpolate_vertices_THB_optimize(interpolate_vertices,
                                     helper_curve,
                                     dk = 4,
                                     dnump = 30,
                                     return_Lspline=True,
                                     make_thb = True,
                                     nose_constraint=True,
                                     Bmax=None,
                                     hullblockextremes=None):
    """NOT USED YET!
    
    At some point we will be FPD'ing the THB representation
    Thoughts:
        -lots of vertices are not active - we just want to optimze
        the active ones... And LEAVE THE OTHERS WHERE THEY ARE
        -because those cannot change - they are part of a 'lower level'
    
    Lspline optimization to pass a new curve through 
        a set of initial vertices
        (quite often this function returns 'lcurvenet' curves)
        

    Parameters
    ----------
    
        interpolate_vertices:   -vertices we want to interpolate
        
        helper_curve:       -THB representation of the actual hull longitudinals
        
        dnump:      curve.nump number of points to compute on the curve
                    (standard data from curve.Bsplines class)
        
        return_Lspline:     -return the curve solver entire
        make_thb:           -return the curve as a THB style curve
        nose_constraint:    -flat end conditions for the bow of a bare hull
        Bmax:               -max breadth 
                            (to be replaced by hullblockextremes)
        hullblockextremes:  -comprehensive block extremes of the bare hull
                            non-symmetric 'half'
        
    Returns
    ----------
        curve   (reparameterized with a 'standard' knot vector)
    
    this is a speciallization of the way
    to match any curve of equal dimension.
    """
    #if verbose:
    print '\n utility function: interpolate_vertices \n'
    # get the vertices to interpolate:
    ukbar = helper_curve.ukbar  #ukbar are the parameters 
                                #at which the interpolation happens
    #print 'the ukbar are'
    #print ukbar
    #   
    if make_thb:
        curve = rbspline(helper_curve.vertices,
                         k=dk,
                         nump=dnump)
    else:
        curve = spline.Bspline(helper_curve.vertices,
                                 k=dk,
                                 nump=dnump)
    
    
    
    inicurve = copy.deepcopy(curve)
    def initial_form(curve):
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve)
        #
        #**************************************************************************
        # Fairness
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .01)
        FPD.add_ArcLengthApprox(kind='LS', 
                                weight = 1.)
        #
        #**************************************************************************
        # interoplation constraints
        #
        # interpolate_vertices : vertices to interpolate
        # ukbar : 
        #
        this_kind = 'equality'
        #this_kind = 'LS'
        for pt,loc in zip(interpolate_vertices[1:-1],ukbar[1:-1]):
            #print 'pt = ', pt
            FPD.add_xPointConstraint(kind= this_kind, 
                                     value=pt[0], 
                                     location=loc )
            FPD.add_yPointConstraint(kind= this_kind, 
                                     value=pt[1], 
                                     location=loc )
            FPD.add_zPointConstraint(kind= this_kind, 
                                     value=pt[2], 
                                     location=loc )
        #
        #**************************************************************************
        #
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                            interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, 
                                         data = interval_data)
        Lspline.curve.pts_M_pts()
        Lspline.curve.compute_arclength()
        #
        #**************************************************************************
        #
        print 'lcurve initial_form solve:'
        Lspline.optimize(stop=30)
        return Lspline
    
    
    
    def simplified_form(Lspline):
        norms = package_norms(Lspline)
        E1_0 = norms[0]
        E2_0 = norms[1]
        E3_0 = norms[2]
        S_0  = norms[3]
        curve = Lspline.curve
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve)
        #
        #**************************************************************************
        # interoplation constraints
        this_kind = 'equality'
        #this_kind = 'LS'
        for pt,loc in zip(interpolate_vertices[1:-1],ukbar[1:-1]):
            #print 'pt = ', pt
            FPD.add_xPointConstraint(kind= this_kind,
                                     value=pt[0], 
                                     location=loc )
            FPD.add_yPointConstraint(kind= this_kind, 
                                     value=pt[1], 
                                     location=loc )
            FPD.add_zPointConstraint(kind= this_kind, 
                                     value=pt[2], 
                                     location=loc )
        #
        #**********************************************************************
        # block extrema?
        #
        # be careful here:
        #  the extrema of the hull are not necessarily the max extents
        #  that the lcurve vertices should be allowed to excurt to
        #  -maybe better to change this to max extents for 
        #   __actual_points_on_the_curves_at_even_intervals__
        #  -for now we just give it a little lee-way
        if hullblockextremes is not None:
            for i,pt in enumerate(curve.vertices[1:-1]):
                this = i+1
                #--------------------------------------x
                FPD.add_xVertexConstraint(kind='min',
                                          value=hullblockextremes.bb,
                                          index=this, 
                                          weight = 1./E1_0)
                FPD.add_xVertexConstraint(kind='max',
                                          value=hullblockextremes.be*1.5,
                                          index=this, 
                                          weight = 1./E1_0)
                #--------------------------------------y
                FPD.add_yVertexConstraint(kind='min',
                                          value=hullblockextremes.db,
                                          index=this, 
                                          weight = 1./E1_0)
                FPD.add_yVertexConstraint(kind='max',
                                          value=hullblockextremes.de*1.5,
                                          index=this, 
                                          weight = 1./E1_0)
                #--------------------------------------z
                FPD.add_zVertexConstraint(kind='min',
                                          value=hullblockextremes.lb,
                                          index=this, 
                                          weight = 1./E1_0)
                FPD.add_zVertexConstraint(kind='max',
                                          value=hullblockextremes.le*1.5,
                                          index=this, 
                                          weight = 1./E1_0)
                #--------------------------------------done
        #
        #**************************************************************************
        # Fairness
        FPD.add_E1(kind='LS', weight = 1.5/E1_0)
        FPD.add_E2(kind='LS', weight = .5/E2_0)
        FPD.add_E3(kind='LS', weight = .05/E3_0)
        FPD.add_ArcLengthApprox(kind='LS', weight = 5.5/S_0)
        #
        #
        #**************************************************************************
        #
        interval_data, small = interval_bounds(curve)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                            interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, 
                                         data = interval_data)
        #
        #**************************************************************************
        # Optimize
        print 'lcurve fixing solve:'
        Lspline.optimize(stop=30)
        return Lspline
    
    
    
    def nose_constraint_form(Lspline):
        """ leading edge of the bare hull shape 
        expressed in the longitudinal curve list
        
        Actually, how about solving the nose constraint
        in the THB refined space only!
        """
        message = 'ERROR in utility_optimization: send Bmax to: \n function <interpolate_vertices>'
        assert(Bmax is not None),message
        #
        norms = package_norms(Lspline)
        E1_0 = norms[0]
        E2_0 = norms[1]
        E3_0 = norms[2]
        S_0  = norms[3]
        curve = Lspline.curve
        #
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve)
        #
        #**************************************************************************
        # make the nose:
        #
        #        FPD.add_xVertexConstraint(kind='equality',
        #                                 index = 1,
        #                                 value = 0.)
        #
        #        FPD.add_relative_xVertexConstraint(
        #                                 kind='equality',
        #                                 location=None,
        #                                 index = 0,
        #                                 index2 = 1,
        #                                 value = 0.,
        #                                 weight=1.)
        FPD.add_xFixity(index=1,
                        value=curve.vertices[0,0])
        #
        FPD.add_yVertexConstraint(kind='min',
                                 index = 1,
                                 value = 0.)
        #        FPD.add_yVertexConstraint(kind='LS',
        #                                 index = 1,
        #                                 value = Bmax)
        #        FPD.add_yVertexConstraint(kind='LS',
        #                                 index = 2,
        #                                 value = Bmax)
        #
        #**********************************************************************
        # interoplation constraints
        for pt,loc in zip(interpolate_vertices[1:-1],ukbar[1:-1]):
            #print 'pt = ', pt
            FPD.add_xPointConstraint(kind= 'equality', 
                                     value=pt[0], 
                                     location=loc )
            FPD.add_yPointConstraint(kind= 'equality', 
                                     value=pt[1], 
                                     location=loc )
            FPD.add_zPointConstraint(kind= 'equality', 
                                     value=pt[2], 
                                     location=loc )   
        #
        #**********************************************************************
        # block extrema?
        #
        # be careful here:
        #  the extrema of the hull are not necessarily the max extents
        #  that the lcurve vertices should be allowed to excurt to
        #  -maybe better to change this to max extents for 
        #   __actual_points_on_the_curves_at_even_intervals__
        #  -for now we just give it a little lee-way
        if hullblockextremes is not None:
            for i,pt in enumerate(curve.vertices[1:-1]):
                this = i+1
                #--------------------------------------x
                FPD.add_xVertexConstraint(kind='min',
                                          weight=.5,
                                          value=hullblockextremes.bb,
                                          index=this)
                FPD.add_xVertexConstraint(kind='max',
                                          weight=.5,
                                          value=hullblockextremes.be*1.5,
                                          index=this)
                #--------------------------------------y
                FPD.add_yVertexConstraint(kind='min',
                                          weight=.5,
                                          value=hullblockextremes.db,
                                          index=this)
                FPD.add_yVertexConstraint(kind='max',
                                          weight=.5,
                                          value=hullblockextremes.de*1.5,
                                          index=this)
                #--------------------------------------z
                FPD.add_zVertexConstraint(kind='min',
                                          weight=.5,
                                          value=hullblockextremes.lb,
                                          index=this)
                FPD.add_zVertexConstraint(kind='max',
                                          weight=.5,
                                          value=hullblockextremes.le*1.5,
                                          index=this)
                #--------------------------------------done
        
        #
        #**************************************************************************
        # Fairness
        FPD.add_E1(kind='LS', weight = .5/E1_0)
        FPD.add_E2(kind='LS', weight = .5/E2_0)
        FPD.add_E3(kind='LS', weight = .5/E3_0)
        FPD.add_ArcLengthApprox(kind='LS', weight = 5./S_0)
        #
        #**************************************************************************
        #
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                            interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, 
                                         data = interval_data)
        #
        #**************************************************************************
        # Optimize
        print 'lcurve stage 2 solve (nose constraint):'
        Lspline.optimize(stop=25)
        #
        #**************************************************************************
        # Done
        return Lspline
    
    
    def getckit(Lspline):
        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
            ckit =False
        else:
            ckit = True 
        return ckit
        
    
    #
    #**********************************************************************
    # initialize
    ckit = False
    #
    #**********************************************************************
    # solve Lcurve stage 1
    Lspline = initial_form(curve)
    Lsave1 = copy.deepcopy(Lspline)
    Lsave_ini = copy.deepcopy(Lspline)
    #
    #**********************************************************************
    # check Lcurve stage 1
    ckit = getckit(Lspline)
    #
    #**********************************************************************
    # Lcurve nose constraint  (LEAVE for THB!)
    if nose_constraint:
        Lspline = nose_constraint_form(Lspline)
        if Lspline.conv > Lspline.tol: #TODO: fix need for this due to stupid delta_z = NAN bail out:
            Lspline.error_code = 'max_iter'
            
            
        if Lspline.error_code is 'NAN':
            if ckit:
                Lspline = Lsave1
                Lspline = simplified_form(Lspline)
                ckit = getckit(Lspline)
                
                
            else: #initial Lspline is bad
                Lspline = Lsave1
                
                
        elif Lspline.error_code is 'max_iter':
            if Lspline.conv>.01:
                if ckit:
                    Lspline = Lsave1
                    Lspline = simplified_form(Lspline)
                    ckit = getckit(Lspline)
                    if not ckit:
                        Lspline = Lsave1
                    else:
                        pass #done!
                else: #initial Lspline is bad
                    Lspline = initial_form(curve) #give up
            else:
                pass #take Lspline from solution 2!
                
    #
    #**********************************************************************
    #  
    #
    #---
    #TODO: only put the nose on THB fine level(s)!
            
            
    
    print 'Lspline interpolated curve done'
            
    if return_Lspline:
        return Lspline
    else:
        return Lspline.curve



##
##*****************************************************************************
##
##
##*****************************************************************************
##
def fast_interpolate_vertices(vertices_list,
                         helper_curve_list,
                         dk = 4,
                         dnump = 30,
                         return_Lspline=False,
                         make_thb = False):
    """
        input:  a curve we'd like to re-parameterize
                (probably with a problamatic knot vector)
            
        output: a curve  
                with a standardized knot vector
        
        this is a speciallization of the way
        to match any curve of equal dimension.

        this is an unfinished start at making a fast
        version of the above function interpolate_vertices

        make it fast method:  short circuit FPD and Lspline
        instantiation
    """
    if verbose:print 'utility function: fast_interpolate_vertices'
    outputlist = []
    #start = True #need to make startup faster
    for initial_vertices, helper_curve in zip(vertices_list,
                                              helper_curve_list):        
        
        ukbar = helper_curve.ukbar
        if make_thb:
            curve = rbspline(helper_curve.vertices,
                             k=dk,
                             nump=dnump)
        else:
            curve = spline.Bspline(helper_curve.vertices,
                                   k=dk,
                                   nump=dnump)
        
        
        interval_data, small = interval_bounds(curve)
        #if start: #need to make startup faster!
        FPD = FormParameterDict(curve)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', 
                                weight = .5)
        for pt,loc in zip(initial_vertices[1:-1],ukbar[1:-1]):
            
            FPD.add_xPointConstraint(kind= 'equality', 
                                     value=pt[0], 
                                     location=loc )
            FPD.add_yPointConstraint(kind= 'equality', 
                                     value=pt[1], 
                                     location=loc )
            FPD.add_zPointConstraint(kind= 'equality', 
                                     value=pt[2], 
                                     location=loc )
            
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, 
                            interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, 
                                         data = interval_data)
        
        Lspline.curve.pts_M_pts()
        Lspline.curve.compute_arclength()
        
        Lspline.optimize(stop=30)
    
    if return_Lspline:
        return Lspline
    else:
        return Lspline.curve


##
##*****************************************************************************
##
##
##*****************************************************************************
##
   
def lift_curve(curve,
               elevation,
               index,
               use_knots=False):
    """
    """
    if verbose:print 'utility function: lift_curve'
    array = curve.vertices
    size,dim = np.shape(array) 
    new_array = np.zeros((size,dim+1),float)
    indices = [ el for el in range(dim+1)]
    indices.pop(index)
    new_array[:,indices] = array[:]
    new_array[:,index] = elevation
    if use_knots:
        return spline.Bspline(vertices = new_array,
                              k     = curve.k,
                              nump  = curve.nump,
                              t     = curve.t)
    else:
        return spline.Bspline(vertices  = new_array,
                              k         = curve.k,
                              nump      = curve.nump)




##
##*****************************************************************************
##
##
##*****************************************************************************
##
def duplicate_shift_curve(curve,
                            shift,
                            index,
                            use_knots=False):
    if verbose:print 'utility function: duplicate_shift_curve'
    new_array = copy.copy(curve.vertices)
    new_array[:,index] = new_array[:,index] + shift
    
    if use_knots:
        return spline.Bspline(vertices  = new_array,
                              k         = curve.k,
                              nump      = curve.nump,
                              t         = curve.t)
    else:
        return spline.Bspline(vertices  = new_array,
                              k         = curve.k,
                              nump      = curve.nump)
        

##
##*****************************************************************************
##
## intend to factor utility functions out of various 
##  user files and into this utility file.
##
##*****************************************************************************
##

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

def arc_vertices(start, end, num):

    start = np.asarray(start)
    end = np.asarray(end)
    D = end - start
    radius = 0.5*D[1]
    rs = radius**2
    dim = len(D)
    vertices = []
    for i in range(dim):
        vi = np.linspace(start[i],end[i],num)
        vertices.append(vi)
    vertices = np.asarray(vertices)
    vertices = np.transpose(vertices)
    for i in range(len(vertices)):
        vertices[i,0] = np.sqrt(rs - (vertices[i,1]-radius)**2)
    return vertices


def interpolate(rise,run,step,ini):
    return ini + step*(rise/run)

def interpolate_vertical(rise,run,step,ini):
    return ini + step*(rise/run)

def push_vertices(vertices,
                  bottom_vertex,bottom_index,
                  top_vertex,top_index):
    return
    

##
##*****************************************************************************
##
##
##
##
##
##   TESTING AND CONVERGENCE CHECKING, etc...
##
##
##
##
##*****************************************************************************
def THBspline_to_Bspline_ckdiff_evaluator(hull,mod_level = None):
    """check the difference between B-spline and THBspline
    suface evaluation for a given THBsurface called hull
    """
    
    print 'Evaluating Fit between THBsurface and its Bspline counterpart'
    #
    P = hull.project_all_to_finest()
    #
    if mod_level is not None:
        surf = hull.get_surface_by_level(mod_level)
        print 'comparing the THB surface to the Bspline at default level',surf.level
        print 'change level with "mod_level" option'
    else:
        surf = hull
        print 'comparing the THB surface to the Bspline at level',surf.level
    print 'This comparison will not match if changes have been made at another level'
    print 'Because the B-spline will have no way to match it'
    
    nn = len(hull.mastercurve_v.s)
    findmax00 = np.zeros((nn,nn),float)
    findmax01 = np.zeros((nn,nn),float)
    for i,u in enumerate(hull.mastercurve_v.s):
        for j,v in enumerate(hull.mastercurve_v.s):
            thb = hull.eval(u,v)[0,0]
            thb2 = hull.eval_smart(u,v,P)[0,0]
            bspline = surf.surface_point(u,v)
            ck = np.linalg.norm(bspline)
            if ck<.001:
                print 'skipping ', thb, bspline, thb-bspline
            else:
                findmax00[i,j] = np.linalg.norm(thb - bspline)#/ck
                findmax01[i,j] = np.linalg.norm(thb2 - bspline)#/ck
    print 'max +/- difference between THB and B-spline evaluation = ', findmax00.max(), findmax00.min()
    print 'max +/- difference between THB and B-spline evaluation = ', findmax01.max(), findmax01.min()
    ui = np.argmax(findmax01.max(axis=0))
    vi = np.argmax(findmax01.max(axis=1))
    print 'indicies of max = {},{}'.format(np.argmax(findmax01.max(axis=0)),
                                           np.argmax(findmax01.max(axis=1)))
    print 'indicies of min = {},{}'.format(ui,vi)
    print 'u = {}, v = {}'.format(hull.u[ui], hull.v[vi])
    
    
    
    findmax10 = np.zeros((nn,nn),float)
    findmax11 = np.zeros((nn,nn),float)
    for i,u in enumerate(hull.mastercurve_v.s):
        for j,v in enumerate(hull.mastercurve_v.s):
            thb = hull.eval_old(u,v)[0,0]
            thb2 = hull.eval_smart(u,v,P)[0,0]
            bspline = surf.surface_point(u,v)
            ck = np.linalg.norm(bspline)
            if ck<.001:
                print 'skipping ', thb, bspline, thb-bspline
            else:
                findmax10[i,j] = np.linalg.norm(thb - bspline)#/ck
                findmax11[i,j] = np.linalg.norm(thb2 - bspline)#/ck
    ui = np.argmax(findmax11.max(axis=0))
    vi = np.argmax(findmax11.max(axis=1))
    print 'max +/- difference between THB and B-spline evaluation = ', findmax10.max(), findmax10.min()
    print 'max +/- difference between THB and B-spline evaluation = ', findmax11.max(), findmax11.min()
    print 'indicies of max = {},{}'.format(np.argmax(findmax11.max(axis=0)),
                                           np.argmax(findmax11.max(axis=1)))
    print 'indicies of min = {},{}'.format(ui,vi)
    print 'u = {}, v = {}'.format(hull.u[ui], hull.v[vi])
  
    return

def THBspline_evaluator_diffs(surf):
    """check the difference between THBspline
    suface evaluation for various eval methods
    """
    print 'Evaluating Fit between THBsurface eval and eval_smart'
    P = surf.project_all_to_finest()
    #
    hull = surf
    nn = len(hull.mastercurve_v.s)
    findmax00 = np.zeros((nn,nn),float)
    for i,u in enumerate(hull.mastercurve_v.s):
        for j,v in enumerate(hull.mastercurve_v.s):
            thb = hull.eval(u,v)[0,0]
            thb2 = hull.eval_smart(u,v,P)[0,0]
            #bspline = hull.surface_point(u,v)
            #ck = np.linalg.norm(thb)
            #if ck<.001:
            #    print 'skipping ', thb, bspline, thb-bspline
            #else:
            #    findmax00[i,j] = np.linalg.norm(thb - thb2)#/ck
            findmax00[i,j] = np.linalg.norm(thb - thb2)#/ck
    print 'max +/- difference between THB and B-spline evaluation = ', findmax00.max(), findmax00.min()
    #print 'max +/- difference between THB and B-spline evaluation = ', findmax01.max(), findmax01.min()
    ui = np.argmax(findmax00.max(axis=0))
    vi = np.argmax(findmax00.max(axis=1))
    print 'indicies of max = {},{}'.format(np.argmax(findmax00.max(axis=0)),
                                           np.argmax(findmax00.max(axis=1)))
    print 'indicies of min = {},{}'.format(ui,vi)
    print 'u = {}, v = {}'.format(hull.u[ui], hull.v[vi])
    
    
    
    findmax11 = np.zeros((nn,nn),float)
    for i,u in enumerate(hull.mastercurve_v.s):
        for j,v in enumerate(hull.mastercurve_v.s):
            thb = hull.eval_old(u,v)[0,0]
            thb2 = hull.eval_smart(u,v,P)[0,0]
            #bspline = surf.surface_point(u,v)
            #ck = np.linalg.norm(thb)
            #if ck<.001:
            #    print 'skipping ', thb, bspline, thb-bspline
            #else:
            #    findmax11[i,j] = np.linalg.norm(thb - thb2)#/ck
            findmax11[i,j] = np.linalg.norm(thb - thb2)#/ck
    ui = np.argmax(findmax11.max(axis=0))
    vi = np.argmax(findmax11.max(axis=1))
    #print 'max +/- difference between THB and B-spline evaluation = ', findmax10.max(), findmax10.min()
    print 'max +/- difference between THB and B-spline evaluation = ', findmax11.max(), findmax11.min()
    print 'indicies of max = {},{}'.format(np.argmax(findmax11.max(axis=0)),
                                           np.argmax(findmax11.max(axis=1)))
    print 'indicies of min = {},{}'.format(ui,vi)
    print 'u = {}, v = {}'.format(hull.u[ui], hull.v[vi])
    return
##
##*****************************************************************************
##
##
##
##
##   Funky real valued array logic for special purposes
##
##
##
##*****************************************************************************
def data_range(outer,inner, nump=30):
    """Not used
    --disjoint version used instead--
    parameters:
        outer = ia(x_inf,x_sup) interval
        inner = ia(y_inf,y_sup) interval
        such that
        {x_inf, y_inf, y_sup, x_sup}
        is well ordered and inner is inside outer
        
    returns:  
        numpy array parameterizing 
        ui = [outer_:inner_]
        ue = [inner^:outer^]
        
    Notes:
        This was added late in 2012..
        used in thbsurface plotcurves method
        to tease out the portions of the parameter space
        belonging to the current and next level
    """
    assert(nump%2 == 0),'Must use even nump for plotting.  \n'+\
                        'Note also that disjoint inner level bounds \n'+\
                        'are not supported in plotting!'
    outer = outer
    inner = inner
    sui = set(np.linspace(outer[0],inner[0],nump/2,endpoint=True))
    sue = set(np.linspace(inner[1],outer[1],nump/2,endpoint=True))
    goodsui = len(sui) > 1 
    goodsue = len(sue) >1
    
    if goodsui and goodsue:
        su = np.asarray(list(sui)+list(sue))
    elif goodsui:
        su = np.asarray(list(sui))
    elif goodsue:
        su = np.asarray(list(sue))
    return su


def data_range_disjoint_vectors(outer,inner, nump=30):
    """Essentially Deprecated in favor of new dijoint methods
    in extended_interval_arithmetic.py
    and THBsurface
    
    parameters:
        outer = ia(x_inf,x_sup) a simple ia type interval
        inner = ia(y_inf,y_sup) a simple ia type interval
        such that
        {x_inf, y_inf, y_sup, x_sup}
        is well ordered and inner is inside outer
        
    returns:  
        numpy array parameterizing 
        ui = [outer_:inner_]
        ue = [inner^:outer^]
    
    Notes:
        This was added late in 2012..
        used in thbsurface plotcurves method
        to tease out the portions of the parameter space
        belonging to the current and next level
    """
    assert(nump%2 == 0),'Must use even nump for plotting.  \n'+\
                        'Note also that disjoint inner level bounds \n'+\
                        'are not supported in plotting!'
    #outer = outer
    #inner = inner
    ui = np.linspace(outer[0],inner[0],nump/2,endpoint=True)
    ue = np.linspace(inner[1],outer[1],nump/2,endpoint=True)
    sui = set(ui)
    sue = set(ue)
    goodsui = len(sui) > 1 
    goodsue = len(sue) >1
    
    if goodsui and goodsue:
        return  ui,ue
    elif goodsui:
        return  ui, np.asarray([])
    elif goodsue:
        return ue, np.asarray([])
    else:
        return np.asarray([]),np.asarray([])
    

    
def package_norms(Lspline, recalc=True):
    if recalc:
        Lspline.curve.pts_M_pts()
        Lspline.curve.compute_arclength()
    try:
        E1 = copy.copy(Lspline.curve.E1.value)
    except:
        E1 = copy.copy(Lspline.curve.E1)
    try:
        E2 = copy.copy(Lspline.curve.E2.value)
    except:
        E2 = copy.copy(Lspline.curve.E2)
    try:
        E3 = copy.copy(Lspline.curve.E3.value)
    except:
        E3 = copy.copy(Lspline.curve.E3)
    
    if Lspline.curve.AL is not None:
        try:
            S = copy.copy(Lspline.curve.AL.value)
        except:
            S = copy.copy(Lspline.curve.AL)
    elif Lspline.curve.ALapprox is not None:
        try:
            S = copy.copy(Lspline.curve.ALapprox.value)
        except:
            S = copy.copy(Lspline.curve.ALapprox)
    else:
        S = 1.0
    norms = [E1,E2,E3,S]
    return norms





#def issue_longitudinals_stepwise(tcurvenet,
#                                 k=4,nump=30):
#    """THB Optimize longitudinal curves
#    
#    #
#    self = thbcurve
#    #
#    
#    dev
#    ----------
#    tvertnet = make_vertex_matrix(tcurvenet)
#    """
#    #
#    tvertnet = make_vertex_matrix(tcurvenet)
#    lcurvenet = []
#    curve = tcurvenet[0]
#    N = curve.n
#    nnew = N-2 #just optimize the interior THB curves?
#    #
#    #
#    i=1 #dev
#    #
#    #
#    for i in range(1      , 2): #      ,   nnew+1): #
#        #**************************************
#        #linear interpolation
#        #
#        #vanilla_spline = spline.interpolatedBspline(
#        #                tvertnet[:,i], k, nump)
#        #**************************************
#        # simple 4 c.v. starting linear curve
#        vtinterp = tvertnet[:,i] 
#        
#        total_len = len(vtinterp)
#        ti = 4
#        #diff = total_len-ti
#        #----------------------------------
#        # -initial vertices and intial curve:
#        lvts = linear_vertices(vtinterp[0],
#                               vtinterp[ti],4)
#        thbcurve = rbspline(lvts,k,nump)
#        #
#        Lspline = THBsolver_stepwise(thbcurve,
#                            vtinterp[0:ti])
#        #----------------------------------
#        for j in range(ti+1,total_len):
#            #
#            Lspline = Lspline.get_finest_Lspline()
#            curve = Lspline.curve
#            vL = list(curve.vertices)
#            vL.append(vtinterp[j])
#            new_vertices = np.asarray(vL)
#            if len(new_vertices)>len(curve.vertices):
#                curve = curve.dyadic_refinement()
#                curve.vertices[-1] = vtinterp[j]
#                #curve.parent=None
#            #
#            Lspline = THBsolver(curve,
#                                vtinterp[0:j],
#                                maxlevel=4,
#                                normalize=False)
#        
#        
#        print 'finished new Lspline'
#        lcurvenet.append(Lspline)
#    return lcurvenet
    



#def THBsolver_stepwise(thbcurve=None,
#                      vtinterp=None,
#                      Lspline=None,
#                      refinemax=0,
#                      normalize=True):
#    """Lspline solver (adaptive) for THB splines
#    
#    intention : demonstrate an adaptive solver
#    for THB-spline optimization of a longitudinal 
#    hull curve!
#    
#    
#    inputs
#    ----------
#        thbcurve    = starting base level THB curve (n=4)
#        vinterp     = the transverse vertices to be interpolated
#                                        (including the endpoints!)
#    process
#    ----------
#        1.  solve one iteration on the base level
#        2.  check the gradient
#        
#    notes
#    ----------
#        11 control vertices and 6 interpolation vertices
#        makes for 51 elements of the gradient
#        11*3 + 6*3 = 51
#    """
#    print '   uopt.THBsolver'
#    num_refined = 0
#    #because we are using equality constraints
#    # we must begin with a soluble number of vertices
#    # to get any information out of our solver!
#    #**************************************************
#    # refine until the ini curve is feasible 
#    # given the number of vertices to 
#    # be interpolated
#    #
#    #**************************************************
#    #
#    if vtinterp is not None and Lspline is None:
#        if thbcurve.n<11:
#            while nthb_vs_nvi(thbcurve,vtinterp) and thbcurve.n<11:
#                print 'pre-refinement'
#                thbcurve = thbcurve.dyadic_refinement()
#        #else: we may get a singular curve until enough vertices 
#        # are there to match the vertices to be interpolated
#        #**************************************************
#        thbcurve.parent = None
#        #
#        print '   setup'
#        Lspline = setup_interpolation_solver(thbcurve,
#                                             vtinterp)
#        #******************************************
#    #    
#    #******************************************
#    this = Lspline.get_finest_Lspline()
#    if this is not Lspline:
#        num_refined +=1
#    #******************************************************
#    cmax = 2 #max num times to iterate
#    count =0
#    setconv = False
#    while this is not None and count<cmax:
#        print '   optimize'
#        #**************************************************
#        print '   optimize eq constraints'
#        this.verbose = True
#        """(1) makes THB hierarchical opti work! (1 of 2)"""
#        save_vertices = copy.copy(this.curve.vertices)
#        #*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
#        this.tol = 1.e-7
#        this.optimize_thb(stop=15,
#                          ceiling = 1.e-10,
#                          relax_grad_bounds=False)
#        this.tol = 1.e-7
#        #**************************************************
#        #normalize
#        if normalize: this.normalize()
#        #**************************************************
#        ckit = ckLspline(this)
#        if ckit and this.conv<this.tol:
#            setconv = True
#            break
#        elif not ckit:
#            print 'restore vertices'
#            """(2) makes THB hierarchical opti work! (2 of 2)"""
#            this.curve.vertices = save_vertices
#            #*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
#        #**************************************************
#        # manually call refine:
#        #**************************************************
#        print 'again = ',this.refine_by_these_grad_indices
#        if num_refined < refinemax:
#            print 'Not converged, refining ...'
#            this.refine()
#            newlevel = this.get_finest_Lspline()
#            if newlevel is not this:
#                num_refined +=1
#            this = newlevel
#        #**************************************************
#        count +=1
#
#    
#    #******************************************************
#    # retore equality constraints
#    print 'restore LS to equality'
#    this.restore_Lagrangian()
#    #******************************************************
#    #print '   optimize'
#    #Lspline.optimize_thb(stop=5,
#    #                     relax_grad_bounds=True)
#    print 'Done with this longitudinal'
#    if setconv:
#        print 'convergence is true'
#    else:
#        print 'convergence is false'
#    return this.get_coarsest_Lspline()