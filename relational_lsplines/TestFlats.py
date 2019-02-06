#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 21:56:50 2018

@author: luke

aft ship hull needs vertical flat near DWL

fwd ship hull needs horizontal flat in a sort of triangle 
from midship to bulb
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




def init_lagrangian(curve,weight = 1.):
    """
    """
    FPD = FormParameterDict(curve)
    #
    FPD.add_E1(kind='LS', weight = weight)
    FPD.add_E2(kind='LS', weight = weight)
    FPD.add_E3(kind='LS', weight = weight)
    FPD.add_ArcLengthApprox(kind='LS',
                          weight = weight)
    return FPD

def fini_lagrangian(FPD,
                    this_kind='equality'):
    curve = FPD.curve
    interval_data, small = interval_bounds(curve)
    L = Lagrangian(FPD)
    interval_data, small = lagrangian_bounds(L, 
                        interval_data, small, 1.e4)
    Lspline = IntervalLagrangeSpline(curve, L, 
                                     data = interval_data)
    return Lspline



if __name__ == """__main__""":
    
    k=4
    nump = 30
    #
    ini = [0.,0.]
    intersect = [6.,0.]
    fini = [12.,12.]
    #**************************************************************************
    # first version:  flat to intersection then flat up
    # what is the area?
    #
    v1 = linear_vertices(ini,
                         intersect,
                         3)
    v2 =  linear_vertices(intersect,
                         fini,
                         3)
    vertices1 = list(v1)+list(np.asarray([[6.,0.]]))+list(v2)
    c1 = rbspline(vertices1,
                  k = k,
                  nump = nump)
    #c2 = c1.dyadic_refinement()
    #c2.bounds = BoundsList( [ ia(.3,.7) ] )
    
    ax = c1.plotcurvehierarchy()
    
    
    FPD = init_lagrangian(c1)
    
    FPD.add_Area_to_y_Constraint(110.)
    
    
    #second version, start at the intersection and go flat to top
    # what is the area?
    vertices01 = linear_vertices(ini,
                                 intersect,
                                 7)
    vertices02 = linear_vertices(intersect,
                               fini,
                               7)
    c2 = rbspline(vertices01,
                  k = k,
                  nump = nump)
    c3 = rbspline(vertices02,
                  k = k,
                  nump = nump)
    
    #ax = c2.plotcurvehierarchy(canvas = ax)
    ax = c3.plotcurvehierarchy(canvas = ax)
    
    print '\n C1:'
    c1.compute_area()
    print c1.area
    c1.compute_area_to_y()
    print c1.area_to_y
    
    print '\n C2:'
    c2.compute_area()
    print c2.area
    c2.compute_area_to_y()
    print c2.area_to_y

    print '\n C3:'
    c3.compute_area()
    print c3.area
    c3.compute_area_to_y()
    print c3.area_to_y
    
    
    
     
    #ini = [6.,0.]
    #fini = [12.,12.]
    vertices2 = linear_vertices(ini,
                               fini,
                               7)
    c00 = rbspline(vertices2,
                  k = k,
                  nump = nump)
    print '\n C00:'
    c00.compute_area()
    print c00.area
    c00.compute_area_to_y()
    print c00.area_to_y