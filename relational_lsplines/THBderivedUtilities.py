#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 09:51:05 2018

@author: luke

Stuff that uses THB surfaces and curves to ... do things...
"""


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




level2n = { 0:4,
            1:5,
            2:7,
            3:11,
            4:19,
            5:35}
n2level = { 4:0,
            5:1,
            7:2,
            11:3,
            19:4,
            35:5}

class DyadicMap(object):
    def __init__(self):
        self.level2n = {0:4,
                        1:5,
                        2:7,
                        3:11,
                        4:19,
                        5:35}
        self.n2level = {4:0,
                        5:1,
                        7:2,
                        11:3,
                        19:4,
                        35:5}
    
def make_surface(lcurvenet,
                 hullcurvenet=None,
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
    
    """
    ul = n2level[len(lcurvenet)]
    vl = n2level[lcurvenet[0].n]
    
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
                              hullcurvelist=hullcurvenet,
                              ulevel = ul,
                              vlevel = vl)
    
    if refine:
        s1 = s0.dyadic_loft_refinement()
    s0.compute_surface()
    return s0