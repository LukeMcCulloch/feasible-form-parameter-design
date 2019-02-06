#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 17:01:40 2017

@author: luke
"""
import numpy as np
import curve             as spline
from   ADILS             import IntervalLagrangeSpline, Lagrangian
from   FormParameter     import FormParameterDict
from   initialValues     import interval_bounds, lagrangian_bounds

import copy   

if __name__ == '__main__':
    
##    vertices=np.asarray([[0.,0.,0.],
##                         [-6.4,0.,0.],
##                         [-9.,0.,2.],
##                         [-6.4,0.,4.],
##                         [-10.,0.,10.]])


    points = np.array([[ 11.,  6.,   1.],
                       [ 30.,  16.,   2.5],
                        [5.,5.,.5]])
    

    vertices = np.array([[0.,0.,0.],
                        [10.,5.,1.],
                        [15.,7.5,1.5],
                        [20.,10.,2.],
                        [20.,10.,2.],
                        [20.,10.,2.],
                        [25.,15.,2.5],
                        [25.,15.,2.5],
                        [25.,15.,2.5],
                        [30.,15.,3.],
                        [40.,20.,4.]])
    
    vertices = np.array([[0.,0.,0.],
                        [10.,5.,1.],
                        [15.,7.5,1.5],
                        [30.,15.,3.],
                        [40.,20.,4.]])
    
    k=4
    nump=50
    curve = spline.interpolatedBspline(vertices, k, nump)
    cc = copy.deepcopy(curve)
    
    interval_data, small = interval_bounds(curve)
    
    FPD = FormParameterDict(curve) 
    FPD.add_E1(kind='LS', weight = .5)
    FPD.add_E2(kind='LS', weight = .5)
    FPD.add_E3(kind='LS', weight = .5)
    FPD.add_ArcLengthApprox(kind='LS', weight = .5)
    
    for i in range(len(curve.vertices[1:-1])):
        j = i+1
        """TODO, TLM OCT 30 2017:
            find where the curve interpolated pts are : 
                done:  curve.ukbar[i]
            and use these as the targets for the
            optimization z_point_constraints!
        
            Done: NOV 1, 2017
        """
        #"""
        thisloc = curve.ukbar[j]
        #thisvertex = curve.FindPointOmnidirectional(findx=vertices[j,2],2)
        
        FPD.add_xPointConstraint(kind='equality',
                                 value = vertices[j,0],
                                 location=thisloc)
        FPD.add_yPointConstraint(kind='equality',
                                 value = vertices[j,1],
                                 location=thisloc)
        FPD.add_zPointConstraint(kind='equality',
                                 value = vertices[j,2],
                                 location=thisloc)
        #"""
    
    L = Lagrangian(FPD)
    interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
    Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
    Lspline.optimize(stop = 50)
    
    
    print curve.CurvePoint(curve.ukbar[1])
    print cc.CurvePoint(curve.ukbar[1])
    print curve.CurvePoint(curve.ukbar[2])
    print cc.CurvePoint(curve.ukbar[2])