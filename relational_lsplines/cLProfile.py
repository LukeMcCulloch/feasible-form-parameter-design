#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 14 17:25:33 2018

@author: luke
"""



import numpy as np
from scipy import ndimage
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.axislines import SubplotZero
import matplotlib.patches as patches
from matplotlib.patches import Ellipse
import matplotlib.path    as path
mpl.rcParams['legend.handlelength'] = 0

import curve as spline





def compute_cLProfile(LocMaxSection=None, 
                      Xc=None,
                      flat=False, 
                      height=None, 
                      alphab=0.,
                      alphae=0.,
                      drop_aft=.25):
    
    elevation = MaxDraft
    k = 4
    nump = 30
    
    xb      = 0.
    yb      = 0.
    xe      = lwl 
    ye      = 0.+MaxDraft*drop_aft
    
    Cab_given = 0.
    Cae_given = 0.
    slope = 'down'
    
    Dmax = MaxDraft
    area = CLarea

    ab = alphab
    ae = alphae

    curve = initial_curve((xb,yb),
                          (xe,ye),
                          num=11, k=k, nump=nump)
    
    
    v = copy.deepcopy(curve.vertices)
    v[4:8,1] = Dmax
    v[3,1] = Dmax*2./3.
    v[8,1] = Dmax*2./3.
    
    
    curve.vertices = v 
    
    
    
    def OSV_bbow(curve):
        """Lspline solver 
        
            for the 
            
            CProfile keel curve
            
         (intended for a transom sterned hull)
        with bulbous bow
        
        dev
        ----------
        
        elevation = self.MaxDraft
        k           = self.k
        nump        = self.nump
        xb          = 0.
        yb          = 0.
        xe          = self.LengthDWL
        ye          = 0.+self.MaxDraft*drop_aft
        #alphab     = 85.##35
        #alphae     = -85.##-35
        Cab_given   = 0.
        Cae_given   = 0.
        slope       = 'down'

        Dmax = self.MaxDraft
        area = self.CLarea

        ab = alphab
        ae = alphae
        """
        print 'OSV_bbow solver Oct 14, 2018'
#            #
#            #**********************************************************************
#            # get weights
#            if Lspline is not None:
#                norms = package_norms(Lspline)
#                E1_0 = norms[0]
#                E2_0 = norms[1]
#                E3_0 = norms[2]
#                S_0  = norms[3]
#                curve = Lspline.curve
#            else:
#                E1_0 = 1.
#                E2_0 = 1.
#                E3_0 = 1.
#                S_0  = 1.
        # CPK
        #**********************************************************************
        # initialize
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve)
        #
        #**********************************************************************
        # area  CPKeelProfile stage 1
        FPD.add_AreaConstraint(kind='equality',
                               value = area)
        #
        #**********************************************************************
        # Xc CPKeelProfile
        if Xc is None:
            print 'CPK OSV bbow, equality LCG'
            FPD.add_XcConstraint(kind='LS', 
                                 value=LCG*.75,
                                 weight = 100.) #cLProfile flat to bulb...
            #FPD.add_XcConstraint(kind='max', 
            #                     value=self.LCG*.95)
            #FPD.add_XcConstraint(kind='min', 
            #                     value=self.LCG*0.8)
        else:
            print 'Centerplane setting Xc constraint to ',Xc
            FPD.add_XcConstraint(kind='LS', 
                                 value=Xc)
        #                                     value = 0.)
        #
        ##**************************************************************
        # CPkeelProfile  
        #print 'Adding x - fixity!'
        #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # TLM, straight flat stem
#            FPD.add_xFixity(index=1,
#                            value=xb)
#            FPD.add_xFixity(index=2,
#                            value=xb)
        #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # Dr. Birk, flat stem
#            FPD.add_AngleConstraint(
#                                    kind        = 'equality',
#                                    location    = 0.,
#                                    value       = 75.)
        FPD.add_CurvatureConstraint(
                                kind        = 'equality',
                                location    = 0.,
                                value       = 0.,
                                weight      = 100.)
        #--------------------------------------
        # +++
        #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # TLM, straight flat stem
        #FPD.add_xVertexConstraint(kind='LS',
        #                         index = 3,
        #                         value = 0.)
        FPD.add_xFixity(index=3, #5,
                        value=0.)
        FPD.add_relative_xVertexConstraint(
                                    kind = 'equality',
                                    value = 0.,
                                    index = 3,
                                    index2 = 4,
                                    seperation=0.)
        #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        #--------------------------------------
        # ++
#            FPD.add_xVertexConstraint(kind='LS',
#                                     index = 8,
#                                     value = xe)
#            FPD.add_xVertexConstraint(kind='LS',
#                                     index = 9,
#                                     value = xe)
        #--------------------------------------
        ##**************************************************************
        # CPkeelProfile  
        #print 'Adding x - fixity!'
        #tlmOct 14 2018 tweak this:
        FPD.add_xFixity(index=5, #5,
                        value=FOCP0)
        FPD.add_xFixity(index=6, #6,
                        value=FOCP1)
        #
        #------------------------------- #tlmOct 14 2018 tweak this:
        #print 'Adding y - fixity!'
        #FPD.add_yFixity(index=3,
        #                value=Dmax)
        FPD.add_yFixity(index=4,
                        value=Dmax)
        FPD.add_yFixity(index=5,
                        value=Dmax)
        FPD.add_yFixity(index=6,
                        value=Dmax)
        FPD.add_yFixity(index=7,
                        value=Dmax)
        #-------------------------------
        FPD.add_yFixity(index=8,
                        value=Dmax)
        #FPD.add_yFixity(index=9,
        #                value=Dmax)
        #?-------------------------------?
        #            FPD.add_yVertexConstraint(kind='max',
        #                                     index = 8,
        #                                     value = Dmax)
        #            FPD.add_yVertexConstraint(kind='max',
        #                                     index = 9,
        #                                     value = Dmax)
        #?-------------------------------?
        #
        #
        #**********************************************************************
        # CPkeelProfile Fairness 
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .001)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.5)
        #
        #**********************************************************************
        #Setup and Solve:
        #
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        print 'Doing first stage CPkeel solve (include bbow weighting)'
        Lspline.optimize(stop=50)
        return Lspline
    
    
    
    def OSV_bbow2(Lspline):
        """Lspline solver 
        
            for the 
            
            CProfile keel curve
            
         (intended for a transom sterned hull)
        with bulbous bow
        
        dev
        ----------
        
        elevation = self.MaxDraft
        k           = self.k
        nump        = self.nump
        xb          = 0.
        yb          = 0.
        xe          = self.LengthDWL
        ye          = 0.+self.MaxDraft*drop_aft
        #alphab     = 85.##35
        #alphae     = -85.##-35
        Cab_given   = 0.
        Cae_given   = 0.
        slope       = 'down'

        Dmax = self.MaxDraft
        area = self.CLarea

        ab = alphab
        ae = alphae
        """
        print 'OSV_bbow solver Oct 14, 2018'
        #
        #**********************************************************************
        # get weights
#        if Lspline is not None:
#            norms = package_norms(Lspline)
#            E1_0 = norms[0]
#            E2_0 = norms[1]
#            E3_0 = norms[2]
#            S_0  = norms[3]
#            curve = Lspline.curve
#        else:
        E1_0 = 1.
        E2_0 = 1.
        E3_0 = 1.
        S_0  = 1.
        # CPK
        #**********************************************************************
        # initialize
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve)
        #
        #**********************************************************************
        # area  CPKeelProfile stage 1
        FPD.add_AreaConstraint(kind='equality',
                               value = area)
        #
        #**********************************************************************
        # Xc CPKeelProfile
        if Xc is None:
            print 'CPK OSV bbow, equality LCG'
            FPD.add_XcConstraint(kind='LS', 
                                 value=LCG*.75,
                                 weight = 100.) #cLProfile flat to bulb...
            #FPD.add_XcConstraint(kind='max', 
            #                     value=self.LCG*.95)
            #FPD.add_XcConstraint(kind='min', 
            #                     value=self.LCG*0.8)
        else:
            print 'Centerplane setting Xc constraint to ',Xc
            FPD.add_XcConstraint(kind='LS', 
                                 value=Xc)
        #                                     value = 0.)
        #
        ##**************************************************************
        # CPkeelProfile  
        #print 'Adding x - fixity!'
        #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # TLM, straight flat stem
#            FPD.add_xFixity(index=1,
#                            value=xb)
#            FPD.add_xFixity(index=2,
#                            value=xb)
        #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # Dr. Birk, flat stem
#            FPD.add_AngleConstraint(
#                                    kind        = 'equality',
#                                    location    = 0.,
#                                    value       = 75.)
        FPD.add_CurvatureConstraint(
                                kind        = 'LS',
                                location    = 0.,
                                value       = 0.,
                                weight      = 100.)
        #--------------------------------------
        # +++
        #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        # TLM, straight flat stem
        #FPD.add_xVertexConstraint(kind='LS',
        #                         index = 3,
        #                         value = 0.)
        FPD.add_xFixity(index=3, #5,
                        value=0.)
        FPD.add_relative_xVertexConstraint(
                                    kind = 'equality',
                                    value = 0.,
                                    index = 3,
                                    index2 = 4,
                                    seperation=0.)
        #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        #--------------------------------------
        # ++
#            FPD.add_xVertexConstraint(kind='LS',
#                                     index = 8,
#                                     value = xe)
#            FPD.add_xVertexConstraint(kind='LS',
#                                     index = 9,
#                                     value = xe)
        #--------------------------------------
        ##**************************************************************
        # CPkeelProfile  
        #print 'Adding x - fixity!'
        #tlmOct 14 2018 tweak this:
        FPD.add_xFixity(index=5, #5,
                        value=FOCP0)
        FPD.add_xFixity(index=6, #6,
                        value=FOCP1)
        #
        #------------------------------- #tlmOct 14 2018 tweak this:
        #print 'Adding y - fixity!'
        #FPD.add_yFixity(index=3,
        #                value=Dmax)
        FPD.add_yFixity(index=4,
                        value=Dmax)
        FPD.add_yFixity(index=5,
                        value=Dmax)
        FPD.add_yFixity(index=6,
                        value=Dmax)
        FPD.add_yFixity(index=7,
                        value=Dmax)
        #-------------------------------
        FPD.add_yFixity(index=8,
                        value=Dmax)
        #FPD.add_yFixity(index=9,
        #                value=Dmax)
        #?-------------------------------?
        #            FPD.add_yVertexConstraint(kind='max',
        #                                     index = 8,
        #                                     value = Dmax)
        #            FPD.add_yVertexConstraint(kind='max',
        #                                     index = 9,
        #                                     value = Dmax)
        #?-------------------------------?
        #
        #
        #**********************************************************************
        # CPK fairness norms
        FPD.add_E1(kind='LS', weight = 1./E1_0)
        FPD.add_E2(kind='LS', weight = .5/E2_0)
        FPD.add_E3(kind='LS', weight = .01/E3_0)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.5/S_0)
        #
        #**********************************************************************
        #Setup and Solve:
        #
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        print 'Doing first stage CPkeel solve (include bbow weighting)'
        Lspline.optimize(stop=50)
        return Lspline
    
    
    
    
    
    #
    #**********************************************************************
    # CPK save init state
    inicurve = copy.deepcopy(curve)
    #func1_for_fix = fix1()
        
    #Lspline = OSV_1(curve)
    Lspline = OSV_bbow(curve)
    #
    #**********************************************************************
    #
    ckit = False
    if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
        Lspline = OSV_bbow2(Lspline)
        ckit = False
    else:    
        Lcopy = copy.deepcopy(Lspline)
        ckit = True
        
        
    curve = Lspline.curve
    return Lspline




def initial_curve(start=(0.,12.), end=(12.,0.), num=7, k=4,nump=50):
    vertices = linear_vertices(start, end, num)
    k=4
    nump_=nump
    curve = spline.Bspline(vertices, k, nump_)
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
        


if __name__ == "__main__":
    # TLM B-Spline Curve class
    import curve             as     spline
    from   initialValues     import InitializeControlPoints, InitializeControlVertices
    import copy
    from   ADILS             import IntervalLagrangeSpline, Lagrangian
    from   FormParameter     import FormParameterDict
    from   initialValues     import interval_bounds, lagrangian_bounds
    
    
    import utility_optimization as uopt
    package_norms = uopt.package_norms
    ##
    ##
    ##
    doall = False
    ##
    ##*****************************************
    ## Define Stations
    
    
    cp = 0.5
    
    ##
    ##*****************************************
    ##
    
    Acp = 1544.25474515
    
    Ccp = 0.87039360383
    
    LCG = 57.2087722367
    
    lwl = 117.625323999
    
    draft = 15.0902796256
    MaxDraft = 15.087279625615295
    
    # SD.hull.hyperparameter_input.aft_drop
    drop_aft = .25
    
    CLarea = 1544.2517451500446 #SD.hull.CLarea
    
    FOCP0 = 52.10086641415609# SD.hull.stations.FOCP[0]
    
    FOCP1 = 65.52145758493774 # SD.hull.stations.FOCP[1]
    ##
    ##*****************************************
    ##
    
    if doall:
       Lspline = compute_cLProfile()
       Lspline.curve.plotcurve_detailed()
    
    
    
    
    
    ##
    ##*****************************************
    ##
    
    # SD.design_space.rgp.get_name_and_value('Acp')[1][0].getpoint(.5)
    Acp = 1562.208645198735
    
    # SD.design_space.rgp.get_name_and_value('Ccp')[1][0].getpoint(.5)
    Ccp = 0.8750076650671084
    
    
    # SD.design_space.rgp.get_name_and_value('LCG')[1][0].getpoint(.5)
    LCG = 60.497434002285
    
    
    # SD.design_space.rgp.get_name_and_value('lwl')[1][0].getpoint(.5)
    lwl = 118.21419378914734
    
    
    # SD.design_space.rgp.get_name_and_value('draft')[1][0].getpoint(.5)
    draft = 15.10280297957117
    
    #SD.hull.MaxDraft
    MaxDraft = 15.10280297957117
    
    # SD.hull.hyperparameter_input.aft_drop
    drop_aft = .25
    
    CLarea = 1562.208645198735   #SD.hull.CLarea
    
    FOCP0 = 54.428736216476295    # SD.hull.stations.FOCP[0]
    
    FOCP1 = 63.78545757267105    # SD.hull.stations.FOCP[1]
    ##
    ##*****************************************
    ##
    
    if doall:
        Lspline = compute_cLProfile()
        Lspline.curve.plotcurve_detailed()
    
    
    
    
    
    
    ##
    ##*****************************************
    ##
    
    """
    print 'Acp =', SD.design_space.rgp.get_name_and_value('Acp')[1][0].getpoint(.5)
    
    print 'Ccp = ', SD.design_space.rgp.get_name_and_value('Ccp')[1][0].getpoint(.5)
    
    print 'LCG = ', SD.design_space.rgp.get_name_and_value('LCG')[1][0].getpoint(.5)    
    
    print 'lwl = ', SD.design_space.rgp.get_name_and_value('lwl')[1][0].getpoint(.5)
    
    print 'draft = ', SD.design_space.rgp.get_name_and_value('draft')[1][0].getpoint(.5)
    
    try:
        print SD.hull
        print 'MaxDraft = ', SD.hull.MaxDraft
        print 'drop_aft = ', SD.hull.hyperparameter_input.aft_drop
        print 'CLarea = ', SD.hull.CLarea
        print 'FOCP0 = ', SD.hull.stations.FOCP[0]
        print 'FOCP1 = ', SD.hull.stations.FOCP[1]
        
        
    except:
        print 'MaxDraft = ', SD.design_space.rgp.get_name_and_value('draft')[1][0].getpoint(.5)
    

        print 'drop_aft = ', .25 #SD.hull.hyperparameter_input.aft_drop
    
        print 'CLarea = ', SD.design_space.rgp.get_name_and_value('Acp')[1][0].getpoint(.5)
    
    
        flat_bottom = SD.design_space.rgp.get_name_and_value('lfcp')[1][0].getpoint(.5)
        args = [flat_bottom,cp]
        len_flat, cp = args
        lwl =  121.288129841
        pflat  = len_flat/lwl
        a = (cp-.5*pflat)
        b = (cp+.5*pflat)
        FOCP =  [a*lwl,b*lwl,a,b]
        print 'FOCP0 = ',FOCP[0]
        print 'FOCP1 = ',FOCP[1]
    #"""
    
    ##
    ##*****************************************
    ##
    Acp = 1549.23787542
    Ccp =  0.840787447695
    LCG =  58.7799605381
    lwl =  121.288129841
    draft =  15.1919529408
    MaxDraft =  15.1919529408
    drop_aft =  0.25
    CLarea =  1549.23787542
    FOCP0 =  53.5087539896
    FOCP1 =  67.7793758514
    ##
    ##*****************************************
    ##
    
    
    Lspline = compute_cLProfile()
    Lspline.curve.plotcurve_detailed()
    
    
    
    
    
    ##
    ##*****************************************
    ##
    Acp = 1657.14867213
    Ccp =  0.868338674979
    LCG =  60.7096661116
    lwl =  120.039549083
    draft =  15.8982012853
    MaxDraft =  15.8982012853
    drop_aft =  0.25
    CLarea =  1657.14867213
    FOCP0 =  56.6168164234
    FOCP1 =  64.6713134176
    ##
    ##*****************************************
    ##
    Lspline = compute_cLProfile()
    Lspline.curve.plotcurve_detailed()