#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 16:42:38 2018

@author: luke

Break out the bare hull curve solvers into separate classes
to keep valuble data on hand locally
-hullrules thin ship design parameters
-heuristic methods of computing each curve (can become quite extensive
there is a need to reformulate to something amenable to hyperparameter search)

Here's an idea of possible relevance to both hyperparameter search
and to the interdependence of constraints within one curve, between curves, etc.
-so of relevance to the optimization of geometry itself, and
-also of relevance to the picking of which constraint methods to use to specify 
a particular geometry 


-interdependence of constraints can be seen geoemtrically, which led to the 
mk rules pre-geometry filtering of the design (configuration) space
-but it is classically measured in econometrics and the correlations
of various parameters with each other and with the error (same as our error)
-Can we use 2 stage least squares and IV techniques to uncorrelate out
conflicting constraints?

TODO:
    constraint choice is an interdependent multi-arm bandit
    problem.  
    -Each constraint is one 'slot machine'.
    -Payoff statistics depend on what other constraints
    have been set, assuming we say a particular choice oblem is 
    always with it's particular curve.
    
    
    

Idea:  given a sub problem in hull design,
    i.e., in this instance, the development of
    the optimum cLProfile Keel curve
    which is actually the leading transverse cufve in the 
    bare hull form,
    provide a set of alternatives which are always the same.
    
    -use RL to determine how to pick from the alternatives
    given the initial curve and the desired Form Parameter 
    Constraint values.
    (alternatives = alternative 'solvers', 
    which are function continuations containing sets of procedures
    to be used in reaching the desired curve, awaiting inicurve to start)
    -(later one can think of tweaking inside the solvers with more fine grained optimization
    but not this go-round)
    -drawbacks:
        Well, RL looks ahead to examine the worth of moves...
        -We can just look ahead at all the moves and pick the best.
        -It would be nice to learn better and better predictors for
        which curve constraint combos will be best
        -but this is more about efficiency than getting a feasible answer.
    -but if the parameter space is infinite, or just large, in any way,
    lookahead is only an approximate search!

"""
import numpy as np


class SAC(object):
    def __init__(self,k,nump,LengthDWL):
        self.k = k
        self.nump = nump
        self.LengthDWL = LengthDWL
    
    def modify_area(self, new_area):
        def setter(new_area, FPD):
            return FPD
        return setter
    def compute_SAC(self, LocMaxSection=.5,
                    Xc=None, flat=True):
        """
            TODO : use bow_fairness_section properties
            to fix the bow fairness curve feasability
        """
        print '\nSAC'
        k=self.k
        nump=self.nump
        xb = 0.
        yb = 0.
        xe = self.LengthDWL
        ye = 0.
        alphab = 15.#-5.
        alphae = -70.#-15.
        #Cab_given = 0.
        #Cae_given = 0.
        #slope = 'down'

        #ae = alphae
        #ab = alphab
        curve = initial_curve((xb,yb),
                              (xe,ye),
                              num=12, k=k,nump=nump)
        v = copy.deepcopy(curve.vertices)
        v[4:8,1] = self.maxSecArea
        v[3,1] = self.maxSecArea*2./3.
        v[8,1] = self.maxSecArea*2./3.

        curve.vertices = v #updates the curve automagically
        
        inicurve = copy.copy(curve)
        
        
        def stage1(curve):
            #
            #******************************************************************
            # Begin Stage 1
            interval_data, small = interval_bounds(curve)
            
            FPD = FormParameterDict(curve)
            #
            #******************************************************************
            #  Integral Parameters:
            FPD.add_AreaConstraint(kind='equality',
                                   value = self.BareHullVolume)
            FPD.add_XcConstraint(kind='equality', 
                                 value=self.LCG)
            #
            #******************************************************************
            #  Flat Portion
            if flat:
                FPD.add_yPointConstraint(kind='equality',
                                     location = 0.3333333333,
                                     value = self.maxSecArea)#*2./3.)
    
                FPD.add_yPointConstraint(kind='equality',
                                     location = 0.6666666666,
                                     value = self.maxSecArea)#*2./3.)
            else:
                print 'SAC stage 1 not enforcing flat section'
            #
            #******************************************************************
            # starting, ending angles
            FPD.add_AngleConstraint(kind='LS',
                                    location = 0.,
                                    value = alphab)
            FPD.add_AngleConstraint(kind='LS',
                                    location = 1.,
                                    value = alphae)
            #
            #******************************************************************
            # Max S.A.
            FPD.add_yPointConstraint(kind='equality',
                                     location = LocMaxSection,
                                     value = self.maxSecArea)
            #
            #******************************************************************
            # Fairness
            FPD.add_E1(kind='LS', weight = .5)
            FPD.add_E2(kind='LS', weight = .5)
            FPD.add_E3(kind='LS', weight = .5) 
            FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            #
            #******************************************************************
            #
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            #
            #******************************************************************
            #OPTIMIZE - 1st pass
            #
            print 'Doing first stage SAC solve'
            Lspline.optimize()
            #
            #******************************************************************
            # Done
            self.Lsplines.SAC=Lspline
            self.SAC = Lspline.curve
            #
            return Lspline,FPD
        
        
        def stage2_elegant(curve,Lspline=None,
                           fi = .3333333333, fe = .6666666666):
            """the original way to get a flat SAC curve
            
                usage of norms:
                    given some Lspline, then:
                    Lspline.curve.compute_arclength()
                    E1_0 = copy.deepcopy(Lspline.curve.E1)
                    E2_0 = copy.deepcopy(Lspline.curve.E2)
                    E3_0 = copy.deepcopy(Lspline.curve.E3)
                    S_0 = copy.deepcopy(Lspline.curve.AL)
                    
                    norms = [E1_0,E2_0,E3_0,S_0]
            """
            #
            #******************************************************************
            # Get scalars
            if Lspline is not None:
                norms = self.package_norms(Lspline)
            #
            #******************************************************************
            # Begin standard procedure
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve)
            #
            #******************************************************************
            #  Integral Parameters:
            FPD.add_AreaConstraint(kind='equality',
                                   value = self.BareHullVolume)
            FPD.add_XcConstraint(kind='equality', 
                                 value=self.LCG)
            #
            #******************************************************************
            #  Flat Portion
            #fi = .3333333333
            #fe = .6666666666
            #
            #******************************************************************
            # flat x
            FPD.add_yPointConstraint(kind='equality',
                                 location = fi,
                                 value = self.maxSecArea)

            FPD.add_yPointConstraint(kind='equality',
                                 location = fe,
                                 value = self.maxSecArea)
            #
            #******************************************************************
            # Flat derivatives
            FPD.add_AngleConstraint(kind='equality',
                                    location = fi,
                                    value = 0.)
            FPD.add_AngleConstraint(kind='equality',
                                    location = fe,
                                    value = 0.)
            
            FPD.add_CurvatureConstraint(kind='equality',
                                    location = fi,
                                    value = 0.)
            FPD.add_CurvatureConstraint(kind='equality',
                                    location = fe,
                                    value = 0.)
            #
            #******************************************************************
            # flat x
            s1 = self.stations.FOSAC[0]
            s2 = self.stations.FOSAC[1]
            FPD.add_xPointConstraint(kind='equality',
                                 location = fi,
                                 value = s1)
            FPD.add_xPointConstraint(kind='equality',
                                 location = fe,
                                 value = s2)
            #
            #******************************************************************
            # Fairness
            if norms is None:
                FPD.add_E1(kind='LS', weight = .5)
                FPD.add_E2(kind='LS', weight = .5)
                FPD.add_E3(kind='LS', weight = .5) 
                FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            else:
                E1_0    = norms[0]
                E2_0    = norms[1]
                E3_0    = norms[2]
                S_0     = norms[3]
                FPD.add_E1(kind='LS', weight = .5/E1_0.value)
                FPD.add_E2(kind='LS', weight = .5/E2_0.value)
                FPD.add_E3(kind='LS', weight = .5/E3_0.value)
                FPD.add_ArcLengthApprox(kind='LS', weight = .5/S_0)
            #
            #******************************************************************
            #
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            #
            #******************************************************************
            #OPTIMIZE - elegant pass
            #
            print 'Trying for elegant SAC solve'
            Lspline.optimize()
            #
            #******************************************************************
            # Done
            self.Lsplines.SACe=Lspline
            self.SACe = Lspline.curve
            #
            return Lspline,FPD
        
        def stage2(Lspline,FPD):
            #
            #******************************************************************
            # NOW Entering Solver Stage 2
            Lspline.curve.pts_M_pts() 
            Lspline.curve.compute_arclength()
            E1_0 = copy.copy(Lspline.curve.E1)
            E2_0 = copy.copy(Lspline.curve.E2)
            E3_0 = copy.copy(Lspline.curve.E3)
            S_0 = copy.copy(Lspline.curve.AL)
            
            curve = Lspline.curve
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve)
            #
            #******************************************************************
            # Integral Parameters
            FPD.add_AreaConstraint(kind='equality',
                                   value = self.BareHullVolume)
            FPD.add_XcConstraint(kind='equality', value=self.LCG)
            #
            #******************************************************************
            # do not use a fine grained curve to drive a coarse one!
            #if self.fairbow:
            #
            #if self.fairnesscurves_drive_hull_curves:
            #    start,end,ao,bo = self.stations.Bow_Fairing
            #    AreaBowFairingCurve = self.Abfc#/2.
            #    where = ao
            #    FPD.add_yPointConstraint(kind='equality',
            #                             location = where,
            #                             value = AreaBowFairingCurve)
            #
            #******************************************************************
            # FLAT Portion
            #
            #FPD.add_AngleConstraint(kind='LS',
            #                        location = 0.,
            #                        value = alphab)
            #FPD.add_AngleConstraint(kind='LS',
            #                        location = 1.,
            #                        value = alphae)
            """
            IDEA:
            instead of points and angle contraint, use...
            
            -three vertices collinear
            -x vertex location same as .3333 .666
            constraints below
            
            This may be the only way to approach
            the performance of the 3 split SAC components
            
            Note that this may be even more parsimoneous than
            that was- it allows the flat spot to 
            be spec'd in 3 vertices when curve order, k, is 4
            """
            #FPD.add_yPointConstraint(kind='LS',
            #                         location = 0.3333333333,
            #                         value = self.maxSecArea)
            #FPD.add_yPointConstraint(kind='LS',
            #                         location = LocMaxSection,
            #                         value = self.maxSecArea)
            #FPD.add_yPointConstraint(kind='LS',
            #                         location = 0.6666666666,
            #                         value = self.maxSecArea)
            
            #s1 = self.stations.FOSAC[0]
            #s2 = self.stations.FOSAC[1]
            #hdx = s2-s1
            #sepvec = np.linspace(0,hdx,3)#,4)#
#            FPD.add_xPointConstraint(kind='equality',
#                                     location = 0.3333333333,
#                                     value = s1)
#            FPD.add_xPointConstraint(kind='equality',
#                                     location = 0.6666666666,
#                                     value = s2)
            
            #
            #******************************************************************
            # X FLAT SAC (4,5,6 + keep 3,7 out)
            FPD.add_xVertexConstraint(kind='max',
                                     index = 3,
                                     value = self.stations.FOSAC[0])
            FPD.add_xVertexConstraint(kind='equality',
                                     index = 4,
                                     value = self.stations.FOSAC[0])
            FPD.add_xVertexConstraint(kind='min',
                                     index = 5,
                                     value = self.stations.FOSAC[0])
            FPD.add_xVertexConstraint(kind='max',
                                     index = 5,
                                     value = self.stations.FOSAC[1])
            FPD.add_xVertexConstraint(kind='equality',
                                     index = 6,
                                     value = self.stations.FOSAC[1])
            FPD.add_xVertexConstraint(kind='min',
                                     index = 7,
                                     value = self.stations.FOSAC[1])
            #
            #******************************************************************
            # Y FLAT SAC (3,4,5,6,7)
            FPD.add_yVertexConstraint(kind = 'equality', 
                                      value=self.maxSecArea, 
                                      index =3)
            FPD.add_yVertexConstraint(kind = 'equality', 
                                      value=self.maxSecArea, 
                                      index =4)
            FPD.add_yVertexConstraint(kind = 'equality', 
                                      value=self.maxSecArea, 
                                      index =5)
            FPD.add_yVertexConstraint(kind = 'equality', 
                                      value=self.maxSecArea, 
                                      index =6)
            FPD.add_yVertexConstraint(kind = 'equality', 
                                      value=self.maxSecArea, 
                                      index =7)
            #
            #******************************************************************
            # Let's Straighten up the Sides of the SACurve - fit within the 'box'
            #
            # X limits:
            #FPD.add_xVertexConstraint(kind='min',
            #                         index = 1,
            #                         value = 0.)
            #
            #******************************************************************
            # not used:
            #            FPD.add_relative_xVertexConstraint(kind = 'min', 
            #                                               location=None, 
            #                                               value = 0., 
            #                                               weight = 1.0, 
            #                                               index = 1, 
            #                                               index2 = 2, 
            #                                               seperation = 0. )
            #            FPD.add_relative_xVertexConstraint(kind = 'min', 
            #                                               location=None, 
            #                                               value = 0., 
            #                                               weight = 1.0, 
            #                                               index = 2, 
            #                                               index2 = 3, 
            #                                               seperation = 0. )
            #
            #******************************************************************
            # MAKE X go ''out''
            #            FPD.add_xVertexConstraint(kind='LS',
            #                                     index = 1,
            #                                     value = 0.,
            #                                     weight=.5/E1_0.value)
            #            FPD.add_xVertexConstraint(kind='LS',
            #                                     index = 2,
            #                                     value = 0.,
            #                                     weight=0.5/E1_0.value)
            #            FPD.add_xVertexConstraint(kind='LS',
            #                                     index = 3,
            #                                     value = 0.,
            #                                     weight=0.5/E1_0.value)
            #            
            #            FPD.add_xVertexConstraint(kind='LS',
            #                                     index = 7,
            #                                     value = self.LengthDWL,
            #                                     weight=0.5/E1_0.value)
            #            FPD.add_xVertexConstraint(kind='LS',
            #                                     index = 8,
            #                                     value = self.LengthDWL,
            #                                     weight=0.5/E1_0.value)
            #            FPD.add_xVertexConstraint(kind='LS',
            #                                     index = 9,
            #                                     value = self.LengthDWL,
            #                                     weight=0.5/E1_0.value)
            #            FPD.add_xVertexConstraint(kind='LS',
            #                                     index = 10,
            #                                     value = self.LengthDWL,
            #                                     weight=.5/E1_0.value)
            #            FPD.add_xVertexConstraint(kind='max',
            #                                     index = 10,
            #                                     value = self.LengthDWL)
            #
            #******************************************************************
            # Finally make y 'want' to go up, instead of pushing x to balloon out.
            # -still getting a weird hook in SAC
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 2,
            #                                     value = self.maxSecArea,
            #                                     weight=1./E1_0.value)
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 3,
            #                                     value = self.maxSecArea,
            #                                     weight=1./E1_0.value)
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 7,
            #                                     value = self.maxSecArea,
            #                                     weight=1./E1_0.value)
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 8,
            #                                     value = self.maxSecArea,
            #                                     weight=1./E1_0.value)
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 9,
            #                                     value = self.maxSecArea,
            #                                     weight=.1/E1_0.value)
            #            FPD.add_yVertexConstraint(kind='max',
            #                                     index = 9,
            #                                     value = self.maxSecArea,
            #                                     weight=1./E1_0.value)
            #            # also, make sure y stays positive:
            #            FPD.add_yVertexConstraint(kind='min',
            #                                     index = 9,
            #                                     value = 0.,
            #                                     weight=1./E1_0.value)
            #            FPD.add_yVertexConstraint(kind='min',
            #                                     index = 10,
            #                                     value = 0.,
            #                                     weight=1./E1_0.value)
            #
            #******************************************************************
            # Fairness
            FPD.add_E1(kind='LS', weight = 1.5/E1_0.value)
            FPD.add_E2(kind='LS', weight = .5/E2_0.value)
            FPD.add_E3(kind='LS', weight = .5/E3_0.value)
            FPD.add_ArcLengthApprox(kind='LS', weight = 10./S_0)
            #
            #******************************************************************
            # 
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            #
            #******************************************************************
            #OPTIMIZE - 2nd pass
            print 'Doing second stage SAC solve'
            Lspline.optimize(stop = 100)
            #
            #******************************************************************
            # Done
            return Lspline, FPD
        
        def stagefix1(curve,Lspline):
            #
            #******************************************************************
            # Begin Stage 1
            interval_data, small = interval_bounds(curve)
            
            FPD = FormParameterDict(curve)
            #
            #******************************************************************
            #
            #**********************************************************************
            # NOW Entering Solver fixup stages - because hyperparameter 
            # auto-tuning is absent
            # lacking to date.  ..Future work: consider Domn.Speci.Languge. expresivity  
            # for hyperparameters 
            Lspline.curve.pts_M_pts() 
            Lspline.curve.compute_arclength()
            #
            E1_0 = copy.copy(Lspline.curve.E1)
            E2_0 = copy.copy(Lspline.curve.E2)
            E3_0 = copy.copy(Lspline.curve.E3)
            S_0 = copy.copy(Lspline.curve.AL)
            #
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve)
            #
            #******************************************************************
            # Integral Parameters
            FPD.add_AreaConstraint(kind='equality',
                                   value = self.BareHullVolume)
            FPD.add_XcConstraint(kind='equality', value=self.LCG)
            #
            #******************************************************************
            # Fairness
            FPD.add_E1(kind='LS', weight = 1.5/E1_0.value)
            FPD.add_E2(kind='LS', weight = .5/E2_0.value)
            FPD.add_E3(kind='LS', weight = .5/E3_0.value)
            FPD.add_ArcLengthApprox(kind='LS', weight = 10./S_0)
            #
            #******************************************************************
            # 
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            #
            #******************************************************************
            #OPTIMIZE - 1st stage fix
            print 'SAC fixit stage1 solve:'
            Lspline.optimize(stop = 50)
            return Lspline
        
        
        def stagefix2(Lspline):
            #
            #**********************************************************************
            # Fixit stage 2: flat
            Lspline.curve.pts_M_pts() 
            Lspline.curve.compute_arclength()
            #
            E1_0 = copy.copy(Lspline.curve.E1)
            E2_0 = copy.copy(Lspline.curve.E2)
            E3_0 = copy.copy(Lspline.curve.E3)
            S_0 = copy.copy(Lspline.curve.AL)
            #
            curve = Lspline.curve
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve)
            #
            #******************************************************************
            # Integral Parameters
            FPD.add_AreaConstraint(kind='equality',
                                   value = self.BareHullVolume)
            FPD.add_XcConstraint(kind='equality', value=self.LCG)
            #
            #******************************************************************
            # y flat
            return
    
class DWL(object):
    def __init__(self):
        pass
    def compute_DWL(self, LocMaxSection=.5, Xc=None,
                    flat=True, height=None, alphab=0.,alphae=0.):
        print '\nDWL'
        if height is None:
            elevation = self.MaxDraft
        else:
            elevation = height
        k       = self.k
        nump    = self.nump
        xb      = 0.
        yb      = 0.
        xe      = self.LengthDWL
        ye      = 0.
        #alphab = 0.#-5.
        #alphae = 0.#-15.
        Cab_given = 0.
        Cae_given = 0.
        slope = 'down'

        Bmax = self.MaxBreadth/2.
        area = self.WLarea/2.

        #ab = alphab
        #ae = alphae

        curve = initial_curve((xb,yb),
                              (xe,ye),
                              num=11, k=k,nump=nump)
        v = copy.deepcopy(curve.vertices)
        v[4:8,1] = Bmax
        v[3,1] = Bmax*2./3.
        v[8,1] = Bmax*2./3.

        curve.vertices = v #updates the curve automagically
        def stage1(curve):
            #
            #******************************************************************
            # Begin Stage 1
            #
            #**********************************************************************
            # initialize
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve)
            #
            #**********************************************************************
            # area
            FPD.add_AreaConstraint(kind='equality',
                                   value = area)
            #
            #**********************************************************************
            #   DWL/2 might be a good starting parameter for a work deck (not DWL)
            if Xc is None:
                FPD.add_XcConstraint(kind='LS', 
                                     value=self.LCG)
            else:
                FPD.add_XcConstraint(kind='equality', value=Xc)
            #
            #******************************************************************
            # nose intercept should be purely transverse 
            #- no longitudinal slope
            
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 1,
            #                                     value = Bmax/4.)
            #
            #******************************************************************
            # x-flat DWL  (4,5,6 + optional:  keep 3,7 out)
            FPD.add_xVertexConstraint(kind='max',
                                     index = 3,
                                     value = self.stations.FOWL[0])
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 4,
                                     value = self.stations.FOWL[0])
            FPD.add_xVertexConstraint(kind='min',
                                     index = 5,
                                     value = self.stations.FOWL[0])
            FPD.add_xVertexConstraint(kind='max',
                                     index = 5,
                                     value = self.stations.FOWL[1])
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 6,
                                     value = self.stations.FOWL[1])
            FPD.add_xVertexConstraint(kind='min',
                                     index = 7,
                                     value = self.stations.FOWL[1])
            #
            #******************************************************************
            # y-flat DWL (3,4,5,6,7)
            # ensure curvature flatness at vertex 4 via vertex 3 (and 5,6,7)
            FPD.add_yVertexConstraint(kind='LS',
                                     index = 3,
                                     value = Bmax)
            # Now the mid section 4,5,6
            FPD.add_yVertexConstraint(kind='LS',
                                     index = 4,
                                     value = Bmax)
            #ensures the solver doesn't put 
            # a bulge right in the middle (to meet the area constraint in backhanded fashion)
            FPD.add_yVertexConstraint(kind='LS',
                                     index = 5,
                                     value = Bmax)
            #
            FPD.add_yVertexConstraint(kind='LS',
                                     index = 6,
                                     value = Bmax)    
            #curvature continuity
            FPD.add_yVertexConstraint(kind='LS',
                                     index = 7,
                                     value = Bmax)  
            #
            #******************************************************************
            #    
            #FPD.add_E1(kind='LS', weight = 1.)
            #FPD.add_E2(kind='LS', weight = .5)
            #FPD.add_E3(kind='LS', weight = .5) #always commented out TLM OCT 27 2017
            #FPD.add_ArcLengthApprox(kind='LS', weight = 1.5)
            
            FPD.add_E1(kind='LS', weight = .5)
            FPD.add_E2(kind='LS', weight = .5)
            FPD.add_E3(kind='LS', weight = .001)
            FPD.add_ArcLengthApprox(kind='LS', weight = 1.5)
            
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            # make curve
            print 'Doing first stage DWL solve'
            Lspline.optimize(stop=50)
            #Lspline.curve.knot_insertion(.5) #get to 11
            self.Lsplines.DWL=Lspline
            self.DWL2D = Lspline.curve
            
            return Lspline
        
        def stage2(Lspline):
            """DWL stage 2
            """
            norms = self.package_norms(Lspline)
            E1_0 = norms[0]
            E2_0 = norms[1]
            E3_0 = norms[2]
            S_0  = norms[3]
            curve = Lspline.curve
            #
            #**********************************************************************
            # initialize
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve)
            #
            #*********************************************************************
            # area
            FPD.add_AreaConstraint(kind='equality',
                                   value = area)
            #
            #*********************************************************************
            # X-centroid
            #   DWL/2 might be a good parameter for a work deck (not DWL)
            if Xc is None:
                FPD.add_XcConstraint(kind='LS', 
                                     value=self.LCG)
            else:
                print 'DWL setting Xc constraint to ',Xc
                FPD.add_XcConstraint(kind='equality', value=Xc)
            #
            #******************************************************************
            # NOSE CONSTRAINTS
            # nose intercept should be purely transverse 
            #- no longitudinal slope
            FPD.add_relative_xVertexConstraint(
                                     kind='equality',
                                     location=None,
                                     index = 0,
                                     index2 = 1,
                                     value = 0.,
                                     weight=1.)
            
            FPD.add_yVertexConstraint(kind='min',
                                     index = 1,
                                     value = 0.)
            # end of nose tangent constraint 
            #******************************************************************
            # pull DWL vertices outwards:
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 2,
            #                                     value = Bmax)#/2.)   
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 3,
            #                                     value = Bmax)#/2.) 
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 7,
            #                                     value = Bmax)
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 8,
            #                                     value = Bmax)
            #
            #**********************************************************************
            # DWL stay within parameters
            FPD.add_xVertexConstraint(kind='min',
                                     index = 2,
                                     value = 0.)
            FPD.add_xVertexConstraint(kind='max',
                                     index = 8,
                                     value = 1.05*self.LengthDWL)
            FPD.add_xVertexConstraint(kind='max',
                                     index = 9,
                                     value = 1.05*self.LengthDWL)
            #
            #******************************************************************
            # DWL ensure no 0 crossing
            FPD.add_yVertexConstraint(kind='min',
                                     index = 2,
                                     value = 0.)
            FPD.add_yVertexConstraint(kind='min',
                                     index = 8,
                                     value = 0.)
            FPD.add_yVertexConstraint(kind='min',
                                     index = 9,
                                     value = 0.)
            #
            #**********************************************************************
            # x-flat DWL  (4,5,6 + optional:  keep 3,7 out)
            """based on n=10 for this curve
            4 free in front
            3 free in back
            """
            FPD.add_xVertexConstraint(kind='max',
                                     index = 3,
                                     value = self.stations.FOWL[0])
            FPD.add_xVertexConstraint(kind='equality',
                                     index = 4,
                                     value = self.stations.FOWL[0])
            FPD.add_xVertexConstraint(kind='min',
                                     index = 5,
                                     value = self.stations.FOWL[0])
            FPD.add_xVertexConstraint(kind='max',
                                     index = 5,
                                     value = self.stations.FOWL[1])
            FPD.add_xVertexConstraint(kind='equality',
                                     index = 6,
                                     value = self.stations.FOWL[1])
            FPD.add_xVertexConstraint(kind='min',
                                     index = 7,
                                     value = self.stations.FOWL[1])
            #
            #**********************************************************************
            # y-flat DWL (3,4,5,6,7)
            FPD.add_yVertexConstraint(kind='equality',
                                     index = 3,
                                     value = Bmax)  
            # Now the mid section 4,5,6
            FPD.add_yVertexConstraint(kind='equality',
                                     index = 4,
                                     value = Bmax)
            FPD.add_yVertexConstraint(kind='equality',
                                     index = 5,
                                     value = Bmax)
            FPD.add_yVertexConstraint(kind='equality',
                                     index = 6,
                                     value = Bmax)  
            FPD.add_yVertexConstraint(kind='equality',
                                     index = 7,
                                     value = Bmax)    
            #
            #**********************************************************************
            # DWL fairness norms
            FPD.add_E1(kind='LS', weight = .5/E1_0.value)
            FPD.add_E2(kind='LS', weight = .5/E2_0.value)
            FPD.add_E3(kind='LS', weight = .5/E3_0.value)
            FPD.add_ArcLengthApprox(kind='LS', weight = .5/S_0)
            #
            #**********************************************************************
            # setup
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            #
            #**********************************************************************
            # solve
            print 'Doing second stage DWL solve'
            Lspline.optimize(stop = 100)
            return Lspline
        
        #
        #**********************************************************************
        # save init state
        inicurve = copy.copy(curve)
        #
        #**********************************************************************
        # First solve:
        Lspline = stage1(curve)
        #
        #**********************************************************************
        # Second solve:
        Lspline = stage2(Lspline)
        #
        #**********************************************************************
        #check solution:
        #        if Lspline.error_code is 'max_iter':
        #            pass 
        #            #Lspline = stage_fix(Lspline)
        #        elif Lspline.error_code is 'NAN':
        #            Lspline = stage1(curve)
        #            Lspline = stage_fix(Lspline)=
        #
        #**********************************************************************
        # save
        self.Lsplines.DWL=Lspline
        self.DWL2D = Lspline.curve
        #
        #**********************************************************************
        # rotate to y-z plane and move x=DWL
        verts = Lspline.curve.vertices
        #
        #**********************************************************************
        #  lift to 3D
        vnew = frame.xy_to_zxy_elevate( verts,
                                            elevation=elevation)
        #
        #**********************************************************************
        # SAVE
        self.DWL = spline.Bspline(vnew, k, nump)
        
        return
    
class CPKeelProfile(object):
    def __init__(self):
        pass
    def compute_cLProfile(self, LocMaxSection=None, Xc=None,
                          flat=False, height=None, alphab=0.,alphae=0.):
        print '\n hull_from_simple_designspace: compute_cLProfile (Keel Profile)'
        if height is None:
            elevation = self.MaxDraft
        else:
            elevation = height
        k       = self.k
        nump    = self.nump
        xb      = 0.
        yb      = 0.
        xe      = self.LengthDWL
        ye      = 0.
        #alphab = 85.##35
        #alphae = -85.##-35
        Cab_given = 0.
        Cae_given = 0.
        slope = 'down'

        Dmax = self.MaxDraft
        area = self.CLarea

        ab = alphab
        ae = alphae

        curve = initial_curve((xb,yb),
                              (xe,ye),
                              num=11, k=k, nump=nump)
        v = copy.deepcopy(curve.vertices)
        v[4:8,1] = Dmax
        v[3,1] = Dmax*2./3.
        v[8,1] = Dmax*2./3.

        curve.vertices = v #updates the curve automagically
        def stage1(curve):
            #
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
                FPD.add_XcConstraint(kind='LS', 
                                     value=self.LCG)
            else:
                print 'Centerplane setting Xc constraint to ',Xc
                FPD.add_XcConstraint(kind='equality', value=Xc)
            #
            #******************************************************************
            #  CPKeelProfile stage 1
            # end condition stays nice:
            #            FPD.add_yVertexConstraint(kind='max',
            #                                     index = 1,
            #                                     value = Dmax*1.3)
            #            FPD.add_yVertexConstraint(kind='min',
            #                                     index = 1,
            #                                     value = 0.)
            #            FPD.add_yVertexConstraint(kind='max',
            #                                     index = 2,
            #                                     value = Dmax)
            #            FPD.add_yVertexConstraint(kind='max',
            #                                     index = 8,
            #                                     value = Dmax)
            #
            #
            #            FPD.add_yVertexConstraint(kind='min',
            #                                     index = 2,
            #                                     value = 0.)
            #            FPD.add_yVertexConstraint(kind='min',
            #                                     index = 8,
            #                                     value = 0.)
            #
            #******************************************************************
            # x-flat CPkeelProfile  (4,5,6 + optional:  keep 3,7 out)
            # stage 1 solve
            #        FPD.add_xVertexConstraint(kind='equality',
            #                                 index = 3,
            #                                 value = self.stations.FOCP[0])
            #        FPD.add_xVertexConstraint(kind='equality',
            #                                 index = 7,
            #                                 value = self.stations.FOCP[1])
            FPD.add_xVertexConstraint(kind='max',
                                     index = 3,
                                     value = self.stations.FOCP[0])
            FPD.add_xVertexConstraint(kind='equality',
                                     index = 4,
                                     value = self.stations.FOCP[0])
            FPD.add_xVertexConstraint(kind='min',
                                     index = 5,
                                     value = self.stations.FOCP[0])
            FPD.add_xVertexConstraint(kind='max',
                                     index = 5,
                                     value = self.stations.FOCP[1])
            FPD.add_xVertexConstraint(kind='equality',
                                     index = 6,
                                     value = self.stations.FOCP[1])
            FPD.add_xVertexConstraint(kind='min',
                                     index = 7,
                                     value = self.stations.FOCP[1])
            #
            #**********************************************************************
            # y-flat CPkeelProfile (3,4,5,6,7) stage 1 solve
            FPD.add_yVertexConstraint(kind='equality',
                                     index = 3,
                                     value = Dmax)
            FPD.add_yVertexConstraint(kind='equality',
                                     index = 4,
                                     value = Dmax)
            FPD.add_yVertexConstraint(kind='equality',
                                     index = 6,
                                     value = Dmax)
            FPD.add_yVertexConstraint(kind='equality',
                                     index = 7,
                                     value = Dmax)
            FPD.add_yVertexConstraint(kind='equality',
                                     index = 5,
                                     value = Dmax)
            #
            #**********************************************************************
            # curvd cL_Keelprofile stern fairness curve attachment:
            #        if self.fairnesscurves_drive_hull_curves:
            #            start,end,ao,bo = self.stations.Stern_Fairing
            #            where = ao
            #            FPD.add_yPointConstraint(kind='equality',
            #                                     location = where,
            #                                     value = self.dsfc)
            #            FPD.add_xPointConstraint(kind='equality',
            #                                     location = where,
            #                                     value = start)
            #
            #----------------------- Fairness Functional:
            #
            #FPD.add_E1(kind='LS', weight = 1.) 
            #FPD.add_E2(kind='LS', weight = .5)
            #FPD.add_E3(kind='LS', weight = .5) #always commented out TLM Oct 28 2017
            #FPD.add_ArcLengthApprox(kind='LS', weight = 1.5)
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
            print 'Doing first stage CPkeel solve'
            Lspline.optimize()
            return Lspline
        
        
        
        def stage2(Lspline):
            """CLKeelProfile stage 2
            """
            norms = self.package_norms(Lspline)
            E1_0 = norms[0]
            E2_0 = norms[1]
            E3_0 = norms[2]
            S_0  = norms[3]
            curve = Lspline.curve
            #
            #**********************************************************************
            # initialize
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve)
            #
            #**********************************************************************
            # area  CPKeelProfile stage 2
            FPD.add_AreaConstraint(kind='equality',
                                   value = area)
            #
            #**********************************************************************
            # Xc CPKeelProfile
            if Xc is None:
                FPD.add_XcConstraint(kind='LS', 
                                     value=self.LCG)
            else:
                print 'Centerplane setting Xc constraint to ',Xc
                FPD.add_XcConstraint(kind='equality', value=Xc)
            #
            #**********************************************************************
            # cL Keel Profile - stay within parameters
            FPD.add_xVertexConstraint(kind='min',
                                     index = 2,
                                     value = 0.)
            FPD.add_xVertexConstraint(kind='max',
                                     index = 8,
                                     value = 1.05*self.LengthDWL)
            FPD.add_xVertexConstraint(kind='max',
                                     index = 9,
                                     value = 1.05*self.LengthDWL)
            #
            #******************************************************************
            #  CPKeelProfile stage 2
            # end condition stays nice:
            #            FPD.add_yVertexConstraint(kind='max',
            #                                     index = 1,
            #                                     value = Dmax*1.3)
            # ensure no 0 crossing
            FPD.add_yVertexConstraint(kind='min',
                                     index = 1,
                                     value = 0.)
            FPD.add_yVertexConstraint(kind='min',
                                     index = 2,
                                     value = 0.)
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 2,
            #                         value = Dmax*1.5)
            FPD.add_yVertexConstraint(kind='min',
                                     index = 8,
                                     value = 0.)
            FPD.add_yVertexConstraint(kind='min',
                                     index = 9,
                                     value = 0.)
            #
            #******************************************************************
            # x-flat CPkeelProfile  (4,5,6 + optional:  keep 3,7 out)
            # stage 1 solve
            #        FPD.add_xVertexConstraint(kind='equality',
            #                                 index = 3,
            #                                 value = self.stations.FOCP[0])
            #        FPD.add_xVertexConstraint(kind='equality',
            #                                 index = 7,
            #                                 value = self.stations.FOCP[1])
            FPD.add_xVertexConstraint(kind='max',
                                     index = 3,
                                     value = self.stations.FOCP[0])
            FPD.add_xVertexConstraint(kind='equality',
                                     index = 4,
                                     value = self.stations.FOCP[0])
            FPD.add_xVertexConstraint(kind='min',
                                     index = 5,
                                     value = self.stations.FOCP[0])
            FPD.add_xVertexConstraint(kind='max',
                                     index = 5,
                                     value = self.stations.FOCP[1])
            FPD.add_xVertexConstraint(kind='equality',
                                     index = 6,
                                     value = self.stations.FOCP[1])
            FPD.add_xVertexConstraint(kind='min',
                                     index = 7,
                                     value = self.stations.FOCP[1])
            #
            #**********************************************************************
            # y-flat CPkeelProfile (3,4,5,6,7) stage 1 solve
            FPD.add_yVertexConstraint(kind='equality',
                                     index = 3,
                                     value = Dmax)
            FPD.add_yVertexConstraint(kind='equality',
                                     index = 4,
                                     value = Dmax)
            FPD.add_yVertexConstraint(kind='equality',
                                     index = 6,
                                     value = Dmax)
            FPD.add_yVertexConstraint(kind='equality',
                                     index = 7,
                                     value = Dmax)
            FPD.add_yVertexConstraint(kind='equality',
                                     index = 5,
                                     value = Dmax)
            #
            #**********************************************************************
            # DWL fairness norms
            FPD.add_E1(kind='LS', weight = .5/E1_0.value)
            FPD.add_E2(kind='LS', weight = .5/E2_0.value)
            FPD.add_E3(kind='LS', weight = .5/E3_0.value)
            FPD.add_ArcLengthApprox(kind='LS', weight = .5/S_0)
            #
            #**********************************************************************
            #Setup and Solve:
            #
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            print 'Doing 2nd stage CPkeel solve'
            Lspline.optimize()
            return Lspline
        
        #
        #**********************************************************************
        # save init state
        inicurve = copy.copy(curve)
        #
        #**********************************************************************
        # First solve:
        Lspline = stage1(curve)
        #
        #**********************************************************************
        # Second solve:
        Lspline = stage2(Lspline)
        #
        #**********************************************************************
        #check solution:
        #        if Lspline.error_code is 'max_iter':
        #            pass 
        #            #Lspline = stage_fix(Lspline)
        #        elif Lspline.error_code is 'NAN':
        #            Lspline = stage1(curve)
        #            Lspline = stage_fix(Lspline)=
        #
        #**********************************************************************
        # save CLKeelProfile
        #Lspline.curve.knot_insertion(.5) #get to 11
        self.Lsplines.CProfile=Lspline
        self.CProfile2D = Lspline.curve
        #
        #**********************************************************************
        # lift to 3D
        self.CProfile = self.lifting_2DCProfile(Lspline.curve.vertices, 
                                                elevation,k,nump)
        return
    
class MidShip(object):
    def __init__(self):
        pass
    
class BowFairness(object):
    def __init__(self):
        pass
    
class SternFairness(object):
    def __init__(self):
        pass
    
    
class CPKfwd(object):
    def __init__(self):
        pass
    
    
class CPKsplitcurves(object):
    """
    """
    def __init__(self, 
                 initial_curve,
                 LazyFPD,
                 ordered_lazy_solvercontainer,
                 hull=None):
        
        self.initial_curve = initial_curve
        self.o_solvers = ordered_lazy_solvercontainer
        self.LazyFPD = LazyFPD
        self.hull = hull
        pass
    
    
    #
    #**********************************************************************
    # FWD flat finding
    #**********************************************************************
    #
    def compute_max_flat_fwd_extents(self, 
                        curve, maxDraft,
                        tol=1.e-7,nump=100,
                        ini=0.,fini=1.):
        """This funny routine is needed to find
        -THE- parametric location of maximum
        foward extents of the flat spot of the cLKeelProfile curve
        (This is the kind of thing that seems needed to -really-
        finish a hull solver which is faithful to the constraints
        which it began with)
        """
        param=.5
        pvector = np.linspace(ini,fini,nump)
        i_p_x_y = []
        
        minflatloc = curve.vertices[-1,0]
        maxflatloc = curve.vertices[0,0]
        
        for i,p in enumerate(pvector):
            val = curve.CurvePoint(p)
            if np.linalg.norm(val[1]-maxDraft)<tol:
                if val[0]<minflatloc:
                    i_p_x_y.append((i,p,val))
                    minflatloc = min(minflatloc,val[0])
                #if val[0]>maxflatloc:
                #    maxflatloc = max(maxflatloc,val[0])
                
        return i_p_x_y, pvector[i_p_x_y[0][0]-1:i_p_x_y[0][0]+1]
    
    
    
    #
    #**********************************************************************
    # FWD flat finding
    #**********************************************************************
    #
    def find_xy_max_flat_fwd_extents(self, curve, maxDraft,
                        tol=1.e-15,nump=100,maxiter=100,
                        ini=0.,fini=1.):
        """This funny routine is needed to find
        -THE- parametric location of maximum
        foward extents of the flat spot of the cLKeelProfile curve
        (This is the kind of thing that seems needed to -really-
        finish a hull solver which is faithful to the constraints
        which it began with)
        """
        #robust(very accurate) starting value:
        i_p_x_y, range_next = self.compute_max_flat_fwd_extents( 
                                                    curve, maxDraft,
                                                    tol=1.e-15, nump=100,
                                                    ini=0., fini=1.)
        
        #binary search for the earliest availible point
        # which is at max draft
        mid       = sum(range_next)*.5
        
        low = range_next[0]
        high = range_next[1]
        
        
        lowcv = curve.CurvePoint(low)[1]
        highcv = curve.CurvePoint(high)[1]
        q = curve.CurvePoint(mid)[1]
        
        ans_tol = np.linalg.norm(maxDraft-lowcv)
        count = 0
        while (ans_tol>=tol) and count<maxiter :
            
            if q<maxDraft:
                low = mid
            else:
                high = mid
            mid = (low+high)*.5
            
            q = curve.CurvePoint(mid)[1]
            lowcv = curve.CurvePoint(low)[1]
            
            ans_tol = np.linalg.norm(q-lowcv)
            
            count+=1
        
        final_param_val = high
        if count >= maxiter:
            print 'error, point not found'
            found = False
            return final_param_val, q, found
        else:
            found = True
            q = curve.CurvePoint(mid)
            return final_param_val, q, found  
    #
    #**********************************************************************
    # AFT flat finding
    #**********************************************************************
    #
    
    def compute_max_flat_aft_extents(self, 
                        curve, maxDraft,
                        tol=1.e-7,nump=100,
                        ini=0.,fini=1.):
        """This funny routine is needed to find
        -THE- parametric location of maximum
        foward extents of the flat spot of the cLKeelProfile curve
        (This is the kind of thing that seems needed to -really-
        finish a hull solver which is faithful to the constraints
        which it began with)
        
        TODO: obviously it is more efficient to use one max
        flat finder for both fore and aft.
        Unfortunately programmer time is at a all time premium.  :-(
        
        
        SLIGHT BUG IS BEING LEFT IN HERE:
            this returns the edge of the array of values meeting
            the criteria.
            -we want the last value meeting the criteria
            and the 'first' value violating it instead
            
        Furthermore, naive application of this funciton 
        on the reduced range+outside bounds
        seems like it might outpreform the binary
        search in testing..
        --at any rate it would be better to do 
        CPKeel as a 3 part curve.
        """
        param=.5
        pvector = np.linspace(ini,fini,nump)
        i_p_x_y = []
        
        minflatloc = curve.vertices[-1,0]
        maxflatloc = curve.vertices[0,0]
        
        for i,p in enumerate(pvector):
            val = curve.CurvePoint(p)
            if np.linalg.norm(val[1]-maxDraft)<tol:
                #if val[0]<minflatloc:
                #    i_p_x_y.append((i,p,val))
                #    minflatloc = min(minflatloc,val[0])
                if val[0]>maxflatloc:
                    i_p_x_y.append((i,p,val))
                    maxflatloc = max(maxflatloc,val[0])
                
        return i_p_x_y,pvector[i_p_x_y[-1][0]-1:i_p_x_y[-1][0]+1]
    
    
    
    def find_xy_max_flat_aft_extents(self, curve, maxDraft,
                        tol=1.e-15,nump=100,maxiter=100,
                        ini=0.,fini=1.):
        """This funny routine is needed to find
        -THE- parametric location of maximum
        foward extents of the flat spot of the cLKeelProfile curve
        (This is the kind of thing that seems needed to -really-
        finish a hull solver which is faithful to the constraints
        which it began with)
        
        
        dev:
            
            curve=self.CProfile2D
            maxDraft = self.MaxDraft
            
          i_p_x_y, range_next   =      cLKP_solver.compute_max_flat_aft_extents( 
                                                    curve, maxDraft,
                                                    tol=1.e-15, nump=100,
                                                    ini=0., fini=1.)
        """
        #robust(very accurate) starting value:
        i_p_x_y, range_next = self.compute_max_flat_aft_extents( 
                                                    curve, maxDraft,
                                                    tol=1.e-15, nump=100,
                                                    ini=0., fini=1.)
        
        #binary search for the earliest availible point
        # which is at max draft
        mid       = sum(range_next)*.5
        
        low = range_next[0]
        high = range_next[1]
        
        
        lowcv = curve.CurvePoint(low)[1]
        highcv = curve.CurvePoint(high)[1]
        q = curve.CurvePoint(mid)[1]
        
        ans_tol = np.linalg.norm(maxDraft-lowcv)
        count = 0
        while (ans_tol>=tol) and count<maxiter :
            
            if q<maxDraft:
                high = mid
            else:
                low = mid
            mid = (low+high)*.5
            
            q = curve.CurvePoint(mid)[1]
            lowcv = curve.CurvePoint(low)[1]
            
            ans_tol = np.linalg.norm(q-lowcv)
            
            count+=1
        
        final_param_val = high
        if count >= maxiter:
            print 'error, point not found'
            found = False
            return final_param_val, q, found
        else:
            found = True
            q = curve.CurvePoint(mid)
            return final_param_val, q, found  
        
        
        
    #
    #**********************************************************************
    # AFT optimize
    #**********************************************************************
    #
    
#    def find_good_solve(self,
#                   curve,
#                   LazyFPD=None,
#                   return_Lspline=True,
#                   make_thb = True,
#                   ncpts=self.mshp_sections[0].n,
#                   wt=10):
#        
#        return
    

class DWLsplitcurves(object):
    def __init__(self, 
                 initial_curve,
                 LazyFPD,
                 ordered_lazy_solvercontainer,
                 hull=None):
        
        self.initial_curve = initial_curve
        self.o_solvers = ordered_lazy_solvercontainer
        self.LazyFPD = LazyFPD
        self.hull = hull
        pass
    
    
    
    #
    #**********************************************************************
    # FWD flat finding
    #**********************************************************************
    #
    def compute_max_flat_fwd_extents(self, 
                        curve, maxDraft,
                        tol=1.e-7,nump=100,
                        ini=0.,fini=1.):
        """This funny routine is needed to find
        -THE- parametric location of maximum
        foward extents of the flat spot of the cLKeelProfile curve
        (This is the kind of thing that seems needed to -really-
        finish a hull solver which is faithful to the constraints
        which it began with)
        """
        param=.5
        pvector = np.linspace(ini,fini,nump)
        i_p_x_y = []
        
        minflatloc = curve.vertices[-1,0]
        maxflatloc = curve.vertices[0,0]
        
        for i,p in enumerate(pvector):
            val = curve.CurvePoint(p)
            if np.linalg.norm(val[1]-maxDraft)<tol:
                if val[0]<minflatloc:
                    i_p_x_y.append((i,p,val))
                    minflatloc = min(minflatloc,val[0])
                #if val[0]>maxflatloc:
                #    maxflatloc = max(maxflatloc,val[0])
                
        return i_p_x_y, pvector[i_p_x_y[0][0]-1:i_p_x_y[0][0]+1]
    
    
    
    #
    #**********************************************************************
    # FWD flat finding
    #**********************************************************************
    #
    def find_xy_max_flat_fwd_extents(self, curve, maxDraft,
                        tol=1.e-15,nump=100,maxiter=100,
                        ini=0.,fini=1.):
        """This funny routine is needed to find
        -THE- parametric location of maximum
        foward extents of the flat spot of the ___DWL___ curve
        (This is the kind of thing that seems needed to -really-
        finish a hull solver which is faithful to the constraints
        which it began with)
        """
        #robust(very accurate) starting value:
        i_p_x_y, range_next = self.compute_max_flat_fwd_extents( 
                                                    curve, maxDraft,
                                                    tol=1.e-15, nump=100,
                                                    ini=0., fini=1.)
        
        #binary search for the earliest availible point
        # which is at max draft
        mid       = sum(range_next)*.5
        
        low = range_next[0]
        high = range_next[1]
        
        
        lowcv = curve.CurvePoint(low)[1]
        highcv = curve.CurvePoint(high)[1]
        q = curve.CurvePoint(mid)[1]
        
        ans_tol = np.linalg.norm(maxDraft-lowcv)
        count = 0
        while (ans_tol>=tol) and count<maxiter :
            
            if q<maxDraft:
                low = mid
            else:
                high = mid
            mid = (low+high)*.5
            
            q = curve.CurvePoint(mid)[1]
            lowcv = curve.CurvePoint(low)[1]
            
            ans_tol = np.linalg.norm(q-lowcv)
            
            count+=1
        
        final_param_val = high
        if count >= maxiter:
            print 'error, point not found'
            found = False
            return final_param_val, q, found
        else:
            found = True
            q = curve.CurvePoint(mid)
            return final_param_val, q, found  
    #
    #**********************************************************************
    # AFT flat finding
    #**********************************************************************
    #
    
    def compute_max_flat_aft_extents(self, 
                        curve, maxDraft,
                        tol=1.e-7,nump=100,
                        ini=0.,fini=1.):
        """This funny routine is needed to find
        -THE- parametric location of maximum
        foward extents of the flat spot of the ___DWL___ curve
        (This is the kind of thing that seems needed to -really-
        finish a hull solver which is faithful to the constraints
        which it began with)
        
        TODO: obviously it is more efficient to use one max
        flat finder for both fore and aft.
        Unfortunately programmer time is at a all time premium.  :-(
        
        
        SLIGHT BUG IS BEING LEFT IN HERE:
            this returns the edge of the array of values meeting
            the criteria.
            -we want the last value meeting the criteria
            and the 'first' value violating it instead
            
        Furthermore, naive application of this funciton 
        on the reduced range+outside bounds
        seems like it might outpreform the binary
        search in testing..
        --at any rate it would be better to do 
        ___DWL___ as a 3 part curve.
        """
        param=.5
        pvector = np.linspace(ini,fini,nump)
        i_p_x_y = []
        
        minflatloc = curve.vertices[-1,0]
        maxflatloc = curve.vertices[0,0]
        
        for i,p in enumerate(pvector):
            val = curve.CurvePoint(p)
            if np.linalg.norm(val[1]-maxDraft)<tol:
                #if val[0]<minflatloc:
                #    i_p_x_y.append((i,p,val))
                #    minflatloc = min(minflatloc,val[0])
                if val[0]>maxflatloc:
                    i_p_x_y.append((i,p,val))
                    maxflatloc = max(maxflatloc,val[0])
                
        return i_p_x_y,pvector[i_p_x_y[-1][0]-1:i_p_x_y[-1][0]+1]
    
    
    
    def find_xy_max_flat_aft_extents(self, curve, maxDraft,
                        tol=1.e-15,nump=100,maxiter=100,
                        ini=0.,fini=1.):
        """This funny routine is needed to find
        -THE- parametric location of maximum
        foward extents of the flat spot of the ___DWL___ curve
        (This is the kind of thing that seems needed to -really-
        finish a hull solver which is faithful to the constraints
        which it began with)
        
        
        dev:
            
            curve=self.CProfile2D
            maxDraft = self.MaxDraft
            
          i_p_x_y, range_next   =      cLKP_solver.compute_max_flat_aft_extents( 
                                                    curve, maxDraft,
                                                    tol=1.e-15, nump=100,
                                                    ini=0., fini=1.)
        """
        #robust(very accurate) starting value:
        i_p_x_y, range_next = self.compute_max_flat_aft_extents( 
                                                    curve, maxDraft,
                                                    tol=1.e-15, nump=100,
                                                    ini=0., fini=1.)
        
        #binary search for the earliest availible point
        # which is at max draft
        mid       = sum(range_next)*.5
        
        low = range_next[0]
        high = range_next[1]
        
        
        lowcv = curve.CurvePoint(low)[1]
        highcv = curve.CurvePoint(high)[1]
        q = curve.CurvePoint(mid)[1]
        
        ans_tol = np.linalg.norm(maxDraft-lowcv)
        count = 0
        while (ans_tol>=tol) and count<maxiter :
            
            if q<maxDraft:
                high = mid
            else:
                low = mid
            mid = (low+high)*.5
            
            q = curve.CurvePoint(mid)[1]
            lowcv = curve.CurvePoint(low)[1]
            
            ans_tol = np.linalg.norm(q-lowcv)
            
            count+=1
        
        final_param_val = high
        if count >= maxiter:
            print 'error, point not found'
            found = False
            return final_param_val, q, found
        else:
            found = True
            q = curve.CurvePoint(mid)
            return final_param_val, q, found  
        
        