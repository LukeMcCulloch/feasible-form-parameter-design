#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 18:17:57 2018

@author: luke

3 Piece Fwd Complex Hull : Bare Fwd Surface + Bulb

Output:  Chiefly the 5 Rhino surfaces

surface0data.txt
surface1data.txt
surface2data.txt
surface3data.txt
surface4data.txt

"""

from THBLSplineHullcomplex import *
#substitute the usual formalities to just pretend everything from there is here.

from iaBox import BoxList, Box
import os #check if file exists

class sketch(object):
    """assumes that you have already 
    *refined the fwd bare hull surface
       to a level that will have the bottom surface tcurve vertices
    line up nicely
       with the bulb tcurve vertices.
    """
    def __init__(self, bbowcurve, hullfwdcurve, hullaftcurve ):
        """
        bbowcurve :     bbow fwd transverse
        
        hullcurffve:    curve- hullfwd transverse fwd
        
        aftcurve:       curve- hullfwd transverse aft
        """
        self.bbowcurve = bbowcurve
        self.hullfwdcurve = hullfwdcurve
        self.hullaftcurve = hullaftcurve
        self.bulbn = bbowcurve.n
        self.hulln = hullfwdcurve.n-bbowcurve.n
        mssg = 'ERROR: fore and aft transverse curves of the fwd refined hull dont match'
        assert(hullfwdcurve.n == hullaftcurve.n),mssg
        vertices = []
        #----------------------------------------
        #bulb vertices at fwd edge
        for el in bbowcurve.vertices:
            vertices.append(el)
        #hull vertices above them, hopefully
        for el in hullfwdcurve.vertices[self.bulbn:]:
            vertices.append(el)
        self.fwd_vertices = vertices
        #----------------------------------------
        self.aft_vertices = self.hullaftcurve.vertices
    
    
    def refine(self):
        """
        """
        bbowcurve = self.bbowcurve.dyadic_refinement()
        hullfwdcurve = self.hullfwdcurve.dyadic_refinement()
        hullaftcurve = self.hullaftcurve.dyadic_refinement()
        return sketch(bbowcurve, hullfwdcurve, hullaftcurve)




class make_upper(object):
    """The trick here is 
    5 transverse vertices are incoming since 7+4 = 11 
    but we have one shared boundary!
    """
    def __init__(self, tcurves, bbow, lower):
        self.k = 4
        self.nump=30
        self.nvi = 4
        self.nvs = 7
        self.tcurves = tcurves
        self.ntcurves = len(self.tcurves)
        self.nvcurves = self.tcurves[0].n
        assert(self.nvcurves == 7),'ERROR: there should be 7 transverse vertices'
        self.bbow = bbow
        self.maxlevel_thbsolver = 3
        self.init_lcurves_solver()
        self.lower = lower
        
        
    def init_lcurves_solver(self):
        """
        class:  make_upper
        
        refine lcurves so as to 
        interpolate each vertex on the original bare hull
        transverse curves
        and each vertex on the bulb transverse curves
        """
        self.lcurves = []
        self.solver_records = {}
        for i in range(self.tcurves[0].n):
            ini = self.tcurves[0].vertices[i]
            fin = self.tcurves[-1].vertices[i]
            vtx = linear_vertices(ini,fin,self.nvi)
            curve = rbspline(vtx, k=self.k, nump=self.nump)
            while curve.n<11:
                curve = curve.dyadic_refinement()
            curve.parent = None
            self.solver_records[i] = records(curve)
            self.solver_records[i].vtinterp = []
            #
            #******************************************************************
            # no bbow in the upper portion
            for curve in self.tcurves[1:self.ntcurves-1]:
                vertex = curve.vertices[i]
                self.solver_records[i].vtinterp.append(vertex)
        return
        
        
        
        
    def make_lcurves(self):
        """
        class:  make_upper
        """
        self.lcurves = []
        for i in range(self.tcurves[0].n):
            print '------------------------------------------'
            print '\n\nBegin Solving New Longitudinal({}):'.format(i)
            print ' On the Upper Hull \n\n'
            print '------------------------------------------'
            curve = copy.deepcopy(self.solver_records[i].curve)
            verts = self.solver_records[i].vtinterp
            
            
            Lspline = self.build_lagrangian(thbcurve = curve,
                                            vtinterp = verts,
                                            ith=i)
            
            self.solver_records[i].Lspline = Lspline
            
                        
            Lspline = uopt.THBsolver(thbcurve = curve,
                                     Lspline = Lspline,
                                     maxlevel = self.maxlevel_thbsolver,
                                     normalize = False)
        
            self.lcurves.append(Lspline.curve)
        return
    
    
    
    def build_lagrangian(self, 
                         thbcurve, 
                         vtinterp,
                         ith,
                         this_kind='equality'):
        """
        class:  make_upper
        TODO:
            -fwd 5 longitudinals fixed
            -fore and aft curveture continuity to bulb and midship
            -interpolate 2 interior bare hull vertices
        """
        #other_kind = 'LS'
        #this_kind = 'equality'
        big = 100000.
        mig = big/2.
        interval_data, small = interval_bounds(thbcurve)
        FPD = FormParameterDict(thbcurve)
        
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .05)
        FPD.add_ArcLengthApprox(kind='LS',
                              weight = mig)
        #FPD.add_ArcLength(kind='LS',
        #                      weight = mig)
        #
        #**********************************************************************
        # specify where interpolation will take place:
        ntinterp = len(vtinterp)
        ukbar = np.linspace(0.,1.,ntinterp+2)[1:-1]
        #
        #**********************************************************************
        # fair ini
        xb = thbcurve.vertices[0,0]
        yb = thbcurve.vertices[0,1]
        zb = thbcurve.vertices[0,2]
        # fair fini
        xe = thbcurve.vertices[-1,0]
        ye = thbcurve.vertices[-1,1]
        ze = thbcurve.vertices[-1,2]
        #**********************************************************************
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
        #*********************************************
        #*********************************************
        # nose tangent - Rhino style?
#        FPD.add_yFixity(index = 1,
#                        value = yb,
#                        track_prolongation=False)
#        #--------------------------
#        FPD.add_zFixity(index = 1,
#                        value =zb,
#                        track_prolongation=False)
        #if ith == 0:# or ith == 1:
        #if ith == 0 or ith == 1: #july 6 2018 attempt to smooth new hull style
        #if ith == 0 or \
        #    ith == 1 or \
        #    ith == 2: #july 6 2018 attempt to smooth new hull style
        if False:
            #then we are on the bottom of the upper - make it like 
            # the lower so as the ease the work of the THB rep later.
            
            #*********************************************
            # C2 continuity to bulb.  (bottom of upper)
            FPD.add_xFixity(index = 1,
                            value = xb,
                            track_prolongation=False)
            #FPD.add_xFixity(index = 2,
            #                value = xb,
            #                track_prolongation=False)
            #--------------------------
            FPD.add_yFixity(index = 1,
                            value = yb,
                            track_prolongation=False)
            #FPD.add_yFixity(index = 2,
            #                value = yb,
            #                track_prolongation=False)
        elif ith == 0:
            #*********************************************
            # nose tangent
            FPD.add_yFixity(index = 1,
                            value = yb,
                            track_prolongation=False)
            #--------------------------
            FPD.add_zFixity(index = 1,
                            value =zb,
                            track_prolongation=False)
            
            FPD.add_xVertexConstraint(kind='LS',
                                      value=xb,
                                      weight=1.,
                                      index=1)
            
        else:
            #*********************************************
            # nose tangent
            FPD.add_yFixity(index = 1,
                            value = yb,
                            track_prolongation=False)
            #--------------------------
            FPD.add_zFixity(index = 1,
                            value =zb,
                            track_prolongation=False)
            #if ith == self.tcurves[0].n-1:
            #    
            #    #*********************************************
            #    # nose curvature
            #    FPD.add_yFixity(index = 2,
            #                    value = yb,
            #                    track_prolongation=False)
            #    #--------------------------
            #    FPD.add_zFixity(index = 2,
            #                    value =zb,
            #                    track_prolongation=False)
        #*********************************************
        # C1 continuity to midship.
        FPD.add_xFixity(index = -2,
                        value = xe,
                        track_prolongation=False)
        FPD.add_yFixity(index = -2,
                        value = ye,
                        track_prolongation=False)
        #
        #*********************************************
        # C2 continuity to midship
        """idea:  activate, say the ye Fixity here
        to try and smooth out the buttock lines...?
        and tiny lump induced by making 
        transition to midship 
        totally C1 ??
        TRY this 5-27-2018 ?
        """
#        FPD.add_xFixity(index = -3,
#                        value = xe,
#                        track_prolongation=False)
#        FPD.add_yFixity(index = -3,
#                        value = ye,
#                        track_prolongation=False)
        
            
#                #*********************************************
#                # nose tangent/curvature
#                #edit: just c1 to nose.
#                FPD.add_yFixity(index = 1,
#                                value = yb,
#                                track_prolongation=False)
#                #FPD.add_yFixity(index = 2,
#                #                value = yb,
#                #                track_prolongation=False)
#                #--------------------------
#                FPD.add_zFixity(index = 1,
#                                value =zb,
#                                track_prolongation=False)
#                #FPD.add_zFixity(index = 2,
#                #                value = zb,
#                #                track_prolongation=False)
#            #elif ith >= self.tcurves[0].n-2:
#            else:
#                #*********************************************
#                # nose tangent
#                FPD.add_yFixity(index = 1,
#                                value = yb,
#                                track_prolongation=False)
#                #--------------------------
#                FPD.add_zFixity(index = 1,
#                                value =zb,
#                                track_prolongation=False)
            
        #*********************************************
        # conflict with aft interpolation point
        # if that point is not on the true dwl flat edge.
        #-- which in a commercial vessel it looks like it should be.
        #        FPD.add_xFixity(index = -3,
        #                        value = xe,
        #                        track_prolongation=False)
        #        FPD.add_yFixity(index = -3,
        #                        value = ye,
        #                        track_prolongation=False)
        
#        FPD.add_CurvatureXoverZConstraint(  kind        = this_kind,
#                                            location    = 1.,
#                                            value       = 0.,
#                                            weight      = mig)
        mssg = 'Error:  possibly wrong fence posts in the Lagrangian setup'
        assert(count == ntinterp),mssg
        
        #
        #**********************************************************************
        #
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                            interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(thbcurve, L, 
                                         data = interval_data)
        return Lspline
    
    
    
    
    def make_surface(self, refine=False):
        self.surface = make_surface(self.lcurves,
                                    self.tcurves,
                                    refine=refine)
        return 
    
    
    
    def plot(self, canvas = None):
        """
        dev
        ----------
        ax = sfh.mklwr.solver_records[0].Lspline.curve.plotcurvehierarchy()
        
        kk = sfh.mklwr.solver_records.keys()
        for i in kk:
            ax = sfh.mklwr.solver_records[i].Lspline.curve.plotcurvehierarchy(
                        canvas=ax)
        """
        if canvas is None:
            ax = self.solver_records[0].Lspline.curve.plotcurvehierarchy()
        else:
            ax = canvas
        kk = self.solver_records.keys()
        for i in kk:
            ax = self.solver_records[i].Lspline.curve.plotcurvehierarchy(
                    canvas=ax)
        return ax
        
        
        

class make_lower(object):
    def __init__(self, tcurves, bbow):
        self.k = 4
        self.nump=30
        self.nvi = 4
        self.nvs = 11
        self.tcurves = tcurves
        self.ntcurves = len(self.tcurves)   #+nlower_tcurve 
                                            #if this name made sense.
        self.nvcurves = self.tcurves[0].n
        assert(self.nvcurves == 7),'ERROR: there should be 7 transverse vertices'
        self.bbow = bbow
        self.maxlevel_thbsolver = 3
        self.init_lcurves_solver()
        
        
    def init_lcurves_solver(self):
        """
        class:  make_lower
        
        refine lcurves so as to 
        interpolate each vertex on the original bare hull
        transverse curves
        and each vertex on the bulb transverse curves
        """
        self.lcurves = []
        self.solver_records = {}
        for i in range(self.tcurves[0].n):
            ini = self.tcurves[0].vertices[i]
            #ini = self.bbow.tcurvelist[-1].vertices[i]
            #ini = self.bbow.tcurvelist[0].vertices[i]
            fin = self.tcurves[-1].vertices[i]
            vtx = linear_vertices(ini,fin,self.nvi)
            curve = rbspline(vtx, k=self.k, nump=self.nump)
            while curve.n<11:
                curve = curve.dyadic_refinement()
            curve.parent = None
            self.solver_records[i] = records(curve)
            self.solver_records[i].vtinterp = []
            
            #turn off bbow approx.
            #for j in range(len(self.bbow.tcurvelist)):
            #    vertex = self.bbow.tcurvelist[-1-j].vertices[i]
            #    self.solver_records[i].vtinterp.append(vertex)
            
            
            #make the top-fwd of the bottom surface interpolate
            # the upper vertex there.
            # no need.  bottom surface
            # goes only to the aft edge of bulb
            #if i == self.bbow.lcurvelist[0].n:
            #    for curve in self.tcurves[0:self.ntcurves-1]:
            #        vertex = curve.vertices[i]
            #        self.solver_records[i].vtinterp.append(vertex)
            #else:
            #    for curve in self.tcurves[1:self.ntcurves-1]:
            #        vertex = curve.vertices[i]
            #        self.solver_records[i].vtinterp.append(vertex)
                    
            #
            for curve in self.tcurves[1:self.ntcurves-1]:
                vertex = curve.vertices[i]
                self.solver_records[i].vtinterp.append(vertex)
            
        return
        
        
    def make_lcurves(self):
        """
        class:  make_lower
        """
        self.lcurves = []
        for i in range(self.tcurves[0].n):
            print '------------------------------------------'
            print '\n\nBegin Solving New Longitudinal({}):'.format(i)
            print ' On the Lower Hull \n\n'
            print '------------------------------------------'
            curve = self.solver_records[i].curve
            verts = self.solver_records[i].vtinterp
            
            
            Lspline = self.build_lagrangian(thbcurve = curve,
                                            vtinterp = verts,
                                            ith=i)
            self.solver_records[i].Lspline = Lspline
            
                        
            Lspline = uopt.THBsolver(thbcurve = curve,
                                     Lspline = Lspline,
                                     maxlevel = self.maxlevel_thbsolver,
                                     normalize = False)
        
            self.lcurves.append(Lspline.curve)
        return
    
    
    
    
    def build_lagrangian(self, 
                         thbcurve, 
                         vtinterp,
                         ith,
                         this_kind='equality'):
        """
        class:  make_lower
        
        TODO:
            -fwd 5 longitudinals fixed
            -fore and aft curveture continuity to bulb and midship
            -interpolate 2 interior bare hull vertices
            
            
            thbcurve = curve
            vtinterp = verts
        """
        #this_kind = 'equality'
        #other_kind = 'LS'
        big = 100000.
        mig = big/2.
        interval_data, small = interval_bounds(thbcurve)
        FPD = FormParameterDict(thbcurve)
        
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .05)
        FPD.add_ArcLengthApprox(kind='LS',
                              weight = mig)
        #FPD.add_ArcLength(kind='LS',
        #                      weight = mig)
        #
        #**********************************************************************
        # specify where interpolation will take place:
        ntinterp = len(vtinterp)
        ukbar = np.linspace(0.,1.,ntinterp+2)[1:-1]
        #
        #**********************************************************************
        #
        #**********************************************************************
        # fair ini
        xb = thbcurve.vertices[0,0]
        yb = thbcurve.vertices[0,1]
        zb = thbcurve.vertices[0,2]
        # fair fini
        xe = thbcurve.vertices[-1,0]
        ye = thbcurve.vertices[-1,1]
        ze = thbcurve.vertices[-1,2]
        #
        #**********************************************************************
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
        #**********************************************************************
                # stuff to maintain continuity from hull to bulb...
                # NOW replaced by fixity below.
        #        FPD.add_AngleConstraint(kind        = this_kind,
        #                                location    = loc,
        #                                value       = 0.,
        #                                weight      = mig)
        #        FPD.add_CurvatureXoverZConstraint(  kind        = this_kind,
        #                                            location    = loc,
        #                                            value       = 0.,
        #                                            weight      = mig)
                
        #        FPD.add_xPointConstraint(kind= this_kind, 
        #                                 value=pt[0], 
        #                                 location=loc+loc/10., 
        #                                 weight = big ) 
        #        FPD.add_yPointConstraint(kind= this_kind, 
        #                                 value=pt[1], 
        #                                 location=loc+loc/10., 
        #                                 weight = big  )
        #        FPD.add_xPointConstraint(kind= this_kind, 
        #                                 value=pt[0], 
        #                                 location=loc+loc/5., 
        #                                 weight = big ) 
        #        FPD.add_yPointConstraint(kind= this_kind, 
        #                                 value=pt[1], 
        #                                 location=loc+loc/5., 
        #                                 weight = big  )
                #
        #**********************************************************************
        #
        #        if ith == self.nvcurves:
        #            #then we are on the top of the bulb - make it like 
        #            # the upper?
        #            # - can't really, as it would be about 5 indices inwards from the nose.
        #            #*********************************************
        #            pass
        #            # nose tangent/curvature
        #            #            FPD.add_yFixity(index = 1,
        #            #                            value = yb,
        #            #                            track_prolongation=False)
        #            #            FPD.add_yFixity(index = 2,
        #            #                            value = yb,
        #            #                            track_prolongation=False)
        #            #            #--------------------------
        #            #            FPD.add_zFixity(index = 1,
        #            #                            value =zb,
        #            #                            track_prolongation=False)
        #            #            FPD.add_zFixity(index = 2,
        #            #                            value = zb,
        #            #                            track_prolongation=False)
        #*********************************************
        # 
#        if ith == self.tcurves[0].n-1:
#            # nose tangent! - rhino attempt?
#            FPD.add_yFixity(index = 1,
#                            value = yb,
#                            track_prolongation=False)
#            #--------------------------
#            FPD.add_zFixity(index = 1,
#                            value =zb,
#                            track_prolongation=False)
#        else:
        #*********************************************
        # C1 or C2 continuity to bulb.
        if ith < self.nvcurves-2:#2 best so far.  #tried 1 and 3 as well
            FPD.add_xFixity(index = 1,
                            value = xb,
                            track_prolongation=False)
            #FPD.add_xFixity(index = 2,
            #                value = xb,
            #                track_prolongation=False)
            #--------------------------
            FPD.add_yFixity(index = 1,
                            value = yb,
                            track_prolongation=False)
            #FPD.add_yFixity(index = 2,
            #                value = yb,
            #                track_prolongation=False)
        else:
            #*********************************************
            # nose tangent - flush with upper surface
            FPD.add_yFixity(index = 1,
                            value = yb,
                            track_prolongation=False)
            #--------------------------
            FPD.add_zFixity(index = 1,
                            value =zb,
                            track_prolongation=False)
            #
            FPD.add_xVertexConstraint(kind='LS',
                                      value=xb,
                                      weight=1.,
                                      index=1)
                
        #*********************************************
        # C1 continuity to midship.
        FPD.add_xFixity(index = -2,
                        value = xe,
                        track_prolongation=False)
        FPD.add_yFixity(index = -2,
                        value = ye,
                        track_prolongation=False)
        #
        #*********************************************
        # C2 continuity to midship
        """idea:  activate, say the ye Fixity here
        to try and smooth out the buttock lines...?
        and tiny lump induced by making 
        transition to midship 
        totally C1 ??
        TRY this 5-27-2018 ?
        """
#        FPD.add_xFixity(index = -3,
#                        value = xe,
#                        track_prolongation=False)
#        FPD.add_yFixity(index = -3,
#                        value = ye,
#                        track_prolongation=False)
        mssg = 'Error:  possibly wrong fence posts in the Lagrangian setup'
        assert(count == ntinterp),mssg
        
        #
        #**********************************************************************
        #
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                            interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(thbcurve, L, 
                                         data = interval_data)
        return Lspline
    
    def make_surface(self, refine=False):
        self.surface = make_surface(self.lcurves,
                                    self.tcurves,
                                    refine=refine)
        return 
    
    def plot(self, canvas = None):
        """
        dev
        ----------
        ax = sfh.mklwr.solver_records[0].Lspline.curve.plotcurvehierarchy()
        
        kk = sfh.mklwr.solver_records.keys()
        for i in kk:
            ax = sfh.mklwr.solver_records[i].Lspline.curve.plotcurvehierarchy(
                        canvas=ax)
        """
        if canvas is None:
            ax = self.solver_records[0].Lspline.curve.plotcurvehierarchy()
        else:
            ax = canvas
        kk = self.solver_records.keys()
        for i in kk:
            ax = self.solver_records[i].Lspline.curve.plotcurvehierarchy(
                    canvas=ax)
        return ax
        
    
    
    
    
    

class subdivide_fwdhull(object):
    """Assuming that an 11 vertex transverse curve is
    a good fit for the 2 surface bulb incorporation.
    
    That is, we are not worying about 'really short'
    or especially 'really tall' bulbs - which might distort 
    the fitting with the 11 vertices on the transverse
    i.e. (7 vertices to fit the bulb, 4 vertices remain for the upper hull)
    """
    def __init__(self, tcurves, bbow):
        self.k = 4
        self.nump=30
        self.critical_num=11
        self.offset = 0.
        newtcurves = []
        for curve in tcurves:
            curve.children = None
            while curve.n<self.critical_num:
                curve = curve.dyadic_refinement()
            curve.parent = None
            newtcurves.append(curve)
        self.tcurves = newtcurves
        self.bbow = bbow
        self.enlarge_tcurves()
        self.move_bulb_to_bare_hull()
        self.split_tcurves()
        #
        self.mklwr = make_lower(self.lower, self.bbow)
        #
        self.mkupr = make_upper(self.upper, self.bbow,
                                self.mklwr.lcurves)
        
    def check_boundary(self):
        """just checks that the boundary top and bottom
        surfaces are using the same tcurves vertices
        at the interface
        """
        for i in [0,1,2,3]:
            ans = np.linalg.norm(self.mkupr.tcurves[0].vertices[0] - \
                                 self.mklwr.tcurves[0].vertices[-1])
            assert(ans<1.e-15),'Boundary Tcurve Difference Violated at {}'.format(i)
            
        return
    
    def enlarge_tcurves(self):
        """
        1 stage transverse refinement 
        to
            11 vertices a piece
        to
            try and ensure ( 7 , 7 ) surfaces ( low , high )
        """
        for curve in self.tcurves:
            curve = curve.dyadic_refinement()
        return
    
    
    def split_tcurves(self, p=.5):
        """p=.5
        gives a 7 and 7 split with a dyadic 7 knot vector to each curve
        
        e.g.
        cv_top.t
            Out[68]: array([0.  , 0.  , 0.  , 0.  , 0.25, 0.5 , 0.75, 1.  , 1.  , 1.  , 1.  ])
            
        cv_bottom.t
            Out[69]: array([0.  , 0.  , 0.  , 0.  , 0.25, 0.5 , 0.75, 1.  , 1.  , 1.  , 1.  ])
        """
        self.upper = []
        self.lower = []
        #
        #*******************************************************
        # fwd hull surface aft of the nose:
        for curve in self.tcurves[1:]:
            cv_bottom, cv_top = curve.CurveSplit(p)
            self.upper.append(cv_top)
            self.lower.append(cv_bottom)
            #
        #*******************************************************
        # first cuve, top fwd surface
        ini = self.bbow.tcurvelist[0].vertices[-1]
        #
        # handle offset of fwd top nose curve:
        ini[2] += self.offset
        fin = self.tcurves[0].vertices[-1]
        # make:
        vertices = linear_vertices(ini,
                                   fin,
                                   cv_top.n)
        topcurve = rbspline(vertices, k=self.k, nump=self.nump)
        self.upper.insert(0,topcurve)
        #
        #*******************************************************
        # first cuve, bottom fwd surface
        vertices = copy.deepcopy(self.bbow.tcurvelist[0].vertices)
        #
        # handle offset of fwd bottom nose curve:
        vertices[:,2] += self.offset
        # make:
        botcurve = rbspline(vertices, k=self.k, nump=self.nump)
        self.lower.insert(0, botcurve)
        #*******************************************************
        return
    
    
    def move_bulb_to_bare_hull(self):
        bbow = self.bbow
        fwd_hull_curve = self.tcurves[0]
        offset = self.offset
        #
        #*******************************************************
        # intercept bare hull fwd tcurve with bbow top aft.
        bbow_top_aft_pt = bbow.tcurvelist[0].vertices[-1]
        yb = bbow_top_aft_pt[1]
        #
        #*******************************************************
        #
        p_intercept = fwd_hull_curve.FindPoint(x=yb,
                                               index=1,
                                               maxiter=200,
                                               tol=1.e-10)
        self.p_intercept = p_intercept
        #
        #*******************************************************
        # real space location of intercept, denoted by 
        # the fwd tcurve
        barehull_real_intercept = fwd_hull_curve.CurvePoint(p_intercept[0])
        #Vector_hull_bulb = barehull_real_intercept - bbow_top_aft_pt
        Vector_hull_bulb = bbow_top_aft_pt - bbow_top_aft_pt
        #
        #handle offset:
        Vector_hull_bulb[2] += offset
        # translate:
        self.translate_bulb(vector = Vector_hull_bulb)
        self.barehull_real_intercept = barehull_real_intercept
        #
        return
    
    def translate_bulb(self, vector):
        """translate the bulbous bow by vector
        vector = Vector_hull_bulb
        """
        self.bbow = self.bbow.translate_surface(
                                            vector[0],vector[1],vector[2])
        self.bbow.compute_surface()
        return
    
    def plot(self):
        ax = self.mklwr.plot()
        ax = self.mkupr.plot(canvas=ax)
        return
    
    
    
    
    
    
    
    
              
class match_upper_to_lower(object):
    def __init__(self, tcurves, bulb):
        self.k = 4
        self.nump=30
        self.tcurves = tcurves
        self.bbow = bbow
        
    def refine_lcurves(self):
        """refine lcurves so as to 
        interpolate each vertex on the original bare hull
        transverse curves
        and each vertex on the bulb transverse curves
        """
        return
    
    def build_lagrangian(self):
        """
        """
        return
    
    
    
    

class ComplexHull(object):
    def __init__(self):
        self.k = 4
        self.nump=30
        self.barehull       = BareHull()
        self.bulb           = Bulb()

        self.lmap = {4:0,
                     5:1,
                     7:2,
                     11:3,
                     19:4,
                     35:5}

        self.hullcurves2    = get_tcurves()
        #
        #**********************************************************************
        # hi res:  solve the bulb longitudinals on the 35x35 level
        self.make_hires = False #  True #  
        self.highres_num = 19 #35
        self.lowres_num = 11 #19
        self.maxlevel_thbsolver = 4
        #
        #
        self.use_both_bare_hull_curves = True
        #
        # how about the starting curves for blending:
        self.makefresh_blends = True #start from scratch, 11 vertices
        # else: use the existing longitudinal
        #
        #*********************************************************************
        #
        #self.refine_cnt_method = 1 #count vertices
        self.refine_cnt_method = 2 #specify by height of bulb
        #
        #**********************************************************************
        # 
        self.num2use    = 5 #number of  transverse bulb curves 
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

    def make_bare_hull_3part(self):
        self.barehull.make()
        return 
    
    
    def make_bulb(self,ul=2,vl=1):
        self.bulb.make_surface(ul=ul,vl=vl)
        return
    
    
    def make_system(self, bulblevels = None):
        if bulblevels is None:
            bulblevels = [0,2]
            
        self.make_bare_hull_3part() #make the bare hull
        self.make_bulb(ul=bulblevels[0],
                       vl=bulblevels[1])
        
        tcurves = get_tcurves()
        sfh = subdivide_fwdhull(tcurves = tcurves[:4],
                                bbow = ch.bulb.hullbulb)
        self.subdivider = sfh
        self.subdivider.mklwr.make_lcurves()
        self.subdivider.mkupr.make_lcurves()
        #
        self.subdivider.mklwr.make_surface()
        self.subdivider.mkupr.make_surface()
        
        
        self.surfaces = {}
        #self.surfaces[0] = self.bulb.hullbulb
        self.surfaces[0] = self.subdivider.bbow
        self.surfaces[1] = self.subdivider.mklwr.surface
        self.surfaces[2] = self.subdivider.mkupr.surface
        self.surfaces[3] = self.barehull.hullmid
        self.surfaces[4] = self.barehull.hullaft
        return
    
    
    
    def export_rhino_surfaces(self):
        
        for key in self.surfaces:
            self.surfaces[key].IssueRhinoSurface(
                    the_filename = 'surface'+str(key)+'data.txt')
        return
    
    
    
    
    def add_bulb(self):
        k = self.barehull.hullfwd.tcurvelist[0].k
        nump = self.barehull.hullfwd.tcurvelist[0].nump
        #
        hullfwd = self.barehull.hullfwd.get_finest_level()
        #
        #**********************************************************************
        # 
        self.move_bulb_to_bare_hull()
        self.get_bulb_curves()
        #
        #**********************************************************************
        # decide how many levels to refine by
        #self.pfactor = self.estimate_bulbfractionalheight(self.p_intercept[0])
        #
        # refine once to get:
        #
        # 7x7 => 11x11 == 11x7 + 11x4
        self.pfactor = 1 
        #
        #**********************************************************************
        # save the base fwd bare hull
        #self.original_barehullfwd = copy.deepcopy(self.barehull.hullfwd) 
        #
        #**********************************************************************
        # dyadically regine the bare hull to prepare for bulb interpolation:
        # using pfactor from estimate_bulbfractionalheight (unless set not to)
        fs = self.setup_addbulb()
        self.fs = fs
        #
        #**********************************************************************
        #
        bbow = self.bulb.hullbulb
        #nbarehullu = hullfwd.tcurvelist[0].n
        nbbowu = self.bulb.hullbulb.tcurvelist[0].n
        #self.fwd_sketch_vertices = sketch(bbowcurve = bbow.tcurvelist[0], 
        #                                  hullfwdcurve = hullfwd.tcurvelist[0],
        #                                  hullaftcurve = hullfwd.tcurvelist[-1])
        
        self.fwd_sketch_vertices = sketch(bbowcurve = bbow.tcurvelist[-1], 
                                          hullfwdcurve = fs.tcurvelist[0],
                                          hullaftcurve = fs.tcurvelist[-1])
        #
        #**********************************************************************
        #
        self.makeblendedsurface(fs)
        #
        #**********************************************************************
        #
        return
    
    
    
    
    
    def makeblendedsurface(self, fs):
        """
        fwd enpoint of the new longitudinals will be:
            outline == self.fwd_sketch_vertices.vertices
            outline is no longer a single curve!
        
        assumption: we have finished refining the bare hull curves
        and are ready to interpolate the bulb.
        """
        hullfwd = fs
        bbow = self.bulb.hullbulb
        #
        #**********************************************************************
        # Bare hull 4 transverse curves - 2 boundary, 2 interior:
        self.barehullfwd_tcurvenet = copy.deepcopy(self.barehull.hullcurves[:4])
        #
        #**********************************************************************
        #
        #  * chief tcurves from the basic bare hull => interpolate these!
        #
        # refine them if you need too!
        #
        mssg = 'ERROR:  I though always, always.. that the incoming tcurve net was 7 vertices a piece.'
        FPcurvenet = [] 
        lvl = fs.tcurvelist[0].level
        for el in self.barehullfwd_tcurvenet[1:-1]:  #-1 so as to leave off the last!
            assert(el.n==7),mssg
            el.level = 2 #level =2 => its a 7 vertex curve.
            lel = el.level
            while lel<lvl:
                el = el.dyadic_refinement()
                lel = el.level
            FPcurvenet.append(el)
        self.FPcurvenet = FPcurvenet
        #
        #**********************************************************************
        #
        self.ntot = len(self.fwd_sketch_vertices.fwd_vertices)
        #
        #**********************************************************************
        # Create list of lists
        # vertices to interpolate for every ith curve to solve:
        interpolation_vertices = []
        bulbinterp_vertices = []
        nose_constraints = []
        longiconstraints = []
        for i, el in enumerate(self.fwd_sketch_vertices.fwd_vertices):
            interpolation_vertices.append([])
            nose_constraints.append([])
            bulbinterp_vertices.append([])
            #longiconstraints.append([])
        #
        #**********************************************************************
        # roll up all the vertices to be interpolated
        # start fwd and work aft
        #
        
        for i in range(self.ntot): #ith transverse column
            nose_constraints[i]=False
            #--------------------------------------------------------
            # bottom of bulb:
            if i == 0: 
                for curve in bbow.tcurvelist[::-1][1:]:
                    #interpolation_vertices[i].append(curve.vertices[i])
                    bulbinterp_vertices[i].append(curve.vertices[i])
                # and hull
                for curve in self.FPcurvenet:
                    interpolation_vertices[i].append(curve.vertices[i])
            #--------------------------------------------------------
            # interior of bulb:             bbow.tcurvelist[0].n-1 = 6 and curve.vertices[6] is the top.
            elif i>0 and i<bbow.tcurvelist[0].n-1:
                for curve in bbow.tcurvelist[::-1][1:]:
                    #interpolation_vertices[i].append(curve.vertices[i])
                    bulbinterp_vertices[i].append(curve.vertices[i])
                # and hull
                for curve in self.FPcurvenet:
                    interpolation_vertices[i].append(curve.vertices[i])
            #--------------------------------------------------------
            # top edge of bulb:
            elif i == bbow.tcurvelist[0].n-1: 
                for curve in bbow.tcurvelist[::-1][1:]:
                    #interpolation_vertices[i].append(curve.vertices[i])
                    bulbinterp_vertices[i].append(curve.vertices[i])
                # and hull
                for curve in self.FPcurvenet:
                    interpolation_vertices[i].append(curve.vertices[i])
            #--------------------------------------------------------
            # above the bulb
            else:
                #nose_constraints[i].append(self.fwd_sketch_vertices.fwd_vertices[i])
                nose_constraints[i] = True
                for curve in self.FPcurvenet:
                    interpolation_vertices[i].append(curve.vertices[i])
        #
        #**********************************************************************
        # 
        self.interpolation_vertices = interpolation_vertices
        self.bulbinterp_vertices = bulbinterp_vertices
        self.longiconstraints = longiconstraints #blank!
        self.nose_constraints = nose_constraints
        #
        #**********************************************************************
        #
        # SOLVE
        #
        #**********************************************************************
        # 
        self.solver_records = {}
        #
        self.complex_fwd_Lnet = []
        self.complex_fwd_lcurves = []
        for i in range(self.ntot): #ith transverse column
            #
            #--------------------------------------------------
            #    INIT NEW FWD CURVE HERE:
            #--------------------------------------------------
            #vini = self.barehull.hullfwd.tcurvelist[0].vertices[i]
            #vini = bhull.tcurvelist[0].vertices[i]
            vini = self.fwd_sketch_vertices.fwd_vertices[i]
            #--------------------------------------------------
            vfini = self.fwd_sketch_vertices.aft_vertices[i]
            #--------------------------------------------------
            linverts = linear_vertices(vini,vfini,4)
            curve = rbspline(linverts,k=4,nump=30)
            curve = curve.dyadic_refinement() #5
            curve = curve.dyadic_refinement() #7
            curve = curve.dyadic_refinement() #11 #users choice.  #allow up to 19 in solver.
            curve.parent = None
            self.solver_records[i] = records(curve)
            #--------------------------------------------------
            #
            print '------------------------------------------'
            print '\n\nBegin Solving New Longitudinal({}):\n\n'.format(i)
            print '------------------------------------------'
            """
            thbcurve    = curve
            
            bareverts   = self.interpolation_vertices[i]
            bulbverts = self.bulbinterp_vertices[i]
            nosetan     = self.nose_constraints[i]
            
            this_kind = 'equality'
            """
            #curve = self.solver_records[i].curve
            Lspline = self.setup_fwd_bulb_interpolation_solver(
                                        thbcurve    = curve,
                                        bareverts   = self.interpolation_vertices[i],
                                        bulbverts = self.bulbinterp_vertices[i],
                                        nosetan     = self.nose_constraints[i],
                                        this_kind = 'equality',
                                        bulbctype = 'LS')
            
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
        ##
        ##****************************************************************
        ## combine the surface
        
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
        
        ##
        ##****************************************************************
        ## plot
        self.complexhull.hullfwd.compute_surface()
        
        self.plotcomplex_hull(fancyish=True,
                              simplefancy=True)
        
        #alt change to fs.set_THB_curves_from_vertex_net(override_u=True, override_v = True)
        
        #fs = self.complexhull.hullfwd.get_finest_surface()
        fs.bounds[0][0].sup = 1. #.75 #.8125 #.75 #.6875 #.46875 #.4375 #.40625#.375#.34375#3125
        fs.bounds[0][1].sup = 1. #.4
        fs.compute_surface()
        self.complexhull.hullmid.numpu = 20
        self.complexhull.hullmid.numpv = 5
        self.complexhull.hullmid.compute_surface()
        self.complexhull.hullfwd.numpu = 30
        self.complexhull.hullfwd.numpv = 30
        self.complexhull.hullfwd.compute_surface()
        
        self.plotcomplex_hull(fancyish=True,
                              simplefancy=True)
        
        return
    
    
    """
    ax = Lspline.curve.plotcurvehierarchy()
    
    for curve in self.bulb.hullbulb.tcurvelist:
        ax = curve.plotcurvehierarchy(canvas = ax)
        
    ax = Lspline.curve.plotcurvehierarchy(canvas = ax)
    
    """
    
    def setup_fwd_bulb_interpolation_solver(self, 
                                           thbcurve, 
                                           bareverts, #interior vertices to interpolate
                                           bulbverts,
                                           nosetan = False,
                                           this_kind = 'equality',
                                           bulbctype = 'LS'): #bare hull vertices
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
        # hull interp
        nbvts = len(bareverts)
        ntot = thbcurve.n
        #*********************************************
        # bulb interp
        bulb_nbvts = len(bulbverts)
        #  thbcurve.n
        #*********************************************
        # bulb and hull to find interp ukbar
        totlen = nbvts+bulb_nbvts
        
        ukbar = np.linspace(0.,1.,totlen+2)[1:-1]
        bulb_ukbar = ukbar[:bulb_nbvts]
        hull_ukbar = ukbar[bulb_nbvts:]
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
        
        
        for pt,loc in zip(bulbverts,bulb_ukbar):
            #*********************************************
            FPD.add_xPointConstraint(kind= bulbctype, 
                                     value=pt[0], 
                                     location=loc, 
                                     weight = big ) 
            FPD.add_yPointConstraint(kind= bulbctype, 
                                     value=pt[1], 
                                     location=loc, 
                                     weight = big  )
            FPD.add_zPointConstraint(kind= bulbctype, 
                                     value=pt[2], 
                                     location=loc, 
                                     weight = big  )
        
        for pt,loc in zip(bareverts,hull_ukbar):
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
    
    
    
    def get_bulb_curves(self):
        bbow = self.bulb.hullbulb.get_finest_surface()
        #
        #**********************************************************************
        # get the curves
        self.add_tcurves = bbow.tcurvelist
        self.add_lcurves = bbow.lcurvelist
        return
    
    
    
    
    
    
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
        
        #        if self.refine_cnt_method == 1:
        #            while hullfwd.nu <2*nu or hullfwd.nv < 2*nv:
        #                print 'refining fwd surface for bulb'
        #                nlevel = hullfwd.dyadic_loft_refinement()
        #                hullfwd = nlevel
                
        #        elif self.refine_cnt_method == 2:
        #            for i in range(self.pfactor):
        #                hullfwd = hullfwd.dyadic_loft_refinement()
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
        
        # similar but superior idea:  2 surfaces.
        # 1 high, 1 low.
        # low surface incorporates bulb.
        # high surface matches low up to bulb corner. C2
        
        #
        # bulb nu == number of curves in the v direction
        # bulb nv == number of curves in the u direction
        active_v_bounds = hullfwd.mastercurve_u.active_bounds_given_cv_index(nu)
        active_u_bounds = hullfwd.mastercurve_v.active_bounds_given_cv_index(nv)
        
        active_v_bounds.inf = 0.
        active_u_bounds.inf = 0.
        box1 = thbspline.Box( active_u_bounds, active_v_bounds )
        hullfwd.bounds = thbspline.BoxList(box1)
        #
        #*****************************************************
        #
        return hullfwd
    
    
    
    def estimate_bulbfractionalheight(self, p_intercept):
        """
        NO LONGER USED
        -though there will be an issue of getting enough 
        resolution on the bottom edge of the upper surf
        in order to water tight transition to the lower surf.
        
        Ensure there are enough vertices to fair in the bulb
        
        * refine_to_fair_bulb_to_bare_hull
            -will use this directly to refine pfactor number of times
            where pfactor is the integer we return here.
            
        PROBLEM:  
            
            -so what if p is .75 of whatever,
            if there are enough dofs down there, 
            then you do not need to refine!
            
            -conversely, if there are to many, you may 
            need to push the upper lcurves `up'
            maybe the solver will do that now
            since we are solving all of them.
            
        Solution:
            2 surfaces 
            11x4 upper
            11x7 lower with bulb
            
        """
        hullfwd = self.barehull.hullfwd.get_finest_level()
        nhu = hullfwd.lcurvelist[0].n
        nu = len(self.bulb.tcurvenet)
        nv = len(self.bulb.lcurvenet)
        #
        #**********************************************************************
        #
        if p_intercept > .75:
            if nhu >= nu:
                return 0
            else:
                return 1
        elif p_intercept > .3:
            if nhu < 2*nu:
                return 1
            else:
                return 0
        else:
            if nhu < nu:
                return 3
            elif nhu < 2*nu:
                return 2
            elif nhu < 3*nu:
                return 1
    
    def move_bulb_to_bare_hull(self, offset='Not Used At This time!'):
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
        #
        barehull_real_intercept = fwd_hull_curve.CurvePoint(p_intercept[0])
        self.barehull_real_intercept = barehull_real_intercept
        #
        Vector_hull_bulb = barehull_real_intercept - bbow_top_aft_pt
        self.bulb.translate(vector = Vector_hull_bulb)
        return
    
    
    
    
    def check_solver(self, fs=None):
        if fs is None:
            fs = self.complexhull.hullfwd.get_finest_surface()
            
        if self.hullcurves2 is None:
            #self.hullcurves2 = get_tcurves()
            
            self.hullcurves2 = self.barehull.hullcurves
        
        # solved curves:
        #----------------------------------------------------------
        #these do not exist:
        #ax = self.aft_boundary_curve.plotcurvehierarchy()
        #ax = self.fwd_boundary_curve.plotcurvehierarchy(canvas=ax)
        #----------------------------------------------------------
        ax = self.fwd_sketch_vertices.bbowcurve.plotcurvehierarchy()
        for curve in self.hullcurves2:
            ax = curve.plotcurvehierarchy(canvas = ax)
        #this gave an odd result.  It looks like the old fwd bow curve:
        #ax = self.fwd_sketch_vertices.hullaftcurve.plotcurvehierarchy(canvas=ax)
        # but the vertices show the correct curve.
        #
        #this one is not true:
        # b/c the fwd vertices are "broken"
        ##ax =!= self.fwd_sketch_vertices.hullfwdcurve.plotcurvehierarchy(canvas=ax)
        #----------------------------------------------------------
        for i in range(self.ntot):
            #fs.lcurvelist[i].compute_curve()
            curve = rbspline(fs.lcurvelist[i].vertices,k=4,nump=30)
            #ax = fs.lcurvelist[i].plotcurvehierarchy(canvas = ax)
            ax = curve.plotcurvehierarchy(canvas = ax)
            
        #for tcurve in self.bulb.tcurvenet:
        #    ax = tcurve.plotcurvehierarchy(canvas=ax)
        
        
        #----------------------------------------------------------
        #these do not exist:
        #ax = self.aft_boundary_curve.plotcurvehierarchy()
        #ax = self.fwd_boundary_curve.plotcurvehierarchy(canvas=ax)
        #----------------------------------------------------------
        ax = self.fwd_sketch_vertices.bbowcurve.plotcurvehierarchy()
        #ax = self.hullcurves2[0].plotcurvehierarchy()
        for lcurve in fs.lcurvelist:
            #lcurve.compute_curve()
            curve = rbspline(lcurve.vertices,k=4,nump=30)
            ax = curve.plotcurvehierarchy(canvas=ax)
        for tcurve in self.hullcurves2:
            ax = tcurve.plotcurvehierarchy(canvas=ax)
        
        
        for tcurve in self.bulb.tcurvenet:
            ax = tcurve.plotcurvehierarchy(canvas=ax)
        
        return
    
    def add_subdivided_fwd_hull(self, lower, upper, bulb):
        self.fwdlower = lower
        self.fwdupper = uppder
        self.fwdbulb = bulb
        return
    
    def plot_initial_situation(self,
                                 fancyish = True,
                                 color=False,
                                 simplefancy=True):
        """
        
        dev
        ----------
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
    
    
    def make_lines(self, canvas = None):
        if canvas is None:
            ax = None
        else:
            ax = canvas
        
        
        #self.hullfwd.fwdlower.make_lines_curves()
        #self.hullfwd.fwdupper.make_lines_curves()
        self.bulb.hullbulb.make_lines_curves()
        
        #self.mklwr.hullfwd.make_lines_curves()
        #self.mkupr.hullfwd.make_lines_curves()
        self.subdivider.mklwr.surface.make_lines_curves()
        self.subdivider.mkupr.surface.make_lines_curves()
        
        self.barehull.hullmid.make_lines_curves()
        self.barehull.hullaft.make_lines_curves()
        
        
        ax = self.bulb.hullbulb.plot_lines(canvas = ax)
        ax = self.subdivider.mklwr.surface.plot_lines(canvas = ax)
        ax = self.subdivider.mkupr.surface.plot_lines(canvas = ax)
        
        ax = self.barehull.hullmid.plot_lines(canvas = ax)
        ax = self.barehull.hullaft.plot_lines(canvas = ax)
        
        return
    
    def plotcomplex_hull(self,
                         fancyish = False,
                         color=False,
                         simplefancy=False,
                         bulbtcurves=True):
        """
        fancyish=True
        color = True
        bulbtcurves=True
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
        self.canvas = ax
        return ax
    
    
    
    def export_rhino_textfiles(self):
        
        tlist = []
        for curve in self.tcurvenet:
            tlist.append([ list(pt)  for pt in curve.vertices ] )
        tknots = []
        for curve in self.tcurvenet:
            tknots.append(list(curve.t[1:-1]))
        
        llist = []
        for curve in self.lcurvenet:
            llist.append([ list(pt)  for pt in curve.vertices ] )
        lknots = []
        for curve in self.lcurvenet:
            lknots.append(list(curve.t[1:-1]))
            
        with open('transverse_curves','w') as f:
            pickle.dump(tlist,f)
        with open('longitudinal_curves','w') as f:
            pickle.dump(llist,f)
        
        with open('transverse_knots','w') as f:
            pickle.dump(tlist,f)
        with open('longitudinal_knots','w') as f:
            pickle.dump(llist,f)
        
        the_filename = 'transverse_curves.txt'
        with open(the_filename, 'w') as f:
            for el in tlist:
                for s in el:
                    st = '{}, {}, {}\n'.format(s[0],s[1],s[2])
                    f.write(st)
                f.write('new curve   \n')
        
        the_filename = 'transverse_knots.txt'
        with open(the_filename, 'w') as f:
            for el in tknots:
                st = ''
                for s in el:
                    st = st+' '+str(s)
                f.write(st)
                
                f.write('\n new knot vector \n')
        
        the_filename = 'longitudinal_curves.txt'
        with open(the_filename, 'w') as f:
            for el in llist:
                for s in el:
                    st = '{}, {}, {}\n'.format(s[0],s[1],s[2])
                    f.write(st)
                f.write('new curve   \n')
        
        the_filename = 'longitudinal_knots.txt'
        with open(the_filename, 'w') as f:
            for el in lknots:
                st = ''
                for s in el:
                    st = st+' '+str(s)
                f.write(st)
                
                f.write('\n new knot vector \n')
        return


if __name__ == '''__main__''':
    ch = ComplexHull()
    #self = ch
    ch.make_system()
    
    #"""
    #-------------------------------------
    #moved within make_system call
    tcurves = get_tcurves()
    #    sfh = subdivide_fwdhull(tcurves = tcurves[:4],
    #                            bbow = ch.bulb.hullbulb)
    #    self = sfh
    #    self.mklwr.make_lcurves()
    #    self.mkupr.make_lcurves()
    #    self.mklwr.make_surface()
    #    self.mkupr.make_surface()
    #-------------------------------------
    #"""
    self = ch.subdivider
    ## forward:
    ax = self.mklwr.surface.plotSurface(fancy=False)
    ax = self.mkupr.surface.plotSurface(fancy=False,
                                        canvas = ax)
    ## bulb:
    ax = self.bbow.plotSurface(fancy=False,
                               canvas = ax)
    ## mid and aft
    ax = ch.barehull.hullmid.plotSurface(fancy=False,
                                           canvas = ax)
    ax = ch.barehull.hullaft.plotSurface(fancy=False,
                                           canvas = ax)
    
    ## forward
    ax = self.mklwr.surface.plotSurface(fancy=True,
                                        simplefancy=True)
    ax = self.mkupr.surface.plotSurface(fancy=True,
                                        simplefancy=True,
                                        canvas = ax)
    ## bulb:
    ax = self.bbow.plotSurface(fancy=True,
                               simplefancy=True,
                               canvas = ax)
    ## mid and aft
    ax = ch.barehull.hullmid.plotSurface(fancy=True,
                                           simplefancy=True,
                                           canvas = ax)
    ax = ch.barehull.hullaft.plotSurface(fancy=True,
                                           simplefancy=True,
                                           canvas = ax)
    ## other:
    #"""
    self.check_boundary()
    
    self.canvas = ax
    
    #ch.make_lines() #no good
    
    ch.export_rhino_surfaces()