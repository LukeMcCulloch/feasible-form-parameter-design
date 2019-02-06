# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 18:55:13 2016

@author: lukemcculloch

notes:
    Form parameter Design of a Bare Hull
    using Hull Inference Constraint Logic Programming

    correct expression for the correctly masked SAC gradient:
        sac = self.Lsplines.SAC
        mask = sac.mask
        sac.f.grad[:].T[mask]
        
File:
    purpose:    This file builds a bare hull
    contains:   methods to build a hull using form parameter design
                    if in possession of a feasible set of form parameters

TODO: 
    The code in this file is horrific.  Refactor with monads:
        https://codon.com/refactoring-ruby-with-monads
        
    
    todo: 
        replace 'bowfairning' with a name
        that is spelled correctly :|
"""
#no support in anaconda python for rhinoscriptsyntax!
#import rhinoscriptsyntax as rs
#
#import pickle #do it the old fashioned way...
import cPickle as pickle #C pickle with no subclassing
#
import numpy as np
import copy
import matplotlib.pyplot as plt
#
from interval_arithmetic import ia
#import sqKanren as lp
from hull_inference_ob_graph import Hull as hullclp #not used in opt_simple
# instead, opt uses simple_hull_rules!
#
#from extended_interval_arithmetic import ia
#import sqKanren as lp
#from design_space import HullSpace as hullclp
#
#from OldHullDesigner import Hullcurve#, linear_vertices, transpose_vertices
import curve             as     spline
from   ADILS             import IntervalLagrangeSpline, \
                                Lagrangian, \
                                LazyFPDConstraintInterface
from   adials_gui3d      import frame, DrawCurveInteractive, hull_frame
from   FormParameter     import FormParameterDict
#from   initialValues     import InitializeControlPoints 
from   initialValues     import InitializeControlVertices
from   initialValues     import interval_bounds, lagrangian_bounds
from   plots             import Plotter

from frames import Frame,DeepVector

import utility_optimization as uopt
from Quaternion import Quaternion
import SpecializedCurveSolvers as scsolvers

#
import mx_plots as spfuncs # circle, gfunc class, etc...
# 

package_norms = uopt.package_norms
#
from  automatic_differentiation import ad as ad #doing this purely to assist in latex table generation
#

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
    

class Hullcurve(object):
    """Some really old code
    -offset the midship section for a nice looking flat of side...
    -worry about 7-11 dyadic refinement later... (DONE)
    -possibly delay worry until O(N) 2nd Gen Wavelet Splines are in effect
    """
    def __init__(self, kind, initial_curve, plane_dim=None, offset_const = None, definition_curve=None):
        self.k                  = initial_curve.k
        self.nump               = initial_curve.nump
        self.plane_dim          = plane_dim
        self.initial_curve      = initial_curve
        self.offset_const       = offset_const
        self.definition_curve   = definition_curve
        self.kind               = kind #plane or offset
        if self.kind == 'plane':
            self.plane()
        elif self.kind =='offset':
            self.offset()
        return
        
    def plane(self):
        return_vertices = self.initial_curve.vertices
        return_vertices[:,self.plane_dim] = 0.0
        self.return_curve = spline.Bspline(return_vertices, self.k, self.nump)
        return 
    
    def offset(self):
        return_vertices = copy.copy(self.initial_curve.vertices)
        if self.offset_const != None:
            print 'offset' 
            return_vertices[:,self.plane_dim] += self.offset_const
        if self.definition_curve != None:
            return_vertices[:,:] += self.definition_curve.vertices[:,:] - self.initial_curve.vertices[:,:]
        self.return_curve = spline.Bspline(return_vertices, self.k, self.nump)
        return            
        

class Stations(object):
    """Stations are defined with respect to hull.LengthDWL"""
    def __init__(self, lwl, cp=.5):
        self.lwl            = lwl
        self.cp             = cp
        return

    def __str__(self):
        keys = self.__dict__.keys()
        del(keys[keys.index('cp')])
        del(keys[keys.index('lwl')])
        for key in keys:
            print key, self.__dict__[key]
        return ''

    @property
    def FOSAC(self):
        return self._FOSAC
    @FOSAC.setter
    def FOSAC(self, args):
        len_flat, cp = args
        pflat  = len_flat/self.lwl
        a = (cp-.5*pflat)
        b = (cp+.5*pflat)
        self._FOSAC = [a*self.lwl,b*self.lwl,a,b]
        return

    @property
    def FOWL(self):
        """flat of design water line"""
        return self._FOWL
    @FOWL.setter
    def FOWL(self, args):
        len_flat, cp = args
        pflat  = len_flat/self.lwl
        a = (cp-.5*pflat)
        b = (cp+.5*pflat)
        self._FOWL = [a*self.lwl,b*self.lwl,a,b]
        return

    @property
    def FOS(self):
        return self._FOS
    @FOS.setter
    def FOS(self, args):
        len_flat, cp = args
        pflat  = len_flat/self.lwl
        a = (cp-.5*pflat)
        b = (cp+.5*pflat)
        self._FOS = [a*self.lwl,b*self.lwl,a,b]
        return

    @property
    def FOB(self):
        return self._FOB
    @FOB.setter
    def FOB(self, args):
        len_flat, cp = args
        pflat  = len_flat/self.lwl
        a = (cp-.5*pflat)
        b = (cp+.5*pflat)
        self._FOB = [a*self.lwl,b*self.lwl,a,b]
        return

    @property
    def FOCP(self):
        """Flat of Center Plane"""
        return self._FOCP
    @FOCP.setter
    def FOCP(self, args):
        len_flat, cp = args
        pflat  = len_flat/self.lwl
        a = (cp-.5*pflat)
        b = (cp+.5*pflat)
        self._FOCP = [a*self.lwl,b*self.lwl,a,b]
        return

    @property
    def Stern_Profile(self):
        return self._Stern_Profile
    @Stern_Profile.setter
    def Stern_Profile(self, args):
        start,end = args
        a = start/self.lwl
        b = end/self.lwl
        self._Stern_Profile = [start,end,a,b]
        return


    @property
    def Stern_Transition(self):
        return self._Stern_Transition
    @Stern_Transition.setter
    def Stern_Transition(self, args):
        start,end = args
        a = start/self.lwl
        b = end/self.lwl
        self._Stern_Transition = [start,end,a,b]
        return
    
    
    @property
    def Stern_Fairing(self):
        return self._Stern_Fairing
    @Stern_Fairing.setter
    def Stern_Fairing(self, args):
        start,end = args
        a = start/self.lwl
        b = end/self.lwl
        self._Stern_Fairing = [start,end,a,b]
        return

    @property
    def Bow_Profile(self):
        return self._Bow_Profile
    @Bow_Profile.setter
    def Bow_Profile(self, args):
        start,end = args
        if not (start ==0.):print "Bow Profile is starting from {}".format(start)
        a = start/self.lwl
        b = end/self.lwl
        self._Bow_Profile = [start,end,a,b]
        return

    @property
    def Bow_Fairing(self):
        return self._Bow_Fairing
    @Bow_Fairing.setter
    def Bow_Fairing(self, args):
        start,end = args
        #if not (start ==0.):print "Bow Fairing is starting from {}".format(start)
        a = start/self.lwl
        b = end/self.lwl
        self._Bow_Fairing = [start,end,a,b]
        return

    @property
    def Bow_Bulb_Transition(self):
        return self._Bow_Bulb_Transition
    @Bow_Bulb_Transition.setter
    def Bow_Bulb_Transition(self, args):
        start,end = args
        a = start/self.lwl
        b = end/self.lwl
        self._Bow_Bulb_Transition = [start,end,a,b]
        return




##
##SAC\FD========================================================
##
##
##SAC\FD========================================================
##
class SAC(object):
    k       = 4
    nump    = 30
    xb      = 0.
    yb      = 0.
    def __init__(self, Amax, v1,v2,v3,
                 l1,l2,l3,
                 Xe, Xm, Xr, Xlcb=None,
                 Abfc=None, Abfc_loc=None, Asfc=None, Asfc_loc=None,
                 Afp=None, Aap=None, Cpe=None, Cpr=None,
                 stations=None):
        #self.V         = vol
        self.Amax       = Amax
        self.Ve         = v1
        self.Vm         = v2
        self.Vr         = v3
        self.Le         = l1 #length of the entry
        self.Lm         = l2 #length of the midsection
        self.Lr         = l3 #length of the run
        self.L          = l1+l2+l3
        self.Xe         = Xe
        self.Xm         = Xm
        self.Xr         = Xr
        self.Xlcb       = Xlcb
        self.Afp        = Afp
        self.Aap        = Aap
        self.Cpe        = Cpe
        self.Cpr        = Cpr
        self.Abfc       = Abfc
        self.Asfc       = Asfc
        self.Abfc_loc   = Abfc_loc
        self.Asfc_loc   = Asfc_loc
        self.fairbow    = False
        self.fairstern  = False
        self.stations   = stations
        self.go()
        #self.aggregate_SAC_curves_Small()
        return

    def plot(self):
        """Plot the composite SAC curve
        with all Form Parameters used in
        its construction
        """
        self.nice_plot = Plotter(x_axis_name = 'x',
                                 y_axis_name = 'z')
        ADLsplines = [self.entrance,
                        self.midbody,
                        self.run]
        self.nice_plot(self.SAC,
                       ADLsplines,
                       mytitle = 'SAC')

        return

    def go(self):
        self.entrance   = self.compute_entrance_SAC()
        self.midbody    = self.compute_midbody_SAC()
        self.run        = self.compute_run_SAC()
        return

    def aggregate_SAC_curves_regular(self):
        """put fwd, mid, and run
        into one SAC curve
        """
        n3 = self.run.curve.n
        n2 = self.midbody.curve.n
        n1 = self.entrance.curve.n
        self.midbody.translate(self.Le+SAC.xb)
        #self.entrance.translate(self.Lm + self.Le+SAC.xb)
        nv = n1 + n2 + n3

        c = sac.run.curve.curve_join

        Cv = np.zeros((nv,2),float)
        Cv[0:n1,:] = self.run.curve.vertices
        Cv[n1:n1+n2,:] = self.midbody.curve.vertices
        Cv[n1+n2:n1+n2+n3,:] = self.entrance.curve.vertices


        self.SAC = spline.Bspline(Cv, SAC.k, nump = SAC.nump*3)
        return

    def aggregate_SAC_curves_Small(self):
        """put fwd, mid, and run
        into one SAC curve
        """
        n1 = self.run.curve.n
        n2 = self.midbody.curve.n
        n3 = self.entrance.curve.n
        #self.midbody.translate(self.Le+SAC.xb)
        #self.run.translate(self.Lm + self.Le+SAC.xb)
        nv = n1 + n2-2 + n3

        Cv = np.zeros((nv,2),float)
        Cv[0:n1,:] = self.entrance.curve.vertices
        Cv[n1:n1+n2-2,:] = self.midbody.curve.vertices[1:-1]
        Cv[n1+n2-2:n1+n2+n3,:] = self.run.curve.vertices

        self.SAC = spline.Bspline(Cv, SAC.k, nump = SAC.nump*3)
        return

    def compute_entrance_SAC(self):
        """
        from hull_from_simple_designspace import SAC
        from   initialValues     import InitializeControlVertices
        import curve             as     spline
        
        import curve             as     spline
        from   ADILS             import IntervalLagrangeSpline, \
                                        Lagrangian, \
                                        LazyFPDConstraintInterface
        from   adials_gui3d      import frame, DrawCurveInteractive, hull_frame
        from   FormParameter     import FormParameterDict
        #from   initialValues     import InitializeControlPoints 
        from   initialValues     import InitializeControlVertices
        from   initialValues     import interval_bounds, lagrangian_bounds
        from   plots             import Plotter
        """
        k       = SAC.k
        nump    = SAC.nump
        nCV     = 7
        curve_area = self.Ve
        xb      = SAC.xb #0.  #  xb = self.xb
        xe      = self.Le #56
        yb      = SAC.yb      #    yb = self.yb
        ye      = self.Amax
        #yfp     = self.Afp
        Xc      = self.Xe
        ab = 0.
        ae = 0.
        ini_v = InitializeControlVertices(xb,yb,xe,ye,
                                          alphae=ae,
                                          alphab=90.,
                                          area = curve_area,
                                          Cab_given=90.,
                                          Cae_given=0.,
                                          nCV = 7,
                                          slope = 'up')
        curve = spline.Bspline(ini_v.vertices, k, nump)
        """
        curve.compute_area()
        curve.compute_moments()
        curve.computeXc()
        """

        interval_data, small = interval_bounds(curve)

        FPD = FormParameterDict(curve)
        FPD.add_AreaConstraint(kind='equality', 
                               value = curve_area)
        #FPD.add_AngleConstraint(kind = 'LS', value = 0., location = 0.)
        FPD.add_AngleConstraint(kind = 'equality', 
                                value = 0., 
                                location = 1., 
                                weight=100.)
        #FPD.add_CurvatureConstraint(kind = 'equality', value = 0., location = 0.)
        #FPD.add_CurvatureConstraint(kind = 'LS', value = 0., location = 1.)
        if False:
            FPD.add_XcConstraint(kind = 'equality', value = Xc)#, weight=100.)
        #
        #-----------------------
        #
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = .5)
        #FPD.add_E3(kind='LS', weight = .1)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        #        Lspline.display(mytitle = 'SAC_run_ini', x_axis_name = 'x',y_axis_name = 'z')
        Lspline.optimize(stop = 50)

        ##
        ##  curve optimization run 2
        ##----------------------- Increased Curve Contraints
        ##
        #FPD = FormParameterDict(Lspline.curve)
        FPD = FormParameterDict(copy.deepcopy(Lspline.curve))
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve)
        FPD.add_AreaConstraint(kind='equality', 
                               value = curve_area)
        FPD.add_AngleConstraint(kind = 'equality', 
                                value = 0., 
                                location = 1., 
                                weight=100.)
        if True:
            FPD.add_XcConstraint(kind = 'equality', 
                                 value = Xc)#, weight=100.)
        #
        #----------------------- SAC value at bow fairness curve location:
        #
        #
        # Back to Basics!
        #   hull curves define station curves
        #   NOT THE OTHER WAY AROUND!!
        #           Feb 5 2017
        #
        if self.fairbow:
            start,end,ao,bo = self.stations.Bow_Fairing
            AreaBowFairingCurve = self.Abfc
            where = ao
            FPD.add_yPointConstraint(kind='equality',
                                     location = where,
                                     value = AreaBowFairingCurve)
            FPD.add_xPointConstraint(kind='equality',
                                     location = where,
                                     value = start)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .1)
        FPD.add_ArcLengthApprox(kind='LS', weight = 10.5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        return Lspline


    def compute_midbody_SAC(self):
        k       = SAC.k
        nump    = SAC.nump
        n       = 7
        curve_area = self.Vm
        xb = self.Le
        xe = xb + self.Lm
        yb = self.Amax
        ye = self.Amax
        Xc = self.Xm
        ab = 0.
        ae = 0.
        vertices = linear_vertices([xb,yb],
                                   [xe,ye],
                                   num=n)
        curve = spline.Bspline(vertices, k, nump)

        interval_data, small = interval_bounds(curve)

        FPD = FormParameterDict(curve)
        #FPD.add_AreaConstraint(kind='equality', value = curve_area)
        FPD.add_AngleConstraint(kind = 'equality', value = 0., location = 0.)
        FPD.add_AngleConstraint(kind = 'equality', value = 0., location = 1.)
        FPD.add_CurvatureConstraint(kind = 'equality', value = 0., location = 0.)
        FPD.add_CurvatureConstraint(kind = 'equality', value = 0., location = 1.)
        #FPD.add_XcConstraint(kind = 'equality', value = Xc)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)

        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        #Lspline.compute_lagrangian()
        #        Lspline.display(mytitle = 'SAC_mid_ini',
        #                        x_axis_name = 'x',
        #                        y_axis_name = 'z')
        Lspline.optimize(stop = 50)
        #        Lspline.display(mytitle = 'SAC_mid',
        #                        x_axis_name = 'x',
        #                        y_axis_name = 'z')
        return  Lspline


    def compute_run_SAC(self):
        """
            TODO: check if you need to compute
            offset AREA!
        """
    #def compute_entrance_SAC(self):
        """
        """
        k       = SAC.k
        nump    = SAC.nump
        nCV     = 7
        curve_area = self.Vr
        xb  = self.Le + self.Lm#+ SAC.xb
        xe  = xb + self.Lr#+
        yb  = self.Amax
        ye  = 0.
        #yap = self.Aap
        Xc  = self.Xr#-(self.Le + self.Lm)
        ab = 0.
        ae = 0.
        ini_v = InitializeControlVertices(xb,yb,xe,ye,
                                          alphae=ae,
                                          alphab=ab,
                                          Cab_given=0.,
                                          Cae_given=0.,
                                          nCV = 7,
                                          slope = 'down')
        curve = spline.Bspline(ini_v.vertices, k, nump)

        interval_data, small = interval_bounds(curve)

        FPD = FormParameterDict(curve)
        FPD.add_AreaConstraint(kind='equality', 
                               value = curve_area)
        #FPD.add_AngleConstraint(kind = 'LS', value = 0., location = 0., weight=100.)
        #FPD.add_AngleConstraint(kind = 'LS', value = 0., location = 1.)
        #FPD.add_CurvatureConstraint(kind = 'LS',value = 0.,location = 0.)
        #FPD.add_CurvatureConstraint(kind = 'equality',value = 0.,location = 1.)
        FPD.add_AngleConstraint(kind = 'equality', 
                                value = 0., 
                                location = 0., 
                                weight=100.)
        if False:
            FPD.add_XcConstraint(kind = 'equality', value = Xc)
        #
        #-----------------------
        #
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = .5)
        #FPD.add_E3(kind='LS', weight = .1)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)

        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        #Lspline.compute_lagrangian()
        #        Lspline.display(mytitle = 'SAC_fwd_ini', x_axis_name = 'x',y_axis_name = 'z')
        Lspline.optimize(stop = 50)
        #        Lspline.display(mytitle = 'SAC_fwd',
        #                        x_axis_name = 'x',
        #                        y_axis_name = 'z')
        ## Lspline.curve.plot()
        ##
        
        
        ##----------------------- Increased Curve Contraints
        ##
        #FPD = FormParameterDict(Lspline.curve)
        FPD = FormParameterDict(copy.deepcopy(Lspline.curve))
        FPD.add_AreaConstraint(kind='equality', 
                               value = curve_area)
        #FPD.add_AngleConstraint(kind = 'LS', value = 0., location = 0., weight=100.)
        #FPD.add_AngleConstraint(kind = 'LS', value = 0., location = 1.)
        #FPD.add_CurvatureConstraint(kind = 'LS',value = 0.,location = 0.)
        FPD.add_AngleConstraint(kind = 'equality', 
                                value = 0., 
                                location = 0., 
                                weight=100.)
        if True:
            FPD.add_XcConstraint(kind = 'equality', 
                                 value = Xc)
        #
        #----------------------- SAC value at stern fairness curve:
        #
        #
        # Back to Basics!
        #   hull curves define station curves
        #   NOT THE OTHER WAY AROUND!!
        #           Feb 5 2017
        #
        if self.fairstern:
            start,end,ao,bo = self.stations.Stern_Fairing
            AreaSternFairingCurve = self.Asfc#/2.
            where = ao
            FPD.add_yPointConstraint(kind='equality',
                                     location = where,
                                     value = AreaSternFairingCurve)
            FPD.add_xPointConstraint(kind='equality',
                                     location = where,
                                     value = start)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .1)
        FPD.add_ArcLengthApprox(kind='LS', 
                                weight = 10.5)
        L = Lagrangian(FPD)

        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        
        return  Lspline
##
##SAC\FD========================================================
##
##
##SAC\FD========================================================
##
class Lsplines(object):
    def __int__(self):
        self.SAC = None
        self.DWL = None
        self.CProfile = None
        self.FOS = None
        self.midSS = None
        self.stern = None
        self.sternfairing = None
        self.bowfairning = None
        self.bowfairing_fwd = None
        self.bowfairing_aft = None
        self.bow = None
        self.bbow = None
        return


class Hull(object):
    """
    Pass States attribute a Kanren style environement
    """
    def __init__(self, 
                 hullconstraints=None, 
                 States=None, 
                 which_state = 0,
                 revise=True,
                 verbose=False,
                 hullformtype='osv',
                 hyperparameter_input=None,
                 bulbconstraints = None,
                 xax = (1.,0.,0.), 
                 yax = (0.,1.,0.), 
                 zax = (0.,0.,1.) ):
        """
          
        Bare Hull Coordinate System:
        
                 y    fwd
                 |
                 | 
                 | 
        port     o ------x   starboard
                /  
               /
              z  
            
            aft
        
        
        x : transverse      (+ starboard)
        y : vertical        (+ up)
        z : longitudinal    (+ aft)
         
         
        TODO: refactor the getters to @property?
         
        use ShipDesigner.complex_hull_design_and_generation
        to pass bulbconstraints in automatically
        
        """
        self.k                      = 4
        self.n                      = 7
        self.nump                   = 30
        self.small                  = 1.e-3
        self.tol                    = 1.e-6
        self.revise                 = revise
        self._verbose               = verbose
        self.hullformtype           = hullformtype
        self.hyperparameter_input   = hyperparameter_input
        self.package_norms          = package_norms
        self.bulbconstraints        = bulbconstraints ##SD.bbow.rgp
        #
        self.results_dict           = {}
        #
        #**********************************************************************
        #**********************************************************************
        #  Hyperparameter (hand) Tuning:
        #
        self.remove_fwdfairness = False
        self.heavy_config = False
        self.remove_aft_fos = False 
        #
        #*********************
        # transition to auto-tune
        self.use_in_linear_interp = []
        self.use_in_longi_solver = []
        self.tcurve_map = {}
        
        #documentaion:
        self.midshipindices = (3,4,5,6,7,8,9,10)
        self.fwdFOSindices = (3,4,5)
        self.FOSACindices = (6,7)
        self.aftFOSindices = (8,9,10)
        self.aftindices = (11,12,13)
        
        
        #indices to make LS in longitudinal solver
        self.use_in_linear_interp = [0,1,2,
                                     4,5,
                                     6,7,
                                     8,9,
                                     11,13]
        #indices to make LS in longitudinal solver
        self.use_in_longi_solver = [0,1,2,
                                    3,4,5,
                                    6,7,
                                    8,9,10,
                                    11,12,13]
        #indices to make LS in longitudinal solver
        self.dropout = [] 
        
        #TEST:
        self.dropout = [3] #indices to make LS in longitudinal solver
        self.ignore = [7,10,12] #indicies to totally ignore in the longitudinal solver
        
        if bulbconstraints is not None:
            self.bbow_from_SAC = False
            self.BulbVolume = bulbconstraints.get_name_and_value(
                                        'BulbVolume')[1][0].getpoint(.5)
            self.A_mid = bulbconstraints.get_name_and_value(
                                        'A_mid')[1][0].getpoint(.5)
            self.A_lateral = bulbconstraints.get_name_and_value(
                                        'A_lateral')[1][0].getpoint(.5)
            self.BulbLength = bulbconstraints.get_name_and_value(
                                        'BulbLength')[1][0].getpoint(.5)
            self.BulbDepth = bulbconstraints.get_name_and_value(
                                        'BulbDepth')[1][0].getpoint(.5)
            self.BulbBeam = bulbconstraints.get_name_and_value(
                                        'BulbBeam')[1][0].getpoint(.5)
        else:
            self.bbow_from_SAC = True
        
        """
        tcurve_map:
            
        0: nose
        1: fwd fairness
        2: fwd transition
        3: fwd FOS outer
        4: fwd FOS middle
        5: fwd FOS inner
        6: FOSAC start
        7: FOSAC end
        8: aft FOS inner
        9: aft FOS middle
        10: aft FOS outer
        11: aft transition
        12: aft fairness
        13: transom
        """
        if self.hyperparameter_input is None:
            
            """
            """
            self.fwd_fairness_location_ratio    = .3   #smaller is further fwd
            #self.fwd_fairness_location_ratio    = .28   #smaller is further fwd
            #self.fwd_fairness_location_ratio    = .35   #smaller is further fwd
            #
            self.fwd_transition_location_ratio  = .5   #larger is closer to midship
            self.fwd_transition_location_ratio  = .6   #larger is closer to midship
            #
            self.aft_transition_location_ratio  = .35    #smaller is closer to midship
            self.aft_transition_location_ratio  = .4    #smaller is closer to midship
            #self.aft_fairness_location_ratio    = .7   #larger is closer to aft edge
            self.aft_fairness_location_ratio    = .65   #larger is closer to aft edge
            #
            #**********************************************************************
            #**********************************************************************
            #
            self.aft_drop = .25
            #**********************************************************************
            #**********************************************************************
            #
            # midship vector was responsible for
            # curves (by index) [3,4,5,6,7,8,9,10]
            #
            # transverse curve indices to use for linear interpolation
            self.use_in_linear_interp = [0,1,2,
                                         4,5,
                                         6,7,
                                         8,9,
                                         11,13]
            # 
            self.use_in_longi_solver = [0,1,2,
                                        3,4,5,
                                        6,7,
                                        8,9,10,
                                        11,12,13]
            self.dropout = [4] #indices to make LS in longitudinal solver
            self.ignore = [3,7,10] #indicies to totally ignore in the longitudinal solver
            
        ##bottom##
        # End of Hyperparameter Search 'By Hand'  
        #**********************************************************************
        #
        else:
            #input from outside the class
            self.change_hyperparameters(self.hyperparameter_input)
            """
            print 'Interface for automated hyperparameter inputs... some day soon'
            print 'Not at all ready!  Not today!'
            hyptuples, hclass = self.hyperparameter_input
            #build rules (not a 'real' rules net) with hyptuples
            self.fwd_fairness_location_ratio    = hclass.fwd_fairness.getpoint(.38)
            self.fwd_transition_location_ratio  = hclass.fwd_transition.getpoint(.78)
            self.aft_transition_location_ratio  = hclass.aft_transition(.45 )
            self.aft_fairness_location_ratio    = hclass.aft_fairness.getpoint(.5) 
            
            self.midship_design_vector = {1:'h',   # h=3, s=2, l=1
                                          2:'h',   # h=2, s=1, l=0
                                          3:'l',}  # h=3, s=2, l=1
            #"""
        
        
        #**********************************************************************
        #**********************************************************************
        #
        """
        self.use_in_linear_interp = []
        self.use_in_longi_solver = []
        
        self.midshipindices = (3,4,5,6,7,8,9,10)
        self.fwdFOSindices = (3,4,5)
        self.FOSACindices = (6,7)
        self.aftFOSindices = (8,9,10)
        self.aftindices = (11,12,13)
        #"""
        
        #
        # Done with hyperparameters
        #**********************************************************************
        print '--------------------------------------------------------------'
        print 'hyperparameter configuration:\n'
        #print 'midship info:'
        #print self.midship_design_vector
        print ''
        print 'longitudinal dropout (which transverse vertices '
        print 'to  least squares approximately interpolate):'
        print self.dropout
        print 'transverse curves to completely ignore'
        print self.ignore
        print 'The other vertices will be exactly interpolated'
        print ''
        print 'There are 14 separate possible '
        print ' transverse curves for use here.'
        print ''
        print 'The longitudinals are always composed of 11 control vertices'
        print ''
        print '--------------------------------------------------------------'
        #**********************************************************************
        #
        self.fairnesscurves_drive_hull_curves = False #put an end to this idea
        if hullconstraints is None and States is not None:
            self.hullconstraints = States#make states.state from state which inherits from Kanren State!
            self.hullconstraints.state = States.states[which_state]
        else:
            # Jan 10 2018:
            # can be like hullrulesnet.rgp.env using the most up to date rules language
            self.hullconstraints    = hullconstraints #inherits from Kanren States!
        #self.hullconstraints.state = hullconstraints.states[which_state]
        self.stations           = {}
        self.Lsplines = Lsplines()
        self.set_list_of_hull_methods()
        self.example_transverse= spline.Bspline(
                linear_vertices([0.,0.],[1.,1.],7),4,30)
        
        
        ##
        ##**************************************************
        ## Curve Coordinates
        Origin = np.asarray([0.,0.,0.])
        x, y, z = np.identity(3)
        self.ccs = Frame( np.asarray([x,y,z]) ) #curve coordinate system
        self.origin = DeepVector(v = Origin,
                                B=self.ccs)
        self.rots = Quaternion.from_v_theta(x, 0.) #null rotation vector
        
        ## axi of rotation
        self._ax_yaw        = yax           # self._ax_LR = yax #(0, -1, 0) #Left Right
        self._ax_alt_yaw_   = (0, 0, 1)     # self._ax_LR_alt = (0, 0, 1)
        self._ax_pitch      = xax           # self._ax_UD = xax #(1,0,0)
        self._ax_roll       = zax
        
        ##
        ##**************************************************
        ##
        

    def __str__(self):
        """print representation
        """
        keys = self.__dict__.keys()
        for key in keys:
            print key, self.__dict__[key]
        return ''
    
    
    def get_value(self, Q):
        if isinstance(Q, ad):
            return Q.value
        else:
            return Q
    
    def get_results(self):
        self.results_dict['Cb'] = self.SAC.area/(self.LengthDWL*\
                                                 self.MaxBreadth*\
                                                 self.MaxDraft)
        try:
            self.results_dict['Cmidshp'] = self.results_dict['Amsh']/(self.MaxBreadth*\
                                                                     self.MaxDraft)
        except:
            self.results_dict['Amsh'] = self.get_value( 
                                            2.*
                                            self.Lsplines.midSS[0].curve.area_to_y 
                                            )
            self.results_dict['Cmidshp'] = self.get_value( 
                                                self.results_dict['Amsh']/(
                                                        self.MaxBreadth*\
                                                        self.MaxDraft)
                                                )
        
        val = self.get_value(self.LengthDWL*\
                             self.MaxBreadth*\
                             self.MaxDraft*\
                             self.results_dict['Cmidshp'] )
        self.results_dict['Cp'] = self.get_value( self.SAC.area )/(val)
        
        self.results_dict['Cdl'] = self.get_value( self.SAC.area/(self.LengthDWL**3) )
        self.results_dict['Clb'] = self.get_value( self.LengthDWL/self.MaxBreadth )
        
        
        
        
        self.results_dict['Cwp'] = self.get_value(
                self.results_dict['Awp'] / 
                        (self.LengthDWL * self.MaxBreadth)
                )
        self.results_dict['Ccp'] = self.get_value(
                self.results_dict['Acp'] / 
                        (self.LengthDWL * self.MaxDraft)
                )
        
        self.results_dict['lwl'] = self.get_value(
                self.LengthDWL
                )
        self.results_dict['bwl'] = self.get_value(
                self.MaxBreadth
                )
        self.results_dict['draft'] = self.get_value(
                self.MaxDraft
                )
        
        return
    
    
    
    
    def check_nvertices(self, nrequired, nactual, mssg,op='=='):
        if op == '==':
            assert(nrequired==nactual),mssg.format(nrequired,nactual)
        elif op == '<':
            assert(nactual<nrequired),mssg.format(nrequired,nactual)
        elif op == '<=':
            assert(nactual<=nrequired),mssg.format(nrequired,nactual)
            
        return
    
    def change_hyperparameters(self, hyperparameters):
        if hyperparameters is not None: #putting this here to show where it can be accessed
            print 'altering hyperparameters'
            #
            self.aft_drop = hyperparameters.aft_drop
            #
            #dropout:
            # #[3]#-[3]is good!##trial #[3,7]#GOOD #[3] #bad #  [5]#bad, #[2,8] 
            #no good transition drop # [1,9] no good fairness drop
            self.dropout            = hyperparameters.dropout
            
            #
            self.remove_fwdfairness = hyperparameters.remove_fwdfairness
            
            #
            #drop aft fairness before linear interp
            self.heavy_config = hyperparameters.heavy_config
            
            #
            #drop aftmst FOSAC (not FOS!) before linear interp
            self.remove_aft_fos = hyperparameters.remove_aft_fos
            
            # [h,s,l] for each of [fwdFOS,FOSAC,aftFOS]
            self.midship_design_vector =  hyperparameters.midship_design_vector
            #
            #**********************************************************************
            #**********************************************************************
            # smaller is further fwd
            self.fwd_fairness_location_ratio    = hyperparameters.fwd_fairness_location_ratio
            
            # larger is closer to midship
            self.fwd_transition_location_ratio  = hyperparameters.fwd_transition_location_ratio
            #
            # smaller is closer to midship
            self.aft_transition_location_ratio  = hyperparameters.aft_transition_location_ratio
            
            #larger is closer to aft edge
            self.aft_fairness_location_ratio    = hyperparameters.aft_fairness_location_ratio
            #
            #
            #**********************************************************************
            #**********************************************************************
            # New Automation Friendly Parameters
            self.ignore = hyperparameters.ignore
            self.dropout = hyperparameters.dropout
            self.use_in_linear_interp = hyperparameters.use_in_linear_interp
            self.use_in_longi_solver = hyperparameters.use_in_longi_solver
        return
    
    def print_hyperparameters(self):
        print 'dropout'
        print self.dropout
        print 'ignore'
        print self.ignore
        #        print 'remove_fwdfairness'
        #        print self.remove_fwdfairness
        #        print 'heavy_config'
        #        print self.heavy_config
        #        print 'remove_aft_fos'
        #        print self.remove_aft_fos
        #        print 'midship_design_vector'
        #        print self.midship_design_vector
        print 'fwd_fairness_location_ratio'
        print self.fwd_fairness_location_ratio
        print 'fwd_transition_location_ratio'
        print self.fwd_transition_location_ratio
        print 'aft_transition_location_ratio'
        print self.aft_transition_location_ratio
        print 'aft_fairness_location_ratio'
        print self.aft_fairness_location_ratio
        return
        
        
    def set_list_of_hull_methods(self):
        """inspect generted list of hull design methods
        
        Auto generate this list with
        something like this:
        
        import inspect
        methods = inspect.getmembers(Hull, predicate=inspect.ismethod)
        
        for el in methods:
            print 'self.'+str(el[0])
        """
        self.methods = [self.get_maxSecArea,
                        self.get_Cb,
                        self.get_Cmidship,
                        self.get_LWL,
                        self.get_Vol,
                        self.get_bwl,
                        self.get_clArea,
                        self.get_draft,
                        self.get_flat_of_center_plane,
                        self.get_wlArea,
                        #self.insert_ctrl_vertices,
                        #self.issue_longitudinals,
                        #self.make_longitudinals,
                        #self.make_transverse_net,
                        #self.plot_primary_curves,
                        self.set_LCG,
                        #
                        # SAC curve , June 2016
#                        self.set_SAC_fwd_Xc,
#                        self.set_SAC_mid_Xc,
#                        self.set_SAC_run_Xc,
#                        self.get_SACareas,
                        self.get_len_flat_SAC,
                        self.get_len_flat_wL,
                        self.get_len_FOS, #uses self.FOS = self.flat_wL - self.flat_SAC directly
                        #
                        # SAC curve
                        #self.compute_SAC,
                        #self.compute_SAC_new,
                        #Jan 2018# self.get_len_entrance_SAC,
                        #Jan 2018# self.get_len_run_SAC,
                        #self.get_len_FOS, #uses self.FOS = self.flat_wL - self.flat_SAC directly
                        #self.get_transverse_start_end_param_loc,
                        #self.get_transverse_start_end_points,
                        #
                        #self.compute_DWL,
                        # Jan 2018# self.compute_FOSAC,
                        #self.compute_bow_fairing_curve,
                        #self.compute_cLProfile,
                        #self.split_hull_curves,
                        #self.compute_midship_section,
                        #self.compute_midship_section_old,
                        #self.compute_stern_bulb_profile,
                        #self.compute_stern_curve,
                        #self.compute_stern_fairing_curve,
                        #new guy:
                        #self.compute_stern_transition_curve()
                        #
                        #self.compute_stern_profile,
                        #self.compute_stern_transverse,
                        #self.correct_hull_curves,
                        #self.correct_longitudinal_match_hullcurves,
                        #self.definestations,
                        #self.export_rhino,
#                        self.get_loc_bfc,
#                        self.bow_fairness_section_Coefficient,
#                        self.bow_fairness_section_area,
#                        self.bow_fairness_section_beam,
#                        self.bow_fairness_section_draft,
#                        #self.computeFlatOfSide_curves_old,
#                        self.get_loc_sfc,
#                        self.stern_fairness_section_Coefficient,
#                        self.stern_fairness_section_area,
#                        self.stern_fairness_section_beam,
#                        self.stern_fairness_section_draft]
                        ]
        return
    
    def define_hull(self):
        """This triggersthe  functions
        which get the hull parameters from the design_space
        database
        """
        for el in self.methods:
            #print el.__doc__
            el()
        return
    
    def generate_list_of_hull_methods(self):
        """inspect generted list of hull design methods
        
        Auto generate this list with
        something like this:
        
        (working stuff (slightly different utility)
        in the simple_hull_rules_langage now)
        """
        import inspect
        methods = inspect.getmembers(Hull, predicate=inspect.ismethod)
        #print 'self.methods = ['
        #for el in methods:
        #    print 'self.'+str(el[0])
        #print ']'
        #eval?  or a safer meta needed?
        return methods

    def get_transverse_start_end_param_loc(self, loc_start, loc_end, 
                                           which=2):
        """
        all locations are longitudinal:
        i.e. get the longitudinal location some transverse's start and end
            loc_start   : real space longitudinal location start
            loc_end     : real lingi location end
            
        TODO:
            make tolerance of curve.findpoint
            less than tolerance of 
            midsection-form curve check tolerance
        """
        assert(isinstance(self.SAC, spline.Bspline))
        if self.SAC.dim ==2:
            sacwhich =0
        else:
            sacwhich = which
        assert(isinstance(self.DWL, spline.Bspline))
        assert(isinstance(self.CProfile, spline.Bspline))
        a = self.SAC.FindPoint(loc_start,sacwhich)       #2D
        #b = self.CProfile.FindPoint(loc_start,which)    #3D
        #c = self.DWL.FindPoint(loc_end,which)           #3D
        #
        #**********************************************************************
        # station finding:  which curve 
        if self.hullformtype == 'osv':
            c = self.DWL_OSV.FindPoint(loc_end,which)       #3D
            b = self.CPK_cdr_b.FindPoint(loc_start,which)    #3D
        else:
            c = self.DWL.FindPoint(loc_end,which)           #3D
            b = self.CPK_cdr_b.FindPoint(loc_start,which)    #3D
        return a,b,c


    def get_transverse_start_end_points(self, pa,pb):
        assert(0.0<=pa<=1.0),"pa !=[0,1] => {}".format(pa)
        pt1 = self.CProfile.CurvePoint(pa)  #start!
        pt2 = self.DWL.CurvePoint(pb)
        return pt1, pt2

    def get_LWL(self, frac=.5):
        """Unify and propogate constraints previously complete
        just set local
        max waterline length
        """
        self.LengthDWL = self.hullconstraints.state(
                                self.hullconstraints.lwl).getpoint(frac)
        if self.revise:
            self.hullconstraints.lwl = ia(self.LengthDWL-self.small,
                                          self.LengthDWL+self.small)#unify and propogate
            #self.definestations() #no -> do this after defining all basics!
            self.hullconstraints.AC_revise(print_=False)
        return

    def get_bwl(self, frac=.5):
        """Unify and propogate constraints previously complete
        just set local
        max waterline breadth
        """
        bwl = self.hullconstraints.state(
                                self.hullconstraints.bwl)
        if isinstance(bwl, ia):
            self.MaxBreadth = bwl.getpoint(frac)
        elif isinstance(bwl, float):
            self.MaxBreadth = bwl
        if self.revise:
            self.hullconstraints.bwl = bwl#ia(self.MaxBreadth-self.small,
                                          #self.MaxBreadth+self.small)
            self.hullconstraints.AC_revise(print_=False)
        else:
            self.MaxBreadth = bwl.getpoint(frac)
        return

    def get_draft(self, frac=.5):
        """Unify and propogate constraints previously complete
        just set local
        design waterline : draft
        """
        self.MaxDraft = self.hullconstraints.state(
                                self.hullconstraints.draft).getpoint(frac)
        if self.revise:
            self.hullconstraints.draft = ia(self.MaxDraft-self.small,
                                            self.MaxDraft+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return

    def get_Vol(self, frac=.5):
        """Unify and propogate constraints previously complete
        just set local
        design waterline : volume
        """
        self.BareHullVolume = self.hullconstraints.state(
                                self.hullconstraints.vol).getpoint(frac)
        if self.revise:
            self.hullconstraints.vol = ia(self.BareHullVolume-self.small ,
                                      self.BareHullVolume+self.small )
            
            self.hullconstraints.AC_revise(print_=False)
        return

    def get_wlArea(self, frac=.5):
        """Unify and propogate constraints previously complete
        just set local
        design waterline : area
        """
        self.WLarea = self.hullconstraints.state(
                                self.hullconstraints.Awp).getpoint(frac)
        if self.revise:
            self.hullconstraints.Awp = ia(self.WLarea-self.small,
                                          self.WLarea+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return

    def get_clArea(self, frac=.5):
        """Unify and propogate constraints previously complete
        just set local
        design waterline : area
        """
        self.CLarea = self.hullconstraints.state(
                                self.hullconstraints.Acp).getpoint(frac)
        if self.revise:
            #        self.hullconstraints.Acp = ia(self.CLarea-self.small,
            #                                      self.CLarea+self.small)
            self.hullconstraints.Acp = ia(0.,
                                          self.CLarea+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return

    def get_maxSecArea(self, frac=.5):
        """Unify and propogate constraints previously complete
        just set local
        design waterline : max section area

        Note: Amsh and maxSecArea are double sided!
        """
        maxSecArea = self.hullconstraints.state(
                                self.hullconstraints.Amsh).getpoint(frac)
        if self.revise:
            self.hullconstraints.Amsh = ia(maxSecArea-self.small,
                                              maxSecArea+self.small)
            #        Cmidshp = self.hullconstraints.state(
            #                                self.hullconstraints.Cmidshp).getpoint(frac)
            #        self.hullconstraints.Cmidshp = ia(Cmidshp-self.small,
            #                                          Cmidshp+self.small)
            self.maxSecArea = maxSecArea #* Cmidshp
            self.hullconstraints.AC_revise(print_=False)
        else:
            self.maxSecArea = maxSecArea
        return
    
    def get_SACareas(self, frac=.5):
        a1 = self.hullconstraints.state(
                            self.hullconstraints.SAC_entrance_area).getpoint(frac)
        a2 = self.hullconstraints.state(
                            self.hullconstraints.SAC_mid_area).getpoint(frac)
        a3 = self.hullconstraints.state(
                            self.hullconstraints.SAC_run_area).getpoint(frac)
        if self.revise:
            pass
        else:
            self.SAC_entrance_area  = a1
            self.SAC_mid_area       = a2
            self.SAC_run_area       = a3
        return

    def get_Cmidship(self, frac=.5):
        """Unify and propogate constraints previously complete
        just set local
        design waterline : max section area
        """
        self.Cmidship = self.hullconstraints.state(
                                self.hullconstraints.Cmidshp).getpoint(frac)
        if self.revise:
            self.hullconstraints.Cmidshp = ia(self.Cmidship-self.small,
                                              self.Cmidship+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return


    def get_Cb(self, frac=.5):
        """Unify and propogate constraints previously complete
        just set local
        design waterline : max section area
        """
        self.Cb = self.hullconstraints.state(
                                self.hullconstraints.Cb).getpoint(frac)
        if self.revise:
            self.hullconstraints.Cb = ia(self.Cb-self.small,
                                         self.Cb+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return


    def get_len_flat_wL(self, frac=.5):
        """Unify and propogate constraints previously complete
        just set local
        flat of wL area curve

        !! issue with inequality:
            flat of wL does not bound wL_area/bwl  !!
        """
        self.flat_wL = self.hullconstraints.state(
                                self.hullconstraints.lfwl).getpoint(frac)
        if self.revise:
            #self.hullconstraints.lfwl = ia(0.,self.flat_wL+self.small)
            self.hullconstraints.lfwl = ia(self.flat_wL-self.small,
                                           self.flat_wL+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return


    def get_len_FOS(self, frac=.5):
        """Unify and propogate constraints previously complete
        just set local
        flat of sac

        !! issue with inequality:
            flat of sac does not bound vol/(bwl*draft)  !!
        """
#        self.FOS = self.hullconstraints.state(
#                                self.hullconstraints.lfos).getpoint(frac)
        self.FOS = self.flat_wL - self.flat_SAC
        if self.revise:
            #self.hullconstraints.lfos = ia(0.,self.FOS+self.small)
            self.hullconstraints.lfos = ia(self.FOS-self.small,
                                           self.FOS+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return
    """
    #    @property
    #    def LCG(self):
    #        return self._LCG
    #    @LCG.setter
    #    def LCG(self, frac=.5):
    #        return self.hullconstraints.get_val(self.hullconstraints.LCG).getpoint(frac)
    #"""
    def set_LCG(self, frac=.5):
        #self.LCG = self.hullconstraints.get_val(
        #                    self.hullconstraints.LCG).getpoint(frac)
        
        self.LCG = self.hullconstraints.state(
                                self.hullconstraints.LCG).getpoint(frac)
        if self.revise:
            self.hullconstraints.LCG = ia(self.LCG-self.small,  
                                           self.LCG+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return
    def set_SAC_fwd_Xc(self, frac=2./3.):
        #self.SAC_fwd_Xc = self.hullconstraints.get_val(
        #                    self.hullconstraints.SAC_fwd_Xc).getpoint(frac)
        self.SAC_fwd_Xc = self.hullconstraints.state(
                            self.hullconstraints.SAC_fwd_Xc).getpoint(frac)
        if self.revise:
            self.hullconstraints.SAC_fwd_Xc = ia(self.SAC_fwd_Xc-self.small,
                                                 self.SAC_fwd_Xc+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return
    def set_SAC_mid_Xc(self, frac=.5):
        #self.SAC_mid_Xc = self.hullconstraints.get_val(
        #                    self.hullconstraints.SAC_mid_Xc).getpoint(frac)
        self.SAC_mid_Xc = self.hullconstraints.state(
                            self.hullconstraints.SAC_mid_Xc).getpoint(frac)
        if self.revise:
            self.hullconstraints.SAC_mid_Xc = ia(self.SAC_mid_Xc-self.small,
                                           self.SAC_mid_Xc+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return
    def set_SAC_run_Xc(self, frac=1./3.):
#        self.SAC_run_Xc = self.hullconstraints.get_val(
#                            self.hullconstraints.SAC_run_Xc).getpoint(frac)
        self.SAC_run_Xc = self.hullconstraints.state(
                            self.hullconstraints.SAC_run_Xc).getpoint(frac)
        if self.revise:
            self.hullconstraints.SAC_run_Xc = ia(self.SAC_run_Xc-self.small,
                                           self.SAC_run_Xc+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return



    def get_len_flat_SAC(self, frac=.5):
        """Unify and propogate constraints
        flat of SAC

        !! FIXED ?
            issue with inequality:
            flat of SAC does not bound vol  !!
        """
        self.flat_SAC = self.hullconstraints.state(
                                self.hullconstraints.lfsac).getpoint(frac)
        if self.revise:
            #self.hullconstraints.lfsac = ia(0.,self.flat_SAC+self.small)
            self.hullconstraints.lfsac = ia(self.flat_SAC-self.small,
                                            self.flat_SAC+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return
    def get_len_entrance_SAC(self, frac=.5):
        """Unify and propogate constraints
        entrance of SAC
        """
        self.SAC_entrance_len = self.hullconstraints.state(
                                self.hullconstraints.SAC_entrance_len).getpoint(frac)
        if self.revise:
            self.hullconstraints.SAC_entrance_len = ia(self.SAC_entrance_len-self.small,
                                                       self.SAC_entrance_len+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return
    def get_len_run_SAC(self, frac=.5):
        """Unify and propogate constraints
        entrance of SAC
        """
        self.SAC_run_len = self.hullconstraints.state(
                                self.hullconstraints.SAC_run_len).getpoint(frac)
        if self.revise:
            self.hullconstraints.SAC_run_len = ia(self.SAC_run_len-self.small,
                                                  self.SAC_run_len+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return

    def get_flat_of_center_plane(self, frac = .5):
        """
        """
        self.flat_bottom = self.hullconstraints.state(
                                self.hullconstraints.lfcp).getpoint(frac)
        if self.revise:
            self.hullconstraints.lfcp = ia(self.flat_bottom-self.small,
                                           self.flat_bottom+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return

    def get_loc_bfc(self, frac = .5):
        """
        """
        self.loc_bfc = self.hullconstraints.state(
                                self.hullconstraints.loc_bfc).getpoint(frac)
        if self.revise:
            self.hullconstraints.loc_bfc = ia(self.flat_bottom-self.small,
                                           self.flat_bottom+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return

    def get_loc_sfc(self, frac = .5):
        """
        """
        self.loc_sfc = self.hullconstraints.state(
                                self.hullconstraints.loc_sfc).getpoint(frac)
        if self.revise:
            self.hullconstraints.loc_sfc = ia(self.flat_bottom-self.small,
                                           self.flat_bottom+self.small)
            self.hullconstraints.AC_revise(print_=False)
        return


    def definestations(self,cp=.5):
        """Setup the hot spots on the hull
        
        About:
        ----------------
            uses the Stations class
            cp : fraction of nonflat length to send fwd
            
            cp split decides the relative length front and rear
            (should make simple rule for this!)
        
        TODO: redo this
        ----------------
            this is the weakest point in everything.
            
            Still needs work!  Harries sin/linear/etc
            other ideas?
            rules?
            how about code clarity??
        """
        if self._verbose: print 'hull_from_simple_designspace: definestations'
        assert(self.LengthDWL != 0.0),"LengthDWL = {}".format(self.LengthDWL)
        #        lfwl    = self.flat_wL
        #        plfwl   = lfwl/self.LengthDWL
        #
        #        lfos    = self.FOS
        #        plfos   = lfos/self.LengthDWL
        #
        #        lfsac   = self.flat_SAC
        #        plfsac  = lfsac/self.LengthDWL
        #cpfrac = (self.hullconstraints.state(self.hullconstraints.SAC_entrance_len)/
        #                self.hullconstraints.state(self.hullconstraints.SAC_run_len))/2.
        #cpfrac = cpfrac.getpoint(.5)
        #
        # January 2018 divy up the SAC without the over-specified 
        # cicvular, etc., rules
        #
        #**********************************************************************
        # length of ship
        lwl = self.LengthDWL
        #
        #**********************************************************************
        # length flat of sac
        lfsac = self.flat_SAC
        #
        #**********************************************************************
        # LCG may come into play again (3 curve looking at you)
        #lcg = self.LCG
        #
        #**********************************************************************
        # front and back allowables
        difflen = lwl - lfsac
        self.SAC_entrance_len = difflen*cp
        self.SAC_run_len =  difflen*(1.-cp)
        #
        #**********************************************************************
        # 
        # better to control this directly:
        #self.LCG = = cp/2.
        #
        #**********************************************************************
        #
        self.stations = Stations(self.LengthDWL, cp)
        #
        #**********************************************************************
        #  use stations
        self.stations.FOSAC = [self.flat_SAC, cp]
        #self.hullconstraints.SAC_entrance_len = ia(self.stations.FOSAC[0]-self.small,
        #                                           self.stations.FOSAC[0]+self.small)

        self.stations.FOWL  = [self.flat_wL, cp]
        self.stations.FOS   = [self.FOS, cp]
        self.stations.FOCP  = [self.flat_bottom,cp]

        #
        #lbow = (1./6.)*(self.stations.FOCP[0])
        #self.stations.Bow_Profile = [0.,lbow] #does nothing 
        #
        
        
        
        #
        #**********************************************************************
        # Grid Generation:  BOW Fairness Curve
        use = min(self.stations.FOCP[0],
                  self.stations.FOSAC[0],
                  self.stations.FOWL[0])
        #
        #**********************************************************************
        # Set fractions for bow fairning and transition to midship
        #
        #  location of the bow fairness curve
        #       ratio is the fraction of the 
        #           fwd section of (CPKeel) 
        #   or 
        #       fwd section of (SAC) 
        #       where the bow fairness curve shall be placed
        #
        #TODO: make a rule for this!
        #
        #ratio = 0.25   #->  larger ratio move the BFC aftwards #0.7
        #ratio = .3
        #ratio = 0.38
        #
        #ratio2 = .78         # second (further aft!) curve to 
                            # smooth the bow-bulb transition
        ratio = self.fwd_fairness_location_ratio
        ratio2 = self.fwd_transition_location_ratio
        #
        #**********************************************************************
        #  location of the bow fairness curve
        #       ratio is the fraction of the 
        #           fwd section of (CPKeel) 
        #
        locbowfair = (ratio)*use
        #
        """
        the trick here is that 
        this location
        does not know anything about where the 
        beginning of the curves section 
        of the nose will be (at keel)
        therefore you may inadvertedly
        pick a spot that is above
        the keel - bad for bulbous bow!
        """
        self.stations.Bow_Fairing = [locbowfair,
                                     locbowfair]
        # 
        # Done with BOW fairness curve
        #**********************************************************************
        # 'Rhino Bow-Bulb FWD Transition'
        locbowinterface = (ratio2)*use
        #
        self.stations.Bow_Bulb_Transition = [locbowinterface,
                                             locbowinterface]
        
        #
        #**********************************************************************
        # DONE with BOW curves
        
        #**********************************************************************
        #  Grid Generation:  STERN Fairness Curve
        use = max(self.stations.FOCP[1],
                  self.stations.FOSAC[1],
                  self.stations.FOWL[1])
        #
        #**********************************************************************
        #   stern fairness and midship transition ratios (station location)
        #
        #
        #ratio = .8        #-> larger ratio move the SFC aftwards #0.3
        #
        #ratio2 = .35         # second (yet further fwd) curve to smooth the stern
                            # OSV transom-like sections
        #
        #ratio = .75
        #ratio2 = .25
        #
        # back to old ratio
        #ratio = .75  #really nice hull seen here
        #ratio2 = .3 # with this pair (forgot to export)
        
        #match fwd for symmetric curves? 
        # this was tough...
        # ratio = .78  #
        # ratio2 = .38 # 
        #
        #ratio = .68  #
        #
        #ratio = .65  #
        #ratio2 = .28 # 
        ratio = self.aft_fairness_location_ratio
        ratio2 = self.aft_transition_location_ratio
        #
        #ratio = .6
        #ratio2 = .3 
        #ratio2 = .25 #bad
        #
        #
        #**********************************************************************
        #   location of the stern fairness curve
        locsternfair = (ratio)*(self.LengthDWL-use)
        locsternfair = use+locsternfair
        #
        self.stations.Stern_Fairing = [locsternfair,
                                       locsternfair]
        #
        #**********************************************************************
        # additional aft station (actually FWD of the Stern_Fairing):
        
        locsterntransition = (ratio2)*(self.LengthDWL-use)
        locsterntransition = use+locsterntransition
        
        self.stations.Stern_Transition = [locsterntransition,
                                          locsterntransition]
        #
        #**********************************************************************
        # Done with STERN curves
        
        #
        #**********************************************************************
        # DONE
        return
    
    



    def compute_midship_section_old(self):
        """For reference only
        """
        if self._verbose: print 'hull_from_simple_designspace: compute_midship_section_old'
        assert(False),'error, do not use this midship design routine!'
        Bmax = self.MaxBreadth/2.
        bottom  = linear_vertices((0.,0.,self.LengthDWL*0.5),(Bmax,0.,self.LengthDWL*0.5),4)
        side    = linear_vertices((Bmax,0.,self.LengthDWL*0.5),(Bmax,self.MaxDraft,self.LengthDWL*0.5),3)
        self.midshipsec = np.zeros((len(bottom)+len(side),3),float)
        self.midshipsec[:len(bottom)]=bottom
        self.midshipsec[len(bottom):]=side
        self.midshipsec = spline.Bspline(self.midshipsec,self.k,self.nump)
        return
    
    
    def compute_midship_section(self):
        """
        Builds midship curves 
        based on the real locations of 
        the bare hull curves of form.
        
        
        notes:
        ----------
        
        
        
        
        DEV
        ----------
        
        best way to use dropout:
            make all 3+2+3 midship curve every time
            -drop the ones you don't want, buy index
            -you can then LS them in the Longitudinals,
            should you choose
        
        worst way to use dropout:
            -maintain 11 transverse curves beforehand at all times
            -to 'move them around' actually make sure
            not to have to dropout any curves that will be
            very close to other curves
            -those tend to cause loops when dropped in the
            solver as the starting case.
        
        import curve as spline
        from ADILS import interval_bounds
        from FormParameter import FormParameterDict
        from ADILS import Lagrangian
        from ADILS import lagrangian_bounds
        
        from ADILS import IntervalLagrangeSpline
        
        from   adials_gui3d      import frame
        
        
        
        curve = SD.hull.SAC
        self = SD.hull
        
        curve = self.SAC
        
        curve.plot3DmultiList(SD.hull.lcurvenet,
                              SD.hull.initial_longitudinals)
        
        
        curve.plot3DmultiList([curve],
                              self.mshp_sections)
        
        
        curve.plot3DmultiList([],
                              self.mshp_sections)
        
        
        
        self.mshp_sections=[]
        self.Lsplines.midSS = []
        self.mshp_sections2D = []
        
        """
        print '\n compute_midship_section'
        if self._verbose: print 'hull_from_simple_designspace: compute_midship_section'
        local_TOL = 1.e-5
        k = self.k
        nump=self.nump
        area = self.maxSecArea/2.
        #*****************************************
        start, end, a, c = self.stations.FOSAC
        #*****************************************
        mid = (start+end)/2.
        b = (a+c)/2.
        stations = [start, mid, end] #using iterstations instead (mid is not the best use of the 3DOFs)
        fractions = [a,b,c]
        #
        # fwd flat of side
        wl_fwd = self.stations.FOWL[0]
        keel_fwd = self.stations.FOCP[0]
        #
        # aft flat of side
        wl_aft = self.stations.FOWL[1]
        keel_aft = self.stations.FOCP[1]
        #
        #********
        # use to make longitudinals more flat - only once you can 
        # come up with 11 total transverses! - got to get there for 
        # the longitudinals to avvoide having 
        # to fix-up the longitudinals to get to 11 with a good 
        # knot vector
        #"""
        #
        # iterstations : the box of 'straight up and down'
        # midship curves
        #
        diff = end-start
        iterstations = [start,
                        end]
        #-yes I was extending the back flat curve def
        # iterstations.insert(1, (end - .1*diff,))
        # iterstations.append(end + .1*diff) 
        #   just a hair - becuase after all it will
        #   not actually be flat there.
        #-meanwhile the fwd sections are going to have some help ;)
        #"""
        
        #iterstations = [start,mid,end]
        self.midship_total_grid = iterstations
        #
        #******
        #
        Bmax = self.MaxBreadth/2.
        Dmax = self.MaxDraft
        #
        self.mshp_sections2D = []
        self.mshp_sections   = []
        self.Lsplines.midSS  = []
        #loc=stations[0]
        #frac = fractions[0]
        ccc=-1
        
        
        def stage1(section,ab,ae):
            interval_data, small = interval_bounds(section)
            FPD = FormParameterDict(section)
            #-----------------------------------
            # AREA
            FPD.add_Area_to_y_Constraint(kind='equality',
                               value = area)
            """TODO: There is something in the CPKeel
            that is leaving the mad depth during
            the midship stations...
            """
            #attempt to force colinearity:
            #-----------------------------------
            ab = max(ab,0.)
            FPD.add_AngleConstraint(kind='LS',
                                    location = 0.,
                                    value = ab)
            FPD.add_AngleConstraint(kind='LS',
                                    location = 0.2,
                                    value = ab)
            #-----------------------------------
            # force convexity at the keel
            FPD.add_AngleConstraint(kind='min',
                                    location = 0.,
                                    value = 0.)
            #-----------------------------------
            FPD.add_verticalAngleConstraint(kind='LS',
                                            location = 1.,
                                            value = 90.-ae)
            # C2 continuity for flat side and bottom
            FPD.add_CurvatureConstraint(kind='equality',
                                    location = 0.,
                                    value = 0.)
            FPD.add_CurvatureConstraint(kind='equality',
                                    location = 1.0,
                                    value = 0.)
            #            FPD.add_relative_yVertexConstraint(
            #                                    kind = 'equality', location=None,
            #                                    value = 1., weight = 1.0,
            #                                    index = 2, index2 = 3,
            #                                    seperation = 0. )
            #            FPD.add_relative_xVertexConstraint(
            #                                    kind = 'equality', location=None,
            #                                    value = 1., weight = 1.0,
            #                                    index = 4, index2 = 5,
            #                                    seperation = 0. )
            #            FPD.add_relative_xVertexConstraint(
            #                                    kind = 'equality', location=None,
            #                                    value = 1., weight = 1.0,
            #                                    index = 3, index2 = 4,
            #                                    seperation = 0. )
            #            FPD.add_relative_yVertexConstraint(
            #                                    kind = 'min', location=None,
            #                                    value = 1., weight = 1.0,
            #                                    index = 4, index2 = 3,
            #                                    seperation = 0. )
            #
            #
            #**************************************************************************
            # fairness:
            #
            FPD.add_E1(kind='LS', weight = .5)
            FPD.add_E2(kind='LS', weight = .5)
            #FPD.add_E3(kind='LS', weight = .01)
            FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data,
                                                     small, 1.e4)
            Lspline = IntervalLagrangeSpline(section, L, data = interval_data)
            #
            #**************************************************************************
            # 
            print 'midship stage 1 solve:'
            Lspline.optimize(stop=35)
            return Lspline
        
        
        def stagefix(curve,ab,ae):
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve)
            #-----------------------------------
            # AREA
            FPD.add_Area_to_y_Constraint(kind='equality',
                               value = area)
            #-----------------------------------
            ab = max(ab,0.)
            FPD.add_AngleConstraint(kind='LS',
                                    location = 0.,
                                    value = ab)
            #-----------------------------------
            FPD.add_verticalAngleConstraint(kind='LS',
                                            location = 1.,
                                            value = 90.-ae)
            #-----------------------------------
            #
            #**************************************************************************
            # Fairness
            # note the hackish use of the E 'old' weights... what might they be?
            #FPD.add_E1(kind='LS', weight = .5/E1_0)
            #FPD.add_E2(kind='LS', weight = .5/E2_0)
            #FPD.add_E3(kind='LS', weight = .01/E3_0)
            #FPD.add_ArcLengthApprox(kind='LS', weight = 5./S_0)
            FPD.add_E1(kind='LS', weight = .5)
            FPD.add_E2(kind='LS', weight = .5)
            FPD.add_E3(kind='LS', weight = .01)
            FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
            #
            #**************************************************************************
            # 
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data,
                                                     small, 1.e4)
            Lspline = IntervalLagrangeSpline(section, L, data = interval_data)
            #
            #**************************************************************************
            # 
            print 'midship stage fixing solve:'
            Lspline.optimize(stop=35)
            #
            #**************************************************************************
            #Done
            return Lspline
        
        
        
        def stage2(Lspline,ab,ae,
                   pt1,pt2):
            #
            #**************************************************************************
            #Stage 2:
            #
            Lspline = copy.deepcopy(Lspline)
            Lspline.curve.pts_M_pts()
            Lspline.curve.compute_arclength()
            try:
                E1_0 = copy.deepcopy(Lspline.curve.E1.value)
            except:
                E1_0 = copy.deepcopy(Lspline.curve.E1)
            try:
                E2_0 = copy.deepcopy(Lspline.curve.E2.value)
            except:
                E2_0 = copy.deepcopy(Lspline.curve.E2)
            #try:
            #    E3_0 = copy.deepcopy(Lspline.curve.E3.value)
            #except:
            #    E3_0 = copy.deepcopy(Lspline.curve.E3)
            try:
                S_0 = copy.deepcopy(Lspline.curve.AL.value)
            except:
                S_0 = copy.deepcopy(Lspline.curve.AL)
            
            curve = Lspline.curve
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve)
            
            
            FPD = FormParameterDict(section)
            FPD.add_Area_to_y_Constraint(kind='equality',
                               value = area)
            #attempt to force colinearity:
            #-----------------------------------
            # bottom
            FPD.add_yFixity(index=1,
                            value=pt2[1])
            
            FPD.add_yFixity(index=2,
                            value=pt2[1])
            
            
            aa = Bmax*Dmax 
            # if aiming for Cmidship~=1, then this fixity is nice to have:
            if (area/aa) > .92:
                print 'sensing Cmship near 1, going for wider collinear bottom'
                FPD.add_yFixity(index=3,
                                value=pt2[1])
            #
            #-----------------------------------
            # side
            FPD.add_xFixity(index=4,
                            value=pt1[0])
            
            FPD.add_xFixity(index=5,
                            value=pt1[0])
            #-----------------------------------
            # force global convexity at the keel
            #-----------------------------------
            # C2 continuity for flat side and bottom
            #
            #******************************************************************
            # Fairness
            FPD.add_E1(kind='LS', weight = .5/E1_0)
            FPD.add_E2(kind='LS', weight = .5/E2_0)
            #FPD.add_E3(kind='LS', weight = .01/E3_0)
            FPD.add_ArcLengthApprox(kind='LS', weight = 1./S_0)
            #
            #******************************************************************
            # 
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data,
                                                     small, 1.e4)
            Lspline = IntervalLagrangeSpline(section, L, data = interval_data)
            #
            #*******************************************************************
            # 
            print 'midship stage 2 solve:'
            Lspline.optimize(stop=35)
            #
            #******************************************************************
            #Done
            return Lspline
        
        
        #
        #**********************************************************************
        #**********************************************************************
        # Loop and generate the midship curves ##lgmc##
        
        loc, frac = stations[0],fractions[0]
        ccc+=1
        #b = self.DWL.FindPoint(loc,2)
        #a = self.CProfile.FindPoint(loc,2) #reversed!
        #
        #
        #o=SAC a=CProfile  b=DWL
        o,a,b = self.get_transverse_start_end_param_loc(loc,loc)
        #
        #
        #o = self.SAC.FindPoint(loc_start,sacwhich)

        #assert(a[2] is True), "error: Could not locate midship match point on Keel Profile Curve"
        #assert(b[2] is True), "error: Could not locate midship match point on DWL"

        #pt2 = self.CProfile.CurvePoint(a[0]) #start
        #pt1 = self.DWL.CurvePoint(b[0]) #end
        pt1 = b[1] #DWL point
        pt2 = a[1] #cLProfile point
        
        #pt2,pt1=self.get_transverse_start_end_points(a[0],b[0])
        assert(abs(Dmax-pt1[1])<local_TOL), "Beam elevation does not match Profile: {} != {}".format(Dmax,pt1[1])
        
        #pt2,pt1=self.get_transverse_start_end_points(a[0],b[0])
        #assert(abs(Bmax-pt1[0])<local_TOL), "Beam width does not match itself: {} != {}".format(Bmax,pt1[1])
        
        if a[2] is False or b[2] is False or (abs(pt1[0]-Bmax)>local_TOL):
            print 'warning, sections may not match global max!' 
            
        #
        # it's trickier to look at CProfile, because the curve 
        #  has been 'flipped upside down' compared to the way it was
        #  generated.
        #
        # trickier still because CProfile is normalized
        #   to be '0' at the Keel!
        #
        #ckProfile = self.CProfile.vertices[0,1]-Dmax
        #assert(abs(ckProfile-pt2[1])<self.tol), "MaxDraft does not match Profile: {} != {}".format(Dmax,pt1[1])
        if a[2] is False or b[2] is False or (abs(pt1[0]-Bmax)>self.tol):
            print 'warning, sections may not match global max!' 
            

        bottom  = linear_vertices((pt2[0],pt2[1],loc),
                                  (pt1[0],pt2[1],loc),4)
        side    = linear_vertices((pt1[0],pt2[1],loc),
                                  (pt1[0],pt1[1],loc),3)#$
        section = np.zeros((len(bottom)+len(side),3),float)
        section[:len(bottom)]=bottom
        section[len(bottom):]=side
        section = spline.Bspline(section[:,0:2],self.k,self.nump)
        ab = np.rad2deg(section.compute_tangent(0.))
        ae = np.rad2deg(section.compute_tangent(1.))
        #
        #**************************************************************************
        # optimize based on DWL, SAC, and Body cL_Profile:
        """TODO: add bilge radius constraint"""
        """consider its automation with lim Amsh -> 1"""
        inicurve = copy.copy(section)
        
        
        #
        #**************************************************************************
        # Stage 1 
        Lspline = stage1(copy.deepcopy(section),
                         ab,ae)
        LSAVE1 = copy.deepcopy(Lspline)
        #**************************************************************************
        # Stage 2
        Lspline = stage2(copy.deepcopy(Lspline),
                         ab,ae,
                         pt1,pt2)
        #**************************************************************************
        #
        #if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
        #    Lspline = stage2(inicurve,ab,ae)
            #if LSAVE1.error_code is 'max_iter' or LSAVE1.error_code is 'NAN':
            #    Lspline = stagefix(inicurve,ab,ae)
            #else:
            #    Lspline = LSAVE1
        #
        #
        #**********************************************************************
        # FOSAC provides the sections
        #for loc, frac in zip(stations,fractions):
        #
        # start = iterstations[0] = fwd FOSAC
        # end   = iterstayions[1] = aft FOSAC
        #
        
        #if self.midship_design_vector[2] == 'h':
        #    #assert(self.remove_aft_fos),'Something must take away the extra curve'
        #    loop_this = iterstations
        #else:
        #    loop_this = [iterstations[0]]
        loop_this = iterstations
        for i,loc in enumerate(loop_this):
            #
            #
            #******************************************************************
            # store
            Lspline = copy.deepcopy(Lspline)
            #
            #******************************************************************
            #
            print 'midship curve done'
            self.Lsplines.midSS.append(Lspline)
            #
            section = copy.deepcopy(Lspline.curve)
            self.mshp_sections2D.append(section)
            #
            #******************************************************************
            # make 3D   
            vts = frame.add1D(section.vertices,2,where=loc)
            vts[:,2] = loc
            #section.vertices = vts #no - smart handler does not handle extra dimensions!
            section = spline.Bspline(vts, k, nump)
            self.mshp_sections.append(section)

            #
            #******************************************************************
            # tcurve map (hyperparameter control, 
            #                   possible graphical access to Lspline)
            self.tcurve_map[6+i] = {'Lspline':Lspline,
                                   'tcurve':section,
                                   'vertices':vts} #vts[:,:2] to input x,y only
        ##*********************************************************************
        # Go ahead and remove the midship curve ?- use only as a template
        #
        #*** tuning problems
        if self.midship_design_vector[2] == 's':
            Template_Midship = self.mshp_sections[0]
        elif self.midship_design_vector[2] == 'l': #light => no section at all here! FOSAC (implicitly fixed by FOS curves)
            Template_Midship = self.mshp_sections.pop()
        elif self.midship_design_vector[2] == 'h':
            assert(len(self.mshp_sections)==2),'ERROR, there should be 2 midship definition curves here'
            Template_Midship = self.mshp_sections[0]
        else:
            print 'ERROR, miship 2 option ',self.midship_design_vector[2], 'not supported, using standard'
            Template_Midship = self.mshp_sections[0]
        #***
        #
        # wl_fwd = self.stations.FOWL[0]
        # keel_fwd = self.stations.FOCP[0]
        #
        #
        #******************************************************************
        # Adding  Flat Of Side Section Stations 
        # -based off the midship curve
        #
        # -using 11 curves to do the basic design layout
        #
        # -but note later we could use THB to really enforce curvature flatness
        #   and anything else to do with local shape
        #
        diff = keel_fwd - wl_fwd
        #
        #******************************************************************
        #  FWD 
        #*** tuning problems
        #        if self.midship_design_vector[1] == 's': # standard => two forward leaning FOS curves
        #            longi_locs = [[keel_fwd, wl_fwd],
        #                          [keel_fwd-.1*diff, wl_fwd-.1*diff]]
        #        elif self.midship_design_vector[1] == 'h': # heavy => 3 fwd leaning curves (FOS)
        #            longi_locs = [[keel_fwd, wl_fwd],
        #                          [keel_fwd-.1*diff, wl_fwd-.1*diff],
        #                          [keel_fwd-.2*diff, wl_fwd-.2*diff]]
        #        elif self.midship_design_vector[1] == 'l': # light => just one fwd curve
        #            longi_locs = [ [keel_fwd, wl_fwd] ]
        #        else:
        #            print 'error, miship [1] option ',self.midship_design_vector[1], 'not supported, using standard'
        #            longi_locs = [[keel_fwd, wl_fwd],
        #                          [keel_fwd-.1*diff, wl_fwd-.1*diff]]
            
        #always make them all:
        longi_locs = [[keel_fwd, wl_fwd],
                      [keel_fwd-.1*diff, wl_fwd-.1*diff],
                      [keel_fwd-.2*diff, wl_fwd-.2*diff]]
        #***
        for i, pair in enumerate(longi_locs):
            keel_fwd, wl_fwd = pair[0],pair[1]
            
            curve = copy.deepcopy(Template_Midship)
            vertices = copy.deepcopy(curve.vertices)
            
            vertices[:,2] = keel_fwd #move the entire curve up to the 'keel_fwd' position
            vertices[-1,2] = wl_fwd  #wl_fwd is gauranteed to be fwd of the keel
                                     # just as the keel_fwd position is gauranteed to be
                                     # fwd of the flat of SAC positions
                                     # and the fwd fairness curve is gauranteed to 
                                     # be fwd of all the most fwd of these
            
            start = vertices[2,2]
            end = vertices[-1,2]
            #
            bottom =  vertices[2,1]
            top =  vertices[-1,1]
            #
            run = end-start
            rise = top-bottom
            for j in range(3,6):
                vstep = vertices[j,1]-vertices[j-1,1]
                vertices[j,2] = interpolate(run,rise,vstep,vertices[j-1,2])
            curve.vertices = vertices
            print 'flat of side curve done'
            self.mshp_sections.insert(0, curve)
            
            
            #
            #******************************************************************
            # tcurve map (hyperparameter control, 
            #                   possible graphical access to Lspline)
            print 'adding fwd sections to hyperparameter possibilities'
            print 'encapsulated in the tcurve_map'
            print 'tcurve_map [',i+3,']'
            self.tcurve_map[i+3] = {'Lspline':Lspline,
                                   'tcurve':curve,
                                   'vertices':vertices}
        #
        #******************************************************************
        #  AFT
        diff = wl_aft - keel_aft
        if diff<1.:
            assert(diff>0.),'ERROR: keel flat aft > wl flat aft location'
            diff = keel_fwd - wl_fwd
        #
        #*** tuning problems
        #        if self.midship_design_vector[3] == 's': # standard => two aft leaning (FOS) curves
        #            longi_locs = [[keel_aft, wl_aft],
        #                          [keel_aft+.1*diff, wl_aft+.1*diff] ] 
        #        elif self.midship_design_vector[3] == 'l': #light => 1 aft leaning (FOS) curves
        #            longi_locs = [ [keel_aft, wl_aft] ] 
        #        elif self.midship_design_vector[3] == 'h': # heavy => 3 
        #            longi_locs = [[keel_aft, wl_aft],
        #                          [keel_aft+.1*diff, wl_aft+.1*diff],
        #                          [keel_aft+.2*diff, wl_aft+.2*diff] ] 
        #        else:
        #            print 'error, miship option 3 ',self.midship_design_vector[1], 'not supported, using standard'
        #            longi_locs = [[keel_aft, wl_aft],
        #                          [keel_aft+.1*diff, wl_aft+.1*diff] ] 
        #***
        #always make them all, it costs almost nothing:
        longi_locs = [[keel_aft, wl_aft],
                      [keel_aft+.1*diff, wl_aft+.1*diff],
                      [keel_aft+.2*diff, wl_aft+.2*diff] ] 
        #
        for i,pair in enumerate(longi_locs):
            keel_aft, wl_aft = pair[0],pair[1]
            
            curve = copy.deepcopy(Template_Midship)
            vertices = copy.deepcopy(curve.vertices)
            
            vertices[:,2] = keel_aft #move the entire curve back to the 'keel_aft' position
            vertices[-1,2] = wl_aft  #wl_aft is gauranteed to be aft of the keel
                                     # just as the keel_fwd position is gauranteed to be
                                     # fwd of the flat of SAC positions
                                     # and the fwd fairness curve is gauranteed to 
                                     # be fwd of all the most fwd of these
            start = vertices[2,2]
            end = vertices[-1,2]
            #
            run = end-start
            rise = top-bottom
            for j in range(3,6):
                vstep = vertices[j,1]-vertices[j-1,1]
                vertices[j,2] = interpolate(run,rise,vstep,vertices[j-1,2])
            curve.vertices = vertices
            print 'flat of side curve done'
            self.mshp_sections.append(curve)
            #
            #******************************************************************
            # tcurve map (hyperparameter control, 
            #                   possible graphical access to Lspline)
            print 'adding fwd sections to hyperparameter possibilities'
            print 'encapsulated in the tcurve_map'
            print 'tcurve_map [',i+8,']'
            self.tcurve_map[8+i] = {'Lspline':Lspline,
                                   'tcurve':section,
                                   'vertices':vertices}
            
        self.results_dict['Amsh'] = Lspline.curve.area_to_y.value*2.
        return

    def computeFlatOfSide_curves_old(self, fore=27., aft = 45.):
        """
            The sectional area at some station where the flat of side curve is active
            must be at least the same as the boxed area
            under this curve at that station
        """
        Bmax = self.MaxBreadth/2.
        #def computeFlatOfSide_curves(self, self.stations.FOS[0], self.stations.FOS[1])
        fwd_extent = fore   #self.LengthDWL*fore
        aft_extent = aft    #self.LengthDWL*aft
        bottom      = linear_vertices((0.,0.,aft_extent),(Bmax,0.,aft_extent),4)
        side        = linear_vertices((Bmax,0.,aft_extent),(Bmax,self.MaxDraft,fwd_extent),3)
        self.flat_bottom = bottom
        self.flatside = {}
        self.flatside['fwd'] = fwd_extent
        #self.flatbottom['extent'] = [.45,.7]
        self.fos1   = np.zeros((len(bottom)+len(side),3),float)
        self.fos1[:len(bottom)]=bottom
        self.fos1[len(bottom):]=side
        self.fos1   = spline.Bspline(self.fos1,self.k,self.nump)
        #self.curve_by_station[str()]
        self.fos2 = Hullcurve('offset', self.fos1,2, -1., None)
        self.fos2 = self.fos2.return_curve
        self.fos3 = Hullcurve('offset', self.fos1,2, -2., None)
        self.fos3 = self.fos3.return_curve
        self.fos_fwd = [self.fos1,self.fos2,self.fos3]
        return


    def compute_FOSAC(self, frac=.5):
        #---- entrance length
        self.SAC_entrance_len = self.hullconstraints.state(
                                self.hullconstraints.SAC_entrance_len).getpoint(frac)
        if self.revise:
            self.hullconstraints.SAC_entrance_len = ia(self.SAC_entrance_len-self.small,
                                         self.SAC_entrance_len+self.small)
            self.hullconstraints.AC_revise(print_=False) #must run this now to ensure consistency fore n aft
        #---- run length
        self.SAC_run_len = self.hullconstraints.state(
                                self.hullconstraints.SAC_run_len).getpoint(frac)
        if self.revise:
            self.hullconstraints.SAC_run_len = ia(self.SAC_run_len-self.small,
                                         self.SAC_run_len+self.small)
    
            self.hullconstraints.AC_revise(print_=False)
        return

    def compute_stern_fairing_curve(self):
        """stations.Stern_Fairing
        Should be more correct than bow
        
        Csfc
        Asfc
        
        explanation of the getter:
            
        get_transverse_start_end_param_loc( loc_start, loc_end, =which=2)
            a = self.SAC.FindPoint(loc_start,sacwhich)
            b = self.CProfile.FindPoint(loc_start,which)
            c = self.DWL.FindPoint(loc_end,which)
            return a,b,c
            
            (a,b,c) <=> (o,a,b)
            stupid complication
        
        getting the z, or longitudinal, location in space
        
        
        NEW DATA to use:
            
        
            
            self.CPK_parametric_details['bow_split_point']
            
            self.CPK_parametric_details['stern_split_point'] 
            
            
            
            self.CProfile.plot3DmultiList([self.CProfile],
                                     [self.CPK_fwd,
                                     self.CPK_cdr_b,
                                     self.CPK_aft])
            
            self.CProfile.plot3DmultiList([self.CPK_cdr_b],
                                     [self.CPK_fwd,
                                     self.CPK_aft])
            
            
            
            self.CProfile.plot3DmultiList(
                                     [self.CPK_fwd,
                                     self.CPK_aft],
                                     [self.CPK_cdr_b])
            
            
            self.CProfile2D.plotcurve_detailed()
            self.CPK_fwd.plotcurve()
            self.CPK_cdr_b.plotcurve()
            self.CPK_aft.plotcurve()
            
            
            self.CProfile2D.t
            self.CPK_cdr_b.t
            
            
            self.CPK_fwd.t
            self.CPK_aft.t
            
        DEV
        ----------
        self = SD.hull
        import curve as spline
        
        
        from ADILS import interval_bounds
        from FormParameter import FormParameterDict
        from ADILS import Lagrangian
        from ADILS import lagrangian_bounds
        
        from ADILS import IntervalLagrangeSpline
        
        from   adials_gui3d      import frame
        
        
        
        if osv:
            DWL_aft
        else:
            CPK_cdr_b
        
        """
        print '\n hull_from_simple_designspace: compute_stern_fairing_curve\n'
        #
        #**********************************************************************
        # station finding:  which curve 
        if self.hullformtype == 'osv':
            dwl_station_curve = self.DWL_OSV
            #cLp_station_curve = self.CPK_cdr_b
        else:
            dwl_station_curve = self.DWL
            #cLp_station_curve = self.CPK_cdr_b
        #
        #**********************************************************************
        # 
        # o : self.SAC.FindPoint(loc_start,sacwhich) => z location
        # a : self.CProfile.FindPoint(loc_start,which) => z location
        # b : self.DWL.FindPoint(loc_end,which) => z location
        start,end,ao,bo = self.stations.Stern_Fairing
        o,a,b = self.get_transverse_start_end_param_loc(start,end)
        
        #stern_split_point = self.CPK_parametric_details['stern_split_point']
        
        #if a[0]==b[0]:
        if abs(a[1][2]-b[1][2])<self.tol:
            pass
        else:
            print 'stern fairning is not planar'
        if self.fairnesscurves_drive_hull_curves:
            print 'This was a bad idea'
            print 'could be an inequality from ia instead'
            area = self.Asfc/2. 
        else:
            if self._verbose: print 'deriving stern fairness curve area from SAC'
            area = self.SAC.CurvePoint(o[0])[1]*0.5
        
         
        #
        #**********************************************************************
        #  cLProfile - 
        #  Get the parametric location on the split curves
        # 
        #if stern_split_point<a[0]:
        #    print 'ERROR: fairness curves should never intersect transverses!!'
        #    a = self.CPK_aft.FindPoint(start,2)
        #elif stern_split_point>a[0]:
        a = self.CPK_cdr_b.FindPoint(start,2)
        
        
        
        #check if area will fit!
        dy = self.CProfile2D.CurvePoint(a[0])[1]
        dx = self.DWL2D.CurvePoint(b[0])[1]
        area_ck = dy*dx
        print 'area checking for stern fairness curve'
        print 'Area from SAC : ',area
        print 'Area convex with DWL, CProfile',area_ck
        print 'convex area? {}'.format(area_ck>area)
        if area_ck/(dy*dx) < 10.:
            small_area = True
        else:
            small_area = False
        
        #assert(a[2] is True), "error: Could not locate stern fairing match point on Keel Profile Curve"
        #assert(b[2] is True), "error: Could not locate stern fairing match point on DWL"
        if a[2] is False or b[2] is False:
            print 'warning Could not locate stern fairing match point(s)'
            print 'start = {}, end = {}'.format(start,end)
            print 'o = ',o
            print 'a = ',a
            print 'b = ',b
            
        
        pt1 = a[1]
        pt2 = b[1]
        bottom  = linear_vertices((pt1[0],pt1[1],pt1[2]),
                                      (pt2[0],pt1[1],pt1[2]),4)
        side    = linear_vertices((pt2[0],pt1[1],pt2[2]),
                                      (pt2[0],pt2[1],pt2[2]),3)#4
        section = np.zeros((len(bottom)+len(side),3),float)
        section[:len(bottom)]=bottom
        section[len(bottom):]=side
        section = spline.Bspline(section[:,0:2],self.k,self.nump) #inicurve
        ab = np.rad2deg(section.compute_tangent(0.))
        ae = np.rad2deg(section.compute_tangent(1.))
        #
        inicurve = copy.deepcopy(section)
        #
        def stage1(section, ab,ae):
            #
            #******************************************************************
            # Begin Solve
            interval_data, small = interval_bounds(section)
            FPD = FormParameterDict(section)
            #
            #******************************************************************
            #  Integral Parameters
            FPD.add_Area_to_y_Constraint(kind='equality',
                                         value = area)
            #
            #******************************************************************
            # Tangents at the ends
            #FPD.add_AngleConstraint(kind='min',
            #                        location = 0.,
            #                        value = ab)
            #FPD.add_verticalAngleConstraint(kind='LS',
            #                                location = 1.,
            #                                value = 90.-ae)
            FPD.add_xFixity(index=5,
                            value=FPD.curve.vertices[-1][0])
            FPD.add_xFixity(index=4,
                            value=FPD.curve.vertices[-1][0])
            FPD.add_yFixity(index=1,
                            value=FPD.curve.vertices[0][1])
            FPD.add_yFixity(index=2,
                            value=FPD.curve.vertices[0][1])
            
            #
            #******************************************************************
            # attempt to force colinearity:
            # force convexity at the keel
            mssg = 'ISSUE: planar area of stern fairness is '+\
                    '> area constraint '+\
                    'fix with a new rule!'
            assert(area_ck>area),mssg
            #if (area_ck>area):
                #for i in range(1,3):
                #    loc = section.greville_abscissa(i)
                #    FPD.add_CurvatureConstraint(kind='equality',
                #                                location = loc,
                #                                value = 0.)
            #
            #******************************************************************
            # Make sure the stern fairness curve does not penetrate the
            # centerplane or waterplane
            #            top = section.vertices[-1][1]
            #            for param in [.1,.2,.3,.7,.8,.9]:
            #                FPD.add_xPointConstraint(kind='min',
            #                                         location=param,
            #                                         value=0.,
            #                                         weight=.5)
            #                FPD.add_yPointConstraint(kind='max',
            #                                         location=param,
            #                                         value=top,
            #                                         weight=.5)
            #
            #******************************************************************
            # Fairness
            if small_area:
                print 'solving stern fairness curve'
                print '   with small cross section in mind'
                FPD.add_E1(kind='LS', weight = .1)
                FPD.add_E2(kind='LS', weight = .1)
                FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            else:
                FPD.add_E1(kind='LS', weight = .5)
                FPD.add_E2(kind='LS', weight = .5)
                #FPD.add_E3(kind='LS', weight = .5)
                FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            #
            #******************************************************************
            # 
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data,
                                                     small, 1.e4)
            Lspline = IntervalLagrangeSpline(section, L, data = interval_data)
            #
            #******************************************************************
            # Solve
            print 'solving stern fairness curve'
            Lspline.optimize()
            #
            #******************************************************************
            # Done
            #
            return Lspline
        
        def stagefix(section, ab,ae):
            #
            #******************************************************************
            # Begin Solve
            interval_data, small = interval_bounds(section)
            FPD = FormParameterDict(section)
            #
            #******************************************************************
            #  Integral Parameters
            FPD.add_Area_to_y_Constraint(kind='equality',
                                         value = area)
            #
            #******************************************************************
            # Tangents at the ends
            #FPD.add_AngleConstraint(kind='LS',
            #                        location = 0.,
            #                        value = ab)
            #FPD.add_verticalAngleConstraint(kind='LS',
            #                                location = 1.,
            #                                value = 90.-ae)
            FPD.add_xFixity(index=5,
                            value=FPD.curve.vertices[-1][0])
            FPD.add_xFixity(index=4,
                            value=FPD.curve.vertices[-1][0])
            FPD.add_yFixity(index=1,
                            value=FPD.curve.vertices[0][1])
            #
            #******************************************************************
            # Fairness
            if small_area:
                print 'solving stern fairness curve'
                print '   with small cross section in mind'
                FPD.add_E1(kind='LS', weight = .1)
                FPD.add_E2(kind='LS', weight = .1)
                FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            else:
                FPD.add_E1(kind='LS', weight = .5)
                FPD.add_E2(kind='LS', weight = .5)
                FPD.add_E3(kind='LS', weight = .1)
                FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
            #
            #******************************************************************
            # 
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data,
                                                     small, 1.e4)
            Lspline = IntervalLagrangeSpline(section, L, data = interval_data)
            #
            #******************************************************************
            # Solve
            print 'solving stern fairness curve'
            Lspline.optimize()
            #
            #******************************************************************
            # Done
            #
            return Lspline
        
        Lspline = stage1(copy.deepcopy(section),ab,ae)  
        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
            Lspline = stagefix(inicurve,ab,ae)  
    
        # save
        self.Lsplines.sternfairing = Lspline
        section = Lspline.curve
        #self.bowfairning.append(section)
        vts = frame.add1D(section.vertices,2,
                          where=pt1[2])
        section = spline.Bspline(vts, self.k, self.nump)
        self.sternfairing = section
        #
        #******************************************************************
        # tcurve map (hyperparameter control, 
        #                   possible graphical access to Lspline)
        self.tcurve_map[12] = {'Lspline':Lspline,
                               'tcurve':section}
        return
    
    
    def compute_stern_transition_curve(self):
        """stations.Stern_Transition
        
        Csfc
        Asfc
        
        explanation of the getter:
            
        get_transverse_start_end_param_loc( loc_start, loc_end, =which=2)
            a = self.SAC.FindPoint(loc_start,sacwhich)
            b = self.CProfile.FindPoint(loc_start,which)
            c = self.DWL.FindPoint(loc_end,which)
            return a,b,c
            
            (a,b,c) <=> (o,a,b)
            stupid complication
        
        getting the z, or longitudinal, location in space
        
        
        NEW DATA to use:
            
        
            
            self.CPK_parametric_details['bow_split_point']
            
            self.CPK_parametric_details['stern_split_point'] 
            
            
            
            self.CProfile.plot3DmultiList([self.CProfile],
                                     [self.CPK_fwd,
                                     self.CPK_cdr_b,
                                     self.CPK_aft])
            
            self.CProfile.plot3DmultiList([self.CPK_cdr_b],
                                     [self.CPK_fwd,
                                     self.CPK_aft])
            
            
            
            self.CProfile.plot3DmultiList(
                                     [self.CPK_fwd,
                                     self.CPK_aft],
                                     [self.CPK_cdr_b])
            
            
            self.CProfile2D.plotcurve_detailed()
            self.CPK_fwd.plotcurve()
            self.CPK_cdr_b.plotcurve()
            self.CPK_aft.plotcurve()
            
            
            self.CProfile2D.t
            self.CPK_cdr_b.t
            
            
            self.CPK_fwd.t
            self.CPK_aft.t
            
        DEV
        ----------
        self = SD.hull
        import curve as spline
        
        
        from ADILS import interval_bounds
        from FormParameter import FormParameterDict
        from ADILS import Lagrangian
        from ADILS import lagrangian_bounds
        
        from ADILS import IntervalLagrangeSpline
        
        from   adials_gui3d      import frame
        
        """
        print '\n hull_from_simple_designspace: compute_stern_transition_curve\n'
        #
        #**********************************************************************
        # 
        # o : self.SAC.FindPoint(loc_start,sacwhich) => z location
        # a : self.CProfile.FindPoint(loc_start,which) => z location
        # b : self.DWL.FindPoint(loc_end,which) => z location
        start,end,ao,bo = self.stations.Stern_Transition
        o,a,b = self.get_transverse_start_end_param_loc(start,end)
        
        #stern_split_point = self.CPK_parametric_details['stern_split_point']
        
        #if a[0]==b[0]:
        if abs(a[1][2]-b[1][2])<self.tol:
            pass
        else:
            print 'stern Transition is not planer'
        if self.fairnesscurves_drive_hull_curves:
            print 'This was a bad idea'
            print 'could be an inequality from ia instead'
            area = self.Asfc/2. 
        else:
            if self._verbose: print 'deriving stern fairness curve area from SAC'
            area = self.SAC.CurvePoint(o[0])[1]*0.5
        
         
        #
        #**********************************************************************
        #  cLProfile - 
        #  Get the parametric location on the split curves
        # 
        #if stern_split_point<a[0]:
        #    print 'ERROR: Transition curve should never intersect other transverses!!'
        #    a = self.CPK_aft.FindPoint(start,2)
        #elif stern_split_point>a[0]:
        a = self.CPK_cdr_b.FindPoint(start,2)
        
        
        
        #check if area will fit!
        dy = self.CProfile2D.CurvePoint(a[0])[1]
        dx = self.DWL2D.CurvePoint(b[0])[1]
        area_ck = dy*dx
        print 'area checking for stern transition curve'
        print 'Area from SAC : ',area
        print 'Area convex with DWL, CProfile',area_ck
        print 'convex area? {}'.format(area_ck>area)
        if area_ck/(dy*dx) < 10.:
            small_area = True
        else:
            small_area = False
        
        #assert(a[2] is True), "error: Could not locate stern fairing match point on Keel Profile Curve"
        #assert(b[2] is True), "error: Could not locate stern fairing match point on DWL"
        if a[2] is False or b[2] is False:
            print 'warning Could not locate stern transition match point(s)'
            print 'start = {}, end = {}'.format(start,end)
            print 'o = ',o
            print 'a = ',a
            print 'b = ',b
            
        
        pt1 = a[1]
        pt2 = b[1]
        bottom  = linear_vertices((pt1[0],pt1[1],pt1[2]),
                                      (pt2[0],pt1[1],pt1[2]),4)
        side    = linear_vertices((pt2[0],pt1[1],pt2[2]),
                                      (pt2[0],pt2[1],pt2[2]),3)#4
        section = np.zeros((len(bottom)+len(side),3),float)
        section[:len(bottom)]=bottom
        section[len(bottom):]=side
        section = spline.Bspline(section[:,0:2],self.k,self.nump) #inicurve
        ab = np.rad2deg(section.compute_tangent(0.))
        ae = np.rad2deg(section.compute_tangent(1.))
        #
        inicurve = copy.deepcopy(section)
        #
        def stage1(section, ab,ae):
            #
            #******************************************************************
            # Begin Solve
            interval_data, small = interval_bounds(section)
            FPD = FormParameterDict(section)
            #
            #******************************************************************
            #  Integral Parameters
            FPD.add_Area_to_y_Constraint(kind='equality',
                                         value = area)
            #
            #******************************************************************
            # Tangents at the ends
            #FPD.add_AngleConstraint(kind='min',
            #                        location = 0.,
            #                        value = ab)
            #FPD.add_verticalAngleConstraint(kind='LS',
            #                                location = 1.,
            #                                value = 90.-ae)
            FPD.add_xFixity(index=5,
                            value=FPD.curve.vertices[-1][0])
            FPD.add_xFixity(index=4,
                            value=FPD.curve.vertices[-1][0])
            FPD.add_yFixity(index=1,
                            value=FPD.curve.vertices[0][1])
            FPD.add_yFixity(index=2,
                            value=FPD.curve.vertices[0][1])
            
            #
            #******************************************************************
            # attempt to force colinearity:
            # force convexity at the keel
            mssg = 'ISSUE: planar area of stern fairness is '+\
                    '> area constraint '+\
                    'fix with a new rule!'
            assert(area_ck>area),mssg
            #if (area_ck>area):
            #    for i in range(1,3):
            #        loc = section.greville_abscissa(i)
            #        FPD.add_CurvatureConstraint(kind='equality',
            #                                    location = loc,
            #                                    value = 0.)
            #
            #******************************************************************
            # Make sure the stern fairness curve does not penetrate the
            # centerplane or waterplane
            top = section.vertices[-1][1]
            #            for param in [.1,.2,.3,.7,.8,.9]:
            #                FPD.add_xPointConstraint(kind='min',
            #                                         location=param,
            #                                         value=0.,
            #                                         weight=.5)
            #                FPD.add_yPointConstraint(kind='max',
            #                                         location=param,
            #                                         value=top,
            #                                         weight=.5)
            #
            #******************************************************************
            # Fairness
            if small_area:
                print 'solving stern transition curve'
                print '   with small cross section in mind'
                FPD.add_E1(kind='LS', weight = .1)
                FPD.add_E2(kind='LS', weight = .1)
                FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            else:
                FPD.add_E1(kind='LS', weight = .5)
                FPD.add_E2(kind='LS', weight = .5)
                #FPD.add_E3(kind='LS', weight = .5)
                FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            #
            #******************************************************************
            # 
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data,
                                                     small, 1.e4)
            Lspline = IntervalLagrangeSpline(section, L, data = interval_data)
            #
            #******************************************************************
            # Solve
            print 'solving stern transition curve'
            Lspline.optimize()
            #
            #******************************************************************
            # Done
            #
            return Lspline
        
        def stagefix(section, ab,ae):
            #
            #******************************************************************
            # Begin Solve
            interval_data, small = interval_bounds(section)
            FPD = FormParameterDict(section)
            #
            #******************************************************************
            #  Integral Parameters
            FPD.add_Area_to_y_Constraint(kind='equality',
                                         value = area)
            #
            #******************************************************************
            # Tangents at the ends
            FPD.add_AngleConstraint(kind='LS',
                                    location = 0.,
                                    value = ab)
            #FPD.add_verticalAngleConstraint(kind='LS',
            #                                location = 1.,
            #                                value = 90.-ae)
            FPD.add_xFixity(index=5,
                            value=FPD.curve.vertices[-1][0])
            FPD.add_xFixity(index=4,
                            value=FPD.curve.vertices[-1][0])
            FPD.add_yFixity(index=1,
                            value=FPD.curve.vertices[0][1])
            #
            #******************************************************************
            # Fairness
            if small_area:
                print 'solving stern transition curve'
                print '   with small cross section in mind'
                FPD.add_E1(kind='LS', weight = .1)
                FPD.add_E2(kind='LS', weight = .1)
                FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            else:
                FPD.add_E1(kind='LS', weight = .5)
                FPD.add_E2(kind='LS', weight = .5)
                FPD.add_E3(kind='LS', weight = .1)
                FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
            #
            #******************************************************************
            # 
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data,
                                                     small, 1.e4)
            Lspline = IntervalLagrangeSpline(section, L, data = interval_data)
            #
            #******************************************************************
            # Solve
            print 'solving stern transition curve'
            Lspline.optimize()
            #
            #******************************************************************
            # Done
            #
            return Lspline
        
        Lspline = stage1(copy.deepcopy(section),ab,ae)  
        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
            Lspline = stagefix(inicurve,ab,ae)  
    
        # save
        self.Lsplines.sterntransition = Lspline
        section = Lspline.curve
        #self.bowfairning.append(section)
        vts = frame.add1D(section.vertices,2,
                          where=pt1[2])
        section = spline.Bspline(vts, self.k, self.nump)
        self.sterntransition = section
        #
        #******************************************************************
        # tcurve map (hyperparameter control, 
        #                   possible graphical access to Lspline)
        self.tcurve_map[11] = {'Lspline':Lspline,
                               'tcurve':section}
        return
    
    
    def compute_bow_fairing_curve(self):
        """stations.Bow_Profile
        
        -now this one can
        be responsible for fairning into the
        bulbous bow.
        
        Cbfc
        Abfc
        
        
        
        DEV
        ----------
        import curve as spline
        from ADILS import interval_bounds
        from FormParameter import FormParameterDict
        from ADILS import Lagrangian
        from ADILS import lagrangian_bounds
        
        from ADILS import IntervalLagrangeSpline
        
        from   adials_gui3d      import frame
        
        
        
        
        explanation of the getter:
            
        get_transverse_start_end_param_loc( loc_start, loc_end, =which=2)
            a = self.SAC.FindPoint(loc_start,sacwhich)
            b = self.CProfile.FindPoint(loc_start,which)
            c = self.DWL.FindPoint(loc_end,which)
            return a,b,c
            
            (a,b,c) <=> (o,a,b)
            stupid complication
        
        
        NEW DATA to use:
            
        
            
            self.CPK_parametric_details['bow_split_point']
            
            self.CPK_parametric_details['stern_split_point'] 
            
            
            
            self.SAC.plot3DmultiList([self.CProfile],
                                     [self.CPK_fwd,
                                     self.CPK_cdr_b,
                                     self.CPK_aft])
        
        """
        print '\n hull_from_simple_designspace: compute_bow_fairing_curve \n'
        #
        #**********************************************************************
        # 
        if self.hullformtype == 'osv':
            dwlcurve= self.DWL_OSV  #3D
            bcpkcurve = self.CPK_cdr_b  #3D
        else:
            dwlcurve = self.DWL          #3D
            bcpkcurve = self.CPK_cdr_b    #3D
        #
        #**********************************************************************
        # 
        # o : self.SAC.FindPoint(loc_start,sacwhich) => z location
        # a : self.CProfile.FindPoint(loc_start,which) => z location
        # b : self.DWL.FindPoint(loc_end,which) => z location
        start,end,ao,bo = self.stations.Bow_Fairing
        o,a,b = self.get_transverse_start_end_param_loc(start,end)
        
        bow_split_point = self.CPK_parametric_details['bow_split_point']
        

        if self.fairnesscurves_drive_hull_curves:
            print 'This was a bad idea'
            print 'could be an inequality from ia instead'
            area = self.Abfc/2.
        else:
            if self._verbose: print 'deriving bow fairness curve area from SAC'
            area = self.SAC.CurvePoint(o[0])[1]*0.5 #ao)/2.#ao doesnt work as SAC has been reparameterized
        
            
        #check if area will fit
        dy = self.CProfile2D.CurvePoint(a[0])[1]
        dx = self.DWL2D.CurvePoint(b[0])[1]
        area_ck = dy*dx
        print 'area checking for bow fairness curve'
        print 'Area from SAC : ',area
        print 'Area convex with DWL, CProfile',area_ck
        print 'convex area? {}'.format(area_ck>area)
        if area_ck/(dy*dx) < 10.:
            small_area = True
        else:
            small_area = False
        
        if a[2] is False or b[2] is False:
            print 'warning Could not locate bow fairing match point(s)'
            print 'start = {}, end = {}'.format(start,end)
            print 'o = ',o
            print 'a = ',a
            print 'b = ',b
        
        a = a[1]
        #a = [a[0],self.MaxDraft,a[2]] #setting it flush to bottom since we will have a bulbous bow in the end.
        a = [a[0],0.,a[2]] 
        pt1 = a
        pt2 = b[1]
        Total_Height = pt2[1]-pt1[1]
        bottom  = linear_vertices((pt1[0],pt1[1],pt1[2]),
                                      (pt2[0],pt1[1],pt1[2]),4)
        side    = linear_vertices((pt2[0],pt1[1],pt2[2]),
                                      (pt2[0],pt2[1],pt2[2]),3)#4
        section = np.zeros((len(bottom)+len(side),3),float)
        section[:len(bottom)]=bottom
        section[len(bottom):]=side
        section = spline.Bspline(section[:,0:2],self.k,self.nump)
        ab = np.rad2deg(section.compute_tangent(0.))
        ae = np.rad2deg(section.compute_tangent(1.))
        #
        #**********************************************************************
        # special functions - circular bulb section of a curve
        circle = spfuncs.circle
        #ellipse = spfuncs.ellipse
        gfunc = spfuncs.gfunc
        #
        #radius = Total_Height*.15  #2.25 # assumed radius of bulbous bow
        radius_bow = np.sqrt(self.A_mid/np.pi) #radius tied to bbow
        radius_frac = np.sqrt(2.*area*self.hyperparameter_input.bulb_volume_fraction/np.pi) #independent, 
                                                #related to bbow by transition to bbow_SAC curve
        
        #Cblend = radius_bow/radius_frac
        if radius_bow*1.25<=radius_frac:
            #radius = radius_bow*1.25
            radius = radius_bow*1.18
        else:
            radius = radius_frac
        self.bowbulb_blending_raidus = radius
            
        self.fwd_fairness_blendingbulb_area     = area
        self.fwd_fairness_blendingbulb_radius   = radius
        #rb_ = 1.8   # elliptic foci
        #ra_ = 2.25  # elliptic foci
        #
        nump=30
        mp = .4 #
        org = pt1[1] + radius
        c1 = circle(origin=[pt1[1],org],
                    radius=radius)
        #        e1 = ellipse(origin=[0.,rb_],
        #                     ra=ra_, 
        #                     rb=rb_)
        circ_p = np.linspace(-np.pi/2., 
                             0., 
                             nump, 
                             endpoint=False)
        gfc3 = gfunc(copy.deepcopy(section), 
                     c1.circle_pt, 0.,mp,
                     circ_p, 10.)
        #
        #**********************************************************************
        # 
        inicurve = copy.copy(section)
        #
        def stage1(section,ab,ae):
            #
            #******************************************************************
            # Begin Solve
            interval_data, small = interval_bounds(section)
            FPD = FormParameterDict(section)
            #
            #******************************************************************
            #  Integral Parameters
            FPD.add_Area_to_y_Constraint(kind='equality',
                                         value = area)
            #
            #******************************************************************
            # Tangents at the ends
#            FPD.add_AngleConstraint(kind='min',
#                                    location = 0.,
#                                    value = ab)
#            FPD.add_verticalAngleConstraint(kind='LS',
#                                            location = 1.,
#                                            value = 90.-ae)
            FPD.add_yFixity(index=1,
                            value=FPD.curve.vertices[0][0])
            FPD.add_yFixity(index=2,
                            value=FPD.curve.vertices[0][0])
            #FPD.add_xFixity(index=5,
            #                value=FPD.curve.vertices[-1][0])
            #FPD.add_xFixity(index=4,
            #                value=FPD.curve.vertices[-1][0])
            #
            #******************************************************************
            # Make sure the bow fairness curve does not penetrate 
            # the centerplane or waterplane
            top = section.vertices[-1][1]
            for param in [.1,.2,.3,.7,.8,.9]:
                FPD.add_xPointConstraint(kind='min',
                                         location=param,
                                         value=0.,
                                         weight=.5)
                FPD.add_yPointConstraint(kind='max',
                                         location=param,
                                         value=top,
                                         weight=.5)
            #
            #******************************************************************
            # Fairness
            if small_area:
                print 'solving bow fairness curve'
                print '   with small cross section in mind'
                FPD.add_E1(kind='LS', weight = .1)
                FPD.add_E2(kind='LS', weight = .1)
                FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            else:
                FPD.add_E1(kind='LS', weight = .5)
                FPD.add_E2(kind='LS', weight = .5)
                #FPD.add_E3(kind='LS', weight = .5)
                FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            #
            #******************************************************************
            # 
            L = Lagrangian(FPD)
            #
            #******************************************************************
            # 
            print 'adding bulb interface objective to the functional'
            L.sp['circ']     = gfc3
            #print 'not adding bulb interface objective to the functional'
            #
            #******************************************************************
            # 
            interval_data, small = lagrangian_bounds(L, interval_data,
                                                     small, 1.e4)
            Lspline = IntervalLagrangeSpline(section, L, data = interval_data)
            #
            #******************************************************************
            # Solve
            Lspline.optimize()
            #
            #******************************************************************
            # Done
            return Lspline
        
        
        
        def stagefix(section, ab,ae):
            #
            #******************************************************************
            # Begin Solve
            interval_data, small = interval_bounds(section)
            FPD = FormParameterDict(section)
            #
            #******************************************************************
            #  Integral Parameters
            FPD.add_Area_to_y_Constraint(kind='equality',
                                         value = area)
            #
            #******************************************************************
            # Tangents at the ends
#            FPD.add_AngleConstraint(kind='LS',
#                                    location = 0.,
#                                    value = ab)
#            FPD.add_verticalAngleConstraint(kind='LS',
#                                            location = 1.,
#                                            value = 90.-ae)
            FPD.add_yFixity(index=1,
                            value=FPD.curve.vertices[0][0])
            FPD.add_yFixity(index=2,
                            value=FPD.curve.vertices[0][0])
            #FPD.add_xFixity(index=5,
            #                value=FPD.curve.vertices[-1][0])
            #FPD.add_xFixity(index=4,
            #                value=FPD.curve.vertices[-1][0])
            #
            #******************************************************************
            # Fairness
            if small_area:
                print 'solving bow fairness curve'
                print '   with small cross section in mind'
                FPD.add_E1(kind='LS', weight = .1)
                FPD.add_E2(kind='LS', weight = .1)
                FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            else:
                FPD.add_E1(kind='LS', weight = .5)
                FPD.add_E2(kind='LS', weight = .5)
                FPD.add_E3(kind='LS', weight = .1)
                FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
            #
            #******************************************************************
            # 
            L = Lagrangian(FPD)
            #
            #******************************************************************
            # 
            print 'adding bulb interface objective to the functional'
            L.sp['circ']     = gfc3
            #print 'not adding bulb interface objective to the functional'
            #
            #******************************************************************
            # 
            interval_data, small = lagrangian_bounds(L, interval_data,
                                                     small, 1.e4)
            Lspline = IntervalLagrangeSpline(section, L, data = interval_data)
            #
            #******************************************************************
            # Solve
            Lspline.optimize()
            #
            #******************************************************************
            # Done
            return Lspline
            
        
        Lspline = stage1(section,ab,ae)  
        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
            Lspline = stagefix(inicurve,ab,ae)  
        #
        # save
        #
        self.Lsplines.bowfairning = Lspline
        section = Lspline.curve
        #self.bowfairning.append(section)
        vts = frame.add1D(section.vertices,2,
                          where=pt1[2])
        section = spline.Bspline(vts, self.k, self.nump)
        self.bowfairning = section
        #
        #******************************************************************
        # tcurve map (hyperparameter control, 
        #                   possible graphical access to Lspline)
        self.tcurve_map[1] = {'Lspline':Lspline,
                               'tcurve':section}
        return



    def compute_bow_transition_curve(self):
        """stations.compute_bow_transition_curve
        
        Cbfc
        Abfc
        
        
        
        DEV
        ----------
        import curve as spline
        from ADILS import interval_bounds
        from FormParameter import FormParameterDict
        from ADILS import Lagrangian
        from ADILS import lagrangian_bounds
        
        from ADILS import IntervalLagrangeSpline
        
        from   adials_gui3d      import frame
        
        
        
        
        explanation of the getter:
            
        get_transverse_start_end_param_loc( loc_start, loc_end, =which=2)
            a = self.SAC.FindPoint(loc_start,sacwhich)
            b = self.CProfile.FindPoint(loc_start,which)
            c = self.DWL.FindPoint(loc_end,which)
            return a,b,c
            
            (a,b,c) <=> (o,a,b)
            stupid complication
        
        
        plotting for test:
            
            
            self.SAC.plot3DmultiList([self.CProfile],
                                     [self.CPK_fwd,
                                     self.CPK_cdr_b,
                                     self.CPK_aft])
                
        idea:
            pull the top forward --somehow--
            (constraint plus... make the hull a bit different)
            
            to smooth the flow to midship
            
            ... Actually this is kinda fixed by making
            the bow-midship transition C1 instead of C2?!
        
        """
        print '\n hull_from_simple_designspace: compute_bow_transition_curve \n'
        #
        #**********************************************************************
        # 
        if self.hullformtype == 'osv':
            dwlcurve= self.DWL_OSV  #3D
            bcpkcurve = self.CPK_cdr_b  #3D
        else:
            dwlcurve = self.DWL          #3D
            bcpkcurve = self.CPK_cdr_b    #3D
        #
        #**********************************************************************
        # 
        # o : self.SAC.FindPoint(loc_start,sacwhich) => z location
        # a : self.CProfile.FindPoint(loc_start,which) => z location
        # b : self.DWL.FindPoint(loc_end,which) => z location
        start,end,ao,bo = self.stations.Bow_Bulb_Transition
        o,a,b = self.get_transverse_start_end_param_loc(start,end)
        
        bow_split_point = self.CPK_parametric_details['bow_split_point']
        
        
        if self.fairnesscurves_drive_hull_curves:
            print 'This was a bad idea'
            print 'could be an inequality from ia instead'
            area = self.Abfc/2.
        else:
            if self._verbose: print 'deriving bow fairness curve area from SAC'
            area = self.SAC.CurvePoint(o[0])[1]*0.5 #ao)/2.#ao doesnt work as SAC has been reparameterized
        
            
        #check if area will fit
        dy = self.CProfile2D.CurvePoint(a[0])[1]
        dx = self.DWL2D.CurvePoint(b[0])[1]
        area_ck = dy*dx
        print 'area checking for bow transition curve'
        print 'Area from SAC : ',area
        print 'Area convex with DWL, CProfile',area_ck
        print 'convex area? {}'.format(area_ck>area)
        if area_ck/(dy*dx) < 10.:
            small_area = True
        else:
            small_area = False
        
        if a[2] is False or b[2] is False:
            print 'warning Could not locate bow fairing match point(s)'
            print 'start = {}, end = {}'.format(start,end)
            print 'o = ',o
            print 'a = ',a
            print 'b = ',b
        
        pt1 = a[1]
        pt2 = b[1]
        bottom  = linear_vertices((pt1[0],pt1[1],pt1[2]),
                                      (pt2[0],pt1[1],pt1[2]),4)
        side    = linear_vertices((pt2[0],pt1[1],pt2[2]),
                                      (pt2[0],pt2[1],pt2[2]),3)#4
        section = np.zeros((len(bottom)+len(side),3),float)
        section[:len(bottom)]=bottom
        section[len(bottom):]=side
        section = spline.Bspline(section[:,0:2],self.k,self.nump)
        ab = np.rad2deg(section.compute_tangent(0.))
        ae = np.rad2deg(section.compute_tangent(1.))
        #
        #**********************************************************************
        # special functions - circular bulb section of a curve
        # 
        #  NOT USED - this curve is further aft
        #  potentially used if bbow makes tunnel.
        #
        #        circle = spfuncs.circle
        #        ellipse = spfuncs.ellipse
        #        gfunc = spfuncs.gfunc
        #        #
        #        radius = 2.25 # assumed radius of bulbous bow
        #        rb_ = 1.8   # elliptic foci
        #        ra_ = 2.25  # elliptic foci
        #        #
        #        nump=30
        #        mp = .4 #
        #        c1 = circle(origin=[0.,2.25],radius=2.25)
        #        e1 = ellipse(origin=[0.,rb_],ra=ra_, rb=rb_)
        #        circ_p = np.linspace(-np.pi/2., 
        #                             np.pi/4., 
        #                             nump, 
        #                             endpoint=False)
        #        gfc3 = gfunc(copy.deepcopy(section), 
        #                     e1.circle_pt, 0., 
        #                     mp,circ_p, 10.)
        #
        #**********************************************************************
        # 
        inicurve = copy.deepcopy(section)
        #
        def stage1(section,ab,ae):
            #
            #******************************************************************
            # Begin Solve
            interval_data, small = interval_bounds(section)
            FPD = FormParameterDict(section)
            #
            #******************************************************************
            #  Integral Parameters
            FPD.add_Area_to_y_Constraint(kind='equality',
                                         value = area)
            #
            #******************************************************************
            # Tangents at the ends
            #            FPD.add_AngleConstraint(kind='min',
            #                                    location = 0.,
            #                                    value = ab)
            #            FPD.add_verticalAngleConstraint(kind='LS',
            #                                            location = 1.,
            #                                            value = 90.-ae)
            FPD.add_yFixity(index=1,
                            value=FPD.curve.vertices[0][0])
            FPD.add_yFixity(index=2,
                            value=FPD.curve.vertices[0][0])
            #FPD.add_xFixity(index=5,
            #                value=FPD.curve.vertices[-1][0])
            #FPD.add_xFixity(index=4,
            #                value=FPD.curve.vertices[-1][0])
            #
            #******************************************************************
            # Make sure the stern fairness curve does not penetrate 
            # the centerplane or waterplane
            top = section.vertices[-1][1]
            for param in [.1,.2,.3,.7,.8,.9]:
                FPD.add_xPointConstraint(kind='min',
                                         location=param,
                                         value=0.,
                                         weight=.5)
                FPD.add_yPointConstraint(kind='max',
                                         location=param,
                                         value=top,
                                         weight=.5)
            #
            #******************************************************************
            # Fairness
            if small_area:
                print 'solving bow transition curve'
                print '   with small cross section in mind'
                FPD.add_E1(kind='LS', weight = .1)
                FPD.add_E2(kind='LS', weight = .1)
                FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            else:
                FPD.add_E1(kind='LS', weight = .5)
                FPD.add_E2(kind='LS', weight = .5)
                #FPD.add_E3(kind='LS', weight = .5)
                FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            #
            #******************************************************************
            # 
            L = Lagrangian(FPD)
            #
            #******************************************************************
            # 
            interval_data, small = lagrangian_bounds(L, interval_data,
                                                     small, 1.e4)
            Lspline = IntervalLagrangeSpline(section, L, data = interval_data)
            #
            #******************************************************************
            # Solve
            print 'solving stage 1 bow transition curve'
            Lspline.optimize()
            #
            #******************************************************************
            # Done
            return Lspline
        
        
        
        def stagefix(section, ab,ae):
            #
            #******************************************************************
            # Begin Solve
            interval_data, small = interval_bounds(section)
            FPD = FormParameterDict(section)
            #
            #******************************************************************
            #  Integral Parameters
            FPD.add_Area_to_y_Constraint(kind='equality',
                                         value = area)
            #
            #******************************************************************
            # Tangents at the ends
            #            FPD.add_AngleConstraint(kind='LS',
            #                                    location = 0.,
            #                                    value = ab)
            #            FPD.add_verticalAngleConstraint(kind='LS',
            #                                            location = 1.,
            #                                            value = 90.-ae)
            FPD.add_yFixity(index=1,
                            value=FPD.curve.vertices[0][0])
            FPD.add_yFixity(index=2,
                            value=FPD.curve.vertices[0][0])
            #FPD.add_xFixity(index=5,
            #                value=FPD.curve.vertices[-1][0])
            #FPD.add_xFixity(index=4,
            #                value=FPD.curve.vertices[-1][0])
            #
            #******************************************************************
            # Fairness
            if small_area:
                print 'solving bow fairness curve'
                print '   with small cross section in mind'
                FPD.add_E1(kind='LS', weight = .1)
                FPD.add_E2(kind='LS', weight = .1)
                FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            else:
                FPD.add_E1(kind='LS', weight = .5)
                FPD.add_E2(kind='LS', weight = .5)
                FPD.add_E3(kind='LS', weight = .1)
                FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
            #
            #******************************************************************
            # 
            L = Lagrangian(FPD)
            #
            #******************************************************************
            # 
            interval_data, small = lagrangian_bounds(L, interval_data,
                                                     small, 1.e4)
            Lspline = IntervalLagrangeSpline(section, L, data = interval_data)
            #
            #******************************************************************
            # Solve
            print 'solving stage-fix bow transition curve'
            Lspline.optimize()
            #
            #******************************************************************
            # Done
            return Lspline
            
        
        Lspline = stage1(section,ab,ae)  
        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
            Lspline = stagefix(inicurve,ab,ae)  
        #
        # save
        #
        self.Lsplines.bowtransition = Lspline
        section = Lspline.curve
        #self.bowfairning.append(section)
        vts = frame.add1D(section.vertices,2,
                          where=pt1[2])
        section = spline.Bspline(vts, self.k, self.nump)
        self.bowtransition = section
        #
        #******************************************************************
        # tcurve map (hyperparameter control, 
        #                   possible graphical access to Lspline)
        self.tcurve_map[2] = {'Lspline':Lspline,
                               'tcurve':section}
        return

    def compute_stern_curve(self):
        Bmax = self.MaxBreadth/2.
        stern_depth = 0.5*self.MaxDraft
        bottom  = linear_vertices((0.,stern_depth,self.LengthDWL),(Bmax,stern_depth,self.LengthDWL),4)
        side    = linear_vertices((Bmax,stern_depth,self.LengthDWL),(Bmax,self.MaxDraft,self.LengthDWL),3)
        self.sternsec = np.zeros((len(bottom)+len(side),3),float)
        self.sternsec[:len(bottom)]=bottom
        self.sternsec[len(bottom):]=side
        self.sternsec = spline.Bspline(self.sternsec,self.k,self.nump)
        return

    def compute_stern_transverse(self):
        k=self.k
        nump=self.nump
        #xb = self.stations.Stern_Profile[0] #probably want to make a stern curve station thingy?
        #xe = self.stations.Stern_Profile[1]
        a =self.stations.Stern_Profile[2]  #fwd end of run, start of Flat SAC sections
        b =self.stations.Stern_Profile[3] #aft end of ship
        xb = self.DWL.CurvePoint(a)
        xe = 0.
        yb = 0.
        ye = self.MaxDraft
        return

    def bow_fairness_section_beam(self, frac=.5):
        """get the beam at the bow fairness control curve
        """
        self.bbfc = self.hullconstraints.state(
                self.hullconstraints.bbfc).getpoint(frac)
        if self.revise:
            self.hullconstraints.bbfc = ia(self.bbfc-self.small,
                                           self.bbfc+self.small)
            self.hullconstraints.AC_revise()
        return

    def bow_fairness_section_draft(self, frac=.5):
        """get the draft at the bow fairness control curve
        """
        self.dbfc = self.hullconstraints.state(
                self.hullconstraints.dbfc).getpoint(frac)
        if self.revise:
            self.hullconstraints.dbfc = ia(self.dbfc-self.small,
                                           self.dbfc+self.small)
            self.hullconstraints.AC_revise()
        return

    def bow_fairness_section_Coefficient(self, frac=.5):
        """get the Coefficient of area
        of the bow fairness control curve
        """
        self.Cbfc = self.hullconstraints.state(
                self.hullconstraints.Cbfc).getpoint(frac)
        if self.revise:
            self.hullconstraints.Cbfc = ia(self.Cbfc-self.small,
                                           self.Cbfc+self.small)
            self.hullconstraints.AC_revise()
        return

    def bow_fairness_section_area(self, frac=.5):
        """get the area of the bow fairness control beam
        """
        
        Abfc = self.hullconstraints.state(
                self.hullconstraints.Abfc).getpoint(frac)
        if self.revise:
            #        Cbfc = self.hullconstraints.state(
            #                self.hullconstraints.Cbfc).getpoint(frac)
            #        self.hullconstraints.Abfc = ia(Abfc-self.small,
            #                                       Abfc+self.small)
            self.Abfc = Abfc #* Cbfc
            self.hullconstraints.AC_revise()
        else:
            self.Abfc = Abfc 
        return

    def stern_fairness_section_beam(self, frac=.5):
        """get the beam at the stern fairness control curve
        """
        self.bsfc = self.hullconstraints.state(
                self.hullconstraints.bsfc).getpoint(frac)
        if self.revise:
            self.hullconstraints.bsfc = ia(self.bsfc-self.small,
                                           self.bsfc+self.small)
            self.hullconstraints.AC_revise()
        return

    def stern_fairness_section_draft(self, frac=.5):
        """get the draft at the stern fairness control curve
        """
        self.dsfc = self.hullconstraints.state(
                self.hullconstraints.dsfc).getpoint(frac)
        if self.revise:
            self.hullconstraints.dsfc = ia(self.dsfc-self.small,
                                           self.dsfc+self.small)
            self.hullconstraints.AC_revise()
        return

    def stern_fairness_section_Coefficient(self, frac=.5):
        """get the Coefficient of area
        of the stern fairness control curve
        """
        self.Csfc = self.hullconstraints.state(
                self.hullconstraints.Csfc).getpoint(frac)
        if self.revise:
            self.hullconstraints.Csfc = ia(self.Csfc-self.small,
                                           self.Csfc+self.small)
            self.hullconstraints.AC_revise()
        return

    def stern_fairness_section_area(self, frac=.5):
        """get the area of the stern fairness control beam
        """
        Asfc = self.hullconstraints.state(
                self.hullconstraints.Asfc).getpoint(frac)
        if self.revise:
            #        self.hullconstraints.Asfc = ia(Asfc-self.small,
            #                                       Asfc+self.small)
            #        Csfc = self.hullconstraints.state(
            #                self.hullconstraints.Csfc).getpoint(frac)
            self.Asfc = Asfc#*Csfc
            self.hullconstraints.AC_revise()
        else:
            self.Asfc = Asfc
        return

    #       from hull_from_simple_designspace import SAC
    def compute_SAC_new(self, LocMaxSection=.5, flat=False):
        sac = SAC(Amax = self.maxSecArea,
                  v1=self.SAC_entrance_area,
                  v2=self.SAC_mid_area,
                  v3=self.SAC_run_area,
                  l1 =self.SAC_entrance_len,
                  l2=self.flat_SAC,
                  l3=self.SAC_run_len,
                  Xe = self.SAC_fwd_Xc,
                  Xm = self.SAC_mid_Xc,
                  Xr = self.SAC_run_Xc,
                  #Abfc = self.Abfc,
                  #Abfc_loc = None,
                  #Asfc = self.Asfc,
                  #Asfc_loc = None,
                  stations = self.stations)
        self.SAC_class = sac
        sac.aggregate_SAC_curves_Small()
        self.SAC = sac.SAC
        self.sac_co = sac
        self.Lsplines.SAC_entrance=sac.entrance
        self.Lsplines.SAC_mid=sac.midbody
        self.Lsplines.SAC_run=sac.run
        return
    
#    def package_norms(self,Lspline):
#        Lspline.curve.pts_M_pts()
#        Lspline.curve.compute_arclength()
#        try:
#            E1 = copy.copy(Lspline.curve.E1.value)
#        except:
#            E1 = copy.copy(Lspline.curve.E1)
#        try:
#            E2 = copy.copy(Lspline.curve.E2.value)
#        except:
#            E2 = copy.copy(Lspline.curve.E2)
#        try:
#            E3 = copy.copy(Lspline.curve.E3.value)
#        except:
#            E3 = copy.copy(Lspline.curve.E3)
#        try:
#            S = copy.copy(Lspline.curve.AL.value)
#        except:
#            S = copy.copy(Lspline.curve.AL)
#        norms = [E1,E2,E3,S]
#        return norms
        

    def compute_SAC(self, LocMaxSection=.5,
                    Xc=None, flat=True,
                    Ab=None,Ae=None):
        """
            TODO : use bow_fairness_section properties
            to fix the bow fairness curve feasability
            
            TODO:
                use SD.hull.stations.Bow_Fairing
                together with SD.hull.A_mid
                to enforce SAC 
                to have enough area to design in the 
                bulbous bow
                
                
                
                
            SD.hull.stations.Bow_Fairing
                Out[6]: [19.224785534211648, 
                           19.224785534211648, 
                           0.15266043165628, 
                           0.15266043165628]
                
                
                Naive SAC computation:
                SD.hull.SAC(.15266043165628)
                
                Out[8]: array([ 14.69608985, 227.49719917])
                
                
                
            SAC required by bbow:
                
                SD.hull.A_mid
                Out[7]: 11.298465400737983
                
                    
                not this, we are using half this per side:
                SD.hull.A_mid*2.
                Out[37]: 22.596930801475967
                
                
                
                
            Actual availible SAC:
                
                start,end,ao,bo = self.stations.Bow_Fairing
                o,a,b = self.get_transverse_start_end_param_loc(start,end)
                
                area = self.SAC.CurvePoint(o[0])[1]*0.5
                
                self.SAC.CurvePoint(o[0])[1]*0.5
                Out[25]: 132.00296824600264
                
                
            1/2 Volume of a Spherical Cap:
                
                V = (1/2)(pi/3)(h**2)(3r-h)
                
        """
        print '\nSAC'
        k=self.k
        nump=self.nump
        xb = 0.
        xe = self.LengthDWL
        #
        if Ab is None:
            yb = 0.
        else:
            yb = Ab
        if Ae is None:
            ye = 0.
        else:
            ye = Ae
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
            #            FPD.add_AngleConstraint(kind='LS',
            #                                    location = 0.,
            #                                    value = alphab)
            #            FPD.add_AngleConstraint(kind='LS',
            #                                    location = 1.,
            #                                    value = alphae)
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
            # SAC  Get scalars
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
                FPD.add_E1(kind='LS', weight = .5/E1_0)
                FPD.add_E2(kind='LS', weight = .5/E2_0)
                FPD.add_E3(kind='LS', weight = .5/E3_0)
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
            # SAC
            # NOW Entering Solver Stage 2
            norms = self.package_norms(Lspline)
            E1_0 = norms[0]
            E2_0 = norms[1]
            E3_0 = norms[2]
            S_0  = norms[3]
            curve = Lspline.curve
            #
            #**********************************************************************
            # 
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve)
            #
            #******************************************************************
            # Integral Parameters
            FPD.add_AreaConstraint(kind='equality',
                                   value = self.BareHullVolume)
            FPD.add_XcConstraint(kind='equality', 
                                 value=self.LCG)
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
            #            FPD.add_xVertexConstraint(kind='max',
            #                                     index = 3,
            #                                     value = self.stations.FOSAC[0])
            #            FPD.add_xVertexConstraint(kind='equality',
            #                                     index = 4,
            #                                     value = self.stations.FOSAC[0])
            #            FPD.add_xVertexConstraint(kind='min',
            #                                     index = 5,
            #                                     value = self.stations.FOSAC[0])
            #            FPD.add_xVertexConstraint(kind='max',
            #                                     index = 5,
            #                                     value = self.stations.FOSAC[1])
            #            FPD.add_xVertexConstraint(kind='equality',
            #                                     index = 6,
            #                                     value = self.stations.FOSAC[1])
            #            FPD.add_xVertexConstraint(kind='min',
            #                                     index = 7,
            #                                     value = self.stations.FOSAC[1])
            #
            #print 'using x - fixity'
            FPD.add_xFixity(index=5,
                            value=self.stations.FOSAC[0])
            FPD.add_xFixity(index=6,
                            value=self.stations.FOSAC[1])
            #
            #******************************************************************
            # Y FLAT SAC (3,4,5,6,7)
            #            FPD.add_yVertexConstraint(kind = 'equality', 
            #                                      value=self.maxSecArea, 
            #                                      index =3)
            #            FPD.add_yVertexConstraint(kind = 'equality', 
            #                                      value=self.maxSecArea, 
            #                                      index =4)
            #            FPD.add_yVertexConstraint(kind = 'equality', 
            #                                      value=self.maxSecArea, 
            #                                      index =5)
            #            FPD.add_yVertexConstraint(kind = 'equality', 
            #                                      value=self.maxSecArea, 
            #                                      index =6)
            #            FPD.add_yVertexConstraint(kind = 'equality', 
            #                                      value=self.maxSecArea, 
            #                                      index =7)
            #-----------------------?
            #print 'using y - fixity'
            #FPD.add_yFixity(index=3,
            #                value=self.maxSecArea)
            FPD.add_yFixity(index=4,
                            value=self.maxSecArea)
            FPD.add_yFixity(index=5,
                            value=self.maxSecArea)
            FPD.add_yFixity(index=6,
                            value=self.maxSecArea)
            FPD.add_yFixity(index=7,
                            value=self.maxSecArea)
            #-----------------------?
            #FPD.add_yFixity(index=8,
            #                value=self.maxSecArea)
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
            FPD.add_E1(kind='LS', weight = 1.5/E1_0)
            FPD.add_E2(kind='LS', weight = .5/E2_0)
            FPD.add_E3(kind='LS', weight = .5/E3_0)
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
            """SAC fix stage 1
            """
            norms = self.package_norms(Lspline)
            E1_0 = norms[0]
            E2_0 = norms[1]
            E3_0 = norms[2]
            S_0  = norms[3]
            curve = Lspline.curve
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
            FPD.add_E1(kind='LS', weight = 1.5/E1_0)
            FPD.add_E2(kind='LS', weight = .5/E2_0)
            FPD.add_E3(kind='LS', weight = .5/E3_0)
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
            # scaling values
            norms = self.package_norms(Lspline)
            E1_0 = norms[0]
            E2_0 = norms[1]
            E3_0 = norms[2]
            S_0  = norms[3]
            curve = Lspline.curve
            #
            #**********************************************************************
            # Fixit stage 2: flat
            # initialize
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
            # x flat
            s1 = self.stations.FOSAC[0]
            s2 = self.stations.FOSAC[1]
            FPD.add_xVertexConstraint(kind = 'equality', 
                                      value=s1, 
                                      index =4)
            FPD.add_xVertexConstraint(kind = 'equality', 
                                      value=s2, 
                                      index =6)
            ##-------------------------------------- well ordered vertices
            FPD.add_xVertexConstraint(kind = 'max', 
                                      value=s1, 
                                      index =3)
            FPD.add_xVertexConstraint(kind = 'min', 
                                      value=s1, 
                                      index =5)
            FPD.add_xVertexConstraint(kind = 'max', 
                                      value=s2, 
                                      index =5)
            FPD.add_xVertexConstraint(kind = 'min', 
                                      value=s1, 
                                      index =7)
            #
            #******************************************************************
            # Fairness
            FPD.add_E1(kind='LS', weight = 1.5/E1_0)
            FPD.add_E2(kind='LS', weight = .5/E2_0)
            FPD.add_E3(kind='LS', weight = .5/E3_0)
            FPD.add_ArcLengthApprox(kind='LS', weight = 10./S_0)
            #
            #******************************************************************
            # 
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            #
            #******************************************************************
            #OPTIMIZE - 2nd stage fix
            print 'SAC fixit stage2 solve:'
            Lspline.optimize(stop = 50)
            return Lspline
            
        self.savebadlsplines = []
        
        inicurve = copy.copy(curve)
        Lspline,FPD = stage1(curve)
        self.savebadlsplines.append(Lspline)
        
        inicurve2 = copy.copy(Lspline.curve)
        Lspline,FPD = stage2(Lspline,curve)
        self.savebadlsplines.append(Lspline)
        
        #check solution:
        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
            Lspline = stagefix1(inicurve,Lspline)
            self.savebadlsplines.append(Lspline)
            Lspline = stagefix2(Lspline)
            self.savebadlsplines.append(Lspline)
        #elegant solution
        #
        # save solution:
        self.Lsplines.SAC=Lspline
        self.SAC = Lspline.curve
        verts = self.SAC.vertices
        vnew = frame.xy_to_zyx_elevate(verts,elevation=0.)
        self.SAC_3dvis = spline.Bspline(vnew, k, nump)
        
        
        
        #
        #******************************************************************
        # Get the total sectional area at the bow fairness curve:
        """
        
            Actual availible SAC:
                
                start,end,ao,bo = self.stations.Bow_Fairing
                o,a,b = self.get_transverse_start_end_param_loc(start,end)
                
                area = self.SAC.CurvePoint(o[0])[1]*0.5
                
                self.SAC.CurvePoint(o[0])[1]*0.5
                Out[25]: 132.00296824600264
                
                
            1/2 Volume of a Spherical Cap:
                
                V = (1/2)(pi/3)(h**2)(3r-h)
        """
        start,end,ao,bo = self.stations.Bow_Fairing
        sacwhich = 2
        if self.SAC.dim ==2:
            sacwhich =0
        #do it this way as the class function for this is not availible yet.
        o = self.SAC.FindPoint(start,sacwhich) 
        #
        xloc_bow_fair = self.SAC.CurvePoint(o[0])[0]
        Abfc = self.SAC.CurvePoint(o[0])[1]*0.5
        
        
        if self.bbow_from_SAC:
        
            self.A_mid = self.hyperparameter_input.bulb_volume_fraction*Abfc
            
            
            self.BulbLength =xloc_bow_fair*\
                                    self.hyperparameter_input.bulb_length_fraction
                                    
            self.Bulb_SAC_intersect = self.SAC.FindPoint(self.A_mid, 1)
            bbow_SAC_fwd, junk = self.SAC.CurveSplit(
                                    self.Bulb_SAC_intersect[0]) #fwd SAC becomes 
            
            transform = 1.5
            newvertices = []
            for vtx in bbow_SAC_fwd.vertices:
                #frac = vtx[0]/xloc_bow_fair
                newvertices.append([vtx[0]*transform,
                                    vtx[1]])
            bbow_vertices = newvertices#list(bbow_SAC_fwd.vertices)
            
            bbow_vertices.append(bbow_SAC_fwd.vertices[-1])
            bbow_vertices.append(bbow_SAC_fwd.vertices[-1])
            bbow_vertices.append(bbow_SAC_fwd.vertices[-1])
            bbow_vertices.append(np.asarray(    [xloc_bow_fair,
                                                 bbow_SAC_fwd.vertices[-1][1]])
                                )
            self.Bulb_SAC = spline.Bspline(bbow_vertices,
                                           k=self.k, nump=self.nump)
        
        else: #info flows from bbow to SAC:
            bbow_vertices = [np.asarray([0.,0.])]
            bbow_vertices.append( self.SAC.vertices[0] )
            bbow_vertices.append( self.SAC.vertices[0] )
            bbow_vertices.append( self.SAC.vertices[0] )
            
            bbow_Xe  = self.BulbLength
            bbow_Ye =  self.SAC.vertices[0,1]
            
            bbow_vertices.append( np.asarray([bbow_Xe,bbow_Ye]) )
            bbow_vertices.append( np.asarray([bbow_Xe,bbow_Ye]) )
            bbow_vertices.append( np.asarray([bbow_Xe,bbow_Ye]) )
            
            
            
            bbow_vertices.append( np.asarray([
                                xloc_bow_fair,
                                Abfc*self.hyperparameter_input.bulb_volume_fraction]) )
            
            #interpolatedBspline
            self.Bulb_SAC = spline.Bspline(
                                            bbow_vertices,
                                            k=self.k, 
                                            nump=self.nump)
            
        
            
        self.results_dict['vol'] = self.SAC.area
        self.results_dict['lfsac'] = self.stations.FOSAC[1]-self.stations.FOSAC[0]
        self.results_dict['lfsac_data'] = 'min flat length'
        
        
        self.results_dict['LCG'] = self.SAC.Xc
        return


    def compute_DWL(self, LocMaxSection=.5, Xc=None,
                    flat=True, height=None, alphab=0.,alphae=0.):
        """
        
        Xc = None
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

        curve.vertices = v
        
        
        design notes
        ----------
        *pull LCG of this one forward
        (and LCG of cLProfile aft?)
        *what about SAC vol. distribution
        I really don't know what to do to make good shapes more likely
        
        dev
        ----------
        import curve as spline
        from ADILS import interval_bounds
        from FormParameter import FormParameterDict
        from ADILS import Lagrangian
        from ADILS import lagrangian_bounds
        
        from ADILS import IntervalLagrangeSpline
        
        from   adials_gui3d      import frame
        
        elevation = self.MaxDraft
        Xc = None
        
        
        """
        print '\nDWL'
        print '\n hull_from_simple_designspace: compute_DWL (Design Waterline)'
        if height is None:
            elevation = self.MaxDraft #location of the DWL above the keel
        else:
            elevation = height
        k       = self.k
        nump    = self.nump
        
        Bmax = self.MaxBreadth/2.
        area = self.WLarea/2.
        
        xb      = 0.
        yb      = 0.
        xe      = self.LengthDWL
        if self.hullformtype == 'osv':
            ye  = Bmax
        else:
            ye  = 0.
        #alphab = 0.#-5.
        #alphae = 0.#-15.
        Cab_given = 0.
        Cae_given = 0.
        slope = 'down'


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
            # DWL
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
            #FPD.add_AreaConstraint(kind='LS',
            #                       value = area)
            #FPD.add_AreaConstraint(kind='min',
            #                       value = area-.25*area)
            #FPD.add_AreaConstraint(kind='max',
            #                       value = area+.25*area)
            FPD.add_AreaConstraint(kind='equality',
                                   value = area)
            #
            #**********************************************************************
            #   DWL/2 might be a good starting parameter for a work deck (not DWL)
            if Xc is None:
                FPD.add_XcConstraint(kind='LS', 
                                     value=self.LCG)
            else:
                FPD.add_XcConstraint(kind='LS', 
                                     value=Xc)
            #
            #******************************************************************
            # nose intercept should be purely transverse 
            #- no longitudinal slope
            #--------------------------------------
            # beyond baseline DWL
            #--------------------------------------
            # +
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 1,
                                     value = 0.)
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 2,
                                     value = 0.)
            #--------------------------------------
            # +++
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 3,
                                     value = 0.)
            #
            #--------------------------------------
            # ++
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 8,
                                     value = xe)
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 9,
                                     value = xe)
            #--------------------------------------
            #
            #FPD.add_yVertexConstraint(kind='LS',
            #                         index = 1,
            #                         value = Bmax/4.)
            #FPD.add_xFixity(index=1,
            #                value=0.)
            #
            #******************************************************************
            # x-flat DWL  (4,5,6 + optional:  keep 3,7 out)
#            FPD.add_xVertexConstraint(kind='max',
#                                     index = 3,
#                                     value = self.stations.FOWL[0])
#            FPD.add_xVertexConstraint(kind='max',
#                                     index = 4,
#                                     value = self.stations.FOWL[0])
            #FPD.add_xVertexConstraint(kind='min',
            #                         index = 5,
            #                         value = self.stations.FOWL[0])
            #FPD.add_xVertexConstraint(kind='min',
            #                         index = 6,
            #                         value = self.stations.FOWL[1])
#            FPD.add_xVertexConstraint(kind='min',
#                                     index = 7,
#                                     value = self.stations.FOWL[1])
            ##**************************************************************
            # DWL
            #print 'Adding x - fixity!'
            FPD.add_xFixity(index=5,
                            value=self.stations.FOWL[0])
            FPD.add_xFixity(index=6,
                            value=self.stations.FOWL[1])
            #
            #******************************************************************
            # y-flat DWL (3,4,5,6,7)
            # ensure curvature flatness at vertex 4 via vertex 3 (and 5,6,7)
            
            #FPD.add_CurvatureConstraint(kind='equality',
            #                            location = 0.25,
            #                            value = 0.)
            #FPD.add_CurvatureConstraint(kind='equality',
            #                            location = 0.75,
            #                            value = 0.)
            #FPD.add_xPointConstraint(kind='equality',
            #                        location = 0.25,
            #                        value = self.stations.FOWL[0])
            #FPD.add_xPointConstraint(kind='equality',
            #                        location = 0.75,
            #                        value = self.stations.FOWL[1])
            
            # curvature continuity----
#            FPD.add_yVertexConstraint(kind='max',
#                                     index = 3,
#                                     value = Bmax)
            #FPD.add_yVertexConstraint(kind='LS',
            #                         index = 3,
            #                         value = Bmax)
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 3,
            #                         value = Bmax+self.tol*.25) 
            
            # Now the mid section 4,5,6
            #FPD.add_yVertexConstraint(kind='equality',
            #                         index = 4,
            #                         value = Bmax)
            #FPD.add_yVertexConstraint(kind='equality',
            #                         index = 5,
            #                         value = Bmax)
            ##
            #FPD.add_yVertexConstraint(kind='equality',
            #                         index = 6,
            #                         value = Bmax)    
            
            #curvature continuity-------------
            #FPD.add_yVertexConstraint(kind='LS',
            #                         index = 7,
            #                         value = Bmax)  
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 7,
            #                         value = Bmax+self.tol*.25) 
            #
            #**********************************************************************
            #  y-flat DWL (use fixity)
            #-----------------------?
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 1,
            #                         value = Bmax)  
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 2,
            #                         value = Bmax)  
            FPD.add_yFixity(index=2,
                            value=Bmax)
            #-----------------------?
            #print 'Adding y - fixity!'
            FPD.add_yFixity(index=3,
                            value=Bmax)
            FPD.add_yFixity(index=4,
                            value=Bmax)
            FPD.add_yFixity(index=5,
                            value=Bmax)
            FPD.add_yFixity(index=6,
                            value=Bmax)
            FPD.add_yFixity(index=7,
                            value=Bmax)
            #-----------------------?
            FPD.add_yFixity(index=8,
                            value=Bmax)
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 8,
            #                         value = Bmax)  
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 9,
            #                         value = Bmax)  
            #-----------------------?
            #
            #******************************************************************
            #    
            #            FPD.add_yVertexConstraint(kind='min',
            #                                     index = 1,
            #                                     value = 0.)
            #            FPD.add_yVertexConstraint(kind='min',
            #                                     index = 2,
            #                                     value = 0.)
            #            FPD.add_yVertexConstraint(kind='min',
            #                                     index = 8,
            #                                     value = 0.)
            #            FPD.add_yVertexConstraint(kind='min',
            #                                     index = 9,
            #                                     value = 0.)
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
            Lspline.optimize(stop=25)
            #Lspline.curve.knot_insertion(.5) #get to 11
            self.Lsplines.DWL=Lspline
            self.DWL2D = Lspline.curve
            
            return Lspline
        
        
        def Transom_DWL(curve):
            """Lspline solver 
            
                for the 
                
                DWL curve
            
            (intended for a transom sterned hull)
            """
            # DWL
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
                                     value=self.LCG*1.1) #DWL square aft
            else:
                FPD.add_XcConstraint(kind='LS', 
                                     value=Xc)
            #
            #******************************************************************
            # nose intercept should be purely transverse 
            #- no longitudinal slope
            #--------------------------------------
            # beyond baseline DWL
            #--------------------------------------
            # +
            FPD.add_xFixity(index=1,
                            value=xb)
            #
            # curvature continuity is "too flat nosed?"
            #FPD.add_xFixity(index=2,
            #                value=xb)
            #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
            # min y so the nose can't go negative
            FPD.add_yVertexConstraint(kind='min',
                                     index = 1,
                                     value = 0.)
            # pull the nose so it tends not to dip past CL
            FPD.add_yVertexConstraint(kind='LS',
                                     index = 1,
                                     value = Bmax)
            FPD.add_yVertexConstraint(kind='LS',
                                     index = 2,
                                     value = Bmax)
            FPD.add_yVertexConstraint(kind='LS',
                                     index = 3,
                                     value = Bmax)
            #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
            #--------------------------------------
            # +++
            FPD.add_xVertexConstraint(kind='min',
                                     index = 3,
                                     value = 0.)
            #
            #--------------------------------------
            # ++
            FPD.add_xVertexConstraint(kind='max',
                                     index = 8,
                                     value = xe)
            ##**************************************************************
            # DWL
            #print 'Adding x - fixity!'
            FPD.add_xFixity(index=5,
                            value=self.stations.FOWL[0])
            FPD.add_xFixity(index=6,
                            value=self.stations.FOWL[1])
            
            #
            #**********************************************************************
            #  Or Square off the DWL...
            # taking this from ++
            FPD.add_xFixity(index=9,
                            value=xe)
            #no longer necessary, since we will make
            # this "really a transom ending' DWL curve:
            #FPD.add_yFixity(index=9,
            #                value=Bmax)
            #
            #
            #**********************************************************************
            #  y-flat DWL (use fixity)
            #-----------------------?
            #FPD.add_yFixity(index=2,
            #                value=Bmax)
            #-----------------------?
            #print 'Adding y - fixity!'
            #FPD.add_yFixity(index=3,
            #                value=Bmax)
            FPD.add_yFixity(index=4,
                            value=Bmax)
            FPD.add_yFixity(index=5,
                            value=Bmax)
            FPD.add_yFixity(index=6,
                            value=Bmax)
            FPD.add_yFixity(index=7,
                            value=Bmax)
            #-----------------------?
            FPD.add_yFixity(index=8,
                            value=Bmax)
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 8,
            #                         value = Bmax)  
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 9,
            #                         value = Bmax)  
            #-----------------------?
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
            print 'Doing first stage Transom DWL solve'
            Lspline.optimize(stop=25)
            self.Lsplines.DWL=Lspline
            self.DWL2D = Lspline.curve
            
            return Lspline
        
        
        
        def Transom_DWL_fixit(curve):
            # DWL
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
                FPD.add_XcConstraint(kind='LS', 
                                     value=Xc)
            #
            #******************************************************************
            # nose intercept should be purely transverse 
            #- no longitudinal slope
            #--------------------------------------
            # beyond baseline DWL
            #--------------------------------------
            # +
            FPD.add_xFixity(index=1,
                            value=xb)
            #--------------------------------------
            # +++
            FPD.add_xVertexConstraint(kind='min',
                                     index = 3,
                                     value = 0.)
            #
            #--------------------------------------
            # ++
            FPD.add_xVertexConstraint(kind='max',
                                     index = 8,
                                     value = xe)
            ##**************************************************************
            # DWL
            #print 'Adding x - fixity!'
            FPD.add_xFixity(index=5,
                            value=self.stations.FOWL[0])
            FPD.add_xFixity(index=6,
                            value=self.stations.FOWL[1])
            
            #
            #**********************************************************************
            #  Or Square off the DWL...
            # taking this from ++
            FPD.add_xFixity(index=9,
                            value=xe)
            #no longer necessary, since we will make
            # this "really a transom ending' DWL curve:
            #FPD.add_yFixity(index=9,
            #                value=Bmax)
            #
            #
            #**********************************************************************
            #  y-flat DWL (use fixity)
            #-----------------------?
            #FPD.add_yFixity(index=2,
            #                value=Bmax)
            #-----------------------?
            #print 'Adding y - fixity!'
            #FPD.add_yFixity(index=3,
            #                value=Bmax)
            FPD.add_yFixity(index=4,
                            value=Bmax)
            FPD.add_yFixity(index=5,
                            value=Bmax)
            FPD.add_yFixity(index=6,
                            value=Bmax)
            FPD.add_yFixity(index=7,
                            value=Bmax)
            #-----------------------?
            FPD.add_yFixity(index=8,
                            value=Bmax)
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 8,
            #                         value = Bmax)  
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 9,
            #                         value = Bmax)  
            #-----------------------?
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
            print 'Doing Fixit stage Transom DWL solve'
            Lspline.optimize(stop=25)
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
            # DWL initialize
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve)
            #
            #*********************************************************************
            # DWL area
            FPD.add_AreaConstraint(kind='LS',
                                   value = area,
                                   weight=1.5/E1_0)
            #FPD.add_AreaConstraint(kind='min',
            #                       value = area-.25*area)
            #FPD.add_AreaConstraint(kind='max',
            #                       value = area+.25*area)
            #FPD.add_AreaConstraint(kind='equality',
            #                       value = area)
            #
            #*********************************************************************
            # X-centroid
            #   DWL/2 might be a good parameter for a work deck (not DWL)
            if Xc is None:
                FPD.add_XcConstraint(kind='LS', 
                                     value=self.LCG,
                                     weight=.5/E1_0)
            else:
                print 'DWL setting Xc constraint to ',Xc
                FPD.add_XcConstraint(kind='LS', 
                                     value=Xc,
                                     weight=.5/E1_0)
            #
            #******************************************************************
            # DWL NOSE CONSTRAINTS
            # nose intercept should be purely transverse 
            #- no longitudinal slope
            # 
            #            FPD.add_relative_xVertexConstraint(
            #                                     kind='equality',
            #                                     location=None,
            #                                     index = 0,
            #                                     index2 = 1,
            #                                     value = 0.,
            #                                     weight=1.)
            #print 'Adding x - nose fixity!'
            FPD.add_xFixity(index=1,
                            value=0.)
            #FPD.add_xVertexConstraint(kind='LS',
            #                         index = 1,
            #                         value = 0.,
            #                         weight=.5/E1_0)
            
            # end of nose tangent constraint 
            #**********************************************************************
            # DWL    - stay within parameters
            # ensure no 0 crossing
            #FPD.add_xVertexConstraint(kind='min',
            #                         index = 1,
            #                         value = 0.)
            FPD.add_xVertexConstraint(kind='min',
                                     index = 2,
                                     value = 0.)
            FPD.add_xVertexConstraint(kind='max',
                                     index = 8,
                                     value = 1.1*self.LengthDWL)
            FPD.add_xVertexConstraint(kind='max',
                                     index = 9,
                                     value = 1.1*self.LengthDWL)
            #FPD.add_xFixity(index=9,
            #                value=self.LengthDWL)
            #
            #******************************************************************
            # DWL stage2
            # end condition stays nice:
            #            FPD.add_yVertexConstraint(kind='max',
            #                                     index = 1,
            #                                     value = Bmax*1.3)
            # pull up on y
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 1,
            #                                     value = Bmax,
            #                                     weight=.5/E1_0)
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 2,
            #                                     value = Bmax)
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 8,
            #                                     value = Bmax)
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 9,
            #                                     value = Bmax)
            # pull in on x
            #            FPD.add_xVertexConstraint(kind='LS',
            #                                     index = 1,
            #                                     value = self.stations.FOWL[0])
            #            FPD.add_xVertexConstraint(kind='LS',
            #                                     index = 2,
            #                                     value = self.stations.FOWL[0])
            #            FPD.add_xVertexConstraint(kind='LS',
            #                                     index = 8,
            #                                     value = self.stations.FOWL[0])
            #            FPD.add_xVertexConstraint(kind='LS',
            #                                     index = 9,
            #                                     value = self.stations.FOWL[0])
            
            #  ensure no 0 crossing
            #FPD.add_yVertexConstraint(kind='min',
            #                         index = 1,
            #                         value = 0.)
            FPD.add_yVertexConstraint(kind='min',
                                     index = 2,
                                     value = 0.)
            FPD.add_yVertexConstraint(kind='min',
                                     index = 8,
                                     value = 0.)
            FPD.add_yVertexConstraint(kind='min',
                                     index = 9,
                                     value = 0.)
            
            # ensure no wild runoff
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 1,
            #                         value = 1.05*Bmax)
            #            FPD.add_yVertexConstraint(kind='max',
            #                                     index = 2,
            #                                     value = 1.05*Bmax)
            #            FPD.add_yVertexConstraint(kind='max',
            #                                     index = 8,
            #                                     value = 1.05*Bmax)
            #            FPD.add_yVertexConstraint(kind='max',
            #                                     index = 9,
            #                                     value = 1.05*Bmax)
            #
            #**********************************************************************
            # DWL 
            # x-flat (4,5,6 + 3,7 )
            # stage 2solve
#            FPD.add_xVertexConstraint(kind='max',
#                                     index = 3,
#                                     value = self.stations.FOWL[0])
            FPD.add_xVertexConstraint(kind='max',
                                     index = 4,
                                     value = self.stations.FOWL[0])
#            #FPD.add_xVertexConstraint(kind='min',
#            #                         index = 5,
#            #                         value = self.stations.FOWL[0])
#            #FPD.add_xVertexConstraint(kind='max',
#            #                         index = 5,
#            #                         value = self.stations.FOWL[1])
#            #FPD.add_xVertexConstraint(kind='equality',
#            #                         index = 6,
#            #                         value = self.stations.FOWL[1])
            FPD.add_xVertexConstraint(kind='min',
                                     index = 7,
                                     value = self.stations.FOWL[1])
            #
            #
            #print 'Adding x - fixity!'
            FPD.add_xFixity(index=5,
                            value=self.stations.FOWL[0])
            FPD.add_xFixity(index=6,
                            value=self.stations.FOWL[1])
            #
            #**********************************************************************
            # DWL
            # y-flat (3,4,5,6,7) stage 2
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 3,
            #                                     value = Bmax)  
            #            FPD.add_yVertexConstraint(kind='max',
            #                                     index = 3,
            #                                     value = Bmax+self.tol*.25)  
            #            FPD.add_yVertexConstraint(kind='equality',
            #                                     index = 4,
            #                                     value = Bmax)
            #            FPD.add_yVertexConstraint(kind='equality',
            #                                     index = 5,
            #                                     value = Bmax)
            #            FPD.add_yVertexConstraint(kind='equality',
            #                                     index = 6,
            #                                     value = Bmax)  
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 7,
            #                                     value = Bmax)    
            #            FPD.add_yVertexConstraint(kind='max',
            #                                     index = 7,
            #                                     value = Bmax+self.tol*.25)  
            #print 'Adding y - fixity!'
            FPD.add_yFixity(index=3,
                            value=Bmax)
            FPD.add_yFixity(index=4,
                            value=Bmax)
            FPD.add_yFixity(index=5,
                            value=Bmax)
            FPD.add_yFixity(index=6,
                            value=Bmax)
            FPD.add_yFixity(index=7,
                            value=Bmax)
            #FPD.add_yFixity(index=8,
            #                value=Bmax)
            #
            #**********************************************************************
            # DWL fairness norms
            FPD.add_E1(kind='LS', weight = .5/E1_0)
            FPD.add_E2(kind='LS', weight = .5/E2_0)
            FPD.add_E3(kind='LS', weight = .5/E3_0)
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
            Lspline.optimize(stop = 25)
            return Lspline
        
        
        def fix1():
            #
            #**********************************************************************
            # continuation passing style curve fixer, DWL
            def doit(Lspline=None,curve=None):
                """DWL fix
                """
                #
                #**********************************************************************
                # get weights
                if Lspline is not None:
                    norms = self.package_norms(Lspline)
                    E1_0 = norms[0]
                    E2_0 = norms[1]
                    E3_0 = norms[2]
                    S_0  = norms[3]
                    curve = Lspline.curve
                else:
                    E1_0 = 1.
                    E2_0 = 1.
                    E3_0 = 1.
                    S_0  = 1.
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
                #if Xc is None:
                #    FPD.add_XcConstraint(kind='LS', 
                #                         value=self.LCG)
                #else:
                #    print 'DWL setting Xc constraint to ',Xc
                #    FPD.add_XcConstraint(kind='LS', value=Xc)
                #
                #**********************************************************************
                # x-flat DWL  (4,5,6 + optional:  keep 3,7 out)
                # stage fix solve
                
                FPD.add_xVertexConstraint(kind='max',
                                         index = 3,
                                         value = self.stations.FOWL[0])
                FPD.add_xVertexConstraint(kind='max',
                                         index = 4,
                                         value = self.stations.FOWL[0])
                FPD.add_xVertexConstraint(kind='min',
                                         index = 7,
                                         value = self.stations.FOWL[1])
                #**************************************************************
                #print 'Adding x - fixity!'
                FPD.add_xFixity(index=5,
                                value=self.stations.FOWL[0])
                FPD.add_xFixity(index=6,
                                value=self.stations.FOWL[1])
                #
                #**********************************************************************
                # y-flat DWL (4,5,6) 
                # stage fix solve
                #-----------------------?
                FPD.add_yVertexConstraint(kind='max',
                                         index = 1,
                                         value = Bmax)  
                FPD.add_yVertexConstraint(kind='max',
                                         index = 2,
                                         value = Bmax)  
                #print 'Adding y - fixity!'
                FPD.add_yFixity(index=3,
                                value=Bmax)
                FPD.add_yFixity(index=4,
                                value=Bmax)
                FPD.add_yFixity(index=5,
                                value=Bmax)
                FPD.add_yFixity(index=6,
                                value=Bmax)
                FPD.add_yFixity(index=7,
                                value=Bmax)
                #-----------------------?
                #FPD.add_yFixity(index=8,
                #                value=Bmax)
                FPD.add_yVertexConstraint(kind='max',
                                         index = 8,
                                         value = Bmax)  
                FPD.add_yVertexConstraint(kind='max',
                                         index = 9,
                                         value = Bmax)  
                #
                #
                #**********************************************************************
                # DWL fairness norms
                FPD.add_E1(kind='LS', weight = 1./E1_0)
                FPD.add_E2(kind='LS', weight = .5/E2_0)
                FPD.add_E3(kind='LS', weight = .01/E3_0)
                FPD.add_ArcLengthApprox(kind='LS', weight = 1.5/S_0)
                #
                #**********************************************************************
                # setup
                L = Lagrangian(FPD)
                interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
                Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
                #
                #**********************************************************************
                # solve
                print 'Doing DWL Repair solve (old fix)'
                Lspline.optimize(stop = 25)
                return Lspline
            return doit
        
        
        
        
        def improve_stage1(Lspline):
            """DWL improve normalization of stage 1
            use only if stage 1 succeeded
            """
            norms = self.package_norms(Lspline)
            E1_0 = norms[0]
            E2_0 = norms[1]
            E3_0 = norms[2]
            S_0  = norms[3]
            curve = Lspline.curve
            # DWL
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
            #FPD.add_AreaConstraint(kind='LS',
            #                       value = area)
            #FPD.add_AreaConstraint(kind='min',
            #                       value = area-.25*area)
            #FPD.add_AreaConstraint(kind='max',
            #                       value = area+.25*area)
            FPD.add_AreaConstraint(kind='equality',
                                   value = area)
            #
            #**********************************************************************
            #   DWL/2 might be a good starting parameter for a work deck (not DWL)
            if Xc is None:
                FPD.add_XcConstraint(kind='LS', 
                                     value=self.LCG, 
                                     weight = 1./E1_0)
            else:
                FPD.add_XcConstraint(kind='LS', 
                                     value=Xc, 
                                     weight = 1./E1_0)
            #
            #******************************************************************
            # nose intercept should be purely transverse 
            #- no longitudinal slope
            #--------------------------------------
            # beyond baseline DWL
            #--------------------------------------
            # +
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 1,
                                     value = 0., 
                                     weight = 1./E1_0)
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 2,
                                     value = 0., 
                                     weight = 1./E1_0)
            #--------------------------------------
            # +++
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 3,
                                     value = 0., 
                                     weight = 1./E1_0)
            #
            #--------------------------------------
            # ++
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 8,
                                     value = xe, 
                                     weight = 1./E1_0)
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 9,
                                     value = xe, 
                                     weight = 1./E1_0)
            #--------------------------------------
            #
            ##**************************************************************
            # DWL
            #print 'Adding x - fixity!'
            FPD.add_xFixity(index=5,
                            value=self.stations.FOWL[0])
            FPD.add_xFixity(index=6,
                            value=self.stations.FOWL[1])
            #
            #******************************************************************
            # y-flat DWL (3,4,5,6,7)
            # ensure curvature flatness at vertex 4 via vertex 3 (and 5,6,7)
            #curvature continuity-------------
            #FPD.add_yVertexConstraint(kind='LS',
            #                         index = 7,
            #                         value = Bmax)  
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 7,
            #                         value = Bmax+self.tol*.25) 
            #
            #**********************************************************************
            #  y-flat DWL (use fixity)
            #-----------------------?
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 1,
            #                         value = Bmax)  
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 2,
            #                         value = Bmax)  
            FPD.add_yFixity(index=2,
                            value=Bmax)
            #-----------------------?
            #print 'Adding y - fixity!'
            FPD.add_yFixity(index=3,
                            value=Bmax)
            FPD.add_yFixity(index=4,
                            value=Bmax)
            FPD.add_yFixity(index=5,
                            value=Bmax)
            FPD.add_yFixity(index=6,
                            value=Bmax)
            FPD.add_yFixity(index=7,
                            value=Bmax)
            #-----------------------?
            FPD.add_yFixity(index=8,
                            value=Bmax)
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 8,
            #                         value = Bmax)  
            #FPD.add_yVertexConstraint(kind='max',
            #                         index = 9,
            #                         value = Bmax)  
            #-----------------------?
            #
            #
            #******************************************************************
            #    
            #FPD.add_E1(kind='LS', weight = 1.)
            #FPD.add_E2(kind='LS', weight = .5)
            #FPD.add_E3(kind='LS', weight = .5) #always commented out TLM OCT 27 2017
            #FPD.add_ArcLengthApprox(kind='LS', weight = 1.5)
            #
            #
            #**********************************************************************
            # DWL fairness norms
            FPD.add_E1(kind='LS', weight = 1.5/E1_0)
            FPD.add_E2(kind='LS', weight = .5/E2_0)
            FPD.add_E3(kind='LS', weight = .05/E3_0)
            FPD.add_ArcLengthApprox(kind='LS', weight = 1.5/S_0)
            #
            #**********************************************************************
            #Setup and Solve:
            #
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            # make curve
            print 'Improving first stage DWL solve with scaling factors'
            Lspline.optimize(stop=25)
            return Lspline
            
            
        
        #
        #**********************************************************************
        # DWL save init state
        #inicurve = copy.deepcopy(curve)
        #func1_for_fix = fix1()
        #
        #**********************************************************************
        # DWL First solve:
        #Lspline = stage1(curve)
        if self.hullformtype == 'osv':
            print  'DWL to end at transom '
            Lspline = Transom_DWL(curve)
            stop = uopt.ckLspline(Lspline)
            if not stop:
                Lspline = Transom_DWL_fixit(curve)
        else:
            Lspline = stage1(curve)
        #
        #ckit = False
        #        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
        #            #func1_for_fix = fix1(curve = inicurve)
        #            Lspline = stage1(curve)
        #        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
        #            #func1_for_fix = fix1(curve = inicurve)
        #            Lspline = func1_for_fix(curve = inicurve)
        #    ckit = False
        #        else:    
        #            Lcopy = copy.deepcopy(Lspline)
        #            ckit = True
        #
        #
        #**********************************************************************
        # DWL Second solve:
        #Lspline = stage2(Lspline)print 'attempt to improve upon the normalization of stage 1'
        #        Lspline = improve_stage1(Lspline)
        #        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
        #            if ckit:
        #                print 'reverting to 1st stage'
        #                Lspline = Lcopy
        #                #Lcopy2 = copy.deepcopy(Lspline)
        #                print 'attempt to improve upon the normalization of stage 1'
        #                #Lspline = improve_stage1(Lspline)
        #                #if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
        #                #    print 'reverting to stage 1 solution'
        #                #    Lspline = Lcopy2
        #                #else:
        #                #    pass
        #            else:
        #                print 'attempt to improve base curve stages 1 and 2 failed'
        #                Lspline = func1_for_fix(curve = inicurve)
                #                #initial curve, alt try
                #                Lspline = improve_stage1(Lspline)
                #                if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
                #                    Lspline = Lcopy2
                #                else:
                #                    
                #                    Lcopy2 = copy.deepcopy(Lspline)
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
        
        
        if self.hullformtype=='osv':
            self.split_DWL_curves()
            
        
        
        self.results_dict['Awp'] = self.DWL2D.area.value * 2.
        self.results_dict['lfwl'] = self.stations.FOWL[1] - \
                                            self.stations.FOWL[0]
        self.results_dict['lfwl_data'] = 'min flat length'
        
        
        #self.results_dict['LCG'] = self.DWL.Xc
        return

    def compute_cLProfile(self, LocMaxSection=None, Xc=None,
                          flat=False, height=None, alphab=0.,alphae=0.,
                          drop_aft=.25):
        print '\ncLProfile'
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
        ye      = 0.+self.MaxDraft*drop_aft
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
                FPD.add_XcConstraint(kind='LS', 
                                     value=self.LCG)
            else:
                print 'Centerplane setting Xc constraint to ',Xc
                FPD.add_XcConstraint(kind='LS', 
                                     value=Xc)
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
            #--------------------------------------
            # beyond baseline 
            # ++
            #--------------------------------------
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 1,
                                     value = 0.)
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 2,
                                     value = 0.)
            #--------------------------------------
            # +++
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 3,
                                     value = 0.)
            #--------------------------------------
            # ++
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 8,
                                     value = xe)
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 9,
                                     value = xe)
            #--------------------------------------
            #
            #******************************************************************
            # x-flat CPkeelProfile  (4,5,6 + optional:  keep 3,7 out)
            # stage 1 solve
#            FPD.add_xVertexConstraint(kind='max',
#                                     index = 3,
#                                     value = self.stations.FOCP[0])
#            FPD.add_xVertexConstraint(kind='max',
#                                     index = 4,
#                                     value = self.stations.FOCP[0])
            #FPD.add_xVertexConstraint(kind='min',
            #                         index = 5,
            #                         value = self.stations.FOCP[0])
            #FPD.add_xVertexConstraint(kind='max',
            #                         index = 5,
            #                         value = self.stations.FOCP[1])
            #FPD.add_xVertexConstraint(kind='equality',
            #                         index = 6,
            #                         value = self.stations.FOCP[1])
#            FPD.add_xVertexConstraint(kind='min',
#                                     index = 7,
#                                     value = self.stations.FOCP[1])
            ##**************************************************************
            # CPkeelProfile  
            #print 'Adding x - fixity!'
            FPD.add_xFixity(index=5,
                            value=self.stations.FOCP[0])
            FPD.add_xFixity(index=6,
                            value=self.stations.FOCP[1])
            #
            #**********************************************************************
            # y-flat CPkeelProfile (3,4,5,6,7) stage 1 solve
#            FPD.add_yVertexConstraint(kind='max',
#                                     index = 3,
#                                     value = Dmax)
            #FPD.add_yVertexConstraint(kind='equality',
            #                         index = 4,
            #                         value = Dmax)
            #FPD.add_yVertexConstraint(kind='equality',
            #                         index = 5,
            #                         value = Dmax)
            #FPD.add_yVertexConstraint(kind='equality',
            #                         index = 6,
            #                         value = Dmax)
            #FPD.add_yVertexConstraint(kind='equality',
            #                         index = 7,
            #                         value = Dmax)
            #
            #**********************************************************************
            # y-flat CPkeelProfile  (use fixity) stage 1 solve 
            #?-------------------------------?
            #            FPD.add_yVertexConstraint(kind='max',
            #                                     index = 2,
            #                                     value = Dmax)
            #?-------------------------------?
            #FPD.add_yFixity(index=2,
            #                value=Dmax)
            #?-------------------------------?
            #print 'Adding y - fixity!'
            FPD.add_yFixity(index=3,
                            value=Dmax)
            FPD.add_yFixity(index=4,
                            value=Dmax)
            FPD.add_yFixity(index=5,
                            value=Dmax)
            FPD.add_yFixity(index=6,
                            value=Dmax)
            FPD.add_yFixity(index=7,
                            value=Dmax)
            #?-------------------------------?
            #FPD.add_yFixity(index=8,
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
            print 'Doing first stage CPkeel solve'
            Lspline.optimize(stop=100)
            return Lspline
        
        
        def OSV_1(curve):
            """Lspline solver 
            
                for the 
                
                CProfile keel curve
                
             (intended for a transom sterned hull)
            
            
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
                FPD.add_XcConstraint(kind='LS', 
                                     value=self.LCG) #cLProfile flat to bulb...
                FPD.add_XcConstraint(kind='max', 
                                     value=self.LCG*1.05)
                FPD.add_XcConstraint(kind='min', 
                                     value=self.LCG*0.95)
            else:
                print 'Centerplane setting Xc constraint to ',Xc
                FPD.add_XcConstraint(kind='LS', 
                                     value=Xc)
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
                                    value       = 0.)
#            FPD.add_AngleConstraint(
#                                    kind        = 'LS',
#                                    location    = 0.1,
#                                    value       = 75.)
#            FPD.add_CurvatureConstraint(
#                                    kind        = 'LS',
#                                    location    = .08,
#                                    value       = 0.)
            #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
            #--------------------------------------
            # beyond baseline 
            # ++
            #--------------------------------------
            #            FPD.add_xVertexConstraint(kind='LS',
            #                                     index = 1,
            #                                     value = 0.)
            #            FPD.add_xVertexConstraint(kind='LS',
            #                                     index = 2,
            #                                     value = 0.)
            #--------------------------------------
            # +++
            #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
            # TLM, straight flat stem
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 3,
                                     value = 0.)
            #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
            #--------------------------------------
            # ++
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 8,
                                     value = xe)
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 9,
                                     value = xe)
            #--------------------------------------
            #
            #******************************************************************
            # x-flat CPkeelProfile  (4,5,6 + optional:  keep 3,7 out)
            # stage 1 solve
#            FPD.add_xVertexConstraint(kind='max',
#                                     index = 3,
#                                     value = self.stations.FOCP[0])
#            FPD.add_xVertexConstraint(kind='max',
#                                     index = 4,
#                                     value = self.stations.FOCP[0])
            #FPD.add_xVertexConstraint(kind='min',
            #                         index = 5,
            #                         value = self.stations.FOCP[0])
            #FPD.add_xVertexConstraint(kind='max',
            #                         index = 5,
            #                         value = self.stations.FOCP[1])
            #FPD.add_xVertexConstraint(kind='equality',
            #                         index = 6,
            #                         value = self.stations.FOCP[1])
#            FPD.add_xVertexConstraint(kind='min',
#                                     index = 7,
#                                     value = self.stations.FOCP[1])
            ##**************************************************************
            # CPkeelProfile  
            #print 'Adding x - fixity!'
            FPD.add_xFixity(index=5,
                            value=self.stations.FOCP[0])
            FPD.add_xFixity(index=6,
                            value=self.stations.FOCP[1])
            #
            #**********************************************************************
            # y-flat CPkeelProfile (3,4,5,6,7) stage 1 solve
#            FPD.add_yVertexConstraint(kind='max',
#                                     index = 3,
#                                     value = Dmax)
            #FPD.add_yVertexConstraint(kind='equality',
            #                         index = 4,
            #                         value = Dmax)
            #FPD.add_yVertexConstraint(kind='equality',
            #                         index = 5,
            #                         value = Dmax)
            #FPD.add_yVertexConstraint(kind='equality',
            #                         index = 6,
            #                         value = Dmax)
            #FPD.add_yVertexConstraint(kind='equality',
            #                         index = 7,
            #                         value = Dmax)
            #
            #**********************************************************************
            # y-flat CPkeelProfile  (use fixity) stage 1 solve 
            #?-------------------------------?
            #            FPD.add_yVertexConstraint(kind='max',
            #                                     index = 2,
            #                                     value = Dmax)
            #-------------------------------
            # for good blending of the bulb:
            # try and keep the bottom flat past
            # where the fwd fairness curve intersects the cLProfile
            
            
            #this is 'just too much?':
            #FPD.add_yFixity(index=1,
            #                value=Dmax)
            
            FPD.add_yFixity(index=2,
                            value=Dmax)
            #            FPD.add_relative_xVertexConstraint(value = 0., 
            #                                               weight = 1.0, 
            #                                               index = 3, 
            #                                               index2 = 2,
            #                                               seperation=0.)
            #-------------------------------
            #print 'Adding y - fixity!'
            FPD.add_yFixity(index=3,
                            value=Dmax)
            FPD.add_yFixity(index=4,
                            value=Dmax)
            FPD.add_yFixity(index=5,
                            value=Dmax)
            FPD.add_yFixity(index=6,
                            value=Dmax)
            FPD.add_yFixity(index=7,
                            value=Dmax)
            #-------------------------------
            #FPD.add_yFixity(index=8,
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
            print 'Doing first stage CPkeel solve'
            Lspline.optimize(stop=100)
            return Lspline
        
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
            if self._verbose: print 'OSV_bbow solver Oct 14, 2018'
            #
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
                if self._verbose: 'CPK OSV bbow, equality LCG'
                FPD.add_XcConstraint(kind='LS', 
                                     value=self.LCG*.75,
                                     weight = 100.) #cLProfile flat to bulb...
                #FPD.add_XcConstraint(kind='max', 
                #                     value=self.LCG*.95)
                #FPD.add_XcConstraint(kind='min', 
                #                     value=self.LCG*0.8)
            else:
                print 'Centerplane setting Xc constraint to ',Xc
                FPD.add_XcConstraint(kind='LS', 
                                     value=Xc)
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
                                    value       = 0.)
            #-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
            #--------------------------------------
            # beyond baseline 
            # ++
            #--------------------------------------
            #            FPD.add_xVertexConstraint(kind='LS',
            #                                     index = 1,
            #                                     value = 0.)
            #            FPD.add_xVertexConstraint(kind='LS',
            #                                     index = 2,
            #                                     value = 0.)
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
            #
            ##**************************************************************
            # CPkeelProfile  
            #print 'Adding x - fixity!'
            #tlmOct 14 2018 tweak this:
            FPD.add_xFixity(index=5, #5,
                            value=self.stations.FOCP[0])
            FPD.add_xFixity(index=6, #6,
                            value=self.stations.FOCP[1])
            #
            #**********************************************************************
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
            if self._verbose: print 'OSV_bbow solver Oct 14, 2018'
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
            #
            # not using weights because this
            # solver is used after a fail
            # so the E data is no good.
            #(alt it could be used after a good solve
            # to improve further, in which case we would
            # want to use real E magnatiude normalization)
            #
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
                if self._verbose: print 'CPK OSV bbow, equality LCG'
                FPD.add_XcConstraint(kind='LS', 
                                     value=self.LCG*.75,
                                     weight = 100.) #cLProfile flat to bulb...
                #FPD.add_XcConstraint(kind='max', 
                #                     value=self.LCG*.95)
                #FPD.add_XcConstraint(kind='min', 
                #                     value=self.LCG*0.8)
            else:
                if self._verbose: print 'Centerplane setting Xc constraint to ',Xc
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
                            value=self.stations.FOCP[0])
            FPD.add_xFixity(index=6, #6,
                            value=self.stations.FOCP[1])
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
            if self._verbose: print 'Doing 2nd Stage CPkeel solve (include bbow weighting)'
            Lspline.optimize(stop=50)
            return Lspline
        
        
        
        
            
        def OSV_THB(curve):
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
                FPD.add_XcConstraint(kind='LS', 
                                     value=self.LCG)
            else:
                print 'Centerplane setting Xc constraint to ',Xc
                FPD.add_XcConstraint(kind='LS', 
                                     value=Xc)
            #
            #
            #**********************************************************************
            # Nose to 90
            FPD.add_verticalAngleConstraint(kind='LS',
                                            location = 0.,
                                            value = 90.)
            FPD.add_CurvatureConstraint(kind='equality',
                                    location = 0.,
                                    value = 0.)
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
            #Lspline.optimize(stop=100)
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
            # CPK
            # initialize
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve)
            #
            #**********************************************************************
            # area  CPKeelProfile stage 2
            FPD.add_AreaConstraint(kind='LS',
                                   value = area,
                                   weight=1.5/E1_0)
            #FPD.add_AreaConstraint(kind='min',
            #                       value = area-.25*area)
            #FPD.add_AreaConstraint(kind='max',
            #                       value = area+.25*area)
            #FPD.add_AreaConstraint(kind='equality',
            #                       value = area)
            #
            #**********************************************************************
            # Xc CPKeelProfile center of waterplane needs to be a rule
            # consider rules with a max and min!
            if Xc is None:
                FPD.add_XcConstraint(kind='LS', 
                                     value=self.LCG,
                                     weight=.5/E1_0)
            else:
                print 'Centerplane setting Xc constraint to ',Xc
                FPD.add_XcConstraint(kind='LS', 
                                     value=Xc,
                                     weight=.5/E1_0)
            #
            #**********************************************************************
            # Nose behavior
            # - weight the slope at the nose?
            #            FPD.add_relative_xVertexConstraint(
            #                                     kind='LS',
            #                                     location=None,
            #                                     index = 0,
            #                                     index2 = 1,
            #                                     value = 0.,
            #                                     weight=1.)
            #print 'Adding x - nose fixity!'
            #FPD.add_xFixity(index=1,
            #                value=0.)
            #
            #**********************************************************************
            # cL Keel Profile - stay within parameters
            FPD.add_xVertexConstraint(kind='min',
                                     index = 1,
                                     value = 0.)
            FPD.add_xVertexConstraint(kind='min',
                                     index = 2,
                                     value = 0.)
            FPD.add_xVertexConstraint(kind='max',
                                     index = 8,
                                     value = 1.1*xe)
            FPD.add_xVertexConstraint(kind='max',
                                     index = 9,
                                     value = 1.1*xe)
            #FPD.add_xFixity(index=9,
            #                value=self.LengthDWL)
            #
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
            FPD.add_yVertexConstraint(kind='min',
                                     index = 8,
                                     value = 0.)
            FPD.add_yVertexConstraint(kind='min',
                                     index = 9,
                                     value = 0.)
#            #
#            # ensure no wild runoff
#            #FPD.add_yVertexConstraint(kind='max',
#            #                         index = 1,
#            #                         value = 1.05*Dmax)
#            FPD.add_yVertexConstraint(kind='max',
#                                     index = 2,
#                                     value = 1.05*Dmax)
#            FPD.add_yVertexConstraint(kind='max',
#                                     index = 8,
#                                     value = 1.05*Dmax)
#            FPD.add_yVertexConstraint(kind='max',
#                                     index = 9,
#                                     value = 1.05*Dmax)
            #
            #******************************************************************
            # CPkeelProfile
            # x-flat (4,5,6 + 3,7 )
            # stage 2solve
#            FPD.add_xVertexConstraint(kind='max',
#                                     index = 3,
#                                     value = self.stations.FOCP[0])
            FPD.add_xVertexConstraint(kind='max',
                                     index = 4,
                                     value = self.stations.FOCP[0])
            #FPD.add_xVertexConstraint(kind='min',
            #                         index = 5,
            #                         value = self.stations.FOCP[0])
            #FPD.add_xVertexConstraint(kind='max',
            #                         index = 5,
            #                         value = self.stations.FOCP[1])
            #FPD.add_xVertexConstraint(kind='equality',
            #                         index = 6,
            #                         value = self.stations.FOCP[1])
            FPD.add_xVertexConstraint(kind='min',
                                     index = 7,
                                     value = self.stations.FOCP[1])
            #
            #
            #print 'Adding x - fixity!'
            FPD.add_xFixity(index=5,
                            value=self.stations.FOCP[0])
            FPD.add_xFixity(index=6,
                            value=self.stations.FOCP[1])
            #
            #**********************************************************************
            # CPkeelProfile
            # y-flat (3,4,5,6,7) stage 2 solve
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 3,
            #                                     value = Dmax)
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 3,
            #                                     value = Dmax+self.tol*.25)  
            #            FPD.add_yVertexConstraint(kind='equality',
            #                                     index = 4,
            #                                     value = Dmax)
            #            FPD.add_yVertexConstraint(kind='equality',
            #                                     index = 6,
            #                                     value = Dmax)
            #            FPD.add_yVertexConstraint(kind='equality',
            #                                     index = 5,
            #                                     value = Dmax)
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 7,
            #                                     value = Dmax)
            #            FPD.add_yVertexConstraint(kind='LS',
            #                                     index = 7,
            #                                     value = Dmax+self.tol*.25)  
            #print 'Adding y - fixity!' stage 2 solve
            FPD.add_yFixity(index=3,
                            value=Dmax)
            FPD.add_yFixity(index=4,
                            value=Dmax)
            FPD.add_yFixity(index=5,
                            value=Dmax)
            FPD.add_yFixity(index=6,
                            value=Dmax)
            FPD.add_yFixity(index=7,
                            value=Dmax)
            #
            #
            #**********************************************************************
            # CPK fairness norms
            FPD.add_E1(kind='LS', weight = .5/E1_0)
            FPD.add_E2(kind='LS', weight = .5/E2_0)
            FPD.add_E3(kind='LS', weight = .5/E3_0)
            FPD.add_ArcLengthApprox(kind='LS', weight = .5/S_0)
            #
            #**********************************************************************
            #Setup and Solve:
            #
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            print 'Doing 2nd stage CPkeel solve'
            Lspline.optimize(stop=25)
            return Lspline
        
        def fix1():
            #
            #**********************************************************************
            # continuation passing style curve fixer, CPkeel
            def doit(Lspline=None,curve=None):
                """CPK fix
                """
                #
                #**********************************************************************
                # get weights
                if Lspline is not None:
                    norms = self.package_norms(Lspline)
                    E1_0 = norms[0]
                    E2_0 = norms[1]
                    E3_0 = norms[2]
                    S_0  = norms[3]
                    curve = Lspline.curve
                else:
                    E1_0 = 1.
                    E2_0 = 1.
                    E3_0 = 1.
                    S_0  = 1.
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
                #if Xc is None:
                #    FPD.add_XcConstraint(kind='LS', 
                #                         value=self.LCG)
                #else:
                #    print 'Centerplane setting Xc constraint to ',Xc
                #    FPD.add_XcConstraint(kind='equality', value=Xc)
                #
                #******************************************************************
                # x-flat CPkeelProfile  (4,5,6 + optional:  keep 3,7 out)
                # stage fix solve
                #                FPD.add_xVertexConstraint(kind='LS',
                #                                         index = 3,
                #                                         value = self.stations.FOCP[0])
                #                FPD.add_xVertexConstraint(kind='equality',
                #                                         index = 4,
                #                                         value = self.stations.FOCP[0])
                #                FPD.add_xVertexConstraint(kind='min',
                #                                         index = 5,
                #                                         value = self.stations.FOCP[0])
                #                FPD.add_xVertexConstraint(kind='max',
                #                                         index = 5,
                #                                         value = self.stations.FOCP[1])
                #                FPD.add_xVertexConstraint(kind='equality',
                #                                         index = 6,
                #                                         value = self.stations.FOCP[1])
                #                FPD.add_xVertexConstraint(kind='LS',
                #                                         index = 7,
                #                                         value = self.stations.FOCP[1])
                FPD.add_xVertexConstraint(kind='max',
                                         index = 3,
                                         value = self.stations.FOCP[0])
                FPD.add_xVertexConstraint(kind='max',
                                         index = 4,
                                         value = self.stations.FOCP[0])
                FPD.add_xVertexConstraint(kind='min',
                                         index = 7,
                                         value = self.stations.FOCP[1])
                #**************************************************************
                #print 'Adding x - fixity!'
                FPD.add_xFixity(index=5,
                                value=self.stations.FOCP[0])
                FPD.add_xFixity(index=6,
                                value=self.stations.FOCP[1])
                #
                #**********************************************************************
                # y-flat CPkeelProfile (4,5,6,) 
                # stage fix solve
                #FPD.add_yVertexConstraint(kind='LS',
                #                         index = 3,
                #                         value = Dmax)  
                #FPD.add_yVertexConstraint(kind='max',
                #                         index = 3,
                #                         value = Dmax+self.tol*.25) 
                #                FPD.add_yVertexConstraint(kind='equality',
                #                                         index = 4,
                #                                         value = Dmax)
                #                FPD.add_yVertexConstraint(kind='equality',
                #                                         index = 5,
                #                                         value = Dmax)
                #                FPD.add_yVertexConstraint(kind='equality',
                #                                         index = 6,
                #                                         value = Dmax)
                #FPD.add_yVertexConstraint(kind='LS',
                #                         index = 7,
                #                         value = Dmax)  
                #FPD.add_yVertexConstraint(kind='max',
                #                         index = 7,
                #                         value = Dmax+self.tol*.25) 
                #print 'Adding y - fixity!'
                FPD.add_yFixity(index=3,
                                value=Dmax)
                FPD.add_yFixity(index=4,
                                value=Dmax)
                FPD.add_yFixity(index=5,
                                value=Dmax)
                FPD.add_yFixity(index=6,
                                value=Dmax)
                FPD.add_yFixity(index=7,
                                value=Dmax)
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
                print 'Repairing cL-Keel-Profile curve (old repair)'
                Lspline.optimize(stop=25)
                return Lspline
            return doit
        
        def fix_improve_stage1(Lspline):
            """Should be passed a successful instance of stage1
            """
            # CPK
            norms = self.package_norms(Lspline)
            E1_0 = norms[0]
            E2_0 = norms[1]
            E3_0 = norms[2]
            S_0  = norms[3]
            curve = Lspline.curve
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
                                     value=self.LCG,
                                     weight=.5/E1_0)
            else:
                print 'Centerplane setting Xc constraint to ',Xc
                FPD.add_XcConstraint(kind='LS', 
                                     value=Xc,
                                     weight=.5/E1_0)
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
            #--------------------------------------
            # beyond baseline 
            # ++
            #--------------------------------------
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 1,
                                     value = 0.,
                                     weight=.5/E1_0)
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 2,
                                     value = 0.,
                                     weight=.5/E1_0)
            #--------------------------------------
            # +++
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 3,
                                     value = 0.,
                                     weight=.5/E1_0)
            #--------------------------------------
            # ++
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 8,
                                     value = xe)
            FPD.add_xVertexConstraint(kind='LS',
                                     index = 9,
                                     value = xe,
                                     weight=.5/E1_0)
            #--------------------------------------
            #
            #******************************************************************
            # x-flat CPkeelProfile  (4,5,6 + optional:  keep 3,7 out)
            # stage 1 solve
#            FPD.add_xVertexConstraint(kind='max',
#                                     index = 3,
#                                     value = self.stations.FOCP[0])
#            FPD.add_xVertexConstraint(kind='max',
#                                     index = 4,
#                                     value = self.stations.FOCP[0])
            #FPD.add_xVertexConstraint(kind='min',
            #                         index = 5,
            #                         value = self.stations.FOCP[0])
            #FPD.add_xVertexConstraint(kind='max',
            #                         index = 5,
            #                         value = self.stations.FOCP[1])
            #FPD.add_xVertexConstraint(kind='equality',
            #                         index = 6,
            #                         value = self.stations.FOCP[1])
#            FPD.add_xVertexConstraint(kind='min',
#                                     index = 7,
#                                     value = self.stations.FOCP[1])
            ##**************************************************************
            # CPkeelProfile  
            #print 'Adding x - fixity!'
            FPD.add_xFixity(index=5,
                            value=self.stations.FOCP[0])
            FPD.add_xFixity(index=6,
                            value=self.stations.FOCP[1])
            #
            #**********************************************************************
            # y-flat CPkeelProfile (3,4,5,6,7) stage 1 solve
#            FPD.add_yVertexConstraint(kind='max',
#                                     index = 3,
#                                     value = Dmax)
            #FPD.add_yVertexConstraint(kind='equality',
            #                         index = 4,
            #                         value = Dmax)
            #FPD.add_yVertexConstraint(kind='equality',
            #                         index = 5,
            #                         value = Dmax)
            #FPD.add_yVertexConstraint(kind='equality',
            #                         index = 6,
            #                         value = Dmax)
            #FPD.add_yVertexConstraint(kind='equality',
            #                         index = 7,
            #                         value = Dmax)
            #
            #**********************************************************************
            # y-flat CPkeelProfile  (use fixity)
            #?-------------------------------?
            #            FPD.add_yVertexConstraint(kind='max',
            #                                     index = 2,
            #                                     value = Dmax)
            #?-------------------------------?
            #FPD.add_yFixity(index=2,
            #                value=Dmax)
            #?-------------------------------?
            #print 'Adding y - fixity!'
            FPD.add_yFixity(index=3,
                            value=Dmax)
            FPD.add_yFixity(index=4,
                            value=Dmax)
            FPD.add_yFixity(index=5,
                            value=Dmax)
            FPD.add_yFixity(index=6,
                            value=Dmax)
            FPD.add_yFixity(index=7,
                            value=Dmax)
            #?-------------------------------?
            #FPD.add_yFixity(index=8,
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
            #
            #bad:--------------------------
            #FPD.add_E1(kind='LS', weight = 1.5/E1_0)
            #FPD.add_E2(kind='LS', weight = .5/E2_0)
            #FPD.add_E3(kind='LS', weight = .05/E3_0)
            #--------------------------------
            #FPD.add_E1(kind='LS', weight = .5/E1_0)
            #FPD.add_E2(kind='LS', weight = .5/E2_0)
            #FPD.add_E3(kind='LS', weight = .001/E3_0)
            #FPD.add_ArcLengthApprox(kind='LS', weight = 1.5/S_0)
            #
            FPD.add_ArcLengthApprox(kind='LS', weight = 1.5)
            #
            #**********************************************************************
            #Setup and Solve:
            #
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            print 'Improving first stage CPkeel solve with scaling data'
            Lspline.optimize(stop=100)
            return Lspline
        
        
        
        #
        #**********************************************************************
        # CPK save init state
        inicurve = copy.deepcopy(curve)
        func1_for_fix = fix1()
        #
        #**********************************************************************
        # First solve:
        #        Lspline = stage1(curve)
        #        #
        #        ckit = False
        #        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
        #            Lspline = func1_for_fix(curve = inicurve)
        #            ckit = False
        #        else:    
        #            Lcopy = copy.deepcopy(Lspline)
        #            ckit = True
        #            
        #            
        #        curve = Lspline.curve
            
        #Lspline = OSV_1(curve)
        Lspline = OSV_bbow(curve)
        #
        #
        #**********************************************************************
        # Second solve:
        ckit = False
        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
            Lspline = OSV_bbow2(Lspline)
        #Lspline = stage2(Lspline)
        #        Lspline = fix_improve_stage1(Lspline)
        #        if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
        #            if ckit:
        #                print 'reverting to 1st stage'
        #                Lspline = Lcopy
        #                print 'attempt to improve upon the normalization of stage 1'
        #                #Lcopy2 = copy.deepcopy(Lspline)
        #                # stage1, improved
        #                #Lspline = fix_improve_stage1(Lspline)
        #                #if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
        #                #    Lspline = Lcopy2
        #                #else:
        #                #    pass
        #                
        #            else:
        #                print 'attempt to improve base curve stages 1 and 2 failed'
        #                Lspline = func1_for_fix(curve = inicurve)
        #                #                #initial curve, alt try
        #                #                Lspline = fix_improve_stage1(Lspline)
        #                #                if Lspline.error_code is 'max_iter' or Lspline.error_code is 'NAN':
        #                #                    Lspline = Lcopy2
        #                #                else:
        #                #                    
        #                #                    Lcopy2 = copy.deepcopy(Lspline)
        #
        #**********************************************************************
        # save CLKeelProfile
        self.Lsplines.CProfile=Lspline
        self.CProfile2D = Lspline.curve
        #
        #**********************************************************************
        # lift to 3D
        self.CProfile = self.lifting_2DCProfile(Lspline.curve.vertices, 
                                                elevation,k,nump)
        
        #        if self.hullformtype == 'osv':
        #            self.split_cLprofile_curves(OSV=True)
        #        else:
        #            self.split_cLprofile_curves(OSV=False)
        # yet more implicit state...
        self.split_cLprofile_curves()
        
        
        self.results_dict['Acp'] = self.CProfile2D.area.value
        self.results_dict['lfcp'] = self.stations.FOCP[1] - \
                                            self.stations.FOCP[0]
        self.results_dict['lfcpdata'] = 'min flat length'
        
        
        #self.results_dict['LCG'] = self.DWL.Xc
        return
    
    def lifting_2DCProfile(self, 
                           verts, elevation,k,nump):
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
        vnew  = frame.xy_to_zyx_elevate_noflp( verts,
                                            elevation=elevation)
        return spline.Bspline(vnew, k, nump)

    def compute_stern_bulb_profile(self):
        """mirrored
        """
        k=self.k
        nump=self.nump
        xb = self.stations.Stern_Profile[0]
        xe = self.stations.Stern_Profile[1]
        yb = self.MaxDraft
        ye = 0.

        return

    def compute_stern_profile(self,area=None):
        """NOT USED
        
        mirrored
        cL plane == z-y plane
        """
        k=self.k
        nump=self.nump
        xb = self.stations.Stern_Profile[0]
        xe = self.stations.Stern_Profile[1]
        a =self.stations.Stern_Profile[2]
        b =self.stations.Stern_Profile[3]
        yb = self.MaxDraft
        ye = 0.
        alphab = 0.#-5.
        alphae = 0.#-15.
        Cab_given = 0.
        Cae_given = 0.
        slope = 'down'

        #Dmax = self.DWL.CurvePoint(a)
        if area is None:
            areafac = (2./3.)
            area = areafac*((xb-xe)*(yb-ye))
        xb3d = self.DWL.CurvePoint(a)
        xb = np.asarray([0.,xb3d[2]])
        xe3d = self.DWL.CurvePoint(b)
        xe = np.asarray([xe3d[0],xe3d[2]])
        return

    def plot_primary_curves(self):
        """
        self.SAC.plot3DmultiList([self.SAC_3dvis,
                                  self.DWL,
                                  self.CProfile],[])
                
                
        self.SAC.plotcurve_detailed()
        self.DWL2D.plotcurve()
        self.CProfile2D.plotcurve()
        
        ax = self.SAC.plotcurve_detailed()
        ax = self.DWL2D.plotcurve_detailed(canvas = ax)
        ax = self.CProfile2D.plotcurve_detailed(canvas = ax)

        #"""
        
        if self.hullformtype == 'osv':
            fwdtcurve = self.CPK_fwd
            afttcurve = self.DWL_aft
        else:
            fwdtcurve = self.CPK_fwd
            afttcurve = self.CPK_aft
        
        if self.hullformtype !='osv':
        
            l2 = [fwdtcurve]
            l2.append(self.bowtransition)
            l2.append(self.bowfairning)
            for cv in self.mshp_sections:
                l2.append(cv)
            l2.append(self.sterntransition)
            l2.append(self.sternfairing)
            l2.append(afttcurve)
    
            self.DWL.plot3DmultiList([self.DWL,
                                      self.CProfile],
                                      l2,
                                      view='x-vertical')
            
#            self.DWL.plot3DmultiList([self.DWL,
#                                      self.CProfile],
#                                      l2,
#                                      view='y-vertical')
#            
#            self.DWL.plot3DmultiList([self.DWL,
#                                      self.CProfile],
#                                      l2,
#                                      view='z-vertical')
        else:
            l2 = [fwdtcurve]
            l2.append(self.bowtransition)
            l2.append(self.bowfairning)
            for cv in self.mshp_sections:
                l2.append(cv)
            l2.append(self.sterntransition)
            l2.append(self.sternfairing)
            l2.append(afttcurve)
    
            self.DWL.plot3DmultiList([self.DWL_OSV,
                                      self.CPK_cdr_b],
                                      l2,
                                      view='x-vertical')
            
#            self.DWL.plot3DmultiList([self.DWL_OSV,
#                                      self.CPK_cdr_b],
#                                      l2,
#                                      view='y-vertical')
#            
#            self.DWL.plot3DmultiList([self.DWL_OSV,
#                                      self.CPK_cdr_b],
#                                      l2,
#                                      view='z-vertical')
            
        return
    
    
    def plot_hull_transverse_curves_and_longitudinal_curves(self):
        """
        self.SAC.plot3DmultiList([self.SAC_3dvis,
                                  self.DWL,
                                  self.CProfile],[])

        #"""
        #        self.DWL.plot3DmultiList([self.DWL,
        #                                  self.CProfile],
        #                                  [],
        #                                  limx=100.,
        #                                  limy=100.,
        #                                  limz=100.)
        #        self.DWL.plot3DmultiList([],
        #                                  self.mshp_sections,
        #                                  limx=100.,
        #                                  limy=100.,
        #                                  limz=100.)
        
        l2 = []
        l2.append(self.bowtransition)
        l2.append(self.bowfairning)
        for cv in self.mshp_sections:
            l2.append(cv)
        l2.append(self.sterntransition)
        l2.append(self.sternfairing)

        self.DWL.plot3DmultiList(self.lcurvenet,
                                 self.tcurvenet,
                                  view='x-vertical')
        
        #        self.DWL.plot3DmultiList(self.lcurvenet,
        #                                 self.tcurvenet,
        #                                  view='y-vertical')
        #        
        #        
        #        self.DWL.plot3DmultiList(self.lcurvenet,
        #                                 self.tcurvenet,
        #                                  view='z-vertical')
        return


    def insert_ctrl_vertices(self, curve, n):
        knots = np.linspace(0.,1.,n+2)[1:-1]
        for knot in knots:
            curve.knot_insertion(knot)
        return curve
    
    def match_curve_knots(self, curve_to_fix, good_curve, tol=.1):
        matched_curve = uopt.match_knots(
                    initial_curve = copy.copy(curve_to_fix),
                    curve_w_good_properties=good_curve,
                    TOL=tol)
        return matched_curve
    
    
    

    def split_cLprofile_curves(self):
        """Split the end(s) of the CP Curve 
        
        Make these fwd and aft ends have == # vertices and same knot vector
        as the mid ship sections 
        
        (*note, if we are making an osv, there is no aft cLProfile curve
        instead the cLProfile longitudinal section runs all the way
        to the standard aft termination point of the cLProfile curve)
        
        and the center section must have 
        the same knot vector as the longitudinals
        
        
        Notes:
        --------------------
            *the CPKeel curve is defined fwd, (starting from the origin at 0.) 
                                      to aft, (running to the max lwl)
            *careful with parameterizations!
        
        developer shortcuts:
        --------------------
        
        
            import utility_optimization as uopt
            import SpecializedCurveSolvers as scsolvers
            from ADILS import LazyFPDConstraintInterface
            import curve as spline
            
            self = SD.hull
            cpk = self.CProfile2D
            
            
            
            
            bb, be = self.stations.Bow_Fairing[2:]
            sb, se = self.stations.Stern_Fairing[2:]
            
            sb = sb + 0.6*(1.-sb)
            
            
            CPK_fwd, CPK_cdr_a = cpk.CurveSplit(bb)
            
            #self.CPK_cdr_b, self.CPK_aft = cpk.CurveSplit(sb)
            
            #commenting out DWL mods.  Nov 4 2017
            #self.DWL_fwd, self.DWL_cdr_a = dwl.CurveSplit(be) #Nov  2017 - this is not used?
            # DWL does not need to be split because no longitudinal curves intersect it.
            #
            
            CPK_cdr_b, CPK_aft = CPK_cdr_a.CurveSplit(sb)
        
            elevation = self.MaxDraft
            
            
            side check, 
            ---------------------
            Question: if we lift a whole CP curve, and lift a split off piece of it,
                            do they perfectly match?
            Answer:  yes, yes they do.
            
            tc = self.lifting_2DCProfile(self.CProfile2D.vertices,
                                    elevation,
                                    4,
                                    30)
            
            
            hh = self.lifting_2DCProfile(CPK_cdr_b.vertices,
                                    elevation,
                                    4,
                                    30)
            ---------------------------
            
            matched_curve = uopt.match_knots(initial_curve = CPK_cdr_b,
                                             curve_w_good_properties=self.CProfile2D,
                                             TOL=100.1)
            
            
            lc = self.lifting_2DCProfile(matched_curve.vertices,
                                         elevation = self.MaxDraft,4,30)
            lc.plot3DmultiList([hh],[lc])
        
        
            self.CPK_aft.plot3DmultiList(
                                    [self.CPK_aft,self.CPK_fwd],
                                    [self.CPK_cdr_b])
                                    
            self.CPK_aft.plot3DmultiList(
                                    [self.CPK_aft,self.CPK_fwd],
                                    [])
                                    
            
            SD.hull.SAC.plot3DmultiList(
            SD.hull.lcurvenet[:1],
            [SD.hull.CProfile])
            
            
            
            
            self.CProfile.plot3DmultiList([self.CProfile,
                                     self.CPK_cdr_b],
                                     [self.bowfairning,
                                     self.CPK_fwd,
                                     self.sternfairing,
                                     self.CPK_aft,])
            
            
            self.CProfile.plot3DmultiList([self.CProfile,
                                         self.CPK_cdr_b],
                                         [self.CPK_fwd,
                                         self.CPK_aft,])
            
            OSV:
            self.CProfile.plot3DmultiList([self.CProfile,
                                         self.CPK_cdr_b],
                                         [self.CPK_fwd])
                                     
            self.CProfile2D.plotcurve_detailed()
            self.CPK_aft.plotcurve(color_='green')
            self.CPK_fwd.plotcurve(color_='green')
            Lspline.curve.plotcurve()
            
            self.CPK_cdr_b.plotcurve_detailed()
            self.CPK_aft.plotcurve(color_='green')
            self.CPK_fwd.plotcurve(color_='green')
        
        """
        mssg = '\n hull_from_simple_designspace : split_cLProfile \n'
        #if self._verbose:print mssg
        print mssg
        #
        #**********************************************************************
        # basics for this project
        k       = self.k
        nump    = self.nump
        #
        #**********************************************************************
        # 
        cpk = self.CProfile2D  # more efficient to use 2D 
                               # when reparameterizing w/ knot removal
        #                      
        # import SpecializedCurveSolvers as scsolvers
        #
        
        self.CPK_parametric_details = {}
        cLKP_solver = scsolvers.CPKsplitcurves(initial_curve=self.CProfile2D,
                                               LazyFPD=None,
                                               ordered_lazy_solvercontainer=None)
        #
        #**********************************************************************
        # 
        fp_index = False #no regard for types, do I seem to have.
        # in my ignorance, I apologize, Vladimir Voevodsky.
        # if this doesn't work, I promise to solve a knot vector
        # reparameterization problem via homotopic deformation of the knots.
        #
        if cpk.dim ==2:
            fp_index = 0
        elif cpk.dim ==3:
            fp_index = 0
        else:
            mssg = '\n ERROR: Fatal; Could not find dimensions '+\
                                        'of CPKeel curve to be split \n'
            assert(fp_index),mssg
        #
        #**********************************************************************
        # Which is outside the other: flat of Cprofile?, 
        #    or intersections with fairness curves?
        #longi_ini,longi_fini = self.stations.FOCP[:2] #these are vertex constraints
        # not curve point constarints!
        
        #
        #**********************************************************************
        # BOW
        #
        # Begin:  find location to split the Bow curve from the CProfile
        #
        #
        #**********************************************************************
        initial_param_val, q, found = cLKP_solver.find_xy_max_flat_fwd_extents(
                                                curve=self.CProfile2D,
                                                maxDraft = self.MaxDraft)
        mssg = '\n ERROR: Fatal; Could not find the fwd edge '+\
                                'of the flat of keel profile\n'
        assert(found),mssg
        bow_split_point = initial_param_val
       
        
        #
        #**********************************************************************
        #  STERN
        final_param_val, q, found = cLKP_solver.find_xy_max_flat_aft_extents(
                                                curve=self.CProfile2D,
                                                tol=1.e-15,
                                                maxDraft = self.MaxDraft)
        mssg = '\n ERROR: Fatal; Could not find the aft edge of the flat of keel profile\n'
        assert(found),mssg
        stern_split_point = final_param_val
        #        dmaxlist = []
        #        minflatloc = self.LengthDWL
        #        maxflatloc = self.CProfile2D.vertices[0,0]
        #        for vert in self.CProfile2D.vertices:
        #            if vert[1] == self.MaxDraft:
        #                minflatloc = min(minflatloc,vert[0])
        #                maxflatloc = max(maxflatloc,vert[0])
        #                dmaxlist.append(vert)
                
                
        #
        #
        #**********************************************************************
        # NOTE: fairness curves must intersect 
        # the longitudinal keel curve
        #
        #**********************************************************************
        # Bow 
        #   split point data gathering from Cprofile and stations
        #
        # stations along water line length:
        bstart,bend = self.stations.Bow_Fairing[2:]
        sstart,send = self.stations.Stern_Fairing[2:]
        
        # corresponding locations on the unsplit CProfile curve:
        bend = self.CProfile2D.FindPoint(self.stations.Bow_Fairing[1],0)[0]
        send = self.CProfile2D.FindPoint(self.stations.Stern_Fairing[1],0)[0]
        
        bend = bend - 0.25*bend
        send = send + 0.25*send
        #
        bow_split_point = min(bend,
                              bow_split_point)
        bow_split_point = max(.05,bow_split_point)
        
        stern_split_point = max(send,
                                stern_split_point)
        stern_split_point = min(.95,stern_split_point)
        
        final_stern_split_pt = self.CProfile2D.CurvePoint(stern_split_point)[0]
        
        
         
        self.CPK_parametric_details['bow_split_point'] = bow_split_point
        
        
        self.CPK_parametric_details['stern_split_point'] = stern_split_point
        
        
        #
        #**********************************************************************
        #
        #
        #**********************************************************************
        # find the location at which to cut the CProfile curve for
        #   the fwd portions of the longitudinal keel curve (and fore and possibly aft transverses)
        #
        #
        #**********************************************************************
        # cut the bow CProfile just fwd of the bow fairness curve
        #
        self.CPK_fwd, self.CPK_cdr_a = self.CProfile2D.CurveSplit(bow_split_point)
        self.CPK_OSV = self.CPK_cdr_a
        #
        #**********************************************************************
        #**********************************************************************
        #
        #
        # Bow split done.  Now Aft split 
        #
        #   (only needed if not making an OSV 
        #    but keeping this way at the moment) 
        #
        #
        #**********************************************************************
        #
        #**********************************************************************
        # CProfile aft split
        #
        # cut the stern CProfile just aft of the stern fairness curve
        #
        junk, self.CPK_aft = self.CProfile2D.CurveSplit(stern_split_point)
        #
        #**********************************************************************
        # CProfile flat section split
        """
        #now  find the fwd and aft flat edges to split on
        # and find the location of the fairness curves as they intersect the 
        #
        # NO:  The fairness curves only ever intersect the center keel longitudinal!
        # this is still true (even easier to see) with the OSV cLProfile
        # which 
        """
        a = self.CPK_cdr_a.FindPoint(final_stern_split_pt,0)[0]
        #
        #**********************************************************************
        # CProfile flat section split
        self.CPK_cdr_b, junk = self.CPK_cdr_a.CurveSplit(a)
        #self.CPK_cdr_b is the initial representation of the keel
        # but the knots are not yet correct.
        #
        #**********************************************************************
        #
        # Aft splitting Done
        #
        #**********************************************************************
        # 
        #
        #**********************************************************************
        #
        # Optimizing and fixing parameterization of the fwd cLProfile curve
        #
        # do this whether or not you are making an OSV
        #
        # 
        #**********************************************************************
        # 
        #**********************************************************************
        #
        # Special Functions:
        #
        # Let's give the bow a bulb shape to help integrate with the
        # fwd fairness curve (all to eventually help integrate with the Kraft bulb)
        #
        #
        """TODO: fix this!
        
        the bulb function is active at the start
        of the curve right now....
        corresponding to the part of the fwd profile curve
        near the dwl
        
        self = SD.hull
        """
        #        # ERROR: do note that this curve runs 'top to bottom'
        #        # when it begins life!
        #        #
        #        # special functions - circular bulb section of a curve
        #        circle = spfuncs.circle
        #        #ellipse = spfuncs.ellipse
        #        gfunc = spfuncs.gfunc
        #        #
        #        radius = 3.5 # assumed radius of bulbous bow
        #        #rb_ = 1.8   # elliptic foci
        #        #ra_ = 2.25  # elliptic foci
        #        #
        #        nump=30
        #        mp = .4 #
        #        pt1 = self.CPK_fwd.vertices[0]
        #        org = pt1[1] + radius
        #        c1 = circle(origin=[pt1[1],org],
        #                    radius=radius)
        #        #e1 = ellipse(origin=[0.,rb_],
        #        #             ra=ra_, 
        #        #             rb=rb_)
        #        circ_p = np.linspace(-np.pi/2., 
        #                             np.pi/2., 
        #                             nump, 
        #                             endpoint=False)
        #        def gfunc_closure(curve):
        #            """this little closure is necessary to wrap 
        #            the mx_plots.gfunc class 
        #            because we do not 
        #            have immediate access to the 
        #            Lagrangian here.
        #            -simply delay until we do have it!
        #            """
        #            return  gfunc(copy.deepcopy(curve), 
        #                             c1.circle_pt, 0., 
        #                             mp,circ_p, 10.)
        
        # 
        #**********************************************************************
        #  BOW
        #  FWD LSPLINE setup start (both transport and OSV get this one)
        # actually I think this one may be universal
        #
        # 
        # fwd goes from keel to dwl :-
        tan_fwd = np.rad2deg( self.CPK_fwd.compute_tangent(1.) )
        curv_fwd = np.rad2deg( self.CPK_fwd.compute_curvature(1.) )
        #
        CI = LazyFPDConstraintInterface()
        #
        #
        #**********************************************************************
        #  integral parameters
        self.CPK_fwd.compute_area()
        area = self.CPK_fwd.area
        CI.Constraint_Area(value=area,
                           kind='equality')
        
        self.CPK_fwd.compute_moments()
        self.CPK_fwd.computeXc()
        Xc = self.CPK_fwd.Xc
        CI.Constraint_Xc(value = Xc,
                         kind = 'equality')
        #
        #**********************************************************************
        # fixity conditions
        # needed for longitudinal flat portion!
        for i,vert in enumerate(self.CPK_fwd.vertices[1:-1]):
            this = i+1
            if np.linalg.norm(vert[1]-self.MaxDraft)<self.tol:
                #print 'setting fixity'
                CI.Constraint_Fix(index_tuple=[this,1],
                                  value=self.MaxDraft)
        #
        #**********************************************************************
        # end condition parameters
        CI.Constraint_Tangent(value = tan_fwd,
                              location=1.,
                              kind = 'equality')
        CI.Constraint_Curvature(value = curv_fwd,
                              location=1.,
                              kind = 'equality')
        #
        #**********************************************************************
        # solve
        print 'solving initial bow transverse (CP split curve)'
        Lspline = uopt.match_curve(
                        initial_curve = copy.deepcopy(self.CPK_fwd),
                        LazyFPD = CI,
                        return_Lspline=True,
                        make_thb = True,
                        ncpts=self.example_transverse.n,
                        special_func = None,#{'bulb':gfunc_closure},
                        nmatch=None,
                        wt=1)
        stop = uopt.ckLspline(Lspline)
        if not stop:
        
            #
            #**********************************************************************
            # FWD intermediate
            
            CI = LazyFPDConstraintInterface()
            self.CPK_fwd.compute_area()
            area = self.CPK_fwd.area
            CI.Constraint_Area(value=area,
                               kind='equality')
            
            self.CPK_fwd.compute_moments()
            self.CPK_fwd.computeXc()
            Xc = self.CPK_fwd.Xc
            CI.Constraint_Xc(value = Xc,
                             kind = 'LS')
            
            
            #
            #**********************************************************************
            # fixity conditions
            # needed for longitudinal flat portion!
            for i,vert in enumerate(self.CPK_fwd.vertices[3:-1]):
                this = i+3
                if np.linalg.norm(vert[1]-self.MaxDraft)<self.tol:
                    #print 'setting fixity'
                    CI.Constraint_Fix(index_tuple=[this,1],
                                      value=self.MaxDraft)
            #
            #**********************************************************************
            # end condition parameters
            CI.Constraint_Tangent(value = tan_fwd,
                                  location=1.,
                                  kind = 'equality')
            CI.Constraint_Curvature(value = curv_fwd,
                                  location=1.,
                                  kind = 'LS')
            
            #CPK_fwd.curve.plotcurve_detailed()
            print 'solving stage 2 bow transverse (CP split curve)'
            Lspline = uopt.match_curve(
                                    initial_curve = copy.deepcopy(self.CPK_fwd),
                                        LazyFPD = CI,
                                        return_Lspline=True,
                                        make_thb = True,
                                        ncpts=self.example_transverse.n,
                                        nmatch=None,
                                        #special_func = {'bulb':gfunc_closure},
                                        wt=10,)
                
        #
        #**********************************************************************
        # SAVE fwd cL Keel Profile (transverse)
        self.Lsplines.fwd_cpk = Lspline
        self.CPK_fwd = Lspline.curve
        
        
        
        #
        #**********************************************************************
        # Lift CPF_fwd (to be transverse curve), TODO: check parameterization
        self.CPK_fwd = self.lifting_2DCProfile(self.CPK_fwd.vertices,
                                                elevation = self.MaxDraft,
                                                k=self.CPK_fwd.k,
                                                nump=self.CPK_fwd.nump)
        #
        #**********************************************************************
        # CPK_fwd came from a "0 to lwl", aka "front to back" hull curve
        # therefor it must be reversed!
        #self.CPK_fwd = spline.Bspline(self.CPK_fwd.vertices[::-1], k,nump)
        self.CPK_fwd = spline.reverse_spline(self.CPK_fwd)
        
        
        #
        #******************************************************************
        # tcurve map (hyperparameter control, 
        #                   possible graphical access to Lspline)
        self.tcurve_map[0] = {'Lspline':None,
                               'tcurve':self.CPK_fwd}
        
        #
        # FWD END
        #**********************************************************************
        #
        if self.hullformtype == 'osv':
            print ''
            print ' CLProfile:  OSV TYPE'
            print 'cLProfile split'
            print 'only splitting fwd portion'
            print 'now fixing parameterizatin of aft portion'
            def fixit(curve):
                CI = LazyFPDConstraintInterface() 
                #
                #**********************************************************************
                # AREA from the curve to the x-axis
                curve.compute_area()
                area = curve.area
                CI.Constraint_Area(value=area,
                                   kind='equality')
                #
                #**********************************************************************
                # X centroid
                curve.compute_moments()
                curve.computeXc()
                Xc = curve.Xc
                CI.Constraint_Xc(value = Xc,
                                 kind = 'equality')
                #
                #**********************************************************************
                # fixity conditions
                # needed for longitudinal flat portion!
                for i,vert in enumerate(curve.vertices[1:-1]):
                    this = i+1
                    if np.linalg.norm(vert[1]-self.MaxDraft)<self.tol:
                        print 'fixing ',this
                        CI.Constraint_Fix([this,1],self.MaxDraft)
                #
                #**********************************************************************
                # locations of primary transverse joins
                #
                # inputs to the tan and curvature constraints
                # are in degrees
                #
                tan_fwd = np.rad2deg( curve.compute_tangent(0.) )
                curv_fwd = np.rad2deg( curve.compute_curvature(0.) )
                #
                tan_aft = np.rad2deg( curve.compute_tangent(1.) )
                curv_aft = np.rad2deg( curve.compute_curvature(1.) )
                
                
                CI.Constraint_Tangent(value = tan_fwd,
                                      location=0.,
                                      kind = 'equality')
                CI.Constraint_Curvature(value = curv_fwd,
                                      location=0.,
                                      kind = 'equality')
                
                CI.Constraint_Tangent(value = tan_aft,
                                      location=1.,
                                      kind = 'equality')
                CI.Constraint_Curvature(value = curv_aft,
                                      location=1.,
                                      kind = 'equality')
                #
                #**********************************************************************
                # solve for the reparameterization, optimization method:
                #
                print 'solving initial aft transverse (CP split curve)'
                Lspline = uopt.match_curve(initial_curve = copy.deepcopy(curve),
                                            LazyFPD = CI,
                                            return_Lspline=True,
                                            make_thb = True,
                                            ncpts=self.CProfile2D.n,
                                            nmatch=None,
                                            wt=10,)
                #
                #**********************************************************************
                # DONE
                return Lspline
            Lspline = fixit(copy.deepcopy(self.CPK_OSV))
            self.Lsplines.OSV_keel = Lspline
            matched_curve = Lspline.curve
            #
            #**********************************************************************
            #  Lift the curve back to 3D (if it was indeed necessary to go 2D)
            matched_curve = self.lifting_2DCProfile(matched_curve.vertices,
                                                    elevation = self.MaxDraft,
                                                    k=matched_curve.k,
                                                    nump=matched_curve.nump)
            matched_curve = uopt.make_thb(matched_curve)
            self.CPK_cdr_b = matched_curve  #this is the new CPKeel after splitting
            self.CPK_OSV = matched_curve
            a = self.CPK_cdr_b.FindPoint(final_stern_split_pt, 2)
            self.CPK_parametric_details['stern_fairness_CLK_FindPoint'] = a
            print 'DONE WITH OSV Center Profile, OSV hull type'
        else:
            print 'cLProfile split'
            print 'making traditional hull'
            print 'splitting aft portion of cLProfile'
            #
            print 'now fixing parameterization of the longitudinal portion'
            
            
            def fixit(curve):
                CI = LazyFPDConstraintInterface() 
                #
                #**********************************************************************
                # AREA from the curve to the x-axis
                curve.compute_area()
                area = curve.area
                CI.Constraint_Area(value=area,
                                   kind='equality')
                #
                #**********************************************************************
                # X centroid
                curve.compute_moments()
                curve.computeXc()
                Xc = curve.Xc
                CI.Constraint_Xc(value = Xc,
                                 kind = 'equality')
                #
                #**********************************************************************
                # fixity conditions
                # needed for longitudinal flat portion!
                for i,vert in enumerate(curve.vertices[1:-1]):
                    this = i+1
                    if np.linalg.norm(vert[1]-self.MaxDraft)<self.tol:
                        print 'fixing ',this
                        CI.Constraint_Fix([this,1],self.MaxDraft)
                #
                #**********************************************************************
                # locations of primary transverse joins
                #
                # inputs to the tan and curvature constraints
                # are in degrees
                #
                tan_fwd = np.rad2deg( curve.compute_tangent(0.) )
                curv_fwd = np.rad2deg( curve.compute_curvature(0.) )
                #
                tan_aft = np.rad2deg( curve.compute_tangent(1.) )
                curv_aft = np.rad2deg( curve.compute_curvature(1.) )
                
                
                CI.Constraint_Tangent(value = tan_fwd,
                                      location=0.,
                                      kind = 'equality')
                CI.Constraint_Curvature(value = curv_fwd,
                                      location=0.,
                                      kind = 'equality')
                
                CI.Constraint_Tangent(value = tan_aft,
                                      location=1.,
                                      kind = 'equality')
                CI.Constraint_Curvature(value = curv_aft,
                                      location=1.,
                                      kind = 'equality')
                #
                #**********************************************************************
                # solve for the reparameterization, optimization method:
                #
                Lspline = uopt.match_curve(initial_curve = copy.deepcopy(curve),
                                            LazyFPD = CI,
                                            return_Lspline=True,
                                            make_thb = True,
                                            ncpts=self.CProfile2D.n,
                                            nmatch=None,
                                            wt=10,)
                #
                #**********************************************************************
                # DONE
                return Lspline
            #
            #**********************************************************************
            #  longitudinal cLProfile returned from reparameterization:
            Lspline = fixit(copy.deepcopy(self.CPK_cdr_b))
            self.Lsplines.long_keel = Lspline
            matched_curve = Lspline.curve
            #
            #**********************************************************************
            #  Lift the curve back to 3D (if it was indeed necessary to go 2D)
            matched_curve = self.lifting_2DCProfile(matched_curve.vertices,
                                                    elevation = self.MaxDraft,
                                                    k=matched_curve.k,
                                                    nump=matched_curve.nump)
            matched_curve = uopt.make_thb(matched_curve)
            self.CPK_cdr_b = matched_curve  #this is the new CPKeel after splitting
            #a = self.CPK_cdr_b.FindPoint(loc_x,2) 
            a = self.CPK_cdr_b.FindPoint(final_stern_split_pt, 2)
            self.CPK_parametric_details['stern_fairness_CLK_FindPoint'] = a
            
            
            print 'now fixing parameterizatin of aft portion'
            #
            #**********************************************************************
            #**********************************************************************
            #
            # 
            #
            #**********************************************************************
            #**********************************************************************
            # start the constraints for this curve, 
            # but delay computation a bit ;) 
            # (because we don't have a curve with the correct knot vector 
            # with which to instantiate the FPD class)
            CI = LazyFPDConstraintInterface() #lazy FPD factory simple magic
            #
            #**********************************************************************
            # AREA from the curve to the x-axis
            self.CPK_cdr_b.compute_area()
            area = self.CPK_cdr_b.area
            CI.Constraint_Area(value=area,
                               kind='equality')
            #
            #**********************************************************************
            # X centroid
            self.CPK_cdr_b.compute_moments()
            self.CPK_cdr_b.computeXc()
            Xc = self.CPK_cdr_b.Xc
            CI.Constraint_Xc(value = Xc,
                             kind = 'equality')
            #
            #**********************************************************************
            # fixity conditions
            # needed for longitudinal flat portion!
            for i,vert in enumerate(self.CPK_cdr_b.vertices[1:-1]):
                this = i+1
                if np.linalg.norm(vert[1]-self.MaxDraft)<self.tol:
                    print 'fixing ',this
                    CI.Constraint_Fix([this,1],self.MaxDraft)
            #
            #**********************************************************************
            # joint center section with fore and aft transverses (if they are both split)
            # (stations.FOSAC, stations.FOCP, stern and bow Fairness Curves)
            #
            # always remember - inputs to the tan and curvature constraints
            # are in degrees!
            #
            tan_fwd = np.rad2deg( self.CPK_cdr_b.compute_tangent(0.) )
            curv_fwd = np.rad2deg( self.CPK_cdr_b.compute_curvature(0.) )
            #
            tan_aft = np.rad2deg( self.CPK_cdr_b.compute_tangent(1.) )
            curv_aft = np.rad2deg( self.CPK_cdr_b.compute_curvature(1.) )
            
            
            CI.Constraint_Tangent(value = tan_fwd,
                                  location=0.,
                                  kind = 'equality')
            CI.Constraint_Curvature(value = curv_fwd,
                                  location=0.,
                                  kind = 'equality')
            
            CI.Constraint_Tangent(value = tan_aft,
                                  location=1.,
                                  kind = 'equality')
            CI.Constraint_Curvature(value = curv_aft,
                                  location=1.,
                                  kind = 'equality')
            
            
            #
            #**********************************************************************
            # solve for the reparameterization, optimization method:
            #
            # note: FPD.curve <=> initial_curve
            #
            Lspline = uopt.match_curve(initial_curve = copy.deepcopy(self.CPK_cdr_b),
                                            LazyFPD = CI,
                                            return_Lspline=True,
                                            make_thb = True,
                                            ncpts=self.CProfile.n,
                                            nmatch=None,
                                            wt=10,)
            self.Lsplines.long_keel = Lspline
            matched_curve = Lspline.curve
            #
            #**********************************************************************
            #  Lift the curve back to 3D (if it was indeed necessary to go 2D)
            matched_curve = self.lifting_2DCProfile(matched_curve.vertices,
                                                    elevation = self.MaxDraft,
                                                    k=matched_curve.k,
                                                    nump=matched_curve.nump)
            matched_curve = uopt.make_thb(matched_curve)
            self.CPK_cdr_b = matched_curve  #this is the new CPKeel after splitting
            #a = self.CPK_cdr_b.FindPoint(loc_x,2) 
            a = self.CPK_cdr_b.FindPoint(final_stern_split_pt, 2)
            self.CPK_parametric_details['stern_fairness_CLK_FindPoint'] = a
            #
            #**********************************************************************
            #(then split the cLProfile curve aft!)
            #
            #**********************************************************************
            #  AFT more FPD
            # could add a stern bulb to the profile here
            # using the Lagrangian special functions
            
            #
            #**********************************************************************
            # 
            # 
            #  Repaameterize : Aft cL Keel Profile 
            #
            #**********************************************************************
            #  AFT start
            
            tan0 = np.rad2deg( self.CPK_aft.compute_tangent(0.) )
            curv0 = np.rad2deg( self.CPK_aft.compute_curvature(0.) )
            
            
            CI = LazyFPDConstraintInterface()
            
            self.CPK_aft.compute_area()
            area = self.CPK_aft.area
            CI.Constraint_Area(value=area,
                               kind='equality')
            
            self.CPK_aft.compute_moments()
            self.CPK_aft.computeXc()
            Xc = self.CPK_aft.Xc
            CI.Constraint_Xc(value = Xc,
                             kind = 'equality',
                             weight=100.)
            
            #self.CPK_aft.pts_M_pts()
            #self.CPK_aft.compute_arclength()
            #CI.Constraint_ArcLength(self.CPK_aft.AL,
            #                        kind = 'equality)
            #CI.Constraint_ArcLengthApprox(value=self.CPK_aft.AL,
            #                              kind = 'equality')
            
            
            #
            #**********************************************************************
            # fixity conditions
            # needed for longitudinal flat portion!
            for i,vert in enumerate(self.CPK_aft.vertices[1:-1]):
                this = i+1
                if np.linalg.norm(vert[1]-self.MaxDraft)<self.tol:
                    print 'setting fixity, ',this
                    CI.Constraint_Fix([this,1],self.MaxDraft)
            #
            
                
            CI.Constraint_Tangent(value = tan0,
                                  location=0.,
                                  kind = 'equality')
            CI.Constraint_Curvature(value = curv0,
                                  location=0.,
                                  kind = 'equality',
                                  weight=100.)
            
            
            Lspline = uopt.match_curve(
                                    initial_curve = copy.deepcopy(self.CPK_aft),
                                        LazyFPD = CI,
                                        return_Lspline=True,
                                        make_thb = True,
                                        ncpts=self.example_transverse.n,
                                        nmatch=None,
                                        wt=10,)
            
            stop = uopt.ckLspline(Lspline)
            if not stop:
                #
                #**********************************************************************
                #  AFT more FPD
                CI = LazyFPDConstraintInterface()
                
                self.CPK_aft.compute_area()
                area = self.CPK_aft.area
                CI.Constraint_Area(value=area,
                                   kind='equality')
                
                self.CPK_aft.compute_moments()
                self.CPK_aft.computeXc()
                Xc = self.CPK_aft.Xc
                CI.Constraint_Xc(value = Xc,
                                 kind = 'LS')
                
                
                
                #
                #**********************************************************************
                # fixity conditions
                # needed for longitudinal flat portion!
                for i,vert in enumerate(self.CPK_aft.vertices[1:3]):
                    this = i+1
                    if np.linalg.norm(vert[1]-self.MaxDraft)<self.tol:
                        print 'setting fixity, ',this
                        CI.Constraint_Fix([this,1],self.MaxDraft)
                #
            
                
                CI.Constraint_Tangent(value = tan0,
                                      location=0.,
                                      kind = 'equality')
                CI.Constraint_Curvature(value = curv0,
                                      location=0.,
                                      kind = 'LS')
                
                
                Lspline = uopt.match_curve(
                                        initial_curve = copy.deepcopy(self.CPK_aft),
                                            LazyFPD = CI,
                                            return_Lspline=True,
                                            make_thb = True,
                                            ncpts=self.example_transverse.n,
                                            nmatch=None,
                                            wt=1,)
                
            
            
            
            stop = uopt.ckLspline(Lspline)
            if not stop:
                #
                #**********************************************************************
                #  AFT Max FPD
                
                CI = LazyFPDConstraintInterface()
                
                self.CPK_aft.compute_area()
                area = self.CPK_aft.area
                CI.Constraint_Area(value=area,
                                   kind='equality')
                
                self.CPK_aft.compute_moments()
                self.CPK_aft.computeXc()
                Xc = self.CPK_aft.Xc
                CI.Constraint_Xc(value = Xc,
                                 kind = 'LS')
                
                
                #
                #**********************************************************************
                # fixity conditions
                # needed for longitudinal flat portion!
                for i,vert in enumerate(self.CPK_aft.vertices[1:3]):
                    this = i+1
                    if np.linalg.norm(vert[1]-self.MaxDraft)<self.tol:
                        CI.Constraint_Fix([this,1],self.MaxDraft)
                #
            
                
                
                #        CI.Constraint_Tangent(value = 0.,
                #                              location=0.,
                #                              kind = 'equality',
                #                              weight=100.)
                #        CI.Constraint_Curvature(value = 0.,
                #                              location=0.,
                #                              kind = 'LS',
                #                              weight=.01)
                
                CI.Constraint_Tangent(value = tan0,
                                      location=0.,
                                      kind = 'equality')
                CI.Constraint_Curvature(value = curv0,
                                      location=0.,
                                      kind = 'LS')
                
                
                Lspline = uopt.match_curve(
                                        initial_curve = copy.deepcopy(self.CPK_aft),
                                            LazyFPD = CI,
                                            return_Lspline=True,
                                            make_thb = True,
                                            ncpts=self.example_transverse.n,
                                            nmatch=None,
                                            wt=10,)
            
            
            stop = uopt.ckLspline(Lspline)
            if not stop:
            
            
                #
                #**********************************************************************
                #  AFT fix
                
                
                CI = LazyFPDConstraintInterface()
                
                self.CPK_aft.compute_area()
                area = self.CPK_aft.area
                CI.Constraint_Area(value=area,
                                   kind='equality')
                
                self.CPK_aft.compute_moments()
                self.CPK_aft.computeXc()
                Xc = self.CPK_aft.Xc
                CI.Constraint_Xc(value = Xc,
                                 kind = 'LS')
                
                
                
                #
                #**********************************************************************
                # fixity conditions
                # needed for longitudinal flat portion!
                for i,vert in enumerate(self.CPK_aft.vertices[1:3]):
                    this = i+1
                    if np.linalg.norm(vert[1]-self.MaxDraft)<self.tol:
                        CI.Constraint_Fix([this,1],self.MaxDraft)
                #
            
                
                CI.Constraint_Tangent(value = tan0,
                                      location=0.,
                                      kind = 'LS')
                CI.Constraint_Tangent(value = curv0,
                                      location=0.,
                                      kind = 'LS')
                
                Lspline = uopt.match_curve(
                                        initial_curve = copy.deepcopy(self.CPK_aft),
                                            LazyFPD = CI,
                                            return_Lspline=True,
                                            make_thb = True,
                                            ncpts=self.example_transverse.n,
                                            nmatch=None,
                                            wt=10,)
                
            self.Lsplines.aft_cpk = Lspline
            self.CPK_aft = Lspline.curve
            #
            #**********************************************************************
            # Lift CPK_aft (to be transverse curve), TODO: check parameterization
            self.CPK_aft = self.lifting_2DCProfile(self.CPK_aft.vertices,
                                                    elevation = self.MaxDraft,
                                                    k=self.CPK_aft.k,
                                                    nump=self.CPK_aft.nump)
            print '\n\n Done with Center Profile, Traditional Stern hull type\n\n'
            #
            #******************************************************************
            #
        
        #
        #******************************************************************
        # tcurve map (hyperparameter control, 
        #                   possible graphical access to Lspline)
        self.tcurve_map[13] = {'Lspline':None,
                               'tcurve':self.CPK_aft}
            
        #
        #**********************************************************************
        # old double check
        if self.hullformtype != 'osv':
            if self.CPK_fwd.n==self.CPK_aft.n:
                pass
            else:
                print "CPK Split Curves Do not match!"
        #
        #**********************************************************************
        # Done
        print 'Done with Center Profile'
        return
    
    
    
    
    
    
    
    def split_DWL_curves(self):
        """
        TODO:  make the drop in the cLProfile curve
        at transom (with associated area (volume) 
        carried into the aft edge)
        incorporated into the SAC.
        
        
        assumptions
        ----------
            -we want an OSV-like set of curves
            
            self.DWL_OSV, self.DWL_aft
            
        therefore
        ----------
            -we do not actually need to split the DWL then.
            
            
        notes
        ----------
            -the DWL now ends wide - with breadth similar to that
            of midship (need to make it a design parameter, by the way)
            and so we do not split the DWL,
            but instead just tack on the aft most transverse
        
            
            import curve             as     spline
            from hull_from_simple_designspace import linear_vertices
        
            import utility_optimization as uopt
            from   adials_gui3d      import frame, DrawCurveInteractive, hull_frame
        
        """
        #self.DWL_parametric_details = {}
        #self.DWL_parametric_details['bow_split_point'] = 0. #no splitting
        #self.DWL_parametric_details['stern_split_point'] = 1.0 #no splitting
        #
        #**********************************************************************
        #
        k = self.k
        nump = self.nump
        #
        #**********************************************************************
        #  Make the DWL 3D longitudinal thb-spline
        self.DWL_OSV = copy.deepcopy(self.DWL2D)
        
        verts = copy.deepcopy(self.DWL_OSV.vertices)
        #
        #**********************************************************************
        #  lift to 3D
        vnew = frame.xy_to_zxy_elevate(verts,
                                       elevation=self.MaxDraft)
        #
        #**********************************************************************
        # SAVE
        self.DWL_OSV = spline.Bspline(vnew, k, nump)
        self.DWL_OSV = uopt.make_thb(self.DWL_OSV)
        print 'Done with primary DWL generation, OSV hull type  \n\n'
        #
        #**********************************************************************
        # 
        
        
        #
        #**********************************************************************
        #  Make the transverse transom thb-spline
        #
        start = self.CProfile.CurvePoint(1.)
        end = self.DWL_OSV.CurvePoint(1.)
        yintercept = start[1] + (end-start)[1]/2.
        intersect = np.asarray([ end[0], yintercept, end[2] ])
        #
        #        verts = linear_vertices(tuple(start),
        #                                tuple(end),
        #                                self.n)
        #
        v1 = linear_vertices(start,
                             intersect,
                             5)
        v2 =  linear_vertices(intersect,
                             end,
                             3)
        #vertices1 = list(v1)+list(v2[0])+list(v2)
        verts = list(v1[:4])+list(v2)
        verts = np.asarray(verts)
        #
        #**********************************************************************
        # 
        #
        #**********************************************************************
        # SAVE (transverse must run from cLProfile (0) to end at DWL (1.))
        print 'new stuff'
        print 'verts = ',
        print verts
        print verts[:,:2]
        curve = spline.Bspline(verts[:,:2], k, nump)
        
        #
        #**********************************************************************
        # TBD: make it match SAC aftmost.
        #
        #**********************************************************************
        # 
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve)
        Aarea = self.SAC(1.)[1]
        FPD.add_Area_to_y_Constraint(kind='equality',
                                         value = Aarea)
        FPD.add_xFixity(index=5,
                        value=FPD.curve.vertices[-1][0])
        FPD.add_xFixity(index=4,
                        value=FPD.curve.vertices[-1][0])
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5) 
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
        #
        vts = []
        for vert,v3 in zip(Lspline.curve.vertices,verts[:,2:]):
            lv = list(vert)
            lv.append(v3)
            vts.append(lv)
        self.DWL_aft = spline.Bspline(np.asarray(vts), k, nump)
        #
        #
        #******************************************************************
        # tcurve map (hyperparameter control, 
        #                   possible graphical access to Lspline)
        self.tcurve_map[13] = {'Lspline':None,
                               'tcurve':self.DWL_aft}
        print 'Done with aft transverse generation, OSV hull type \n\n'
        return
    
    
    
    
    
    
    
    


    def correct_hull_curves(self):
        if self._verbose: print 'hull_from_simple_designspace: correct_hull_curves'
        
        diff = self.CProfile.n-self.DWL.n
        assert(diff==0),"Complication, DWL.n != CProfile.n"
        #        if not (diff==0):
        #            #if self.CProfile.n>self.DWL.n:            
        #            assert(self.CProfile.n==self.DWL.n),'self.DWL /= self.CProfile'
        #                #self.DWL = self.insert_ctrl_vertices(self.DWL,diff)
        #            #elif self.CProfile.n<self.DWL.n:
        #                #self.CProfile = self.insert_ctrl_vertices(self.CProfile,diff)
        assert(diff==0),"Error, DWL.n => CProfile.n correction failed!"
        return
    

    def make_transverse_net(self):
                            #,dropout = None):
        """
        
        Set up the transverse net with the picture
        that longitudinals are defined starting fore and running aft

        
        input
        ----------
            dropout : 11 transverse curves have been solved for
                -this is, heuristically, a lot of definition for the global surface
                -as informed by looking at various 'nice designs' from industry software
                examples.
                -more to the piont,
                this is totally constraining to the 11 control vertex longitudinals
                that we have on hand.
                -dropout designates one transverse curve to be forgotten
                so as to free up 1 R**3 vertex for fairinig
                per longitudinal curve.
                
            idea : do dropout at the last possible second, 
                    encourageing smoother curves which coform 
                    more to 'what we want"?
                ergo, do not dropout here.
                
                dropout 1 would drop the first transverse aft of the nose transverse.
        
        Transverse Curve Network,
        curve definitions by index:
        0: nose curve (do not drop!)
        1: bowfairning
        2: bowtransition
        
        [3-7] it varies!
        3: mshp_sections[0] (flat of side)
        4: mshp_sections[0] (flat of side)
        5: mshp_sections[0] (midship for real)   # (flat of side)
        6: mshp_sections[0] (flat of side)       #   (midship for real)
        7: mshp_sections[0] (flat of side)       #   (midship for real)
        
        8: sterntransition
        9: sternfairing
        10: aftcurve
        
        """
        if self._verbose: print 'hull_from_simple_designspace: make_transverse_net'
        self.lfixities = [[],[],[]]
        
        if self.hullformtype == 'osv':
            fwdtcurve = self.CPK_fwd
            afttcurve = self.DWL_aft    #vestigal name - DWL is now one piece, 
                                        #non split in OSV case
        else: #trad hull - CProfile was split fore and aft to make transverse ends
            fwdtcurve = self.CPK_fwd
            afttcurve = self.CPK_aft
        
        
        #
        #**********************************************************************
        # 
        
        #        if self.heavy_config:
        #            self.tcurvenet = [afttcurve,
        #                              self.sterntransition]
        #        else:
        
        # aft
        self.tcurvenet = [afttcurve,
                          self.sternfairing,
                          self.sterntransition]
            
        #from left to right, first clip one curve off the net (end is aft), 
        # then flip the order to start aft
        #        if self.remove_aft_fos:
        #            for c in self.mshp_sections[:-1][::-1]: 
        #                self.tcurvenet.append(c)
        #        else: #do not clip any curves out
        
        
        #FOS fwd + FOSAC + FOS aft
        for c in self.mshp_sections[::-1]: #loop mship aft to fore
            self.tcurvenet.append(c)
            self.lfixities[0].append(c.vertices)
        
        #        if self.remove_fwdfairness:
        #            self.tcurvenet.append(self.bowtransition)
        #            self.tcurvenet.append(fwdtcurve)
        #        else:
        
        #bow
        self.tcurvenet.append(self.bowtransition)
        self.tcurvenet.append(self.bowfairning)
        self.tcurvenet.append(fwdtcurve)
        #
        #**********************************************************************
        # 
        
        
    
        
        #   self.tcurvenet
        #   is built aft to front
        self.tcurvenet.reverse() 
        self.lfixities[0].reverse()
        #
        #**********************************************************************
        # dropout: eliminate extra curve from the system
        # idea: leave this in at the moment,
        # this might help
        # put the DOF where they are needed most??
        #if dropout is not None:
        #    self.dropped_transverse = self.tcurvenet.pop(dropout)
        
        
        #
        #**********************************************************************
        # 
        #
        #*****************************
        # to use heavy_config, you must
        # have asked the FPD routines to build in 1 extra curve
        #        if self.heavy_config:
        #            tvertnet = [afttcurve.vertices,
        #                        self.sterntransition.vertices]
        #        else:
        
        
        #aft
        tvertnet = [afttcurve.vertices,
                    self.sternfairing.vertices,
                    self.sterntransition.vertices]
        
        
        #
        #*****************************
        # to use remove_aft_fos, you must
        # have asked the FPD routines to build in 1 extra curve
#        if self.remove_aft_fos:
#            for c in self.mshp_sections[:-1][::-1]: #first truncate the list, then reverse
#                tvertnet.append(c.vertices)
#        else:
        
        #FOS fwd + FOSAC + FOS aft
        for c in self.mshp_sections[::-1]:
            tvertnet.append(c.vertices)
        
        #        if self.remove_fwdfairness:
        #            tvertnet.append(self.bowtransition.vertices)
        #            tvertnet.append(fwdtcurve.vertices)
        #        else:
        
        #bow
        tvertnet.append(self.bowtransition.vertices)
        tvertnet.append(self.bowfairning.vertices)
        tvertnet.append(fwdtcurve.vertices)
        
        
        #
        #**********************************************************************
        # 
                     
        tvertnet.reverse() 
        #
        #**********************************************************************
        # dropout: eliminate extra curve from the system
        # idea: leave this in at the moment,
        # this might help
        # put the DOF where they are needed most??
        #if dropout is not None:
        #    hide_me = tvertnet.pop(dropout) #yay garbage collectors ;)
        tvertnet = np.asarray(tvertnet)
        
        
        #now select from all curves and vertices
        # those curves and their vertices to use in the 
        # linear solver (a)
        # nonlinear solver (b)
        a = self.use_in_linear_interp
        b = self.use_in_longi_solver
        
        m1 = 'ERROR: Must have {} transverse curves in the linear longitudinal solver. n= {}'
        m2 = 'ERROR: Must have {} transverse curves in the nonlinear longitudinal solver. n= {}'
        this = len(a)
        that = len(b) - len(self.dropout) - len(self.ignore)
        #self.check_nvertices(11,this,m1)
        #self.check_nvertices(11,that,m2,op='<=')
        
        self.solve_linear_list_curves = [ self.tcurvenet[i] for i in a]
        self.solve_linear_list_vertices = np.asarray([ tvertnet[i] for i in a])
        self.solve_nonlinear_list_curves = [ self.tcurvenet[i] for i in b]
        self.solve_nonlinear_list_vertices = np.asarray([ tvertnet[i] for i in b])
        ##return##
        return tvertnet
    
    
    def issue_longitudinals(self, 
                            refine=True,
                            extrema_enforcement = None,
                            dropout = None):
        """
        TODO: make dropout use a list instead of a single index
        dropout somewhere with less shape (presently dropping the
        first traverse aft of the nose!)
        
        issue_longitudinals is part of the overall 
            set of functions called by the function make_longitudinals
            
            
            extrema_enforcement : 
            
            note: this will compute the transverse curve net as well.
            
        inputs
        ----------
            refine: use optimization to hit vertex targets.  
                        straighten out the curves
                        but hit bow end conditions for flattness
                        and parameterize nicely
            
            extrema_enforcement: object to provide constraints to the solver
                                for
                                enforcing bare hull extremes 
                                as inequality constraints
                                
                                from ShipDesigner class
                                this is the extrema found
                                after a thin hull design is generated
                                    
            dropout: 11 transverse curves have been solved for
                        -this is, heuristically, a lot of definition for the global surface
                        -as informed by looking at various 'nice designs' from industry software
                        examples.
                        -more to the piont,
                        this is totally constraining to the 11 control vertex longitudinals
                        that we have on hand.
                        -dropout designates one transverse curve to be forgotten
                        so as to free up 1 R**3 vertex for fairinig
                        per longitudinal curve.
                        
                    idea - do dropout at the last possible second, 
                            encourageing smoother curves which coform 
                            more to 'what we want"?
                        ergo, dropout in the optimization algorithm itself
                        
                    Important definition:
                        dropout 1 would drop the first transverse aft of the nose transverse.
                
                
        notes
        ----------
            nose constraint : moved to THB as it is a local constraint
            
            
            dropout, which curve am I dropping?:
            ---
            dropout drops from the actual existing list of
            tcurves that make it to the longitudinal solver
            -- 
            this is not great because 
            what the index refers to depends on 
            what has been put into the list
            ---
            the total number of tcurves must be 11
            but which 11? - that depends
            ---
            here is a possible list:
            0: nose curve (do not drop!)
            1: bowfairning
            2: bowtransition
            3: mshp_sections[0] (flat of side)
            4: mshp_sections[0] (flat of side)
            5: mshp_sections[0] (midship for real)   #?# (flat of side)
            6: mshp_sections[0] (flat of side)       #?# (midship for real)
            7: mshp_sections[0] (flat of side)       # (midship for real)
            8: sterntransition
            9: sternfairing
            10: aftcurve
            
        dev:
        -----------
            self = SD.hull
            
            import curve as spline
            import utility_optimization as uopt
        """
        if self._verbose: print 'hull_from_simple_designspace: issue_longitudinals'
        k=self.k
        nump = self.nump
        if self.hullformtype != 'osv':
            N = self.CPK_aft.n
        else:
            N = self.DWL_aft.n
        dropout = self.dropout
        nnew = N-2
        #m = self.CProfile.n
        #lvertices = np.zeros((nnew,m),float)
        tvertnet = self.make_transverse_net()
        lcurvenet = []
        #N = self.tcurvenet[0].n
        #ptslist = []
        #n=7:
        #for i in range(N):
        # n=11:
        self.save_junk = []
        
        
        
        self.dropped_transverse = []
        #
        totl = self.dropout + self.ignore 
        dlist = self.use_in_longi_solver
        dlist = [ el for el in dlist if el not in totl]
        d2list = [ self.tcurvenet[i] for i in dlist]
        self.tcurvenet = d2list
        
        #**********************************************************************
        #**********************************************************************
        #
        #
        #           uopt.issue_longitudinals:
        #
        lc = uopt.issue_longitudinals(self.tcurvenet,
                                      k,nump,
                                      refine1st=True)
        #
        #
        #**********************************************************************
        #**********************************************************************
        lcurvenet3 = [Lspline.curve.get_coarsest_curve() for Lspline in lc]
        self.lcurvenet = lcurvenet3
        self.Lsplines.lcurves = [Lspline for Lspline in lc]
        
        
        #    for el in totl[::-1]:
        #        if el in self.dropout:
        #            self.dropped_transverse.append(el)
        #            
        #                elif el in self.ignore:
        #                    self.ignored_transverse.append(el)
        #                    self.tcurvenet.remove(el)
        #                elif el not in self.use_in_longi_solver:
        #                    self.not_in_longi.append(el)
        #                    self.tcurvenet.remove(el)
        #self.tcurvenet = self.solve_nonlinear_list_curves
        print '------------------------------------------'
        print 'correction finished'
        print 'Longitudinal Bare Hull Curve Solving is Done'
        print '------------------------------------------'
        return

    def issue_longitudinals_old(self, refine=True,
                            extrema_enforcement=None,
                            dropout = None):
        """
        TODO: make dropout use a list instead of a single index
        dropout somewhere with less shape (presently dropping the
        first traverse aft of the nose!)
        
        issue_longitudinals is part of the overall 
            set of functions called by the function make_longitudinals
            
            
            extrema_enforcement : 
            
            note: this will compute the transverse curve net as well.
            
        inputs
        ----------
            refine: use optimization to hit vertex targets.  
                        straighten out the curves
                        but hit bow end conditions for flattness
                        and parameterize nicely
            
            extrema_enforcement: object to provide constraints to the solver
                                for
                                enforcing bare hull extremes 
                                as inequality constraints
                                
                                from ShipDesigner class
                                this is the extrema found
                                after a thin hull design is generated
                                    
            dropout: 11 transverse curves have been solved for
                        -this is, heuristically, a lot of definition for the global surface
                        -as informed by looking at various 'nice designs' from industry software
                        examples.
                        -more to the piont,
                        this is totally constraining to the 11 control vertex longitudinals
                        that we have on hand.
                        -dropout designates one transverse curve to be forgotten
                        so as to free up 1 R**3 vertex for fairinig
                        per longitudinal curve.
                        
                    idea - do dropout at the last possible second, 
                            encourageing smoother curves which coform 
                            more to 'what we want"?
                        ergo, dropout in the optimization algorithm itself
                        
                    Important definition:
                        dropout 1 would drop the first transverse aft of the nose transverse.
                
                
        notes
        ----------
            nose constraint : moved to THB as it is a local constraint
            
            
            dropout, which curve am I dropping?:
            ---
            dropout drops from the actual existing list of
            tcurves that make it to the longitudinal solver
            -- 
            this is not great because 
            what the index refers to depends on 
            what has been put into the list
            ---
            the total number of tcurves must be 11
            but which 11? - that depends
            ---
            here is a possible list:
            0: nose curve (do not drop!)
            1: bowfairning
            2: bowtransition
            3: mshp_sections[0] (flat of side)
            4: mshp_sections[0] (flat of side)
            5: mshp_sections[0] (midship for real)   #?# (flat of side)
            6: mshp_sections[0] (flat of side)       #?# (midship for real)
            7: mshp_sections[0] (flat of side)       # (midship for real)
            8: sterntransition
            9: sternfairing
            10: aftcurve
            
        dev:
        -----------
            self = SD.hull
            
            import curve as spline
            import utility_optimization as uopt
        """
        if self._verbose: print 'hull_from_simple_designspace: issue_longitudinals'
        k=self.k
        nump = self.nump
        if self.hullformtype != 'osv':
            N = self.CPK_aft.n
        else:
            N = self.DWL_aft.n
        dropout = self.dropout
        nnew = N-2
        #m = self.CProfile.n
        #lvertices = np.zeros((nnew,m),float)
        tvertnet = self.make_transverse_net()
        lcurvenet = []
        #N = self.tcurvenet[0].n
        #ptslist = []
        #n=7:
        #for i in range(N):
        # n=11:
        self.save_junk = []
        
        self.interpolation_schemes = {0:'equal spaced ukbar',
                                      1:'centripetal ukbar',
                                      2:'chord length ukbar'}
        
        #go_further = False #we have to go on - to get the right parameterization
        for i in range(1,nnew+1): #[1,2,3,4,5]
            #ptslist.append(tvertnet[:,i])
            """
            TODO: change this back! 
            but modify so that in the event the optimization
            fails, 
            then and only then do we fall back to
            these secondary, inferior methods of linear interpolation
            such that the knot vector is okay
            """
            # centripetal ukbar spacing, average knot spacing:
            try:
                #vanilla_spline = spline.interpolatedBspline(
                #                        tvertnet[:,i], k, nump)
                vanilla_spline = spline.interpolatedBspline(
                        self.solve_linear_list_vertices[:,i], k, nump)
                dm1 = uopt.make_thb(vanilla_spline)
                dm1.ukbar =vanilla_spline.ukbar
                lcurvenet.append(dm1)
                #            try:
                #                # uniform ukbar spacing, uniform knot spacing:
                #                try:
                #                    # uniform ukbar spacing, uniform knot spacing:
                #                    v0 = copy.deepcopy(vanilla_spline)
                #                    v0.GlobalCurveInterp_gen1_THB_compatible_ver0(tvertnet[:,i])
                #                    sm0 = v0.smart_extremes()
                #                except:
                #                    sm0 = [(1.e10,1.e10),(1.e10,1.e10),(1.e10,1.e10)]
                #                
                #                try:
                #                    # centripetal ukbar spacing, uniform knot spacing:
                #                    v1 = copy.deepcopy(vanilla_spline)
                #                    v1.GlobalCurveInterp_gen1_THB_compatible_ver1(tvertnet[:,i])
                #                    sm1 = v0.smart_extremes()
                #                except:
                #                    sm1 = [(1.e10,1.e10),(1.e10,1.e10),(1.e10,1.e10)]
                #                
                #                try:
                #                    # chord length ukbar spacing, uniform knot spacing:
                #                    v2 = copy.deepcopy(vanilla_spline)
                #                    v2.GlobalCurveInterp_gen1_THB_compatible_ver2(tvertnet[:,i])
                #                    sm2 = v0.smart_extremes()
                #                except:
                #                    sm2 = [(1.e10,1.e10),(1.e10,1.e10),(1.e10,1.e10)]
                #                
                #                pick = {0:v0,1:v1,2:v2}
                #                cksum = {
                #                            0:sum(np.linalg.norm(el) for el in sm0),
                #                            1:sum(np.linalg.norm(el) for el in sm1),
                #                            2:sum(np.linalg.norm(el) for el in sm2)
                #                        }
                #                key = min(cksum, key=cksum.get)
                #                #key = min(cksum, key=lambda k: cksum[k])  #equiv code
                #                print 'curve ',i,' interpolation type: ',self.interpolation_schemes[key]
                #                dm1 = uopt.make_thb(pick[key])
                #                dm1.ukbar = pick[key].ukbar
                #                lcurvenet.append(dm1)
            except:
                print '-----------------------------------------------------'
                print 'ERROR: warning'
                print '  Bare Hull Longitudinals'
                print '  Interpolation Failed at curve {}, '.format(i)
                print '  Switching to Bspline starting curve for this longidutinal.'
                print '-----------------------------------------------------'
                #inicurve = spline.Bspline(tvertnet[:,i], k, nump)
                inicurve = spline.Bspline(
                        self.solve_linear_list_vertices[:,i], k, nump)
                
                #nu = len(tvertnet[:,i])
                #nu = len(self.solve_linear_list_vertices[:,i])
                
                #stupid trick to guess some interpolation points
                #a good program would have much less control flow
                #alas, there is no time to add 'another layer'
                # to optimize boat construction hyperperameters
                
                #inicurve.ukbar = np.linspace(0.,1.,nu+2,
                #                             endpoint=True)#[1:-1]  #don't take the inner part here.  
                                                                    #This is done in the function
                                                                    
                ukbar = []
                for i in range(inicurve.n):
                    ukbar.append(inicurve.greville_abscissa(i))
                inicurve.ukbar = ukbar
                    
                lcurvenet.append(inicurve)
                self.save_junk.append(
                        lcurvenet[-1])
                #go_further = True
                
        # store the basic lcurve net created via interpolation 
        # for comparison
        # with the 'optimized' version we are about to create
        self.initial_longitudinals = copy.deepcopy(lcurvenet)
        """Switching to more efficient longitudinal lofting
        via Lspline interpolation of the transverses
        with the interpolated curve as starting point 
        (i.e. we have a good answer)
        point being just to re-parameterize our good answer
        If we make init ADILSpline fast,
        this might beat some other re-parameterization method.
        
        see also geometry_generator.py
        make_longitudinals
        where this issue also arrises.
        #        nnet = []
        #        for curve in lcurvenet:
        #            nnet.append(uopt.match_curve(initial_curve = curve,
        #                                   return_Lspline = False))
        #"""
        print 'correcting longitudinal parameterization'
        nnet = []
        longLsplines = []
        #for i,curve in zip(range(N),lcurvenet):
        #
        # -must change for 11 vs 7 curves:
        #
        if refine:# and go_further: 
            #   it is assumed that the DWL and CPKeel are the 
            #   first and last longitudinal, respectively,
            #   and that they are dyadic-cubic-compatible i.e. 
            #   n = n in {4,5,7,11,19,35...}
            #   Generator of this set (for a cubic curve only):
            #       N(g) = 2*(N(g-1) - N(g-2) ) + N(g-1)
            #   with
            #       N(1) = 4
            #       N(2) = 5
            #
            #  n=11:
            print '-correcting longitudinal parameterization'
            print '-looping over the interior of the lcurvenet'
            if extrema_enforcement is not None:
                print '-enforcing bare hull extremes as inequality constraints'
                for i,curve in zip(range(1,nnew+1),lcurvenet):
                    #pts =  tvertnet[:,i]
                    pts = self.solve_nonlinear_list_vertices[:,i]
                    print '------------------------------------------'
                    Lspline = uopt.interpolate_vertices(
                                initial_vertices = pts,
                                helper_curve     = curve,
                                return_Lspline   = True,
                                make_thb         = True,
                                nose_constraint  = True,
                                Bmax = self.MaxBreadth/2.,
                                hullblockextremes=extrema_enforcement,
                                exact_fit = True,
                                ignore = self.ignore,
                                dropout = self.dropout) 
                    #Lspline.curve = Lspline.curve.dyadic_refinement()
                    if curve.n == 7:
                        new_curve = curve.dyadic_refinement()
                        nnet.append(new_curve)
                    else:
                        assert(curve.n==11),'Oops, wrong number of transverse curves used. n={}'.format(curve.n)
                        nnet.append(Lspline.curve)
                    longLsplines.append(Lspline)
            else:
                print 'ERROR warning,: not enforcing extrema '
                print '             of Longitudinal bare hull form Curves'
                for i,curve in zip(range(1,nnew+1),lcurvenet):
                    #pts =  tvertnet[:,i]
                    pts = self.solve_nonlinear_list_vertices[:,i]
                    print '------------------------------------------'
                    Lspline = uopt.interpolate_vertices(
                                initial_vertices = pts,
                                helper_curve     = curve,
                                return_Lspline   = True,
                                make_thb         = True,
                                nose_constraint  = True,
                                Bmax = self.MaxBreadth/2.,
                                exact_fit = True,
                                ignore = self.ignore,
                                dropout = self.dropout) 
                    if curve.n == 7:
                        new_curve = curve.dyadic_refinement()
                        nnet.append(new_curve)
                    else:
                        assert(curve.n==11),'Oops, wrong number of transverse curves used. n={}'.format(curve.n)
                        nnet.append(Lspline.curve)
                    longLsplines.append(Lspline)
                        
        else:
            #   n=7:
            print 'ERROR: warning, lcurves will not match DWL and CProfile'
            print '-looping over the entirety of the lcurvenet'
            for i,curve in zip(range(N),lcurvenet):
                #pts =  tvertnet[:,i]
                pts = self.solve_nonlinear_list_vertices[:,i]
                print '------------------------------------------'
                curve = uopt.interpolate_vertices(
                            initial_vertices=pts,
                            helper_curve = curve,
                            return_Lspline = False,
                            make_thb=True,
                            nose_constraint = True,
                            Bmax = self.MaxBreadth/2.,
                            exact_fit = True,
                            ignore = self.ignore,
                            dropout = self.dropout) 
                nnet.append(curve)
        
        self.lcurvenet = nnet
        self.Lsplines.lcurves = longLsplines
        self.dropped_transverse = []
        #        self.ignored_transverse = []
        #        self.not_in_longi = []
        #if self.dropout is not None or self.ignore is not None:
        totl = self.dropout + self.ignore 
        dlist = self.use_in_longi_solver
        dlist = [ el for el in dlist if el not in totl]
        d2list = [ self.tcurvenet[i] for i in dlist]
        self.tcurvenet = d2list
        #    for el in totl[::-1]:
        #        if el in self.dropout:
        #            self.dropped_transverse.append(el)
        #            
        #                elif el in self.ignore:
        #                    self.ignored_transverse.append(el)
        #                    self.tcurvenet.remove(el)
        #                elif el not in self.use_in_longi_solver:
        #                    self.not_in_longi.append(el)
        #                    self.tcurvenet.remove(el)
        #self.tcurvenet = self.solve_nonlinear_list_curves
        print '------------------------------------------'
        print 'correction finished'
        print 'Longitudinal Bare Hull Curve Solving is Done'
        print '------------------------------------------'
        return

    def correct_longitudinal_match_hullcurves(self, refine=True):
        """DWL_cdr_b
        CPK_cdr_b
        
        import utility_optimization as uopt
        """
        if self._verbose: print 'hull_from_simple_designspace: correct_longitudinal_match_hullcurves'
        
        
        if self.hullformtype == 'osv':
            dwl = self.DWL_OSV
            clp = self.CPK_cdr_b
        else:
            self.DWL = uopt.make_thb(self.DWL)
            dwl = self.DWL
            clp = self.CPK_cdr_b
            
        #if refine:
            
        assert(self.DWL.n == 11),'DWL is not dyadic, n={}'.format(self.DWL.n)
        assert(self.CProfile.n == 11),'CProfile is not dyadic, n={}'.format(self.CProfile.n)

        self.lcurvenet.insert(0,clp)  #lcurve boundary CProfile gets the split version

        self.lcurvenet.append(dwl)
        return

    def make_longitudinals(self,
                           refine=True,
                           extrema_data=None):
        """interpolate the transverse vertices
        via Piegel and Tiller
        
        Notes:  
        ----------
            
            this function will compute the transverse curve net as well.
        
        
            
            dropout, which curve am I dropping?:
            0: nose curve (do not drop!)
            1: bowfairning
            2: bowtransition
            3: mshp_sections[0] (flat of side)
            4: mshp_sections[0] (flat of side)
            5: mshp_sections[0] (flat of side)
            6: mshp_sections[0] (midship for real)
            7: mshp_sections[0] (midship for real)
            8: sterntransition
            9: sternfairing
            10: aftcurve
        
        """
        #if self._verbose: print 'hull_from_simple_designspace: make_longitudinals'
        print '\n hull_from_simple_designspace: ',' make_longitudinals \n'
        assert(self.sternfairing.n==self.mshp_sections[0].n),"Transvere Curves Do not match!"
        assert(self.sternfairing.n==self.bowfairning.n),"Transvere Curves Do not match!"
        #self.split_hull_curves() #put this at the end of the CPKeel curve generation!
        #assert(self.sternfairing.n==self.CPK_aft.n),"Transvere Curves Do not match!"
        self.correct_hull_curves() #must come after split because fairness curves assume fixed parameter loc!
        self.issue_longitudinals(refine,
                                 extrema_data, #extrema_data sends block hull data to lcurvenet lspline-solver wrapper function
                                 dropout=self.dropout) #dropout to allow free vertices to smooth the longitudinals
        self.correct_longitudinal_match_hullcurves(refine)
        for curve in self.tcurvenet:
            print curve.n
            #assert (curve.n == 7),'transverse curve n = {}'.format(curve.n)
            
        for curve in self.lcurvenet:
            print curve.n
            #assert (curve.n == 11),'longtudinal curve n = {}'.format(curve.n)
            
        self.hullsurf = spline.BsplineSurface(self.tcurvenet,
                                              self.lcurvenet)
        return
    
    def issue_flats_outline(self):
        vertices = []
        for curve in self.tcurvenet[:3]:
            vertices.append(curve.vertices[2])
        for curve in self.tcurvenet[3:6]:
            vertices.append(curve.vertices[3])
        for curve in self.tcurvenet[6:]:
            vertices.append(curve.vertices[2])
        vertices = np.asarray(vertices)
        self.flats_outline = spline.Bspline(vertices, 
                                            k=self.k,
                                            nump=self.nump)
        return
    def issue_detail_curves(self):
        self.issue_flats_outline()
        return
    
    def swap_longitudinal(self, option='direct'):
        """TODO: match initial_longitudinals with 
        curve that have the same Xc and area, etc.....
        BUT the correct knot vector.
        
        How?:- dyadic refinement! 
        
        """
        if option is 'direct':
            self.store_lcurves = copy.deepcopy(self.lcurvenet)
            self.lcurvenet = self.initial_longitudinals
            if (self.initial_longitudinals[0].n == 7):
                for cv in self.lcurvenet:
                    cv.dyadic_refinement()
            else:
                assert(self.initial_longitudinals[0].n == 11),'Oops wrong'
            self.lcurvenet.insert(0,self.CPK_cdr_b)
            self.DWL = uopt.make_thb(self.DWL)
        else:
            print 'option = ',option,' not implemented.'  
            print ' Try option = ','direct'
        return
    
    def make_bulb(self, bbow):
        """
        Deprecated in favor of making this a native
        part of the bbow class itself
        """
        print 'Deprecated bulbous bow curve and surface gen'
#        def make3D(curve,dpt,index):
#            """
#            assumes knot vectors are the 
#            *standard* vector 
#            for a given # of control pts
#            """
#            vv = [0,1,2]
#            vv.pop(index)
#            olv = curve.vertices
#            nlv = np.zeros((curve.n,curve.dim+1),float)
#            nlv[:,vv] = olv[:]
#            nlv[:,index] = dpt
#            curve.vertices = nlv
#            return spline.Bspline(curve.vertices,
#                                  curve.k,
#                                  curve.nump)
        
        def liftvertex(vertex,index):
            vv = [0,1,2]
            vv.pop(index)
            nlv = np.zeros((3),float)
            nlv[vv] = vertex[:]
            nlv[index] = 0.
            return nlv
        
        #bowc = split_nose()
        """
        tcurvenet = [bbow.idealdict['inib_ideal'].Lspline.curve,
                          bbow.idealdict['mb_ideal'].Lspline.curve,
                          bbow.idealdict['fsb_ideal'].Lspline.curve,
                          bowc]
        #"""
        #bvt = DeepVector(liftvertex(bowc.vertices[0],2),
        #                 B=bbow.idealdict['inib_ideal'].base.B)
        
        idnet = [(bbow.idealdict['inib_ideal'].Lspline.curve,
                  bbow.idealdict['inib_ideal'].base.v[0],
                  2),
                 (bbow.idealdict['mb_ideal'].Lspline.curve,
                  bbow.idealdict['mb_ideal'].base.v[0],
                  2),
                 (bbow.idealdict['fsb_ideal'].Lspline.curve,
                  bbow.idealdict['fsb_ideal'].base.v[0],
                  2),
                 (bbow.idealdict['cp_nose'].Lspline.curve,
                  0.,
                  1)]
                 
        """
        Done, Nov 3 2017:
            make the bow curves 3D
            make the bow surface
        TODO, Nov 3 2017:
            THB...
        #"""
        d3net = []
        for tpl in idnet:
            curve   = tpl[0]
            dv      = tpl[1]
            index   = tpl[2]
            #d3net.append(make3D(curve,dv,index))
            d3net.append(uopt.lift_curve(curve,
                                         elevation=dv,
                                         index=index))
            
#        tvertnet = np.asarray([d3net[0].vertices,
#                               d3net[1].vertices,
#                               d3net[2].vertices,
#                               d3net[3].vertices])
        
        bbow.tcurvenet = d3net
        
        
        bbow.lcurvenet = self.make_other_longitudinals(bbow.tcurvenet)
        
        bbow.hullsurf = spline.BsplineSurface(bbow.tcurvenet,
                                              bbow.lcurvenet)
        return  #bbow
    
    

    def make_other_longitudinals(self, tcurvenet):
        """tcurvenet = bbow.tcurvenet
        
        import curve as spline
        """
        k=tcurvenet[0].k
        nump = tcurvenet[0].nump
        N = tcurvenet[0].n
        nnew = N
        #m = self.CProfile.n
        #lvertices = np.zeros((nnew,m),float)
        tvertnet = np.asarray([el.vertices for el in tcurvenet])
        
        
        lcurvenet = []
        for i in range(N):
            lcurvenet.append(
                spline.interpolatedBspline(
                tvertnet[:,i], k, nump))
            
        """TODO, TLM OCT 30 2017:
            make this use 3D optimization
            to optimize the longitudinals
            and hit the vertex targets!
            (as point interpolations of course)
        """
        
        
        return lcurvenet
    
    def export_curves(self):
        """Note that this is a standard B-spline hull.
        Here the tcurvenet is the same thing
        as the THB hullcurve net
        """
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
        return
    
    def export_rhino(self):
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

#class controls(self):
#    def __init__(self, hull):
#        self.hull = hull

def compute_hull(hulldesign, flat_=True):
    """TODO:
        all hull curves can be made in PARALLEL!
        all section curves can be made in parallel.
    """
    print '------------------------------------------'
    #hulldesign.compute_SAC_new(flat=True)
    hulldesign.compute_SAC(flat=flat_) #single piece
    print '------------------------------------------'

    hulldesign.compute_cLProfile(flat=flat_)
    print '------------------------------------------'
    np.rad2deg(hulldesign.CProfile2D.compute_tangent(0.))
    np.rad2deg(hulldesign.CProfile2D.compute_tangent(1.))

    print '------------------------------------------'
    hulldesign.compute_DWL(flat=True, alphab=30.,alphae=-30.)#35.,35. works #25.,25. works#15., -15. works# 0.,0. works#

    print '------------------------------------------'
    hulldesign.compute_midship_section()

    print '------------------------------------------'
    hulldesign.compute_bow_fairing_curve()
    print '------------------------------------------'
    hulldesign.compute_bow_transition_curve()
    print '------------------------------------------'
    hulldesign.compute_stern_fairing_curve()
    print '------------------------------------------'
    hulldesign.compute_stern_transition_curve()
    print '------------------------------------------'

    
    
    #    the_filename = 'transverse_knots.txt'
    #    with open(the_filename, 'r') as f:
    #        tknots = []
    #        ci = []
    #        for line in f:
    #            if not 'new' in line:
    #                ci.append([float(el) for el in line.split()])
    #            else:
    #                tknots.append(ci[0])
    #                ci = []
    return

def compute_hull_and_transverse_curves(hulldesign, 
                                       sacflat_=True):
    """
    
    Efficiency idea:
    ----------------
        all hull curves can be made in parallel
        all section curves can be made in parallel.
        
    this function 
    gets the basic hull together up to but not
    including the longitudinal net
    nor compiling the transverse net
    """
        
    ##
    ##******************************************************************
    ## Handle aft SAC edge this way:
    DROP        = hulldesign.MaxDraft*hulldesign.aft_drop
    # aft section area must be 
    #               -greater than a triangle
    #               -less than a box.
    Abox = DROP*hulldesign.MaxBreadth/2.
    Atri = .5*Abox
    Aaft = Atri+.75*(Abox-Atri)
    ##
    ##******************************************************************
    ## SAC entrance...
    Atube = hulldesign.A_mid #not A_lateral!  This is the Max SAC of the bulb
    ##
    ##******************************************************************
    ## 
    #hulldesign.compute_SAC_new(flat=True)
    print '------------------------------------------'
    """
    TDB:  issue - SAC really does need to gracefully go to zero 
    at the fore of the entrance.
    
    -idea for SAC specs:
        yb = 0.
        
        FPD.add_yPointConstraint(kind = 'equality',
                                 location= 0.1,
                                 value = Atube,
                                 )
        
    Q: But what location for a local spot to incorporate the bow
        i.e. - how to make sure there
        is enough area for the bbow + deadrise to deck?
    A: No.  Do not go that route.
    Design SAC without regard to it.
        add something to tradeoff sectional area
        between deadrise and deck
        and design both deck curve and bbow to meet this trade
        
    """
    hulldesign.compute_SAC(flat=sacflat_,
                           Ab = Atube, 
                           Ae = Aaft) #single piece
                           #Ab = 0., 
                           #Ae = Aaft) #single piece
    print '------------------------------------------'
    hulldesign.compute_cLProfile(flat=True, 
                                 alphab=80.,
                                 alphae=-80.,
                                 drop_aft = hulldesign.aft_drop)
    np.rad2deg(hulldesign.CProfile2D.compute_tangent(0.))
    np.rad2deg(hulldesign.CProfile2D.compute_tangent(1.))
    print '------------------------------------------'
    hulldesign.compute_DWL(flat=True, 
                           alphab=80.,alphae=-80.)#35.,35. works #25.,25. works#15., -15. works# 0.,0. works#
    print '------------------------------------------'
    hulldesign.compute_midship_section()
    print '------------------------------------------'
    hulldesign.compute_bow_fairing_curve()
    print '------------------------------------------'
    hulldesign.compute_bow_transition_curve()
    print '------------------------------------------'
    hulldesign.compute_stern_fairing_curve()
    print '------------------------------------------'
    hulldesign.compute_stern_transition_curve()
    print '------------------------------------------'
    return

def loft_hull(hulldesign, 
              refine=True,
              extrema=None):
    hulldesign.make_longitudinals(refine,
                                  extrema)
    print '------------------------------------------'
    return 



def hull_gui(hulldesign):
    class wrapper(object):
        def __init__(self, curve):
            self.curve = curve

    HF = hull_frame()

    HF.add_section(hulldesign.Lsplines.bowfairning,
                   hulldesign.bowfairning.vertices[0,2])


    for el, key in zip(hulldesign.Lsplines.midSS,
                       hulldesign.mshp_sections):
        HF.add_section(el, key.vertices[0,2])



    HF.add_section(hulldesign.Lsplines.sternfairing,
                   hulldesign.sternfairing.vertices[0,2])

    #HF.add_bow(hulldesign.Lsplines.DWL,0.)
    #HF.add_bow(hulldesign.Lsplines.CProfile,0.)

    HF.add_bow(wrapper(hulldesign.DWL),0.)
    HF.add_bow(wrapper(hulldesign.CProfile),0.)

    #for lcurve in hulldesign.lcurvenet:
    #    HF.add_bow(wrapper(lcurve),0.)

    fig = plt.figure(figsize=(4, 4))
    ax = DrawCurveInteractive(fig, sections=HF, xkcd=False)
    fig.add_axes(ax)
    ax.draw_curves()
    #hulldesign = ax
    
    hull_gui_no_Lsplines(hulldesign)
    curves_gui_(hulldesign)
    return ax



def hull_gui_no_Lsplines(hulldesign):
    """Attempt to provide full on 3D curves in the 'real GUI'
    """
    class wrapper(object):
        def __init__(self, curve):
            self.curve = curve

    HF = hull_frame()
    #
    #**********************************************************************
    # Fairness
    HF.add_section(hulldesign.Lsplines.bowfairning,
                   hulldesign.bowfairning.vertices[0,2])



    HF.add_section(hulldesign.Lsplines.sternfairing,
                   hulldesign.sternfairing.vertices[0,2])

    #HF.add_bow(hulldesign.Lsplines.DWL,0.)
    #HF.add_bow(hulldesign.Lsplines.CProfile,0.)

    HF.add_bow(wrapper(hulldesign.DWL),0.)
    HF.add_bow(wrapper(hulldesign.CProfile),0.)

    #for lcurve in hulldesign.lcurvenet:
    #    HF.add_bow(wrapper(lcurve),0.)
    
    

    for curve in hulldesign.mshp_sections:
        HF.add_3Dsection(wrapper(curve))
        
        
    for curve in hulldesign.tcurvenet:
        HF.add_3Dsection(wrapper(curve))

    fig = plt.figure(figsize=(4, 4))
    ax = DrawCurveInteractive(fig, 
                              sections=HF, 
                              xkcd=False)
    fig.add_axes(ax)
    ax.draw_curves()
    #hulldesign = ax
    return ax

def curves_gui_(hulldesign):
    """Use full on 3D curve nets
    
    
    
    dev
    ----------
    
    from hull_from_designspace import hull_frame
    
    hulldesign=SD.hull
    
    import matplotlib.pyplot as plt
    
    from adials_gui3d import DrawCurveInteractive
    
    """
    class wrapper(object):
        def __init__(self, curve):
            self.curve = curve

    #
    #**********************************************************************
    #  Init markers
    HF = hull_frame()
    #
    #**********************************************************************
    # Fairness
    HF.add_section(hulldesign.Lsplines.bowfairning,
                   hulldesign.bowfairning.vertices[0,2])



    HF.add_section(hulldesign.Lsplines.sternfairing,
                   hulldesign.sternfairing.vertices[0,2])
    
    #
    #**********************************************************************
    # mid ship curves
    
    for curve in hulldesign.mshp_sections:
        HF.add_3Dsection(wrapper(curve))
        
        
    #
    #**********************************************************************
    # hull curves!
    
    for curve in hulldesign.lcurvenet:
        HF.add_3Dsection(wrapper(curve))
    for curve in hulldesign.tcurvenet:
        HF.add_3Dsection(wrapper(curve))
        
    
    #
    #**********************************************************************
    # make GUI

    fig = plt.figure(figsize=(4, 4))
    ax = DrawCurveInteractive(fig, 
                              sections=HF, 
                              xkcd=False)
    fig.add_axes(ax)
    ax.draw_curves()
    #
    #**********************************************************************
    # attach gui
    hulldesign.ax = ax
    return ax


def quicktest():
    k       = self.k
    nump    = self.nump
    #
    #**********************************************************************
    # 
    cpk = self.CProfile2D  # more efficient to use 2D 
                           # when reparameterizing w/ knot removal
    #
    #**********************************************************************
    # points we absolutely must hit
    #pwamh = [] #doing knot removal... things seem ... improved... xluckx
    #
    #**********************************************************************
    # 
    fp_index = False #no regard for types, do I seem to have.
    # in my ignorance, I apologize, Vladimir Voevodsky.
    # if this doesn't work, I promise to solve a knot vector
    # reparameterization problem via homotopic deformation of the knots.
    #..... which did not work either... win some lose some.
    #
    if cpk.dim ==2:
        fp_index = 0
    elif cpk.dim ==3:
        fp_index = 0
    else:
        mssg = '\n ERROR: Fatal; Could not find dimensions of CPKeel curve to be split \n'
        assert(fp_index),mssg
    #
    #**********************************************************************
    # Which is outside the other: flat of Cprofile?, 
    #    or intersections with fairness curves?
    #longi_ini,longi_fini = self.stations.FOCP[:2] #these are vertex constraints
    # not curve point constarints!
    #
    bstart,bend = self.stations.Bow_Fairing[:2]
    sstart,send = self.stations.Stern_Fairing[:2]
    #
    #bow_split_point = min(bstart,sstart)
    #stern_split_point = max(bend,send)
    bow_split_point = min(bstart,bend)
    stern_split_point = max(sstart,send)
    #
    #**********************************************************************
    #
    #
    #**********************************************************************
    # find the location at which to cut the CProfile curve for
    #   the fwd portions of the longitudinal keel curve (and fore and aft transverses)
    #
    #start,end,ao,bo = self.stations.Bow_Fairing
    b1 = cpk.FindPoint(bow_split_point,
                       fp_index) 
    mssg = '\n ERROR: Fatal; Could not find bow fairing intercept with CProfile curve\n'
    assert(b1[2]),mssg
    #
    #**********************************************************************
    # cut the bow CProfile just fwd of the bow fairness curve
    #
    bb = b1[0] - 0.01*b1[0]
    self.CPK_fwd, self.CPK_cdr_a = cpk.CurveSplit(bb) #split fwd end
    #
    #**********************************************************************
    # Bow split done.  Now Aft split
    #**********************************************************************
    # find the location at which to cut the CProfile curve for
    #   the aft parts of the longitudinal keel curve (and fore and aft transverses)
    #
    #start,end,ao,bo = self.stations.Stern_Fairing
    b2 = self.CPK_cdr_a.FindPoint(stern_split_point,
                                  fp_index) 
    mssg = '\n ERROR: Fatal; Could not find bow fairing intercept with CProfile curve\n'
    assert(b2[2]),mssg
    #
    #**********************************************************************
    # CProfile aft split
    # cut the stern CProfile just aft of the stern fairness curve
    sb = b2[0] + 0.01*(1.-b2[0])
    self.CPK_cdr_b, self.CPK_aft = self.CPK_cdr_a.CurveSplit(sb)
    #
    #**********************************************************************
    # start the constraints for this curve, 
    # but delay computation a bit ;) 
    # (because we don't have a curve with the correct knot vector 
    # with which to instantiate the FPD class)
    CI = LazyFPDConstraintInterface() #lazy FPD factory simple magic
    #
    """Note that at present
    we are __NOT__ solving an optimziation problem
    for the good parameterization.
    Instead, we are doing knot insertion and removal to 
    fix up the curve whilst changing it as little as possible
    'CI' is here though, as example of 
    'where I think it fruitful to go'
    with abstracting out FPD setup
    """
    #
    #**********************************************************************
    # checks in this block are
    #      no longer relevant (and causing a great deal of headache):
    #
    #      but there is a lesson here:  split should be both
    #          -'outside' of fairness curves 
    #          (so that they meet the longitudinal keel, obviously)
    #          - and 'outside' of the flat of CProfile 
    #          (to arrive at a 'best' set of longitudinals on the interior)
    #
    #longi_ini,longi_fini = self.stations.FOCP[:2]
    #param_a,point_a,sucess_a=self.CPK_cdr_b.FindPoint(longi_ini,0) #not there
    #param_b,point_b,sucess_b=self.CPK_cdr_b.FindPoint(longi_fini,0) #this one is there
    #mssg = 'could not find aft point to match on the CPKeel after split'
    #assert(sucess_b),mssg
    #
    #**********************************************************************
    # area from the curve to the x-axis (tinkering with translation to axis)
    #        minx,maxx = self.CPK_cdr_b.extremes1d(0)
    #        miny,maxy = self.CPK_cdr_b.extremes1d(1)
    #        minz,maxz = self.CPK_cdr_b.extremes1d(2)
    #        normalized_vertices = spline.move_copy_vertices([-minx,-miny],
    #                                                        self.CPK_cdr_b.vertices)
    #        test_curve = copy.copy(self.CPK_cdr_b)
    #        test_curve.vertices = normalized_vertices
    #
    #**********************************************************************
    # AREA from the curve to the x-axis
    self.CPK_cdr_b.compute_area()
    area = self.CPK_cdr_b.area
    CI.Constraint_Area(value=area,
                       kind='equality')
    #
    #**********************************************************************
    # X centroid
    self.CPK_cdr_b.compute_moments()
    self.CPK_cdr_b.computeXc()
    Xc = self.CPK_cdr_b.Xc
    CI.Constraint_Xc(value = Xc,
                     kind = 'equality')
    return

def test():
    self.define_hull()
    #
    self.compute_SAC_new()
    print '------------------------------------------'
    self.compute_cLProfile()
    print '------------------------------------------'
    self.compute_DWL()
    print '------------------------------------------'
    self.compute_midship_section()
    print '------------------------------------------'
    self.compute_bow_fairing_curve()
    print '------------------------------------------'
    self.compute_bow_transition_curve()
    print '------------------------------------------'
    self.compute_stern_fairing_curve()
    print '------------------------------------------'
    self.compute_stern_transition_curve()
    print '------------------------------------------'
    self.make_longitudinals()
    print '------------------------------------------'
    #self.split_hull_curves() #moved within CLProfile Keel
    #quickdirty()
    return

if __name__=='__main__':
    frac=.5
    h1 = hullclp()

    Vsmall = 2000.
    Vbig = 300000.
    Len = 120. #120.
    BLEN = 20. #10.
    DLEN = 20. #20.

    h1.vol =ia(Vsmall,Vbig)
    h1.lwl = ia(120.,Len)
    h1.lfwl = ia(0.,Len)
    h1.lfsac = ia(0.,Len)
    h1.lfos = ia(0.,Len)
    h1.bwl = ia(15.,BLEN) #20.#
    
    h1.draft = ia(DLEN,DLEN)
    #h1.draft = ia(DLEN-h1.SMALL,DLEN)

    #    h1.vol =ia(2000.,300000.)
    #    h1.lwl = ia(120.,120.)
    #    h1.lfwl = ia(0.,120.)
    #    h1.lfsac = ia(0.,120.)
    #    h1.lfos = ia(0.,120.)
    #    h1.bwl = ia(20.,20.) #20.#

    #self.Cb = ia(.5,.7)
    
    h1.Cwp = ia(.2,.8)  #thin boat cannot have much FOS!!  Constraints need to reflect this
    h1.Cmidshp =  ia(.7,.9999)
    h1.dsmax =  ia(DLEN,DLEN)
    h1.compute_prismatic_coefficient()
    h1.compute_Cb()
    h1.compute_midship_coefficient()
    h1.compute_prismatic_coefficient()
    h1.print_state()
    h1.Cp  = ia(0.,1.)
    h1.Cb  = ia(0.,1.)
    h1.Ccp = ia(.5,.8)
    h1.lfcp = ia(0.,100.)

    h1.AC_revise(print_=False)
    h1.print_state()

    #self = Hull(h1) #crucial for pointers
    self = Hull(hullconstraints = h1)
    self.get_flat_of_center_plane(frac=.5)#.01 is good#
    #self.get_LWL(frac=.5)

    #self.get_bwl(frac=.5)
    #self.get_draft(frac=.5)

    #self.get_Vol(frac=.5)

    #self.hullconstraints.AC_revise(print_=False)

    #self.hullconstraints.Cb.arcs[0]()
    #self.hullconstraints.Cb.arcs[1]()
    #self.get_wlArea(frac=.5)
    #self.hullconstraints.AC_revise(print_=True)
    #self.arc_set = set([])

    #self.get_maxSecArea(frac=.999)
    #self.hullconstraints.AC_revise(print_=True)

    #self.get_Cmidship(frac=.5)
    #self.hullconstraints.AC_revise(print_=True)

    self.hullconstraints.print_state()

    self.get_Vol(frac=.7)#.7#.8#.6

    self.get_LWL(frac=.5)
    self.get_bwl(frac=.5)
    self.get_draft(frac=.5)
    self.hullconstraints.AC_revise()
    #self.hullconstraints.print_state()

    self.get_Cb(frac=.7)#.75#.85#.6
    self.hullconstraints.print_state()

    self.get_Cmidship(frac=1.)#frac=.5 is good as well#
    self.hullconstraints.print_state()

    self.get_maxSecArea(frac=.5)
    self.hullconstraints.print_state()

    self.hullconstraints.AC_revise()
    self.get_wlArea(frac=.5)#.75 is good#.9 is good#.5 from .75
    self.get_clArea(frac=.5)#.42 is good# .75 is good#.9 is good#.5 from .75
    #self.hullconstraints.print_state()

    self.hullconstraints.AC_revise()
    self.hullconstraints.print_state()
    #    self.get_LWL(frac=.5)
    #
    #    self.get_bwl(frac=.5)
    #    self.get_draft(frac=.5)

    #self.get_len_flat_wL(frac=.5)
    #self.get_len_flat_SAC(frac=.5)
    #self.get_len_FOS(.5)

    self.get_len_flat_SAC(frac=.1)
    self.get_len_flat_wL(frac=.1)
    self.get_len_FOS(.3)#.2 #doesn't do anything yet
    self.hullconstraints.AC_revise()
    #self.hullconstraints.print_state()


    ##SAC sections:
    x=self.hullconstraints.state(self.hullconstraints.lfsac).getpoint(.9)
    re = (Len - x)/2.
    self.hullconstraints.SAC_entrance_len = ia(0.,re)

    v = self.hullconstraints.state(self.BareHullVolume)
    sma = self.hullconstraints.state(self.hullconstraints.SAC_mid_area) #SAC area == Volume
    re = (v-sma)/2.
    self.hullconstraints.SAC_entrance_area = re
    #self.get_len_FOS(.9)
    self.hullconstraints.AC_revise()

    #self.get_len_FOS(.3)
    #self.get_len_flat_wL(frac=.5)
    #self.get_len_flat_SAC(frac=.49)
    """
    TODO:  think about the following:
    # For the bow fairness curve sectional area info:
    #   the GA should choose a value here,
    # just like the rest of the self.hullconstraints:
    #"""
    self.hullconstraints.Cbfc = ia(.1,1.)

    #self.hullconstraints.Abfc = self.hullconstraints.state.value_of(
    #    self.hullconstraints.Amsh)/(7.)

    #    self.hullconstraints.bbfc = self.hullconstraints.state.value_of(
    #        self.hullconstraints.bwl)/(5.)
    #    self.hullconstraints.print_state()
    #

    self.hullconstraints.Csfc = ia(.1,1.) #Coeff stern fairness curve
    #self.hullconstraints.Asfc = self.hullconstraints.state.value_of(
    #        self.hullconstraints.Amsh)/(7.)

    """THERE is an unspecified relation between Asfc Abfc and
    A cL Keel Profile and A DWL!!!*******

    These hullconstraints.Asfc sets were too restrictive.
    """
    self.hullconstraints.AC_revise()
    #self.hullconstraints.print_state()

    self.set_LCG()#cb!
    self.set_SAC_fwd_Xc()
    self.set_SAC_mid_Xc()
    self.set_SAC_run_Xc()
    self.hullconstraints.AC_revise()
    #self.hullconstraints.print_state()
    #self.hullconstraints.bsfc = self.state.value_of(self.bwl)/(10.)
    """#  With the above though, the issue is
    # self.hullconstraints are 'one level more abstract'
    # than the actuall hull parmeters
    #"""

    #move up...?
    self.get_len_entrance_SAC(frac=.5)
    self.get_len_run_SAC(frac=.5)
    self.hullconstraints.AC_revise()
    #self.hullconstraints.print_state()

    #if False:
    self.bow_fairness_section_Coefficient(frac=.7)
    self.bow_fairness_section_beam(frac=.4)
    self.bow_fairness_section_draft(frac=.7)
    self.bow_fairness_section_area(frac=.5)
    self.hullconstraints.AC_revise()
    #self.hullconstraints.print_state()

    self.stern_fairness_section_Coefficient(frac=.7)
    self.stern_fairness_section_beam(frac=.7)
    self.stern_fairness_section_draft(frac=.4)
    self.stern_fairness_section_area(frac=.5)
    self.hullconstraints.AC_revise()
    #self.hullconstraints.print_state()



    self.hullconstraints.print_state()
    self.definestations()

    #"""
    #    s=Stations(lwl=self.LengthDWL,cp=.5)
    #    s.FOSAC = [self.flat_SAC, .5]
    #    s.FOWL = [self.flat_wL, .5]
    #    s.FOS = [self.FOS, .5]

    #self.compute_FOSAC(frac=.5) #same frac for both at the moment




    """
    compute_hull(self)
    hull_gui(self)
    #"""
    """
    self.Lsplines.sternfairing.curve.compute_area_to_y()
    self.Lsplines.sternfairing.curve.area_to_y
    self.Lsplines.bowfairning.curve.compute_area_to_y()
    self.Lsplines.bowfairning.curve.area_to_y

    self.sternfairing.vertices
    start,end,ao,bo = self.stations.Stern_Fairing
    self.SAC.CurvePoint(ao)/2.
   # self.SAC.


    #self.plot_primary_curves()

    #self.SAC2D.plotcurve_detailed()
    #self.DWL2D.plotcurve_detailed()
    #"""
    #"""
    #"""