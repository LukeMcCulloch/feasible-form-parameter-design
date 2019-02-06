# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 18:55:13 2016

@author: lukemcculloch

probably best called:
    Constraint Interval Logic Programming on Objects for Form parameter Design

    using Hull Inference Constraint Logic Programming

    correct expression for the correctly masked SAC gradient:
        sac = self.Lsplines.SAC
        mask = sac.mask
        sac.f.grad[:].T[mask]
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
#from interval_arithmetic import ia
#import uKanren as lp
#from hull_inference_ob_graph import Hull as hullclp
from extended_interval_arithmetic import ia
import sqKanren as lp
from design_space import HullSpace as hullclp
#
#from OldHullDesigner import Hullcurve#, linear_vertices, transpose_vertices
import curve             as spline
from   ADILS             import IntervalLagrangeSpline, Lagrangian
from   adials_gui3d      import frame, DrawCurveInteractive, hull_frame
from   FormParameter     import FormParameterDict
from   initialValues     import InitializeControlPoints, InitializeControlVertices
from   initialValues     import interval_bounds, lagrangian_bounds
from   plots             import Plotter
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

class Stations(object):
    """Stations are defined with respect to Hull.LWL"""
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




class Hullcurve(object):

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
    def __init__(self, Amax, #v1,v2,v3,
                 l1,l2,l3,
                 Xe, Xm, Xr, Xlcb=None,
                 Abfc=None, Abfc_loc=None, Asfc=None, Asfc_loc=None,
                 Afp=None, Aap=None, Cpe=None, Cpr=None,
                 stations=None):
        #self.V         = vol
        self.Amax       = Amax
        #self.Ve        = v1
        #self.Vm        = v2
        #self.Vr        = v3
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
        self.stations   = stations
        self.go()
        self.aggregate_SAC_curves_Small()
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
    #def compute_run_SAC(self):
        k       = SAC.k
        nump    = SAC.nump
        nCV     = 7
        xb      = SAC.xb #0.
        xe      = self.Le #56
        yb      = SAC.yb
        ye      = self.Amax
        #yfp     = self.Afp
        Xc      = self.Xe
        ab = 0.
        ae = 0.
        ini_v = InitializeControlVertices(xb,yb,xe,ye,
                                          alphae=ae,
                                          alphab=ab,
                                          Cab_given=0.,
                                          Cae_given=0.,
                                          nCV = 7,
                                          slope = 'up')
        curve = spline.Bspline(ini_v.vertices, k, nump)

        interval_data, small = interval_bounds(curve)

        FPD = FormParameterDict(curve)
        #FPD.add_AreaConstraint(kind='equality', value = curve_area)
        FPD.add_AngleConstraint(kind = 'LS', value = 0., location = 0.)
        FPD.add_AngleConstraint(kind = 'equality', value = 0., location = 1.)
        #FPD.add_CurvatureConstraint(kind = 'equality', value = 0., location = 0.)
        FPD.add_CurvatureConstraint(kind = 'equality', value = 0., location = 1.)
        #FPD.add_XcConstraint(kind = 'equality', value = Xc)
        #
        #----------------------- SAC value at bow fairness curve location:
        #
        start,end,ao,bo = self.stations.Bow_Fairing
        AreaBowFairingCurve = self.Abfc
        where = ao
        FPD.add_yPointConstraint(kind='equality',
                                 location = where,
                                 value = AreaBowFairingCurve)
        FPD.add_xPointConstraint(kind='equality',
                                 location = where,
                                 value = start)
        #
        #-----------------------
        #
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        #FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        #        Lspline.display(mytitle = 'SAC_run_ini', x_axis_name = 'x',y_axis_name = 'z')
        Lspline.optimize(stop = 50)
        #        Lspline.display(mytitle = 'SAC_run',
        #                        x_axis_name = 'x',
        #                        y_axis_name = 'z')
        return Lspline


    def compute_midbody_SAC(self):
        k       = SAC.k
        nump    = SAC.nump
        n       = 7
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
        #FPD.add_XcConstraint(kind = 'equality', value = 6.)
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
    #def compute_entrance_SAC(self):
        """
        """
        k       = SAC.k
        nump    = SAC.nump
        nCV     = 7
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
        #FPD.add_AreaConstraint(kind='equality', value = curve_area)
        FPD.add_AngleConstraint(kind = 'equality', value = 0., location = 0.)
        FPD.add_AngleConstraint(kind = 'LS', value = 0., location = 1.)
        FPD.add_CurvatureConstraint(kind = 'equality',
                                    value = 0.,
                                    location = 0.)
        #FPD.add_CurvatureConstraint(kind = 'equality',
        #                            value = 0.,
        #                            location = 1.)
        #
        #----------------------- SAC value at stern fairness curve:
        #
        start,end,ao,bo = self.stations.Stern_Fairing
        AreaSternFairingCurve = self.Asfc#/2.
        where = ao
        FPD.add_yPointConstraint(kind='equality',
                                 location = where,
                                 value = AreaSternFairingCurve)
        FPD.add_xPointConstraint(kind='equality',
                                 location = where,
                                 value = start)
        #
        #-----------------------
        #
        #FPD.add_XcConstraint(kind = 'equality', value = Xc)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        #FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)

        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        #Lspline.compute_lagrangian()
        #        Lspline.display(mytitle = 'SAC_fwd_ini', x_axis_name = 'x',y_axis_name = 'z')
        Lspline.optimize(stop = 50)
        #        Lspline.display(mytitle = 'SAC_fwd',
        #                        x_axis_name = 'x',
        #                        y_axis_name = 'z')
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
    def __init__(self, hullconstraints):
        """
            TODO: refactor the gets to @property
        """
        self.k                  = 4
        self.n                  = 7
        self.nump               = 30
        self.small              = 1.e-3
        self.tol                = 1.e-6
        self.hullconstraints    = hullconstraints
        self.stations           = {}
        self.Lsplines = Lsplines()
        return

    def __str__(self):
        """print representation
        """
        keys = self.__dict__.keys()
        for key in keys:
            print key, self.__dict__[key]
        return ''

    def get_transverse_start_end_param_loc(self, loc_start, loc_end, which=2):
        """
        all locations are longitudinal:
        i.e. get the longitudinal location some transverse's start and end
            loc_start   : real longi location start
            loc_end     : real lingi location end
        """
        assert(isinstance(self.SAC, spline.Bspline))
        if self.SAC.dim ==2:
            sacwhich =0
        else:
            sacwhich = which
        assert(isinstance(self.DWL, spline.Bspline))
        assert(isinstance(self.CProfile, spline.Bspline))
        a = self.SAC.FindPoint(loc_start,sacwhich)
        b = self.CProfile.FindPoint(loc_start,which)
        c = self.DWL.FindPoint(loc_end,which)
        return a,b,c

    def get_transverse_start_end_points(self, pa,pb):
        assert(0.0<=pa<=1.0),"pa !=[0,1] => {}".format(pa)
        pt1 = self.CProfile.CurvePoint(pa)  #start!
        pt2 = self.DWL.CurvePoint(pb)
        return pt1, pt2

    def get_LWL(self, frac=.5):
        """Unify and propogate constraints
        max waterline length
        """
        self.LengthDWL = self.hullconstraints.state(
                                self.hullconstraints.lwl).getpoint(frac)
        self.hullconstraints.lwl = ia(self.LengthDWL-self.small,
                                      self.LengthDWL+self.small)#unify and propogate
        #self.definestations() #no -> do this after defining all basics!
        self.hullconstraints.AC_revise(print_=False)
        return

    def get_bwl(self, frac=.5):
        """Unify and propogate constraints
        max waterline breadth
        """
        bwl = self.hullconstraints.state(
                                self.hullconstraints.bwl)
        if isinstance(bwl, ia):
            self.MaxBreadth = bwl.getpoint(frac)
        elif isinstance(bwl, float):
            self.MaxBreadth = bwl

        self.hullconstraints.bwl = bwl#ia(self.MaxBreadth-self.small,
                                      #self.MaxBreadth+self.small)
        self.hullconstraints.AC_revise(print_=False)
        return

    def get_draft(self, frac=.5):
        """Unify and propogate constraints
        design waterline : draft
        """
        self.MaxDraft = self.hullconstraints.state(
                                self.hullconstraints.draft).getpoint(frac)
        self.hullconstraints.draft = ia(self.MaxDraft-self.small,
                                        self.MaxDraft+self.small)
        self.hullconstraints.AC_revise(print_=False)
        return

    def get_Vol(self, frac=.5):
        """Unify and propogate constraints
        design waterline : volume
        """
        self.BareHullVolume = self.hullconstraints.state(
                                self.hullconstraints.vol).getpoint(frac)
        self.hullconstraints.vol = ia(self.BareHullVolume-self.small ,
                                      self.BareHullVolume+self.small )
        self.hullconstraints.AC_revise(print_=False)
        return

    def get_wlArea(self, frac=.5):
        """Unify and propogate constraints
        design waterline : area
        """
        self.WLarea = self.hullconstraints.state(
                                self.hullconstraints.Awp).getpoint(frac)
        self.hullconstraints.Awp = ia(self.WLarea-self.small,
                                      self.WLarea+self.small)
        self.hullconstraints.AC_revise(print_=False)
        return

    def get_clArea(self, frac=.5):
        """Unify and propogate constraints
        design waterline : area
        """
        self.CLarea = self.hullconstraints.state(
                                self.hullconstraints.Acp).getpoint(frac)
        #        self.hullconstraints.Acp = ia(self.CLarea-self.small,
        #                                      self.CLarea+self.small)
        self.hullconstraints.Acp = ia(0.,
                                      self.CLarea+self.small)
        self.hullconstraints.AC_revise(print_=False)
        return

    def get_maxSecArea(self, frac=.5):
        """Unify and propogate constraints
        design waterline : max section area

        Note: Amsh and maxSecArea are double sided!
        """
        maxSecArea = self.hullconstraints.state(
                                self.hullconstraints.Amsh).getpoint(frac)
        self.hullconstraints.Amsh = ia(maxSecArea-self.small,
                                          maxSecArea+self.small)
#        Cmidshp = self.hullconstraints.state(
#                                self.hullconstraints.Cmidshp).getpoint(frac)
#        self.hullconstraints.Cmidshp = ia(Cmidshp-self.small,
#                                          Cmidshp+self.small)
        self.maxSecArea = maxSecArea #* Cmidshp
        self.hullconstraints.AC_revise(print_=False)
        return

    def get_Cmidship(self, frac=.5):
        """Unify and propogate constraints
        design waterline : max section area
        """
        self.Cmidship = self.hullconstraints.state(
                                self.hullconstraints.Cmidshp).getpoint(frac)
        self.hullconstraints.Cmidshp = ia(self.Cmidship-self.small,
                                          self.Cmidship+self.small)
        self.hullconstraints.AC_revise(print_=False)
        return


    def get_Cb(self, frac=.5):
        """Unify and propogate constraints
        design waterline : max section area
        """
        self.Cb = self.hullconstraints.state(
                                self.hullconstraints.Cb).getpoint(frac)
        self.hullconstraints.Cb = ia(self.Cb-self.small,
                                     self.Cb+self.small)
        self.hullconstraints.AC_revise(print_=False)
        return


    def get_len_flat_wL(self, frac=.5):
        """Unify and propogate constraints
        flat of wL area curve

        !! issue with inequality:
            flat of wL does not bound wL_area/bwl  !!
        """
        self.flat_wL = self.hullconstraints.state(
                                self.hullconstraints.lfwl).getpoint(frac)
        #self.hullconstraints.lfwl = ia(0.,self.flat_wL+self.small)
        self.hullconstraints.lfwl = ia(self.flat_wL-self.small,
                                       self.flat_wL+self.small)
        self.hullconstraints.AC_revise(print_=False)
        return


    def get_len_FOS(self, frac=.5):
        """Unify and propogate constraints
        flat of sac

        !! issue with inequality:
            flat of sac does not bound vol/(bwl*draft)  !!
        """
        self.FOS = self.hullconstraints.state(
                                self.hullconstraints.lfos).getpoint(frac)
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
        self.LCG = self.hullconstraints.get_val(self.hullconstraints.LCG).getpoint(frac)
        self.hullconstraints.LCG = ia(self.LCG-self.small,
                                       self.LCG+self.small)
        self.hullconstraints.AC_revise(print_=False)
        return
    def set_SAC_fwd_Xc(self, frac=2./3.):
        self.SAC_fwd_Xc = self.hullconstraints.get_val(
                            self.hullconstraints.SAC_fwd_Xc).getpoint(frac)
        self.hullconstraints.SAC_fwd_Xc = ia(self.SAC_fwd_Xc-self.small,
                                             self.SAC_fwd_Xc+self.small)
        self.hullconstraints.AC_revise(print_=False)
        return
    def set_SAC_mid_Xc(self, frac=.5):
        self.SAC_mid_Xc = self.hullconstraints.get_val(
                            self.hullconstraints.SAC_mid_Xc).getpoint(frac)
        self.hullconstraints.SAC_mid_Xc = ia(self.SAC_mid_Xc-self.small,
                                       self.SAC_mid_Xc+self.small)
        self.hullconstraints.AC_revise(print_=False)
        return
    def set_SAC_run_Xc(self, frac=1./3.):
        self.SAC_run_Xc = self.hullconstraints.get_val(
                            self.hullconstraints.SAC_run_Xc).getpoint(frac)
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
        self.hullconstraints.SAC_run_len = ia(self.SAC_run_len-self.small,
                                              self.SAC_run_len+self.small)
        self.hullconstraints.AC_revise(print_=False)
        return

    def get_flat_of_center_plane(self, frac = .5):
        """
        """
        self.flat_bottom = self.hullconstraints.state(
                                self.hullconstraints.lfcp).getpoint(frac)
        self.hullconstraints.lfcp = ia(self.flat_bottom-self.small,
                                       self.flat_bottom+self.small)
        self.hullconstraints.AC_revise(print_=False)
        return


    def definestations(self,cp=.5):
        """
            start and end longitudinal vectors on the body
            start and end longitudinal parameters on ???
        """
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
        cpfrac = 1.-.5*(self.SAC_entrance_len/self.SAC_run_len)
        #self.LCG = = cpfrac/2.
        self.stations       = Stations(lwl=self.LengthDWL,cp=cpfrac)
        self.stations.FOSAC = [self.flat_SAC, cp]
        #self.hullconstraints.SAC_entrance_len = ia(self.stations.FOSAC[0]-self.small,
        #                                           self.stations.FOSAC[0]+self.small)

        self.stations.FOWL  = [self.flat_wL, cp]
        self.stations.FOS   = [self.FOS, cp]
        self.stations.FOCP  = [self.flat_bottom,cp]


        lbow = (1./6.)*(self.stations.FOCP[0])
        self.stations.Bow_Profile = [0.,lbow]

        locbowfair  = (1./3.)*(self.LengthDWL - self.stations.FOCP[1])
        self.stations.Bow_Fairing = [locbowfair,locbowfair]

        lstern = (1./6.)*(self.LengthDWL - self.stations.FOCP[1])
        self.stations.Stern_Profile = [self.LengthDWL-lstern,
                                       self.LengthDWL]
        locsternfair = (1./3.)*(self.stations.FOCP[0])
        self.stations.Stern_Fairing = [self.LengthDWL-locsternfair,
                                       self.LengthDWL-locsternfair]
        return



    def compute_midship_section_old(self):
        """For reference only
        """
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
        print 'compute_midship_section'
        k = self.k
        nump=self.nump
        area = self.maxSecArea/2.
        start, end, a, c = self.stations.FOSAC
        mid = (start+end)/2.
        b = (a+c)/2.
        stations = [start, mid, end]
        fractions = [a,b,c]
        Bmax = self.MaxBreadth/2.
        Dmax = self.MaxDraft
        self.mshp_sections2D = []
        self.mshp_sections   = []
        self.Lsplines.midSS  = []
        #loc=stations[0]
        #frac = fractions[0]
        for loc, frac in zip(stations,fractions):
            #b = self.DWL.FindPoint(loc,2)
            #a = self.CProfile.FindPoint(loc,2) #reversed!
            o,a,b = self.get_transverse_start_end_param_loc(loc,loc)
            #o = self.SAC.FindPoint(loc_start,sacwhich)

            assert(a[2] is True), "error: Could not locate midship match point on Keel Profile Curve"
            assert(b[2] is True), "error: Could not locate midship match point on DWL"

            #pt2 = self.CProfile.CurvePoint(a[0]) #start
            #pt1 = self.DWL.CurvePoint(b[0]) #end
            pt1 = b[1]
            pt2 = a[1]
            #pt2,pt1=self.get_transverse_start_end_points(a[0],b[0])
            assert(abs(Dmax-pt1[1])<self.tol), "MaxDraft does not match Profile: {} != {}".format(Dmax,pt1[1])
            if a[2] is False or b[2] is False or (abs(pt1[0]-Bmax)>self.tol):
                print 'warning, sections may not match global max!'
#                bottom  = linear_vertices((0.,0.,loc),
#                                          (Bmax,0.,loc),4)
#                side    = linear_vertices((Bmax,0.,loc),
#                                          (Bmax,self.MaxDraft,loc),3)
            bottom  = linear_vertices((pt2[0],pt2[1],loc),
                                      (pt1[0],pt2[1],loc),4)
            side    = linear_vertices((pt1[0],pt2[1],loc),
                                      (pt1[0],pt1[1],loc),4)
            section = np.zeros((len(bottom)+len(side),3),float)
            section[:len(bottom)]=bottom
            section[len(bottom):]=side
            section = spline.Bspline(section[:,0:2],self.k,self.nump)
            ab = np.rad2deg(section.compute_tangent(0.))
            ae = np.rad2deg(section.compute_tangent(1.))
            #
            #optimize based on DWL, SAC, and Body cL_Profile:
            #
            """TODO: add bilge radius constraint"""
            """consider its automation with lim Amsh -> 1"""
            interval_data, small = interval_bounds(section)
            FPD = FormParameterDict(section)
            FPD.add_Area_to_y_Constraint(kind='equality',
                               value = area)
            """TODO: There is something in the CPKeel
            that is leaving the mad depth during
            the midship stations...
            """
            FPD.add_AngleConstraint(kind='LS',
                                    location = 0.,
                                    value = ab)
            FPD.add_verticalAngleConstraint(kind='equality',
                                            location = 1.,
                                            value = 90.-ae)
            # C2 continuity for flat side and bottom
            FPD.add_CurvatureConstraint(kind='equality',
                                    location = 0.,
                                    value = 0.)
            FPD.add_CurvatureConstraint(kind='equality',
                                    location = 1.0,
                                    value = 0.)
            FPD.add_relative_yVertexConstraint(
                                    kind = 'equality', location=None,
                                    value = 1., weight = 1.0,
                                    index = 2, index2 = 3,
                                    seperation = 0. )
            FPD.add_relative_xVertexConstraint(
                                    kind = 'equality', location=None,
                                    value = 1., weight = 1.0,
                                    index = 4, index2 = 5,
                                    seperation = 0. )
            FPD.add_relative_xVertexConstraint(
                                    kind = 'equality', location=None,
                                    value = 1., weight = 1.0,
                                    index = 3, index2 = 4,
                                    seperation = 0. )
            #            FPD.add_relative_yVertexConstraint(
            #                                    kind = 'min', location=None,
            #                                    value = 1., weight = 1.0,
            #                                    index = 4, index2 = 3,
            #                                    seperation = 0. )
            #
            #----------------------- fairness:
            #
            FPD.add_E1(kind='LS', weight = .5)
            FPD.add_E2(kind='LS', weight = .5)
            #FPD.add_E3(kind='LS', weight = .5)
            FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data,
                                                     small, 1.e4)
            Lspline = IntervalLagrangeSpline(section, L, data = interval_data)
            Lspline.optimize()
            #
            # ----Done
            #
            self.Lsplines.midSS.append(Lspline)
            #
            #
            section = Lspline.curve
            self.mshp_sections2D.append(section)
            vts = frame.add1D(section.vertices,2,where=loc)
            section = spline.Bspline(vts, k, nump)
            self.mshp_sections.append(section)
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
        self.SAC_entrance_len = self.hullconstraints.state(
                                self.hullconstraints.SAC_entrance_len).getpoint(frac)
        self.hullconstraints.SAC_entrance_len = ia(self.SAC_entrance_len-self.small,
                                     self.SAC_entrance_len+self.small)
        self.hullconstraints.AC_revise(print_=False) #must run this now to ensure consistency fore n aft
        self.SAC_run_len = self.hullconstraints.state(
                                self.hullconstraints.SAC_run_len).getpoint(frac)
        self.hullconstraints.SAC_run_len = ia(self.SAC_run_len-self.small,
                                     self.SAC_run_len+self.small)

        self.hullconstraints.AC_revise(print_=False)
        return

    def compute_stern_fairing_curve(self):
        """stations.Stern_Fairing
        Should be more correct than bow
        """
        print 'compute_stern_fairing_curve'
        start,end,ao,bo = self.stations.Stern_Fairing

        o,a,b = self.get_transverse_start_end_param_loc(start,end)
        #if a[0]==b[0]:
        if abs(a[1][2]-b[1][2])<self.tol:
            pass
        else:
            print 'stern fairning is not plainer'
        area = self.Asfc/2. #self.SAC.CurvePoint(o[0])[1]/2.
        assert(a[2] is True), "error: Could not locate stern fairing match point on Keel Profile Curve"
        assert(b[2] is True), "error: Could not locate stern fairing match point on DWL"
        pt1 = a[1]
        pt2 = b[1]
        bottom  = linear_vertices((pt1[0],pt1[1],pt1[2]),
                                      (pt2[0],pt1[1],pt1[2]),4)
        side    = linear_vertices((pt2[0],pt1[1],pt2[2]),
                                      (pt2[0],pt2[1],pt2[2]),4)
        section = np.zeros((len(bottom)+len(side),3),float)
        section[:len(bottom)]=bottom
        section[len(bottom):]=side
        section = spline.Bspline(section[:,0:2],self.k,self.nump)
        ab = np.rad2deg(section.compute_tangent(0.))
        ae = np.rad2deg(section.compute_tangent(1.))
        interval_data, small = interval_bounds(section)
        FPD = FormParameterDict(section)
        FPD.add_Area_to_y_Constraint(kind='equality',
                                     value = area)
        FPD.add_AngleConstraint(kind='equality',
                                location = 0.,
                                value = ab)
        FPD.add_verticalAngleConstraint(kind='equality',
                                        location = 1.,
                                        value = 90.-ae)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        #FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data,
                                                 small, 1.e4)
        Lspline = IntervalLagrangeSpline(section, L, data = interval_data)
        Lspline.optimize()
        #
        # ----Done
        #
        self.Lsplines.sternfairing = Lspline
        #
        #
        section = Lspline.curve
        #self.bowfairning.append(section)
        vts = frame.add1D(section.vertices,2,
                          where=pt1[2])
        section = spline.Bspline(vts, self.k, self.nump)
        self.sternfairing = section
        return
    def compute_bow_fairing_curve(self):
        """stations.Bow_Profile"""

        print 'compute_bow_fairing_curve'
        start,end,ao,bo = self.stations.Bow_Fairing
        o,a,b = self.get_transverse_start_end_param_loc(start,end)

        #o = self.SAC.FindPoint(start,0)
        #assert(abs(o[1][1]-self.Abfc)<self.tol)

        if abs(a[1][2]-b[1][2])<self.tol:
            pass
        else:
            print 'bow fairning is not plainer'
        area = self.Abfc/2.#self.SAC.CurvePoint(o[0])[1]/2. #ao)/2.#ao doesnt work as SAC has been reparameterized
        assert(a[2] is True), "error: Could not locate bow fairing match point on Keel Profile Curve"
        assert(b[2] is True), "error: Could not locate bow fairing match point on DWL"
        pt1 = a[1]
        pt2 = b[1]
        bottom  = linear_vertices((pt1[0],pt1[1],pt1[2]),
                                      (pt2[0],pt1[1],pt1[2]),4)
        side    = linear_vertices((pt2[0],pt1[1],pt2[2]),
                                      (pt2[0],pt2[1],pt2[2]),4)
        section = np.zeros((len(bottom)+len(side),3),float)
        section[:len(bottom)]=bottom
        section[len(bottom):]=side
        section = spline.Bspline(section[:,0:2],self.k,self.nump)
        ab = np.rad2deg(section.compute_tangent(0.))
        ae = np.rad2deg(section.compute_tangent(1.))
        interval_data, small = interval_bounds(section)
        FPD = FormParameterDict(section)
        FPD.add_Area_to_y_Constraint(kind='equality',
                                     value = area)
        FPD.add_AngleConstraint(kind='equality',
                                location = 0.,
                                value = ab)
        FPD.add_verticalAngleConstraint(kind='equality',
                                        location = 1.,
                                        value = 90.-ae)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        #FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data,
                                                 small, 1.e4)
        Lspline = IntervalLagrangeSpline(section, L, data = interval_data)
        Lspline.optimize()
        #
        # ----Done
        #
        self.Lsplines.bowfairning = Lspline
        #
        #
        section = Lspline.curve
        #self.bowfairning.append(section)
        vts = frame.add1D(section.vertices,2,
                          where=pt1[2])
        section = spline.Bspline(vts, self.k, self.nump)
        self.bowfairning = section
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
        self.hullconstraints.bbfc = ia(self.bbfc-self.small,
                                       self.bbfc+self.small)
        self.hullconstraints.AC_revise()
        return

    def bow_fairness_section_draft(self, frac=.5):
        """get the draft at the bow fairness control curve
        """
        self.dbfc = self.hullconstraints.state(
                self.hullconstraints.dbfc).getpoint(frac)
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
        self.hullconstraints.Cbfc = ia(self.Cbfc-self.small,
                                       self.Cbfc+self.small)
        self.hullconstraints.AC_revise()
        return

    def bow_fairness_section_area(self, frac=.5):
        """get the area of the bow fairness control beam
        """
        Abfc = self.hullconstraints.state(
                self.hullconstraints.Abfc).getpoint(frac)
#        Cbfc = self.hullconstraints.state(
#                self.hullconstraints.Cbfc).getpoint(frac)
#        self.hullconstraints.Abfc = ia(Abfc-self.small,
#                                       Abfc+self.small)
        self.Abfc = Abfc #* Cbfc
        self.hullconstraints.AC_revise()
        return

    def stern_fairness_section_beam(self, frac=.5):
        """get the beam at the stern fairness control curve
        """
        self.bsfc = self.hullconstraints.state(
                self.hullconstraints.bsfc).getpoint(frac)
        self.hullconstraints.bsfc = ia(self.bsfc-self.small,
                                       self.bsfc+self.small)
        self.hullconstraints.AC_revise()
        return

    def stern_fairness_section_draft(self, frac=.5):
        """get the draft at the stern fairness control curve
        """
        self.dsfc = self.hullconstraints.state(
                self.hullconstraints.dsfc).getpoint(frac)
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
        self.hullconstraints.Csfc = ia(self.Csfc-self.small,
                                       self.Csfc+self.small)
        self.hullconstraints.AC_revise()
        return

    def stern_fairness_section_area(self, frac=.5):
        """get the area of the stern fairness control beam
        """
        Asfc = self.hullconstraints.state(
                self.hullconstraints.Asfc).getpoint(frac)
#        self.hullconstraints.Asfc = ia(Asfc-self.small,
#                                       Asfc+self.small)
#        Csfc = self.hullconstraints.state(
#                self.hullconstraints.Csfc).getpoint(frac)
        self.Asfc = Asfc#*Csfc
        self.hullconstraints.AC_revise()
        return

    def compute_SAC_new(self, LocMaxSection=.5, flat=False):
        sac = SAC(Amax = self.maxSecArea,
                  l1 =self.SAC_entrance_len,
                  l2=self.flat_SAC,
                  l3=self.SAC_run_len,
                  Xe = self.SAC_fwd_Xc,
                  Xm = self.SAC_mid_Xc,
                  Xr = self.SAC_run_Xc,
                  Abfc = self.Abfc,
                  Abfc_loc = None,
                  Asfc = self.Asfc,
                  Asfc_loc = None,
                  stations = self.stations)
        self.SAC = sac.SAC
        self.Lsplines.SAC_entrance=sac.entrance
        self.Lsplines.SAC_entrance=sac.midbody
        self.Lsplines.SAC_entrance=sac.run
        return

    def compute_SAC(self, LocMaxSection=.5,
                    Xc=None, flat=False):
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
        Cab_given = 0.
        Cae_given = 0.
        slope = 'down'

        ae = alphae
        ab = alphab
        curve = initial_curve((xb,yb),
                              (xe,ye),
                              num=12, k=k,nump=nump)
        v = copy.deepcopy(curve.vertices)
        v[4:8,1] = self.maxSecArea
        v[3,1] = self.maxSecArea*2./3.
        v[8,1] = self.maxSecArea*2./3.

        curve.vertices = v #updates the curve automagically
        interval_data, small = interval_bounds(curve)



        FPD = FormParameterDict(curve)
        FPD.add_AreaConstraint(kind='equality',
                               value = self.BareHullVolume)
#        if Xc is None:
#            FPD.add_XcConstraint(kind='equality', value=self.LengthDWL/2.)
#        else:
#            FPD.add_XcConstraint(kind='equality', value=Xc)
        if flat:
            FPD.add_yPointConstraint(kind='equality',
                                 location = 0.3333333333,
                                 value = self.maxSecArea)#*2./3.)

            FPD.add_yPointConstraint(kind='equality',
                                 location = 0.6666666666,
                                 value = self.maxSecArea)#*2./3.)
        #if self.fairbow:
        start,end,ao,bo = self.stations.Bow_Fairing
        AreaBowFairingCurve = self.Abfc#/2.
        where = ao
        FPD.add_yPointConstraint(kind='equality',
                                 location = where,
                                 value = AreaBowFairingCurve)
        #else:
        #    pass
        FPD.add_AngleConstraint(kind='LS',
                                location = 0.,
                                value = alphab)
        FPD.add_AngleConstraint(kind='LS',
                                location = 1.,
                                value = alphae)

        FPD.add_yPointConstraint(kind='equality',
                                 location = LocMaxSection,
                                 value = self.maxSecArea)

        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        #FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
        self.Lsplines.SAC=Lspline
        self.SAC = Lspline.curve

        if flat:
            #Lspline.curve = self.Lsplines.SAC
            curve = Lspline.curve
            #curve.knot_insertion(.35)
            #curve.knot_insertion(.65)
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve)
            FPD.add_AreaConstraint(kind='equality',
                                   value = self.BareHullVolume)
            #            if Xc is None:
            #                FPD.add_XcConstraint(kind='equality', value=self.LengthDWL/2.)
            #            else:
            #                FPD.add_XcConstraint(kind='equality', value=Xc)

            #            FPD.add_AngleConstraint(kind='LS',
            #                                    location = 0.,
            #                                    value = alphab)
            #            FPD.add_AngleConstraint(kind='LS',
            #                                    location = 1.,
            #                                    value = alphae)


            #if self.fairbow:
            start,end,ao,bo = self.stations.Bow_Fairing
            AreaBowFairingCurve = self.Abfc#/2.
            where = ao
            FPD.add_yPointConstraint(kind='equality',
                                     location = where,
                                     value = AreaBowFairingCurve)
            #else:
            #    pass
            FPD.add_yPointConstraint(kind='equality',
                                     location = 0.3333333333,
                                     value = self.maxSecArea)
            FPD.add_yPointConstraint(kind='equality',
                                     location = LocMaxSection,
                                     value = self.maxSecArea)
            FPD.add_yPointConstraint(kind='equality',
                                     location = 0.6666666666,
                                     value = self.maxSecArea)

            FPD.add_AngleConstraint(kind='equality',
                                    location = .3333333333,
                                    value = 0.)
            FPD.add_AngleConstraint(kind='equality',
                                    location = .6666666666,
                                    value = 0.)
            s1 = self.stations.FOSAC[0]
            s2 = self.stations.FOSAC[1]
            FPD.add_xPointConstraint(kind='equality',
                                     location = 0.3333333333,
                                     value = s1)
            FPD.add_xPointConstraint(kind='equality',
                                     location = 0.6666666666,
                                     value = s2)

            FPD.add_E1(kind='LS', weight = .5)
            FPD.add_E2(kind='LS', weight = .5)
            #FPD.add_E3(kind='LS', weight = .5)
            FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            Lspline.optimize(stop = 50)
            self.Lsplines.SAC=Lspline
            self.SAC = Lspline.curve
        verts = self.SAC.vertices
        vnew = frame.xy_to_zyx_elevate(verts,elevation=0.)
        self.SAC_3dvis = spline.Bspline(vnew, k, nump)
        return


    def compute_DWL(self, LocMaxSection=.5, Xc=None,
                    flat=False, height=None, alphab=0.,alphae=0.):
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
                              num=13, k=k,nump=nump)
        v = copy.deepcopy(curve.vertices)
        v[4:8,1] = Bmax
        v[3,1] = Bmax*2./3.
        v[8,1] = Bmax*2./3.

        curve.vertices = v #updates the curve automagically
        interval_data, small = interval_bounds(curve)

        FPD = FormParameterDict(curve)
        #
        #----------------------- area
        #
        FPD.add_AreaConstraint(kind='LS',
                               value = area)
        #
        #----------------------- Xc
        #
        #if Xc is None:
        #    FPD.add_XcConstraint(kind='equality', value=self.LengthDWL/2.)
        #else:
        #    FPD.add_XcConstraint(kind='equality', value=Xc)
        #
        #----------------------- prepare for flat of DWL:
        #
        #start,end,ao,bo = self.stations.FOWL
        start, end, ao,bo= self.stations.FOSAC
        if flat:
            FPD.add_yPointConstraint(kind='equality',
                     location = ao,
                     value = Bmax)
            FPD.add_yPointConstraint(kind='equality',
                     location = bo,
                     value = Bmax)
        #problomatic center bulge:---------------
        mid = (start+end)/2.
        b = (ao+bo)/2.
        #FPD.add_yPointConstraint(kind='equality',
        #                         location = LocMaxSection,
        #                         value = Bmax)
        FPD.add_yPointConstraint(kind='equality',
                                 location = b,
                                 value = Bmax)
        #----------------------------------------
        #
        #----------------------- curved dwl bow fairness curve attachment:
        #
        start,end,ao,bo = self.stations.Bow_Fairing
        where = bo
        FPD.add_yPointConstraint(kind='equality',
                                 location = where,
                                 value = self.bbfc/2.)
        FPD.add_xPointConstraint(kind='equality',
                                 location = where,
                                 value = end)
        #
        #----------------------- curved dwl stern fairness curve attachment:
        #
        start,end,ao,bo = self.stations.Stern_Fairing
        where = bo
        FPD.add_yPointConstraint(kind='equality',
                                 location = where,
                                 value = self.bsfc/2.)
        FPD.add_xPointConstraint(kind='equality',
                                 location = where,
                                 value = end)
        #
        #-----------------------  Angles at the ends
        #
        FPD.add_AngleConstraint(kind='LS',
                                location = 0.,
                                value = alphab)
        FPD.add_AngleConstraint(kind='LS',
                                location = 1.,
                                value = alphae)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        #FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        # make curve
        Lspline.optimize()
        self.Lsplines.DWL=Lspline
        self.DWL2D = Lspline.curve
#
#---------------------------------------------- Make Flat DWL:
#
        if flat:
            curve = Lspline.curve
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve)
            #
            #-----------------------  area
            #
            FPD.add_AreaConstraint(kind='LS',
                                   value = area)
            #
            #-----------------------  Xc
            #
            #if Xc is None:
            #    FPD.add_XcConstraint(kind='equality', value=self.LengthDWL/2.)
            #else:
            #   wrong to expect this: FPD.add_XcConstraint(kind='equality', value=self.LCG)
            #
            #-----------------------  Angles at the ends
            FPD.add_AngleConstraint(kind='LS',
                                    location = 0.,
                                    value = alphab)#,
                                    #weight = 10.)
            FPD.add_AngleConstraint(kind='LS',
                                    location = 1.,
                                    value = alphae)#,
                                    #weight = 10.)
            #
            #----------------------- flat dwl bow fairness curve:
            #
            start,end,ao,bo = self.stations.Bow_Fairing
            where = bo
            FPD.add_yPointConstraint(kind='equality',
                         location = where,
                         value = self.bbfc/2.)
            FPD.add_xPointConstraint(kind='equality',
                         location = where,
                         value = end)
            #
            #----------------------- Design the DWL Flat Section:
            #"""more economical to split the DWL"""
            #start,end,ao,bo = self.stations.FOWL
            #st_sac, en_sac, ao_sac,bo_sac= self.stations.FOSAC
            start,end,ao,bo = self.stations.FOSAC
            mid = (start+end)/2.
            b = (ao+bo)/2.
            #xpts
            FPD.add_xPointConstraint(kind='equality',
                                     location = ao,
                                     value = start)
            FPD.add_xPointConstraint(kind='equality',
                                     location = bo,
                                     value = end)
            #ypts
            FPD.add_yPointConstraint(kind='equality',
                                     location = ao,
                                     value = Bmax)
            #problomatic center bulge:---------------
            #
            #FPD.add_yPointConstraint(kind='equality',
            #                         location = LocMaxSection,
            #                         value = Bmax)
            FPD.add_yPointConstraint(kind='equality',
                                     location = b,
                                     value = Bmax)
            #----------------------------------------
            FPD.add_yPointConstraint(kind='equality',
                                     location = bo,
                                     value = Bmax)
            #angles
            FPD.add_AngleConstraint(kind='equality',
                                    location = ao,
                                    value = 0.)
            FPD.add_AngleConstraint(kind='equality',
                                    location = bo,
                                    value = 0.)
            #curvature, could replace this with points
            FPD.add_CurvatureConstraint(kind='LS',
                                    location = ao,
                                    value = 0.)
            FPD.add_CurvatureConstraint(kind='LS',
                                    location = bo,
                                    value = 0.)
            #
            #----------------------- curved dwl stern fairness curve attachment:
            #
            start,end,ao,bo = self.stations.Stern_Fairing
            where = bo
            FPD.add_yPointConstraint(kind='equality',
                                     location = where,
                                     value = self.bsfc/2.)
            FPD.add_xPointConstraint(kind='equality',
                                     location = where,
                                     value = end)
            #
            #-----------------------
            #
            FPD.add_E1(kind='LS', weight = .5)
            FPD.add_E2(kind='LS', weight = .5)
            #FPD.add_E3(kind='LS', weight = .5)
            FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            Lspline.optimize(stop = 50)
            self.Lsplines.DWL=Lspline
            self.DWL2D = Lspline.curve
        #---------------------------------
        #rotate to y-z plane and move x=DWL
        verts = Lspline.curve.vertices
        vnew = frame.xy_to_zxy_elevate( verts,
                                            elevation=elevation)
        self.DWL = spline.Bspline(vnew, k, nump)
        return

    def compute_cLProfile(self, LocMaxSection=.5, Xc=None,
                          flat=False, height=None):
        print '\nKeel Profile'
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
        alphab = 85.##35
        alphae = -85.##-35
        Cab_given = 0.
        Cae_given = 0.
        slope = 'down'

        Dmax = self.MaxDraft
        area = self.CLarea

        ab = alphab
        ae = alphae

        curve = initial_curve((xb,yb),
                              (xe,ye),
                              num=13, k=k, nump=nump)
        v = copy.deepcopy(curve.vertices)
        v[4:8,1] = Dmax
        v[3,1] = Dmax*2./3.
        v[8,1] = Dmax*2./3.

        curve.vertices = v #updates the curve automagically
        interval_data, small = interval_bounds(curve)


        FPD = FormParameterDict(curve)
        #
        #
        #----------------------- area
        #
        FPD.add_AreaConstraint(kind='LS',
                               value = area)
        #
        #----------------------- Xc
        #
        #if Xc is None:
        #    FPD.add_XcConstraint(kind='equality', value=self.LengthDWL/2.)
        #else:
        #    FPD.add_XcConstraint(kind='equality', value=Xc)
        #
        #----------------------- cvd cL_Keelprofile bow fairness:
        #
        start,end,ao,bo = self.stations.Bow_Fairing
        where = ao
        FPD.add_yPointConstraint(kind='equality',
                                 location = where,
                                 value = self.dbfc)
        FPD.add_xPointConstraint(kind='equality',
                                 location = where,
                                 value = start)
        #
        #----------------------- end point angles
        #

        FPD.add_AngleConstraint(kind='LS',
                                location = 0.,
                                value = ab)
        #        FPD.add_verticalAngleConstraint(kind='equality',
        #                                location = 0.,
        #                                value = 0.)
        FPD.add_AngleConstraint(kind='LS',
                                location = 1.,
                                value = ae)
        #        FPD.add_verticalAngleConstraint(kind='equality',
        #                                location = 1.,
        #                                value = -2.)
        #
        #------------------------if doing flat of center plane
        #
        if True:
            #start,end,ao,bo = self.stations.FOCP
            start, end, ao, bo = self.stations.FOSAC
            FPD.add_xPointConstraint(kind='equality',
                                     location = ao,
                                     value = start)
            FPD.add_xPointConstraint(kind='equality',
                                     location = bo,
                                     value = end)
            FPD.add_yPointConstraint(kind='equality',
                                 location = ao,
                                 value = Dmax)
            FPD.add_yPointConstraint(kind='equality',
                                 location = bo,
                                 value = Dmax)
        else:
            pass
        #
        #----------------------- Max section
        #
        FPD.add_yPointConstraint(kind='equality',
                                 location = LocMaxSection,
                                 value = Dmax)
        #
        #----------------------- curvd cL_Keelprofile stern fairness curve attachment:
        #
        start,end,ao,bo = self.stations.Stern_Fairing
        where = ao
        FPD.add_yPointConstraint(kind='equality',
                                 location = where,
                                 value = self.dsfc)
        FPD.add_xPointConstraint(kind='equality',
                                 location = where,
                                 value = start)
        #
        #----------------------- Fairness Functional:
        #
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        #FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        #
        #----------------------- Setup and Solve:
        #
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
        #save
        self.Lsplines.CProfile=Lspline
        self.CProfile2D = Lspline.curve
#
#---------------------------------------------- Make Flat Keel Cp Curve:
#
        if flat:
            curve = Lspline.curve
            interval_data, small = interval_bounds(curve)
            FPD = FormParameterDict(curve)
            #
            #
            #----------------------- area
            #
            FPD.add_AreaConstraint(kind='LS',
                                   value = area)
            #
            #----------------------- Xc
            #
            #if Xc is None:
            #    FPD.add_XcConstraint(kind='equality', value=self.LengthDWL/2.)
            #else:
            #    FPD.add_XcConstraint(kind='equality', value=Xc)
            #
            #----------------------- flat cLprofile bow fairness:
            #
            start,end,ao,bo = self.stations.Bow_Fairing
            where = ao
            FPD.add_yPointConstraint(kind='equality',
                         location = where,
                         value = self.dbfc)
            FPD.add_xPointConstraint(kind='equality',
                                 location = where,
                                 value = end)
            #
            #------------------- end point angles
            #
            FPD.add_AngleConstraint(kind='LS',
                                    location = 0.,
                                    value = ab)
            #            FPD.add_verticalAngleConstraint(kind='equality',
            #                                    location = 0.,
            #                                    value = alphab)
            FPD.add_AngleConstraint(kind='LS',
                                    location = 1.,
                                    value = ae)
            #            FPD.add_verticalAngleConstraint(kind='equality',
            #                                    location = 1.,
            #                                    value = -2.)
            #
            #------------------- Flat of Center Plane
            #
            #start,end,ao,bo = self.stations.FOCP
            start, end, ao,bo = self.stations.FOSAC
            FPD.add_xPointConstraint(kind='equality',
                                     location = ao,
                                     value = start)
            FPD.add_xPointConstraint(kind='equality',
                                     location = bo,
                                     value = end)
            FPD.add_yPointConstraint(kind='equality',
                                     location = ao,
                                     value = Dmax)
            FPD.add_yPointConstraint(kind='equality',
                                     location = bo,
                                     value = Dmax)
            #
            #----------------------- Max section
            #
            FPD.add_yPointConstraint(kind='equality',
                                     location = LocMaxSection,
                                     value = Dmax)
            # C1 continuity to join the midsection
            FPD.add_AngleConstraint(kind='equality',
                                    location = ao,
                                    value = 0.)
            FPD.add_AngleConstraint(kind='equality',
                                    location = bo,
                                    value = 0.)
            # C2 continuity to join the midsection
            FPD.add_CurvatureConstraint(kind='equality',
                                    location = ao,
                                    value = 0.)
            FPD.add_CurvatureConstraint(kind='equality',
                                    location = bo,
                                    value = 0.)
            #
            #----------------------- cvd cL_Keelprofile stern fairness curve attachment:
            #
            start,end,ao,bo = self.stations.Stern_Fairing
            where = ao
            FPD.add_yPointConstraint(kind='equality',
                                     location = where,
                                     value = self.dsfc)
            FPD.add_xPointConstraint(kind='equality',
                                     location = where,
                                     value = start)
            #
            #----------------------- Fairness Functional:
            #
            FPD.add_E1(kind='LS', weight = .5)
            FPD.add_E2(kind='LS', weight = .5)
            FPD.add_E3(kind='LS', weight = .5)
            FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            #
            #----------------------- Setup and Solve:
            #
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            Lspline.optimize(stop = 50)
            self.Lsplines.CProfile=Lspline
            self.CProfile2D = Lspline.curve

        #rotate to y-z plane and move x=DWL
        verts = Lspline.curve.vertices
        #vnew = frame.xy_to_zyx_elevate( verts,
        #                                    elevation=elevation)
        vnew = frame.xy_to_zyx_elevate_noflp( verts,
                                            elevation=elevation)
        self.CProfile = spline.Bspline(vnew, k, nump)
        return

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
        """mirrored
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

        #"""
        self.DWL.plot3DmultiList([self.DWL,
                                  self.CProfile],
                                  [])
        self.DWL.plot3DmultiList([],
                                  self.mshp_sections)

        self.DWL.plot3DmultiList([self.DWL,
                                  self.CProfile],
                                  self.mshp_sections)
        return


    def insert_ctrl_vertices(self, curve, n):
        knots = np.linspace(0.,1.,n+2)[1:-1]
        for knot in knots:
            curve.knot_insertion(knot)
        return curve

    def split_hull_curves(self):
        """Split the ends of the CP Curve
        and DWL
        Make these fwd and aft ends have == # vertices
        as the mid ship sections
        """
        k=self.k
        nump = self.nump
        #
        dwl = self.DWL
        cpk = self.CProfile
        #
        bb, be = self.stations.Bow_Fairing[2:]
        sb, se = self.stations.Stern_Fairing[2:]
        #
        #split:
        self.CPK_fwd, self.CPK_cdr_a = cpk.CurveSplit(bb)
        self.DWL_fwd, self.DWL_cdr_a = dwl.CurveSplit(be)
        #
        self.CPK_cdr_b, self.CPK_aft = self.CPK_cdr_a.CurveSplit(sb)
        self.DWL_cdr_b, self.DWL_aft = self.DWL_cdr_a.CurveSplit(se)
        #
        view = False
        if view:
            self.CPK_aft.plot3DmultiList([self.CPK_fwd,
                                          self.CPK_cdr_b,self.CPK_aft],
                                         [self.DWL_fwd,
                                          self.DWL_cdr_b,self.DWL_aft])
            self.CPK_aft.plot3DmultiList([self.DWL_fwd,
                                          self.DWL_cdr_b,self.DWL_aft],
                                          [self.CPK_fwd,
                                          self.CPK_cdr_b,self.CPK_aft])
            self.CPK_aft.plot3DmultiList([self.CPK_fwd,self.CPK_aft],
                                         [self.DWL_fwd,self.DWL_aft])
        #
        lM = self.mshp_sections[0].n
        lA = self.CPK_aft.n
        lF = self.CPK_fwd.n
        self.CPK_aft = self.insert_ctrl_vertices(self.CPK_aft, lM-lF)
        self.CPK_fwd = self.insert_ctrl_vertices(self.CPK_fwd, lM-lA)
        self.CPK_fwd = spline.Bspline(self.CPK_fwd.vertices[::-1],
                                      k,nump)
        assert(self.CPK_fwd.n==self.CPK_aft.n),"CPK Split Curves Do not match!"
        return

    def correct_hull_curves(self):
        
        diff = self.CProfile.n-self.DWL.n
        assert(diff==0),"Complication, DWL.n != CProfile.n"
        if not (diff==0):
            if self.CProfile.n>self.DWL.n:
                self.DWL = self.insert_ctrl_vertices(self.DWL,diff)
            elif self.CProfile.n<self.DWL.n:
                self.CProfile = self.insert_ctrl_vertices(self.CProfile,diff)
        assert(diff==0),"Error, DWL.n => CProfile.n correction failed!"
        return

    def make_transverse_net(self):
        """
        self.tvertnet = [self.CPK_aft,
                     self.mshp_sections[2],
                     self.mshp_sections[1],
                     self.mshp_sections[0],
                     self.CPK_fwd]
        """
        self.tcurvenet = [self.CPK_aft,
                         self.sternfairing,
                         self.mshp_sections[2],
                         self.mshp_sections[1],
                         self.mshp_sections[0],
                         self.bowfairning,
                         self.CPK_fwd]
        tvertnet = np.asarray([self.CPK_aft.vertices,
                         self.sternfairing.vertices,
                         self.mshp_sections[2].vertices,
                         self.mshp_sections[1].vertices,
                         self.mshp_sections[0].vertices,
                         self.bowfairning.vertices,
                     self.CPK_fwd.vertices])
        return tvertnet

    def issue_longitudinals(self):
        """
        """
        k=self.k
        nump = self.nump
        N = self.CPK_aft.n
        nnew = N-2
        #m = self.CProfile.n
        #lvertices = np.zeros((nnew,m),float)
        tvertnet = self.make_transverse_net()
        lcurvenet = []
        for i in range(1,nnew+1):
            lcurvenet.append(
                spline.interpolatedBspline(
                tvertnet[:,i], k, nump))
        self.lcurvenet = lcurvenet
        #plot = False
        #if plot:
        #    self.CPK_aft.plot3DmultiList(self.tcurvenet,self.lcurvenet)
        return

    def correct_longitudinal_match_hullcurves(self):
        """DWL_cdr_b
        CPK_cdr_b
        """
        ln = self.lcurvenet[0].n
        dwln = self.DWL.n
        assert(dwln == self.CProfile.n),"earlier dwl <=> CPKeel n failed"
        diff = dwln - ln
        assert(diff>0),"DWL.n < lcurve[0].n  ??"
        for curve in self.lcurvenet:
            curve = self.insert_ctrl_vertices(curve, diff)
        self.lcurvenet.insert(0,self.CPK_cdr_b)
        #self.lcurvenet.append(self.DWL_cdr_b)
        self.lcurvenet.append(self.DWL)
        return

    def make_longitudinals(self):
        """interpolate the transverse vertices
        via Piegel and Tiller
        """
        assert(self.sternfairing.n==self.mshp_sections[0].n),"Transvere Curves Do not match!"
        assert(self.sternfairing.n==self.bowfairning.n),"Transvere Curves Do not match!"
        self.split_hull_curves()
        assert(self.sternfairing.n==self.CPK_aft.n),"Transvere Curves Do not match!"
        self.correct_hull_curves() #must come after split because fairness curves assume fixed parameter loc!
        self.issue_longitudinals()
        self.correct_longitudinal_match_hullcurves()
        self.hullsurf = spline.BsplineSurface(self.tcurvenet,
                                              self.lcurvenet)
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

def compute_hull(hulldesign):
    hulldesign.compute_SAC_new(flat=True)
    #hulldesign.compute_SAC(flat=True) #single piece

    hulldesign.compute_cLProfile(flat=True)
    np.rad2deg(hulldesign.CProfile2D.compute_tangent(0.))
    np.rad2deg(hulldesign.CProfile2D.compute_tangent(1.))

    hulldesign.compute_DWL(flat=True, alphab=30.,alphae=-30.)#35.,35. works #25.,25. works#15., -15. works# 0.,0. works#

    self.compute_midship_section()

    hulldesign.compute_bow_fairing_curve()
    hulldesign.compute_stern_fairing_curve()
    hulldesign.make_longitudinals()

    
    
    
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



def hull_gui(hulldesign):
    class wrapper(object):
        def __init__(self, curve):
            self.curve = curve

    HF = hull_frame()

    HF.add_section(hulldesign.Lsplines.bowfairning,
                   hulldesign.bowfairning.vertices[0,2])


    for el, key in zip(hulldesign.Lsplines.midSS,hulldesign.mshp_sections):
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
    
    
    
    h1.draft = ia(DLEN-h1.SMALL,DLEN)

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

    h1.AC_revise(print_=False)
    h1.print_state()

    self = Hull(h1) #crucial for pointers
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

    self.get_len_FOS(.3)#.2 #doesn't do anything yet
    self.get_len_flat_SAC(frac=.1)
    self.get_len_flat_wL(frac=.1)
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




    #"""
    #compute_hull(self)
    #hull_gui(self)
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
