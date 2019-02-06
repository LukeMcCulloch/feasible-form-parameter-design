# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 18:45:16 2015

@author: lukemcculloch

souces:
    http://www.randalolson.com/2014/06/28/
        how-to-make-beautiful-data-visualizations-in-python-with-matplotlib/
        
        
propositions:
    -bulb area can be peeled from within the general SAC area
        *perhaps with provisions made for making sure
        a good shape fit is availible within the basic SAC curve form
    -start and end SAC y heights can map to tractable rule systems
        *including relations with waterline and cp curves
"""
import numpy as np
#from scipy import ndimage
import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid.axislines import SubplotZero
#import matplotlib.patches as patches
#from matplotlib.patches import Ellipse
#import matplotlib.path    as path
mpl.rcParams['legend.handlelength'] = 0

import curve as spline
from interval_arithmetic import ia
from utility_optimization import package_norms
#
from ADILS import interval_bounds
from ADILS import Lagrangian
from ADILS import lagrangian_bounds
from ADILS import IntervalLagrangeSpline
from FormParameter import FormParameterDict
#
import mx_plots as spfuncs # circle, gfunc class, etc...
# 
def default_input( message, defaultVal ):
    """http://stackoverflow.com/
    questions/5403138/how-to-set-a-default-string-for-raw-input
    """
    if defaultVal:
        return raw_input( "%s [%s]:" % (message,defaultVal) ) or defaultVal
    else:
        return raw_input( "%s " % (message) )
        
def check_collinear(vertices):
    v1 = list(vertices[0])
    v2 = list(vertices[1])
    v3 = list(vertices[2])
    v1.insert(0,1.)
    v2.insert(0,1.)
    v3.insert(0,1.)
    test = np.asarray([v1,v2,v3])
    if np.linalg.det(test) ==0:
        return True
    else:
        return False


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

class HullCurve(object):
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

class SAC(object):
    k       = 4
    nump    = 30
    xb      = 0.
    yb      = 0.
    def __init__(self, Amax, #v1,v2,v3, 
                 l1,l2,l3, 
                 Xe, Xm, Xr, Xlcb, 
                 Afp, Aap, Cpe, Cpr):
        #self.V      = vol
        self.Amax   = Amax
        #self.Ve     = v1
        #self.Vm     = v2
        #self.Vr     = v3
        self.Le     = l1
        self.Lm     = l2
        self.Lr     = l3
        self.L      = l1+l2+l3
        self.Xe     = Xe
        self.Xm     = Xm
        self.Xr     = Xr
        self.Xlcb   = Xlcb
        self.Afp    = Afp
        self.Aap    = Aap
        self.Cpe    = Cpe
        self.Cpr    = Cpr
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
        n1 = self.run.curve.n
        n2 = self.midbody.curve.n
        n3 = self.entrance.curve.n
        self.midbody.translate(self.Le+SAC.xb)
        self.entrance.translate(self.Lm + self.Le+SAC.xb)
        nv = n1 + n2 + n3
        
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
        self.midbody.translate(self.Le+SAC.xb)
        self.entrance.translate(self.Lm + self.Le+SAC.xb)
        nv = n1 + n2-2 + n3
        
        Cv = np.zeros((nv,2),float)
        Cv[0:n1,:] = self.run.curve.vertices
        Cv[n1:n1+n2-2,:] = self.midbody.curve.vertices[1:-1]
        Cv[n1+n2-2:n1+n2+n3,:] = self.entrance.curve.vertices
        
        self.SAC = spline.Bspline(Cv, SAC.k, nump = SAC.nump*3)
        return
    
    #def compute_entrance_SAC(self):
    def compute_run_SAC(self):
        k       = SAC.k
        nump    = SAC.nump
        nCV     = 7
        xb      = SAC.xb
        xe      = xb + self.Le
        yb      = SAC.yb
        ye      = self.Amax
        yfp     = self.Afp
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
        FPD.add_AngleConstraint(kind = 'equality', value = 0., location = 0.)
        FPD.add_AngleConstraint(kind = 'equality', value = 0., location = 1.)
        FPD.add_CurvatureConstraint(kind = 'equality', value = 0., location = 0.)
        FPD.add_CurvatureConstraint(kind = 'equality', value = 0., location = 1.)
        FPD.add_XcConstraint(kind = 'equality', value = 8.5)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        
        Lspline.compute_lagrangian()
        Lspline.display(mytitle = 'SAC_run_ini', 
                        x_axis_name = 'x',
                        y_axis_name = 'z')
        Lspline.optimize()
        Lspline.display(mytitle = 'SAC_run', 
                        x_axis_name = 'x',
                        y_axis_name = 'z')
        return Lspline
        
        
    def compute_midbody_SAC(self):
        k       = SAC.k
        nump    = SAC.nump
        n       = 7
        xb = SAC.xb
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
        
        Lspline.compute_lagrangian()
        Lspline.display(mytitle = 'SAC_mid_ini', 
                        x_axis_name = 'x',
                        y_axis_name = 'z')
        Lspline.optimize()
        Lspline.display(mytitle = 'SAC_mid', 
                        x_axis_name = 'x',
                        y_axis_name = 'z')
        return  Lspline
        
    
    #def compute_run_SAC(self):
    def compute_entrance_SAC(self):
        k       = SAC.k
        nump    = SAC.nump
        nCV     = 7
        xb  = SAC.xb
        xe  = xb + self.Lr
        yb  = self.Amax
        ye  = 0.
        yap = self.Aap
        Xc  = self.Xr
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
        FPD.add_AngleConstraint(kind = 'equality', value = 0., location = 1.)
        FPD.add_CurvatureConstraint(kind = 'equality', value = 0., location = 0.)
        FPD.add_CurvatureConstraint(kind = 'equality', value = 0., location = 1.)
        FPD.add_XcConstraint(kind = 'equality', value = 3.5)
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = 1.)
        FPD.add_E3(kind='LS', weight = 1.)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        L = Lagrangian(FPD)
        
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        
        Lspline.compute_lagrangian()
        Lspline.display(mytitle = 'SAC_fwd_ini', 
                        x_axis_name = 'x',
                        y_axis_name = 'z')
        Lspline.optimize()
        Lspline.display(mytitle = 'SAC_fwd', 
                        x_axis_name = 'x',
                        y_axis_name = 'z')
        return  Lspline
        
        
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


class SAC_complex(object):
    def __init__(self, 
                 voldisp, Amsh, lwl, 
                 lfsac, LCG, lfwl, lfcp,
                 Ae,Ab):
        self.k=4
        self.nump=60
        
        self.xb = 0.
        self.xe = lwl
        
        self.lfsac = lfsac
        
        self.alphab = 0.
        self.alphae = 0.
        self.Cab = 0.
        self.Cae = 0.
        
        self.LCG = LCG
        self.yc = 4.
        
        self.maxSecArea = Amsh
        self.Ab = Ab
        self.yb = Ab#0.
        self.ye = Ae
        
        self.BareHullVolume = voldisp
        
        self.slope = 'down'
        
        self.x_offset1 = 3.
        self.x_offset2 = -2.
        self.y_offset = 2.
        
        self.ae = self.alphae
        self.ab = self.alphab
        
        self.LocMaxSection = .5
        
        self.stations = Stations(lwl,.5)
        self.stations.FOSAC = [self.lfsac,.5]
    
    def make(self):
        curve = initial_curve([self.xb,self.yb],
                              [self.xe,self.ye],
                              12,
                              nump=100) 
        interval_data, small = interval_bounds(curve)
        self.Lspline,self.FPD = self.stage1(curve)
        self.Lspline,self.FPD = self.stage2(self.Lspline,
                                            self.FPD)
        return     
        
    def stage1(self, curve, flat=False):
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
                                 location = self.LocMaxSection,
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
        #
        return Lspline,FPD
        
    
        
    def stage2(self, Lspline,FPD):
        #
        #******************************************************************
        # SAC
        # NOW Entering Solver Stage 2
        norms = package_norms(Lspline)
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
        
        #FPD.add_xPointConstraint(kind='max',
        #                         value=0.1,
        #                         location=.001)
        #FPD.add_yPointConstraint(kind='min',
        #                         value=10.,
        #                         location=.001)
        #
        #******************************************************************
        #  FLATS
        #
        #******************************************************************
        # X FLAT SAC 
        #
        #print 'using x - fixity'
        #        FPD.add_xFixity(index=5,
        #                        value=self.stations.FOSAC[0])
        #        FPD.add_xFixity(index=6,
        #                        value=self.stations.FOSAC[1])
        #
        FPD.add_relative_xVertexConstraint(kind = 'equality', location=None, 
                                           value = self.lfsac, weight = 1.0, 
                                           index = 5, index2 = 6, 
                                           seperation = -self.lfsac )
        #        FPD.add_relative_xVertexConstraint(kind = 'min', location=None, 
        #                                           value = 0., weight = 1.0, 
        #                                           index = 5, index2 = 6, 
        #                                           seperation = self.lfsac )
        #        FPD.add_relative_xVertexConstraint(kind = 'max', location=None, 
        #                                           value = self.lfsac+1., weight = 1.0, 
        #                                           index = 5, index2 = 6, 
        #                                           seperation = self.lfsac )
        
        #
        #******************************************************************
        # Y FLAT SAC (3,4,5,6,7)
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
        # 
        #        FPD.add_yPointConstraint(kind = 'equality',
        #                                 location=.01,
        #                                 value = self.Ab,
        #                                 )
        #
        #******************************************************************
        #
        #
        #******************************************************************
        # 
        print 'adding bulb to the SAC functional'
        A_mid = 10.080798436245335
        BulbLength = 4.341062137105493
        #gfunc = spfuncs.gfunc
        #bbow_sac = spfuncs.BBowSAC(Xb=0.,Xe = BulbLength,
        #                           height=A_mid)
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
        """
        gfunc_bbow_SAC  = spfuncs.gfunc(curve,
                                         bbow_sac,
                                         0.,
                                         .00001,
                                         np.linspace(0.,1.,
                                                     bbow_sac.nump*5,
                                                     endpoint=True),
                                        10000.)
        gfunc_bbow_SAC.analysis_type = 'min'
        L.sp['bbow_sac'] = gfunc_bbow_SAC
        #"""
        #
        #******************************************************************
        # 
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
    
if __name__ == "__main__":
    # TLM B-Spline Curve class
    import curve             as     spline
    
    from   initialValues     import InitializeControlPoints, InitializeControlVertices
    import copy
    from   ADILS             import IntervalLagrangeSpline, Lagrangian
    from   FormParameter     import FormParameterDict
    from   initialValues     import interval_bounds, lagrangian_bounds
    option  = 'ini_test'
    test    = None
    test    = 'SAC'  # this one makes the run plot with free stern tube
    do_all  = False
    
    
        
    voldisp = ia(34196.3123868, 34196.3174777).getpoint(.5)
    Amsh    = ia(440.441841483, 440.447841483).getpoint(.5)
    lwl     = ia(113.978485473, 113.984485473).getpoint(.5)
    #
    lfsac   = ia(10.5364794686, 10.5424794686).getpoint(.5)
    LCG     = ia(56.1860734432, 56.1920734432).getpoint(.5)
    #
    lfwl    = ia(14.0303648233, 14.0363648233).getpoint(.5)
    lfcp    = ia(11.3864366677, 11.3924366677).getpoint(.5)
    #
    draft   = ia(18.4490277827, 18.4550277827).getpoint(.5)
    beam    = ia(27.092071177, 27.098071177).getpoint(.5)
    
    BulbVolume  =  ia(7.39104381671, 7.39110381671).getpoint(.5)
    A_mid       =  ia(20.4224589012, 20.4225189012).getpoint(.5)
    A_lateral   =  ia(24.2218639047, 24.2219239047).getpoint(.5)
    BulbLength  =  ia(4.18207557456, 4.18212145051).getpoint(.5)
    BulbDepth   =  ia(6.28644007111, 6.28650007111).getpoint(.5)
    BulbBeam    =  ia(3.26128501745, 3.26133238331).getpoint(.5)
    ##
    ##******************************************************************
    ## Handle aft SAC edge this way:
    drop_aft    = 0.1 #.25
    DROP        = draft*drop_aft
    # aft section area must be 
    #               -greater than a triangle
    #               -less than a box.
    Abox = DROP*beam/2.
    Atri = .5*Abox
    Aaft = Atri+.1*(Abox-Atri)
    ##
    ##******************************************************************
    ## SAC entrance...
    Atube = A_lateral
    
    ##
    ##******************************************************************
    ## 
    if test == 'SAC':
        self = SAC_complex(voldisp, Amsh, lwl, 
                 lfsac, LCG, lfwl, lfcp,
                 Ab=Atube, Ae = Aaft)
        
        self.make()
        
        #
        normit = True
        ax = self.Lspline.curve.plotcurve_detailed(normalize=normit)
        """
        #
        self.Lspline.Lagrangian.sp[
                'bbow_sac'].func.SAC.plotcurve_detailed(
                                normalize=normit,
                                scale_curve=self.Lspline.curve,
                                canvas=ax)
        #
        #Out[24]: <matplotlib.axes._subplots.AxesSubplot at 0x7f780159eb90>
        
        self.Lspline.Lagrangian.sp['bbow_sac'].func
        #Out[25]: <mx_plots.BBowSAC at 0x7f77fe6e87d0>
        
        self.Lspline.Lagrangian.sp['bbow_sac'].func.Lflat
        #Out[26]: 3.6899028165396692
        
        self.Lspline.Lagrangian.sp['bbow_sac'].func.height
        #Out[27]: 10.080798436245335
        
        self.Lspline.Lagrangian.sp['bbow_sac'].func.SAC.vertices
        print 'area = ',self.Lspline.Lagrangian.sp['bbow_sac'].func.area 
        print 'area = ',self.Lspline.Lagrangian.sp['bbow_sac'].func.SAC.area 
        #
        #"""
        """
        #
        #
        args = [self.Lspline.vertices]
        verts   = args[0]
        xpts    = verts[0]
        ypts    = verts[1]
        oldself = self
        self = self.Lspline.Lagrangian.sp['bbow_sac']
        #
        print self.result
        #
        #
        self = oldself
        #
        #"""