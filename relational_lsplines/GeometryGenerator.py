#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 18:09:31 2017

@author: luke

Notes:
    
    
File:
    purpose:
    contains:
"""

import copy
import sqKanren as lp
np = lp.np
import sobol
ia = lp.ia #use ia import from lp instead to match isinstance
#
from  automatic_differentiation import ad as ad #doing this purely to assist in latex table generation
#
from initialValues import InitializeControlVertices
from initialValues import interval_bounds, lagrangian_bounds
from plots         import Plotter
from adials_gui3d  import Quaternion#, frame, DrawCurveInteractive, hull_frame
#
from frames import Frame, DeepVector
import curve as spline
#
from   ADILS         import IntervalLagrangeSpline, Lagrangian
from   FormParameter import FormParameterDict

import utility_optimization as uopt

GraphScaleNode = lp.GraphScaleNode
RulesGraphProcessor = lp.RulesGraphProcessor
#
import random
#
import cPickle as pickle

class LocalSystem(object):
    """Get the local 2D curve system
    with respect to the local 3D curve system
    with respect to the global system
    """
    def __init__(self):
        return
    

class BoundBulbCurves(object):
    def __init__(self,
                 BulbBeam,BulbDepth,BulbLength):
        """use global coordinate directions 
        translate the origin
        """
        self.coords2D = {}
        self.coords2D['x'] = 0
        self.coords2D['y'] = 1
        self.coords2D['z'] = 2
        self.def_rotate_about_x_axis()
        self.def_rotate_about_y_axis()
        self.def_rotate_about_z_axis()
        self.setuprotation()
        
        
        #bulb common coordinate system
        Origin = np.asarray([0.,0.,0.])
        x, y, z = np.identity(3)
        #        self.ccs = Frame( np.asarray([[1.,0.,0.],
        #                                      [0.,1.,0.],
        #                                      [0.,0.,1.]]) )
        self.ccs = Frame( np.asarray([x,y,z]) )
        self.origin = DeepVector(v = Origin,
                                 B=self.ccs)
        
        
        #ccs = self.ccs
    
        #vector: interface point (0.,0.,0.)
        self.v_itfc = DeepVector(v = Origin, 
                                  B = self.ccs)
        q0 = self.v_itfc[0]
        q1 = self.v_itfc[1]
        q2 = self.v_itfc[2]
        
        """
        # INI-BULB CURVE:
        # interface curve:
        """
        #curve defining the interface from the bare hull to the bulb:
        #vector: ini bulb curve start (0.,cL,bL)
        #
        self.v_inibc_cL_bL = DeepVector(
                v = np.asarray([q0,q1,q2]), 
                                  B = self.ccs)
        #vector: mid bulb curve end (mid_length,cL,Top)
        self.v_inibc_cL_tp = DeepVector(
                v = np.asarray([q0,q1,BulbDepth]), 
                                  B = self.ccs)
    
        
        """
        # MID-BULB CURVE, MB(C):
        """
        #curve defining the mid-section of the bulb:
        #vector: mid bulb curve start (mid_length,cL,bL)
        #
        # mbc : mid bulb curve (like mid ship curve)
        #
        self.v_mbc_cL_bL = DeepVector(
                v = np.asarray([BulbLength/2.,q1,q2]), 
                                  B = self.ccs)
        #vector: mid bulb curve end (mid_length,cL,Top)
        self.v_mbc_cL_tp = DeepVector(
                v = np.asarray([BulbLength/2.,q1,BulbDepth]), 
                                  B = self.ccs)
    
        
        """FWD-BULB SHAPE CURVE, FSB(C)
        """
        #curve defining the fwd-shape of the bulb:
        #
        self.v_fsbc_cL_bL = DeepVector(
                v = np.asarray([3.*BulbLength/4.,q1,q2]), 
                                  B = self.ccs)
        #vector: mid bulb curve end (mid_length,cL,Top)
        self.v_fsbc_cL_tp = DeepVector(
                v = np.asarray([3.*BulbLength/4.,q1,BulbDepth]), 
                                  B = self.ccs)
        
        
        """
        #HORIZONTAL "WATER" PLANE CURVE, BWL(C):
        """
        # curve defining the horizontal mid 'water' plane of the bulb:
        #(mid elevation curve)
        #vector: midL bulb curve start (inipt, cL, mid_depth)
        #
        # hbc : horizontal bulb curve (pseudo water plane curve)
        #
        self.v_hbc_cL_bL = DeepVector(
                v = np.asarray([q0,BulbBeam,BulbDepth/2.]), 
                                  B = self.ccs)
        #vector: midL bulb curve end (BulbLength, cL, mid_depth)
        self.v_hbc_cL_tp = DeepVector(
                v = np.asarray([BulbLength,q1,BulbDepth/2.]), 
                                  B = self.ccs)
        
        
        """
        #VERTICAL CENTER PLANE CURVE, CP(C):
        """
        # curve defining the vertical center plane of the bulb:
        #(center plane curve)
        #vector: midL bulb curve start (inipt, cL, mid_depth)
        #
        # cpbc : center plane bulb curve, cL (centerline) to top
        #
        self.v_cpbc_cL_bL = DeepVector(
                v = np.asarray([q0,q1,q2]), 
                                  B = self.ccs)
        #vector: midL bulb curve end (BulbLength, cL, mid_depth)
        self.v_cpbc_cL_tp = DeepVector(
                v = np.asarray([q0,q1,BulbDepth]), 
                                  B = self.ccs)
        
        
        
        
        
        
        self.x_interface    = Origin# [0.,0.,0.]
        self.cLb            = self.x_interface
        self.cLe            = [BulbLength,0.,BulbDepth]
        self.bLb            = [BulbLength/2.,0.,0.]
        self.bLe            = [BulbLength/2.,0.,0.]
    
    def simple_projection(self, r1='x',r2='y'):
        
        return self.coords2D[r1], self.coords2D[r2]
    
    
    def cardinal_rotations(self):
        return
    
    def def_rotate_about_x_axis(self):
        def neg_sin(radians):
            return -np.sin(radians)
        
        def evalx(r):
            r0 = [1.,    0.,       0.]
            r1 = [0., np.cos(r), neg_sin(r)]
            r2 = [0., np.sin(r), np.cos(r)]
            return np.asarray([r0,r1,r2])
        
        self.rotx = evalx
        return
    
    
    
    def def_rotate_about_y_axis(self):
        def neg_sin(radians):
            return -np.sin(radians)
        
        def evaly(r):
            r0 = [np.cos(r), 0., np.sin(r),]
            r1 = [0.,     1.,    0.  ]
            r2 = [neg_sin(r), 0., np.cos(r)]
            return np.asarray([r0,r1,r2])
        
        self.roty = evaly
        return
    
    
    def def_rotate_about_z_axis(self):
        def neg_sin(radians):
            return -np.sin(radians)
        
        def evalz(r):
            r0 = [np.cos(r),  np.sin(r), 0.]
            r1 = [neg_sin(r), np.cos(r), 0.]
            r2 = [0.,     0.,   1.   ]
            return np.asarray([r0,r1,r2])
        
        self.rotz = evalz
        return
    
    def setuprotation(self):
        
        self.alias = {'x':0,'y':1,'z':2}
        
        self.rdptch = {0:self.rotx,
                       1:self.roty,
                       2:self.rotz}
        return
    
    
    def rotation(self, radians, index):
        
        rotation = np.identity((3),float)
        
        if isinstance(index, list):
            for i, rad in zip(index,radians):
                if isinstance(i,str):
                    i = self.alias[i]
                rotation = np.dot(rotation,self.rdptch[i](rad))
        else:
            i = index
            rad = radians
            if isinstance(index,str):
                i = self.alias[i]
            rotation = np.dot(rotation,self.rdptch[i](rad))
                
        return rotation
        
        
    def rot_bulb_geometry(self, tcurvelist,lcurvelist):
        
        return


class FPDIdeal(object):
    """FPD data, rule, and curve class
        input: 
            -form parameters for a curve
            -rules defining its validity
    """
    def __init__(self, xb=0.,yb=12.,xe=12.,ye=0.,
                 alphab=None,alphae=None,
                 valphab=None,valphae=None,
                 Cab_given=None,Cae_given=None,
                 area=None,xc=None,yc=None,
                 ymax=None,
                 nCV = 7, 
                 stationloc=None, base=None,
                 nump=30, k=4,
                 tol = 1.e-7):
        self.xb         = xb
        self.yb         = yb
        self.xe         = xe
        self.ye         = ye
        self.alphab     = alphab
        self.alphae     = alphae
        self.valphab    = valphab
        self.valphae    = valphae
        self.Cab_given  = Cab_given
        self.Cae_given  = Cae_given
        self.area       = area
        self.ymax       = ymax
        self.xc         = xc
        self.yc         = yc
        self.nCV        = nCV
        self.nump       = nump
        self.stationloc = stationloc,
        self.base       = base
        self.k          = k
        self.tol        = tol
        
        if self.xc  is 'default':
            self.xc = self.xe/2.
        
    def make_curve(self):
        self.vertices = InitializeControlVertices(self.xb,
                                                  self.yb,
                                                  self.xe,
                                                  self.ye,
                                                  None,
                                                  None,
                                                  None,
                                                  None,
                                                  self.area,
                                                  self.xc,
                                                  self.yc,
                                                  self.nCV)
        #print 'got vertices:'
        #print self.vertices.vertices
        
        self.curve = spline.Bspline(self.vertices.vertices,
                                    self.k,self.nump)
        return #self.curve
    
    def solve_curve_ini(self):
        """
            STAGE 1
        """
        interval_data, small = interval_bounds(self.curve)
        FPD = self.setup_basic()
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                                    interval_data, small, 1.e4)
        self.Lspline = IntervalLagrangeSpline(self.curve,
                                         L, data = interval_data)
        
        self.Lspline.curve.pts_M_pts()
        self.Lspline.curve.compute_arclength()
        self.Lspline.optimize(stop=30)
        return
    
    def solve_curve_fini(self,Xc=None):
        """
            STAGE 2
        """
        interval_data, small = interval_bounds(self.Lspline.curve)
        
        if Xc is None:
            FPD = self.setup_stage2()
        else:
            FPD = self.setup_stage2(Xc)
        
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                                    interval_data, small, 1.e4)
        self.Lspline = IntervalLagrangeSpline(self.curve,
                                         L, data = interval_data)
        
        self.Lspline.curve.pts_M_pts()
        self.Lspline.curve.compute_arclength()
        
        self.Lspline.optimize(stop=30)
        return
    
        
    def setup_basic(self):
        FPD = FormParameterDict(self.curve)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', 
                                weight = .5)
        if self.area is not None:
            FPD.add_AreaConstraint(kind='equality', 
                               value = self.area)
        
        return FPD
        
    def setup_stage2(self, Xc = None):
        """
        self.alphab
        self.alphae
        self.Cab_given
        self.Cae_given
        """
        interval_data, small = interval_bounds(
                self.Lspline.curve)
        
        E1_0 = copy.deepcopy(self.Lspline.curve.E1)
        E2_0 = copy.deepcopy(self.Lspline.curve.E2)
        E3_0 = copy.deepcopy(self.Lspline.curve.E3)
        S_0  = copy.deepcopy(self.Lspline.curve.AL)
        
        FPD = FormParameterDict(self.curve)
        if Xc is None and self.xc is None:
            pass
        else:
            if Xc is None: Xc = self.xc
            FPD.add_XcConstraint(kind = 'equality', 
                                 value = Xc)
        FPD.add_E1(kind='LS', weight = .5/E1_0.value)
        FPD.add_E2(kind='LS', weight = .5/E2_0.value)
        FPD.add_E3(kind='LS', weight = .5/E3_0.value)
        FPD.add_ArcLengthApprox(kind='LS', 
                                weight = .5/S_0)
        if self.area is not None:
            FPD.add_AreaConstraint(kind='equality', 
                               value = self.area)
        
        """
        TODO: match form curves
        """
        FPD.add_yPointConstraint(kind= 'max', 
                                 value=self.ymax, 
                                 location=0.5 )
        
        if (np.sqrt((self.ymax-self.yb)**2)<self.tol):
            
            dx = (self.xe-self.xb)/2.
            
            FPD.add_yVertexConstraint(kind = 'equality', 
                                      value=self.ymax, 
                                      index =1)
            FPD.add_yVertexConstraint(kind = 'equality', 
                                      value=self.ymax, 
                                      index =2)
            FPD.add_xVertexConstraint(kind = 'equality', 
                                      value=dx, 
                                      index =2)
        if self.alphab is not None:
            FPD.add_AngleConstraint(kind = 'equality',
                                    value = self.alphab,
                                    location = 0.)
            #FPD.add_xFixity
        if self.alphae is not None:
            FPD.add_verticalAngleConstraint(kind = 'equality',
                                    value = self.alphae,
                                    location = 1.)
        ##-------------------------------------------------
        if self.valphab is not None:
            #FPD.add_AngleConstraint(kind = 'equality',
            #                        value = self.valphab,
            #                        location = 0.)
            FPD.add_xFixity(index = 1,
                            value = self.xb)
        if self.valphae is not None:
            #FPD.add_verticalAngleConstraint(kind = 'equality',
            #                        value = self.valphae,
            #                        location = 1.)
            FPD.add_xFixity(index = self.nCV-2,
                            value = self.xe)
        ##-------------------------------------------------
        if self.Cab_given is not None:
            FPD.add_CurvatureConstraint(kind = 'equality',
                                        value = self.Cab_given,
                                        location = 0.)
        if self.Cae_given is not None:
            FPD.add_CurvatureConstraint(kind = 'equality',
                                        value = self.Cae_given,
                                        location = 1.)
            
        
        return FPD
    
    def check_area_to_box_relation(self):
        dx = self.xe-self.xb
        dy1 = self.ye - self.yb
        dy2 = self.ymax - self.yb
        dy = max(dy1,dy2)
        Amax = dx*dy
        assert Amax>self.area, \
             "Area {} > {}feasibility fail".format(Amax,
                                            self.area)
        return  
    
    
    
        
def ini_coeffs(bbobj):#clist, rgp):
    """initialize coefficients to ia(0.,1.)
    """
    print 'initializing coefficients'
    for co in bbobj.Coefficients:
        co   = co == ia(0.,1.)
    
        bbobj.rgp.add_one_rule(co,co.name)
    
        bbobj.rgp.compute_fresh_rules_graph()
    return bbobj


class GeometryGenerator(object):
    """Rules for the Feasible Design
    of a Bulbous Bow
    
    
    bulb coordinate system:
    
        z
        |
        |
        o -----y
       /  
      /
     x  
     
     x : longitudinal
     y : transverse
     z : vertical
     
     
    Notes
    ----------
        New style language based rules and AC_revise
    """
    def __init__(self,
                 ship_beam=None,
                 ship_depth=None,
                 ship_Lpp=None,
                 ship_Amshp=None,
                 ship_Acp=None,
                 ship_Awp=None,
                 ship_Vol=None,
                 bare_hull=None,
                 given_A_mid=None,
                 given_length=None):
        self.nCV = 7
        self.idealdict       = {} #this will eventually store primary bulb curves
        self.init_ship_beam  = ship_beam
        self.init_ship_depth = ship_depth
        self.init_ship_Lpp   = ship_Lpp
        self.init_ship_Amshp = ship_Amshp
        self.init_ship_Acp   = ship_Acp
        self.init_ship_Awp   = ship_Awp
        self.init_ship_Vol   = ship_Vol
        self.bare_hull       = bare_hull # sqKanren States() processed by opt_simple_hull
        
        self.given_A_mid = given_A_mid
        self.given_length = given_length
        #
        self.results_dict = {}
        #
        self.rgp = RulesGraphProcessor(verbose=False)
        self.sobol_seq = sobol.sobolSeq([1,1],[1,1])
        self.tol = 1.e-4
        self.verbose = True
        if bare_hull is None:
            self = setup_dummy_bare_hull(self)
            self.rgp.compute_fresh_rules_graph()
        else:
            self.initialize_ship_parameters_and_values()
        
        self.initialize_bulb_parameters()
        self.coupling_constants()
        
        self.basic_bulb_rules() #cp_bulb = 1.0
        
        self.linear_parameters()
        self.nonlinear_parameters()
        
        
        self.initialize_lists()
        self = ini_coeffs(self)#.Coefficients, self.rgp)
        
        self.experiment_bulb_parameters() #add new rules here
        
        
        
    
    def get_value(self, name):
        """
            put in this level of indirection
            because the hull is technically a list of sates
            it's not so wise to [0] everywhere!
        """
        if isinstance(name, lp.Variable):
            val = self.bare_hull(name)[0]
            return val
        else:
            print 'GeometryGenerator'
            print 'name ',name,' not found'
            return None
    
    def get_value_from_ad(self, Q):
        if isinstance(Q, ad):
            return Q.value
        else:
            return Q
    
    def get_results(self):
        
        
        BulbBeam    = self.get_for_design(self.BulbBeam)
        BulbDepth   = self.get_for_design(self.BulbDepth)
        BulbLength  = self.get_for_design(self.BulbLength)
        
        
        ship_beam    = self.get_for_design(self.ship_beam)
        ship_depth   = self.get_for_design(self.ship_depth)
        ship_length  = self.get_for_design(self.ship_Lpp)
        ship_Acp     = self.get_for_design(self.ship_Acp)
        ship_Amshp   = self.get_for_design(self.ship_Amshp)
        ship_Awp     = self.get_for_design(self.ship_Awp)
        ship_Vol     = self.get_for_design(self.ship_Vol)
        
        area = self.get_value_from_ad( self.mb_ideal.Lspline.curve.area )
        boxvol = self.assumed_vol #= len*width*beam - should make tubular instead...
        alt_vol = BulbLength * area*2. #true max
        
        Cbbulb = alt_vol/boxvol #bulb block coefficient
        Cpbulb = alt_vol/alt_vol
        
        self.results_dict['BulbBeam' ] = BulbBeam
        self.results_dict['BulbDepth' ] = BulbDepth
        self.results_dict['BulbLength' ] = BulbLength
        
        # hull basics
        self.results_dict['BulbVolume' ] = self.get_value_from_ad(
                                                        alt_vol)
        self.results_dict['CBulbBlock'] = self.get_value_from_ad(
                                                        Cbbulb)
        
        self.results_dict['CBulbPrismatic'] = self.get_value_from_ad(
                                                        Cpbulb)
        
        ## cplane
        """
        #VERTICAL CENTER PLANE CURVE, CP(C):
        """
        # bbow.ideal.v_cpbc_cL_bL
        area = self.get_value_from_ad( self.cp_ideal.Lspline.curve.area )
        self.results_dict['A_lateral'] = area
        self.results_dict['CBulbCtrPln'] = self.get_value_from_ad(
                                        area/
                                        (BulbLength*BulbDepth)
                                        )
        
        #self.CBulbCtrPln = self.CBulbCtrPln == self.A_lateral/(self.BulbLength*self.BulbDepth)
        
        
        
        ## bwl  TIMES 2  (!?) XX ?
        """
        #HORIZONTAL "WATER" PLANE CURVE, BWL(C):
        """
        #bbow.ideal.v_hbc_cL_bL
        area = self.get_value_from_ad( self.bwl_ideal.Lspline.curve.area )
        area = area#*2.
        self.results_dict['A_BBwl'] = area
        self.results_dict['CBulbWtrPln'] = self.get_value_from_ad(
                                        area/
                                        (BulbBeam*BulbLength)
                                        )
        
        
        
        ## midship    TIMES 2  (!?) XX ?
        """
        # MID-BULB CURVE, MB(C):
        """
        #bbow.ideal.v_mbc_cL_bL
        area = self.get_value_from_ad( self.mb_ideal.Lspline.curve.area )
        area = area*2.
        self.results_dict['A_mid' ] = area
        self.results_dict['CBulbMid'] = self.get_value_from_ad(
                                        area/
                                        (BulbBeam*BulbDepth)
                                        )
        
        #self.CBulbMid = self.CBulbMid ==  self.A_mid/(self.BulbBeam*self.BulbDepth)
        
        
        
        
        """
        # INI-BULB CURVE:
        # interface curve:
        """
        #bbow.ideal.v_inibc_cL_tp
        #bbow.inib_ideal.Lspline
        """FWD-BULB SHAPE CURVE, FSB(C)
        """
        area = self.get_value_from_ad( self.inib_ideal.Lspline.curve.area )
        
        
        
        
        self.results_dict['Cbb'] = self.get_value_from_ad( 
                BulbBeam/ship_beam
                )
        self.results_dict['Clpr'] = self.get_value_from_ad(
                BulbLength/ship_length 
                ) 
        self.results_dict['Czb'] = self.get_value_from_ad( 
                BulbDepth/ship_depth 
                )
        
        self.results_dict['Cabt'] = self.get_value_from_ad(
                self.results_dict['A_mid' ] /ship_Amshp 
                ) 
        self.results_dict['Cabl'] = self.get_value_from_ad(
                self.results_dict['A_lateral'] /ship_Acp 
                ) 
        self.results_dict['Cvpr'] = self.get_value_from_ad(
                alt_vol/ship_Vol 
                ) 
        return



    def get_parameter(self, name):
        """Not used so far
        """
        param = self.rgp(name)
        return param
    
    
    def get_from_rgp(self, name):
        """
        -its odd that this works,
        yet dropping the treednode straight into the env
         does nothing
        yet construct returns the right answers
         supposedly without adding a new varialb to env
         
         
        cur_el = self.rgp.env.bind(self.BulbBeam.name)
        >>> self.rgp.env(cur_el)
        [ia(5.0, 8.0)]
        """
        mk_name, mk_value = self.rgp.get_name_and_value(name)
        return mk_name, mk_value
    
    def get_for_design(self, var):
        """TODO: return a thin value for each FPD parameter
        """
        mk_name,mk_value = self.get_from_rgp(var.name)
        #assert(len(mk_value)==1),"Design space is multivalued"
        x = .5 #self.sobol_seq.next()
        return  mk_value[0].getpoint(x)
    
    
    
    def get_curves(self):
        """package the various bulb curves
        with curve and rule data together.
        
        bulb coordinate system:
        
            z
            |
            |
            o -----y
           /  
          /
         x  
         
         x : longitudinal
         y : transverse
         z : vertical
        
        TODO: make transverse curves
        match bwl curve at midline!
        -probably want a bwl curve which 
        is more cylindrical and less conical...
        --working now?
        TODO: include a check.
        
        
        TODO:
            annotate each 2D ideal with it's
            3D location, highlight the missing
            coordinate
        """
        BulbBeam    = self.get_for_design(self.BulbBeam)
        BulbDepth   = self.get_for_design(self.BulbDepth)
        BulbLength  = self.get_for_design(self.BulbLength)
        
        self.assumed_vol = BulbLength*BulbBeam*BulbDepth
        
        ideal = BoundBulbCurves(BulbBeam,
                                BulbDepth,
                                BulbLength)
        self.ideal = ideal
        
        """
        # INI-BULB CURVE:
        # interface curve:
        """
        r1m,r2m = ideal.simple_projection('z','y')
        zordinate = ideal.coords2D['x']
        xb = ideal.v_inibc_cL_bL[r1m]
        xe = ideal.v_inibc_cL_tp[r1m]
        xc = .525*(xe-xb)
        self.inib_ideal = FPDIdeal(
                xb = ideal.v_inibc_cL_bL[r1m],
                xe = ideal.v_inibc_cL_tp[r1m],
                yb = ideal.v_inibc_cL_bL[r2m],
                ye = ideal.v_inibc_cL_tp[r2m],
                valphab = 0.,
                valphae = 0.,
                #Cab_given = 0.,
                #Cae_given = 0.,
                ymax = BulbBeam,
                area = self.get_for_design(self.A_mid)/2.,
                xc   = xc,
                stationloc = zordinate,
                base = ideal.v_inibc_cL_bL,
                nCV = self.nCV
                )
        self.inib_ideal.make_curve()
        self.idealdict['inib_ideal'] = self.inib_ideal
        
        
        
        """
        # MID-BULB CURVE, MB(C):
        -xb : start at the centerline in x
        -xe: run to the centerline again
        -yb: start at the baseline
        -ye: run to the top = BulbDepth
        -ymax: input the maximum width (BulbBeam) encountered along the way
        -area: A_mid, the midbulb area is the total area for this curve.
        """
        #        self.midbulbcurce_varlst = [self.BulbDepth,
        #                                    self.BulbBeam,
        #                                    self.BulbVolume,
        #                                    self.CBulbMid,
        #                                    self.CBulbBlock]
        
        r1m,r2m = ideal.simple_projection('z','y')
        zordinate = ideal.coords2D['x']
        xb = ideal.v_mbc_cL_bL[r1m]
        xe = ideal.v_mbc_cL_tp[r1m]
        xc = .525*(xe-xb)
        self.mb_ideal = FPDIdeal(
                        xb = ideal.v_mbc_cL_bL[r1m],
                        xe = ideal.v_mbc_cL_tp[r1m],
                        yb = ideal.v_mbc_cL_bL[r2m],
                        ye = ideal.v_mbc_cL_tp[r2m],
                        valphab = 0.,
                        valphae = 0.,
                        #alphab = 90.,
                        #alphae = 90.,
                        #Cab_given = 0.,
                        #Cae_given = 0.,
                        ymax = BulbBeam,
                        area = self.get_for_design(self.A_mid)/2.,
                        xc   = xc,
                        stationloc = zordinate,
                        base = ideal.v_mbc_cL_bL,
                        nCV = self.nCV
                        )
        self.mb_ideal.make_curve()
        self.idealdict['mb_ideal'] = self.mb_ideal
        
        
        
        """
        #HORIZONTAL "WATER" PLANE CURVE, BWL(C):
        -xb, the center line of the bulb
        -to xe, the BulbBeam
        -from yb, the 
        -area:  the planar area of the ~water plane~
                which is really the 
                area of elevation of widest breadth
        """
        r1w,r2w = ideal.simple_projection('x','y')
        zordinate = ideal.coords2D['z']
        self.bwl_ideal = FPDIdeal(
                xb = ideal.v_hbc_cL_bL[r1w],
                xe = ideal.v_hbc_cL_tp[r1w],
                yb = ideal.v_hbc_cL_bL[r2w],
                ye = ideal.v_hbc_cL_tp[r2w],
                ymax = BulbBeam,
                area = self.get_for_design(self.A_BBwl),
                stationloc = zordinate,
                base = ideal.v_hbc_cL_bL,
                nCV = self.nCV#,
                #xc = .3*ideal.v_hbc_cL_tp[r1w]
                )
        self.bwl_ideal.make_curve()
        self.idealdict['bwl_ideal'] = self.bwl_ideal
        
        
        
        """
        #VERTICAL CENTER PLANE CURVE, CP(C):
        -runs from xb = the interface between the bulb and the hull
        -to xe = the bulblength
        -from yb = ~0~ (by fiat) the bottom of the bulb
        -to ye =  Bulbdepth, the top of the bulb
        -with area:  the planar area of the center line
                of the bulb
        """
        r1p,r2p = ideal.simple_projection('z','x')
        zordinate = ideal.coords2D['y']
        self.cp_ideal =  FPDIdeal(
                xb = ideal.v_cpbc_cL_bL[r1p],
                xe = ideal.v_cpbc_cL_tp[r1p],
                yb = ideal.v_cpbc_cL_bL[r2p],
                ye = ideal.v_cpbc_cL_tp[r2p],
                ymax = BulbLength,
                area = self.get_for_design(self.A_lateral),
                stationloc = zordinate,
                base = ideal.v_cpbc_cL_bL,
                nCV = self.nCV
                )
        self.cp_ideal.make_curve()
        self.idealdict['cp_ideal'] = self.cp_ideal
        
        
        
        
        """FWD-BULB SHAPE CURVE, FSB(C) (not the nose!)
        """
        """TODO
            parameterize everything nicely so
            things like this are easy:
                
            -find the location of the fwd bulb shape curve.
            
            Question:  where is this curve supposed to sit, exactly?
            -I know it's 'somewhere near the bulb front end'...
            
            -how about wherever bwl curve is 
                at .75 of the total length of the bulb
        """
        sq = self.bwl_ideal.curve.FindPoint(.75*BulbLength,0)
        s = sq[0]
        #q = sq[1] #pt location of at the 
        thisymax = self.bwl_ideal.curve.CurvePoint(s)[1]
        
        """
        dev
        ----------
        
        r1m,r2m = ideal.simple_projection('z','y')
        
        zordinate = ideal.coords2D['x']
        
        xb = ideal.v_fsbc_cL_bL[r1m]
        xe = ideal.v_fsbc_cL_tp[r1m]  #full bulb depth
        yb = ideal.v_fsbc_cL_bL[r2m]
        ye = ideal.v_fsbc_cL_tp[r2m]
        
        """
        r1m,r2m = ideal.simple_projection('z','y')
        zordinate = ideal.coords2D['x']
        self.fsb_ideal = FPDIdeal(
                            xb = ideal.v_fsbc_cL_bL[r1m],
                            xe = ideal.v_fsbc_cL_tp[r1m],
                            yb = ideal.v_fsbc_cL_bL[r2m],
                            ye = ideal.v_fsbc_cL_tp[r2m],
                            #alphab = 0.,
                            #alphae = 0.,
                            valphab = 0.,
                            valphae = 0.,
                            #Cab_given = 0.,
                            #Cae_given = 0.,
                            #ymax = thisymax,
                            ymax = BulbBeam,
                            ##not used #area = np.sqrt(self.get_for_design(self.A_mid)/np.pi)
                            #area = 0.5*(np.pi*(thisymax**2)), #was using this!
                            area = self.get_for_design(self.A_mid)/2.,
                            stationloc = zordinate,
                            base = ideal.v_fsbc_cL_bL,
                            nCV = self.nCV
                            )
        self.fsb_ideal.make_curve()
        self.idealdict['fsb_ideal'] = self.fsb_ideal
        
        return
    
    
    def specialized_bow_development_nose(self):
        """
            self.cp_nose
        """
        def split_nose():
            
            zloc = self.ideal.v_fsbc_cL_bL.v[0]
            ini = self.idealdict['fsb_ideal'].Lspline.curve.CurvePoint(0.03)
            fini = self.idealdict['fsb_ideal'].Lspline.curve.CurvePoint(0.97)
            #ini = self.idealdict['fsb_ideal'].Lspline.curve.CurvePoint(0.1)
            #fini = self.idealdict['fsb_ideal'].Lspline.curve.CurvePoint(0.9)
    
            a = self.idealdict['cp_ideal'].Lspline.curve.FindNearestPoint([ini[0],zloc])
            a = 1.25*a
            
            #can't search here, the parameter will make no sense:
            #b = self.idealdict['cp_ideal'].Lspline.curve.FindNearestPoint([fini[0],zloc])
            #b=0.75*b
            
            bowc = self.idealdict['cp_ideal'].Lspline.curve.CurveSplit(a)
            #bowc = self.idealdict['cp_ideal'].Lspline.curve.CurveSplit(.25)
            b = bowc[1].FindNearestPoint([fini[0],zloc])
            b=0.75*b
            bowc = bowc[1].CurveSplit(b)[0] #final bow curve
            #bowc = bowc[1].CurveSplit(.75)[0]
            
            if self.nCV == 7:
                bowc.knot_insertion(.3333333)
                bowc.knot_insertion(.6666666)
            
            return bowc
        
        
        def make_nose(yloc = None, Wratio=.75, Aratio=.7):
            """from GeometryGenerator import FPDIdeal
            """
            ymax = self.idealdict['cp_ideal'].ymax
            if yloc is None:
                yloc = .95*ymax
            xb = self.idealdict['cp_ideal'].xb*Wratio
            xe = self.idealdict['cp_ideal'].xe*Wratio
            #            ye = self.idealdict['cp_ideal'].ye
            #            yb = self.idealdict['cp_ideal'].yb
            
            self.nose_ideal =  FPDIdeal(
                    xb = xb,
                    xe = xe,
                    yb = yloc,
                    ye = yloc,
                    ymax = ymax,
                    area = Aratio*(xb-xe)*(ymax-yloc) + \
                                yloc*(xe-xb),
                    stationloc = self.ideal.coords2D['y'],
                    base = self.ideal.v_cpbc_cL_bL, #unused DeepVector thing
                    nCV = self.nCV
                    )
            self.nose_ideal.make_curve()
            self.nose_ideal.solve_curve_ini()
            self.nose_ideal.solve_curve_fini()
            return self.nose_ideal.Lspline.curve
        
        #bowc = split_nose()
        bowc = make_nose()
        
        #more curvature at the nose of the bbow:
        Lspline = uopt.match_curve(initial_curve = bowc,
                                   return_Lspline = True,
                                   wt=1000.)
        
        r1p,r2p = self.ideal.simple_projection('z','x')
        zordinate = self.ideal.coords2D['y']
        bb = copy.deepcopy(self.ideal.v_inibc_cL_bL)
        bb.v[0] = Lspline.curve.vertices[0,1]
        bb.v[2] = Lspline.curve.vertices[0,0]
        self.cp_nose = FPDIdeal(
                xb = Lspline.curve.vertices[0,0],
                xe = Lspline.curve.vertices[-1,0],
                yb = Lspline.curve.vertices[0,1],
                ye = Lspline.curve.vertices[-1,1],
                stationloc = zordinate,
                base = bb )
        
        self.cp_nose.Lspline = Lspline
        
        self.idealdict['cp_nose'] = self.cp_nose
        
        return 
    
    
        
    #def make_nose(self, yloc = None, Wratio=.8,Aratio=.7):
    def make_nose(self, yloc = None, Wratio=.5,Aratio=.7):
        """from GeometryGenerator import FPDIdeal
        """
        ytot = self.idealdict['cp_ideal'].ymax
        xb = self.idealdict['cp_ideal'].xb
        xe = self.idealdict['cp_ideal'].xe
        xdiff = xe-xb
        ymax = self.idealdict['cp_ideal'].ymax
        if yloc is None:
            yloc = .95*ymax
        xb = xb+.4*(1.-Wratio)*xdiff
        xe = xe-.25*(Wratio)*xdiff
        #            ye = self.idealdict['cp_ideal'].ye
        #            yb = self.idealdict['cp_ideal'].yb
        
        self.nose_ideal =  FPDIdeal(
                xb = xb,
                xe = xe,
                yb = yloc,
                ye = yloc,
                ymax = ymax,
                area = Aratio*(xe-xb)*(ymax-yloc) + \
                            yloc*(xe-xb),
                stationloc = self.ideal.coords2D['y'],
                base = self.ideal.v_cpbc_cL_bL, #unused DeepVector thing
                nCV = self.nCV
                )
        self.nose_ideal.make_curve()
        self.nose_ideal.solve_curve_ini()
        self.nose_ideal.solve_curve_fini()
        return self.nose_ideal.Lspline.curve
    
    def make_bulb(self):
        
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
        bowc = self.make_nose()
        """
        tcurvenet = [self.idealdict['inib_ideal'].Lspline.curve,
                          self.idealdict['mb_ideal'].Lspline.curve,
                          self.idealdict['fsb_ideal'].Lspline.curve,
                          bowc]
        #"""
        #bvt = DeepVector(liftvertex(bowc.vertices[0],2),
        #                 B=self.idealdict['inib_ideal'].base.B)
        
        nose_curve = self.make_nose()
        
        idnet = [(self.idealdict['inib_ideal'].Lspline.curve,
                  self.idealdict['inib_ideal'].base.v[0],
                  2),
                 (self.idealdict['mb_ideal'].Lspline.curve,
                  self.idealdict['mb_ideal'].base.v[0],
                  2),
                    #                 (self.idealdict['fsb_ideal'].Lspline.curve,
                    #                  self.idealdict['fsb_ideal'].base.v[0],
                    #                  2),
                (nose_curve,
                 0.,
                 1)
                 #(self.idealdict['cp_nose'].Lspline.curve,
                 # 0.,
                 # 1)
                 ]
                 
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
            d3net.append(
                    uopt.lift_curve(curve,
                                    elevation=dv,
                                    index=index))
            
        u0 = d3net[0].vertices[0,2]
        u1 = d3net[1].vertices[0,2]
        u01 = (0.25*(u1-u0)) #taking this out to reduce # of vertices April 2018
        u12 = (0.5*(u1-u0)) #shift the copy of the 1st curve to 
        # a location in between the 2nd and third curves:
        # this assumes the second curve is just like the first one
        # and that you want a third curve that is just like it agin.
        # note that duplicate-shift will take care of the shift
        assert(u12<d3net[1].vertices[0,2]),'curves are in the wrong place'
        #print 'shift U01 by ',u01
        print 'shift U12 by ',u12
        intermediary_curve1 = uopt.duplicate_shift_curve(curve = d3net[0],
                                                        shift = u01,
                                                        index = 2)
        #note the trick!  
        # this will shift based off the original d3net[1] curve,
        #which is what was used for computation as well
        # there is thus a bit going on to pay attention to
        # in these "simple" shifts.
        intermediary_curve2 = uopt.duplicate_shift_curve(curve = d3net[0],
                                                        shift = u12,
                                                        index = 2)
        
        d3net.insert(1,intermediary_curve1) #removed April 2018
        d3net.insert(2,intermediary_curve2)
        
        
        
        self.tcurvenet = d3net
        
        # incorrect method:  longitudinals will arc out to interp both!
        #for i in range(self.tcurvenet[0].n):
        #    self.tcurvenet[-2].vertices[i,2] = self.tcurvenet[-1].vertices[i,2]
            
            
        
        
        self.lcurvenet = self.make_longitudinals()
        
        
        self.hullsurf = spline.BsplineSurface(self.tcurvenet,
                                              self.lcurvenet)
        return  
    
    
    def make_longitudinals(self):
        """tcurvenet = bbow.tcurvenet
        
        import curve as spline
        """
        tcurvenet = self.tcurvenet
        k=tcurvenet[0].k
        nump = tcurvenet[0].nump
        N = tcurvenet[0].n
        #m = self.CProfile.n
        #lvertices = np.zeros((nnew,m),float)
        tvertnet = np.asarray([el.vertices for el in tcurvenet])
        
        
        
        """
        Ship Designer acess:
            
        SD.thin_bbow_design.tcurvenet
        """
        try:
            lcurvenet = []
            for i in range(N):
                lcurvenet.append(
                    spline.interpolatedBspline(
                    tvertnet[:,i], k, nump))
        except:
            print '\n\nException thrown with tvernet in function make_longitudinals'
            print 'Data, tcurvenet[i].n:'
            print [ el.n for el in tcurvenet]
            print 'printing tvertnet...'
            for i in range(N):
                print tvertnet[:,i]
            assert(False),'ERROR: tvetnet failed.  See above'
            
           
        """Switching to more efficient longitudinal lofting
        via Lspline interpolation of the transverses
        with the interpolated curve as starting point 
        (i.e. we have a good answer)
        point being just to re-parameterize our good answer
        If we make init ADILSpline fast,
        this might beat some other re-parameterization method.
        
        see also hull_from_design_space.py
        def issue_longitudinals
        where this same issue arrises.
        
        #                nnet = []
        #                for curve in lcurvenet:
        #                    nnet.append(uopt.match_curve(initial_curve = curve,
        #                                           return_Lspline = False))
        #"""
        nnet = []
        for i,curve in zip(range(N),lcurvenet):
            pts =  tvertnet[:,i]
            # could add fixities here:
            nnet.append(uopt.interpolate_vertices(
                    initial_vertices=pts,
                    helper_curve = curve,
                    return_Lspline = False,
                    make_thb=False) )
        
        lcurvenet = nnet
        for i in range(self.tcurvenet[0].n):
            lcurvenet[i].vertices[-2,2] = lcurvenet[i].vertices[-1,2]
        return lcurvenet
    
    
    def export_curves(self, tcurvenet=None, lcurvenet=None):
        """
        """
        if tcurvenet is None: tcurvenet = self.tcurvenet
        if lcurvenet is None: lcurvenet = self.lcurvenet
        
        tlist = []
        for curve in tcurvenet:
            tlist.append([ list(pt)  for pt in curve.vertices ] )
        tknots = []
        for curve in tcurvenet:
            tknots.append(list(curve.t[1:-1]))
        
        llist = []
        for curve in lcurvenet:
            llist.append([ list(pt)  for pt in curve.vertices ] )
        lknots = []
        for curve in lcurvenet:
            lknots.append(list(curve.t[1:-1]))
            
        with open('transverse_bulb_curves','w') as f:
            pickle.dump(tlist,f)
        with open('longitudinal_bulb_curves','w') as f:
            pickle.dump(llist,f)
        
        with open('transverse_bulb_knots','w') as f:
            pickle.dump(tlist,f)
        with open('longitudinal_bulb_knots','w') as f:
            pickle.dump(llist,f)
        
        the_filename = 'transverse_bulb_curves.txt'
        with open(the_filename, 'w') as f:
            for el in tlist:
                for s in el:
                    st = '{}, {}, {}\n'.format(s[0],s[1],s[2])
                    f.write(st)
                f.write('new curve   \n')
        
        the_filename = 'transverse_bulb_knots.txt'
        with open(the_filename, 'w') as f:
            for el in tknots:
                st = ''
                for s in el:
                    st = st+' '+str(s)
                f.write(st)
                
                f.write('\n new knot vector \n')
        
        the_filename = 'longitudinal_bulb_curves.txt'
        with open(the_filename, 'w') as f:
            for el in llist:
                for s in el:
                    st = '{}, {}, {}\n'.format(s[0],s[1],s[2])
                    f.write(st)
                f.write('new curve   \n')
        
        the_filename = 'longitudinal_bulb_knots.txt'
        with open(the_filename, 'w') as f:
            for el in lknots:
                st = ''
                for s in el:
                    st = st+' '+str(s)
                f.write(st)
                
                f.write('\n new knot vector \n')
        return
    
    
    def tree_search(self,
                    sobol_seq = None):
        """States tracks it's self history
        This may cause objects never to be 
        garbage collected?
        Anyway thats how design_sq was working in opt_simple_hull
        
        
        Notes
        ----------
        http://www.ps.uni-saarland.de/~duchier/python/continuations.html
        We could search states in other ways:
            
            class Foo:
              def search(self,info,yes,no):
                if self.check(info):
                  return yes(info,no)
                else:
                  return no()
                  
              sucess could return the narrowed state
              
              failure could return the old state
              
              Then we just keep searching as before
              no logic required at the `higher level'
              -we already ask for thick variables
              and randomly thin them.
              So we'd expect fresh attempts
              until success.
        
        Honestly though, let's get back to paver style asap.
            -not to pave it all, but to retain all the info
            we actually compute about the space.
        """
        if sobol_seq is None:sobol_seq=self.sobol_seq
        #if initial_state is None: initial_state = self.rgp
        
        #design = initial_state
        
        k1 = self.Primary
        k2 = self.Coefficients
        k3 = self.Areas
        
        self.wrapped_search(k1,sobol_seq)
        self.wrapped_search(k3,sobol_seq) #sept 21
        self.wrapped_search(k2,sobol_seq)
        #self.wrapped_search(k1,sobol_seq)
        #self.wrapped_search(k3,sobol_seq)
        
        """
        redo = self.search_loop(k1,sobol_seq)
        redo = self.search_loop(k1,sobol_seq,
                            maxinner = 2)
        #"""
        """
        redo = self.search_loop([BulbLength],
                                sobol_seq)
        #"""
        #self.search_loop(k2,sobol_seq)
        #self.search_loop(k3,sobol_seq)
        return
    
    def wrapped_search(self,
                       klist,
                       sobol_seq = None,
                       maxinner = 1):
        if sobol_seq is None:sobol_seq=self.sobol_seq
        loop = True
        inner = 0
        redolist = klist
        while loop and inner<maxinner:
            redolist = self.search_loop(redolist,
                                        sobol_seq)
            if len(redolist) == 0: #now this is redundant
                loop=False
                print '*'
            else:
                print 'redo list contents: ',redolist
            inner +=1
        return
    
    
        

    def search_loop(self, 
                    inlist,
                    generator=None,
                    this_iter = 0,
                    maxinner = 1,
                    ith = -1,
                    operator='Default'):
        if generator is None: generator = self.sobol_seq
        redo = []
        for i in range(len(inlist)):
            var = inlist[i]
            #x = generator.next() 
            x = random.random()
            #x = .5
            #checker = self.set_this(var,x)
            if operator == 'Default':
                checker = self.set_this(var,x)
            elif operator == 'Narrow':
                checker = self.narrow_this(var)
                
            if not checker:
                redo.append(var)
        if len(redo)==0:
            return redo
        
        else:
            if this_iter>maxinner:
                return redo
            else:
                this_iter += 1
                print 'redo',redo
                if this_iter %2 == 1:
                    return self.search_loop(redo,
                                            generator,
                                            this_iter,
                                            operator='Narrow')
                else:
                    return self.search_loop(redo,
                                            generator,
                                            this_iter,
                                            operator='Default')
    
    

    
    def set_this(self, var, x, ith=-1):
        """
            var = key
            ith=-1
            
            ?Do this the old fashioned way 
            to go back one state in the event
            that an equality rule returns the null state
            
            or, -since that would require adding 
                  def equal()
                  to the code,
                -Which would then mean future users
                  building similar classes would have to
                  be smart about internals...
            
            ?how about copy.copy(self.rgp.env)
            try the equality
            accept new env only if it works.
            
            (only) issue: sometimes backing up
            further is nice
            
        """
        mk_name, mk_val = self.get_from_rgp(var) #
        if mk_val[ith] is None:
            return True
        val = mk_val[ith].getpoint(x)
        #------ performance boost ;)
        space = self.tol*.3
        value = ia(val-space,val+space)
        #------
        # set new value:
        var = var == value
        #-----------------------------------------
        ret_state = copy.copy(self.rgp.env)
        self.rgp.add_one_rule(var,var.name)
        self.rgp.compute_fresh_rules_graph()#maxiter=8)
        #-----------------------------------------
        #
        #now check if ok or not
        #(checking here instead of 
          #RulesGraphProcessor's AC_revise function )
        if len(self.rgp.env.states)==0:
            self.rgp.reset(ret_state,var)
            #v1i = mk_val[ith].inf
            #v2s = mk_val[ith].sup
            #v1 = ia(v1i,val)
            #v2 = ia(val,v2s)
            print ' sending ',mk_name, 'to the redo list'
            #print 'v1 = ',v1
            #print 'v2 = ',v2
            #self.rgp.env = self.rgp.env.split( (mk_val[ith],v1,v2) )
            return False #, ret_state
        else:
            self.rgp.env.parent = ret_state #sept 30 
            ret_state.children.append(self.rgp.env) #sept 30 
            if self.rgp._verbose: print 'done ', mk_name,'=>', value
            self.rgp.varsmap[var.name].flist.pop()
            return True #,  self.rgp.env
    
    

    
    def narrow_this(self, var, ith=-1):
        """
            var = key
            ith=-1
            
        """
        x1=.25
        x2=.75
        mk_name, mk_val = self.get_from_rgp(var) #
        if mk_val[ith] is None:
            return True
        vali = mk_val[ith].getpoint(x1)
        vals = mk_val[ith].getpoint(x2)
        value = ia(vali,vals)
        var = var == value
        #-----------------------------------------
        ret_state = copy.copy(self.rgp.env)
        self.rgp.add_one_rule(var,var.name)
        self.rgp.compute_fresh_rules_graph()#maxiter=8)
        #-----------------------------------------
        if len(self.rgp.env.states)==0:
            #----------------------------
            self.rgp.reset(ret_state,var)
            #----------------------------
            v1i = mk_val[ith].inf
            v2s = mk_val[ith].sup
            v1 = ia(v1i,vali)
            v2 = ia(vals,v2s)
            if self.verbose:
                #print 'not splitting ',mk_name
                print ' splitting wide! ',mk_name
                print 'v1 = ',v1
                print 'v2 = ',v2
            self.rgp.env = self.rgp.env.split( (mk_val[ith],v1,v2) )
            return False #, ret_state
        else:
            if self.rgp._verbose: print 'narrowed ', mk_name,'=>', value
            self.rgp.varsmap[var.name].flist.pop()
            return True #,  self.rgp.env
    
    
    def get_thick_list(self):
        todo_list = []
        for var in self.rgp.vars:
            try:
                attr = self.__getattribute__(var.name)
                mk_name, mk_val = self.get_from_rgp(attr)
                #if multivalued, pick first
                iaval = mk_val[0]
                diff = iaval.sup - iaval.inf
                if np.linalg.norm(diff)<self.tol:
                    pass
                else:
                    todo_list.append(attr)
            except:
                pass
        return todo_list
    
    
    def get_node_list(self, option=''):
        if option.lower() == 'thick':
            return self.get_thick_list()
        else:
            todo_list = []
            for var in self.rgp.vars:
                try:
                    attr = self.__getattribute__(var.name)
                    mk_name, mk_val = self.get_from_rgp(attr)
                    # if multivalued, simply pick the first.
                    iaval = mk_val[0]
                    diff = iaval.sup - iaval.inf
                    todo_list.append(attr)
                except:
                    pass
            return todo_list
    
    
    def coupling_constants(self, rgp =None):
        if rgp is None:
            rgp = self.rgp
        #
        #linear
        #-----------------------------------------
        self.Cbb    = lp.PStates(name='Cbb')
        self.Clpr   = lp.PStates(name='Clpr')
        self.Czb    = lp.PStates(name='Czb')
        #
        #nonlinear
        #-----------------------------------------
        self.Cabt   = lp.PStates(name='Cabt')
        self.Cabl   = lp.PStates(name='Cabl')
        self.Cvpr   = lp.PStates(name='Cvpr')
        #
        return rgp
    
    
    def initialize_ship_parameters_and_values(self):#, rgp=None):
        #if rgp is None:
        #    rgp = self.rgp
        #-----------------------------------------
        # quantities of m**1
        self.ship_beam = lp.PStates(name='ship_beam') 
        # note: could use bare_hull var names instead. 
        #       e.g. lp.PStates(name=self.init_ship_beam.name)
        self.ship_depth = lp.PStates(name='ship_depth')
        self.ship_Lpp = lp.PStates(name='ship_Lpp')
        #
        #quantities of m**2
        self.ship_Amshp = lp.PStates(name='ship_Amshp')
        self.ship_Acp = lp.PStates(name='ship_Acp')
        self.ship_Awp = lp.PStates(name='ship_Awp')
        #
        #quantities of m**3
        self.ship_Vol = lp.PStates(name='ship_Vol')
        #-----------------------------------------
        #
        #existing ship:
        #-----------------------------------------
        # set the ship values in the bulb environement:
        self.ship_beam   = self.ship_beam == self.get_value(
                                    self.init_ship_beam)
        self.ship_depth  = self.ship_depth == self.get_value(
                                    self.init_ship_depth)
        self.ship_Lpp    = self.ship_Lpp == self.get_value(
                                    self.init_ship_Lpp)
        self.ship_Amshp  = self.ship_Amshp == self.get_value(
                                    self.init_ship_Amshp)
        self.ship_Acp    = self.ship_Acp == self.get_value(
                                    self.init_ship_Acp)
        self.ship_Awp    = self.ship_Awp == self.get_value(
                                    self.init_ship_Awp)
        self.ship_Vol    = self.ship_Vol == self.get_value(
                                    self.init_ship_Vol)
        #-----------------------------------------
        
        #-----------------------------------------
        self.rgp.add_one_rule(self.ship_beam,self.ship_beam.name)
        self.rgp.add_one_rule(self.ship_depth,self.ship_depth.name)
        self.rgp.add_one_rule(self.ship_Lpp,self.ship_Lpp.name)
        self.rgp.add_one_rule(self.ship_Amshp,self.ship_Amshp.name)
        self.rgp.add_one_rule(self.ship_Acp,self.ship_Acp.name)
        self.rgp.add_one_rule(self.ship_Awp,self.ship_Awp.name)
        self.rgp.add_one_rule(self.ship_Vol,self.ship_Vol.name)
        #-----------------------------------------
        self.rgp.compute_fresh_rules_graph()
        return #rgp
    
    
    
    def initialize_bulb_parameters(self):#, rgp=None):
        #if rgp is None:
        #    rgp = self.rgp
        
        #bulb areas
        self.A_mid      = lp.PStates(name='A_mid')
        self.A_lateral  = lp.PStates(name='A_lateral')
        #self.A_flat     = lp.PStates(name='A_flat') #A_BBwl instead
        self.A_BBwl     = lp.PStates(name='A_BBwl')
        
        
        #use the principle of the minimum square enclosing box (cartesian b/c ia is cartesian)
        self.BulbBeam    = lp.PStates(name = 'BulbBeam')     #Bulb half beam
        self.BulbDepth   = lp.PStates(name = 'BulbDepth')    #Bulb depth
        self.BulbLength  = lp.PStates(name = 'BulbLength')   #Bulb max length (min square enclosing box length)
        
        self.BulbVolume  = lp.PStates(name = 'BulbVolume') 
        
        self.CBulbMid    = lp.PStates(name = 'CBulbMid') #Bulb midsection coefficient
        self.CBulbCtrPln = lp.PStates(name = 'CBulbCtrPln') #Bulb centerplane profile area coefficient
        self.CBulbWtrPln = lp.PStates(name = 'CBulbWtrPln') #Bulb waterplane area coefficient
        
        self.CBulbBlock     = lp.PStates(name = 'CBulbBlock') #Bulb block coefficient 
        self.CBulbPrismatic = lp.PStates(name = 'CBulbPrismatic') #Bulb prismatic coefficient 
        
        return #rgp
    
    
    def initialize_lists(self):
        
        self.Primary = [#self.CBulbBlock,
                        self.BulbLength,
                        self.BulbDepth,
                        self.BulbBeam,
                        self.BulbVolume]
        
        self.Areas = [self.A_mid,
                      self.A_lateral,
                      self.A_BBwl]
        
        
        
        self.Coefficients = [self.Cbb,
                             self.Clpr,
                             self.Czb,
                             self.Cabt,
                             self.Cabl,
                             self.Cvpr,
                             self.CBulbMid,
                             self.CBulbCtrPln,
                             self.CBulbWtrPln,
                             self.CBulbBlock,
                             self.CBulbPrismatic]
        return 
    
    def print_bulb_state(self):
        self.print_list(self.Primary)
        self.print_list(self.Coefficients)
        self.print_list(self.Areas)
        return
        
    def print_list(self, list_):
        for key in list_:
            print self.get_from_rgp(key)
        return
    
    def print_thick_list(self):
        hl = self.get_thick_list()
        self.print_list(hl)
        return
    
    
    
    
    def basic_bulb_rules(self):#, rgp=None):
        #if rgp is None:
        #    rgp = self.rgp
        """-----------------------------------------------
        RULES DERIVED FROM THE BARE HULL SOLVER (not bare hull CSP/CLP)
        """
        """-----------------------------------------------
        1.)If A_mid comes from the bare hull, 
        incorporate it first
        """
        if self.given_A_mid is not None:
            frac = self.tol*2.
            self.A_mid = self.A_mid == ia(self.given_A_mid-frac, 
                                          self.given_A_mid+frac)
            self.rgp.add_one_rule(self.A_mid, 'A_mid')
            self.rgp.compute_fresh_rules_graph()
        
        """-----------------------------------------------
        1.)If BulbLength comes from the bare hull, 
        incorporate it first
        """
        if self.given_length is not None:
            frac = self.tol*2.
            self.BulbLength = self.BulbLength == ia(self.given_length-frac, 
                                                    self.given_length+frac)
            self.rgp.add_one_rule(self.BulbLength, 'BulbLength')
            self.rgp.compute_fresh_rules_graph()
            
        
        """-----------------------------------------------
        # cross section parameter:
        Rule:  Midbulb_Area < max_Beam * max_Depth
        CBulbMid -> [0.,1.]
        CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
        """
        self.CBulbMid = self.CBulbMid ==  self.A_mid/(self.BulbBeam*self.BulbDepth)
        self.rgp.add_one_rule(self.CBulbMid, self.CBulbMid.name  )
        self.rgp.compute_fresh_rules_graph()
        #
        """-----------------------------------------------
        Rule: z-y_area < max_length * max_Depth
        """
        self.CBulbCtrPln = self.CBulbCtrPln == self.A_lateral/(self.BulbLength*self.BulbDepth)
        self.rgp.add_one_rule(self.CBulbCtrPln, 'CBulbCtrPln')
        self.rgp.compute_fresh_rules_graph()
        #
        """-----------------------------------------------
        Rule: wL_area < max_length * max half-Beam
        """
        self.CBulbWtrPln = self.CBulbWtrPln == self.A_BBwl/(self.BulbLength*self.BulbBeam)
        self.rgp.add_one_rule(self.CBulbWtrPln, 'CBulbWtrPln')
        self.rgp.compute_fresh_rules_graph()
        #
        """-----------------------------------------------
        Rule: Bulb_vol < max_length * max_Depth * max *BulbBeam
        """
        self.CBulbBlock = self.CBulbBlock == self.BulbVolume/(self.BulbLength*self.BulbDepth*self.BulbBeam)
        self.rgp.add_one_rule(self.CBulbBlock, 'CBulbBlock')
        self.rgp.compute_fresh_rules_graph()
        mkvar,mkval = self.get_from_rgp('CBulbBlock')
        print 'CBulbBlock value = ',mkval
        #
        
        """-----------------------------------------------
        Rule: CBulbPrismatic ~= 1.0 for a tube!
        """
#        self.CBulbPrismatic = self.CBulbPrismatic ==  ia(.99,1.00001)
#        self.rgp.add_one_rule(self.CBulbPrismatic, self.CBulbPrismatic.name)
#        self.rgp.compute_fresh_rules_graph()
        
        self.CBulbPrismatic = self.CBulbPrismatic == self.BulbVolume/(self.BulbLength*self.A_mid)
        self.rgp.add_one_rule(self.CBulbPrismatic, 'CBulbPrismatic')
        self.rgp.compute_fresh_rules_graph()
        #
        self.BulbVolume = self.BulbVolume == self.BulbLength*self.A_mid
        self.rgp.add_one_rule(self.BulbVolume, self.BulbVolume.name)
        self.rgp.compute_fresh_rules_graph()
        """-----------------------------------------------
        Rule: Bulb_vol < max_length * mid_bulb_area
        """
        #
        #self.rgp = rgp
        
        return #rgp
    
    
    def linear_parameters(self):#, rgp=None):
        """Kracht Design of Bulbous Bows, 1978
        """
        #if rgp is None:
        #    rgp = self.rgp
        
        #Kracht linear relations:
        
        # breadth parameter:
        self.Cbb     = self.Cbb  ==  self.BulbBeam/self.ship_beam
        self.rgp.add_one_rule( self.Cbb)
        self.rgp.compute_fresh_rules_graph()
        
        # length parameter:
        self.Clpr    = self.Clpr ==  self.BulbLength/self.ship_Lpp
        self.rgp.add_one_rule(self.Clpr)
        self.rgp.compute_fresh_rules_graph()
        
        # depth parameter:
        #self.Czb     = self.Czb  ==  self.BulbLength/self.ship_depth
        self.Czb     = self.Czb  ==  self.BulbDepth/self.ship_depth
        self.rgp.add_one_rule( self.Czb)
        self.rgp.compute_fresh_rules_graph()
        #
        #
        #
        #self.rgp = rgp
        return #rgp
    
    
    
    
    def nonlinear_parameters(self):#, rgp=None):
        """Kracht Design of Bulbous Bows, 1978
        """
        #Kracht nonlinear relations:
        
        # cross section parameter:
        self.Cabt    = self.Cabt ==  self.A_mid/self.ship_Amshp
        
        # lateral section parameter:
        self.Cabl    = self.Cabl ==  self.A_lateral/self.ship_Acp #departure from Kracht
        
        # volumetric parameter:
        self.Cvpr    = self.Cvpr ==  self.BulbVolume/self.ship_Vol
        #
        self.rgp.add_one_rule( self.Cabt)
        self.rgp.compute_fresh_rules_graph()
        #
        self.rgp.add_one_rule(self.Cabl)
        self.rgp.compute_fresh_rules_graph()
        #
        self.rgp.add_one_rule( self.Cvpr)
        self.rgp.compute_fresh_rules_graph()
        #
        #self.rgp = rgp
        return #rgp
    
    
    
    def lowhigh(self, lowc,highc,this):
        """automated sensible setter for the bulb geometry
        starting points
        
        usage:
            
        self.BulbLength     = self.BulbLength ==  self.lowhigh(plow,phigh, 
                                                               this = self.BulbLength)'''
        """
        plow = lowc*self.get_from_rgp(this)[1][0].getpoint(.5)
        phigh = highc*self.get_from_rgp(this)[1][0].getpoint(.5)
        #use'''self.BulbLength     = self.BulbLength == ia(plow,phigh)'''
        return ia(plow,phigh)
    def fraction(self, lowc,highc,this):
        """automated sensible setter for the bulb geometry
        starting points
        
        usage:
            
        self.BulbLength     = self.BulbLength ==  self.lowhigh(plow,phigh, 
                                                               this = self.BulbLength)'''
        """
        plow = lowc*self.get_from_rgp(this)[1][0].getpoint(1.)
        phigh = highc*self.get_from_rgp(this)[1][0].getpoint(1.)
        #use'''self.BulbLength     = self.BulbLength == ia(plow,phigh)'''
        return ia(plow,phigh)
    
    def experiment_bulb_parameters(self):#, rgp=None):
        #if rgp is None:
        #    rgp = self.rgp
        
        print 'testing experiment_bulb_parameters'
        
        #self.A_mid     = self.A_mid     == self.lowhigh(.05**2,0.1**2,self.ship_Amshp)
        #self.BulbBeam  = self.BulbBeam  == self.lowhigh(.05,.4,self.ship_beam)
        
        #self.A_mid = self.A_mid == ia(10.,15.)
        #self.A_mid = self.A_mid == ia(5.,15.)
        
        # smaller bulbs make for easier fitting:
        #self.A_mid = self.A_mid == ia(2.,15.) #removed May 22 2018
        
        #self.BulbBeam = self.BulbBeam == ia(1.,2.5)
        #self.BulbBeam  = self.BulbBeam  == self.lowhigh(.1,.14,self.ship_beam)  #removed May 22 2018
        # and replaced with this:
        
        if self.given_A_mid is None:
            self.Cbb = self.Cbb == ia(.15,.24) #beam fraction
            self.rgp.add_one_rule( self.Cbb)
            self.rgp.compute_fresh_rules_graph()
        else:
            self.Cbb = self.Cbb == ia(.1,.5) #beam fraction
            self.rgp.add_one_rule( self.Cbb)
            self.rgp.compute_fresh_rules_graph()
        
        if self.given_A_mid is None:
            self.CBulbMid = self.CBulbMid == ia(.35,1.) #Cx area to area availible
            self.rgp.add_one_rule(self.CBulbMid,'CBulbMid')
        else:
            self.CBulbMid = self.CBulbMid == ia(.1,1.) #Cx area to area availible
            #BulbDepth = BulbDepth == ia(-10.1,10.)
            #self.rgp.add_one_rule(self.A_mid,'A_mid')
            #self.rgp.add_one_rule(self.BulbBeam,'BulbBeam')
            self.rgp.add_one_rule(self.CBulbMid,'CBulbMid')
        self.rgp.compute_fresh_rules_graph()
        
        
        
        if self.given_A_mid is None:
            self.CBulbWtrPln = self.CBulbWtrPln == ia(.7,.97)
            #
            self.rgp.add_one_rule(self.CBulbWtrPln,'CBulbWtrPln')
            
            self.rgp.compute_fresh_rules_graph()
        else:
            self.CBulbWtrPln = self.CBulbWtrPln == ia(.5,.97)
            #
            self.rgp.add_one_rule(self.CBulbWtrPln,'CBulbWtrPln')
            
            self.rgp.compute_fresh_rules_graph()
            
        
        
        
        
        #self.BulbLength     = self.BulbLength == ia(10.,15.)
        #self.BulbLength     = self.BulbLength == self.lowhigh(0.03,0.1,self.ship_Lpp)
        #self.Clpr = self.Clpr == ia(.03,.1)
        
        if self.given_length is None:
            self.Clpr = self.Clpr == ia(.03,.045) #length fraction
            self.rgp.add_one_rule( self.Clpr)
            self.rgp.compute_fresh_rules_graph()
        else:
            self.Clpr = self.Clpr == ia(.03,.5) #length fraction
            self.rgp.add_one_rule( self.Clpr)
            self.rgp.compute_fresh_rules_graph()
        
        
        
        if self.given_A_mid is None:
            #
            #self.BulbDepth      = self.BulbDepth == ia(2.5,10.)
            #self.BulbDepth  = self.BulbDepth  == self.lowhigh(.3,.4,self.ship_depth) #may 22 2018
            #replaced with this:
            self.Czb = self.Czb == ia(.28,.35)  #depth fraction
            self.rgp.add_one_rule( self.Czb)
            self.rgp.compute_fresh_rules_graph()
        else:
            self.Czb = self.Czb == ia(.1,.45)  #depth fraction
            self.rgp.add_one_rule( self.Czb)
            self.rgp.compute_fresh_rules_graph()
            
        
        #
        
        if self.given_A_mid is None:
            self.CBulbCtrPln    = self.CBulbCtrPln == ia(.5,1.) #Cx ctr plane fraction
            #now add them
            #self.rgp.add_one_rule(self.BulbLength,'BulbLength')#may 22
            #self.rgp.add_one_rule(self.BulbDepth,'BulbDepth') #may 22
            self.rgp.add_one_rule(self.CBulbCtrPln)
            
            self.rgp.compute_fresh_rules_graph()
        else:
            self.CBulbCtrPln    = self.CBulbCtrPln == ia(.5,1.) #Cx ctr plane fraction
            #now add them
            #self.rgp.add_one_rule(self.BulbLength,'BulbLength')#may 22
            #self.rgp.add_one_rule(self.BulbDepth,'BulbDepth') #may 22
            self.rgp.add_one_rule(self.CBulbCtrPln)
            
            self.rgp.compute_fresh_rules_graph()
            
        
        #
        ##
        ##************************************* new bulb blending rules
        ##
        
        #keep the bulb short:
        #self.Clpr    = self.Clpr ==  ia(.01,.045)
        #self.rgp.add_one_rule(self.Clpr, 'Clpr')
        #self.rgp.compute_fresh_rules_graph()
        
        #
        ##
        ##************************************* end rules
        ##
        
        self.rgp.compute_fresh_rules_graph()
        
        return #rgp
    

#def generate_a_thin_bulb_design():
#    return
    
def setup_dummy_bare_hull(bbobj):#, rgp=None):
    #if rgp is None:
    #    rgp = bbobj.rgp
    #-----------------------------------------
    # quantities of m**1
    bbobj.ship_beam = lp.PStates(name='ship_beam') 
    # note: could use bare_hull var names instead. 
    #       e.g. lp.PStates(name=self.init_ship_beam.name)
    bbobj.ship_depth = lp.PStates(name='ship_depth')
    bbobj.ship_Lpp = lp.PStates(name='ship_Lpp')
    #
    #quantities of m**2
    bbobj.ship_Amshp = lp.PStates(name='ship_Amshp')
    bbobj.ship_Acp = lp.PStates(name='ship_Acp')
    bbobj.ship_Awp = lp.PStates(name='ship_Awp')
    #
    #quantities of m**3
    bbobj.ship_Vol = lp.PStates(name='ship_Vol')
    #-----------------------------------------
    #
    #existing ship:
    #-----------------------------------------
    # quantities of m**1
    #    ship_beam   = bbobj.ship_beam 
    #    ship_depth  = bbobj.ship_depth
    #    ship_Lpp    = bbobj.ship_Lpp 
    #    #
    #    #quantities of m**2
    #    ship_Amshp  = bbobj.ship_Amshp 
    #    ship_Acp    = bbobj.ship_Acp 
    #    ship_Awp    = bbobj.ship_Awp 
    #    #
    #    #quantities of m**3
    #    ship_Vol    = bbobj.ship_Vol 
    #-----------------------------------------
    #
    #-----------------------------------------
    # set the ship values in the bulb environement:
    bbobj.ship_beam   = bbobj.ship_beam == ia(17.4663142374, 17.4663142374)
    
    bbobj.ship_depth  = bbobj.ship_depth == ia(16.2051841085, 16.2051841085)
    
    bbobj.ship_Lpp    = bbobj.ship_Lpp == ia(111.099919763, 111.099919763)
    
    bbobj.ship_Amshp  = bbobj.ship_Amshp == ia(261.639572047, 261.639572047)
    
    bbobj.ship_Acp    = bbobj.ship_Acp == ia(1656.36308186, 1656.36308186)
    
    bbobj.ship_Awp    = bbobj.ship_Awp == ia(1736.75296874, 1736.75296874)
    
    bbobj.ship_Vol    = bbobj.ship_Vol == ia(27043.7825521, 27043.7825521)
    #-----------------------------------------
    
    #-----------------------------------------
    bbobj.rgp.add_one_rule(bbobj.ship_beam,'ship_beam')
    bbobj.rgp.add_one_rule(bbobj.ship_depth,'ship_depth')
    bbobj.rgp.add_one_rule(bbobj.ship_Lpp,'ship_Lpp')
    bbobj.rgp.add_one_rule(bbobj.ship_Amshp,'ship_Amshp')
    bbobj.rgp.add_one_rule(bbobj.ship_Acp,'ship_Acp')
    bbobj.rgp.add_one_rule(bbobj.ship_Awp,'ship_Awp')
    bbobj.rgp.add_one_rule(bbobj.ship_Vol,'ship_Vol')
    #-----------------------------------------
    #rgp.compute_fresh_rules_graph()
    return bbobj#, rgp



    
def generate_quasirandom_sequence(self,
                                  allow_range):   
    """Not used
    """     
    s1 = sobol.sobolSeq([1,1],[1,1])
    inf = allow_range[0]
    sup = allow_range[1]
    r = sup-inf
    x = np.asarray(
            [ next(s1) for _ in range(self.N) ]
        )
    return r*x+inf



def old_test():
    print 'Kraft Bulbous Bow Parameters'
    print'\n Linear Parameters:'
    
 
    
    A_mid = lp.PStates(name='A_mid')

    
    #use the principle of the minimum square enclosing box (cartesian b/c ia is cartesian)
    BulbBeam    = lp.PStates(name = 'BulbBeam')     #Bulb half beam
    BulbDepth   = lp.PStates(name = 'BulbDepth')    #Bulb depth

    
    CBulbMid = lp.PStates(name = 'CBulbMid') #Bulb midsection coefficient


    
    
    #TODO: fix this with class to construct rules graph!
    
    A_mid = A_mid == ia(10.,20.) #issue: order of assignment "matters"
    BulbBeam = BulbBeam == ia(5.,10.)
    CBulbMid = CBulbMid == ia(.5,1.)
    
    
    rgp = RulesGraphProcessor()
    
    rgp.add_one_rule(A_mid,A_mid.name)
    rgp.add_one_rule(BulbBeam,BulbBeam.name)
    rgp.add_one_rule(CBulbMid,CBulbMid.name)
    
    
    rgp.compute_fresh_rules_graph()
    
    
    """-----------------------------------------------
    Rule:  Midbulb_Area < max_Beam * max_Depth
    CBulbMid -> [0.,1.]
    CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
    """
    CBulbMid = CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
    
    #rgp.add_one_rule(BulbDepth,BulbDepth.name)
    rgp.add_one_rule(CBulbMid,CBulbMid.name)
    
    rgp.compute_fresh_rules_graph()
    #"""
    
    #rgp.compute_fresh_rules_graph()
    
    #BulbDepth = BulbDepth == ia(5.,10.)
    #rgp.add_one_rule(BulbDepth,BulbDepth.name)
    #rgp.compute_fresh_rules_graph()
    
    
    
    
    
    print '\n\n state after:'
    #print st
    
    print rgp.env
    
    
    
    ellist = [BulbDepth,
              CBulbMid,
              BulbBeam,
              A_mid]
    
    ith = -1
    x=.5
    redo = []
    for i in range(len(ellist)):
        var = ellist[i]
        mk_name, mk_val = rgp.get_name_and_value(var)
        print mk_name, mk_val
        val = mk_val[ith].getpoint(x)
        value = ia(val,val)
        var = var == value
        #-----------------------------------------
        ret_state = copy.copy(rgp.env)
        rgp.add_one_rule(var,var.name)
        rgp.compute_fresh_rules_graph()
        #-----------------------------------------
        
        if len(rgp.env.states)==0:
            #rgp.env = ret_state
            rgp.reset(ret_state,var)
            redo.append(var)
        else:
            if rgp._verbose: print 'done ', mk_name,'=>', value
        
    
    return



def small_test():
    print 'Kraft Bulbous Bow Parameters'
    print'\n Linear Parameters:'
    
    
    
    
    ##
    ##************************************* bulb
    ##
    A_mid = lp.PStates(name='A_mid')
    A_lateral = lp.PStates(name='A_lateral')
    A_BBwl = lp.PStates(name='A_BBwl')
  
    
    #use the principle of the minimum square enclosing box (cartesian b/c ia is cartesian)
    BulbBeam    = lp.PStates(name = 'BulbBeam')     #Bulb half beam
    BulbDepth   = lp.PStates(name = 'BulbDepth')    #Bulb depth
    BulbLength  = lp.PStates(name = 'BulbLength')   #Bulb max length (min square enclosing box length)
    
    
    BulbVolume  = lp.PStates(name = 'BulbVolume') 
        
    
    #    CBulbMid = lp.PStates(name = 'CBulbMid') #Bulb midsection coefficient
    #
    #    CBulbBlock     = lp.PStates(name = 'CBulbBlock')
    #    CBulbPrismatic = lp.PStates(name = 'CBulbPrismatic') 
    #    
    #    
    #    CBulbCtrPln = lp.PStates(name = 'CBulbCtrPln') #Bulb centerplane profile area coefficient
    #    CBulbWtrPln = lp.PStates(name = 'CBulbWtrPln') #Bulb waterplane area coefficient
    
        
    CBulbMid = lp.PStates(name = 'CBulbMid') #Bulb midsection coefficient
    CBulbCtrPln = lp.PStates(name = 'CBulbCtrPln') #Bulb centerplane profile area coefficient
    CBulbWtrPln = lp.PStates(name = 'CBulbWtrPln') #Bulb waterplane area coefficient
    
    CBulbBlock     = lp.PStates(name = 'CBulbBlock')
    CBulbPrismatic = lp.PStates(name = 'CBulbPrismatic') 
    
    
    
    #TODO: fix this with class to construct rules graph!
    
    A_mid = A_mid == ia(10.,20.) #issue: order of assignment "matters"
    BulbBeam = BulbBeam == ia(5.,10.)
    #BulbDepth = BulbDepth == ia(-10.1,10.)
    
    CBulbMid = CBulbMid == ia(.1,1.)
    CBulbWtrPln = CBulbWtrPln == ia(.1,1.)
    
    
    rgp = RulesGraphProcessor()
    #"""
    
    #rgp.add_one_rule(BulbDepth,'BulbDepth')
    rgp.add_one_rule(A_mid,'A_mid')
    rgp.add_one_rule(BulbBeam,'BulbBeam')
    rgp.add_one_rule(CBulbMid,'CBulbMid')
    rgp.add_one_rule(CBulbWtrPln,'CBulbWtrPln')
    #rgp.add_one_rule(BulbDepth,'BulbDepth') #should not be needed!
    
    rgp.compute_fresh_rules_graph()
    
    #rgp.AC_revise()
    #rgp.env
    
    ##
    ##*************************CBulbCtrPln************ end bulb
    ##
    
    ##
    ##************************************* ship
    ##
    ship_beam = lp.PStates(name='ship_beam') 
    # note: could use bare_hull var names instead. 
    #       e.g. lp.PStates(name=self.init_ship_beam.name)
    ship_depth = lp.PStates(name='ship_depth')
    ship_Lpp = lp.PStates(name='ship_Lpp')
    #
    #quantities of m**2
    ship_Amshp = lp.PStates(name='ship_Amshp')
    ship_Acp = lp.PStates(name='ship_Acp')
    ship_Awp = lp.PStates(name='ship_Awp')
    #
    #quantities of m**3
    ship_Vol = lp.PStates(name='ship_Vol')
    
    
    #-----------------------------------------
    # set the ship values in the bulb environement:
    ship_beam   = ship_beam == ia(17.4663142374, 17.4663142374)
    
    ship_depth  = ship_depth == ia(16.2051841085, 16.2051841085)
    
    ship_Lpp    = ship_Lpp == ia(111.099919763, 111.099919763)
    
    ship_Amshp  = ship_Amshp == ia(261.639572047, 261.639572047)
    
    ship_Acp    = ship_Acp == ia(1656.36308186, 1656.36308186)
    
    ship_Awp    = ship_Awp == ia(1736.75296874, 1736.75296874)
    
    ship_Vol    = ship_Vol == ia(27043.7825521, 27043.7825521)
    #-----------------------------------------
    
    
    
    #-----------------------------------------
    rgp.add_one_rule(ship_beam,'ship_beam')
    rgp.add_one_rule(ship_depth,'ship_depth')
    rgp.add_one_rule(ship_Lpp,'ship_Lpp')
    rgp.add_one_rule(ship_Amshp,'ship_Amshp')
    rgp.add_one_rule(ship_Acp,'ship_Acp')
    rgp.add_one_rule(ship_Awp,'ship_Awp')
    rgp.add_one_rule(ship_Vol,'ship_Vol')
    #-----------------------------------------
    ##
    ##************************************* end ship
    ##
    rgp.compute_fresh_rules_graph()
    
    
    
    ##
    ##************************************* bulb rules
    ##
    """-----------------------------------------------
    Rule:  Midbulb_Area < max_Beam * max_Depth
    CBulbMid -> [0.,1.]
    CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
    """
    CBulbMid = CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
    
    """-----------------------------------------------
    Rule: z-y_area < max_length * max_Depth
    """
    CBulbCtrPln = CBulbCtrPln == A_lateral/(BulbLength*BulbDepth)
    
    """-----------------------------------------------
    Rule: wL_area < max_length * max_Depth
    """
    #CBulbWtrPln = CBulbWtrPln == A_BBwl/(BulbLength*BulbDepth)
    CBulbWtrPln = CBulbWtrPln == A_BBwl/(BulbLength*BulbBeam)
    
    
    
    #2 more rules!    
    """-----------------------------------------------
    Rule: Bulb_vol < max_length * max_Depth * max *BulbBeam
    """
    CBulbBlock = CBulbBlock == BulbVolume/(BulbLength*BulbDepth*BulbBeam)
    #
    """-----------------------------------------------
    Rule: Bulb_vol < max_length * mid_bulb_area
    """
    CBulbPrismatic = CBulbPrismatic == BulbVolume/(BulbLength*A_mid)    
    #


    
    rgp.add_one_rule(CBulbMid,'CBulbMid')
    rgp.add_one_rule(CBulbCtrPln,'CBulbCtrPln')
    rgp.add_one_rule(CBulbWtrPln,'CBulbWtrPln')
    
    
    #
    # Check for interval splitting:
    #
    #CBulbMid = CBulbMid == ia(-1.,1.)
    #rgp.add_one_rule(CBulbMid,'CBulbMid')
    
    rgp.add_one_rule(CBulbBlock, 'CBulbBlock')
    rgp.add_one_rule(CBulbPrismatic, 'CBulbPrismatic')
    
    rgp.compute_fresh_rules_graph()
    #"""
    
    #BulbLength  = BulbLength == ia(10.,15.)
    BulbLength  = BulbLength == ia(1.,15.)
    #CBulbCtrPln = CBulbCtrPln == ia(.5,1.)
    CBulbCtrPln = CBulbCtrPln == ia(0.,1.)
    rgp.add_one_rule(BulbLength,'BulbLength')
    rgp.add_one_rule(CBulbCtrPln,'CBulbCtrPln')
    
    
    #still good--------------------------------------------------
    
    rgp.compute_fresh_rules_graph()
    
    #BulbDepth = BulbDepth == ia(1.,10.)
    BulbDepth = BulbDepth == ia(5.,10.)
    rgp.add_one_rule(BulbDepth,'BulbDepth')
    
    
    #single space left, but still good---------------------------
    
    
    
    CBulbBlock  = CBulbBlock == ia(0.1,1.)
    CBulbPrismatic = CBulbPrismatic == ia(0.1,1.)
    rgp.add_one_rule(BulbLength,'CBulbBlock')
    rgp.add_one_rule(CBulbPrismatic,'CBulbPrismatic')
    
    rgp.compute_fresh_rules_graph()
    
    
    #single space left, but still good ------------------------
    
    
    ##
    ##************************************* end bulb rules
    ##
    
    
    ##
    ##************************************* Ship to Bulb Coefficients
    ##
    #
    #linear
    #-----------------------------------------
    Cbb    = lp.PStates(name='Cbb')
    Clpr   = lp.PStates(name='Clpr')
    Czb    = lp.PStates(name='Czb')
    #
    #nonlinear
    #-----------------------------------------
    Cabt   = lp.PStates(name='Cabt')
    Cabl   = lp.PStates(name='Cabl')
    Cvpr   = lp.PStates(name='Cvpr')
    
    #still goood-------------------------------------------------
    
    
    
    #
    ##
    ##************************************* end Ship to Bulb Coefficients
    ##
    
    
    ##
    ##************************************* nonlinear relations
    ##
    #Kracht nonlinear relations:
    Cabt    = Cabt ==  A_mid/ship_Amshp
    Cabl    = Cabl ==  A_lateral/ship_Acp #departure from Kracht
    Cvpr    = Cvpr ==  BulbVolume/ship_Vol
    #
    rgp.add_one_rule( Cbb, 'Cbb')
    rgp.add_one_rule(Clpr, 'Clpr')
    rgp.add_one_rule( Czb, 'Czb')
    #"""
    #
    rgp.compute_fresh_rules_graph()
    #    
    #    
    #still goood!------------------------------------
    
    #
    ##
    ##************************************* end nonlinear relations
    ##
    
    ##
    ##************************************* linear relations
    ##
    #Kracht linear relations:
    Cbb     = Cbb  ==  BulbBeam/ship_beam
    Clpr    = Clpr ==  BulbLength/ship_Lpp
    Czb     = Czb  ==  BulbLength/ship_depth
    #
    #"""
    rgp.add_one_rule( Cbb, 'Cbb')
    rgp.add_one_rule(Clpr, 'Clpr')
    rgp.add_one_rule( Czb, 'Czb')
    #"""
    #
    rgp.compute_fresh_rules_graph()
    #
    ##
    ##************************************* end linear relations
    ##

    g
    
    def doit(ellist=None,
             sobol_seq=None):
        """
        A_mid = A_mid == ia(14.376,19.99)
        A_mid = A_mid == ia(14.4,19.8)
        A_mid = A_mid == ia(14.8,19.4)
        rgp.add_one_rule(A_mid)
        """
        
        if ellist is None:
            ellist = [BulbDepth,
                      CBulbMid,
                      CBulbBlock,
                      BulbBeam,
                      BulbLength,
                      BulbVolume,
                      #A_mid,
                      #A_lateral,
                      #A_BBwl,
                      CBulbCtrPln,
                      CBulbWtrPln,
                      CBulbPrismatic]
            
        if sobol_seq is None:
            sobol_seq = sobol.sobolSeq([1,1],[1,1])
            
        ith = -1
        #x=.5
        redo = []
        for i in range(len(ellist)):
            #x = sobol_seq.next()
            x = random.random()
            var = ellist[i]
            mk_name, mk_val = rgp.get_name_and_value(var)
            print mk_name, mk_val
            val = mk_val[ith].getpoint(x)
            value = ia(val,val)
            var = var == value
            #-----------------------------------------
            ret_state = copy.copy(rgp.env)
            rgp.add_one_rule(var,var.name)
            rgp.compute_fresh_rules_graph()
            #-----------------------------------------
            
            if len(rgp.env.states)==0:
                rgp.reset(ret_state,var)
                redo.append(var)
            else:
                if rgp._verbose: print 'done ', mk_name,'=>', value
        return redo, sobol_seq
    
    def checker():
        print '\n CBulbMid'
        print gv(A_mid) / (gv(BulbBeam)*gv(BulbDepth))
        print gv(CBulbMid)
        
        print '\n CBulbCtrPln'
        print gv(A_lateral)/(gv(BulbLength)*gv(BulbDepth))
        print gv(CBulbCtrPln)
        
        print '\n CBulbWtrPln'
        print gv(A_BBwl)/(gv(BulbLength)*gv(BulbBeam))
        print gv(CBulbWtrPln)
        
        print '\n CBulbBlock'
        print gv(BulbVolume)/(gv(BulbLength)*gv(BulbBeam)*gv(BulbDepth))
        print gv(CBulbBlock)
        
        print '\n CBulbPrismatic'
        print gv(BulbVolume)/(gv(BulbLength)*gv(A_mid))
        print gv(CBulbPrismatic)
        return
    
    
    def gv(var):
        return rgp.get_name_and_value(var)[1][0]
    
    
    
    redo, sbs = doit()
    mm = 0
    while len(redo)>0 and mm<1:
        redo, sbs = doit(redo, sbs)
        mm+=1
    
    
    checker()
    print 'mm = ',mm
    return redo, rgp

def older_test():
    print 'Kraft Bulbous Bow Parameters'
    print'\n Linear Parameters:'
    
    
    
    
    ##
    ##************************************* bulb
    ##
    A_mid = lp.PStates(name='A_mid')
    A_lateral = lp.PStates(name='A_lateral')
    #A_flat = lp.PStates(name='A_flat')  #BBwl instead
    A_BBwl = lp.PStates(name='A_BBwl')
    
    
    #use the principle of the minimum square enclosing box (cartesian b/c ia is cartesian)
    BulbBeam    = lp.PStates(name = 'BulbBeam')     #Bulb half beam
    BulbDepth   = lp.PStates(name = 'BulbDepth')    #Bulb depth
    BulbLength  = lp.PStates(name = 'BulbLength')   #Bulb max length (min square enclosing box length)
    
    
    BulbVolume  = lp.PStates(name = 'BulbVolume') 
        
    
    CBulbMid = lp.PStates(name = 'CBulbMid') #Bulb midsection coefficient
    CBulbCtrPln = lp.PStates(name = 'CBulbCtrPln') #Bulb centerplane profile area coefficient
    CBulbWtrPln = lp.PStates(name = 'CBulbWtrPln') #Bulb waterplane area coefficient
    
    CBulbBlock     = lp.PStates(name = 'CBulbBlock')
    CBulbPrismatic = lp.PStates(name = 'CBulbPrismatic') 
    
    
    
    
    #TODO: fix this with class to construct rules graph!
    
    A_mid = A_mid == ia(10.,20.) #issue: order of assignment "matters"
    BulbBeam = BulbBeam == ia(5.,10.)
    CBulbMid = CBulbMid == ia(.5,1.)
    CBulbWtrPln = CBulbWtrPln == ia(.5,1.)
    #BulbDepth = BulbDepth == ia(-10.1,10.)
    
    
    rgp = RulesGraphProcessor()
    #"""
    
    #rgp.add_one_rule(BulbDepth,'BulbDepth')
    rgp.add_one_rule(A_mid,'A_mid')
    rgp.add_one_rule(BulbBeam,'BulbBeam')
    rgp.add_one_rule(CBulbMid,'CBulbMid')
    rgp.add_one_rule(CBulbWtrPln,'CBulbWtrPln')
    #rgp.add_one_rule(BulbDepth,'BulbDepth') #should not older_testbe needed!
    
    rgp.compute_fresh_rules_graph()
    
    #rgp.AC_revise()
    #rgp.env
    
    ##
    ##************************************* end bulb
    ##
    
    ##
    ##************************************* ship
    ##
    ship_beam = lp.PStates(name='ship_beam') 
    # note: could use bare_hull var names instead. 
    #       e.g. lp.PStates(name=self.init_ship_beam.name)
    ship_depth = lp.PStates(name='ship_depth')
    ship_Lpp = lp.PStates(name='ship_Lpp')
    #
    #quantities of m**2
    ship_Amshp = lp.PStates(name='ship_Amshp')
    ship_Acp = lp.PStates(name='ship_Acp')
    ship_Awp = lp.PStates(name='ship_Awp')
    #
    #quantities of m**3
    ship_Vol = lp.PStates(name='ship_Vol')
    
    
    #-----------------------------------------
    # set the ship values in the bulb environement:
    ship_beam   = ship_beam == ia(17.4663142374, 17.4663142374)
    
    ship_depth  = ship_depth == ia(16.2051841085, 16.2051841085)
    
    ship_Lpp    = ship_Lpp == ia(111.099919763, 111.099919763)
    
    ship_Amshp  = ship_Amshp == ia(261.639572047, 261.639572047)
    
    ship_Acp    = ship_Acp == ia(1656.36308186, 1656.36308186)
    
    ship_Awp    = ship_Awp == ia(1736.75296874, 1736.75296874)
    
    ship_Vol    = ship_Vol == ia(27043.7825521, 27043.7825521)
    #-----------------------------------------
    
    
    
    #-----------------------------------------
    rgp.add_one_rule(ship_beam,'ship_beam')
    rgp.add_one_rule(ship_depth,'ship_depth')
    rgp.add_one_rule(ship_Lpp,'ship_Lpp')
    rgp.add_one_rule(ship_Amshp,'ship_Amshp')
    rgp.add_one_rule(ship_Acp,'ship_Acp')
    rgp.add_one_rule(ship_Awp,'ship_Awp')
    rgp.add_one_rule(ship_Vol,'ship_Vol')
    #-----------------------------------------
    ##
    ##************************************* end ship
    ##
    rgp.compute_fresh_rules_graph()
    
    
    
    ##
    ##************************************* bulb rules
    ##
    """-----------------------------------------------
    Rule:  Midbulb_Area < max_Beam * max_Depth
    CBulbMid -> [0.,1.]
    CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
    """
    CBulbMid = CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
    
    """-----------------------------------------------
    Rule: z-y_area < max_length * max_Depth
    """
    CBulbCtrPln = CBulbCtrPln == A_lateral/(BulbLength*BulbDepth)
    
    """-----------------------------------------------
    Rule: wL_area < max_length * max_Depth
    """
    #CBulbWtrPln = CBulbWtrPln == A_BBwl/(BulbLength*BulbDepth)
    CBulbWtrPln = CBulbWtrPln == A_BBwl/(BulbLength*BulbBeam)
    
    
    
    #2 more rules!    
    """-----------------------------------------------
    Rule: Bulb_vol < max_length * max_Depth * max *BulbBeam
    """
    CBulbBlock = CBulbBlock == BulbVolume/(BulbLength*BulbDepth*BulbBeam)
    #
    """-----------------------------------------------
    Rule: Bulb_vol < max_length * mid_bulb_area
    """
    CBulbPrismatic = CBulbPrismatic == BulbVolume/(BulbLength*A_mid)    
    #


    
    rgp.add_one_rule(CBulbMid,'CBulbMid')
    rgp.add_one_rule(CBulbCtrPln,'CBulbCtrPln')
    rgp.add_one_rule(CBulbWtrPln,'CBulbWtrPln')
    
    
    rgp.add_one_rule(CBulbBlock, 'CBulbBlock')
    rgp.add_one_rule(CBulbPrismatic, 'CBulbPrismatic')
    
    rgp.compute_fresh_rules_graph()
    #"""
    
    #BulbLength  = BulbLength == ia(10.,15.)
    BulbLength  = BulbLength == ia(1.,15.)
    #CBulbCtrPln = CBulbCtrPln == ia(.5,1.)
    CBulbCtrPln = CBulbCtrPln == ia(0.,1.)
    rgp.add_one_rule(BulbLength,'BulbLength')
    rgp.add_one_rule(CBulbCtrPln,'CBulbCtrPln')
    
    rgp.compute_fresh_rules_graph()
    
    #BulbDepth = BulbDepth == ia(1.,10.)
    BulbDepth = BulbDepth == ia(5.,10.)
    rgp.add_one_rule(BulbDepth,'BulbDepth')
    
    
    
    CBulbBlock  = CBulbBlock == ia(0.,1.)
    CBulbPrismatic = CBulbPrismatic == ia(0.,1.)
    rgp.add_one_rule(BulbLength,'CBulbBlock')
    rgp.add_one_rule(CBulbPrismatic,'CBulbPrismatic')
    
    rgp.compute_fresh_rules_graph()
    
    ##
    ##************************************* end bulb rules
    ##
    
    
    ##
    ##************************************* Ship to Bulb Coefficients
    ##
    #
    #linear
    #-----------------------------------------
    Cbb    = lp.PStates(name='Cbb')
    Clpr   = lp.PStates(name='Clpr')
    Czb    = lp.PStates(name='Czb')
    #
    #nonlinear
    #-----------------------------------------
    Cabt   = lp.PStates(name='Cabt')
    Cabl   = lp.PStates(name='Cabl')
    Cvpr   = lp.PStates(name='Cvpr')
    #
    ##
    ##************************************* end Ship to Bulb Coefficients
    ##
    
    
    ##
    ##************************************* nonlinear relations
    ##
    #Kracht nonlinear relations:
    Cabt    = Cabt ==  A_mid/ship_Amshp
    Cabl    = Cabl ==  A_lateral/ship_Acp #departure from Kracht
    Cvpr    = Cvpr ==  BulbVolume/ship_Vol
    #
    #"""
    rgp.add_one_rule( Cbb, 'Cbb')
    rgp.add_one_rule(Clpr, 'Clpr')
    rgp.add_one_rule( Czb, 'Czb')
    #"""
    #
    rgp.compute_fresh_rules_graph()
    #
    ##
    ##************************************* end nonlinear relations
    ##
    
    ##
    ##************************************* linear relations
    ##
    #Kracht linear relations:
    Cbb     = Cbb  ==  BulbBeam/ship_beam
    Clpr    = Clpr ==  BulbLength/ship_Lpp
    Czb     = Czb  ==  BulbLength/ship_depth
    #
    #"""
    rgp.add_one_rule( Cbb, 'Cbb')
    rgp.add_one_rule(Clpr, 'Clpr')
    rgp.add_one_rule( Czb, 'Czb')
    #"""
    #
    rgp.compute_fresh_rules_graph()
    #
    ##
    ##************************************* end linear relations
    ##
    print '\n\n state after:'
    #print st
    
    print rgp.env
    
    
    
    def doit(ellist=None,
             sobol_seq=None):
        
        if ellist is None:
            ellist = [BulbDepth,
                      CBulbMid,
                      CBulbBlock,
                      BulbBeam,
                      BulbLength,
                      BulbVolume,
                      #A_mid,
                      #A_lateral,
                      #A_BBwl,
                      CBulbCtrPln,
                      CBulbWtrPln,
                      CBulbPrismatic]
        if sobol_seq is None:
            sobol_seq = sobol.sobolSeq([1,1],[1,1])
            
        ith = -1
        #x=.5
        redo = []
        for i in range(len(ellist)):
            x = sobol_seq.next()
            var = ellist[i]
            mk_name, mk_val = rgp.get_name_and_value(var)
            print mk_name, mk_val
            val = mk_val[ith].getpoint(x)
            value = ia(val,val)
            var = var == value
            #-----------------------------------------
            ret_state = copy.copy(rgp.env)
            rgp.add_one_rule(var,var.name)
            rgp.compute_fresh_rules_graph()
            #-----------------------------------------
            
            if len(rgp.env.states)==0:
                rgp.reset(ret_state,var)
                redo.append(var)
            else:
                if rgp._verbose: print 'done ', mk_name,'=>', value
        return redo, sobol_seq
    
    def checker():
        print '\n CBulbMid'
        print gv(A_mid) / (gv(BulbBeam)*gv(BulbDepth))
        print gv(CBulbMid)
        
        print '\n CBulbCtrPln'
        print gv(A_lateral)/(gv(BulbLength)*gv(BulbDepth))
        print gv(CBulbCtrPln)
        
        print '\n CBulbWtrPln'
        print gv(A_BBwl)/(gv(BulbLength)*gv(BulbBeam))
        print gv(CBulbWtrPln)
        
        print '\n CBulbBlock'
        print gv(BulbVolume)/(gv(BulbLength)*gv(BulbBeam)*gv(BulbDepth))
        print gv(CBulbBlock)
        
        print '\n CBulbPrismatic'
        print gv(BulbVolume)/(gv(BulbLength)*gv(A_mid))
        print gv(CBulbPrismatic)
        return
    
    
    def gv(var):
        return rgp.get_name_and_value(var)[1][0]
    
    
    redo, sbs = doit()
    mm = 0
    while len(redo)>0 and mm<10:
        redo, sbs = doit(redo, sbs)
        mm+=1
    
    
    checker()
    print 'mm = ',mm
    return redo, rgp



    

def get_curve(xb,yb,
              xe,ye,
              alphab=None,alphae=None,
              cb=None,ce=None,
              area=None,
              xc=None,yc=None,
              nCV=None):
    """
        given : the basic form parameters
        output : a curve satisfying them
    """
    data_store = FPDIdeal(xb,
                             yb,
                             xe,
                             ye,
                             alphab,
                             alphae,
                             cb,
                             ce,
                             area,
                             xc,
                             yc,
                             nCV)
    
    return data_store.make_curve()




if __name__ == "__main__":
    #pass
    #redo, rgp = older_test()
    #redo1, rgp1 = small_test()
    
    #"""
    self = GeometryGenerator()
    sobol_seq = self.sobol_seq
    self.rgp._verbose = True
    #    k1 = self.Primary
    #    k2 = self.Coefficients
    #    k3 = self.Areas
    self.tree_search()
    #"""
    #print self.rgp.env
    
    #type(self.rgp.nodes[0].vars['ship_beam'])
    
    
    """
    
    self.CBulbMid = self.CBulbMid ==  self.A_mid/(
            self.BulbBeam*self.BulbDepth)
    self.rgp.add_one_rule(self.CBulbMid,
                          self.CBulbMid.name)
    self.rgp.compute_fresh_rules_graph()
    #"""
    
    self.get_curves()
    
    self.inib_ideal.solve_curve_ini()
    self.inib_ideal.solve_curve_fini()#self.inib_ideal.xc)
    self.inib_ideal.Lspline.curve.plotcurve_detailed()
    
    self.mb_ideal.solve_curve_ini()
    self.mb_ideal.solve_curve_fini()#self.mb_ideal.xc)
    self.mb_ideal.Lspline.curve.plotcurve_detailed()


    self.fsb_ideal.solve_curve_ini()
    self.fsb_ideal.solve_curve_fini()#self.fsb_ideal.xc)
    self.fsb_ideal.Lspline.curve.plotcurve_detailed()
    
    self.bwl_ideal.solve_curve_ini()
    self.bwl_ideal.solve_curve_fini()
    self.bwl_ideal.Lspline.curve.plotcurve_detailed()
    
    self.cp_ideal.solve_curve_ini()
    self.cp_ideal.solve_curve_fini()#self.cp_ideal.xc)
    self.cp_ideal.Lspline.curve.plotcurve_detailed()
    
    #self.specialized_bow_development_nose()
    self.make_bulb()
    
    
    print 'Xc = ',self.mb_ideal.Lspline.curve.Xc
    print 'Xc = ',self.bwl_ideal.Lspline.curve.Xc
    print 'Xc = ',self.cp_ideal.Lspline.curve.Xc
    
    
    
    