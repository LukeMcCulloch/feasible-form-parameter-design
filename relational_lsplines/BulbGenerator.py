#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 12 02:05:36 2018

@author: luke
"""

import copy
import sqKanren as lp
np = lp.np
import sobol
ia = lp.ia #use ia import from lp instead to match isinstance
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
                 bare_hull=None):
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
        
        self.linear_parameters()
        self.nonlinear_parameters()
        
        
        self.initialize_lists()
        self = ini_coeffs(self)#.Coefficients, self.rgp)
        
        self.experiment_bulb_parameters() #add new rules here
        self.basic_bulb_rules()
        
        
        
    
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
        
        ideal = BoundBulbCurves(BulbBeam,
                                BulbDepth,
                                BulbLength)
        self.ideal = ideal
        
        """
        # INI-BULB CURVE:
        # interface curve:
        TODO: redo these in a managable way!
        """
        r1m,r2m = ideal.simple_projection('z','y')
        zordinate = ideal.coords2D['x']
        self.inib_ideal = FPDIdeal(
                xb = ideal.v_inibc_cL_bL[r1m],
                xe = ideal.v_inibc_cL_tp[r1m],
                yb = ideal.v_inibc_cL_bL[r2m],
                ye = ideal.v_inibc_cL_tp[r2m],
                #alphab = 0.,
                #alphae = 0.,
                #Cab_given = 0.,
                #Cae_given = 0.,
                ymax = BulbBeam,
                area = self.get_for_design(self.A_mid),
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
        self.mb_ideal = FPDIdeal(
                        xb = ideal.v_mbc_cL_bL[r1m],
                        xe = ideal.v_mbc_cL_tp[r1m],
                        yb = ideal.v_mbc_cL_bL[r2m],
                        ye = ideal.v_mbc_cL_tp[r2m],
                        #alphab = 0.,
                        #alphae = 0.,
                        #Cab_given = 0.,
                        #Cae_given = 0.,
                        ymax = BulbBeam,
                        area = self.get_for_design(self.A_mid),
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
                            #Cab_given = 0.,
                            #Cae_given = 0.,
                            #ymax = thisymax,
                            ymax = BulbBeam,
                            ##not used #area = np.sqrt(self.get_for_design(self.A_mid)/np.pi)
                            #area = 0.5*(np.pi*(thisymax**2)), #was using this!
                            area = self.get_for_design(self.A_mid),
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
    
    
        
    def make_nose(self, yloc = None, Wratio=.8,Aratio=.7):
        """from GeometryGenerator import FPDIdeal
        """
        ytot = self.idealdict['cp_ideal'].ymax
        xb = self.idealdict['cp_ideal'].xb
        xe = self.idealdict['cp_ideal'].xe
        xdiff = xe-xb
        ymax = self.idealdict['cp_ideal'].ymax
        if yloc is None:
            yloc = .95*ymax
        xb = xb+.5*(1.-Wratio)*xdiff
        xe = xe-.5*(Wratio)*xdiff
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
            nnet.append(uopt.interpolate_vertices(
                    initial_vertices=pts,
                    helper_curve = curve,
                    return_Lspline = False,
                    make_thb=False) )
        
        lcurvenet = nnet
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
        
        self.wrapped_search(k2,sobol_seq)
        self.wrapped_search(k1,sobol_seq)
        self.wrapped_search(k3,sobol_seq)
        
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
            if self.rgp._verbose: print 'done ', mk_name,'=>', value
            self.rgp.varsmap[var.name].flist.pop()
            return True #,  self.rgp.env
    
    

    
    def narrow_this(self, var, ith=-1):
        """
            var = key
            ith=-1
            
        """
        x1=.2
        x2=.8
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
                assert(len(mk_val)==1),'ERROR, singularities still present'
                iaval = mk_val[0]
                diff = iaval.sup - iaval.inf
                if np.linalg.norm(diff)<self.tol:
                    pass
                else:
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
        Rule:  Midbulb_Area < max_Beam * max_Depth
        CBulbMid -> [0.,1.]
        CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
        """
        
        self.CBulbMid = self.CBulbMid ==  self.A_mid/(self.BulbBeam*self.BulbDepth)
        #
        """-----------------------------------------------
        Rule: z-y_area < max_length * max_Depth
        """
        self.CBulbCtrPln = self.CBulbCtrPln == self.A_lateral/(self.BulbLength*self.BulbDepth)
        #
        """-----------------------------------------------
        Rule: wL_area < max_length * max half-Beam
        """
        self.CBulbWtrPln = self.CBulbWtrPln == self.A_BBwl/(self.BulbLength*self.BulbBeam)
        #
        """-----------------------------------------------
        Rule: Bulb_vol < max_length * max_Depth * max *BulbBeam
        """
        self.CBulbBlock = self.CBulbBlock == self.BulbVolume/(self.BulbLength*self.BulbDepth*self.BulbBeam)
        #
        """-----------------------------------------------
        Rule: Bulb_vol < max_length * mid_bulb_area
        """
        self.CBulbPrismatic = self.CBulbPrismatic == self.BulbVolume/(self.BulbLength*self.A_mid)
        #
        self.rgp.add_one_rule(   self.CBulbMid, self.CBulbMid.name  )
        self.rgp.add_one_rule(self.CBulbCtrPln, 'CBulbCtrPln')
        self.rgp.add_one_rule(self.CBulbWtrPln, 'CBulbWtrPln')
        self.rgp.add_one_rule(self.CBulbBlock, 'CBulbBlock')
        self.rgp.add_one_rule(self.CBulbPrismatic, 'CBulbPrismatic')
        #
        self.rgp.compute_fresh_rules_graph()
        #
        #self.rgp = rgp
        return #rgp
    
    
    def linear_parameters(self):#, rgp=None):
        """Kracht Design of Bulbous Bows, 1978
        """
        #if rgp is None:
        #    rgp = self.rgp
        
        #Kracht linear relations:
        self.Cbb     = self.Cbb  ==  self.BulbBeam/self.ship_beam
        self.Clpr    = self.Clpr ==  self.BulbLength/self.ship_Lpp
        self.Czb     = self.Czb  ==  self.BulbLength/self.ship_depth
        #
        self.rgp.add_one_rule( self.Cbb)
        self.rgp.add_one_rule(self.Clpr)
        self.rgp.add_one_rule( self.Czb)
        #
        self.rgp.compute_fresh_rules_graph()
        #
        #self.rgp = rgp
        return #rgp
    
    
    
    
    def nonlinear_parameters(self):#, rgp=None):
        """Kracht Design of Bulbous Bows, 1978
        """
        #Kracht nonlinear relations:
        self.Cabt    = self.Cabt ==  self.A_mid/self.ship_Amshp
        self.Cabl    = self.Cabl ==  self.A_lateral/self.ship_Acp #departure from Kracht
        self.Cvpr    = self.Cvpr ==  self.BulbVolume/self.ship_Vol
        #
        self.rgp.add_one_rule( self.Cabt)
        self.rgp.add_one_rule(self.Cabl)
        self.rgp.add_one_rule( self.Cvpr)
        #
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
        self.A_mid = self.A_mid == ia(2.,15.)
        
        #self.BulbBeam = self.BulbBeam == ia(1.,2.5)
        self.BulbBeam  = self.BulbBeam  == self.lowhigh(.1,.14,self.ship_beam)
        #
        self.CBulbMid = self.CBulbMid == ia(.5,1.)
        self.CBulbWtrPln = self.CBulbWtrPln == ia(.7,.97)
        #
        #BulbDepth = BulbDepth == ia(-10.1,10.)
        self.rgp.add_one_rule(self.A_mid,'A_mid')
        self.rgp.add_one_rule(self.BulbBeam,'BulbBeam')
        self.rgp.add_one_rule(self.CBulbMid,'CBulbMid')
        self.rgp.add_one_rule(self.CBulbWtrPln,'CBulbWtrPln')
        
        self.rgp.compute_fresh_rules_graph()
        
        
        
        
        #self.BulbLength     = self.BulbLength == ia(10.,15.)
        self.BulbLength     = self.BulbLength == self.lowhigh(0.03,0.1,self.ship_Lpp)
        #
        #self.BulbDepth      = self.BulbDepth == ia(2.5,10.)
        self.BulbDepth  = self.BulbDepth  == self.lowhigh(.3,.4,self.ship_depth)
        #
        self.CBulbCtrPln    = self.CBulbCtrPln == ia(.5,1.)
        #now add them
        self.rgp.add_one_rule(self.BulbLength,'BulbLength')
        self.rgp.add_one_rule(self.BulbDepth,'BulbDepth')
        self.rgp.add_one_rule(self.CBulbCtrPln,'CBulbCtrPln')
        
        self.rgp.compute_fresh_rules_graph()
        
        #
        ##
        ##************************************* new bulb blending rules
        ##
        
        #keep the bulb short:
        self.Clpr    = self.Clpr ==  ia(.01,.045)
        self.rgp.add_one_rule(self.Clpr, 'Clpr')
        
        
        self.rgp.compute_fresh_rules_graph()
        
        #
        ##
        ##************************************* end rules
        ##
        
        self.rgp.compute_fresh_rules_graph()
        
        return #rgp
    
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



if __name__ == "__main__":
    pass