# -*- coding: utf-8 -*-
"""
Created on Sat May 31 13:05:42 2014

@author: lukemcculloch

-----------------------------------
  Hull Class
-----------------------------------

June 17, 2014
Things to do:
 - fore body builder
 - interface to hold form parameters to match mid longitudinals with fore body longitudinals
     * if the transverse curve vertices at the interface match, then the surface will be continuous. 
         ' given that all longitudinal curves inerpolate those transverse vertices.
         ' similarly:
     * if the tangents of the longitudinal match, then the slope of the surface will be continuous.
     * if the curvature of the longitudinals matches, then the curveature of the surface will match
 - form parameter curves that specify properties across hull segments - so as to ensure continuity.


Agents:

    Type            Performance Measure     Environment                     actuators               Sensors

    Bulbous Bow     Resistance              Ship Hull adjacent curve(s)     vertex placement        Nonlinear Resistance Code
                                            (interface)                     Kracht 7 parameters     Ship Hull adjacent curve(s)

    Main Hull       Resistance              Water                           vertices                WAMIT, NLpanel code
                    Seakeeping              Free Surface                    Form Curves             producability algorithms
                    Stability               canals, dry docks etc.          transverse curves       inequality constraints?
                    Space                   Operations                      form parameters
                    Structural Force
                    Expense
    
    Above wL        space, wave behavior    Free Surface                    interface (adjacent cv) inequality constraints
                                            Operations
    
    Transom         resistance              Free Surface                    vertices                NLpanel code
                    seakeeping?             Froude Number                   tangency
                    etc?                                                    interface


Agents of Agents:

    Form Curves     Heuristics?          <- Design requirements             form parameters        data (collection?) on goodness? 
                                                                            curve vertices


Propogation of Information:
    TODO

"""
#python standards:
import numpy as np
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=300)
import matplotlib.pyplot as plt
import copy
#from mayavi.mlab import * 

#import TLM B-Spline Classes
from curve import Bspline               as Bspline     
from curve import interpolatedBspline   as interpolatedBspline
from curve import curveSet
from curve import BsplineSurface
from initialValues          import InitializeControlPoints, InitializeControlVertices
from Equalityconstraints    import FormParameter
from InequalityConstraints  import InequalityConstraint
from interfaceGUItoOpt      import interface_ADLsplineOpt
from ADL_LS_OptBspline      import LagrangeSpline, ADLsplineOpt
from CustomRange            import customrange
# import TLM 2D optimization with Least Squares Min included as Optimum Criteria:


#point interpolation routines:
from minDistance import naive_minimizer

#import TLM fuzzy logic
from fich import get_TSK_SAC



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

def transpose_vertices(vertices, dim, scalar):
    vertices[:,dim] += scalar
    return vertices

def FormParameters():
   return  

class parametric_curve(Bspline):
    """Not used - over thought it - may revisit..."""
    def __init__(self, NDP = 1.0, FormParameters={}):
        self.NDP     = NDP #non-dimensional parameter
        self.FPnames = {0:'Vertices',
                        1:'TangentAngles',
                        2:'Area',
                        3:'XCentroidValue',
                        4:'YCentroidValue',
                        5:'Curvature'}
        self.form_parameters()
    
    def form_parameters(self):
        self.vblDict = {}
        return
        
    def ini_curve():
        return
    
    def make_curve():
        interface_ADLsplineOpt(self.vblDict)
        return
        
    def HarriesControlPoints(self,   xb=0., 
                                     yb=0., 
                                     alphab=0., 
                                     Cab_given=0.15,
                                     xe=12.,
                                     ye=12.,
                                     alphae  = 0.,
                                     Cae_given=0.1,
                                     area    = 55.,
                                     xc      = 8.,
                                     yc      = 4.):
        
        """
        # We will have a 7 control point curve, Order 4, as recommended by Harries, 1999.
        # Given The Following Form Parameters:
        
        ## Curve Differential parameters:
        
        xb=0.               #1st vertex: x value
        yb=0.               #1st vertex: y value
        alphab  = 0.        #tangent at the start of the curve
        Cab_given=0.15      #curvature desired at the beggining
        
        xe=12.              #last vertex: x value
        ye=12.              #last vertex: y value
        alphae  = 0.        #tangent at the end of the curve
        Cae_given=0.1       #curvature desired at the end
        
        ## Curve Integral parameters:
        
        area    = 72.
        xc      = 8.
        yc      = 4.
        """
        
        #vertices = InitializeControlPoints(xb,yb,alphab,Cab_given,xe,ye,alphae,Cae_given,xc,yc,area)
        vertices = InitializeControlVertices(xb,yb,xe,ye,alphab,alphae,Cab_given,Cae_given,area,xc,yc,7)
        vertices = vertices.vertices
        return vertices

    
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
        self.return_curve = Bspline(return_vertices, self.k, self.nump)
        return 
    
    def offset(self):
        return_vertices = copy.copy(self.initial_curve.vertices)
        if self.offset_const != None:
            print 'offset' 
            return_vertices[:,self.plane_dim] += self.offset_const
        if self.definition_curve != None:
            return_vertices[:,:] += self.definition_curve.vertices[:,:] - self.initial_curve.vertices[:,:]
        self.return_curve = Bspline(return_vertices, self.k, self.nump)
        return            
        


class Hull(object):

    
    def __init__(self, LDWL, MaxBreadth, MaxDraft, Depth, BareHullVolume, bbow=True, tstern=True):
        self.k      = 4
        self.n      = 7
        self.nump   = 50
        self.tol    = 1.e-3
        self.bbow   = bbow
        self.tstern = tstern
        self.principle_particulars(LDWL, MaxBreadth, MaxDraft, Depth, BareHullVolume)
        self.definestations()
        self.derived_particulars()
        self.block_coefficient = self.BareHullVolume/(self.LengthDWL * self.MaxBreadth * self.MaxDraft)
        self.basic_curves() #TLMflag
        self.derived_curves()
        self.make_hull()
        
    def definestationsOLD(self):
        self.stations = {}
        self.stations['bow_shape']      = [.10*self.LengthDWL,.15*self.LengthDWL]
        self.stations['flat_bottom']    = [.43*self.LengthDWL,.71*self.LengthDWL]
        self.stations['flat_side']      = [.27*self.LengthDWL,.45*self.LengthDWL]
        self.stations['flat_DWL']       = [.43*self.LengthDWL,.71*self.LengthDWL]
        return
    def definestations(self):
        self.stations = {}
        self.stations['bow_shape']      = [.10*self.LengthDWL,.15*self.LengthDWL]
        self.stations['flat_bottom']    = [.25*self.LengthDWL,.71*self.LengthDWL]
        self.stations['flat_side']      = [.27*self.LengthDWL,.45*self.LengthDWL]
        self.stations['flat_DWL']       = [.43*self.LengthDWL,.71*self.LengthDWL]
        return
        
    def make_hull(self):
        #num = 5
        self.buildtransversecurves(5)
        self.buildlongitudinalcurves()
        self.build_primary_hull_surface()
        return
        
    def principle_particulars(self, LDWL, MaxBreadth, MaxDraft, Depth, BareHullVolume):
        self.LengthDWL      = LDWL
        self.MaxBreadth     = MaxBreadth
        self.MaxDraft       = MaxDraft
        self.Depth          = Depth
        self.BareHullVolume = BareHullVolume
        return
    
    def derived_particulars(self):
        self.WLarea     = .9*self.LengthDWL*self.MaxBreadth
        self.maxSecArea = 1.5*self.BareHullVolume/self.LengthDWL#1.*(.1*self.LengthDWL)**2
        return
        
    def basic_curves(self):
        self.defineSAC()
        self.defineDWL()
        self.augDWL()
        return
    
    def derived_curves(self):
        self.bow_curve()
        self.compute_midship_section()
        self.computeFlatOfSide_curves(self.stations['flat_side'][0], self.stations['flat_side'][1])
        self.compute_stern_curve()
        return
        
    def bow_curve(self):
        vxyz = linear_vertices((0.,0.,0.),(0.0,self.Depth,-.1*self.LengthDWL),7)
        self.bowcurves = [Bspline(vxyz, self.k, self.nump)]
        self.bulbcurves = []
        return
    
    def add_bulbous_bow(self, bbow):
        #copy bbow aft (i.e. hull interface) verices 3 times [C2 coninuity]
        curve1 = copy.copy(bbow.tCurveNet3D[-1].vertices)
        curve2 = copy.copy(bbow.tCurveNet3D[-1].vertices)
        curve3 = copy.copy(bbow.tCurveNet3D[-1].vertices)
        curve1[:,-1]=0.
        curve2[:,-1]=0.
        curve3[:,-1]=0.
        curve3 = transpose_vertices(curve3, 2,.1*bbow.length)
        curve2 = transpose_vertices(curve2, 2,.15*bbow.length)
        curve1 = transpose_vertices(curve1, 2,.2*bbow.length)
        curve3 = Bspline(curve3, self.k, self.nump)
        curve2 = Bspline(curve2, self.k, self.nump)
        curve1 = Bspline(curve1, self.k, self.nump)
        self.bowcurves = [curve1,curve2,curve3]
        #self.bowcurves = [curve2,curve3]
        self.bulbcurves = bbow.tCurveNet3D
        self.make_hull()
        return
        
    def sub_bulbous_bow(self):
        self.bow_curve()
        self.make_hull()
        return
        
    def compute_midship_section(self):
        bottom  = linear_vertices((0.,0.,self.LengthDWL*0.5),(self.MaxBreadth,0.,self.LengthDWL*0.5),4)
        side    = linear_vertices((self.MaxBreadth,0.,self.LengthDWL*0.5),(self.MaxBreadth,self.MaxDraft,self.LengthDWL*0.5),3)
        self.midshipsec = np.zeros((len(bottom)+len(side),3),float)
        self.midshipsec[:len(bottom)]=bottom
        self.midshipsec[len(bottom):]=side
        self.midshipsec = Bspline(self.midshipsec,self.k,self.nump)
        return
        
    def computeFlatOfSide_curves(self, fore=27., aft = 45.):
        #def computeFlatOfSide_curves(self, self.stations['flat_side'][0], self.stations['flat_side'][1])
        fwd_extent = fore   #self.LengthDWL*fore
        aft_extent = aft    #self.LengthDWL*aft
        bottom      = linear_vertices((0.,0.,aft_extent),(self.MaxBreadth,0.,aft_extent),4)
        side        = linear_vertices((self.MaxBreadth,0.,aft_extent),(self.MaxBreadth,self.MaxDraft,fwd_extent),3)
        self.flat_bottom = bottom
        self.flatside = {}
        self.flatside['fwd'] = fwd_extent
        #self.flatbottom['extent'] = [.45,.7]
        self.fos1   = np.zeros((len(bottom)+len(side),3),float)
        self.fos1[:len(bottom)]=bottom
        self.fos1[len(bottom):]=side
        self.fos1   = Bspline(self.fos1,self.k,self.nump)
        #self.curve_by_station[str()]
        self.fos2 = Hullcurve('offset', self.fos1,2, -1., None)
        self.fos2 = self.fos2.return_curve
        self.fos3 = Hullcurve('offset', self.fos1,2, -2., None)
        self.fos3 = self.fos3.return_curve
        self.fos_fwd = [self.fos1,self.fos2,self.fos3]
        return
    
    def compute_stern_curve(self):
        stern_depth = 0.5*self.MaxDraft
        bottom  = linear_vertices((0.,stern_depth,self.LengthDWL),(self.MaxBreadth,stern_depth,self.LengthDWL),4)
        side    = linear_vertices((self.MaxBreadth,stern_depth,self.LengthDWL),(self.MaxBreadth,self.MaxDraft,self.LengthDWL),3)
        self.sternsec = np.zeros((len(bottom)+len(side),3),float)
        self.sternsec[:len(bottom)]=bottom
        self.sternsec[len(bottom):]=side
        self.sternsec = Bspline(self.sternsec,self.k,self.nump)
        return
    
    def volume_fwd_curve(self):
        return
    def volume_aft_curve(self):
        return
        
    def template_ship_curve(self):
        return
        
    def shiftSAC(self, fore, aft):
        self.SACrange[1] +=fore
        self.SACrange[2] +=aft
        return
    
    def defineSAC(self):
        num_segments = 3
        areas = [.2,.6,.2]
        #self.SACrange = np.linspace(0.,self.LengthDWL,num_segments+1)
        #self.shiftSAC(10.,5.)
        self.SACrange = np.asarray([0.,self.stations['flat_bottom'][0], self.stations['flat_bottom'][1], self.LengthDWL])
        
        k               = 4
        nump            = 50
        #nondim_factor   = self.LengthDWL
        
        #Aggregate SAC parameters:        
        #area_factor     = 1./10.
        #bulk_factor     = 0.333
        #xm_factor       = 0.5
        #ym_factor       = 0.333
        
    
        if num_segments ==3:
            segment_span = {'fwd':self.SACrange[:2],'mid':self.SACrange[1:3],'aft':self.SACrange[2:]}
            #SAC_segments = []
            
            ## fwd segment:
            length = segment_span['fwd'][1] - segment_span['fwd'][0]
            
            xm_factor = .7
            ym_factor = .2
            xb = segment_span['fwd'][0]
            yb = 0.
            alphab = 0.
            Cab_given = 0.
            xe = segment_span['fwd'][1]
            ye = self.maxSecArea
            alphae = 0.
            Cae_given = 0.
            
            height = ye-yb
            xc = xm_factor*length
            yc = ym_factor*height
            area = areas[0]*self.BareHullVolume #####################EEEEERRORRR
            
            spec = {'SAC':{'maxwidth':self.maxSecArea}}
            spec['xmin'] = xb
            spec['xmax'] = xe
            spec['ymin'] = yb
            spec['ymax'] = self.maxSecArea
            #vertices = InitializeControlPoints(xb,yb,alphab,Cab_given,xe,ye,alphae,Cae_given,xc,yc,area,spec)
            #TLMflag            
            temp_object = InitializeControlVertices(xb,yb,xe,ye,
                                                    alphab,alphae,
                                                    Cab_given,Cae_given,
                                                    area,xc,yc,7,'up',spec)
            #xb=0.,yb=12.,xe=12.,ye=0.,alphab=None,alphae=None,
            #Cab_given=None,Cae_given=None,area=None,
            #xc=None,yc=None,nCV = 7, slope= 'up', options = {}
            
            vertices = temp_object.vertices
            self.SAC_fwd = Bspline(vertices,k,nump)
            self.SAC_fwd.plotcurve()
            self.fwd = Bspline(vertices,k,nump)
            
            FP = FormParameter(self.SAC_fwd)
            IC = InequalityConstraint(self.SAC_fwd)
            

            #FP.add_AngleConstraint(0.,0.) 
            #FP.add_VerticalConstraint(0.,0.)
            
            
            #FP.add_CurvatureConstraint(Cab_given,0.)
            #FP.add_CurvatureConstraint(Cae_given,1.)
            #area = self.SAC_fwd.area.value
            
            #FP.add_XcConstraint(xc)
            
            FP.add_AngleConstraint(0.,1.)
            FP.add_AreaConstraint(area)
            self.SAC_fwd.FormParam = FP.FormParam
            self.SAC_fwd.IneqConstraints = IC.InequalityParam
            data = {'area':area}
            lspline = LagrangeSpline(self.SAC_fwd, 
                                                 {'Inumber':25,
                                                   'ws':10.,
                                                   'w1':.1,
                                                   'w2':.1,
                                                   'w3':.1})#, 10., 'area', None, data)
            """
            lspline = LagrangeSpline(self.SAC_fwd, 
                                                 {'Inumber':25,
                                                   'ws':10.,
                                                   'w1':.1,
                                                   'w2':.1,
                                                   'w3':.1}, None, None, None)
            #"""
                                   
            self.SAC_fwd =lspline.curve
            print 'area of fwd SAC = {}, desired area = {}'.format(self.SAC_fwd.area.value, area)
            self.SAC_fwd.plotcurve()
            
            ##
            ## mid-body SAC:
            ##
            length = segment_span['mid'][1] - segment_span['mid'][0]
            
            xm_factor = .5
            ym_factor = .5
            xb = segment_span['mid'][0]
            yb = self.maxSecArea
            alphab = 0.
            Cab_given = 0.
            xe = segment_span['mid'][1]
            ye = self.maxSecArea
            alphae = 0.
            Cae_given = 0.
            
            height = yb
            xc = xm_factor*length
            yc = ym_factor*height
            area = areas[1]*self.BareHullVolume
            
            #spec = {'SAC':{'maxwidth':self.maxSecArea}}
            spec['xmin'] = xb
            spec['xmax'] = xe
            spec['ymin'] = yb
            spec['ymax'] = self.maxSecArea
            #self.SAC_mid = InitializeControlPoints(xb,yb,alphab,Cab_given,xe,ye,alphae,Cae_given,xc,yc,area,spec)
            temp_object2 = InitializeControlVertices(xb,yb,xe,ye,
                                                     alphab,alphae,
                                                     Cab_given,Cae_given,
                                                     area,xc,yc,7,'up',spec)
            vertices2 = temp_object2.vertices
            self.SAC_mid = Bspline(vertices2,k,nump)
            self.SAC_mid.plotcurve()
            FP = FormParameter(self.SAC_mid)
            IC = InequalityConstraint(self.SAC_mid)
            
            
            self.SAC_mid.FormParam = FP.FormParam
            self.SAC_mid.IneqConstraints = IC.InequalityParam
            
    
            lspline = LagrangeSpline(self.SAC_mid, 
                                                 {'Inumber':25,
                                                   'ws':10.,
                                                   'w1':.1,
                                                   'w2':.1,
                                                   'w3':.1}, None, None, None)
                                   
            self.SAC_mid = lspline.curve
            self.SAC_mid.plotcurve()
            
            ##
            ## aft section:
            ##
            length = segment_span['aft'][1] - segment_span['aft'][0]
            
            xm_factor = .3
            ym_factor = .2
            xb = segment_span['aft'][0]
            yb = self.maxSecArea
            alphab = 0.
            Cab_given = 0.
            xe = segment_span['aft'][1]
            ye = 0.
            alphae = 0.
            Cae_given = 0.
            
            height = yb-ye
            xc = xm_factor*length + segment_span['mid'][1] - segment_span['mid'][0] + segment_span['fwd'][1] - segment_span['fwd'][0]
            yc = ym_factor*height
            area = areas[2]*self.BareHullVolume #####################EEEEERRORRR
            
            #spec = {'SAC':{'maxwidth':self.maxSecArea}}
            spec['xmin'] = xb
            spec['xmax'] = xe
            spec['ymin'] = ye
            spec['ymax'] = yb
            #vertices = InitializeControlPoints(xb,yb,alphab,Cab_given,xe,ye,alphae,Cae_given,xc,yc,area,spec)
            temp_object3 = InitializeControlVertices(xb,yb,xe,ye,
                                                     alphab,alphae,
                                                     Cab_given,Cae_given,
                                                     area,xc,yc,7,'down',spec)
            vertices = temp_object3.vertices
            self.aftHarries = Bspline(vertices,k,nump)
            #vertices = linear_vertices((xb,yb),(xe,ye),7)
            #self.aftLinear = Bspline(vertices,k,nump)
            self.SAC_aft = Bspline(vertices,k,nump)
            self.SAC_aft.plotcurve()

            
            FP = FormParameter(self.SAC_aft)
            IC = InequalityConstraint(self.SAC_aft)
            
            FP.add_AngleConstraint(0.,0.) 
            #FP.add_VerticalConstraint(0.,0.)
            FP.add_AngleConstraint(0.,1.)
            
            #FP.add_CurvatureConstraint(Cab_given,0.)
            #FP.add_CurvatureConstraint(Cae_given,1.)
            #area = self.SAC_aft.area.value
            FP.add_AreaConstraint(area)
            #FP.add_XcConstraint(xc)
            
            self.SAC_aft.FormParam = FP.FormParam
            self.SAC_aft.IneqConstraints = IC.InequalityParam
            
    
            lspline = LagrangeSpline(self.SAC_aft, 
                                                 {'Inumber':25,
                                                   'ws':10.,
                                                   'w1':.1,
                                                   'w2':.1,
                                                   'w3':.1}, None, None, None)
                                   
            self.SAC_aft =lspline.curve
            
            self.SAC_aft.plotcurve()
            
            SAC_vertices = np.zeros((self.SAC_fwd.n+self.SAC_mid.n+self.SAC_aft.n,2),float)
            SAC_vertices[:self.SAC_fwd.n] = self.SAC_fwd.vertices
            SAC_vertices[self.SAC_fwd.n:self.SAC_fwd.n+self.SAC_mid.n] = self.SAC_mid.vertices
            SAC_vertices[self.SAC_fwd.n+self.SAC_mid.n:] = self.SAC_aft.vertices
            self.SAC = Bspline(SAC_vertices, k, nump)
            self.SAC.segment_span = segment_span
            self.f = lspline.f
        self.SAC.plotcurve()
        #self.SAC = parametric_curve(nondim_factor)
        #self.get_SAC_from_TSK()
        return
        
    def get_SAC_from_TSK(self):
        """
            TODO:
                implement a 2d ADLopt 
                that incoprporates the sqrt normed distance
                as a objective
                just like in 3d opt.
                Then use that here to create a b-spline curve
                which interpolates the SAC Curve.
                
                -or do that in the SAC stuff itself..?
        """
        self.SACvertices = get_TSK_SAC()
        num = len(self.SACvertices)
        dumbarray = np.asarray([self.SACvertices[0],self.SACvertices[num/5],self.SACvertices[2*num/5],self.SACvertices[3*num/5],self.SACvertices[4*num/5],self.SACvertices[-1]])
        self.curve = interpolatedBspline(dumbarray, 4, 50)
        plt.plot(self.SACvertices[:,0],self.SACvertices[:,1])
        plt.plot(self.curve.r[:,0],self.curve.r[:,1])
        plt.show()
        
        
        #now get sac from opt:
        points = []
        s_vector = []
        for i in range(0,len(self.SACvertices),20):
            points.append(self.SACvertices[i])
            s_vector.append(naive_minimizer(self.curve, points[-1]))
        
        #print s_vector
        #print self.curve.vertices
        #svector = []
        FP=FormParameter(self.curve)
        self.curve.FormParam=FP.FormParam
        self.curve,f = ADLsplineOpt(self.curve, 
                               {'Inumber':15,
                               'ws':10.,
                               'w1':.1,
                               'w2':.1,
                               'w3':.1}, 5000., points, s_vector)
        plt.plot(self.SACvertices[:,0],self.SACvertices[:,1])
        plt.plot(self.curve.r[:,0],self.curve.r[:,1])
        plt.show()
        
        return
        
    def defineDWL(self):
        k               = 4
        nump            = 50
        nondim_factor   = self.LengthDWL

        s1=.4
        s2=.7  ##---what about SAC and the flat of Side - these don't have to match that- but they might want too sometimes.
        self.s1=s1
        self.s2=s2
        #Aggregate DwL parameters:        
        area_factor     = 2./10.
        bulk_factor     = 0.333
        
        centroid_factor = 0.5        
        start           = [0.,0.]
        end             = [self.LengthDWL, 0.]
        max_area        = area_factor*nondim_factor
        Xc              = centroid_factor*nondim_factor
        
        # In Form Parameter terminology:
        xb = start[0]
        yb = start[1]
        alphab = 45.
        Cab_given = 0.
        xe = end[0]
        ye = end[1]
        alphae = 45.
        Cae_given = 0.
        xc = Xc
        yc = max_area*bulk_factor
        area = self.WLarea
        
        #DWL specific initializier:
        spec = {'DWL':{'maxwidth':max_area}}
        spec['xmin'] = xb
        spec['xmax'] = xe
        spec['ymin'] = yb
        spec['ymax'] = self.maxSecArea
        vertices = InitializeControlPoints(xb,yb,alphab,Cab_given,xe,ye,alphae,Cae_given,xc,yc,area,spec)
        #temp_object = InitializeControlVertices(xb,yb,xe,ye,alphab,alphae,Cab_given,Cae_given,area,xc,yc,7,spec)
        #vertices = temp_object.vertices
        #vertices = linear_vertices((xb,yb),(xe,yb),7)        
        self.DWL = Bspline(vertices,k,nump)
        self.DWL.plotcurve()
        
        self.DWL.basis_matrix()
        self.DWL.MomentMatrices()
        self.DWL.curvearea()
        
        
        
        print 'Setting up DWL form parameters'
        FP = FormParameter(self.DWL)
        IC = InequalityConstraint(self.DWL)
        
        #s_b = self.DWL.CurvePoint(self.SAC.segment_span['fwd'][1])[0]
        FP.add_AngleConstraint(0.,s1) 
        FP.add_AngleConstraint(0.,s2)
        #FP.add_VerticalConstraint(0.,0.)
        #FP.add_CurvatureConstraint(0.,s1)
        #FP.add_CurvatureConstraint(0.,s2)
        FP.add_AreaConstraint(area)
        #FP.add_XcConstraint(xc)
        
        self.DWL.FormParam = FP.FormParam
        self.DWL.IneqConstraints = IC.InequalityParam
        

        lspline = LagrangeSpline(self.DWL, 
                                             {'Inumber':25,
                                               'ws':10.,
                                               'w1':.1,
                                               'w2':.1,
                                               'w3':.1}, None, None, None)
                                               
        FP = FormParameter(self.DWL)
        IC = InequalityConstraint(self.DWL)
        
        #s_b = self.DWL.CurvePoint(self.SAC.segment_span['fwd'][1])[0]
        FP.add_AngleConstraint(0.,s1) 
        FP.add_AngleConstraint(0.,s2)
        #FP.add_VerticalConstraint(0.,0.)
        FP.add_CurvatureConstraint(0.,s1)
        FP.add_CurvatureConstraint(0.,s2)
        FP.add_AreaConstraint(area)
        #FP.add_XcConstraint(xc)
        
        self.DWL.FormParam = FP.FormParam
        self.DWL.IneqConstraints = IC.InequalityParam
        

        lspline = LagrangeSpline(self.DWL, 
                                             {'Inumber':25,
                                               'ws':10.,
                                               'w1':.1,
                                               'w2':.1,
                                               'w3':.1}, None, None, None)
                               
        self.DWL =lspline.curve
        
        print 'DWL WaterPlane area = {}'.format(self.DWL.area.value)
        return
    
    def augDWL(self):
        s1=self.s1
        s2=self.s2
        FP = FormParameter(self.DWL)
        IC = InequalityConstraint(self.DWL)
        
        #s_b = self.DWL.CurvePoint(self.SAC.segment_span['fwd'][1])[0]
        FP.add_AngleConstraint(0.,s1) 
        FP.add_AngleConstraint(0.,s2)
        #FP.add_VerticalConstraint(0.,0.)
        FP.add_CurvatureConstraint(0.,s1)
        FP.add_CurvatureConstraint(0.,s2)
        #FP.add_AreaConstraint(area)
        #FP.add_XcConstraint(xc)
        FP.add_yPointConstraint(self.MaxBreadth,s1)        
        FP.add_yPointConstraint(self.MaxBreadth,s2) 
        FP.add_xPointConstraint(self.stations['flat_DWL'][0],s1)
        FP.add_xPointConstraint(self.stations['flat_DWL'][1],s2)
        
        self.DWL.FormParam = FP.FormParam
        self.DWL.IneqConstraints = IC.InequalityParam
        

        lspline = LagrangeSpline(self.DWL, 
                                             {'Inumber':25,
                                               'ws':10.,
                                               'w1':.1,
                                               'w2':.1,
                                               'w3':.1}, None, None, None)
                               
        self.DWL =lspline.curve
        return

#    def buildbulbousbow(self, halfbreadth=.5, length=1., depth=1., volume=np.pi*(.5**2), num_stations=4):
#        """
#            Kracht: Bulbous bows
#               - 3 main bulb types (center of area determined):
#                   1. delta
#                   2. O
#                   3. nabla
#               - 6 (more) parameters are sufficient for design
#                   linear:
#                       1. breadth
#                       2. length
#                       3. depth
#                   nonlinear:
#                       4. cross sectional area
#                       5. lateral projected area (profile)
#                       6. volume
#                       0. vertical center of area of each cross section (already mentioned)
#                       
#                - for cnst vol & depth, interference depends on bulb length / hull length
#                - wave breaking depends on
#        """
#        num_stations                = num_stations
#        self.bulbousbow = None
#        self.bulbousbow.halfbreadth       = halfbreadth
#        self.bulbousbow.length            = length
#        self.bulbousbow.depth             = depth
#        self.bulbousbow.SAC               = None
#        self.bulbousbow.vSAC              = None
#        self.bulbousbow.profile           = None
#        self.bulbousbow.volume            = volume
#        
#        xb=0.
#        yb=0.5*np.pi*( (self.bbow.depth**2)*.25 + (self.bbow.halfbreadth**2)*.5 )
#        alphab=0.
#        Cab_given=0.
#        xe=self.bbow.length
#        ye = 0.
#        alphae = 90.
#        Cae_given = 0.
#        xc = self.bbow.length*0.5
#        yc = yb*0.5
#        area = self.bbow.volume
#        bb_max_sec_area =   InitializeControlPoints(xb,yb,alphab,Cab_given,xe,ye,alphae,Cae_given,xc,yc,area)   
#        #bb_SAC_vertices = linear_vertices(start, end, num)
#        return
        
        
#    def buildtransomstern(self):
#        return
        
    def buildCenterPlane(self):
        self.centershpcontrol ={'deck_aft':[(self.LengthDWL-20.,self.Depth),    (self.LengthDWL,self.Depth)     ],
                                'deck_mid':[(self.LengthDWL-80.,self.Depth),    (self.LengthDWL-20.,self.Depth) ],
                                'deck_fwd':[(0.,self.Depth),                    (self.LengthDWL-80.,self.Depth) ],
                                'baseline_aft':[(self.LengthDWL-20.,0.),        (self.LengthDWL,4.)             ],
                                'baseline_mid':[(self.LengthDWL-80.,0.),        (self.LengthDWL-20.,0.)         ],
                                'baseline_fwd':[(0.,0.),                        (self.LengthDWL-80.,0.)         ]}

        base_fwd = linear_vertices(self.centershpcontrol['baseline_fwd'][0], self.centershpcontrol['baseline_fwd'][1], self.n)
        base_mid = linear_vertices(self.centershpcontrol['baseline_mid'][0], self.centershpcontrol['baseline_mid'][1], self.n)
        base_aft = linear_vertices(self.centershpcontrol['baseline_aft'][0], self.centershpcontrol['baseline_aft'][1], self.n)
        
        deck_fwd =  linear_vertices(self.centershpcontrol['deck_fwd'][0], self.centershpcontrol['deck_fwd'][1], self.n)  
        deck_mid =  linear_vertices(self.centershpcontrol['deck_mid'][0], self.centershpcontrol['deck_mid'][1], self.n)  
        deck_aft =  linear_vertices(self.centershpcontrol['deck_aft'][0], self.centershpcontrol['deck_aft'][1], self.n)  
        
        self.centerplane = []
        
        self.centerplane_fwd_base = Bspline(base_fwd, self.k, self.nump)
        self.centerplane_mid_base = Bspline(base_mid, self.k, self.nump)
        self.centerplane_aft_base = Bspline(base_aft, self.k, self.nump)
        
        self.centerplane_fwd_deck = Bspline(deck_fwd, self.k, self.nump)
        self.centerplane_mid_deck = Bspline(deck_mid, self.k, self.nump)
        self.centerplane_aft_deck = Bspline(deck_aft, self.k, self.nump)
        
        self.centerplane_base = [self.centerplane_fwd_base, self.centerplane_mid_base, self.centerplane_aft_base]
        self.centerplane_deck = [self.centerplane_fwd_deck, self.centerplane_mid_deck, self.centerplane_aft_deck]
        
        return
        
    def MainDeck(self):
        return
        
    def buildtransversecurves(self, num):
        
        n               = 4
        nump            = 50
        nondim_factor   = self.LengthDWL

        # Idefaults:
        xb = 0.#start[0]
        yb = 0.#start[1]
        alphab = 0.
        Cab_given = 0.
        xe = None #end[0]
        ye = None #end[1]
        alphae = 45.
        Cae_given = 0.
        xc = None#Xc
        yc = None#max_area*bulk_factor
        area = None#self.BareHullVolume
        
        #transverse specific initializier:
        self.tCurveNet = []
        
        #self.trange = np.linspace(20.,self.LengthDWL-20.,num)
        fltst3 = self.stations['flat_side'][0] -2.
        fltst2 = self.stations['flat_side'][0] -1.
        fltst1 = self.stations['flat_side'][0] 
        
        fltend3 = self.stations['flat_side'][1]-2.
        fltend2 = self.stations['flat_side'][1]-1.
        fltend1 = self.stations['flat_side'][1]
        
        self.trange = np.asarray([self.stations['bow_shape'][0], fltst3, fltst2, fltst1, fltend3, fltend2, fltend1  ])
        LocalFormCurve = ['None',self.fos3,self.fos2,self.fos1]
        for position in self.trange:
            
            s_SAC  = self.SAC.FindPoint(position, 0)[0]#position/self.LengthDWL
            s_dwl  = self.DWL.FindPoint(position, 0)[0]
            #s_dwl   = position/self.
            
            
            xe = self.DWL.CurvePoint(s_dwl)[1]#self.DWL.CurvePoint(s_SAC)[1]
            ye = self.MaxDraft
            
            area    = self.SAC.CurvePoint(s_SAC)[1]#abs((xe-xb)*(ye-yb)-self.SAC.CurvePoint(s_SAC)[1])
            print 'tartget area = {}'.format(area)
            xc = (xe-xb)/2.
            yc = (ye-yb)/2.
            spec = {'name':'transverse'}
            spec['xmin'] = xb
            spec['xmax'] = xe
            spec['ymin'] = yb
            spec['ymax'] = ye
            vertices = InitializeControlPoints(xb,yb,alphab,Cab_given,xe,ye,alphae,Cae_given,xc,yc,area)
            temp_object = InitializeControlVertices(xb,yb,xe,ye,
                                                    alphab,alphae,
                                                    Cab_given,Cae_given,
                                                    area,xc,yc,7,'up')
            vertices = temp_object.vertices
            #vertices = linear_vertices((xb,yb),(xe,ye),7)
            curve = Bspline(vertices, n, nump)
            
 
            FP=FormParameter(curve)
            IC=InequalityConstraint(curve)
            #if abs(area - self.maxSecArea) < self.tol:
                #FP.add_CurvatureConstraint(0.,0.)
                #FP.add_AngleConstraint(0.,0.)
                #IC.set_max_tan_value(25.,0.)
            #else:
            #FP.add_AreaConstraint(area)
            FP.add_AreaConstraint_y_axis(area)
            curve.FormParam=FP.FormParam
            curve.IneqConstraints = IC.InequalityParam
            lspline = LagrangeSpline(curve, 
                                      {'Inumber':20,
                                       'ws':.0001,
                                       'w1':.0001,
                                       'w2':.0001,
                                       'w3':.0001}, 
                                       None,None,None)
            self.tCurveNet.append(lspline.curve)
        
    
        return
        
#    def buildbow(self):
#        # Idefaults:
#        xb = 0.#start[0]
#        yb = 0.#start[1]
#        alphab = 0.
#        Cab_given = 0.
#        xe = None #end[0]
#        ye = None #end[1]
#        alphae = 45.
#        Cae_given = 0.
#        xc = None#Xc
#        yc = None#max_area*bulk_factor
#        area = None#self.BareHullVolume
#        xe = self.DWL.CurvePoint(s_dwl)[1]#self.DWL.CurvePoint(s_SAC)[1]
#        ye = self.MaxDraft
#        
#        area    = abs((xe-xb)*(ye-yb)-self.SAC.CurvePoint(s_SAC)[1])
#        print 'tartget area = {}'.format(area)
#        xc = (xe-xb)/2.
#        yc = (ye-yb)/2.
#        vertices = InitializeControlPoints(xb,yb,alphab,Cab_given,xe,ye,alphae,Cae_given,xc,yc,area)
#        curve = Bspline(vertices, n, nump)
#        return
        
    
        
    def build_3D_transverse_curves(self):
        self.tCurveNet3D = []
        self.tnet = np.zeros((len(self.tCurveNet),7,3,),float)
        for i, curve in enumerate(self.tCurveNet):
            self.tnet[i,:,0:2]   = curve.vertices[:]
            self.tnet[i,:,2]     = self.trange[i]
            self.tCurveNet3D.append(Bspline(self.tnet[i],self.k,self.nump))
            
        
            
        for bowcurve in self.bowcurves:
            self.tCurveNet3D.insert(0,bowcurve)
        for bulbcurve in reversed(self.bulbcurves):
            self.tCurveNet3D.insert(0,bulbcurve)
        
        offset = len(self.bowcurves)+len(self.bulbcurves)
        #self.store = self.tCurveNet3D.pop(offset+1)
        """
        self.tCurveNet3D.pop(offset)
        self.tCurveNet3D.pop(offset)
        self.tCurveNet3D.pop(offset)
        self.tCurveNet3D.pop(offset)
        self.tCurveNet3D.insert(offset,self.fos1)
        self.tCurveNet3D.insert(offset,self.fos2)
        self.tCurveNet3D.insert(offset,self.fos3)
        self.tCurveNet3D.insert(offset+3,self.midshipsec)
        #"""
        #self.tCurveNet3D.pop(offset)  ## this one requires knot insertion on Lspline fitting on the lCurveNet
        try:
            self.tCurveNet3D.pop(-1)
            self.tCurveNet3D.insert(len(self.tCurveNet3D),self.sternsec)
            print 'added transom stern'
        except:
            print 'standard stern'
        self.tnet = np.zeros((len(self.tCurveNet3D),7,3),float)  # TLM deleted a comma after the '3' Aug7_2014
        for i, curve in enumerate(self.tCurveNet3D):
            self.tnet[i] = self.tCurveNet3D[i].vertices
        return
    def buildlongitudinalcurves(self):
        k               = 4
        nump            = 50
        nondim_factor   = self.LengthDWL
        
        self.build_3D_transverse_curves()
        
        ## --- --- --- Make Longitudinal curve Net --- --- --- ##
        ##
        lnet    = []
        #tray = np.asarray(tnet)
        for i in range(len(self.tnet[0])):
            lnet.append(self.tnet[:,i])
            
        self.lCurveNet = []
        for ii,vertices in enumerate(lnet):
            tempCurve = interpolatedBspline(vertices, k, nump)
    
            self.lCurveNet.append(tempCurve)
        ##
        ##  TODO:
            ## append a curveback into tCurveNet3d
            ## OR
            ## insert a knot in each curve
            ## then do a curve Lspline fit into the original curves
            ## fixing the tcurve points and maybe other properties
        
        return
        
    def build_primary_hull_surface(self):
        self.main_hull = BsplineSurface(self.tCurveNet3D,self.lCurveNet)
        return
 
 
 
#class bulbousbow_data(object):
#    """
#        Observer Pattern, for an agent:
#        
#        If ever the bulbousbow_data changes,
#        then update the bulbousbow
#    """       
#    def __init__(self, kind, breadth, length, depth, max_sec_area, sec_cx, sec_cy, max_lateral_projection_area, lat_cx, lat_cy, volume):
#        self.kind = kind #delta, O, nabla
#        self.breadth = breadth
#        self.length = length
#        self.depth = depth
#        self.max_sec_area = max_sec_area
#        self.sec_cx
#        self.sec_cy
#        self.max_lateral_projection_area
#        self.lat_cx
#        self.lat_cy
#        self.volume = volume
#        
#        return
    

class bulbousbow(Hull):
    """
            Kracht: Bulbous bows
               - 3 main bulb types (center of area determined):
                   1. delta
                   2. O
                   3. nabla
               - 6 (more) parameters are sufficient for design
                   linear:
                       1. breadth
                       2. length
                       3. depth
                   nonlinear:
                       4. cross sectional area
                       5. lateral projected area (profile)
                       6. volume
                       0. vertical center of area of each cross section (already mentioned)
                       
                - for cnst vol & depth, interference depends on bulb length / hull length
                - wave breaking depends on
    """
    def __init__(self, length=1., halfbreadth=.5, depth=1., volume=np.pi*.5, num_stations=4):
        self.k              = 4
        self.n              = 7
        self.nump           = 50
        self.SAC            = None
        self.vSAC           = None
        self.profile        = None
        self.num_stations   = num_stations
        #self.precept(boundary_curve)
        self.principle_particulars(length, halfbreadth, depth, volume)
        self.derived_particulars()
        self.block_coefficient = self.volume/(self.length * self.halfbreadth * self.depth)
        self.basic_curves()
        #self.derived_curves()
        self.buildtransversecurves()
        self.definebowcurve()
        self.buildlongitudinalcurves()
        self.build_primary_hull_surface()
        
        
    #def precept(self, boundary_curve):
    #    self.boundary_curve = boundary_curve
                
        return
    def principle_particulars(self, length, halfbreadth, Depth, volume):
        self.length         = length
        self.LengthDWL      = length
        self.halfbreadth    = halfbreadth
        self.depth          = Depth
        self.volume         = volume
        self.BareHullVolume = volume
        return
    
    def derived_particulars(self):
        self.maxSecArea = 0.4*np.pi*( (self.depth**2)*.25 + (self.halfbreadth**2)*.5 )
        return
        
    def basic_curves(self):
        self.defineSAC()
        return
        
    def derived_curves(self):
        pass
        return
    
    def defineSAC(self):
        xb=0.
        yb=0.5*np.pi*( (self.depth**2)*.25 + (self.halfbreadth**2)*.5 )
        alphab=0.
        Cab_given=0.
        xe=self.length
        ye = 0.
        alphae = 90.
        Cae_given = 0.
        xc = self.length*0.5
        yc = yb*0.5
        area = self.volume
        #bb_max_sec_area =   InitializeControlPoints(xb,yb,alphab,Cab_given,xe,ye,alphae,Cae_given,xc,yc,area)
        bb_max_sec_area = InitializeControlVertices(xb,yb,xe,ye,
                                                    alphab,alphae,
                                                    Cab_given,Cae_given,
                                                    area,xc,yc,7,'up')
        bb_max_sec_area = bb_max_sec_area.vertices
        
        
        #bb_SAC_vertices = linear_vertices(start, end, num)   
        
        num_segments = 2
        self.SACrange = np.linspace(0.,self.length,num_segments+1)
        
        segment_span = {'fwd':self.SACrange[:2],'mid':self.SACrange[1:3]}
        #SAC_segments = []
        
        ## fwd segment:
        length = segment_span['fwd'][1] - segment_span['fwd'][0]
        
        
        xm_factor = .7
        ym_factor = .2
        xb = segment_span['fwd'][0]
        yb = 0.
        alphab = 80.
        Cab_given = 0.
        xe = segment_span['fwd'][1]
        ye = self.maxSecArea
        alphae = 0.
        Cae_given = 0.
        
        height = ye-yb
        xc = xm_factor*length
        yc = ym_factor*height
        area = 0.4*self.volume #####################EEEEERRORRR
        
        spec = {'SAC':{'maxwidth':self.maxSecArea}}
        spec = {'name':'SAC'}
        spec['xmin'] = xb
        spec['xmax'] = xe
        spec['ymin'] = yb
        spec['ymax'] = ye
        #self.SAC_fwd = InitializeControlPoints(xb,yb,alphab,Cab_given,xe,ye,alphae,Cae_given,xc,yc,area,spec)
        self.SAC_fwd = InitializeControlVertices(xb,yb,xe,ye,
                                                 alphab,alphae,
                                                 Cab_given,Cae_given,
                                                 area,xc,yc,7,'down',spec)
        self.SAC_fwd = self.SAC_fwd.vertices
        self.SAC_fwd = Bspline(self.SAC_fwd,self.k,self.nump)
        #"""
        FP = FormParameter(self.SAC_fwd)
        IC = InequalityConstraint(self.SAC_fwd)
        
        
        self.SAC_fwd.FormParam = FP.FormParam
        self.SAC_fwd.IneqConstraints = IC.InequalityParam
        
        
        lspline = LagrangeSpline(self.SAC_fwd, 
                                             {'Inumber':25,
                                               'ws':10.,
                                               'w1':.1,
                                               'w2':.1,
                                               'w3':.1}, None, None, None)
                               
        self.SAC_fwd = lspline.curve
        #"""
        self.SAC_fwd.basis_matrix()
        self.SAC_fwd.MomentMatrices()
        self.SAC_fwd.curvearea()
        area = self.SAC_fwd.area
        FP = FormParameter(self.SAC_fwd)
        IC = InequalityConstraint(self.SAC_fwd)
        
        #FP.add_AngleConstraint(0.,0.) 
        #FP.add_VerticalConstraint(0.,0.)
        #FP.add_AngleConstraint(0.,1.)
        
        #FP.add_CurvatureConstraint(Cab_given,0.)
        #FP.add_CurvatureConstraint(Cae_given,1.)
        FP.add_AreaConstraint(area)
        #FP.add_XcConstraint(xc)
        
        self.SAC_fwd.FormParam = FP.FormParam
        self.SAC_fwd.IneqConstraints = IC.InequalityParam
        

        lspline = LagrangeSpline(self.SAC_fwd, 
                                             {'Inumber':25,
                                               'ws':10.,
                                               'w1':.1,
                                               'w2':.1,
                                               'w3':.1}, None, None, None)
                               
        self.SAC_fwd =lspline.curve
        
        ## mid-body SAC:
        length = segment_span['mid'][1] - segment_span['mid'][0]
        
        xm_factor = .5
        ym_factor = .5
        xb = segment_span['mid'][0]
        yb = self.maxSecArea
        alphab = 0.
        Cab_given = 0.
        xe = segment_span['mid'][1]
        ye = self.maxSecArea
        alphae = 0.
        Cae_given = 0.
        
        height = yb
        xc = xm_factor*length
        yc = ym_factor*height
        area = 0.6*self.volume
        
        spec = {'SAC':{'maxwidth':self.maxSecArea}}
        spec = {'name':'SAC'}
        spec['xmin'] = xb
        spec['xmax'] = xe
        spec['ymin'] = yb
        spec['ymax'] = ye
        self.SAC_mid = linear_vertices((xb,yb),(xe,ye),7)
        self.SAC_mid = Bspline(self.SAC_mid,self.k,self.nump)
        """
        FP = FormParameter(self.SAC_mid)
        IC = InequalityConstraint(self.SAC_mid)
        

        
        self.SAC_mid.FormParam = FP.FormParam
        self.SAC_mid.IneqConstraints = IC.InequalityParam
        lspline = LagrangeSpline(self.SAC_mid, 
                                             {'Inumber':25,
                                               'ws':10.,
                                               'w1':.1,
                                               'w2':.1,
                                               'w3':.1}, None, None, None)
                               
        self.SAC_mid = lspline.curve
        self.f = lspline.f
        #"""
        
        
        
        SAC_vertices = np.zeros((self.SAC_fwd.n+self.SAC_mid.n,2),float)
        SAC_vertices[:self.SAC_fwd.n] = self.SAC_fwd.vertices
        SAC_vertices[self.SAC_fwd.n:] = self.SAC_mid.vertices
        self.SAC = Bspline(SAC_vertices, self.k, self.nump)
        self.SAC.segment_span = segment_span
        
        return  
        
    def buildtransversecurves(self):
        num = self.num_stations
        self.tCurveNet = []
        xb = 0.
        yb = 0.
        xe = 0.
        ye=self.depth
        
        vertices = arc_vertices((xb,yb),(xe,ye),7)
        self.trange = np.linspace(0.,0.8*self.length,num)
        
        for position in self.trange:
            
            s_SAC  = self.SAC.FindPoint(position, 0)[0]
        
            area = self.SAC.CurvePoint(s_SAC)[1]
    
            dummycurve = Bspline(vertices,self.k,self.nump)
            self.SAC_fwd.basis_matrix()
            self.SAC_fwd.MomentMatrices()
            self.SAC_fwd.curvearea()
            FP = FormParameter(dummycurve)
            IC = InequalityConstraint(dummycurve)
            
            #FP.add_AngleConstraint(0.,0.) 
            #FP.add_VerticalConstraint(0.,0.)
            #FP.add_AngleConstraint(0.,1.)
            #FP.add_xPointConstraint()
            #FP.add_CurvatureConstraint(Cab_given,0.)
            #FP.add_CurvatureConstraint(Cae_given,1.)
            FP.add_AreaConstraint_y_axis(area)
            #FP.add
            #FP.add_XcConstraint(xc)
            #IC.set_min_vertex_x_value(0.,1)
            IC.set_min_vertex_x_value(0.,2)
            #IC.set_min_vertex_x_value(0.,3)
            #IC.set_min_vertex_x_value(0.,4)
            IC.set_min_vertex_x_value(0.,5)
            #IC.set_min_vertex_x_value(0.,6)
            dummycurve.FormParam = FP.FormParam
            dummycurve.IneqConstraints = IC.InequalityParam
            
    
            lspline = LagrangeSpline(dummycurve, 
                                     {'Inumber':25,
                                       'ws':10.,
                                       'w1':.1,
                                       'w2':.1,
                                       'w3':.1}, None, None, None)
                                   
            dummycurve =lspline.curve
            self.tCurveNet.append(dummycurve) 
            
        self.definebowcurve()
        return

    def build_3D_transverse_curves(self):
        #Bulbous Bow "Transverse" 3D curves
        self.tCurveNet3D = []
        self.tnet = np.zeros((len(self.tCurveNet),7,3,),float)
        for i, curve in enumerate(self.tCurveNet):
            self.tnet[i,:,0:2]   = curve.vertices[:]
            self.tnet[i,:,2]     = self.trange[i]
            self.tnet[i] = transpose_vertices(self.tnet[i], 2, -self.length)
            self.tCurveNet3D.append(Bspline(self.tnet[i],self.k,self.nump))
            
            
        self.tCurveNet3D.insert(0,self.bowcurve)
        
        self.tnet = np.zeros((len(self.tCurveNet3D),7,3,),float)
        for i, curve in enumerate(self.tCurveNet3D):
            self.tnet[i] = self.tCurveNet3D[i].vertices
        return
    
        
    def buildlongitudinalcurves(self):
        k               = self.k
        nump            = self.nump
        nondim_factor   = self.length
        
        self.build_3D_transverse_curves()
        
        ## --- --- --- Make Longitudinal curve Net --- --- --- ##
        ##
        lnet    = []
        #tray = np.asarray(tnet)
        for i in range(len(self.tnet[0])):
            lnet.append(self.tnet[:,i])
            
        self.lCurveNet = []
        for ii,vertices in enumerate(lnet):
            tempCurve = interpolatedBspline(vertices, k, nump)
    
            self.lCurveNet.append(tempCurve)
        ##
        ##  TODO:
            ## append a curveback into tCurveNet3d
            ## OR
            ## insert a knot in each curve
            ## then do a curve Lspline fit into the original curves
            ## fixing the tcurve points and maybe other properties
        
        return
        
    def definebowcurve(self):
        num = self.num_stations
        xb = 0.
        yb = 0.
        xe = 0.
        ye = self.depth
        
        #length = 0.2*self.length
        #area = np.pi*
        
        vertices = linear_vertices((xb,yb),(xe,ye),7)
        self.trange = np.linspace(0.,0.8*self.LengthDWL,num)
                
        s_SAC  = self.SAC.FindPoint(self.trange[2], 0)[0]   ### not general!!!
        area = 0.5*self.SAC.CurvePoint(s_SAC)[1]  ### fixed artifice
        
        dummycurve = Bspline(vertices,self.k,self.nump)
        dummycurve.basis_matrix()
        dummycurve.MomentMatrices()
        dummycurve.curvearea()
        FP = FormParameter(dummycurve)
        IC = InequalityConstraint(dummycurve)
        FP.add_AreaConstraint(area)
        
        dummycurve.FormParam = FP.FormParam
        dummycurve.IneqConstraints = IC.InequalityParam
        

        lspline = LagrangeSpline(dummycurve, 
                                 {'Inumber':25,
                                   'ws':10.,
                                   'w1':.1,
                                   'w2':.1,
                                   'w3':.1}, None, None, None)
                               
        dummycurve  =lspline.curve
        bnet        = np.zeros((7,3),float)
        bnet[:,1]  = dummycurve.vertices[:,1]
        bnet[:,2]  = dummycurve.vertices[:,0]
        
        bnet = transpose_vertices(bnet, 2, -self.length-2.)
        self.bowcurve = Bspline(bnet,self.k,self.nump)
        #self.tCurveNet.insert(0,dummycurve) 
        
        return
        


class Vessel(Hull): 
    
    def __init__(self, LDWL, MaxBreadth, MaxDraft, Depth, BareHullVolume):
        self.principle_particulars(LDWL, MaxBreadth, MaxDraft, Depth, BareHullVolume)
        self.block_coefficient = self.BareHullVolume/(self.LDWL * self.MaxBreadth * self.MaxDraft)
   
class OSV(Hull):
    def __init__(self, LDWL, MaxBreadth, MaxDraft, Depth, BareHullVolume):
        self.principle_particulars(LDWL, MaxBreadth, MaxDraft, Depth, BareHullVolume)
        self.block_coefficient = self.BareHullVolume/(self.LDWL * self.MaxBreadth * self.MaxDraft)
    
    
if __name__ =="__main__":
    #"""
    LDWL            = 100.
    MaxBreadth      = 10.
    MaxDraft        = 15.
    Depth           = 15.
    BareHullVolume  = 10000.
    self = Hull(LDWL, MaxBreadth, MaxDraft, Depth, BareHullVolume)
    #"""
    #"""
    scale = 0.4
    bblength = scale*15.
    bbhalfbreadth=scale*7.5 
    bbdepth=Depth
    bbvolume=scale*bbdepth*bblength*bbhalfbreadth
    bbnum_stations=4 
    bbinterface = {}#{'hull_boundar_curve':self.tCurveNet[0]}
    bbow = bulbousbow(bblength,bbhalfbreadth,bbdepth,bbvolume,bbnum_stations)
    #self.main_hull.plotMultiSurface([self.main_hull],LDWL,LDWL,LDWL)
    #"""
    s0 = mesh(self.main_hull.surface[:,:,0],self.main_hull.surface[:,:,1],self.main_hull.surface[:,:,2], colormap='YlGnBu')
    s1 = mesh(-self.main_hull.surface[:,:,0],self.main_hull.surface[:,:,1],self.main_hull.surface[:,:,2], colormap='YlGnBu')
    show()
    self.fos1.plot3DmultiList(self.tCurveNet3D,[self.fos1,self.fos2,self.fos3])
    #"""
    self.add_bulbous_bow(bbow)
    s0 = mesh(self.main_hull.surface[:,:,0],self.main_hull.surface[:,:,1],self.main_hull.surface[:,:,2], colormap='YlGnBu')
    s1 = mesh(-self.main_hull.surface[:,:,0],self.main_hull.surface[:,:,1],self.main_hull.surface[:,:,2], colormap='YlGnBu')
    show()

    self.sub_bulbous_bow()
    s0 = mesh(self.main_hull.surface[:,:,0],self.main_hull.surface[:,:,1],self.main_hull.surface[:,:,2], colormap='YlGnBu')
    s1 = mesh(-self.main_hull.surface[:,:,0],self.main_hull.surface[:,:,1],self.main_hull.surface[:,:,2], colormap='YlGnBu')
    #"""
