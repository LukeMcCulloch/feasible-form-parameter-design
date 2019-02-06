# 20082013
# Python module to compute the Frenet Frame
# of B-Spline curves
# with automatic differentiation

"""
CURRENTLY THIS IS OF HISTORICAL USE ONLY!

TODO: integrate this into the current framework
"""

import numpy as np

import copy
import time


from AD import curveAD, surfAD, surf_curveAD    # Automatic differentiation of B-spline curves & surfaces
import AD
from routines import ADcrossproduct   # LB cross product routine

from curve import *                 # bspline class and methods
from Equalityconstraints  import FormParameter 
from InequalityConstraints import InequalityConstraint
import SurfaceFormParameters as sfp

from initialValues import InitializeControlPoints, InitializeControlVertices

import ADL_LS_OptBspline as ADLS

#from ADL_LS_OptBspline import LagrangeSpline, LagrangeSurface
#from ADLOptBspline3D import ADLsplineOpt

#from mayavi.mlab import *


## singlton pattern:
import DimensionTools

#LB - more good ideas:
MINFLOAT = 1.E-12
nonZero = .00000001
lb3dFLOAT = float


def ADnormalize(vector):
    """
       Normalize ADvectors to unit length
    """
    vector = np.asarray(vector)
    
    return vector/np.sqrt(sum(vector**2)) # this will preserve all AD operations!
    
    

def ADtanVector(curve,s):
    """ length preserving tangent
        AD form

        inputs:
            curve
            point

        outputs:
            tangent vector"""
    
    p = curve.p
    dim = curve.dim
    localBasis = np.zeros((curve.n,curve.n),float)
    span = curve.FindSpan(s)
    curve.DersBasisFunc(span,s,localBasis[:,span-p:span+1])

    ADvertices = AD.curveAD(curve)

    tg=[]
    for i in range(dim):
        tg.append(np.dot(ADvertices[i],localBasis[1]))
    
    return tg


def ADbinormalVector(curve, s):
    """ length preserving tangent
        AD form

        inputs:
            curve
            point

        outputs:
            binormal vector"""
    p = curve.p
    dim = curve.dim
    localBasis = np.zeros((curve.n,curve.n),float)
    span = curve.FindSpan(s)
    curve.DersBasisFunc(span,s,localBasis[:,span-p:span+1])

    rp = ADtanVector(curve,s)

    ADvertices = curveAD(curve)
    
    rpp=[]
    for i in range(dim):
        rpp.append(np.dot(ADvertices[i],localBasis[2]))

    binormal = ADcrossproduct(rp,rpp)
    return binormal



def ADnormalVector(curve, s):
    """ length preserving tangent
        AD form

        inputs:
            curve
            point

        outputs:
            normal vector"""

    tg = ADtanVector(curve,s)
    binormal = ADbinormalVector(curve,s)
    normal = ADcrossproduct(binormal,tg)
        

    return normal




def utangent(self, u,v):
    """Automatic Differentiation
        of the u-tangent of a surface
        at a point"""


    dim = self.dim
    V   = self.surfNetAD
    
    knotsu  = self.uknots
    uorder  = self.uorder
    up      = uorder-1
    uspan   = self.uFindSpan(u)
    
    knotsv  = self.vknots
    vorder  = self.vorder
    vp      = vorder-1
    vspan   = self.vFindSpan(v)

    ulocalBasis = np.zeros((self.un,self.un),float)
    vlocalBasis = np.zeros((self.vn,self.vn),float)
    
    self.uDersBasisFunc(uspan,u,ulocalBasis[:,uspan-up:uspan+1])
    self.vDersBasisFunc(vspan,v,vlocalBasis[:,vspan-vp:vspan+1])
    
    Nu = ulocalBasis[1]
    Nv = vlocalBasis[0]
    M = np.outer(Nu,Nv)
    #ru = 0.
    ru=[]
    for d in range(dim):
        store=0.
        for j in range(len(knotsv)-vorder):
            for i in range(len(knotsu)-uorder):
                    store = V[i][j][d]*M[i,j] + store
        ru.append(store)

    return ru

def surf_uu(self, u,v):
    """AD surface uu"""
    dim = self.dim
    V   = self.surfNetAD
    
    knotsu  = self.uknots
    uorder  = self.uorder
    up      = uorder-1
    uspan   = self.uFindSpan(u)
    
    knotsv  = self.vknots
    vorder  = self.vorder
    vp      = vorder-1
    vspan   = self.vFindSpan(v)

    ulocalBasis = np.zeros((self.un,self.un),float)
    vlocalBasis = np.zeros((self.vn,self.vn),float)
    
    self.uDersBasisFunc(uspan,u,ulocalBasis[:,uspan-up:uspan+1])
    self.vDersBasisFunc(vspan,v,vlocalBasis[:,vspan-vp:vspan+1])

    Nu = ulocalBasis[2]
    Nv = vlocalBasis[0]
    M = np.outer(Nu,Nv)
    
    ruu=[]
    for d in range(dim):
        store=0.
        for j in range(len(knotsv)-vorder):
            for i in range(len(knotsu)-uorder):
                    store = V[i][j][d]*M[i,j] + store
        ruu.append(store)

    return ruu


def vtangent(self, u,v):
    """Automatic Differentiation
        of the v-tangent of a surface
        at a point"""


    dim = self.dim
    V   = self.surfNetAD
    
    knotsu  = self.uknots
    uorder  = self.uorder
    up      = uorder-1
    uspan   = self.uFindSpan(u)
    
    knotsv  = self.vknots
    vorder  = self.vorder
    vp      = vorder-1
    vspan   = self.vFindSpan(v)

    ulocalBasis = np.zeros((self.un,self.un),float)
    vlocalBasis = np.zeros((self.vn,self.vn),float)
    
    self.uDersBasisFunc(uspan,u,ulocalBasis[:,uspan-up:uspan+1])
    self.vDersBasisFunc(vspan,v,vlocalBasis[:,vspan-vp:vspan+1])
    
    Nu = ulocalBasis[0]
    Nv = vlocalBasis[1]
    M = np.outer(Nu,Nv)
    #ru = 0.
    rv=[]
    for d in range(dim):
        store=0.
        for j in range(len(knotsv)-vorder):
            for i in range(len(knotsu)-uorder):
                    store = V[i][j][d]*M[i,j] + store
        rv.append(store)

    return rv

def surf_vv(self, u,v):
    """AD surface vv"""
    dim = self.dim
    V   = self.surfNetAD
    
    knotsu  = self.uknots
    uorder  = self.uorder
    up      = uorder-1
    uspan   = self.uFindSpan(u)
    
    knotsv  = self.vknots
    vorder  = self.vorder
    vp      = vorder-1
    vspan   = self.vFindSpan(v)

    ulocalBasis = np.zeros((self.un,self.un),float)
    vlocalBasis = np.zeros((self.vn,self.vn),float)
    
    self.uDersBasisFunc(uspan,u,ulocalBasis[:,uspan-up:uspan+1])
    self.vDersBasisFunc(vspan,v,vlocalBasis[:,vspan-vp:vspan+1])

    Nu = ulocalBasis[0]
    Nv = vlocalBasis[2]
    M = np.outer(Nu,Nv)
    
    rvv=[]
    for d in range(dim):
        store=0.
        for j in range(len(knotsv)-vorder):
            for i in range(len(knotsu)-uorder):
                    store = V[i][j][d]*M[i,j] + store
        rvv.append(store)

    return rvv


def surf_uv(self, u,v):
    """AD surface uv= =vu"""
    dim = self.dim
    V   = self.surfNetAD
    
    knotsu  = self.uknots
    uorder  = self.uorder
    up      = uorder-1
    uspan   = self.uFindSpan(u)
    
    knotsv  = self.vknots
    vorder  = self.vorder
    vp      = vorder-1
    vspan   = self.vFindSpan(v)

    ulocalBasis = np.zeros((self.un,self.un),float)
    vlocalBasis = np.zeros((self.vn,self.vn),float)
    
    self.uDersBasisFunc(uspan,u,ulocalBasis[:,uspan-up:uspan+1])
    self.vDersBasisFunc(vspan,v,vlocalBasis[:,vspan-vp:vspan+1])

    Nu = ulocalBasis[1]
    Nv = vlocalBasis[1]
    M = np.outer(Nu,Nv)
    
    ruv=[]
    for d in range(dim):
        store=0.
        for j in range(len(knotsv)-vorder):
            for i in range(len(knotsu)-uorder):
                    store = V[i][j][d]*M[i,j] + store
        ruv.append(store)

    return ruv


def surf_normal(self, u,v):
    ru = utangent(self, u,v)
    rv = vtangent(self, u,v)
    nn = ADcrossproduct(ru,rv)
    return nn


def surf_gaussian_curvature(self,u,v):
    dim = self.dim

    ru  = utangent(self,u,v)
    rv  = vtangent(self,u,v)
    ruu = surf_uu(self, u,v)
    ruv = surf_uv(self, u,v)
    rvv = surf_vv(self, u,v)

    nn = ADcrossproduct(ru,rv)
    nn = ADnormalize(nn)
    E = np.dot(ru,ru)
    F = np.dot(ru,rv)
    G = np.dot(rv,rv)
    e = np.dot(ruu,nn)
    f = np.dot(ruv,nn)
    g = np.dot(rvv,nn)

    test = E*G - F**2
    if abs(test.value) > MINFLOAT:
        Kcurvature = (e*g - f**2) / (E*G - F**2)
    else:
        Kcurvature = 0.
        
    
    return Kcurvature

def surf_Christoffel(self,u,v):
    dim = self.dim

    ru  = utangent(self,u,v)
    rv  = vtangent(self,u,v)
    ruu = surf_uu(self, u,v)
    ruv = surf_uv(self, u,v)
    rvv = surf_vv(self, u,v)

    nn = ADcrossproduct(ru,rv)
    nn = ADnormalize(nn)
    E = np.dot(ru,ru)
    F = np.dot(ru,rv)
    G = np.dot(rv,rv)
    e = np.dot(ruu,nn)
    f = np.dot(ruv,nn)
    g = np.dot(rvv,nn)
    
    dscrmnt = E*G - F**2
    
    Gamma1 = ()/()
    return 

def surf_gaussian_curvature_quick(self,u,v):
    dim = self.dim
    
    dim = self.dim
    V   = self.surfNetAD
    
    knotsu  = self.uknots
    uorder  = self.uorder
    up      = uorder-1
    uspan   = self.uFindSpan(u)
    
    knotsv  = self.vknots
    vorder  = self.vorder
    vp      = vorder-1
    vspan   = self.vFindSpan(v)

    ulocalBasis = np.zeros((self.un,self.un),float)
    vlocalBasis = np.zeros((self.vn,self.vn),float)
    
    self.uDersBasisFunc(uspan,u,ulocalBasis[:,uspan-up:uspan+1])
    self.vDersBasisFunc(vspan,v,vlocalBasis[:,vspan-vp:vspan+1])
    
    Nu = ulocalBasis[1]
    Nv = vlocalBasis[0]
    M = np.outer(Nu,Nv)
    ru=[]
    for d in range(dim):
        store=0.
        for j in range(len(knotsv)-vorder):
            for i in range(len(knotsu)-uorder):
                    store = V[i,d,j]*M[i,j] + store
        ru.append(store)
    #ru  = utangent(self,u,v)
        
    Nu = ulocalBasis[0]
    Nv = vlocalBasis[1]
    M = np.outer(Nu,Nv)
    rv=[]
    for d in range(dim):
        store=0.
        for j in range(len(knotsv)-vorder):
            for i in range(len(knotsu)-uorder):
                    store = V[i,d,j]*M[i,j] + store
        rv.append(store)
    #rv  = vtangent(self,u,v)
     
    Nu = ulocalBasis[2]
    Nv = vlocalBasis[0]
    M = np.outer(Nu,Nv)
    ruu=[]
    for d in range(dim):
        store=0.
        for j in range(len(knotsv)-vorder):
            for i in range(len(knotsu)-uorder):
                    store = V[i,d,j]*M[i,j] + store
        ruu.append(store)
    #ruu = surf_uu(self, u,v)
        
    Nu = ulocalBasis[1]
    Nv = vlocalBasis[1]
    M = np.outer(Nu,Nv)
    
    ruv=[]
    for d in range(dim):
        store=0.
        for j in range(len(knotsv)-vorder):
            for i in range(len(knotsu)-uorder):
                    store = V[i,d,j]*M[i,j] + store
        ruv.append(store)
    #ruv = surf_uv(self, u,v)
        
    Nu = ulocalBasis[0]
    Nv = vlocalBasis[2]
    M = np.outer(Nu,Nv)
    rvv=[]
    for d in range(dim):
        store=0.
        for j in range(len(knotsv)-vorder):
            for i in range(len(knotsu)-uorder):
                    store = V[i,d,j]*M[i,j] + store
        rvv.append(store)
    #rvv = surf_vv(self, u,v)

    nn = ADcrossproduct(ru,rv)
    nn = ADnormalize(nn)
    E = np.dot(ru,ru)
    F = np.dot(ru,rv)
    G = np.dot(rv,rv)
    e = np.dot(ruu,nn)
    f = np.dot(ruv,nn)
    g = np.dot(rvv,nn)

    test = E*G - F**2
    if abs(test.value) > MINFLOAT:
        Kcurvature = (e*g - f**2) / (E*G - F**2)
    else:
        Kcurvature = 0.
        
    
    return Kcurvature

def surf_mean_curvature(self,u,v):
    dim = self.dim
    
    ru  = utangent(self,u,v)
    rv  = vtangent(self,u,v)
    ruu = surf_uu(self, u,v)
    ruv = surf_uv(self, u,v)
    rvv = surf_vv(self, u,v)

    nn = ADcrossproduct(ru,rv)
    nn = ADnormalize(nn)
    E = np.dot(ru,ru)
    F = np.dot(ru,rv)
    G = np.dot(rv,rv)
    e = np.dot(ruu,nn)
    f = np.dot(ruv,nn)
    g = np.dot(rvv,nn)

    test = E*G - F**2
    if abs(test.value) > MINFLOAT:
        Mcurvature = (E*g - f*F*2.0 +e*G)*0.5 / (E*G - F**2)
    else:
        Mcurvature = 0. 
    
    return Mcurvature



def main():

    v1 = np.array([[0.,0.,0.],
                   [1.,0.,0.],
                    [2.,0.,0.],
                    [3.,0.,0.],
                    [4.,0.,0.]])
                         
    v2 = np.array([[0.,1.,0.],
                   [1.,1.,0.],
                    [2.,1.,0.],
                    [3.,1.,0.],
                    [4.,1.,0.]])
                    
    v3 = np.array([[0.,2.,0.],
                   [1.,2.,0.],
                    [2.,2.,0.],
                    [3.,2.,0.],
                    [4.,2.,0.]])
    v4 = np.array([[0.,3.,0.],
                   [1.,3.,0.],
                    [2.,3.,0.],
                    [3.,3.,0.],
                    [4.,3.,0.]])
    v5 = np.array([[0.,4.,0.],
                   [1.,4.,0.],
                    [2.,4.,0.],
                    [3.,4.,0.],
                    [4.,4.,0.]])

    

    return #curve1
    


if __name__ == '__main__':
    
    """
    Procedure:
        strip off the 2D part of the curve:
            - transverse very in x and z
            - longitudinal vary in y and z
        optimize - make longitudinals hit the transverse vertices 
        - could do in 2D by fixing y at the loc of the transverse curve - then just make it hit (y,z)
                which looks like x,y to the 2D solver
                
        Need a geometry engine to make such frame of reference changes easy...
        
    """
    #curve = main()
    v1 = np.array([[0.,0.,.1],
                   [1.,0.,.1],
                    [2.,0.,.1],
                    [3.,0.,.1],
                    [4.,0.,.1]])
                         
    v2 = np.array([[0.,1.,.1],
                   [1.,1.,.1],
                    [2.,1.,.1],
                    [3.,1.,.1],
                    [4.,1.,.1]])
                    
    v3 = np.array([[0.,2.,.1],
                   [1.,2.,.1],
                    [2.,2.,.2],
                    [3.,2.,.1],
                    [4.,2.,.1]])
                    
    v4 = np.array([[0.,3.,.1],
                   [1.,3.,.1],
                    [2.,3.,.1],
                    [3.,3.,.1],
                    [4.,3.,.1]])
                    
    v5 = np.array([[0.,4.,.1],
                   [1.,4.,.1],
                    [2.,4.,.1],
                    [3.,4.,.1],
                    [4.,4.,.1]])
    k=4
    nump = 30
    tools = DimensionTools.ChangeDimension(v1)
    vlista = [v1,v2,v3,v4,v5]
    vlistb = []
    for vs in vlista:
        vlistb.append(tools.scaleVertices(vs, 2, 10.))
        
    tnet = np.asarray(vlistb)
    tCurveNet = []
    for points in tnet:
        tCurveNet.append(Bspline(points, k, nump))
        

        
    #surf = ADLS.LagrangeSurfaceMaker(Bspline, interpolatedBspline, BsplineSurface, tCurveNet,  lCurveNetOptimize = True, k=4, nump=30)
    surf = ADLS.LagrangeSurfaceMaker(tCurveNet,  lCurveNetOptimize = True, k=4, nump=30)

    self = surf
    SP = sfp.SurfaceFormParameter(self.surf)
    
    gs = .5
    gt = .5
    SP.add_GussianContraint(2.33,gs,gt)
    self.surf.FormParameters = SP.SurfFormParam
    
    adsurf2 = surf_curveAD(self.surf)
    self.surf.define_surfNetAD(adsurf2)
    gi = surf_gaussian_curvature_quick(self.surf,gs,gt)
    print 'initial gaussian curvature value = {}'.format(gi.value)


    n=17
    h=AD.scalarAD(0., np.matrix(np.zeros((n),float)), np.matrix(np.zeros((n,n),float)) )    
    
    self.optimize_surface()
    
    #"""
    check_loft = True
    if check_loft == True:
        
        
        
        
        adsurf1 = surfAD(self.surf)
        self.surf.define_surfNetAD(adsurf1)
        a = surf_gaussian_curvature(self.surf,gs,gt)
        
        adsurf2 = surf_curveAD(self.surf)
        self.surf.define_surfNetAD(adsurf2)
        b = surf_gaussian_curvature_quick(self.surf,gs,gt)
        
        #self.calc_meanCurvature() 
        
        print 'checking loft quality:'
        self.check_loft_quality()
        for el in surf.loft_quality:
            print el
        print ''
        
            
        pt=.5
        lpt=.5
        i=2
        print self.surf.evalSurface(pt,lpt)
        print self.surf.tcurveNet[i].CurvePoint(pt)
        for curve in self.surf.tcurveNet:
            points3d(curve.vertices[:,0],curve.vertices[:,1],curve.vertices[:,2], scale_factor=.03) #,colormap='greens')
        for curve in self.surf.lcurveNet:
            plot3d(curve.r[:,0],curve.r[:,1],curve.r[:,2], tube_radius=.004)
        for curve in self.surf.lcurveNet:
            points3d(curve.vertices[:,0],curve.vertices[:,1],curve.vertices[:,2], scale_factor=.01, colormap='OrRd')  #,colormap='oranges')
        
            
        for curve in self.surf.tcurveNet:
            plot3d(curve.r[:,0],curve.r[:,1],curve.r[:,2], tube_radius=.005)
        s0 = mesh(self.surf.surface[:,:,0],self.surf.surface[:,:,1],self.surf.surface[:,:,2], colormap='YlGnBu', opacity=.75)
    

        print 'match curve 4? ', self.surf.lcurveNet[4].CurvePoint(.25), self.surf.lcurveNet[4].CurvePoint(.75)
        print 'curve 4 vertices :'
        print self.surf.tcurveNet[4].vertices
    #"""
