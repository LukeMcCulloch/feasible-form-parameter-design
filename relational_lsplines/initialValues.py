# 20130731
# Lagrange Optimization of B-Splines
# Starting Values
# based on Form Parameters
# TLM python module


import numpy as np
import matplotlib.pyplot as plt
from functools import partial
import copy

from automatic_differentiation import ad
#from AD import *


small = 0.0000001



def interval_bounds(curve, small = 1.e-6):
    interval_data = {   'xmin':[],
                        'xmax':[],
                        'ymin':[],
                        'ymax':[]}
    xmin,xmax,ymin,ymax = curve.extreme_C0()
    xi = curve.vertices[ 0,0]
    xe = curve.vertices[-1,0]
    yi = curve.vertices[ 0,1]
    ye = curve.vertices[-1,1]
    for i in range(curve.n):
        interval_data['xmin'].append(xmin)
        interval_data['xmax'].append(xmax)
        interval_data['ymin'].append(ymin)
        interval_data['ymax'].append(ymax)
        
    interval_data['xmin'][ 0] = xi - small
    interval_data['xmax'][ 0] = xi + small
    interval_data['ymin'][ 0] = yi - small
    interval_data['ymax'][ 0] = yi + small
    
    #interval_data['xmin'][ 1] = interval_data['xmin'][ 1] +small
    #interval_data['xmax'][ 1] = interval_data['xmax'][ 1] -small
    #interval_data['xmin'][ 2] = interval_data['xmin'][ 2] +small
    #interval_data['xmax'][ 2] = interval_data['xmax'][ 2] -small
    
    #interval_data['xmin'][-2] = interval_data['xmin'][-2] +small
    #interval_data['xmax'][-2] = interval_data['xmax'][-2] -small
    #interval_data['xmin'][-3] = interval_data['xmin'][-3] +small
    #interval_data['xmax'][-3] = interval_data['xmax'][-3] -small
    
    interval_data['xmin'][-1] = xe - small
    interval_data['xmax'][-1] = xe + small
    interval_data['ymin'][-1] = ye - small
    interval_data['ymax'][-1] = ye + small
    return interval_data, small

def lagrangian_bounds(Lagrangian, interval_data, small = 1.e-6, large = 2.e5):
    interval_data['lmin'] = []
    interval_data['lmax'] = []
    for el in Lagrangian.equality:
        interval_data['lmin'].append(-large)
        interval_data['lmax'].append(large)
    return interval_data, small

def solveForTheta(ab_e, alpha, xe, Cab_given, sign, theta):
    """Newton Raphson routine to solve for a theta angle"""
    conv = 1.0
    count = 0
    #if abs(ab_e) < small:
    #    print 'error, ab_e = ',ab_e
    #    ab_e = 1.
    print 'ab_e = ',ab_e
    #store_denominator = (1+ab_e)**2+(1-ab_e**2)*np.cos(2.*np.arctan2(np.tan(alpha/ab_e)))
    #store_denominator = (1+ab_e)**2+(1-ab_e**2)*np.cos(2.*(alpha/ab_e))
    store_denominator = (1+ab_e**2)+(1-ab_e**2)*np.cos(2.*np.arctan2(np.tan(alpha),ab_e))
    store_scalar = sign*(25.*np.sqrt(2.))/(8.*xe)
    while (abs(conv) >.0000001 and count <50):

        Cab_found =  (  ( ((theta-alpha).sin()) *(-1.)*(1.+ab_e) )   +  (theta+alpha).sin()*(1.-ab_e)  ) *store_scalar/store_denominator
        
        f = Cab_found-Cab_given

        theta -= f.value/f.grad
        
        conv = f.value
        count +=1
    print count, f.value, theta.value
    return theta


def calc_curve_aspect_ratio():
    return

def calc_curve_points_for_tangent():
    return


    
class InitializeControlVertices(object):
    """Neglected Documentation, sorry.
    In general this implements Harries starting control
    vertices for form parameter generation of B-spline curves
    """
    def __init__(self, 
                 xb=0.,yb=12.,
                 xe=12.,ye=0.,
                 alphab=None,alphae=None,
                 Cab_given=None,Cae_given=None,
                 area=None,xc=None,yc=None,nCV = 7, 
                 slope= 'up', options = {}):
                     
        self.dim                = 2
        self.small              = 1.e-6
        self.verysmall          = 1.e-12
        self.geosmall           = 1.e-1
        self.maxit              = 50
        self.n                  = nCV
        self.chklst             = self.n*[False]
        self.slope              = slope
        self.bounds             = {}
        
        self.xb                 = xb
        self.yb                 = yb
        self.xe                 = xe
        self.ye                 = ye
        self.xwidth             = xe-xb 
        
        self.alphab = alphab
        self.alphae = alphae
        self.Cab_given = Cab_given
        self.Cae_given = Cae_given
        
        
        
        #assert abs(xb) < self.verysmall
        #assert yb > 0.        
        #assert abs(ye) < self.verysmall
        #assert xe > 0.
        
        self.endtangentsdefined = False
        if alphab is not None:
            self.alphab          = np.radians(alphab)
            self.endtangentsdefined = True 
        else:
            self.alphab = alphab
        if alphae is not None:
            self.alphae          = np.radians(alphae)
            self.endtangentsdefined = True
        else:
            self.alphae = alphae
            
        self.EndCurvaturesDefined = False
        if Cab_given != None:
            self.Cab_given          = Cab_given
            self.EndCurvaturesDefined = True
        if Cae_given != None:
            self.Cae_given          = Cae_given
            self.EndCurvaturesDefined = True
            
        self.XCentroidDefined = False
        if xc != None:
            self.xc                 = xc
            self.XCentroidDefined   = True
        self.YCentroidDefined = False
        if yc != None:
            self.yc                 = yc
            self.YCentroidDefined = False
            
        
        self.vertices           = np.zeros([self.n,self.dim],float)
        self.vertices[0,0]      = xb
        self.vertices[0,1]      = yb
        self.vertices[nCV-1,0]  = xe
        self.vertices[nCV-1,1]  = ye
        
            
        self.options                = options
        
        self.buildvertices()
        return
        
    def __call__(self):
        return
        
    def set_vertex(self, index, value, dim):
        return
        
    def buildvertices(self):
        self.vertices = self.linear_vertices((self.xb,self.yb),(self.xe,self.ye), self.n)
        self.chklst[0]          = True
        self.chklst[-1]         = True        
        self.ab_e = self.compute_AR()
        if self.endtangentsdefined:
            self.vertices[1],self.vertices[-2] = self.use_end_tangents()
            self.chklst[1]  = True
            self.chklst[-2] = True
        if self.EndCurvaturesDefined:
            self.vertices[2], self.vertices[-3] = self.use_end_curvature()
            self.chklst[2]  = True
            self.chklst[-3] = True
        if self.XCentroidDefined and self.YCentroidDefined:
            self.set_interior_vertices()
        for i, el in enumerate(self.chklst):
            vertices_to_fix = []
            if el==False:
                vertices_to_fix.append(i)
            
        return
        
    def linear_vertices(self, start, end, num):
        start = np.asarray(start)
        end = np.asarray(end)
        #dim = len(start)
        dim = self.dim 
        vertices = []
        for i in range(dim):
            vi = np.linspace(start[i],end[i],num)
            vertices.append(vi)
        vertices = np.asarray(vertices)
        vertices = np.transpose(vertices)
        return vertices
        
    def compute_AR(self):
        xb = self.vertices[0,0]
        yb = self.vertices[0,1]
        xe = self.vertices[-1,0]
        ye = self.vertices[-1,1]
        if 'SAC' in self.options:
            ab_e = self.options['SAC']['maxwidth']/abs(xe-xb)
        elif 'DWL' in self.options:
            ab_e = self.options['DWL']['maxwidth']/abs(xe-xb)
        elif 'maxwidth' in self.options:
            ab_e = self.options['maxwidth']/abs(xe-xb)
        else:
            ab_e = abs(ye-yb)/abs(xe-xb)
        return ab_e
    
    def use_end_tangents(self):
        xb = self.vertices[0,0]
        yb = self.vertices[0,1]
        xe = self.vertices[-1,0]
        ye = self.vertices[-1,1]
        #L_polygon = (xe-xb)*np.sqrt(1.+ab_e)
        ab_e    = self.ab_e
        xwidth  = self.xwidth
        save1   = np.arctan2(np.tan(self.alphab),ab_e)
        save2   = np.arctan2(np.tan(self.alphae),ab_e)
        db      = (2.*np.sqrt(2.)/25.) * xwidth * np.sqrt( ((np.cos(save1))**2) + (ab_e**2)*((np.sin(save1))**2) )
        de      = (2.*np.sqrt(2.)/25.) * xwidth * np.sqrt( ((np.cos(save2))**2) + (ab_e**2)*((np.sin(save2))**2) )
        
        #recomended coordinates for the immediate interior points on the curve:
        print 'self.slope = ',self.slope
        if self.slope == 'up':
            x1 = xb + db*np.cos(self.alphab)
            x5 = xe - de*np.cos(self.alphae)
            y1 = yb + db*np.sin(self.alphab)    
            y5 = ye - de*np.sin(self.alphae)
        elif self.slope == 'down':
            x1 = xb + db*np.cos(self.alphab)
            x5 = xe - de*np.cos(self.alphae)
            y1 = yb - db*np.sin(self.alphab)    
            y5 = ye + de*np.sin(self.alphae)
        return np.asarray([x1,y1]), np.asarray([x5,y5])
        
    def use_end_curvature(self):
        xb = self.vertices[0,0]
        x1 = self.vertices[1,0]
        x2 = self.vertices[2,0]
        x4 = self.vertices[-3,0]
        x5 = self.vertices[-2,0]
        xe = self.vertices[-1,1]
        yb = self.vertices[0,1]
        y1 = self.vertices[1,1]
        y2 = self.vertices[2,1]
        y4 = self.vertices[-3,1]
        y5 = self.vertices[-2,1]
        ye = self.vertices[-1,1]
        ab_e = self.ab_e
        Cab_given = self.Cab_given
        Cae_given = self.Cae_given

        if self.alphab is None:
            alphab = 0.
        else:
            alphab = self.alphab
    
        if self.alphae is None:
            alphae = 0.
        else:
            alphae = self.alphae
        
        xb = self.vertices[0,0]
        yb = self.vertices[0,1]
        xe = self.vertices[-1,0]
        ye = self.vertices[-1,1]
        
        # curvature at the start:
        theta1 = ad(0.,1.,0.) #initialize alpha1 - we need alpha1 such that x2 sets the curvature = to the desired curvature
        theta1 = self.solveForTheta(ab_e, alphab, (xe-xb), Cab_given, -1., theta1)
        self.theta1 = theta1.value
        
        x2 = theta1.cos()*(3.*np.sqrt(2.)/25.)*(xe-xb) + x1
        y2 = theta1.sin()*(3.*np.sqrt(2.)/25.)*(xe-xb)*ab_e + y1
        
        x2 = x2.value
        y2 = y2.value
        
        # curvature at the end:
        theta4 = ad(0.,1.,0.)
        theta4 = solveForTheta(ab_e, alphae, (xe-xb), Cae_given, 1., theta4)
        self.theta4 = theta4.value
        
        x4 = theta4.cos()*(-1.)*(3.*np.sqrt(2.)/25.)*(xe-xb) + x5
        y4 = theta4.sin()*(-1.)*(3.*np.sqrt(2.)/25.)*(xe-xb)*ab_e + y5
    
        x4 = x4.value
        y4 = y4.value
        return np.asarray([x2,y2]),np.asarray([x4,y4])
        
    def set_interior_vertices(self):
        xb = self.vertices[0,0]
        x1 = self.vertices[1,0]
        x2 = self.vertices[2,0]
        x4 = self.vertices[-3,0]
        x5 = self.vertices[-2,0]
        xe = self.vertices[-1,1]
        yb = self.vertices[0,1]
        y1 = self.vertices[1,1]
        y2 = self.vertices[2,1]
        y4 = self.vertices[-3,1]
        y5 = self.vertices[-2,1]
        ye = self.vertices[-1,1]
        if self.n ==7:
            term1 = (1./(6.*area))
            term2 = 2.*area*(x2+x4) + x1**2*yb - x1*x2*yb - x1*x4*yb
            term3 = x1*x2*y1 - x2*x4*y1 - x1**2*y2 + x1*x4*y2 - x2*x5*y4 + x5**2*y4 + xe**2*y5
            term4 = -xe*x2*y5 - xe*x4*y5 + x2*x4*y5 + xe*x5*y5 - x4*x5*y5
            termA = term1*(term2+term3+term4)
            
            termB = (1./(6.*area))*(2.*area - x1*yb -x2*y1 + x1*y2 - x4*y2 + x2*y4 - x5*y4 - xe*y5 +x4*y5 )
            
            x3 = (xc - termA)/termB
            if x3>x4 or x3<x2:
                print 'avoiding singularity'
                x3 = (x2+x4)/2.
    
            y3  = ((-2.*area + x1*yb + x2*y1 -x1*y2 + xe*y5 + x5*y4 -x4*y5 )/(x2-x4)) 
            y3 += x3*(y2-y4)/(x2-x4)
            self.vertices[3] = np.asarray([x3,y3])
            self.chklst[3] = True
            return  
            
        elif self.n ==6:
            return
            
        elif self.n==8:
            dummy = self.linear_vertices( (x2,y2), (x4,y4), 4)
            self.vertices[3]    = dummy[1]
            self.vertices[4]    = dummy[2]
            self.chklst[3]      = True
            self.chklst[4]      = True
            return 
            
    def bound_x_midpoints(self):
        """
            Given area (e.g. from the SAC curve)
            Use Equation D.33 on page 166 to bound allowable constraints on xc
        """
        xb = self.vertices[0,0]
        x1 = self.vertices[1,0]
        x2 = self.vertices[2,0]
        x4 = self.vertices[-3,0]
        x5 = self.vertices[-2,0]
        xe = self.vertices[-1,1]
        yb = self.vertices[0,1]
        y1 = self.vertices[1,1]
        y2 = self.vertices[2,1]
        y4 = self.vertices[-3,1]
        y5 = self.vertices[-2,1]
        ye = self.vertices[-1,1]
        #area = self.area
        
        x3min = x2
        x3max = x4
        y3min = min(y2,y4)
        y3max = max(y2,y4)
        
        areamin = self.calc_area(x3min,y3min)
        areamax = self.calc_area(x3max,y3max)
        
        qmin = partial(self.calc_xc, x3min)
        qmax = partial(self.calc_xc, x3max)
        
        xcmin = qmin( areamin)
        xcmax = qmax( areamax)
        
        print 'bounds'
        print 'min area = ', areamin
        print 'max area = ', areamax
        print 'min xc = ', xcmin
        print 'max xc = ',xcmax
        
        self.bounds['min_area'] = areamin
        self.bounds['max_area'] = areamax
        self.bounds['min_xc'] = xcmin
        self.bounds['max_xc'] = xcmax
        
        
        #xcmin = self.calc_xc(x3min, areamin)
        #xcmax = self.calc_xc(x3max, areamax)
        return xcmin, xcmax
    
    def calc_curvature(self):
        """
            tbd:  give a range for the curvature parameters
                contingent on the aspect ratio of the curve
                
                
        Curvature Flow - equation:
            
            d_s Q(s) = -Grad E(Q(s))
            
            where d_s is the derivative w/r/t parameterization or is it 'time'?
            Grad is the spatial gradient of the immersion in ambient space
            E is some energy (dirichlet, wilmore, bending, etc.)
        """
        xb = self.vertices[0,0]
        x1 = self.vertices[1,0]
        x2 = self.vertices[2,0]
        x4 = self.vertices[-3,0]
        x5 = self.vertices[-2,0]
        xe = self.vertices[-1,1]
        yb = self.vertices[0,1]
        y1 = self.vertices[1,1]
        y2 = self.vertices[2,1]
        y4 = self.vertices[-3,1]
        y5 = self.vertices[-2,1]
        ye = self.vertices[-1,1]
        
        theta1 = self.theta1
        theta4 = self.theta4
        
        
        
        return
    
    def calc_area(self, x3, y3):
        xb = self.vertices[0,0]
        x1 = self.vertices[1,0]
        x2 = self.vertices[2,0]
        x4 = self.vertices[-3,0]
        x5 = self.vertices[-2,0]
        xe = self.vertices[-1,1]
        yb = self.vertices[0,1]
        y1 = self.vertices[1,1]
        y2 = self.vertices[2,1]
        y4 = self.vertices[-3,1]
        y5 = self.vertices[-2,1]
        ye = self.vertices[-1,1]
        
        # assume:  xb = ye = 0
        area = 0.5*(x1*yb + x2*y1 - x1*y2 + xe*y5 + x5*y4 - x4*y5 + (x4-x2)*y3 + (y2-y4)*x3)
        return area
    
    def calc_xc(self, x3, area):
        xb = self.vertices[0,0]
        x1 = self.vertices[1,0]
        x2 = self.vertices[2,0]
        x4 = self.vertices[-3,0]
        x5 = self.vertices[-2,0]
        xe = self.vertices[-1,1]
        #
        yb = self.vertices[0,1]
        y1 = self.vertices[1,1]
        y2 = self.vertices[2,1]
        y4 = self.vertices[-3,1]
        y5 = self.vertices[-2,1]
        ye = self.vertices[-1,1]

        #term1 = (1./(area*6.))
        term1 = (area**-1) * (1./6.)
        if isinstance(term1,list):
            term1 = term1[0]
        term2 = 2.*area*(x2+x4) + x1**2*yb - x1*x2*yb - x1*x4*yb
        term3 = x1*x2*y1 - x2*x4*y1 - x1**2*y2 + x1*x4*y2 - x2*x5*y4 + x5**2*y4 + xe**2*y5
        term4 = -xe*x2*y5 - xe*x4*y5 + x2*x4*y5 + xe*x5*y5 - x4*x5*y5
        termA = term1*(term2+term3+term4)
        
        #termB = (1./(6.*area))*(2.*area - x1*yb -x2*y1 + x1*y2 - x4*y2 + x2*y4 - x5*y4 - xe*y5 +x4*y5 )
        
        termB = (2.*area - x1*yb -x2*y1 + \
                     x1*y2 - x4*y2 + x2*y4 - \
                     x5*y4 - xe*y5 +x4*y5 )*term1
        
        xc = termA + termB*x3
        return xc
    
    
    def solveForTheta(self,ab_e, alpha, xe, Cab_given, sign, theta):
        """Newton Raphson routine to solve for a theta angle"""
        if alpha is None:
            alpha = np.deg2rad(45.)
        conv = 1.0
        count = 0
        #if abs(ab_e) < small:
        #    print 'error, ab_e = ',ab_e
        #    ab_e = 1.
        print 'ab_e = ',ab_e
        #store_denominator = (1+ab_e)**2+(1-ab_e**2)*np.cos(2.*np.arctan2(np.tan(alpha/ab_e)))
        #store_denominator = (1+ab_e)**2+(1-ab_e**2)*np.cos(2.*(alpha/ab_e))
        store_denominator = (1+ab_e**2)+(1-ab_e**2)*np.cos(2.*np.arctan2(np.tan(alpha),ab_e))
        store_scalar = sign*(25.*np.sqrt(2.))/(8.*xe)
        while (abs(conv) >self.small and count <self.maxit):
    
            Cab_found =  (  ( ((theta-alpha).sin()) *(-1.)*(1.+ab_e) )   +  (theta+alpha).sin()*(1.-ab_e)  ) *store_scalar/store_denominator
            
            f = Cab_found-Cab_given
    
            theta -= f.value/f.grad
            
            conv = f.value
            count +=1
        print count, f.value, theta.value
        return theta

class InitializeControlVerticesDownSlope(object):
    """Neglected Documentation, sorry.
    In general this implements Harries starting control
    vertices for form parameter generation of B-spline curves
    """
    def __init__(self, 
                 xb=0.,yb=12.,xe=12.,ye=0.,
                 alphab=None,alphae=None,
                 Cab_given=None,Cae_given=None,
                 area=None,xc=None,yc=None,
                 nCV = 7, 
                 slope= 'down', options = {}):
        self.dim                = 2
        self.small              = 1.e-6
        self.verysmall          = 1.e-12
        self.geosmall           = 1.e-1
        self.maxit              = 50
        self.n                  = nCV
        self.chklst             = self.n*[False]
        self.slope              = slope
        self.bounds             = {}
        
        self.xb                 = xb
        self.yb                 = yb
        self.xe                 = xe
        self.ye                 = ye
        self.xwidth             = xe-xb 
        self.endtangentsdefined = False
        if alphab != None:
            self.alphab          = np.radians(alphab)
            self.endtangentsdefined = True            
        if alphae != None:
            self.alphae          = np.radians(alphae)
            self.endtangentsdefined = True
            
        self.EndCurvaturesDefined = False
        if Cab_given != None:
            self.Cab_given          = Cab_given
            self.EndCurvaturesDefined = True
        if Cae_given != None:
            self.Cae_given          = Cae_given
            self.EndCurvaturesDefined = True
            
        self.XCentroidDefined = False
        if xc != None:
            self.xc                 = xc
            self.XCentroidDefined   = True
        self.YCentroidDefined = False
        if yc != None:
            self.yc                 = yc
            self.YCentroidDefined = False
            
        
        self.vertices           = np.zeros([self.n,self.dim],float)
        self.vertices[0,0]      = xb
        self.vertices[0,1]      = yb
        self.vertices[nCV-1,0]  = xe
        self.vertices[nCV-1,1]  = ye
        
            
        #self.options                = options
        #self.buildvertices()
        return
        
def InitializeControlPoints(
        xb,yb,alphab,Cab_given,xe,ye,alphae,Cae_given,
        xc,yc,area,type_specifics = {}):
    """Neglected Documentation, sorry.
    In general this implements Harries starting control
    vertices for form parameter generation of B-spline curves
    """
    xwidth = xe-xb    
    # Main program to find an Initial Guess of Vertex Positions

    alphab = np.radians(alphab)
    alphae = np.radians(alphae)
    
    # Immediate Interior control points:
    """Given the first and last vertex, the 2nd and 2nd to last vertices are
        linked to the known vertices as follows:"""
        
    # aspect ratio of the enclosing box:
    if 'SAC' in type_specifics:
        ab_e = type_specifics['SAC']['maxwidth']/abs(xe-xb)
    elif 'DWL' in type_specifics:
        ab_e = type_specifics['DWL']['maxwidth']/abs(xe-xb)
    else:
        ab_e = abs(ye-yb)/abs(xe-xb) # aspect ration of enclosing box: 
                                     #   ~1. for an approximately quadratic region

    L_polygon = (xe-xb)*np.sqrt(1.+ab_e)
    
    save1   = np.arctan2(np.tan(alphab),ab_e)
    save2   = np.arctan2(np.tan(alphae),ab_e)
    db      = (2.*np.sqrt(2.)/25.) * xwidth * np.sqrt( ((np.cos(save1))**2) + (ab_e**2)*((np.sin(save1))**2) )
    de      = (2.*np.sqrt(2.)/25.) * xwidth * np.sqrt( ((np.cos(save2))**2) + (ab_e**2)*((np.sin(save2))**2) )

    #recomended coordinates for the immediate interior points on the curve:
    x1 = xb + db*np.cos(alphab)
    x5 = xe - de*np.cos(alphae)
    y1 = yb + db*np.sin(alphab)    
    y5 = ye - de*np.sin(alphae)



    ## Solve for the 3rd vertex from each end:
        # Curvature at the ends (if specified there) 
        # determines the location of the 3rd and 3rd to last points.
        # Solve the non-linear equation for curvature via newton's method. 
        #    using   function   "solveForTheta()"

    # curvature at the start:
    theta1 = ad(0.,1.,0.) #initialize alpha1 - we need alpha1 such that x2 sets the curvature = to the desired curvature
    theta1 = solveForTheta(ab_e, alphab, (xe-xb), Cab_given, -1., theta1)

    x2 = theta1.cos()*(3.*np.sqrt(2.)/25.)*(xe-xb) + x1
    y2 = theta1.sin()*(3.*np.sqrt(2.)/25.)*(xe-xb)*ab_e + y1

    x2 = x2.value
    y2 = y2.value


    # curvature at the end:
    theta4 = ad(0.,1.,0.)
    theta4 = solveForTheta(ab_e, alphae, (xe-xb), Cae_given, 1., theta4)
        
    x4 = theta4.cos()*(-1.)*(3.*np.sqrt(2.)/25.)*(xe-xb) + x5
    y4 = theta4.sin()*(-1.)*(3.*np.sqrt(2.)/25.)*(xe-xb)*ab_e + y5

    x4 = x4.value
    y4 = y4.value



    #Now solve for the mid point using the integral parameters
    #    denominator = (2.*area -x1*yb -x2*y1 +x1*y2 -x4*y2 +x2*y4 +x4*y5-x5*y4-xe*y5 )
    #    store1=6.*xc*area/denominator
    #    store2=(2.*area*(x2+x4)+x1*x1*yb-x1*x2*yb-x1*x4*yb+x1*x2*y1)/denominator
    #    store3=(-x2*x4*y1-x1*x1*y2+x1*x4*y2-x2*x5*y4+x5*x5*y4+xe*xe*y5)/denominator
    #    store4=(-xe*x2*y5-xe*x4*y5+x2*x4*y5+xe*x5*y5-x4*x5*y5)/denominator
    #    x3 = store1-store2-store3-store4
    term1 = (1./(6.*area))
    term2 = 2.*area*(x2+x4) + x1**2*yb - x1*x2*yb - x1*x4*yb
    term3 = x1*x2*y1 - x2*x4*y1 - x1**2*y2 + x1*x4*y2 - x2*x5*y4 + x5**2*y4 + xe**2*y5
    term4 = -xe*x2*y5 - xe*x4*y5 + x2*x4*y5 + xe*x5*y5 - x4*x5*y5
    termA = term1*(term2+term3+term4)
    
    termB = (1./(6.*area))*(2.*area - x1*yb -x2*y1 + x1*y2 - x4*y2 + x2*y4 - x5*y4 - xe*y5 +x4*y5 )
    
    
    x3 = (xc - termA)/termB
    if x3>x4 or x3<x2:
        x3 = (x2+x4)/2.
    

    #y3 = ((-2.*area+x1*yb+x2*y1-x1*y2+xe*y5+x5*y4-x4*y5)/(x2-x4)) + x3*(y2-y4)/(x2-x4)
    y3  = ((-2.*area + x1*yb + x2*y1 -x1*y2 + xe*y5 + x5*y4 -x4*y5 )/(x2-x4)) 
    y3 += x3*(y2-y4)/(x2-x4)


    vertices = np.zeros((7,2),float)

    vertices[0,0]=xb
    vertices[1,0]=x1
    vertices[2,0]=x2
    vertices[3,0]=x3
    vertices[4,0]=x4
    vertices[5,0]=x5
    vertices[6,0]=xe

    vertices[0,1]=yb
    vertices[1,1]=y1
    vertices[2,1]=y2
    vertices[3,1]=y3
    vertices[4,1]=y4
    vertices[5,1]=y5
    vertices[6,1]=ye

    #"""
    try:
        print 'checking starting vertices'
        xmax = type_specifics['xmax']
        xmin = type_specifics['xmin']
        ymax = type_specifics['ymax']
        ymin = type_specifics['ymin']
        for vert in vertices:
            if  vert[0]<xmin:
                vert[0]=xmin
            if  vert[0]>xmax:
                vert[0]=xmax
            if  vert[1]<ymin:
                vert[1]=ymin
            if  vert[1]>ymax:
                vert[1]=ymax
        print 'augmented starting vertices'
    except:
        pass
    #"""
    #"""
    for i in range(1,len(vertices)-1):
        if vertices[i,0] <vertices[i-1,0]:
            vertices[i,0] = (vertices[i-1,0]+vertices[i+1,0])/2.
    #"""
    return vertices

if __name__ == '__main__':
    from curve import Bspline
    # We will have a 7 control point curve, Order 4, as recommended by Harries, 1998.
    # Given The Following Form Parameters:

    ## Curve Differential parameters:

    xb=0.               #1st vertex: x value
    yb=12.               #1st vertex: y value
    alphab  = 0.        #tangent at the start of the curve
    Cab_given=0.      #curvature desired at the beggining

    xe=12.              #last vertex: x value
    ye=0.              #last vertex: y value
    alphae  = 0.        #tangent at the end of the curve
    Cae_given=0.       #curvature desired at the end

    ## Curve Integral parameters:

    area    = 72.
    xc      = 4.
    yc      = 4.
    
    ini_CV = InitializeControlVertices(xb,yb,xe,ye,alphab,alphae,Cab_given,Cae_given,area,xc,yc)
    curve = Bspline(ini_CV.vertices, 4, 50)
    print curve.vertices
    curve.plotcurve()
    plt.plot(curve.r[:,0],curve.r[:,1])
    plt.axis('equal')
    plt.plot(curve.vertices[:,0],curve.vertices[:,1])
    
    c1 = InitializeControlVertices(xb,yb,xe,ye,alphab,alphae,Cab_given,Cae_given,area)
    curve1 = Bspline(c1.vertices, 4, 50)
    curve1.plotcurve()
    plt.plot(curve1.r[:,0],curve1.r[:,1])
    plt.axis('equal')
    plt.plot(curve1.vertices[:,0],curve1.vertices[:,1])
    
    c2 = InitializeControlVertices(xb,yb,xe,ye,alphab,alphae,Cab_given,Cae_given)
    curve2 = Bspline(c2.vertices, 4, 50)
    curve2.plotcurve()
    plt.plot(curve2.r[:,0],curve2.r[:,1])
    plt.axis('equal')
    plt.plot(curve2.vertices[:,0],curve2.vertices[:,1])
    c2.bound_x_midpoints()
    self = c2

    
    c4 = InitializeControlVertices(xb,yb,xe,ye)
    curve4 = Bspline(c4.vertices, 4, 50)
    curve4.plotcurve()
    plt.plot(curve4.r[:,0],curve4.r[:,1])
    plt.axis('equal')
    plt.plot(curve4.vertices[:,0],curve4.vertices[:,1])
    
    
    for angle in [0.,30.,60.,90.,120.]:
        alphae = angle
        alphab = angle
        c5 = InitializeControlVertices(xb,yb,xe,ye,alphab,alphae)
        curve5 = Bspline(c5.vertices, 4, 50)
        curve5.plotcurve()
        plt.plot(curve5.r[:,0],curve5.r[:,1])
        plt.axis('equal')
        plt.plot(curve5.vertices[:,0],curve5.vertices[:,1])
    
    """
    for angle in [0.,20.,40.,60.,80.,100.]:
        alphae = angle
        alphab = angle
    
        vertices = InitializeControlPoints(xb,yb,alphab,Cab_given,xe,ye,alphae,Cae_given,xc,yc,area)
        curve = Bspline(vertices, 4, 50)
        print curve.vertices
        #curve.plotcurve()
        plt.plot(curve.r[:,0],curve.r[:,1])
        plt.axis('equal')
        plt.plot(vertices[:,0],vertices[:,1])
    
        
    for Curvature in [0.,20.,40.,60.,80.,100.]:
        Cab_given = Curvature
        Cae_given = Curvature
    
        vertices = InitializeControlPoints(xb,yb,alphab,Cab_given,xe,ye,alphae,Cae_given,xc,yc,area)
        curve = Bspline(vertices, 4, 50)
        print curve.vertices
        #curve.plotcurve()
        plt.plot(curve.r[:,0],curve.r[:,1])
        plt.axis('equal')
        plt.plot(vertices[:,0],vertices[:,1])
        
    for Centroid in [3.,3.5,4.,4.5,5.]:
        xc = Centroid
        yc = Centroid
        vertices = InitializeControlPoints(xb,yb,alphab,Cab_given,xe,ye,alphae,Cae_given,xc,yc,area)
        curve = Bspline(vertices, 4, 50)
        print curve.vertices
        #curve.plotcurve()
        plt.plot(curve.r[:,0],curve.r[:,1])
        plt.axis('equal')
        plt.plot(vertices[:,0],vertices[:,1])
        
    for area in [50.,55.,65.,70.,72.,75.,80.]:
    
        vertices = InitializeControlPoints(xb,yb,alphab,Cab_given,xe,ye,alphae,Cae_given,xc,yc,area)
        curve = Bspline(vertices, 4, 50)
        print curve.vertices
        #curve.plotcurve()
        plt.plot(curve.r[:,0],curve.r[:,1])
        plt.axis('equal')
        plt.plot(vertices[:,0],vertices[:,1])

    #"""