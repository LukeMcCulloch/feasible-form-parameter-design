#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 15:24:34 2018

@author: luke

lessons learned so far:
    1.)  the extended ia module, 
            so helpful for Hull consistency and relational programming, 
            is failing to give correct hessians for simple multilication
            or at least when powers are involved... hmmm
    2.) box consistency is not giving the results 
            stated in Hansen.  expansion point work needed to replicate max accuaracy
            though my stuff gets there in 6 iterations vs 10 for them.  (example book page 228)
            it uses three expansion pts per iteration
            
            
ah... https://stackoverflow.com/questions/2350072/custom-data-types-in-numpy-arrays
"""

from initialValues import InitializeControlVertices
import matplotlib.pyplot as plt

#from interval_arithmetic import ia as ia
from extended_interval_arithmetic import ia

from automatic_differentiation import ad
from automatic_differentiation import IntervalAnalysis


import numpy as np
import copy

        
#def hansen_gradient(f):
#    """
#    """
#    hansen_grad = []  #page 186, page 254
#    for i in range(np.size(f.grad)):
#        hansen_grad.append(ia(f.grad[0,i],f.grad[0,i]))
#    hansen_grad = np.asarray(hansen_grad)
#    return hansen_grad


def plot(ini_CV):
    curve = Bspline(ini_CV.vertices, 4, 50)
    print curve.vertices
    curve.plotcurve()
    plt.plot(curve.r[:,0],curve.r[:,1])
    plt.axis('equal')
    plt.plot(curve.vertices[:,0],curve.vertices[:,1])
    return

def montonic_compute_scalar_newton(func,
                                   c,
                                   x, 
                                   vertex_tpl, 
                                   itmax = 50, 
                                   tol = 1.e-15):
    """The scalar valued interval Newton's method
    TODO: make this exploit monotonicity
    
    
        monotonicity checking box consistency
        (monotonic interval newton iteration on a single DOF)
        
    inputs
    ----------
    func = function to be minimized/extremized
    c = constraints in f(x,c)
    x = free variables in f(x,c)
    vertex_tpl = 
    max_it
    """
    xi = copy.deepcopy(x)
    i=0
    not_converged = True
    xyz = vertex_tpl[0]     # which dimension are we in [0,1,2]?
    index = vertex_tpl[1]   # which vertex are we attempting to contract?
    n = 1##len(x[0])
    gn = xyz*n + index      #pick the correct element of the gradient
    while(i<itmax and not_converged):
        #midx                    = xi[xyz][index].value.getpoint(.5)
        midx                    = xi.value.getpoint(.5)
        xm                      = copy.deepcopy(xi)
        #xm[xyz][index].value    = midx
        xm.value    = midx
        f                       = func(xm,c) #could change to
        F                       = func(xi,c) # func(Xci,area,Xcgiven)
        # but maybe this is off the track enough to fix a 
        # new Xc_given 
        # via function closure in the definition
        # before every time we come here?
        
        nx = xm - f.value/F.grad[0,gn]
        
        xiold = copy.deepcopy(xi)
        
        #xi[xyz][index].value = xi[xyz][index].value & nx.value
        xi.value = xi.value & nx.value
        
        #not_converged = IntervalAnalysis.some_convergence_criteria(
        #                        xiold[xyz][index],xi[xyz][index],tol)
        not_converged = IntervalAnalysis.some_convergence_criteria(
                                xiold,xi,tol)
        
        i+=1
    return xi, i



def BC_newton(func,
               x, 
               vertex_tpl, 
               itmax = 50, 
               tol = 1.e-15):
    """The scalar valued interval Newton's method
    TODO: make this exploit monotonicity
    
    
        monotonicity checking box consistency
        (monotonic interval newton iteration on a single DOF)
        
    inputs
    ----------
    func = function to be minimized/extremized
    c = constraints in f(x,c)
    x = free variables in f(x,c)
    vertex_tpl = [0,0]
    max_it
    """
    xi = copy.deepcopy(x)
    i=0
    not_converged = True
    xyz = vertex_tpl[0]     # which dimension are we in [0,1,2]?
    #index = vertex_tpl[1]   # which vertex are we attempting to contract?
    n = 1##len(x[0])
    gn = xyz*n      #pick the correct element of the gradient
    while(i<itmax and not_converged):
        #for xpt in [.00001,.5,.99999]:
        for xpt in [.5]:
            xm                      = copy.deepcopy(xi)
            midx                    = xi[xyz].value.getpoint(xpt)
            xm[xyz].value           = midx
            f                       = func(xm) #
            F                       = func(xi) #
            
            if False:
                inner_grad      = ia(f.grad[0,gn].inf,f.grad[0,gn].sup)
                inner_nx        = xm[gn].value - f.value/inner_grad
                test_x_inner    = xi[xyz].value & inner_nx
                inner_mid       = test_x_inner.getpoint(.5)
                inner_midscalar = copy.deepcopy(xi)
                inner_midscalar[xyz].value = inner_mid
                inner_midvec    = copy.deepcopy(xi)
                f_innerscalar   = func(inner_midscalar)
                f_inner         = func(inner_midvec)
            #newgrad = ia()
            
            nx = xm[gn].value - f.value/F.grad[0,gn]
            #nx = xm[gn].value - F.value/F.grad[0,gn]
            
            xiold = copy.deepcopy(xi)
            
            xi[xyz].value = xi[xyz].value & nx
        
        not_converged = IntervalAnalysis.some_convergence_criteria(
                                xiold[xyz],xi[xyz],tol)
        
        i+=1
    return xi, i, F


def var_to_ad_ia(var_inf,var_sup,N,dim, name):
    return ad( ia(var_inf,var_sup), 
              name = name, 
              N=N, 
              dim = dim)



def F_Xc(self, area):
    """fix Xc given 
    via function closure
    """
    def closure(x3,Xc_constraint):
        return self.calc_xc(x3,area)-Xc_constraint
    return closure

def interval_lagrange_tester(X):
    """
        Accepts a vector of ia numbers
        returns a vector of ad ia numbers 
        set up for the AD class 
        and Lagrangian methods
    """
    n = np.size(X)
    interval_ad_vector = []
    for i, x in enumerate(X):
        number = ad( x, dim = i, N = n)
        interval_ad_vector.append(number)
    return interval_ad_vector

def test_HC_BC_2(V):
    """
        constraint:
        f(x,y) = x^3 + 100x + 10y = 0
        
        
    #"""
    x = V[0]
    y = V[1]
    return 10.*y + x**3 + 100.*x

def test_2():
    """
    
       func = F
       x = X
       vertex_tpl = [0,0]
       itmax = 50
       tol = 1.e-15
       xpt = .5
    """
    #xo = ia(-100.,100.)
    #yo = ia(-100.,100.)
    
    xo = ad( ia(-100., 100.), name = 'x', N=2, dim=0)
    yo = ad( ia(-100., 100.), name = 'y', N=2, dim=1)
    
    #X = interval_lagrange_tester([xo,yo])
    X = [xo,yo]
    F = test_HC_BC_2
    g = F(X)
    
    ans, i, F = BC_newton(F,X,[0,0],tol=1.e-20)
    print 'test 2'
    print 'ans = ',ans
    print 'i iterations = ',i
    print 'F = ',F
    print F.grad
    print F.hess
    return


def test(xshift = 0.,
         xb     = 0.,  #1st vertex: x value
         xe     = 12., #last vertex: x value
         yb     = 12.,
         ye     = 0.,
         area   = 72.,
         xc     = 4.,
         yc     = 4.):
    
    xb=xb+xshift        
    xe=xe+xshift        
    
    alphab  = 0.        #tangent at the start of the curve
    Cab_given=0.        #curvature desired at the beggining

    
    alphae  = 0.        #tangent at the end of the curve
    Cae_given=0.       #curvature desired at the end

    ## Curve Integral parameters:

    
    ini_CV = InitializeControlVertices(xb,yb,
                                       xe,ye,
                                       alphab,alphae,
                                       Cab_given,Cae_given,
                                       area,
                                       xc,yc)
    
    curve = Bspline(ini_CV.vertices, 4, 50)
    plot(ini_CV)
    
    
    
    curve.compute_area()
    curve.MomentMatrices()
    curve.compute_moments()
    curve.computeXc()
    
    print '\n\nArea interval Finding....'
    print 'area: ',curve.area
    
    
    
    
    print '\nXC interval Finding....'
    
    print 'using x3=6.+{} and actual area:'.format(xshift)
    eXc = ini_CV.calc_xc(6.+xshift,curve.area)
    print 'est Xc = ',eXc
    print 'actual Xc = ', curve.Xc
    print 'diff = ',curve.Xc-eXc
    
    
    
    #"""
    #self = ini_CV
    xcini = ia(5.,9.)
    area = ia(71.,73.)
    print '\nNow use ...'
    print 'area = 72., ',xcini,' Xc = ',ini_CV.calc_xc(xcini,curve.area)
    print 'area = IA(), ',xcini,' Xc = ',ini_CV.calc_xc(ia(0.,12.),area)
    
    area = ia(50.,90.)
    print 'area = IA(big), ',xcini,' Xc = ',ini_CV.calc_xc(xcini,area)
    #"""
    print '\n\n'
    
    return curve, ini_CV

if __name__ == """__main__""":
    from curve import Bspline
    # We will have a 7 control point curve, Order 4, as recommended by Harries, 1998.
    # Given The Following Form Parameters:



    curve, ini_CV = test()
    
    #ivar = var_to_ad_ia(1.,12.,N=2,dim=0,name='x')
    ivar = var_to_ad_ia(2.,4.,N=2,dim=0,name='x')
    x3 = ivar
    #area = ia(50.,90.)
    area = var_to_ad_ia(50.,90.,N=2,dim=1,name='area')
    
    funclos = F_Xc(ini_CV,area)
    ff = funclos(ivar,area)
    """
    func=funclos
    c=ia(3.,9.)
    x=ivar
    vertex_tpl=[0,0]
    itmax = 5
    tol = 1.e-2
    """
#    ans = montonic_compute_scalar_newton(func=funclos,
#                                         c=ia(3.,9.),
#                                         x=ivar, 
#                                         vertex_tpl=[0,0], 
#                                         itmax = 5, 
#                                         tol = 1.e-2)
    
    #"""
    curve, ini_CV = test(xshift=1.)
    
    
    curve, ini_CV = test(xshift =-12.,
                         xb     = 12.,
                         xe     = 24.)
    
    

    xb=0.               #1st vertex: x value
    yb=0.               #1st vertex: y value
    xe=12.              #last vertex: x value
    ye=12.              #last vertex: y value
    
    curve,ini_CV = test(xb=0.,
                         yb = 0.,
                         xe = 12.,
                         ye = 12.,
                         area = 72.)
            
    
    
    #"""
    
    
    """
    self = ini_CV


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
    
    x3 = ia(0.,12.)
    
    area=72.

    """
    
    
    def test_HC_BC_2(V):
        """
            constraint:
            f(x,y) = x^3 + 100x + 10y = 0
            
            
        #"""
        x = V[0]
        y = V[1]
        return 10.*y + x**3 + 100.*x
    
    
    xo = ad( ia(-100., 100.), name = 'x', N=2, dim=0)
    yo = ad( ia(-100., 100.), name = 'y', N=2, dim=1)
    X = [xo,yo]
    g = test_HC_BC_2(X)
    print g.value
    print ''
    print g.grad
    print ''
    print g.hess
    print ''
    print 10.*yo + xo**3 + 100.*xo
    print ''
    print 3.*(xo**2) + 100.
    print 3.*2.*(xo)
    
    Fthis = F_Xc(ini_CV,
         area   = ia(50.,90.))
    
    
    xc = ini_CV.calc_xc(area   = ia(50.,90.),
                                   x3     = ia(0.,12.))