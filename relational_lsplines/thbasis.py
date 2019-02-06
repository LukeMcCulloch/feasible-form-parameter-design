# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 09:31:55 2016

@author: Luke.McCulloch

efficient every other array scatter slicing:
http://stackoverflow.com/questions/25876640/subsampling-every-nth-entry-in-a-numpy-array


Sources, verified:
    
    (1)
    Ma, Cripps, Li, Lin:
        A subdivision-based implementation of non-uniform
        local refinement with THB-splines
        
    (2)
    Sabin, Dodgson, Cashman:
        A symmetric, non-uniform, refine and smooth subdivision algorithm for general degreeB-spline3
        
    (3)
    Schaefer, Goldman:
        Non-uniform Subdivision for B-splines of Arbitrary Degree
    

Sources?

Sabin, Romani:
    The conversion matrix between uniformB-spline and Bezier representations
    
Juttler, Speleers, Giannelli:
    Strongly Stable bases for adaptively refined multilevel spline spaces
    
Jansen:
    Non-equispaced B-spline wavelets

Juttler, Giannelli, Kiss: *
    Algorithms and Data Structures for truncated hierarchical b-splines
    

Giannelli, Juttler, Spelleers:
    THB-splines:  The truncated basis for hierarchical splines    


Giannelli, Kiss, Zore, Juttler, Grobmann, Barner: *
    Adaptive CAD model (re-)construction with THB-splines
    
Bornemann, Cirak:
    A subdivision-based implementation of the hierarchical b-spline finite element method
    
Buffa, Eduardo, Giannelli, Sangalli:
    On quasi-interpolation operators in spline spaces
    
Bertram:
    Single-knot wavelets for non-uniform B-splines 
"""
import numpy as np
import curve as spline
from HierarchicalBspline import LevelSpline as lespline
from nbspbasis import TwoScale
import time


def knot_union(tau,tp):
    return np.sort(np.insert(tau,-1,tp),kind='quicksort')

def knot_union2(tau,tp):
    taulist = list(tau)
    tplist = []
    maxtau = taulist[-1]
    mintau = taulist[0]
    for el in tp:
        if mintau<el and el<maxtau:
            tplist.append(el)
    tp = np.asarray(tplist)
    return np.sort(np.insert(tau,-1,tp),kind='quicksort')
      
    
def affine1(a,b,g,e):
    if a <= b and b<=g and a != e:
        return (g-b)/(e-a)
    else:
        return 0.0
        
def affine2(a,b,e):
    if a <= b and b <= e and a != e:
        return (b-a)/(e-a)
    else:
        return 0.0

def affine3(a,g,e):
    if a <= g and g <= e and a!= e:
        return (e-g)/(e-a)
    else:
        return 1.0


def subdivide_basis(d,P,tau,tp):
    """
    Sabin algorithm
        d   : degree
        P   : control points
        tau : initial knot vector
        tp  : knots to insert
    output
        Pid : refined conrol points
    """
    n = len(P)-1
    dim = len(P[0])
    Pm = np.zeros((d,2*(n+1),dim),float)
    u = knot_union(tau,tp)
    if d%2==1:
        #Refine
        Pm[1,0::2,:] = P[:]
        Pm[1,-2,:] = 0.
        #for i in range(0,n):
        #    Pm[1,2*i] = P[i]
        for i in range(0,n): #altered from n-1 => changed to get 00 on the end
            j = 2*i + 1
            Pm[1,j] = affine3(u[j-d],u[j],u[j+d])*P[i] + affine2(u[j-d],u[j],u[j+d])*P[i+1]
        #Smooth
        for lm in range(1,((d-1)/2)):
            for i in range(lm,2*n-lm):
                if lm%2 == i%2:
                    Pm[2*lm+1,i] = Pm[2*lm-1,i]
                else:
                    c0 = affine3(u[i-d+lm]        ,u[i+lm],u[i+d-lm]) * Pm[2*lm-1,i-1]
                    c1 = affine1(u[i-d+lm],u[i-lm],u[i+lm],u[i+d-lm]) * Pm[2*lm-1,i  ]
                    c2 = affine2(u[i-d+lm],u[i-lm]        ,u[i+d-lm]) * Pm[2*lm-1,i+1]
                    #
                    #co = (u[i+d-lm]-u[i+lm])/(u[i+d-lm]-u[i-d+lm])
                    #c1 = (u[i+lm]-u[i-lm])/(u[i+d-lm]-u[i-d+lm])
                    #c2 = (u[i-lm]-u[i-d+lm])/(u[i+d-lm]-u[i-d+lm])
                    #
                    Pm[2*lm+1,i] = c0+c1+c2
    else:
        print 'even degree not implemented'
    return Pm,u

def subdivide_basis1(d,P,tau,tp):
    """
    Sabin algorithm
        d   : degree
        P   : control points
        tau : initial knot vector
        tp  : knots to insert
    output
        Pid : refined conrol points
    """
    n = len(P)-1
    dim = len(P[0])
    Pm = np.zeros((d,2*(n+1),dim),float)
    u = knot_union(tau,tp)
    assert( d%2==1)#print  'spline must be odd degree'
    #Refine
    Pm[1,0::2,:] = P[:]
    #Pm[1,-2,:] = 0.  #n => changed to get 00 on the end
    #for i in range(0,n):
    #    Pm[1,2*i] = P[i]
    for i in range(0,n): #n-1
        j = 2*i + 1
        Pm[1,j] = affine3(u[j-d],u[j],u[j+d])*P[i] + affine2(u[j-d],u[j],u[j+d])*P[i+1]
    #Smooth
    for lm in range(1,((d)/2)): #added 1 to (d-1)
        for i in range(lm-1,2*n-lm+1):
            if lm%2 == i%2:
                Pm[2*lm+1,i] = Pm[2*lm-1,i]
            else:
                c0 = affine3(u[i-d+lm]        ,u[i+lm],u[i+d-lm]) * Pm[2*lm-1,i-1]
                c1 = affine1(u[i-d+lm],u[i-lm],u[i+lm],u[i+d-lm]) * Pm[2*lm-1,i  ]
                c2 = affine2(u[i-d+lm],u[i-lm]        ,u[i+d-lm]) * Pm[2*lm-1,i+1]
                Pm[2*lm+1,i] = c0+c1+c2
    return Pm,u



        
def subdivide_basis_not_SABIN(p,tau,tp):
    """
    non Sabin algorithm  (source: (1))
        p   : degree (odd? TLM June 2017)
        tau : initial knot vector active for the particular 
                basis function being refined
        tp  : knots to insert
        
    from myspline:
        p = self.p
        tau = tau
        tp = new_knots
    """
    #t = knot_union2(tau,tp)
    t = knot_union(tau,tp)
    
    W = np.zeros((len(tp)),float)
    
    w = np.zeros((len(tp),p),float)  #assign initial weight
    w[:,0] = 1.0
    

    for lm in range(0,(p+1)/2):
        print 'lm =', lm
        #for i in range(p-lm,p+lm+2):
        #for i in range(p-lm,p+lm+2):
        for i in range(p-lm-1,p+lm+2-1):
            j = i-(p+1)/2
            print 'i,j = ', i,j
            if lm%2 == i%2:
                w[j,lm+1] = w[j,lm]
            else:
                check = (t[i+p-lm]-t[i-p+lm])
                if check==0.:
                    c0 = 0.
                else:
                    c0 = (t[i+p-lm]-t[i+lm])/(t[i+p-lm]-t[i-p+lm])
                check = (t[i+p-lm]-t[i-p+lm])
                if check == 0.:
                    c1 = 0.
                else:
                    c1 = (t[i+lm]-t[i-lm])/(t[i+p-lm]-t[i-p+lm])
                
                check = (t[i+p-lm]-t[i-p+lm])
                if check == 0.:
                    c2 = 0.
                else:
                    c2 = (t[i-lm]-t[i-p+lm])/(t[i+p-lm]-t[i-p+lm])
                    
                if i==(p-lm):
                    w[j,lm+1] = c2*w[j+1,lm]
                elif i==(p+lm+2):
                    w[j,lm+1] = c0*w[j-1,lm]
                else:
                    w[j,lm+1] = c0*w[j-1,lm] + c1*w[j,lm] + c2*w[j+1,lm]
            #w0[:] = w1[:]
            #w1[:] = 0.
    for j in range(0,p+1):
        W[j] = w[j,(p+1)/2]
    return W

def subdivide_basis_not_SABIN__old(p,tau,tp):
    """
    non Sabin algorithm  (source: (1))
        p   : degree (odd? TLM June 2017)
        tau : initial knot vector active for the particular 
                basis function being refined
        tp  : knots to insert
        
    from myspline:
        p = self.p
        tau = tau
        tp = new_knots
    """
    t = knot_union2(tau,tp)
    
    W = np.zeros((len(tp)),float)
    
    w = np.zeros((len(tp),p),float)  #assign initial weight
    w[:,0] = 1.0
    

    for lm in range(0,(p+1)/2):
        print 'lm =', lm
        #for i in range(p-lm,p+lm+2):
        for i in range(p-lm,p+lm+2):
            j = i-(p+1)/2
            print 'i,j = ', i,j
            if lm%2 == i%2:
                w[j,lm+1] = w[j,lm]
            else:
                c0 = (t[i+p-lm]-t[i+lm])/(t[i+p-lm]-t[i-p+lm])
                c1 = (t[i+lm]-t[i-lm])/(t[i+p-lm]-t[i-p+lm])
                c2 = (t[i-lm]-t[i-p+lm])/(t[i+p-lm]-t[i-p+lm])
                if i==(p-lm):
                    w[j,lm+1] = c2*w[j+1,lm]
                elif i==(p+lm+2):
                    w[j,lm+1] = c0*w[j-1,lm]
                else:
                    w[j,lm+1] = c0*w[j-1,lm] + c1*w[j,lm] + c2*w[j+1,lm]
            #w0[:] = w1[:]
            #w1[:] = 0.
    for j in range(0,p+1):
        W[j] = w[j,(p+1)/2]
    return W
    
def test_algo():
    k = 4
    p = k-1
    nump = 60
    correction = -1
    cfac = correction
    verts = spline.linear_vertices([0.,0.],[12.,12.],7)
    verts[:,1] = np.cos(verts[:,1])
    self = spline.Bspline(verts,k,nump)
    lself = lespline(verts,k,nump)
    p = self.p
    TAU = self.t
    tpi = np.sort(list(set(TAU)))
    TP = np.zeros((p+2+cfac),float)
    for i in range(len(TP)):
        TP[i] = (tpi[i]+tpi[i+1])*.5
        
    ith = 0
    actv_knots = self.active_knots(ith)
    #tp_insert = knot_union2(actv_knots,TP)  #intersection instead??
    tau = actv_knots
    tp = TP#tp_insert
    W = subdivide_basis_not_SABIN(p,tau,tp)
    
    for i in range(0,self.n-1):
        tau = self.active_knots(i)
        print  subdivide_basis_not_SABIN(p,tau,tp)
    return

def subdivide_uniform(p):
    """
        cite?
    """
    W = np.zeros((p+2,1),float)
    i=0
    while (i<p+2):
        W[i] = TwoScale(i,p+1)
        i+=1
    return W
    
def subdivide_uniform_extra(p):
    """
        mod to match ~k~
    """
    #p=p-1
    W = np.zeros((p+2,1),float)
    i=0
    while (i<p+2):
        W[i] = TwoScale(i,p+1)
        i+=1
    return W

def factorial (n):
    if not isinstance(n, int):
        print 'Factorial is only defined for integers.'
        return None
    elif n < 0:
        print 'Factorial is not defined for negative integers.'
        return None
    elif n == 0:
        return 1
    else:
        return n * factorial(n-1)

def hbj(p,j):
    return factorial(p+1)/(factorial(j)*factorial(p+1-j))

def emulate_basis_old(i,pt,S,coarse,fine):
    """
        i : ith basis function
        S : Two Scale Subdivision Matrix
        B_fine : Fine (L+1) level Basis functions genertor
    """
    import scipy.special as sp
    binom = sp.binom
    
    p = coarse.p
    k = p+1
    B = np.zeros((p+2,1),float)  #p subdivision
    #B = np.zeros((p+1,1),float) #p-1 subdivision
    p = coarse.p
    #for j in range(i,i+k):
    span_fine = fine.FindSpan(pt)
    #B[0:k,0]  = fine.eval_BasisPt(2*i+k,pt).T
    B[0:k,0]  = fine.eval_BasisPt(span_fine,pt).T
    lenB = len(B)
    Bs = B.T.dot(S)[0][0]

      
    
    oB = 0.
    for j in range(lenB):
        oB += B[i]*S[i]
    
    span_coarse = coarse.FindSpan(pt)
    #Bc = coarse.eval_BasisPt(i+k,pt)
    Bc = coarse.eval_BasisPt(span_coarse,pt)
    
    
    ###
    ###***experiments
    ###
    #pt = .8
    Bf = fine.TLMBasisFuns(pt)
    Bc = coarse.TLMBasisFuns(pt)
    span_fine = fine.FindSpan(pt)
    span_coarse = coarse.FindSpan(pt)
    
    #i=3
    ans1 = 0.
    ans2 = 0.
    for j in range(0,k+1):
        print binom(p+1,j)
        print S[j]
        print Bf[span_fine-p+j]
        ans1 += binom(p+1,j)*Bf[span_fine-p+j]
        ans2 += S[j]*Bf[span_fine-p+j]
    ans1 = ans1/(2**(p))
    ans2 = ans2/(2**(p))
    
    print 'coarse func = ',Bc[span_coarse-p+i]
    ans3=0.
    ar = range(0,p+2)
    br = range(0,p+2)
    br.reverse()
    for j,rj in zip(ar,br):
        print j, rj
        #print Bf[span_fine-p+rj]
        print hbj(p,j)
        print hbj(p,j)/(2.**p)
        ans3 += hbj(p,j)*Bf[span_fine-rj]
    ans3 = ans3/(2**p)
    print 'fine func = ',ans3
    
    ans4 = 0.0
    for j in range(p+1):
        ans4 += Bc[span_coarse-p+i]*Bf[span_fine-p+j]
    
    
    #try to match... 
    print Bc[span_coarse-p+i], ans4
    return ans1

def emulate_basis(i,pt,S,coarse,fine):
    """
        i : ith basis function
        S : Two Scale Subdivision Matrix
        B_fine : Fine (L+1) level Basis functions genertor
    """
    import scipy.special as sp
    binom = sp.binom
    
    p = coarse.p
    k = p+1
    
    #pt = .8
    Bf = fine.TLMBasisFuns(pt)
    Bc = coarse.TLMBasisFuns(pt)
    span_fine = fine.FindSpan(pt)
    span_coarse = coarse.FindSpan(pt)
    
    ans = 0.0
    for j in range(p+1):
        ans += Bc[span_coarse-p+i]*Bf[span_fine-p+j]
    
    print Bc[span_coarse-p+i], ans
    return ans

def check_emulated_basis(self,lself):
    """
        Use the coarsesubdivision 
        to compute  
        the fine subdivision
    """
    S = subdivide_uniform(p)
    param = np.linspace(0.,1.,10,endpoint=True)
    for pt in param:
        for i in range(4):
            coarse = self
            fine = lself
            B = emulate_basis(i,pt,S,coarse,fine)
    return


def EvalBlossom(d,k,u,t):
    """siggraph course 1996
    intro to curves and surfaces page 63
    
    not setup correctly for P&T B-splines?
    """
    b = np.zeros(d+1)
    for k in range(d):
        for j in range(d-k-1):
            Beta = (u[k + 1] - t[i + k + j]) / (t[i + d + j] - t[i + k + j])
            Alpha = 1 - Beta
            b[j] = Alpha * b[j] + Beta * b[j + 1]
    return b
    
    
if __name__ == '__main__':
    green = 'green'
    red = 'red'
    blue = 'blue'
    k = 4
    nump = 60
    correction = -1
    cfac = correction
    verts = spline.linear_vertices([0.,0.],[12.,12.],7)
    verts[:,1] = np.cos(verts[:,1])
    self = spline.Bspline(verts,k,nump)
    lself = lespline(verts,k,nump)
    
    time0=time.clock()
    lself = lself.dyadic_refinement()
    time1=time.clock()
    print 'refined Basis calc time method 1  = {} s'.format(time1-time0)

    #self.plotcurve_detailed()
    
    p = self.p
    tau = self.t
    tpi = np.sort(list(set(tau)))
    tp = np.zeros((p+2+cfac),float)
    for i in range(len(tp)):
        tp[i] = (tpi[i]+tpi[i+1])*.5

    d = self.p
    P=verts
    
    Pl,ul = subdivide_basis1(d,P,tau,tp)
    v2 = spline.Bspline(Pl[1,1:-2],k,nump,t=ul)
    time0=time.clock()
    Pm,um = subdivide_basis(d,P,tau,tp)    
    time1=time.clock()
    print 'refined Basis calc time method 2  = {} s'.format(time1-time0)
    print 'but method 2 remains incorrect!'
    v1 = spline.Bspline(Pm[1,1:-2],k,nump,t=um)
    
    nump=100
    interval = [0.,.24]
    offset = 0.
    bs = self.plot_basis_interval(interval,
                                  nump)
    interval = [.76,1.]
    offset = 0.
    bs = self.plot_basis_interval(interval,
                                  nump)
    ## level 2               
    interval = [.24,.3]
    offset += 1.5
    ls = lself.plot_basis_interval(interval,
                                  nump,
                                  offset)
    interval = [.7,.76]
    ls = lself.plot_basis_interval(interval,
                                  nump,
                                  offset)
    
    #level 3        
    interval = [.3,.42]            
    offset += 1.5
    llself = lself.dyadic_refinement()
    lls = llself.plot_basis_interval(interval,
                                  nump,
                                  offset) 
    interval = [.58,.7]
    llself = lself.dyadic_refinement()
    lls = llself.plot_basis_interval(interval,
                                  nump,
                                  offset)
    #level 4      
    interval = [.42,.58] 
    offset += 1.5
    lllself = llself.dyadic_refinement()
    llls = lllself.plot_basis_interval(interval,
                                  nump,
                                  offset)
     
    
    
    verts = spline.linear_vertices([0.,0.],[12.,12.],4)
    verts[:,1] = verts[:,1]*np.cos(verts[:,1])
    self = spline.Bspline(verts,k,nump)
    self = lespline(verts,k,nump)
    
    s2 = self.dyadic_refinement()
    s3 = self.dyadic_refinement().dyadic_refinement()
    s4 = self.dyadic_refinement().dyadic_refinement().dyadic_refinement()
    
    f1 = .25/2.
    
    tnew = np.asarray([f1,.25+f1,.5+f1,.75+f1 ])
    tblossom = np.asarray([0.,1.,1.,1.,1.,0.])  
    
    v = subdivide_basis_not_SABIN(self.p,
                                  self.t,
                                  np.asarray([0.,1.,0.,1.,0.,]) 
                                  )
    