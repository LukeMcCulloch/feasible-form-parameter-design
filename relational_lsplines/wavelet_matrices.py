#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 15:02:35 2017

@author: TLM

nullspace of a matrix:
    http://scipy-cookbook.readthedocs.io/items/RankNullspace.html
  is a better answer to:
    https://stackoverflow.com/questions/5889142/python-numpy-scipy-finding-the-null-space-of-a-matrix
"""

import numpy as np
import scipy as sp
import scipy.misc
#import sympy as sy #LCD used for refinement.. add cite here.

import copy

import HierarchicalBspline as hb
import thbasis as thb
from curve import linear_vertices
#factorial = sp.math.factorial
factorial = sp.misc.factorial
rank = np.linalg.matrix_rank

try:
    import sympy as sy
    rat = sy.Rational
    frac = sy.fraction
except:
    print 'WARNING: sympy not found.. '
    print 'bypassing some simple number theory (LCD) related wavelet stuff...'
    
    def rationalUNO():
        return 'not implemented at UNO'
    def fractionUNO():
        return 'not implemented at UNO'
    rat = rationalUNO
    frac = fractionUNO

svd = np.linalg.svd


#
#***********************************************************
#
def Factorial(m):
    try:
        if len(m) > 1:
            tisarray = 1
        else:
            tisarray = 0   
    except:
        tisarray = 0
        
    if tisarray == 1:
        r,c = np.shape(m)
        f = np.zeros((r,c))
        for i in range(0,r):
            for j in range(0,c):
                f[i,j] = np.prod( np.arange(2,m[i,j]+1) )
        return f
    else:
        if m<0.:
            return -1.
        else:
            return factorial(m)

#
#***********************************************************
#              

def matlab_arii(i,e):
    return np.arange(i,e+1)

#
#***********************************************************
#
def Knots(d,j):
    """
        x = Knots(d, j) returns a vector 
        of knot values for B-spline scaling
        functions of degree d, level j.
    """
    aa = matlab_arii(0.,2.**j-1)/(2.**j)
    x = np.asarray([0. for el in range(d-1)] + list(aa) + [1. for el in range(d)])
    return x
    
#
#***********************************************************
#

def Greville(d,u):
    """verified
        x = Greville(d, u) 
        returns the vector of Greville abscissa values
        corresponding to degree d and knot vector u.
    """
    l = len(u)
    x = u[0:1-d]
    for k in range(2,d+1):
        x = x + u[k-1:l-d+k]
    return x / d

#
#***********************************************************
#
#def Choose(n,r):
#    return Factorial(n) / (Factorial(r) * Factorial(n-r))
def Choose(i,d):
    return np.divide( Factorial(i) , np.multiply( Factorial(d) , Factorial(i-d) ) )
#def hbj(p,j):
#    return factorial(p+1)/(factorial(j)*factorial(p+1-j))
#def Choose(p,j):
#    return factorial(p+1)/(factorial(j)*factorial(p+1-j))

def BernsteinInner(d):
    """
        I = BernsteinInner(d) returns the 
        matrix of inner products of Bernstein
        polynomials of degree d.
    """
    i = np.ones((d+1, 1),int)*np.arange(0,d+1,1,int)
    j = i.T
    I = np.divide( np.multiply( Choose(d, i) , Choose(d, j) ) , (Choose(2*d, i+j)*(2*d + 1)) )
    return I

def BernsteinWeights(d,j):
    w = np.identity(2**j + d)
    if d==0:
        return w
    u = Knots(d,j)
    g = Greville(d,u)
    
    for i in range(0,2**j-1):
        for r in range(0,d):
            u,g,w = InsertKnot(d,u,g,w,(i+1.)/(2.**j))
    return w

def Inner(d,j):
    I0 = BernsteinInner(d)
    n = 2**j + d
    I = np.zeros((n,n))
    w = BernsteinWeights(d,j)
    for k in range(0,n):
        #w1 = np.reshape(w[:,k],d+1,2**j)
        #w1 = w[:,k].reshape(d+1,2**j)
        w1 = w[:,k].reshape(2**j,-1)
        w1 = w1.reshape(2**j,d+1).T
        #setshape??
        for l in range(0,n):
            #w2 = w[:,l].reshape(d+1,2**j)
            
            w2 = w[:,l].reshape(2**j,-1)
            w2 = w2.reshape(2**j,d+1).T
            
            #I[k,l] = np.trace(w1.T*I0*w2)
            I[k,l] = np.matmul(np.matmul(w1.T,I0),w2).trace()
            I[l,k] = I[k,l]
    I = I / 2.**j
    return I

###----------------------------------------------------------------------------
### 2 Scale Relation for Basis Functions  TLMTLMTLMTLM below!
###----------------------------------------------------------------------------

"""
NOT USED!
#import scipy as sp
"""

def binomial(n_,i_):
    """
        P&T : i is scipy k
        (n,i) <=> (n,k) so i<=>k in the literature
        where
            n is `top`
            i is `bottom`
        
    """
    return sp.special.binom( n_,i_)
    
def TwoScale(i,k):
    return binomial(k,i)/(2.**(k-1))

def checkstein(i,n,u):
    return binomial(n,i)*(u**i)*((1.-u)**(n-i))

def Bernstein(i,n,u):
    """return the ith (n+1)th order, (n)th degree Bernstein polynomial
        only if nonzero at fixed u
        i : ith conrol point goes with ith Basis function (ith span index?)
        n : degree
        u : knot location
        B : Polynomial value
    Piegel and Tiller, page 20
    """
    K = n+1
    B      = np.zeros((K),float)
    B[n-i] = 1.0
    u1     = 1.0-u
    for k in range(1,K):
        for j in range(n,k-1,-1): #careful - sly index!
            B[j] = u1*B[j] + u*B[j-1]
    return B[n]

def AllBernstein(n,u):
    """return all of the ith (n+1)th degree Bernstein polynomial
        only compute if nonzero at fixed u
        n : degree
        u : knot location
        B : array of 
    Piegel and Tiller, page 21
    """
    K       = n+1
    B       = np.zeros((K),float)
    B[0]    = 1.0
    u1      = 1.0-u
    for j in range(1,K):
        saved = 0.
        for k in range(0,j):
            temp = B[k]
            B[k] = saved + u1*temp
            saved = u*temp
        B[j] = saved
    return B

#
#***********************************************************
#

def PolyEval(g, p, gnew):
    """
    % pret = PolyEval(g, p, gnew) returns the values of a control polygon
    % defined by abscissas g and ordinates p, evaluated at gnew.
    """
    m, n = np.shape(p)
    assert(np.size(g) == m),'PolyEval: Length of g and rows of p must be the same.'
    
    lgn = np.size(gnew)
    pret = np.zeros((lgn,n)) #TLM GUESS!
    #
    #*******
    #COMMON MISTAKE:
    #
    #for i in range(1,len(gnew)+1):
    #
    #****
    #Correction:
    #
    for i in range(0,len(gnew)):
      #row = max(find(g <= gnew(i)))
      row = max( (g <= gnew[i]).nonzero() )[-1]
      #row = (g <= gnew[i]).nonzero().max()
      #aaa = g <= gnew[i]
      #ara = aaa.nonzero()
      # when i=0 row = 1
      # corresponds to matlab
      # when i=1, row=1
      if row == m-1:
        pret[i,:] = p[m-1,:]
      else:
        frac = (g[row+1] - gnew[i])/(g[row+1] - g[row])
        pret[i,:] = frac*p[row,:] + (1 - frac)*p[row+1,:]
    return pret

#
#***********************************************************
#

def InsertKnot(d, u, g, p, unew):
    """
    % [uret, gret, pret] = InsertKnot(d, u, g, p, unew) inserts a new knot at
    % unew for B-spline scaling functions of degree d, thereby modifying knot
    % vector u, Greville abscissas g, and synthesis matrix p.
    uret : modified knot vector
    gret : modified Greville abscissas
    pret : modified control synthesis matrix, P
    """
    uret = np.sort(np.concatenate( (u , [unew]),axis=0))
    gret = Greville(d, uret)
    pret = PolyEval(g, p, gret)
    return uret, gret, pret
#
#***********************************************************
#
"""
NOT USED!
"""
#def matlab_cat1_old(Vec, Array):
#    Array[:,0] = Array[:,0]+Vec
#    return Array
#def matlab_cat2_old(Array, Vec):
#    Array[:,1] = Array[:,1]+Vec
#    return Array
#
#def matlab_cat1(Vec, Array):
#    return np.concatenate((Vec,Array), axis=1)
#def matlab_cat2(Array, Vec):
#    return np.concatenate((Array,Vec), axis=1)
#
#***********************************************************
#
def FindP(d,j):
    """
        returns the P matrix for B-spline scaling functions
        given a degree and dyadic refinement level
        
    parameters
    ----------
        degree  : d 
        level   : j
    """
    d = int(np.fix(d))
    assert(d>=0),'Error, FindP:  Must have d >= 0.'
    #assert(j >= 1),'Error, FindP:  Must have j >= 1.'
    if d == 0:
        #P = np.asarray([[1],[1]])
        P = np.asarray([[1.,1.]]).T
        
        #mm = np.zeros(shape = [2**j,2**(j-1)])
        for i in range(2,j+1):
            print i
            sp = np.shape(P)
            p1 = np.concatenate((P,np.zeros(sp)), axis=1)
            p2 = np.concatenate((np.zeros(sp),P), axis=1)
            P = np.array([list(p1) + list(p2)])[0]
            
            """
            P = np.array([ [ list(el) for el in p1]+
                            [list(el) for el in p2] ])
            #"""
    else:
        u = Knots(d,j-1)
        g = Greville(d,u)
        P = np.identity(2**(j-1) + d)
        for k in range(0, 2**(j-1)-1+1  ):
            u,g,P = InsertKnot(d, u, g, P, (2*k+1.)/2**j  )
    
    return P
    
def FindP_one(d,j,b):
    """
        returns the P matrix for B-spline scaling functions of 
        degree      : d 
        level       : j
        basis index : b
    """
    d = int(np.fix(d))
    assert(d>=0),'Error, FindP:  Must have d >= 0.'
    assert(j >= 1),'Error, FindP:  Must have j >= 1.'
    if d == 0:
        #P = np.asarray([[1],[1]])
        P = np.asarray([[1.,1.]]).T
        
        #mm = np.zeros(shape = [2**j,2**(j-1)])
        for i in range(2,j+1):
            print i
            sp = np.shape(P)
            p1 = np.concatenate((P,np.zeros(sp)), axis=1)
            p2 = np.concatenate((np.zeros(sp),P), axis=1)
            P = np.array([list(p1) + list(p2)])[0]
            
            """
            P = np.array([ [ list(el) for el in p1]+
                            [list(el) for el in p2] ])
            #"""
    else:
        u = Knots(d,j-1)
        g = Greville(d,u)
        P = np.identity(2**(j-1) + d)
        #for k in range(0, 2**(j-1)-1+1  ):
        k = b
        u,g,P = InsertKnot(d, u, g, P, (2*k+1.)/2**j  )
    
    return P[b]




#
#***********************************************************
#

#def null(A, eps=1e-15):
#    u, s, vh = sp.linalg.svd(A)
#    null_mask = (s <= eps)
#    null_space = sp.compress(null_mask, vh, axis=0)
#    return sp.transpose(null_space)

def nullspace(A, atol=1e-13, rtol=0):
    A = np.atleast_2d(A)
    u, s, vh = svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns
null = nullspace

def gcd(*numbers):
    """Return the greatest common divisor of the given integers"""
    from fractions import gcd
    return reduce(gcd, numbers)

def lcm(*numbers):
    """Return lowest common multiple."""    
    def lcm(a, b):
        return (a * b) // gcd(a, b)
    return reduce(lcm, numbers, 1)

#
#***********************************************************
#
vfrac = np.vectorize(frac)
vrat = np.vectorize(rat)
def LCD(m):
    num, denom = vrat(m)
    return d
#LCD = lcm
#
#***********************************************************
#
normalization='min'#'L2'
def FindQ(d, j, normalization='min'):
    """
        returns the Q matrix for B-spline wavelet functions
        given a degree and dyadic refinement level
        
    parameters
    ----------
        degree  : d 
        level   : j
    """
    P = FindP(d,j)
    I = Inner(d,j)
    #
    M = np.matmul(P.T,I)
    m1,m2 = np.shape(M)
    n = m2 - rank(M) #M.ndim #np.rank(M)
    Q = np.zeros((m2,n))
    found = 0
    start_col = 0
    while ( (found < n/2.) and (start_col < m2) ):
        #beware the matlab indices!  (used verbatum here)
        start_col =  start_col + 1 + int(found > d)
        width = 0
        rank_def = 0
        #while(  rank_def == 0 and (width < m2 - start_col +1) ):
        while(  not rank_def and (width < m2 - start_col +1) ):
            width = width + 1
            #submatrix = M[:,start_col:start_col+width-1]
            submatrix = M[:,start_col-1:start_col+width-1] #adjusted indices here!
            rank_def = width - rank(submatrix)
        #if rank_def != 0:
        if rank_def:
            #print 'width = ', width
            q_col = null(submatrix)
            #--------------------------------------------------------
            if normalization == 'min':
                q_col = q_col/min(abs(q_col + 1e38*(abs(q_col) < 1e-10)))
            elif normalization == 'max':
                q_col = q_col/max(abs(q_col))
            elif normalization == 'lcd':
                print 'error LCD not implemented yet'
                assert(normalization != 'lcd'),'Error, LDC norm not implemented... sympy?'
                pass
                q_col = q_col/min(abs(q_col + 1e38*(abs(q_col) < 1e-10)))
                q_col = q_col*LCD(q_col)
            #--------------------------------------------------------
            # change sign to give consistent orientation
            q_col = q_col*(-1)**(start_col + np.floor((d+1.)/2.) + (q_col[0,0] > 0))
            #print q_col
            # correct any slight error for answers that should be integers
            #if np.all(abs(submatrix*np.round(q_col)) < 1e-10) and np.any(np.round(q_col) != 0):
            if np.all(abs(np.matmul(submatrix,np.round(q_col)) ) < 1e-10) and np.any(np.round(q_col) != 0):
                q_col = np.round(q_col)
        # put column into strcmpleft half of Q
        found = found + 1
        #Q[start_col:start_col+width-1+1,found] = q_col[:,0]
        Q[start_col-1:start_col-1+width-1+1,found-1] = q_col[:,0]
        #Q[start_col-1:start_col-1+width-1+1,found] = q_col[:,0]
        
        # use symmetry to put column into right half of Q in reverse order
        # and negated if degree is even
        Q[:,n-found] = np.flipud(Q[:,found-1])*(-1.)**(d+1.)
        #print rank_def,  (width < m2 - start_col +1) ,width, m2,start_col
    if normalization=='L2':
        ip = np.matmul(Q.T,np.matmul(I,Q))
        Q = np.matmul(Q,np.diag(1./np.sqrt(np.diag(ip))))
    return Q

#
#***********************************************************
#
def LUPQI(P,Q,vertices):
    """Perform LU decomposition
    and solve for the restriction of
    a vector
    
    in linear time.
    
    Q = c0c.rmi
    P = c0c.rm
    """
    Qshape = np.shape(Q)
    Pshape = np.shape(P)
    ncol = Qshape[1] + Pshape[1] 
    PQM1 = np.zeros((Pshape[0],ncol))
    PQM1[:,:Pshape[1]] = P
    PQM1[:,Pshape[1]:] = Q
    #AB1 = np.linalg.inv(PQM1)
    
    
    sp.linalg.lu(PQM1)
    return

#
#***********************************************************
#
if __name__ == '__main__':
    u = np.asarray([1,1,1,0,0,0,1,1,1],float)
    print np.matrix(Greville(3,u)).T
    """
        matlab:
            ans =

                   1.00000
                   0.66667
                   0.33333
                   0.00000
                   0.33333
                   0.66667
                   1.00000
    """
    u = np.asarray([1,1,1,1,0,0,0,1,1,1,1],float)
    print np.matrix(Greville(4,u)).T
    """
        matlab:
            ans =

                   1.00000
                   0.75000
                   0.50000
                   0.25000
                   0.25000
                   0.50000
                   0.75000
                   1.00000
    """
    
    u = np.asarray([1.,2.,3.,4.,5.,6.])
    
    #degree:
    d = 2
    
    #level:
    j = 2
    
    t = Knots(d,j)
    """
        >> t = Knots(d,j);
        >> t
        t =
        
           0.00000   0.00000   0.25000   0.50000   0.75000   1.00000   1.00000
    #"""
    
    """
            >>
        >> Knots(d, j - 1)
        ans =
        
           0.00000   0.00000   0.50000   1.00000   1.00000
           >>> Knots(d, j - 1)
           array([ 0. ,  0. ,  0.5,  1. ,  1. ])
        
        >>
        >>
        >> Greville(d, u)
        ans =
        
           1.5000   2.5000   3.5000   4.5000   5.5000
           
          >>> Greville(d, u)
          array([ 1.5,  2.5,  3.5,  4.5,  5.5])

    #"""
    
#    u = Knots(d,j-1)
#    g = Greville(d,u)
#    P = np.identity(2**(j-1) + d)
#    k=0 #k=1
#    p = P
#    unew = (2*k+1.)/(2**j)
#    uret = np.sort(np.concatenate( (u , [unew]),axis=0))
#    gret = Greville(d, uret)
#    
#    gnew = gret
#    
#    m, n = np.shape(p)
#    
#    
#    u,g,P = InsertKnot(d, u, g, P, (2*k+1.)/(2**j)  )
#    
#    
    
    #"""
    #level = 3
    P = FindP(3,4)
    Q = FindQ(3,4,'min')
    n=11
    Qshape = np.shape(Q)
    Pshape = np.shape(P)
    ncol = Qshape[1] + Pshape[1] 
    PQM1 = np.zeros((Pshape[0],ncol))
    PQM1[:,:Pshape[1]] = P
    PQM1[:,Pshape[1]:] = Q
    AB1 = np.linalg.inv(PQM1)
    #"""
    """
    *********************************
    Wavelet's for Computer Graphics (WFCG)
    
    PQM1    < = P and Q adjoined
    
    AB1     < =  inverse(PQM1) 
                    
                This is analysis (WFCG p.83), that is, to
                get a coarse level from 
                a level one 'finer' via matmul:
                    
    i.e. coarse level control vertices, c1.vertices
    
    can be recovered from the higher level 
    via projection using the PQ composition matrix
                    
        c1.vertices ==  np.matmul(c2.vertices[:].T,AB1.T).T
        
    
    >>> true #the equivalence is true, to numerical zero accuracy
    #(expressed in pseudocode, above)
    
    (note that this dense inverse n**2 time
    can be reduced, already, to linear time
    by taking advantage of a reordering of the PQM1
    matrix and using LU factorization
    to solve the resulting banded system
    in linear time)
    
    see Wavelet's for Computer Graphics (WFCG), p.[95-96]
    and they cite [Numerical Recipes] for details on the
    sparse inversion of the banded system.
    
    
    *********************************
    """
    
    u=.5
    P1 = FindP(2,2)
    
    c1 = linear_vertices((0.,12.),(12.,0.),4)
    c1 = hb.LevelSpline(c1,k=4,nump=30)
    c1 = c1.dyadic_refinement()
    c1 = c1.dyadic_refinement()
    c1 = c1.dyadic_refinement()
    c1.level = 3
    #"""
    P = FindP(3,c1.level+1)
    Q = FindQ(3,c1.level+1,'min')
    n=c1.n
    #"""
    
    Qshape = np.shape(Q)
    Pshape = np.shape(P)
    ncol = Qshape[1] + Pshape[1] 
    PQM1 = np.zeros((Pshape[0],ncol))
    PQM1[:,:Pshape[1]] = P
    PQM1[:,Pshape[1]:] = Q
    #PQM1 = np.zeros((c1.n+1,c1.n+1))
    #PQM1 = np.zeros((4+1,4+1))
    #PQM1 = np.zeros((4+1,4+1))
    #PQM1[:,0:4] = P
    #PQM1[:,4] = Q[:,0]
    #PQM1[:,4] = Q.T
    
    
    AB1 = np.linalg.inv(PQM1)
    A = AB1[:c1.n,]
    B = AB1[c1.n:]
    #page 84 :
    PHI2 = np.matmul(PQM1,AB1) #(? forgot what this is):
    
    #page 82, the big mystery:
    #see answer below!
    # key is two go up a level first!
    
    
    """
    Possible explanation: [CHECK!]
        you need to start at level 2
        and work down to level 1
    """
    print '\n\n\n---------------------------------------'
    print 'Start off:  dyadic refinement of level 1 to get level 2: '
    c2 = c1.dyadic_refinement()
    c2.level = 4
    c2b = c2.TLMBasisFuns(.5)
    c1b = c1.TLMBasisFuns(.5)
    print 'level 2 vertices are:'
    print c2.vertices
    print '\n\n\n---------------------------------------'
    print 'Result:  '
    print ' '
    print '  getting level 1 vertices '
    print '  and level 1 wavelet coefficient'
    print '  from level 2 vertices:'
    print ' '
    print 'The level 1 vertices + the wavelet scaling coefficient:'
    print np.matmul(c2.vertices[:,0].T,AB1.T)
    #print np.matmul(c2.vertices[:,0],PQM1) #unknown TLM made up
    print 'The level 1 vertices + the wavelet scaling coefficient:'
    detail_coeff1 = np.matmul(c2.vertices[:,0].T,AB1.T)[-1]
    print 'matches with the level 1 x vertices:'
    print c1.vertices
    print 'and the level 1 detail coefficient is = ',detail_coeff1
    print ''
    """
    np.dot  :: good!
    np.matmul :: good!
    np.inner :: bad!
    """
    c1_cktwoscale = np.matmul(c2b.T,PQM1)
    """
    yep dyadic refinement shows
    working two scale relations!
    
    diff is very small,
    as hoped from page 82:
    """
    #n=1
    
    print '\n\n\n'
    print '---------------------------------------'
    print 'Result: '
    print ' '
    print ' getting scale 1 basis funs(u) '
    print ' and scale 1 wavelet funs(u)'
    print ' from scale 2 basis funs(u):'
    print ' '
    c1_cktwoscale = np.matmul(c2b.T,PQM1) # recomputed from above to be
    c1_cktwoscale = np.matmul(PQM1.T,c2b)
    #also this guy:
    print 'or this way!',   np.matmul(P.T,c2b)
    # verbose as heck.  (And redundant!)
    #
    diff = c1_cktwoscale[0:c1.n] - c1b
    print 'scale 1 basis(u=.5):'
    print c1b
    print 'matmul ( c2_basis , P(3,1)|Q(3,1) )'
    print 'i.e.'
    print 'np.matmul(c2b.T,PQM1)  to get scale 2 basis funs'
    print c1_cktwoscale.T
    print '(final fun(u) in the vector above is wavelet function at u)'
    wave1 = c1_cktwoscale[4]
    print 'So the wavelet detail function value is = ',wave1
    print '\n'
    print 'checking the difference:'
    print diff
    print '\n\n\n'
    
    #
    #---------------------------------
    #
    #page 83 synthesis:
    P2 = FindP(3,c2.level+1)
    Q2 = FindQ(3,c2.level+1)
    #TODO:... figure it out!
    
    n=11
    Qshape = np.shape(Q2)
    Pshape = np.shape(P2)
    #PQM2 = np.zeros((?,?))
    #PQM2[:,0:c2.n] = P2
    #PQM2[:,5:] = Q2
    
    ncol = Qshape[1] + Pshape[1] 
    PQM2 = np.zeros((Pshape[0],ncol))
    PQM2[:,:Pshape[1]] = P2
    PQM2[:,Pshape[1]:] = Q2
    #    for e1, e2 in zip(PQM2[:,5:], Q2):
    #        print e1
    #        print e2
    #    print ''
    AB2 = np.linalg.inv(PQM2)
    
    """
    The B-spline filter bank
    is described on page 95
    """
    
    print '\nAn oddity:'
    #not wrong!:
    #x2 = np.matmul(P,c1.vertices[:,0].T)
    #y2 = np.matmul(P,c1.vertices[:,1].T)
    print 'Subdivision in Point space using P:'
    print 'where P = P(c1.n,level+1)'
    print c2.vertices - np.matmul(P,c1.vertices[:])
    print 'success, '
    print 'c2.vertices = np.matmul(P,c1.vertices[:])'
    c1.vertices - np.matmul(c2.vertices[:].T,AB1.T).T[:c1.n]
    ###c1.vertices - np.matmul(c2.vertices[:],AB1.T)[:c1.n]
    print 'go the other way:'
    print c1.vertices - np.matmul(A,c2.vertices)
    
    #this is wrong??:
    print c2.TLMBasisFuns(.5) - np.matmul(P,c1.TLMBasisFuns(.5))
    
    print c2.TLMBasisFuns(.5) - np.matmul(A.T,c1.TLMBasisFuns(.5) )
    print c2.TLMBasisFuns(.5) - np.matmul(c1.TLMBasisFuns(.5).T, A)
    
    print c1.TLMBasisFuns(.5) - np.matmul(P.T,c2.TLMBasisFuns(.5))
    
    print c2.TLMBasisFuns(.5) - np.matmul(A.T,c1.TLMBasisFuns(.5) )
    print c2.TLMBasisFuns(.5) - np.matmul(c1.TLMBasisFuns(.5).T,A )
    
    
    
    phi_wave1 = np.zeros((c2.n,1))
    phi_wave1[0:c1.n] = c1.TLMBasisFuns(.5)
    print c2.TLMBasisFuns(.5) - np.matmul(PQM1.T,phi_wave1)
    print c2.TLMBasisFuns(.5) - np.matmul(phi_wave1.T,PQM1.T)
    print c2.TLMBasisFuns(.5) - np.matmul(P,phi_wave1[0:c1.n])
    
    print c2.TLMBasisFuns(.5) - np.matmul(phi_wave1.T,PQM1.T)
    print c2.TLMBasisFuns(.5) - np.matmul(phi_wave1.T,AB1.T)
    
    
    
    
    print '\n\n\n'
    print '---------------------------------------'
    print 'Result: '
    print ' '
    print '  Get LEVEL 2 Basis Functions'
    print '  from level 1 basis functions'
    #phi_wave1 = np.zeros((c1.n+np.shape(Q)[1],1))
    #better looking:
    phi_wave1 = np.zeros((c2.n,1))
    phi_wave1[0:c1.n] = c1b
    
    #c1_cktwoscale = np.matmul(c2b.T,PQM1) #good
    #better looking:
    c1_cktwoscale = np.matmul(PQM1.T,c2b) #good
    #phi_wave1[c1.n:] = c1_cktwoscale.T[c1.n:]
    phi_wave1[c1.n:] = c1_cktwoscale[c1.n:]
    
    
    print c1b - np.matmul(c2b.T,PQM1).T[:c1.n]
    print c1b - np.matmul(PQM1.T,c2b)[:c1.n]
    print c1b - np.matmul(P.T,c2b)
    cc2b = c2b.copy()
    cc2b[8]=0.
    print c1b - np.matmul(P.T,cc2b) #not quite
    print c2b - np.matmul(P,c1b)
    print c2b - np.matmul(A.T,c1b)
    print c2b - np.matmul(c1b.T,A)
    """
    Analysis
    """
    cj_minus = np.matmul(A,c2b)
    dj_minus = np.matmul(B,c2b)
    """
    Synthesis
    """
    print c2b - (np.matmul(P,cj_minus) + np.matmul(Q,dj_minus) )
    
    """
    screwing around:
    """
    try_phi = np.zeros((c2.n,1))
    try_phi[0:c1.n] = c1b
    print c2b - (np.matmul(np.matmul(A.T,c1b).T,P) + np.matmul(Q,dj_minus) )
    cj_minus = np.matmul(A,c2b)
    dj_minus = np.matmul(B,c2b)
    
    wavelets1 = np.matmul(Q.T,c2b)
    ss = np.shape(c1b)[0]+np.shape(wavelets1)[0]
    store = np.zeros((c2.n,1),float)
    store[:c1.n] = c1b
    store[c1.n:] = wavelets1
    print c2b - np.matmul(AB1.T,store)
    store[4] = 0.
    print np.matmul(AB1.T,store)
    print c2b - np.matmul(AB1.T,store)
    print np.matmul(P.T,store)
    print np.matmul(P,store[:c1.n]) + np.matmul(B.T,dj_minus)#np.matmul(B,store[c1.n:])
    print np.matmul(P,c1b) + np.matmul(B.T,dj_minus)
    
    wavelets1 = np.matmul(Q.T,c2b)
    ss = np.shape(c1b)[0]+np.shape(wavelets1)[0]
    fine = np.zeros((c2.n,1),float)
    fine[4] = -c1b[4]
    print 'the change at Level L-1 via change at L is:',np.matmul(P.T,fine)
    
    coarse = np.zeros((c1.n,1),float)
    coarse[4] = -c1b[4]
    print 'the change at Level L via change at L-1 is:',np.matmul(P,coarse)
    
    
    
    
    
    #print c2b - np.matmul(c1b.T,PQM1)
    #print c2b - np.matmul(B,c2b)
    """Coarse to fine:  DO IT!:"""
    print c1b - np.matmul(c2b.T,PQM1).T[:c1.n]
    print c1b - np.matmul(P.T,phi_wave1)
    print c2b - np.matmul(phi_wave1.T,P)
    print c2b - np.matmul(AB1.T,phi_wave1)
    print c2b - np.matmul(phi_wave1.T,AB1)
    print c2b - np.matmul(phi_wave1.T,P)
    #print c2b - np.matmul(AB1,phi_wave1) #bad
    #print c2b - np.matmul(AB1.T,phi_wave1) #good
    print  c2b - np.matmul(AB1.T,c1_cktwoscale) #good
    print  c2b - np.matmul(PQM1,c1_cktwoscale) #bad
    print  c2b - np.matmul(PQM1.T,c1_cktwoscale) #bad
    #c2_cktwoscale = np.matmul(phi_wave1.T,AB1) #good
    #cb2_from_c1b = np.matmul(phi_wave1.T,AB1)
    #better looking:
    cb2_from_c1b = np.matmul(AB1.T,phi_wave1)
    print 'level 2 basis funs from level 1 basis funs:'
    print cb2_from_c1b.T
    print 'difference between computed and actual:'
    print cb2_from_c1b.T - c2b
    
    ##huh???
    
    phi_wave1 = np.zeros((c1.n+np.shape(Q)[1],1))
    phi_wave1[0:c1.n] = c1b
    
    c1_cktwoscale = np.matmul(c2b.T,PQM1)
    phi_wave1[c1.n:] = c1_cktwoscale.T[c1.n:]
    
    cb2_from_c1b = np.matmul(phi_wave1.T,AB1)
    
    
    print '\n\n\n'
    print '---------------------------------------'
    print 'Result: '
    print ' '
    # RIGHT:
    #********************
    # 1x4  dot 4x5 =>   1x5
    x2 = np.matmul(c1.vertices[:,0],P.T)
    aka_x2 = np.matmul(P,c1.vertices[:,0]) #match eq 7.6 page 83
    y2 = np.matmul(c1.vertices[:,1],P.T)
    simple_x2_from_x1 = np.matmul(P,c1.vertices)
    print '  Getting level 2 vertices '
    print '  from level 1 vertices'
    print ''
    print 'check x at scale 2 from x at scale 1:'
    print 'x at scale 2 computed via '
    print 'scale 1 vertices and P(3,1):'
    print x2
    print 'x at scale 2, via TLM dyadic refinement:'
    print c2.vertices[:,0]
    print ''
    print 'y at scale 2 computed via'
    print 'scale 1 y vertices and P(3,1):'
    print y2
    print 'y at scale 2, via TLM dyadic refinement:'
    print c2.vertices[:,1]
    print '\n\n\n'
    
    d0 = c1_cktwoscale[0,4]
    QD0 = d0*Q.T
    
    #
    #*************************************************
    #
    u=.5
    c1tlm = c1.TLMBasisFuns(.5)
    c2tlm = c2.TLMBasisFuns(.5)
    
    Pr = FindP(c1.p,c1.level)
    
    c1_cktwoscale = np.matmul(c2tlm.T,PQM1).T
    
    wrongforthis = np.matmul(PQM1,c2tlm)  #bad for this...
    
    
    #ck = np.matmul(c1tlm.T,AB1)
    
    
    
    cP = FindP(3,3)
    cQ = FindQ(3,3,'min')
    ###cn=11
    cQshape = np.shape(cQ)
    cPshape = np.shape(cP)
    cncol = cQshape[1] + cPshape[1] 
    cPQM1 = np.zeros((cPshape[0],cncol))
    cPQM1[:,:cPshape[1]] = cP
    cPQM1[:,cPshape[1]:] = cQ
    cAB1 = np.linalg.inv(cPQM1)
    
    
    wrongforthis = np.matmul(c1tlm.T,cAB1)
    ck = np.matmul(c1tlm.T,cPQM1).T
    ck = np.matmul(P,c1tlm)
    #
    #*************************************************
    #
    #print x0 + QD0
    #print y0 + QD0
    
    #quick dyadic refinement
    x3 = np.matmul(P2,c2.vertices[:,0].T)
    y3 = np.matmul(P2,c2.vertices[:,1].T)
    c3_vertices = np.matmul(c2.vertices.T,P2.T).T
    
    #d1 = c1_cktwoscale[0,4]
    #QD1 =d1*Q2.T
    
    
    
    print '\n\n\n'
    print '---------------------------------------'
    print 'Result: '
    print ' '
    print '  Get LEVEL 2 Basis Functions'
    print '  from level 1 basis functions'
    #c2basis = np.matmul(c1b.T,P.T)
    phi_wave1 = np.zeros((c1.n+np.shape(Q)[1],1))
    phi_wave1[0:c1.n] = c1b
    
    c1_cktwoscale = np.matmul(c2b.T,PQM1)
    phi_wave1[c1.n:] = c1_cktwoscale.T[c1.n:]
    
    cb2_from_c1b = np.matmul(phi_wave1.T,AB1)
    print 'level 2 basis funs from level 1 basis funs:'
    print cb2_from_c1b.T
    print 'difference between computed and actual:'
    print cb2_from_c1b.T - c2b
    
    
    ##
    ##*************************************************************
    ##
    print '  '
    print '  TODO:  MODIFY Level 2, '
    #**************************************************************
    print '  compute wavelet differences for level 1'
    print '  And reconstruct level 2 from level 1'
    #
    c2_copy = c1.dyadic_refinement()
    #
    print 'make an alternation to c2_copy!'
    vertsc = copy.copy(c2_copy.vertices)
    vertsc[3] = [8.8,5.5]
    c2_copy.vertices = vertsc
    
    print '\n\n\n---------------------------------------'
    print 'modify the level 2 scaling coefficients '
    print 'and propogate that to the level 1 coefficients'
    print 'plus wavelets????'
    print 'Result:  WRONG??'
    
    print ' '
    print ' getting scale 1 basis funs(u) '
    print ' and scale 1 wavelet funs(u)'
    print ' from scale 2 basis funs(u):'
    print ' '
    
    c2cb = c2_copy.TLMBasisFuns(.5)
    c1_alt_basis = np.matmul(c2cb.T,PQM1) #same
    c1_alt_basis = np.matmul(PQM1.T,c2cb) #better looking
    print ' '
    print '  getting level 1 vertices '
    print '  and level 1 detail coefficient'
    print '  from level 2 vertices:'
    print ' '
    print 'The level 1 vertices + the wavelet scaling coefficient:'
    vw = np.matmul(c2_copy.vertices[:,0].T , AB1.T)
    print vw
    detail_coeff1 = vw[-1]
    print 'matches with the level 1 x vertices:'
    print c1.vertices
    print 'and the level 1 detail coefficient is = ',detail_coeff1
    print ''
    
    
    print '\n\n\n---------------------------------------'
    print 'modify the detail coefficients of level 1'
    print 'Result:  '
