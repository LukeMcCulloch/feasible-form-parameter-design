# Luke McCulloch
# Updating the polyintegration notes

# Side Note:
#   Turns out that the nodes of gaussian quadrature are the
#   eigenvelues
#   of a
#   tridiagonal matrix

import numpy as np
# Derivative Functions:
#from dBasis import dinfluence #edit, nope this is passed in
    # because in general we can integrate any old function...
#from bSplineBasis import *
#from curve          import *
#from bSplineBasis   import influence, oldinfluence
#from dBasis         import dinfluence
from gaussxw        import gaussxw, gaussxwab

from KV import make_knot_vector
#from curve import *
#from quickspline import make_knot_vector


## easiest to check the integration of basis functions from the curve module!

## verified : dbasis_quad_Check1 (1 basis function)
## verified : dbasis_quad_Check  (2 basis functions)

def lamda1(x):
    """
    Test function
    """
    return x*x

def lamda2(x):
    """
    Test function
    """
    return x

def lamda3(x):
    """
    Test function
    """
    return np.sqrt(x)

def Deltlamda1(x):
    """
        Test function:
        function to test discontinuous integration
        using the gaussian function
    """
    if x<.5:
        return 1.
    else:
        return 0.

def linlamda1(x):
    if x<=(1/3.):
        return -6. + 18.*x
    else:
        return 0.


def FindSpan(n,p,U,u):
    #n=self.n
    #p=self.p
    #U=self.t
    """
        Determines the knot Span Index
        Needs testing/editing for multiplicative knots.
        B/c it gets stuck when there are multiples
            of the same knot value on the interior
            of the knot vector.
        Input:
                n = number of vertices (contrast with P&T, who use # vertices -1)
                p = degree of curve
                u = knot position
                U = knot vector

        Output:
                mid = the knot span index"""

    
    if u==U[n+1]:
        return n-1 #special case n = m-p-1 #TLM mod 6-10-2013

    #Do a Binary Search
    low  = p
    high = n+1

    mid = (low + high)/2
    while (u < U[mid] or u >= U[mid+1]):
        if u < U[mid]:
            high = mid
        else:
            low = mid
        mid = (low+high)/2
    return mid

def basis(i,k,U,u,N):
    """
    P&T module to compute the
    Cox De-Boor Basis Functions
    Efficiently at a point

    Inputs:
        i = span index
        u = point on the parameterization
        k = order of the basis function (passed for ease of code update)
        U = knot vector
        N = EmptyArray of Basis functions from [0...p], i.e. size k, at u.
        
    Output:
        N = Array of all Basis Function Values at u.  """
    left=np.zeros((k),float)
    right=np.zeros((k),float)
    p=k-1
    N[0]=1.0
    for j in range(1,k): 
        left[j]  =u-U[i+1-j]
        right[j] =U[i+j]-u
        saved=0.0
        for r in range(0,j):
            temp  = N[r]/(right[r+1]+left[j-r])
            N[r]  = saved+right[r+1]*temp
            saved = left[j-r]*temp
        
        N[j]=saved
    return N


def DersBasisFunc(i,u,k,U,ders):
    """
    P&T module to compute all the
        Non Vanishing Derivatives of 
        Cox De-Boor Basis Functions
        Efficiently at a point

    Inputs:
        i = span index
        u = point on the parameterization
        k = order of the basis function (passed for ease of code update)
        U = knot vector
        ders = EmptyArray of Basis functions from [0...p], i.e. size k, at u.
        
    Output:
        ders = Array of all Basis Function Values at u.   """
    
    #k=self.k
    p=k-1
    n=k #for n<=p on page 72.  so not n=self.n.
    U#=t

    left=np.zeros((k),float)
    right=np.zeros((k),float)
    
    ndu = np.zeros((k,k),float) #store basis funcs and knot differences
    a   = np.zeros((2,k),float) #store two most recently computed fows of ak,j ak-1,j
    
    ndu[0,0]=1.0
    for j in range(1,p+1):
        left[j]  = u-U[i+1-j]
        right[j] = U[i+j]-u
        saved = 0.0
        for r in range(j):
                                                #upper tri
            ndu[j,r]    = right[r+1]+left[j-r]
            temp        = ndu[r,j-1]/ndu[j,r]
                                                #lower tri
            ndu[r,j]    = saved+right[r+1]*temp
            saved = left[j-r]*temp
        ndu[j,j] = saved
        
    for j in range(p+1): #LOAD THE BASIS FUNCTIONS
        ders[0,j]=ndu[j,p]
        
    for r in range(p+1):
        s1=0
        s2=1
        a[0,0]=1.0
        for k in range(1,n):
            d=0.0
            rk = r-k
            pk = p-k
            if r>=k:
                a[s2,0] = a[s1,0]/ndu[pk+1,rk]
                d       = a[s2,0]*ndu[rk,pk]
            if rk>=-1:
                j1=1
            else:
                j1=-rk
            if ((r-1)<=pk):
                j2=k-1
            else:
                j2=p-r
            for j in range(j1,j2+1):
                a[s2,j]=(a[s1,j]-a[s1,j-1])/ndu[pk+1,rk+j]
                d += a[s2,j]*ndu[rk+j,pk]
            if r<=pk:
                a[s2,k]=-a[s1,k-1]/ndu[pk+1,r]
                d+=a[s2,k]*ndu[r,pk]
            ders[k,r] = d
            j=s1
            s1=s2
            s2=j
    r=p
    for k in range(1,n):
        for j in range(p+1):
            ders[k,j] *= r
        r*=(p-k)
    return #ders


##
## ---- ---- ---- Start of Any Order Quadrature Methods ---- ---- ---- ##
##

def ibasis_quad(curve, N, func1, t):
    """
        Modified to use correct starting and end vertices for these product basis
            functions.
        N   = order of accuracy of the integration - should match the order of the new D
                basis Dbasis function
        func1 = oldbasisfunc
        n   = number of vertices
        k   = order of the curve
        t   = knot vector
    """

    
    # !Watch your implementation, it appears you passed k, not p, to p in dbasis_quad:  TLM 12-11-2012

    
    #itegration interval
    # based where the function actually exists -
    # use the knot vector
    # based on ...
    # this domain should be based on the min set of active knots...
    # there is a range dependent on each functon which is passed...
    

    ##  ESTABLISH THE LIMITS OF INTEGRATION:
    ##      a=min, b=max
    #a=curve.t[i]
    #b=curve.t[k+i]
    #a=curve.t[0]
    #b=curve.t[-1]

    
    # Calculate the sample vertices and weights, then map them
    # to the required integration domain
    # or, if the domain is from 0 to 1 or remapped later:
    # use: x,w = gaussxw(N)
    
    #x,w = gaussxwab(N,a,b)
    #print 'x={},w={}'.format(x,w)
    """
        returns integration vertices, x
        and
        integration weights, w
        mapped to the integral of a to b
    """

    # these are only used if you use the naive algorithm "gaussxw(N)"
    #xp = 0.5*(b-a)*x + 0.5*(b+a)
    #wp = 0.5*(b-a)*w
    #print 'x = {}; w = {}'.format(x,w)
    
    # Perform the integration


    


    for i in range(curve.n):
        #if (len(t)<2*curve.n):
        a=t[i]
        #b=t[curve.k+i]
        b=t[N+i]
##        else:
##            a=t[i+1]
##            b=t[curve.k+i+1]
        #print 'a = {}'.format(a)
        #print 'b = {}'.format(b)
        x,w = gaussxwab(N,a,b)  #N is the order of the basis function
        #print 'N = {}'.format(N)
        #print 'x = {}'.format(x)
        #print 'w = {}'.format(w)

        ## Gaussian Quadrature
        s = 0.0
        for pt in range(N):
            store = func1(i,N,t,x[pt]) 
            s += w[pt]*store        #sum over all N

        ## Check integrals with retangular approximation
        #dx=curve.s[1]-curve.s[0]
        #for pt in range(1,len(curve.s),1):
        #    store = dx*func1(i,N,t,curve.s[pt])
        #    s += store
            #print 'Basis function value = {} with w({})={}'.format(store,pt,w[pt])

        curve.ibasis[i]=s
        
        #print 'integral is ={}'.format(s)
    return





def ibasis_quad2(curve, N, i, j, k,t,d):
    """
        Modified to use correct starting and end vertices for these product basis
            functions.
        N   = order of accuracy of the integration - should match the order of the new D
                basis Dbasis function
        i   = the ith function in the Mij matrix
        j   = the jth function im the Mij matrix
        p   = degree (actually passing order) of the curve
        t   = knot vector
        d   = derivative order
    """
    p=curve.p
    ##  ESTABLISH THE LIMITS OF INTEGRATION:
    ij_min = max(i,j) #could use this always assuming assending knot order...
    ij_max = min(k+i,k+j)

    ##
    ##
    ##
    #be very careful - use integrate between discontinuities
    #Now we need to divide the span into continuous segments, if not already continuous...
    rel_knot_set = list(set(t[ ij_min : ij_max+1 ]))  # oh boy... am I intoducing something weird here?
    rel_knot_set=sorted(rel_knot_set)
    #print 'i={},j={},ij_min={},ij_max={},rel_knot_set={}'.format(i,j,ij_min,ij_max,rel_knot_set)
    
    # Calculate the sample vertices and weights, then map them
    # to the required integration domain
    s       = 0.0
    #s_check = 0.0
    
    if(len(rel_knot_set)<2):
            #print 'pass'
            pass
    else:
        for ii in range(len(rel_knot_set)-1):

            a=rel_knot_set[ii]
            b=rel_knot_set[ii+1]
            #a=t[ij_min]
            #b=t[ij_max]

            
            # Calculate the sample vertices and weights, then map them
            x,w = gaussxwab(N,a,b)
            #print 'x={},w={}'.format(x,w)
            """
                returns integration vertices, x
                and
                integration weights, w
                mapped to the integral of a to b
            """

            
            # Perform the integration
            ch1 = np.zeros((curve.n,curve.n),float)
            for pt in range(N):
                ch1[:,:]=0.0
                span = curve.FindSpan(x[pt])
                curve.DersBasisFunc(span,x[pt],ch1[:,span-p:span+1])
                s += w[pt]*ch1[0,i]*ch1[1,j]
    return s

def ibasis_quad2_test(curve, N, i, j, k,t,d):
    """
        Modified to use correct starting and end vertices for these product basis
            functions.
        N   = order of accuracy of the integration - should match the order of the new D
                basis Dbasis function
        i   = the ith function in the Mij matrix
        j   = the jth function im the Mij matrix
        p   = degree (actually passing order) of the curve
        t   = knot vector
        d   = derivative order
    """
    p=curve.p
    ##  ESTABLISH THE LIMITS OF INTEGRATION:
    ij_min = max(i,j) #could use this always assuming assending knot order...
    ij_max = min(k+i,k+j)


    rel_knot_set = list(set(t[ ij_min : ij_max+1 ]))  #visual check with curve.plotbasis(num)
    rel_knot_set=sorted(rel_knot_set)
    print 'i={},j={},ij_min={},ij_max={},rel_knot_set={}'.format(i,j,ij_min,ij_max,rel_knot_set)
    
    # Calculate the sample vertices and weights, then map them
    # to the required integration domain
    s       = 0.0
    s_check = 0.0
    if(len(rel_knot_set)<2):
            #print 'pass'
            pass
    else:
        for ii in range(len(rel_knot_set)-1):

            a=rel_knot_set[ii]
            b=rel_knot_set[ii+1]
            #a=t[ij_min]
            #b=t[ij_max]

            
            # Calculate the sample vertices and weights, then map them
            x,w = gaussxwab(N,a,b)
            print 'x={},w={}'.format(x,w)
            """
                returns integration vertices, x
                and
                integration weights, w
                mapped to the integral of a to b
            """

            
            # Perform the integration
            ch1 = np.zeros((curve.n,curve.n),float)
            
            for pt in range(N):
                ch1[:,:]=0.0
                span = curve.FindSpan(x[pt])
                curve.DersBasisFunc(span,x[pt],ch1[:,span-p:span+1])
                s += w[pt]*ch1[0,i]*ch1[1,j]
                #""" simple check program:
                
        rs=0.
        b2 = np.zeros((curve.n,curve.n),float)
        lena=3000
        param = np.linspace(0,1,lena,endpoint=True)
        for ipt in range(lena):
            b2[:,:]=0.
            span = curve.FindSpan(param[ipt])
            curve.DersBasisFunc(span,param[ipt],b2[:,span-p:span+1])
            rs += b2[0,i]*b2[1,j]
        s_check = rs/lena
        #"""
    return s, s_check

def Mbasis_quad(curve, N, i, j, l, k,t,d):
    """
        Moments:
        Multiply 3 Basis Function Together and Integrate
        
        Modified to use correct starting and end vertices for these product basis
            functions.
        N   = order of accuracy of the integration - should match the order of the new D
                basis Dbasis function
        i   = the ith function in the Mij matrix
        j   = the jth function im the Mij matrix
        l   = the lth function im the Mij matrix
        k   = order of the curve
        t   = knot vector
        d   = derivative order
    """
    p=curve.p

    ##  ESTABLISH THE LIMITS OF INTEGRATION:
    ##      a=min, b=max
    #-----------------------
    ijl_min = max(i,j,l)
    ijl_max = min(k+i,k+j,k+l)

    ##
    ##
    ##
    #be very careful - use integrate between discontinuities
    #Now we need to divide the span into continuous segments, if not already continuous...
    rel_knot_set = list(set(t[ ijl_min : ijl_max+1 ]))  # oh boy... am I intoducing something weird here?
    rel_knot_set=sorted(rel_knot_set)
    #print 'i={},j={},ij_min={},ij_max={},rel_knot_set={}'.format(i,j,ij_min,ij_max,rel_knot_set)
    
    # Calculate the sample vertices and weights, then map them
    # to the required integration domain
    s       = 0.0
    s_check = 0.0
    if(len(rel_knot_set)<2):
            #print 'pass'
            pass
    else:
        for ii in range(len(rel_knot_set)-1):
    ##
    ##
    ##
            a=rel_knot_set[ii]
            b=rel_knot_set[ii+1]
            #a=t[ijl_min]
            #b=t[ijl_max]
            #-----------------------
            
            # Calculate the sample vertices and weights, then map them
            x,w = gaussxwab(N,a,b)

            
            # Perform the integration
            ch1 = np.zeros((curve.n,curve.n),float)
            for pt in range(N):
                ch1[:,:]=0.0
                span = curve.FindSpan(x[pt])
                curve.DersBasisFunc(span,x[pt],ch1[:,span-p:span+1])
                s += w[pt]*ch1[0,l]*ch1[0,i]*ch1[1,j]
    return s
    

def dbasis_quad(curve,N, i, j, k,t,d ):
    """
        For E norms:
        
        Multiply 2 Derivative Basis Functions and Integrate
        
        Modified to use correct starting and end vertices for these product basis
            functions.
        N   = order of accuracy of the integration - should match the order of the new D
                basis Dbasis function
        i   = the ith function in the Mij matrix
        j   = the jth function im the Mij matrix
        p   = degree (actually passing order) of the curve
        t   = knot vector
        d   = derivative order
    """
    p=curve.p

    ##  ESTABLISH THE LIMITS OF INTEGRATION:
    ##      a=min, b=max
    ij_min = max(i,j) #could use this always assuming assending knot order...
    ij_max = min(k+i,k+j)

    
    #Now we need to divide the span into continuous segments, if not already continuous...
    rel_knot_set = list(set(t[ ij_min : ij_max+1 ]))  # oh boy... am I intoducing something weird here?
    rel_knot_set=sorted(rel_knot_set) #cruical!

  
    # Calculate the sample vertices and weights, then map them
    # to the required integration domain
    s = 0.0
    if(len(rel_knot_set)<2):
            pass
    else:
        for ii in range(len(rel_knot_set)-1):
            a=rel_knot_set[ii]
            b=rel_knot_set[ii+1]
            x,w = gaussxwab(N,a,b)
            """
                returns integration points, x
                and
                integration weights, w
                mapped to the integral of a to b
            """
            
            # Perform the integration
            b1 = np.zeros((curve.n,curve.n),float)
            for pt in range(N):
                b1[:,:]=0.0
                span = curve.FindSpan(x[pt])
                curve.DersBasisFunc(span,x[pt],b1[:,span-p:span+1])
                s += w[pt]*b1[d,i]*b1[d,j] #sum over all N

    return s


def dbasis_quad_Check(curve,N, i, j, k,t,d ):
    """
        For E norms:
        
        Multiply 2 Derivative Basis Functions and Integrate
        
        Modified to use correct starting and end vertices for these product basis
            functions.
        N   = order of accuracy of the integration - should match the order of the new D
                basis Dbasis function (product)
        i   = the ith function in the Mij matrix
        j   = the jth function im the Mij matrix
        p   = degree (actually passing order) of the curve
        t   = knot vector
        d   = derivative order
    """
    p=curve.p

    ##  ESTABLISH THE LIMITS OF INTEGRATION:
    ij_min = max(i,j) #could use this always assuming assending knot order...
    ij_max = min(k+i,k+j)
    
    #Now we need to divide the span into continuous segments, if not already continuous...
    rel_knot_set = list(set(t[ ij_min : ij_max+1 ]))  # oh boy... am I intoducing something weird here?
    rel_knot_set=sorted(rel_knot_set)
    print 'i={},j={},ij_min={},ij_max={},rel_knot_set={}'.format(i,j,ij_min,ij_max,rel_knot_set)
    
    # Calculate the sample vertices and weights, then map them
    # to the required integration domain
    s       = 0.0
    s_check = 0.0
    if(len(rel_knot_set)<2):
            #print 'pass'
            pass
    else:
        for ii in range(len(rel_knot_set)-1):
            a=rel_knot_set[ii]
            b=rel_knot_set[ii+1]
            x,w = gaussxwab(N,a,b)
            """
                returns integration points, x
                and
                integration weights, w
                mapped to the integral of a to b
            """
 
            # Perform the integration
            b1 = np.zeros((curve.n,curve.n),float)
            for pt in range(N):
                b1[:,:]=0.
                span = curve.FindSpan(x[pt])
                curve.DersBasisFunc(span,x[pt],b1[:,span-p:span+1])
                s += w[pt]*b1[d,i]*b1[d,j] #sum over all N
 

                
        #""" simple check program:
        rs=0.
        b2 = np.zeros((curve.n,curve.n),float)
        lena=15000
        param = np.linspace(0,1,lena,endpoint=True)
        for ipt in range(lena):
            b2[:,:]=0.
            span = curve.FindSpan(param[ipt])
            curve.DersBasisFunc(span,param[ipt],b2[:,span-p:span+1])
            rs += b2[d,i]*b2[d,j]
        s_check = rs/lena
        #"""
    return s, s_check

def dbasis_quad_Check1(curve,N, i, k,t,d ):
    """
        For single basis function integration:
        
        Multiply 2 Derivative Basis Functions and Integrate
        
        Modified to use correct starting and end vertices for these product basis
            functions.
        N   = order of accuracy of the integration - should match the order of the new
                Dbasis Dbasis function (product)
        i   = the ith function in the Mij matrix
        j   = the jth function im the Mij matrix
        p   = degree (actually passing order) of the curve
        t   = knot vector
        d   = derivative order
    """
    p=curve.p

    ##  ESTABLISH THE LIMITS OF INTEGRATION:
    ij_min = i #could use this always assuming assending knot order...
    ij_max = i+k
    
    #Now we need to divide the span into continuous segments, if not already continuous...
    rel_knot_set = list(set(t[ ij_min : ij_max+1 ]))  # oh boy... am I intoducing something weird here?
    rel_knot_set=sorted(rel_knot_set)
    print 'i={},ij_min={},ij_max={},rel_knot_set={}'.format(i,ij_min,ij_max,rel_knot_set)
    
    # Calculate the sample vertices and weights, then map them
    # to the required integration domain
    s       = 0.0
    s_check = 0.0
    if(len(rel_knot_set)<2):
            pass
    else:
        for ii in range(len(rel_knot_set)-1):
            a=rel_knot_set[ii]
            b=rel_knot_set[ii+1]
            x,w = gaussxwab(N,a,b)
            """
                returns integration points, x
                and
                integration weights, w
                mapped to the integral of a to b
            """
 
            # Perform the integration
            b1 = np.zeros((curve.n,curve.n),float)
            for pt in range(N):
                b1[:,:]=0.0
                span = curve.FindSpan(x[pt])
                curve.DersBasisFunc(span,x[pt],b1[:,span-p:span+1])
                s += w[pt]*b1[d,i] #sum over all N
 

                
    #""" simple check program:
    rs=0.
    b2 = np.zeros((curve.n,curve.n),float)
    lena=10000
    param = np.linspace(0,1,lena,endpoint=True)
    for ipt in range(lena):
        b2[:,:]=0.
        span = curve.FindSpan(param[ipt])
        curve.DersBasisFunc(span,param[ipt],b2[:,span-p:span+1])
        #print b2[d,i]
        rs += b2[d,i]
    s_check = rs/float(lena-1)
    #"""
    return s, s_check




def gauss_int(N, func1, func2 ):
    """
        Testing simple functions
        using gaussian integration:
        
        N   = order of accuracy of the integration - should match the order of the new D
                basis Dbasis function, or be higher - this shouldn't hurt but should
                be slower

        func1 = a function
        func2 = a function
    """

    #itegration interval
    a=0.
    b=1.


    # Calculate the sample vertices and weights, then map them
    # to the required integration domain
    #x,w = gaussxw(N)
    x,w = gaussxwab(N,a,b)
    #print 'x = {}'.format(x)
    #print 'w = {}'.format(w)
    """
        input:
            N = 
            a = min
            b = max
        returns integration vertices, x
        and
        integration weights, w
        mapped to the integral of a to b
    """

    # these are only used if you use the naive algorithm "gaussxw(N)"
    #xp = 0.5*(b-a)*x + 0.5*(b+a)
    #wp = 0.5*(b-a)*w
    #print 'x = {}; w = {}'.format(x,w)
    
    # Perform the integration
    s = 0.0
    for k in range(N):
        s += w[k]*func1(x[k])*func2(x[k])
        #print 'k={}; x[{}] = {}, w[{}]={}'.format(k,k,x[k],k,w[k])
        #print 's={}'.format(s)
    #print 'for k = N={}, Final s={}'.format(k, s)
    return s
    
##def gauss_intbasis(i,k,t,N,n,D):
##    """
##        Testing simple functions
##        using gaussian integration:
##
##        inputs:
##            i = which basis function to integrate
##            k = order of the curve
##            t = knot vector        
##            N = order of accuracy of the integration - should match the order of the new D
##                    basis Dbasis function, or be higher - this shouldn't hurt but should
##                    be slower
##            n = number of vertices : needed for FindSpan routine.
##            D = order of basis function derivative
##
##    """
##    p=k-1
##
##    # store the basis evaluations:
##    locBasis=np.zeros((k,n))
##
##    #itegration interval : local support spans the set [i,i+k]
##    a=t[i]
##    b=t[k+i]
##    print 'a = {}'.format(a)
##    print 'b = {}'.format(b)
##    
##    # Calculate the sample vertices and weights, then map them
##    # to the required integration domain
##    x,w = gaussxwab(N,a,b)
##
##    print 'x = {}'.format(x)
##    print 'w = {}'.format(w)
##    """
##        input:
##            N = 
##            a = min
##            b = max
##        returns integration vertices, x
##        and
##        integration weights, w
##        mapped to the integral of a to b
##    """
##
##    # these are only used if you use the naive algorithm "gaussxw(N)"
##    #xp = 0.5*(b-a)*x + 0.5*(b+a)
##    #wp = 0.5*(b-a)*w
##    #print 'x = {}; w = {}'.format(x,w)
##    
##    # Perform the integration
##    s = 0.0
##    store = 0.
##    for pt in range(N):
##        locBasis[:,:]=0.
##        #store = basis(i,k,t,x[pt])
##        print 'integrate using value:',x[pt]
##        spanloc=FindSpan(n,p,t,x[pt])
##        DersBasisFunc(spanloc,x[pt],k,t,locBasis[:,spanloc-p:spanloc+1])
##        store=locBasis[D,:]
##        s += w[pt]*store[i]
##        print 'basis function value = {} with w({})={}'.format(store,pt,w[pt])
##    print 'integral is ={}'.format(s)
##
##    check_s=0.0
##
##    return s, check_s

def VolKernelQuatdrature(surf, N, ix,jx,iy,jy,iz,jz ):
    """

    
        surf = the object whose Volume Kernel is to be evaluated
        N   = list of orders of accuracy of the integration - should match the order of the function
        i#  = | 
        j#  = |---> these denote which basis function is being computed for which dimension
        uk  = order of the surface in u
        vk  = order of the surface in v
        t   = knot vector
    """
    uk = surf.uorder
    ut = surf.uknots
    un = surf.un
    
    vk = surf.vorder
    vt = surf.vknots
    vn = surf.vn
    
    up=uk-1
    vp=vk-1
    ##  ESTABLISH THE LIMITS OF INTEGRATION:
    ##      we need seperate limits
    ##      for the "u-moving stuff"
    ##      and the "v-moving stuff"
    #-----------------------
    #u_ijl_min = max(ix,iy,iz,jx,jy,jz)
    #u_ijl_max = min(uk+ix, uk+iy, uk+iz,vk+jx, vk+jy, vk+jz)
    u_ijl_min = max(ix,iy,iz)
    u_ijl_max = min(uk+ix, uk+iy, uk+iz)

    v_ijl_min = max(jx,jy,jz)
    v_ijl_max = min(vk+jx, vk+jy, vk+jz)

    #integrate continuous segments piecewise
    u_rel_knot_set = list(set(ut[ u_ijl_min : u_ijl_max+1 ])) 
    u_rel_knot_set=sorted(u_rel_knot_set)

    v_rel_knot_set = list(set(vt[ v_ijl_min : v_ijl_max+1 ])) 
    v_rel_knot_set=sorted(v_rel_knot_set)


    
    # Calculate the sample vertices and weights,
    # map them to the required integration domain
    s       = 0.0
    #s_check = 0.0



    ## U PART:
    if(len(u_rel_knot_set)<2 or len(v_rel_knot_set)<2):
            pass #integral is zero
    else:
        for jj in range(len(v_rel_knot_set)-1):
            av=v_rel_knot_set[jj] 
            bv=v_rel_knot_set[jj+1]

            xv,wv = gaussxwab(N[1],av,bv)
            
            for ii in range(len(u_rel_knot_set)-1):
                au=u_rel_knot_set[ii] 
                bu=u_rel_knot_set[ii+1]
    
                xu,wu= gaussxwab(N[0],au,bu)
    
                ch1 = np.zeros((un,un),float)
                ch2 = np.zeros((vn,vn),float)
                for xpt in range(N[0]):
                    for ypt in range(N[1]):
                    
                        ch1[:,:]=0.0
                        uspan = surf.uFindSpan(xu[xpt])
                        surf.uDersBasisFunc(uspan,xu[xpt],ch1[:,uspan-up:uspan+1])
                        
                        ch2[:,:]=0.0
                        vspan = surf.vFindSpan(xv[ypt])
                        surf.vDersBasisFunc(vspan,xv[ypt],ch2[:,vspan-vp:vspan+1])
                        
                        s += wv[ypt]*wu[xpt]*ch1[0,iz]*ch2[0,jz]*\
                                    (ch1[1,ix]*ch1[0,iy]*ch2[0,jx]*ch2[1,jy] - \
                                    ch1[0,ix]*ch1[1,iy]*ch2[1,jx]*ch2[0,jy])


    ## V PART:
#    if(len(v_rel_knot_set)<2):
#            pass #integral is zero
#    else:
#        for ii in range(len(v_rel_knot_set)-1):
#            a=v_rel_knot_set[ii] 
#            b=v_rel_knot_set[ii+1]
#
#            x,w = gaussxwab(N,a,b)
#
#            ch2 = np.zeros((vn,vn),float)
#            for pt in range(N):
#                ch2[:,:]=0.0
#                span = surf.vFindSpan(x[pt])
#                surf.vDersBasisFunc(span,x[pt],ch2[:,span-vp:span+1])
#                s += w[pt]*(ch2[0,jz]*ch2[0,jx]*ch2[1,jy] - ch2[0,jz]*ch2[1,jx]*ch2[0,jy])


    return s
    
    
def VolKernelQuatdrature_u(surf, N, ix,jx,iy,jy,iz,jz ):
    """

    
        surf = the object whose Volume Kernel is to be evaluated
        N   = order of accuracy of the integration - should match the order of the function
        i#  = | 
        j#  = |---> these denote which basis function is being computed for which dimension
        uk  = order of the surface in u
        uv  = order of the surface in v
        t   = knot vector
    """
    uk = surf.uorder
    ut = surf.uknots
    un = surf.un
    
    vk = surf.vorder
    vt = surf.vknots
    vn = surf.vn
    
    up=uk-1
    vp=vk-1
    ##  ESTABLISH THE LIMITS OF INTEGRATION:
    ##      we need seperate limits
    ##      for the "u-moving stuff"
    ##      and the "v-moving stuff"
    #-----------------------
    u_ijl_min = max(ix,iy,iz)
    u_ijl_max = min(uk+ix, uk+iy, uk+iz)

    v_ijl_min = max(jx,jy,jz)
    v_ijl_max = min(vk+jx, vk+jy, vk+jz)

    #integrate continuous segments piecewise
    u_rel_knot_set = list(set(ut[ u_ijl_min : u_ijl_max+1 ])) 
    u_rel_knot_set=sorted(u_rel_knot_set)

    v_rel_knot_set = list(set(vt[ v_ijl_min : v_ijl_max+1 ])) 
    v_rel_knot_set=sorted(v_rel_knot_set)


    
    # Calculate the sample vertices and weights,
    # map them to the required integration domain
    s       = 0.0
    s_check = 0.0



    ## U PART:
    if(len(u_rel_knot_set)<2):
            pass #integral is zero
    else:
        for ii in range(len(u_rel_knot_set)-1):
            a=u_rel_knot_set[ii] 
            b=u_rel_knot_set[ii+1]

            x,w = gaussxwab(N,a,b)

            ch1 = np.zeros((un,un),float)
            for pt in range(N):
                ch1[:,:]=0.0
                span = surf.uFindSpan(x[pt])
                surf.uDersBasisFunc(span,x[pt],ch1[:,span-up:span+1])
                s += w[pt]*(ch1[0,iz]*ch1[1,ix]*ch1[1,iy] - ch1[0,iz]*ch1[0,ix]*ch1[1,iy])


    ## V PART:
    if(len(v_rel_knot_set)<2):
            pass #integral is zero
    else:
        for ii in range(len(v_rel_knot_set)-1):
            a=v_rel_knot_set[ii] 
            b=v_rel_knot_set[ii+1]

            x,w = gaussxwab(N,a,b)

            ch2 = np.zeros((vn,vn),float)
            for pt in range(N):
                ch2[:,:]=0.0
                span = surf.vFindSpan(x[pt])
                surf.vDersBasisFunc(span,x[pt],ch2[:,span-vp:span+1])
                s += w[pt]*(ch2[0,jz]*ch2[0,jx]*ch2[1,jy] - ch2[0,jz]*ch2[1,jx]*ch2[0,jy])


    return s
    
def VolKernelQuatdrature_v(surf, N, ix,jx,iy,jy,iz,jz ):
    """

    
        surf = the object whose Volume Kernel is to be evaluated
        N   = order of accuracy of the integration - should match the order of the function
        i#  = | 
        j#  = |---> these denote which basis function is being computed for which dimension
        uk  = order of the surface in u
        uv  = order of the surface in v
        t   = knot vector
    """
    uk = surf.uorder
    ut = surf.uknots
    un = surf.un
    
    vk = surf.vorder
    vt = surf.vknots
    vn = surf.vn
    
    up=uk-1
    vp=vk-1
    ##  ESTABLISH THE LIMITS OF INTEGRATION:
    ##      we need seperate limits
    ##      for the "u-moving stuff"
    ##      and the "v-moving stuff"
    #-----------------------
    u_ijl_min = max(ix,iy,iz)
    u_ijl_max = min(uk+ix, uk+iy, uk+iz)

    v_ijl_min = max(jx,jy,jz)
    v_ijl_max = min(vk+jx, vk+jy, vk+jz)

    #integrate continuous segments piecewise
    u_rel_knot_set = list(set(ut[ u_ijl_min : u_ijl_max+1 ])) 
    u_rel_knot_set=sorted(u_rel_knot_set)

    v_rel_knot_set = list(set(vt[ v_ijl_min : v_ijl_max+1 ])) 
    v_rel_knot_set=sorted(v_rel_knot_set)


    
    # Calculate the sample vertices and weights,
    # map them to the required integration domain
    s       = 0.0
    s_check = 0.0



    ## U PART:
    if(len(u_rel_knot_set)<2):
            pass #integral is zero
    else:
        for ii in range(len(u_rel_knot_set)-1):
            a=u_rel_knot_set[ii] 
            b=u_rel_knot_set[ii+1]

            x,w = gaussxwab(N,a,b)

            ch1 = np.zeros((un,un),float)
            for pt in range(N):
                ch1[:,:]=0.0
                span = surf.uFindSpan(x[pt])
                surf.uDersBasisFunc(span,x[pt],ch1[:,span-up:span+1])
                s += w[pt]*(ch1[0,iz]*ch1[1,ix]*ch1[1,iy] - ch1[0,iz]*ch1[0,ix]*ch1[1,iy])


    ## V PART:
    if(len(v_rel_knot_set)<2):
            pass #integral is zero
    else:
        for ii in range(len(v_rel_knot_set)-1):
            a=v_rel_knot_set[ii] 
            b=v_rel_knot_set[ii+1]

            x,w = gaussxwab(N,a,b)

            ch2 = np.zeros((vn,vn),float)
            for pt in range(N):
                ch2[:,:]=0.0
                span = surf.vFindSpan(x[pt])
                surf.vDersBasisFunc(span,x[pt],ch2[:,span-vp:span+1])
                s += w[pt]*(ch2[0,jz]*ch2[0,jx]*ch2[1,jy] - ch2[0,jz]*ch2[1,jx]*ch2[0,jy])


    return s

def main():
    """
        Testing functions
    """


    N=1 #integrate a straight line
    value = gauss_int(N, lamda3, lamda3)
    print 'integrate x from 0 to 1 =',value
    print ''

    N=2 #integrate a quadratic, x^2
    value = gauss_int(N, lamda2, lamda2)
    print 'integrate x^2 dx from 0 to 1 = {}'.format(value)
    print ''

    N=4 #integrate x^4 4th power
    value = gauss_int(N, lamda1, lamda1)
    print 'integrate(x^4)dx from 0 to 1 = {}'.format(value)
    print ''

    N=9 #integrate a disontinuous straight line
    value = gauss_int(N, Deltlamda1, Deltlamda1)
    print value
    print ''

    N=9 #integrate a disontinuous straight line
    value = gauss_int(N, linlamda1, linlamda1)
    print value
    print ''


    return 
if __name__ == '__main__':
    main()


