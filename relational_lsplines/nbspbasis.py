#!/usr/bin/env python
##
## Routines for creation and evaluation of B-spline curves and
## surfaces
##
##
## $RCSfile$   $Revision$
## Copyright: Lothar Birk
## Last update: 17.11.2005
##
## Last update: $Date$ by $Author$
##
## $State$
##
## $Log$
##

#from numarray import *
#from matplotlib.numerix import *
from numpy import *
import sys
import string


shift = 0.000001
MINFLOAT = 1.E-12
lb3dFLOAT = float

def make_knot_vector(k, n, aarray, type='open'):
    """
        make_knot_vector

        Routine to make knot vectors for a B-spline basis

        k = order of Basis (or curve, surface), i.e. degree is k-1
            (example: fourth order curve is cubic)
        n = highest vertex index. There are n+1 vertices in total
            Numbers starting with 0.
    """

    #print aarray
    if len(aarray) == 2:
        return make_uniform_knot_vector(k, n, aarray[0], aarray[1], type)
    else:
        if type == 'open':
            Nknots    = n+k+1
            Nintknots = n-k+1
            Nseg      = n-k+2
            if len(aarray) != Nseg+1:
                print """*** Error *** make_knot_vector:
                 Length of array (open knot vector) subdividing the parameter
                 range does match order k of curve and number of vertices n+1.
              """
                return -2
            knots = zeros((Nknots), lb3dFLOAT)
            knots[0:k] = aarray[0]
            for j in range(k,n+1):
                knots[j] = aarray[j-k+1]
            knots[n+1:] = aarray[-1]
        elif type == 'periodic':
            Nknots    = n+k+1
            Nintknots = n-k+1
            Nseg      = n-k+2
            if len(aarray) != Nseg+1:
                print """*** Error *** make_knot_vector:
                 Length of array (periodic knot vector) subdividing the
                 parameter range does match order k of curve and number
                 of vertices n+1.
              """
                return -3
            knots = zeros((Nknots), lb3dFLOAT)
            for j in range(k-1,n+2):
                knots[j] = aarray[j-k+1]
            for i in range(k-1):
                j = k-2-i
                knots[j] = knots[k-1] - (knots[2*k-2-j] - knots[k-1])
                j = n+2+i
                knots[j] = knots[n+1] + (knots[n+1] - knots[2*n+2-j])
        elif type == 'closed':
            Nknots    = n+2*k
            Nintknots = n
            Nseg      = n+1
            if len(aarray) != Nseg+1:
                print """*** Error *** make_knot_vector:
                 Length of array (closed knot vector) subdividing the
                 parameter range does match order k of curve and number
                 of vertices n+1.
              """
                return -4
            knots = zeros((Nknots), lb3dFLOAT)
            for j in range(k-1,n+k+1):
                knots[j] = aarray[j-(k-1)]
                for i in range(k-1):
                    j = k-2-i
                    knots[j] = knots[k-1] - (knots[2*k-2-j] - knots[k-1])
                    j = n+k+i+1
                    knots[j] = knots[n+k] + (knots[n+k] - knots[2*(n+k)-j])
        else:
            print '*** Error *** make_knot_vector: unknown B-spline basis type'
            return -5
        return knots
       
    return -1

def make_uniform_knot_vector(k, n, a0, ap, type='open'):
    """
        make_knot_vector

        Routine to make knot vectors for a B-spline basis

        k = order of Basis (or curve, surface), i.e. degree is k-1
            (example: fourth order curve is cubic)
        n = highest vertex index. There are n+1 vertices in total
            Numbers starting with 0.
    """

    if type == 'open':
        Nknots    = n+k+1
        Nintknots = n-k+1
        Nseg      = n-k+2
        knots = zeros((Nknots), lb3dFLOAT)
        knots[0:k] = a0
        for j in range(k,n+1):
            knots[j] = knots[j-1] + (ap-a0)/float(Nseg)
        knots[n+1:] = ap
    elif type == 'periodic':
        Nknots    = n+k+1
        Nintknots = n-k+1
        Nseg      = n-k+2
        knots = zeros((Nknots), lb3dFLOAT)
        knots[k-1] = a0
        for j in range(k,n+2):
            knots[j] = knots[j-1] + (ap-a0)/float(Nseg)
        for i in range(k-1):
            j = k-2-i
            knots[j] = knots[k-1] - (knots[2*k-2-j] - knots[k-1])
            j = n+2+i
            knots[j] = knots[n+1] + (knots[n+1] - knots[2*n+2-j])
    elif type == 'closed':
        Nknots    = n+2*k
        Nintknots = n-k+1
        Nseg      = n+1
        knots = zeros((Nknots), lb3dFLOAT)
        knots[k-1] = a0
        for j in range(k,n+k+1):
            knots[j] = knots[j-1] + (ap-a0)/float(Nseg)
        for i in range(k-1):
            j = k-2-i
            knots[j] = knots[k-1] - (knots[2*k-2-j] - knots[k-1])
            j = n+k+i+1
            knots[j] = knots[n+k] + (knots[n+k] - knots[2*(n+k)-j])
    else:
        print '*** Error *** make_knot_vector: unknown B-spline basis type'
        sys.exit()
    return knots
    
def Nik(knots, i, k, t):
    """
        Bernstein Basis Polynom

        This is a recursive function following the Cox, de Boor
        algorithm.

        knots   list or array containing knot vector.
        i       no of span  
        k       order (degree is k-1)
        t       parameter value (scalar or array)
    """

    try:
        # Works if t is an array or list
        if len(t) > 1:
            tisarray = 1
        else:
            tisarray = 0   # should not be reached, but to make sure
                           # tisarray is initialized
    except:
        # t hs no length attribute , it must be a scalar
        tisarray = 0

    # Operate on arrays
    #if tisarray:
    t = asarray(t)
    #else:
    #    t = array([t])
    if tisarray:
        f = zeros((len(t)),lb3dFLOAT)
    else:
        #print 'fifi'
        f = zeros((1),lb3dFLOAT)
        #print f, f[0] 
    #print "t=", t
    
    # first order Bernstein polynomial
    if k == 1:
        #if tisarray:
            if knots[i+1] == knots[-k]: # make sure that the last segment gets
                                        # a one for the last point
                f=where(array(knots[i]<=t,int)&array(t<=knots[i+1],int),1.,0.)
            else:
                f=where(array(knots[i]<=t,int)&array(t<knots[i+1],int),1.,0.)
        #else:   ## THIS HAVE SEEMED TO BE WRONG !!!! 7.2.03, lb
        #    f = asarray(1.)
    else:
        deltakni  = knots[i+k-1] - knots[i]
        deltakni1 = knots[i+k] - knots[i+1]
         # Check for identical knots
        if abs(deltakni) > 0.000001:
            f = f + (t-knots[i]) / deltakni * Nik(knots, i, k-1, t)
        if abs(deltakni1) > 0.000001:
            f = f + (knots[i+k]-t) / deltakni1 * Nik(knots, i+1, k-1, t)

    #print asarray(f)
    # Correct any rounding errors to keep positivity of Nik
    f = where(f < 0.,0.,f)

    # Return array of values if t came in as array otherwise return
    # scalar
    #print 'Nik ', f
    return f
#     if tisarray:
#         return f
#     else:
#         return f[0]

def dmNikdtm(knots, i, k, t, m=1):
    """
        mth derivative of Bernstein Basis Polynom

        This is a recursive function following the Cox, de Boor
        algorithm.

        knots   list or array containing knot vector.
        i       no of span  
        k       order (degree is k-1)
        t       parameter value (scalar or array)
    """

    try:
        # Works if t is an array or list
        if len(t) > 1:
            tisarray = 1
        else:
            tisarray = 0   # should not be reached, but to make sure
                           # tisarray is initialized
    except:
        # t hs no length attribute , it must be a scalar
        tisarray = 0

    # Operate on arrays
    t = asarray(t)
    f = zeros((len(t)),lb3dFLOAT)

    if m > k-1: # All these higher order derivatives are zero
        if tisarray:
            return f
        else:
            return f[0]
    elif 0 == m: # No derivative at all
        #print 'hier', Nik(knots,i,k,t)
        return Nik(knots,i,k,t)
    elif m < 0:
        print 'stupid error'
        return -1

    # second order Bernstein polynomial derivative
    # has two support intervals and a positive and negative constant value
    deltakni  = knots[i+k-1] - knots[i]
    deltakni1 = knots[i+k] - knots[i+1]
    if k == 1:
        if tisarray:
            if knots[i+1] == knots[-k]: # make sure that the last segment gets
                                        # a one for the last point
                f=where(array(knots[i]<=t,int)&array(t<=knots[i+1],int),0.,0.)
            else:
                f=where(array(knots[i]<=t,int)&array(t<knots[i+1],int),0.,0.)
        else:
            f = asarray(0.)
    else: # higher orders
         # Check for identical knots
        if abs(deltakni) > 0.000001:
            f = f \
              + float(m) /deltakni * dmNikdtm(knots, i, k-1, t, m-1) \
               + (t-knots[i]) / deltakni * dmNikdtm(knots, i, k-1, t, m)

        if abs(deltakni1) > 0.000001:
            f = f \
              - float(m) /deltakni1 * dmNikdtm(knots, i+1, k-1, t, m-1) \
               + (knots[i+k]-t) / deltakni1 * dmNikdtm(knots, i+1, k-1, t, m)

    # Return array of values if t came in as array otherwise return
    # scalar
    ##print 'k= ',k, ' i= ',i, f
    return f
#     if tisarray:
#         return f
#     else:
#         return f[0]

def crossproduct(v1, v2):

    """
       Make cross product of 2D or 3D vectors
    """
    v1 = asarray(v1)
    dim1 = v1.shape[-1]     # dimension of vector
    n1   = len(v1.shape)   # How many vectors
    
    v2 = asarray(v2)
    dim2 = v2.shape[-1]     # dimension of vector
    n2   = len(v2.shape)   # How many vectors
    vector = zeros(v1.shape,lb3dFLOAT)
    if dim1 != dim2:
        print " *** Error **** crossproduct: Vector dimensions do not match!"
        return -1
    if n1 != n2:
        print " *** Error **** crossproduct: confused array length!"
        return -2
    if dim1 != 3 or dim2 != 3:
         print " *** Error **** crossproduct: only 3D vectors!"
       
    vector[...,0] = v1[...,1]*v2[...,2] - v1[...,2]*v2[...,1]
    vector[...,1] = v1[...,2]*v2[...,0] - v1[...,0]*v2[...,2]
    vector[...,2] = v1[...,0]*v2[...,1] - v1[...,1]*v2[...,0]

    return vector

def normalize(vector):

    """
       Normalize vectors to unit length
    """
    vector=asarray(vector)
    n = len(vector.shape)
    dim = vector.shape[-1]
    length = vector_length(vector)
    #print "length"
    #print length
    #print min(length), max(length)
    #eps = array((length.shape[:-1]), lb3dFLOAT)
    #eps[...] = 0.000001
    #print eps
    if n == 1:
        if length > 0.0000001:
            vector = vector/length
    else:
        #for l in length[...]:
        #
        if len(length.shape) == 2:
            for i in range(length.shape[0]):
                for j in range(length.shape[0]):
                    if length[i,j] > 0.0000001:
                        vector[i,j,:] = vector[i,j,:]/length[i,j]
                    else:
                        vector[i,j,:] = 0.
        else:
            for d in range(dim):
                vector[...,d] = where(length[...] > 0.00001, \
                              vector[...,d]/length[...],0. )
    return vector

def vector_length(vector):

    """
       Compute length of vectors, ||v||2
    """
    vector=asarray(vector)
    n = len(vector.shape)
    dim = vector.shape[-1]
    if n == 1:
        # one vector only
        length = 0.
        for d in range(dim):
            length = length + vector[d]*vector[d]
    else:
        length = zeros(vector.shape[:-1], lb3dFLOAT)
        for d in range(dim):
            length = length + vector[...,d]*vector[...,d]
    return sqrt(length)

def euclidean_length(vector):

    """
       Compute length of vectors, ||v||2 using the first three dimensions only
    """
    vector=asarray(vector)
    n = len(vector.shape)
    dim = 3
    if n == 1:
        # one vector only
        length = 0.
        for d in range(dim):
            length = length + vector[d]*vector[d]
    else:
        length = zeros(vector.shape[:-1], lb3dFLOAT)
        for d in range(dim):
            length = length + vector[...,d]*vector[...,d]
    return sqrt(length)

def nurbsmap(vertices,weights=array([])):
    dim = vertices.shape[-1]
    #print "nurbsmap", vertices.shape
    if len(weights) == 0:
        NURBS=0
    else:
        NURBS=1
        #print "nurbsmap",weights.shape
    if NURBS:
        dim = dim+1
        if len(vertices.shape) == 2:
            V = zeros((vertices.shape[0],dim), lb3dFLOAT)
        elif len(vertices.shape) == 3:
            V = zeros((vertices.shape[0],vertices.shape[1],dim), lb3dFLOAT)
        else:
            print "*** Error *** nurbsmap: " \
                  + "Cannot convert to homogenous coordinates!"
            sys.exit(-1)
        #print weights.shape, vertices.shape
        for d in range(dim-1):
            V[...,d] = weights*vertices[...,d]
        V[...,-1] = weights
        #print " V =", V
    else:
        V = vertices

    return NURBS, dim, V

def nurbsremap(dim, r):

    dim = dim-1
    if r.shape > 1:
        rmap = r[...,0:dim] # copy r without the last dimension containing the
                        # weights
        #print r
        #print rmap
        # remapping operation
        for i in range(dim):
            rmap[...,i] = r[...,i]/r[...,-1]
            #rmap[...,i] = where(abs(r[...,-1]) > 1e-12, r[...,i]/r[...,-1], \
            #                 r[...,i])
    else:
        rmap = r[0:dim]
        for i in range(dim):
            rmap[i] = where(abs(r[-1]) > 1e-12, r[i]/r[-1], r[i])

    return dim, rmap

def compute_curve(order, knots, vertices, t,  \
                  weights=array([]), type='open', remap=1):

    #print "compute_curve"
    #print t
    t   = asarray(t)
    NURBS, dim, V = nurbsmap(vertices,weights)
    #print 'dim ', dim
    #dim = vertices.shape[-1]
    #print t
    #print len(t)
    r   = zeros((len(t),dim), float)
    #print '**', r
    
    n   = len(V[:,0]) - 1  # highest vertex index ##TLM removed
    
    #print "weights=",weights
    #denom = zeros((len(t),1), lb3dFLOAT)
    #print "Curve 1"
    if type=='open':
        #print 'Huhu'
        for i in range(len(knots)-order): # 0 to n
            #print "Curve", i, t, Nik(knots,i,order,t)
            N = Nik(knots,i,order,t)
            #denom[...,0] = denom[...,0] + N*weights[i]
            for d in range(dim):
                r[:,d] = V[i,d] * N  + r[:,d] # * weights[i]
    elif type=='closed' or type=='periodic':
        for i in range(n+order): # 0 to n+k-1
            N = Nik(knots,i,order,t)
            #denom[...,0] = denom[...,0] + N*weights[i]
            for d in range(dim):
                r[:,d] = V[i%(n+1),d] * N + r[:,d] # * weights[i]
    else:
        print "ERROR"

    if NURBS and remap:
        dim, r = nurbsremap(dim, r)

    if len(t) > 1:
        return r[...,0:dim]# / denom
    else:
        return r[0,0:dim]# / denom

    
def compute_curve_tangent(order, knots, vertices, t, weights=array([]), \
                          type='open', remap=1, norm=0):

    t   = asarray(t)
    NURBS, dim, V = nurbsmap(vertices,weights)
    #dim = vertices.shape[-1]
    tg  = zeros((len(t),dim), lb3dFLOAT)
    r  = compute_curve(order, knots, vertices, t, weights=weights, \
                       type=type, remap=0)
    n   = len(V[:,0]) - 1  # highest vertex index
    #print "Curve 1"
    #denom = zeros((len(t),1), lb3dFLOAT)

    if type=='open':
        for i in range(len(knots)-order): # 0 to n
            #print "Curve", i, t, Nik(knots,i,order,t)
            #N = Nik(knots,i,order,t)
            Nm = dmNikdtm(knots,i,order,t,m=1)
            #print i, N
            #denom[...,0] = denom[...,0] + N*weights[i]
            for d in range(dim):
                tg[:,d] = tg[:,d] + V[i,d] * Nm
                #r[:,d] = r[:,d] + V[i,d] * N
    elif type=='closed' or type=='periodic':
        for i in range(n+order): # 0 to n+k-1
            #N = Nik(knots,i,order,t)
            Nm = dmNikdtm(knots,i,order,t,m=1)
            for d in range(dim):
                tg[:,d] = tg[:,d] + V[i%(n+1),d] * Nm
                #r[:,d] = r[:,d] + V[i%(n+1),d] * N
    #tg = tg / denom
    # remap r and subtract from tg
    if NURBS:
        for d in range(dim-1):
            r[...,d] = where(abs(r[...,-1]) > 1e-12, r[...,d]/r[...,-1], \
                      r[...,d])
            tg[:,d] = (tg[:,d] - tg[:,-1]*r[:,d])#/r[:,-1]
        tg[:,-1] = r[:,-1]
    #print tg
    if NURBS and remap:
        dim, tg = nurbsremap(dim, tg)

    #print "tgmap = ", tgmap
    if norm:            
        tg = normalize(tg)
    if len(t) > 1:
        return tg[...,0:dim]
    else:
        return tg[0,0:dim]

def compute_curve_binormal(order, knots, vertices, t, weights=array([]), \
                           type='open',\
                           norm=0):

    t      = asarray(t)
    NURBS, dim, V = nurbsmap(vertices,weights)
    r      = compute_curve(order, knots, vertices, t, weights, type, remap=0)
    rp     = compute_curve_tangent(order, knots, vertices, t, weights, type)
    rpp    = zeros((len(t),dim), lb3dFLOAT)
    binormal = zeros((len(t),dim), lb3dFLOAT)

    n   = len(V[:,0]) - 1  # highest vertex index
    #print "Curve 1"
    
    if type=='open':
        for i in range(len(knots)-order): # 0 to n
            #print "Curve", i, t, Nik(knots,i,order,t)
            N = dmNikdtm(knots,i,order,t,m=2)
            for d in range(dim):
                rpp[:,d] = rpp[:,d] + V[i,d] * N
    elif type=='closed' or type=='periodic':
        for i in range(n+order): # 0 to n+k-1
            N = dmNikdtm(knots,i,order,t,m=2)
            for d in range(dim):
                rpp[:,d] = rpp[:,d] + V[i%(n+1),d] * N
    if NURBS:
        for d in range(dim-1):
            r[...,d] = where(abs(r[...,-1]) > 1e-12, r[...,d]/r[...,-1], \
                      r[...,d])
            rpp[:,d] = rpp[:,d] - 2.*rp[:,d] - rpp[:,-1]*r[:,d]
        rpp[:,-1] = r[:,-1]
        dim, rpp = nurbsremap(dim, rpp)
        
    binormal = crossproduct( rp, rpp )
        
    if norm:
        binormal = normalize(binormal)
    
    if len(t) > 1:
        return binormal
    else:
        return binormal[0,:]
    
def compute_curve_normal(order, knots, vertices, t, weights=array([]),\
                         type='open', \
                         norm=0):

    t      = asarray(t)

    tg     = compute_curve_tangent(order, knots, vertices, t, weights, type)
    binormal = compute_curve_binormal(order, knots, vertices, t, weights, type)

    n   = len(vertices[:,0]) - 1  # highest vertex index
    #print "Curve 1"
    
    normal = crossproduct( binormal, tg )
    if norm:
        normal = normalize(normal)
    
    if len(t) > 1:
        return normal
    else:
        return normal[0,:]

def compute_curve_curvature(order, knots, vertices, t, weights=array([]), \
                            type='open'):
    """
       compute_curve_curvature
       Computes curvature of a curve
       
    """

    t     = asarray(t)
    dim = vertices.shape[-1]

    tg    = compute_curve_tangent(order, knots, vertices, t, weights, type)
    speed = vector_length(tg)

    if dim == 2: # 2D curve 
        NURBS, dim, V = nurbsmap(vertices,weights)
        rpp    = zeros((len(t),dim), lb3dFLOAT)

        if NURBS:
            r = compute_curve(order, knots, vertices, t, weights, type,\
                              remap=0)
        n   = len(V[:,0]) - 1  # highest vertex index
    
        if type=='open':
            for i in range(len(knots)-order): # 0 to n
                #print "Curve", i, t, Nik(knots,i,order,t)
                N = dmNikdtm(knots,i,order,t,m=2)
                for d in range(dim):
                    rpp[:,d] = rpp[:,d] + V[i,d] * N
        elif type=='closed' or type=='periodic':
            for i in range(n+order): # 0 to n+k-1
                N = dmNikdtm(knots,i,order,t,m=2)
                for d in range(dim):
                    rpp[:,d] = rpp[:,d] + V[i%(n+1),d] * N
        if NURBS:
            for d in range(dim-1):
                r[...,d] = where(abs(r[...,-1]) > 1e-12, r[...,d]/r[...,-1], \
                      r[...,d])
                rpp[:,d] = rpp[:,d] - 2.*tg[:,d] - rpp[:,-1]*r[:,d]
            rpp[:,-1] = r[:,-1]
            dim, rpp = nurbsremap(dim, rpp)
        if len(t) > 1:
            curvature = where(speed > 0.000001, \
                  (tg[...,0]*rpp[...,1]-tg[...,1]*rpp[...,0])/speed**3, \
                              0. )
        else:
            pass
    else: # 3D
        binormal = compute_curve_binormal(order, knots, vertices, t, \
                                          weights, type)
        binormallength = vector_length(binormal)
        if len(t) > 1:
            curvature = where(speed > 0.000001, \
                          binormallength / speed**3, 0.)
        else:
            if speed > 0.000001:
                curvature = binormallength / speed**3
            else:
                curvature = 0. # zero means not defined
    # all done
    return curvature

def compute_curve_torsion(order, knots, vertices, t, \
                          weights=array([]), type='open'):
    """
       compute_curve_torsion
       Computes torsion of a curve
       
    """
    dim    = vertices.shape[-1]

    if dim == 2: # 2D curve has no torsion
        if len(t) > 1:
            torsion = zeros((len(t)), lb3dFLOAT)
        else:
            torsion = 0.
    else:
        binormal = compute_curve_binormal(order, knots, vertices, t, \
                                          weights, type)
        binormallength = vector_length(binormal)
        NURBS, dim, V = nurbsmap(vertices,weights)
        if NURBS:
            r      = compute_curve(order, knots, vertices, t, weights, \
                               type, remap=0)
            rp     = compute_curve_tangent(order, knots, vertices, t, \
                                           weights, type, remap=0)
            rpp    = zeros((len(t),dim), lb3dFLOAT)
        rppp   = zeros((len(t),dim), lb3dFLOAT)
        if type=='open':
            for i in range(len(knots)-order): # 0 to n
                #print "Curve", i, t, Nik(knots,i,order,t)
                if NURBS:
                    Nm2 = dmNikdtm(knots,i,order,t,m=2)
                Nm3 = dmNikdtm(knots,i,order,t,m=3)
                #print N
                for d in range(dim):
                    if NURBS:
                        rpp[:,d] = rpp[:,d] + V[i,d] * Nm2
                    rppp[:,d] = rppp[:,d] + V[i,d] * Nm3
        elif type=='closed' or type=='periodic':
            for i in range(n+order): # 0 to n+k-1
                if NURBS:
                    Nm2 = dmNikdtm(knots,i,order,t,m=2)
                Nm3 = dmNikdtm(knots,i,order,t,m=3)
                for d in range(dim):
                    if NURBS:
                        rpp[:,d] = rpp[:,d] + V[i%(n+1),d] * Nm2
                    rppp[:,d] = rppp[:,d] + V[i%(n+1),d] * Nm3

        if NURBS:
            for d in range(dim-1):
                r[...,d] = where(abs(r[...,-1]) > 1e-12, r[...,d]/r[...,-1], \
                                 r[...,d])
                rp[...,d] = where(abs(rp[...,-1]) > 1e-12, \
                                  rp[...,d]/r[...,-1], \
                                  rp[...,d])
                rpp[...,d] = where(abs(rpp[...,-1]) > 1e-12, \
                                  rpp[...,d]/r[...,-1], \
                                  rpp[...,d])
                rppp[:,d] = rppp[:,d] - 3.* rp[:,-1]*rpp[:,d] \
                            - 3.* rpp[:,-1]*rp[:,d] - rppp[:,-1]*r[:,d]
            rppp[:,-1] = r[:,-1]
            dim, rppp = nurbsremap(dim, rppp)

        #print 'rppp = ',rppp
        if len(t) > 1:
            M = zeros(binormal.shape[:-1], lb3dFLOAT)
            #print binormal.shape, transpose(rppp).shape
            for i in range(len(M)):
                M[i] = dot(binormal[i,0:dim],rppp[i,0:dim])
            #print M
            torsion = where(binormallength > 0.000001, \
                      M / binormallength**2,0.)
        else:
            # only one vector
            if binormallength > 0.000001:
                torsion = dot(binormal,rppp) / binormallength**2
            else:
                torsion = 0.  # zero means not defined
                
    # all done
    return torsion


    
def compute_surf(uorder,vorder, knotsu, knotsv, vertices, u, v, \
                 weights=array([]), remap=1):
    """
      remap=1 return 2D/3D coordinates
      remap=0 return homogenous coordinates
    """

    try:
        _itmp = len(u)
        u = asarray(u)   # u is a list or array
    except:
        u = asarray([u]) # u is a scalar
    try:
        _itmp = len(v)
        v = asarray(v)   # v is a list or array
    except:
        v = asarray([v]) # v is a scalar
    #print "u,v=",u,v
    #dim = vertices.shape[-1]
    NURBS, dim, V = nurbsmap(vertices,weights)
    r = zeros((len(u),len(v),dim), lb3dFLOAT)
        
    if len(shape(knotsu))==2 and len(shape(knotsv))==2:
        for j in range(len(knotsv[0,:])-vorder):
            for i in range(len(knotsu[0,:])-uorder):
                Nu = Nik(knotsu[j,:],i,uorder,u)
                Nv = Nik(knotsv[i,:],j,vorder,v)
                M = outer(Nu,Nv)
                for d in range(dim):
                    r[:,:,d] = r[:,:,d]+matrixmultiply(V[i,j,d],M)
    elif len(shape(knotsu))==2 and len(shape(knotsv))==1:
        for j in range(len(knotsv)-vorder):
            Nv = Nik(knotsv,j,vorder,v)
            for i in range(len(knotsu[0,:])-uorder):
                Nu = Nik(knotsu[j,:],i,uorder,u)
                M = outer(Nu,Nv)
                for d in range(dim):
                    r[:,:,d] = r[:,:,d]+matrixmultiply(V[i,j,d],M)
    elif len(shape(knotsv))==2 and len(shape(knotsu))==1:
        for j in range(len(knotsv[0,:])-vorder):
            for i in range(len(knotsu)-uorder):
                Nu = Nik(knotsu,i,uorder,u)
                Nv = Nik(knotsv[i,:],j,vorder,v)
                M = outer(Nu,Nv)
                for d in range(dim):
                    r[:,:,d] = r[:,:,d]+matrixmultiply(V[i,j,d],M)
    else:
        for j in range(len(knotsv)-vorder):
            Nv = Nik(knotsv,j,vorder,v)
            for i in range(len(knotsu)-uorder):
                Nu = Nik(knotsu,i,uorder,u)
                M = outer(Nu,Nv)
                #print M.getshape()
                #print r.getshape()
                #print V.getshape()
                #MM = zeros((Nu.getshape()[0], Nv.getshape()[0],dim),"f")
                #for d in range(dim):
                #    MM[:,:,d] = M
                #print "Nu, Nv =", Nu, Nv

                # Problem here: caused by change from Numeric to numarray
                
                for d in range(dim):
                #    r[:,:,d] = r[:,:,d]+matrixmultiply(V[i,j,d],M)
                    r[:,:,d] = r[:,:,d] + V[i,j,d]*M
                
    if NURBS and remap:
        dim, r = nurbsremap(dim, r)
        
    if len(u) == 1 and len(v) == 1:
        return r[0,0,:]
    else:
        return r
    #return r

def compute_surf_utangent(uorder, vorder, knotsu, knotsv, \
                          vertices, u, v, weights=array([]), norm=0, \
                          remap=1):

    try:
        _itmp = len(u)
        u = asarray(u)   # u is a list or array
    except:
        u = asarray([u]) # u is a scalar
    try:
        _itmp = len(v)
        v = asarray(v)   # v is a list or array
    except:
        v = asarray([v]) # v is a scalar

    #dim = vertices.shape[-1]
    NURBS, dim, V = nurbsmap(vertices,weights)
    if NURBS:
        r = compute_surf(uorder, vorder, knotsu, knotsv, \
                          vertices, u, v, weights=weights, remap=0)
    ru = zeros((len(u),len(v),dim), lb3dFLOAT)
        
    if len(shape(knotsu))==2 and len(shape(knotsv))==2:
        for j in range(len(knotsv[0,:])-vorder):
            for i in range(len(knotsu[0,:])-uorder):
                Nu = dmNikdtm(knotsu[j,:],i,uorder,u,m=1)
                Nv = Nik(knotsv[i,:],j,vorder,v)
                M = outer(Nu,Nv)
                for d in range(dim):
                    ru[:,:,d] = ru[:,:,d]+matrixmultiply(V[i,j,d],M)
    elif len(shape(knotsu))==2 and len(shape(knotsv))==1:
        for j in range(len(knotsv)-vorder):
            Nv = Nik(knotsv,j,vorder,v)
            for i in range(len(knotsu[0,:])-uorder):
                Nu = dmNikdtm(knotsu[j,:],i,uorder,u,m=1)
                M = outer(Nu,Nv)
                for d in range(dim):
                    ru[:,:,d] = ru[:,:,d]+matrixmultiply(V[i,j,d],M)
    elif len(shape(knotsv))==2 and len(shape(knotsu))==1:
        for j in range(len(knotsv[0,:])-vorder):
            for i in range(len(knotsu)-uorder):
                Nu = dmNikdtm(knotsu,i,uorder,u,m=1)
                Nv = Nik(knotsv[i,:],j,vorder,v)
                M = outer(Nu,Nv)
                for d in range(dim):
                    ru[:,:,d] = ru[:,:,d]+matrixmultiply(V[i,j,d],M)
    else:
        for j in range(len(knotsv)-vorder):
            Nv = Nik(knotsv,j,vorder,v)
            for i in range(len(knotsu)-uorder):
                Nu = dmNikdtm(knotsu,i,uorder,u,m=1)
                M = outer(Nu,Nv)
                for d in range(dim):
                    #ru[:,:,d] = ru[:,:,d]+matrixmultiply(V[i,j,d],M)
                    ru[:,:,d] = ru[:,:,d] + V[i,j,d]*M

    #print r
    if NURBS:
        for d in range(dim-1):
            r[...,d] = r[...,d]/r[...,-1]
            #r[...,d] = where(abs(r[...,-1]) > MINFLOAT, r[...,d]/r[...,-1],
            #          MINFLOAT ) #r[...,d])
            ru[...,d] = ru[...,d] - ru[...,-1]*r[...,d]
        ru[...,-1] = r[...,-1]

    #print ru
    if NURBS and remap:
        dim, ru = nurbsremap(dim, ru)
        
    if norm and remap:
        ru = normalize(ru)
        
#     # Check for degenerated vectors
#     index = zeros((len(u),len(v)), int)
#     index = where(euclidean_length(ru) < MINFLOAT, 1, 0)
#     #print index
#     #shift = 0.0000001 # !!!!!!!!!!!!!!!!!!!!!! Should be controlled somehow
#     for i in range(len(u)):
#         for j in range(len(v)):
#             if index[i,j]==1:
#                 vs = v[j] - shift
#                 if  vs < 0.0:
#                     vs = v[j] + shift
#                 ru[i,j,:] = compute_surf_utangent(uorder, vorder, \
#                                                   knotsu, knotsv, \
#                           vertices, u[i], vs, weights=weights,\
#                                                   norm=norm, remap=remap)

    if len(u) == 1 and len(v) == 1:
        return ru[0,0,:]
    else:
        return ru
    #return ru

def compute_surf_uu(uorder, vorder, knotsu, knotsv, \
                    vertices, u, v, weights=array([]), \
                    norm=0, remap=1 ):

    try:
        _itmp = len(u)
        u = asarray(u)   # u is a list or array
    except:
        u = asarray([u]) # u is a scalar
    try:
        _itmp = len(v)
        v = asarray(v)   # v is a list or array
    except:
        v = asarray([v]) # v is a scalar

    #dim = vertices.shape[-1]
    NURBS, dim, V = nurbsmap(vertices,weights)
    if NURBS:
        r = compute_surf(uorder, vorder, knotsu, knotsv, \
                          vertices, u, v, weights, remap=0)
        ru = compute_surf_utangent(uorder, vorder, knotsu, knotsv, \
                          vertices, u, v, weights, norm=0, remap=0)
    ruu = zeros((len(u),len(v),dim), lb3dFLOAT)
        
    if len(shape(knotsu))==2 and len(shape(knotsv))==2:
        for j in range(len(knotsv[0,:])-vorder):
            for i in range(len(knotsu[0,:])-uorder):
                Nu = dmNikdtm(knotsu[j,:],i,uorder,u,m=2)
                Nv = Nik(knotsv[i,:],j,vorder,v)
                M = outer(Nu,Nv)
                for d in range(dim):
                    ruu[:,:,d] = ruu[:,:,d]+matrixmultiply(V[i,j,d],M)
    elif len(shape(knotsu))==2 and len(shape(knotsv))==1:
        for j in range(len(knotsv)-vorder):
            Nv = Nik(knotsv,j,vorder,v)
            for i in range(len(knotsu[0,:])-uorder):
                Nu = dmNikdtm(knotsu[j,:],i,uorder,u,m=2)
                M = outer(Nu,Nv)
                for d in range(dim):
                     ruu[:,:,d] = ruu[:,:,d]+matrixmultiply(V[i,j,d],M)
    elif len(shape(knotsv))==2 and len(shape(knotsu))==1:
        for j in range(len(knotsv[0,:])-vorder):
            for i in range(len(knotsu)-uorder):
                Nu = dmNikdtm(knotsu,i,uorder,u,m=2)
                Nv = Nik(knotsv[i,:],j,vorder,v)
                M = outer(Nu,Nv)
                for d in range(dim):
                    ruu[:,:,d] = ruu[:,:,d]+matrixmultiply(V[i,j,d],M)
    else:
        for j in range(len(knotsv)-vorder):
            Nv = Nik(knotsv,j,vorder,v)
            for i in range(len(knotsu)-uorder):
                Nu = dmNikdtm(knotsu,i,uorder,u,m=2)
                M = outer(Nu,Nv)
                for d in range(dim):
                    #ruu[:,:,d] = ruu[:,:,d]+matrixmultiply(V[i,j,d],M)
                    ruu[:,:,d] = ruu[:,:,d] + V[i,j,d]*M

    if len(u) == 1 and len(v) == 1:
        ruu = ruu[0,0,:]

    if NURBS:
        for d in range(dim-1):
            r[...,d] = r[...,d]/r[...,-1]
            #r[...,d] = where(abs(r[...,-1]) > 1e-12, r[...,d]/r[...,-1], \
            #          r[...,d])
            #print r.shape
            #print ru.shape
            #print ruu.shape
            if len(r.shape) > 1:
                ruu[...,d] = ruu[...,d] - 2.*ru[...,d] - ruu[...,-1]*r[...,d]
            else:
                ruu[d] = ruu[d] - 2.*ru[d] - ruu[-1]*r[d]
                
        ruu[...,-1] = r[...,-1]        
    if NURBS and remap:
        dim, ruu = nurbsremap(dim, ruu)
        
    if norm and remap:
        ruu = normalize(ruu)

#     # Check for degenerated vectors
#     index = zeros((len(u),len(v)), int)
#     index = where(euclidean_length(ruu) < MINFLOAT, 1, 0)
#     print "index ", index.getshape()
#     #print ruu
#     #shift = 0.0000001 # !!!!!!!!!!!!!!!!!!!!!! Should be controlled somehow
#     if len(index.getshape()) > 1:
#         for i in range(len(u)):
#             for j in range(len(v)):
#                 if index[i,j]==1:
#                     vs = v[j] - shift
#                     if  vs < 0.0:
#                         vs = v[j] + shift
#                         ruu[i,j,:] = compute_surf_uu(uorder, vorder, \
#                                                   knotsu, knotsv, \
#                           vertices, u[i], vs, weights=weights,\
#                                                   norm=norm, remap=remap)
#     else:
#         vs = v[0] - shift
#         if vs < 0.0:
#             vs = v[0] + shift
#         ruu[0,0,:] = compute_surf_uu(uorder, vorder, \
#                                                   knotsu, knotsv, \
#                           vertices, u[0], vs, weights=weights,\
#                                                   norm=norm, remap=remap)
        
    if len(u) == 1 and len(v) == 1:
        return ruu[0,0,:]
    else:
        return ruu
#    return ruu

def compute_surf_uv(uorder, vorder, knotsu, knotsv, \
                    vertices, u, v, weights=array([]), \
                    norm=0, remap=1 ):

    try:
        _itmp = len(u)
        u = asarray(u)   # u is a list or array
    except:
        u = asarray([u]) # u is a scalar
    try:
        _itmp = len(v)
        v = asarray(v)   # v is a list or array
    except:
        v = asarray([v]) # v is a scalar

    #dim = vertices.shape[-1]
    NURBS, dim, V = nurbsmap(vertices,weights)
    if NURBS:
        r = compute_surf(uorder, vorder, knotsu, knotsv, \
                          vertices, u, v, weights, remap=0)
        ru = compute_surf_utangent(uorder, vorder, knotsu, knotsv, \
                          vertices, u, v, weights, norm=0, remap=0)
        rv = compute_surf_vtangent(uorder, vorder, knotsu, knotsv, \
                          vertices, u, v, weights, norm=0, remap=0)
    ruv = zeros((len(u),len(v),dim), lb3dFLOAT)
        
    if len(shape(knotsu))==2 and len(shape(knotsv))==2:
        for j in range(len(knotsv[0,:])-vorder):
            for i in range(len(knotsu[0,:])-uorder):
                Nu = dmNikdtm(knotsu[j,:],i,uorder,u,m=1)
                Nv = dmNikdtm(knotsv[i,:],j,vorder,v,m=1)
                M = outer(Nu,Nv)
                for d in range(dim):
                    ruv[:,:,d] = ruv[:,:,d]+matrixmultiply(V[i,j,d],M)
    elif len(shape(knotsu))==2 and len(shape(knotsv))==1:
        for j in range(len(knotsv)-vorder):
            Nv = dmNikdtm(knotsv,j,vorder,v,m=1)
            for i in range(len(knotsu[0,:])-uorder):
                Nu = dmNikdtm(knotsu[j,:],i,uorder,u,m=1)
                M = outer(Nu,Nv)
                for d in range(dim):
                    ruv[:,:,d] = ruv[:,:,d]+matrixmultiply(V[i,j,d],M)
    elif len(shape(knotsv))==2 and len(shape(knotsu))==1:
        for j in range(len(knotsv[0,:])-vorder):
            for i in range(len(knotsu)-uorder):
                Nu = dmNikdtm(knotsu,i,uorder,u,m=1)
                Nv = dmNikdtm(knotsv[i,:],j,vorder,v,m=1)
                M = outer(Nu,Nv)
                for d in range(dim):
                    ruv[:,:,d] = ruv[:,:,d]+matrixmultiply(V[i,j,d],M)
    else:
        for j in range(len(knotsv)-vorder):
            Nv = dmNikdtm(knotsv,j,vorder,v,m=1)
            for i in range(len(knotsu)-uorder):
                Nu = dmNikdtm(knotsu,i,uorder,u,m=1)
                M = outer(Nu,Nv)
                for d in range(dim):
                    #ruv[:,:,d] = ruv[:,:,d]+matrixmultiply(V[i,j,d],M)
                    ruv[:,:,d] = ruv[:,:,d] + V[i,j,d]*M

    if NURBS:
        for d in range(dim-1):
            r[...,d] = where(abs(r[...,-1]) > MINFLOAT, \
                             r[...,d]/r[...,-1], r[...,d])
            ru[...,d] = where(abs(ru[...,-1]) > MINFLOAT, \
                              ru[...,d]/r[...,-1], ru[...,d])
            rv[...,d] = where(abs(rv[...,-1]) > MINFLOAT, \
                              rv[...,d]/r[...,-1], rv[...,d])
            ruv[...,d] = ruv[...,d] - rv[...,-1]*ru[...,d] \
                                         - ru[...,-1]*rv[...,d] \
                         - ruv[...,-1]*r[...,d]
        ruv[...,-1] = r[...,-1]
        
    if NURBS and remap:
        dim, ruv = nurbsremap(dim, ruv)

    if norm and remap:
        ruv = normalize(ruv)
#     # Check for degenerated vectors
#     index = zeros((len(u),len(v)), int)
#     index = where(euclidean_length(ruv) < MINFLOAT, 1, 0)
#     #print index
#     #shift = 0.0000001 # !!!!!!!!!!!!!!!!!!!!!! Should be controlled somehow
#     for i in range(len(u)):
#         for j in range(len(v)):
#             if index[i,j]==1:
#                 us = u[i] - shift
#                 vs = v[j] - shift
#                 if  vs < 0.0:
#                     vs = v[j] + shift
#                 if  us < 0.0:
#                     us = u[i] + shift
#                 ruv[i,j,:] = compute_surf_uv(uorder, vorder, \
#                                                   knotsu, knotsv, \
#                           vertices, us, vs, weights=weights,\
#                                                   norm=norm, remap=remap)

    if len(u) == 1 and len(v) == 1:
        return ruv[0,0,:]
    else:
        return ruv
#    return ruv

def compute_surf_vtangent(uorder,vorder, knotsu, knotsv, \
                          vertices, u, v, weights=array([]), \
                          norm=0, remap=1 ):

    try:
        _itmp = len(u)
        u = asarray(u)   # u is a list or array
    except:
        u = asarray([u]) # u is a scalar
    try:
        _itmp = len(v)
        v = asarray(v)   # v is a list or array
    except:
        v = asarray([v]) # v is a scalar

    #dim = vertices.shape[-1]
    NURBS, dim, V = nurbsmap(vertices,weights)
    if NURBS:
        r = compute_surf(uorder, vorder, knotsu, knotsv, \
                          vertices, u, v, weights=weights, remap=0)
    rv = zeros((len(u),len(v),dim), lb3dFLOAT)
        
    if len(shape(knotsu))==2 and len(shape(knotsv))==2:
        for j in range(len(knotsv[0,:])-vorder):
            for i in range(len(knotsu[0,:])-uorder):
                Nu = Nik(knotsu[j,:],i,uorder,u)
                Nv = dmNikdtm(knotsv[i,:],j,vorder,v,m=1)
                M = outer(Nu,Nv)
                for d in range(dim):
                    rv[:,:,d] = rv[:,:,d]+matrixmultiply(V[i,j,d],M)
    elif len(shape(knotsu))==2 and len(shape(knotsv))==1:
        for j in range(len(knotsv)-vorder):
            Nv = dmNikdtm(knotsv,j,vorder,v,m=1)
            for i in range(len(knotsu[0,:])-uorder):
                Nu = Nik(knotsu[j,:],i,uorder,u)
                M = outer(Nu,Nv)
                for d in range(dim):
                    rv[:,:,d] = rv[:,:,d]+matrixmultiply(V[i,j,d],M)
    elif len(shape(knotsv))==2 and len(shape(knotsu))==1:
        for j in range(len(knotsv[0,:])-vorder):
            for i in range(len(knotsu)-uorder):
                Nu = Nik(knotsu,i,uorder,u)
                Nv = dmNikdtm(knotsv[i,:],j,vorder,v,m=1)
                M = outer(Nu,Nv)
                for d in range(dim):
                    rv[:,:,d] = rv[:,:,d]+matrixmultiply(V[i,j,d],M)
    else:
        for j in range(len(knotsv)-vorder):
            Nv = dmNikdtm(knotsv,j,vorder,v,m=1)
            for i in range(len(knotsu)-uorder):
                Nu = Nik(knotsu,i,uorder,u)
                M = outer(Nu,Nv)
                for d in range(dim):
                    #rv[:,:,d] = rv[:,:,d]+matrixmultiply(V[i,j,d],M)
                    rv[:,:,d] = rv[:,:,d] + V[i,j,d]*M
                    
    if NURBS:
        for d in range(dim-1):
            r[...,d] = where(abs(r[...,-1]) > 1e-12, r[...,d]/r[...,-1], \
                      r[...,d])
            rv[...,d] = rv[...,d] - rv[...,-1]*r[...,d]
        rv[...,-1] = r[...,-1]

    if NURBS and remap:
        dim, rv = nurbsremap(dim, rv)
        
    if norm and remap:
        rv = normalize(rv)
        
#     # Check for degenerated vectors
#     index = zeros((len(u),len(v)), int)
#     index = where(euclidean_length(rv) < MINFLOAT, 1, 0)
#     #print index
#     #shift = 0.0000001 # !!!!!!!!!!!!!!!!!!!!!! Should be controlled somehow
#     for i in range(len(u)):
#         for j in range(len(v)):
#             if index[i,j]==1:
#                 us = u[i] - shift
#                 if  us < 0.0:
#                     us = u[i] + shift
#                 rv[i,j,:] = compute_surf_vtangent(uorder, vorder, \
#                                                   knotsu, knotsv, \
#                           vertices, us, v[j], weights=weights,\
#                                                   norm=norm, remap=remap)
    if len(u) == 1 and len(v) == 1:
        return rv[0,0,:]
    else:
        return rv
#    return rv

def compute_surf_vv(uorder,vorder, knotsu, knotsv, \
                    vertices, u, v, weights=array([]), \
                    norm=0, remap=1 ):

    try:
        _itmp = len(u)
        u = asarray(u)   # u is a list or array
    except:
        u = asarray([u]) # u is a scalar
    try:
        _itmp = len(v)
        v = asarray(v)   # v is a list or array
    except:
        v = asarray([v]) # v is a scalar

    #dim = vertices.shape[-1]
    NURBS, dim, V = nurbsmap(vertices,weights)
    if NURBS:
        r = compute_surf(uorder, vorder, knotsu, knotsv, \
                          vertices, u, v, weights, remap=0)
        rv = compute_surf_vtangent(uorder, vorder, knotsu, knotsv, \
                          vertices, u, v, weights, norm=0, remap=0)
    rvv = zeros((len(u),len(v),dim), lb3dFLOAT)
        
    if len(shape(knotsu))==2 and len(shape(knotsv))==2:
        for j in range(len(knotsv[0,:])-vorder):
            for i in range(len(knotsu[0,:])-uorder):
                Nu = Nik(knotsu[j,:],i,uorder,u)
                Nv = dmNikdtm(knotsv[i,:],j,vorder,v,m=2)
                M = outer(Nu,Nv)
                for d in range(dim):
                    rvv[:,:,d] = rvv[:,:,d]+matrixmultiply(V[i,j,d],M)
    elif len(shape(knotsu))==2 and len(shape(knotsv))==1:
        for j in range(len(knotsv)-vorder):
            Nv = dmNikdtm(knotsv,j,vorder,v,m=2)
            for i in range(len(knotsu[0,:])-uorder):
                Nu = Nik(knotsu[j,:],i,uorder,u)
                M = outer(Nu,Nv)
                for d in range(dim):
                    rvv[:,:,d] = rvv[:,:,d]+matrixmultiply(V[i,j,d],M)
    elif len(shape(knotsv))==2 and len(shape(knotsu))==1:
        for j in range(len(knotsv[0,:])-vorder):
            for i in range(len(knotsu)-uorder):
                Nu = Nik(knotsu,i,uorder,u)
                Nv = dmNikdtm(knotsv[i,:],j,vorder,v,m=2)
                M = outer(Nu,Nv)
                for d in range(dim):
                    rvv[:,:,d] = rvv[:,:,d]+matrixmultiply(V[i,j,d],M)
    else:
        for j in range(len(knotsv)-vorder):
            Nv = dmNikdtm(knotsv,j,vorder,v,m=2)
            for i in range(len(knotsu)-uorder):
                Nu = Nik(knotsu,i,uorder,u)
                M = outer(Nu,Nv)
                for d in range(dim):
                    #rvv[:,:,d] = rvv[:,:,d]+matrixmultiply(V[i,j,d],M)
                    rvv[:,:,d] = rvv[:,:,d] + V[i,j,d]*M
    if NURBS:
        for d in range(dim-1):
            r[...,d] = where(abs(r[...,-1]) > MINFLOAT, r[...,d]/r[...,-1], \
                      r[...,d])
            rvv[...,d] = rvv[...,d] - 2.*rv[...,d] - rvv[...,-1]*r[...,d]
        rvv[...,-1] = r[...,-1]
        
    if NURBS and remap:
        dim, rvv = nurbsremap(dim, rvv)
        
    if norm and remap:
        rvv = normalize(rvv)
        
#     # Check for degenerated vectors
#     index = zeros((len(u),len(v)), int)
#     index = where(euclidean_length(rvv) < MINFLOAT, 1, 0)
#     #print index
#     #shift = 0.0000001 # !!!!!!!!!!!!!!!!!!!!!! Should be controlled somehow
#     for i in range(len(u)):
#         for j in range(len(v)):
#             if index[i,j]==1:
#                 us = u[i] - shift
#                 if  us < 0.0:
#                     us = u[i] + shift
#                 rvv[i,j,:] = compute_surf_vv(uorder, vorder, \
#                                                   knotsu, knotsv, \
#                           vertices, us, v[j], weights=weights,\
#                                                   norm=norm, remap=remap)
    if len(u) == 1 and len(v) == 1:
        return rvv[0,0,:]
    else:
        return rvv
#    return rvv

def compute_surf_normal(uorder, vorder, uknots, vknots, \
                        vertices, u, v, weights=array([]), norm=0 ):

    try:
        _itmp = len(u)
        u = asarray(u)   # u is a list or array
    except:
        u = asarray([u]) # u is a scalar
    try:
        _itmp = len(v)
        v = asarray(v)   # v is a list or array
    except:
        v = asarray([v]) # v is a scalar

    ru = compute_surf_utangent(uorder, vorder, \
                         uknots, vknots, \
                         vertices, u, v, weights, norm=0)
    rv = compute_surf_vtangent(uorder, vorder, \
                         uknots, vknots, \
                         vertices, u, v, weights, norm=0)
    nn = crossproduct(ru,rv)

    if norm:
        nn = normalize(nn)

    if len(u) == 1 and len(v) == 1:
        return nn[0,0,:]
    else:
        return nn

def compute_gaussian_curvature(uorder, vorder, uknots, vknots, \
                        vertices, u, v, weights=array([])):

    ## TESTED with mbag_hood, 24.10.02, lb and sphere
    
    try:
        _itmp = len(u)
        u = asarray(u)   # u is a list or array
    except:
        u = asarray([u]) # u is a scalar
    try:
        _itmp = len(v)
        v = asarray(v)   # v is a list or array
    except:
        v = asarray([v]) # v is a scalar

    dim = vertices.shape[-1]
    ru = compute_surf_utangent(uorder, vorder, \
                         uknots, vknots, \
                         vertices, u, v, weights, norm=0)
    #print "ru"
    #print ru
    rv = compute_surf_vtangent(uorder, vorder, \
                         uknots, vknots, \
                         vertices, u, v, weights, norm=0)
    #print "rv"
    #print rv
    ruu = compute_surf_uu(uorder, vorder, \
                         uknots, vknots, \
                         vertices, u, v, weights, norm=0)
    #print "ruu"
    #print ruu
    ruv = compute_surf_uv(uorder, vorder, \
                         uknots, vknots, \
                         vertices, u, v, weights, norm=0)
    #print "ruv"
    #print ruv
    rvv = compute_surf_vv(uorder, vorder, \
                         uknots, vknots, \
                         vertices, u, v, weights, norm=0)
    #print "rvv"
    #print rvv
    #print rvv
    nn = crossproduct(ru,rv)
    nn = normalize(nn)   # This seems to be the correct way
    
    Kcurvature = zeros((len(u),len(v)), lb3dFLOAT)
    for i in range(len(u)):
        for j in range(len(v)):
             E = dot(ru[i,j,0:dim],ru[i,j,0:dim])
             F = dot(ru[i,j,0:dim],rv[i,j,0:dim])
             G = dot(rv[i,j,0:dim],rv[i,j,0:dim])
             e = dot(ruu[i,j,0:dim],nn[i,j,0:dim])
             f = dot(ruv[i,j,0:dim],nn[i,j,0:dim])
             g = dot(rvv[i,j,0:dim],nn[i,j,0:dim])
             #print E,F,G, e,f,g
             #print i,j, e,f,g, e*g - f**2, E*G - F**2
             if abs(E*G - F**2) > MINFLOAT:
                 Kcurvature[i,j] = (e*g - f**2) / (E*G - F**2)
             else:
                 Kcurvature[i,j] = 0.
    if len(u) == 1 and len(v) == 1:
        return Kcurvature[0,0]
    else:
        return Kcurvature
    #return Kcurvature
    
def compute_mean_curvature(uorder, vorder, uknots, vknots, \
                        vertices, u, v, weights=array([])):

    ## CORRECT FOR SPHERE, 14.03.02, lb
    
    try:
        _itmp = len(u)
        u = asarray(u)   # u is a list or array
    except:
        u = asarray([u]) # u is a scalar
    try:
        _itmp = len(v)
        v = asarray(v)   # v is a list or array
    except:
        v = asarray([v]) # v is a scalar

    dim = vertices.shape[-1]
    ru = compute_surf_utangent(uorder, vorder, \
                         uknots, vknots, \
                         vertices, u, v, weights, norm=0)
    rv = compute_surf_vtangent(uorder, vorder, \
                         uknots, vknots, \
                         vertices, u, v, weights, norm=0)
    ruu = compute_surf_uu(uorder, vorder, \
                         uknots, vknots, \
                         vertices, u, v, weights, norm=0)
    ruv = compute_surf_uv(uorder, vorder, \
                         uknots, vknots, \
                         vertices, u, v, weights, norm=0)
    rvv = compute_surf_vv(uorder, vorder, \
                         uknots, vknots, \
                         vertices, u, v, weights, norm=0)
    #print rvv
    nn = crossproduct(ru,rv)
    nn = normalize(nn)
    
    Hcurvature = zeros((len(u),len(v)), lb3dFLOAT)
    for i in range(len(u)):
        for j in range(len(v)):
             E = dot(ru[i,j,0:dim],ru[i,j,0:dim])
             F = dot(ru[i,j,0:dim],rv[i,j,0:dim])
             G = dot(rv[i,j,0:dim],rv[i,j,0:dim])
             e = dot(ruu[i,j,0:dim],nn[i,j,0:dim])
             f = dot(ruv[i,j,0:dim],nn[i,j,0:dim])
             g = dot(rvv[i,j,0:dim],nn[i,j,0:dim])
             #print E,F,G, e,f,g
             if abs(E*G - F**2) > 0.000001:
                 # The formula was corrected Oct 28, 2006, There is a
                 # print error in Nowacki's book.
                 Hcurvature[i,j] = 0.5 * (E*g - 2.0*f*F +e*G) / (E*G - F**2)
             else:
                 Hcurvature[i,j] = 0.
    if len(u) == 1 and len(v) == 1:
        return Hcurvature[0,0]
    else:
        return Hcurvature
    #return Hcurvature
        

###----------------------------------------------------------------------------
### 2 Scale Relation for Basis Functions
###----------------------------------------------------------------------------
import scipy as sp
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
        
###----------------------------------------------------------------------------
### Tests
###----------------------------------------------------------------------------
      
def run_test():
    from time import clock, time
    t = arrayrange(0.,1.01,0.1) #error
    #t = 0.5
    knots=(0.,0.,0.,0.5,1.,1.,1.)
    order=3
    a =array(knots[3] <= t,int)
    c= array(t< knots[4],int)
    print a
    print c
    print a & c
    print equal(t,knots[-1])
    b = Nik(knots,3,order,1.0)
    print b
    return

def test_Bernstein():
    n=3
    ts=[]
    
    ti = clock()
    for u in np.linspace(0.,1.1,11,endpoint=False):
        print AllBernstein(n,u)
    te = clock()
    ts.append(te-ti)
    
    ti = clock()
    for u in np.linspace(0.,1.1,11,endpoint=False):
        print '[{},{},{},{}]'.format(Bernstein(0,n,u),
               Bernstein(1,n,u),
                Bernstein(2,n,u),
                Bernstein(3,n,u))
    te = clock()
    ts.append(te-ti)
    
    ti = clock()
    for u in np.linspace(0.,1.1,11,endpoint=False):
        print '[{},{},{},{}]'.format(checkstein(0,n,u),
               checkstein(1,n,u),
                checkstein(2,n,u),
                checkstein(3,n,u))
    te = clock()
    ts.append(te-ti)
           
    print 'times'
    print 'P&T All: {} sec'.format(ts[0])
    print 'P&T indi: {} sec'.format(ts[1])
    print 'scipy binomial: {} sec'.format(ts[2])
    
    return
    
###----------------------------------------------------------------------------
### End
###----------------------------------------------------------------------------

if __name__ == '__main__':
    test_Bernstein()

    
    
