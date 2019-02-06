
# Luke McCulloch
# Using Dr. Birk's code
# to write a uniform knot vecot
import numpy as np
from numpy import *
import sys
import string


#######################################################################

def make_knot_vector(k, n, aarray, type='open'):
    
    """
        make_knot_vector

        Routine to make knot vectors for a B-spline basis

        k = order of Basis (or curve, surface), i.e. degree is k-1
            (example: fourth order curve is cubic)
        n = highest vertex index. There are n+1 vertices in total
            Numbers starting with 0.
    """
    #aarray = np.asarray([0.,1.])
    #print aarray
    if len(aarray) == 2:
        return make_uniform_knot_vector(k, n, aarray[0], aarray[1], type)
    else:
        print 'calling make knot vector'
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
            knots = zeros((Nknots), float)
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
            knots = zeros((Nknots), float)
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
            knots = zeros((Nknots), float)
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


######################################################################

def make_uniform_knot_vector(k, n, a0, ap, type='open'):
    
    """
        make_knot_vector

        Routine to make knot vectors for a B-spline basis

        k = order of Basis (or curve, surface), i.e. degree is k-1
            (example: fourth order curve is cubic)
        n = highest vertex index. There are n+1 vertices in total
            Numbers starting with 0.
    """
    #print 'calling make uniform knot vector'
    if type == 'open':
        Nknots    = n+k+1
        Nintknots = n-k+1
        Nseg      = n-k+2
        #knots = zeros((Nknots), lb3dFLOAT)
        knots = zeros((Nknots), float)
        knots[0:k] = a0
        for j in range(k,n+1):
            knots[j] = knots[j-1] + (ap-a0)/float(Nseg)
        knots[n+1:] = ap
    elif type == 'periodic':
        Nknots    = n+k+1
        Nintknots = n-k+1
        Nseg      = n-k+2
        #knots = zeros((Nknots), lb3FLOAT)
        knots = zeros((Nknots), float)
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
        #knots = zeros((Nknots), lb3dFLOAT)
        knots = zeros((Nknots), float)
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

