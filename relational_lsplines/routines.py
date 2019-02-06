# 20130628
# TLM resturctured codes
# routines for ADLspline
##
## contains:
##           crossproduct
##           ADcrossproduct

import numpy as np
from knot_operations import *
import copy
lb3dFLOAT = float

def crossproduct(v1, v2):

    """
    Dr. Birk's Cross Product Function
    """
    v1 = np.asarray(v1)
    dim1 = v1.shape[-1]     # dimension of vector
    n1   = len(v1.shape)   # How many vectors
    
    v2 = np.asarray(v2)
    dim2 = v2.shape[-1]     # dimension of vector
    n2   = len(v2.shape)   # How many vectors
    vector = np.zeros(v1.shape,lb3dFLOAT)
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



def ADcrossproduct(v1, v2):

    """
    Dr. Birk's Cross Product Function
       take cross product of 3D vectors
    """
    v1 = np.asarray(v1)
    dim1 = v1.shape[-1]     # dimension of vector
    n1   = len(v1.shape)   # How many vectors
    
    v2 = np.asarray(v2)
    dim2 = v2.shape[-1]     # dimension of vector
    n2   = len(v2.shape)   # How many vectors
    vector = []
    if dim1 != dim2:
        print " *** Error **** crossproduct: Vector dimensions do not match!"
        return -1
    if n1 != n2:
        print " *** Error **** crossproduct: confused array length!"
        return -2
    if dim1 != 3 or dim2 != 3:
         print " *** Error **** crossproduct: only 3D vectors!"
       
    vector.append(v1[1]*v2[2] - v1[2]*v2[1])
    vector.append( v1[2]*v2[0] - v1[0]*v2[2])
    vector.append( v1[0]*v2[1] - v1[1]*v2[0])

    return vector

def linear_pts(start, end, number):
    """
        start   : vector in R2 or R3
        end     : vector in R2 or R3
        number  : total number of points in the linear sequence
    """
    change_vector = end-start
    #dist = np.linalg.norm(change_vector)
    #nump = number-2
    pointlist = [start]
    for i in range(1,number-1,1):
        newpoint = start+i*change_vector/(number-1)
        pointlist.append(newpoint)
    pointlist.append(end)
    curve_array = np.asarray(pointlist)
    return curve_array
    
    
def set_curve_knots(curveA, curveB):
    """
        Routine to change curve
        to a specified knot vector
        
        INCOMPLETE! - not developed, untested, and generally wrong.
    """
    curve1 = copy.deepcopy(curveA)
    curve2 = copy.deepcopy(curveB)
    #check curve order:
    if curve1.k == curve2.k:
        #check curve knot vector:
        if curve1.t == curve2.t:
            pass
        else:
            for knot in curve1.t:
                if knot in curve2.t:
                    pass
                else:
                    curve2.t , curve2.vertices  = knot_insertion(curve2.k,curve2.vertices,curve2.t,knot)
            for knot in curve2.t:
                if knot in curve1.t:
                    pass
                else:
                    curve1.t , curve1.vertices  = knot_insertion(curve1.k,curve1.vertices,curve1.t,knot)

        curve1.n = len(curve1.vertices)
        curve2.n = len(curve2.vertices)
    else:
        print 'WARNING: Curves have different orders'
        print 'implement a curve order change algorithm to proceed'
        ## then move knot insertion out of the other branch - to do it for this casetoo.
        ## since this case is not yet implemented, just exit with failure.
    
    
    return curve1, curve2
        