#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 17:12:08 2017

@author: luke  

created from adials_gui3d.py:
    
Created on Wed Apr 29 18:54:51 2015




@author: 

primarily Jake Vanderplas:
https://jakevdp.github.io/blog/2012/11/24/simple-3d-visualization-in-matplotlib/
    cube example for most of the 3D plotting, moving, zorder trick, and zooming.
    and quaternion class.



Jake Vanderplas:
https://jakevdp.github.io/blog/2012/11/24/simple-3d-visualization-in-matplotlib/
    cube example for most of the 3D plotting, moving, zorder trick, and zooming.
    and quaternion class.
    
Luke McCulloch:  added translation to this and generalized it from the cube 
     to B-splines.

Matplotlib docs and examples:
    for point picking and zooming

About Quaternions:

    Transformation to and from Axis-Angle Representation:
    
        Schneider and Eberly 
        Geometric Tools for Computer Graphics, (GTCG)
    
    For cool ideas that expand massively on Quaterions 
    and normed division algebras generlly, see
        Dorst, Frontijne, Mann
        Geometric Algebra for Computer Science, 112182
    
    Axis-Angle representation to and from Quaternion:
        GTCG p 897
        
        a = (xo,yo,zo)  : axis of rotation
        theta           : angle of rotation
        
        (1) Axis-Angle to Quaternion = w + xi + yj + zk 
            w = cos(theta/2)
            x = a[0]*sin(theta/2)
            y = a[1]*sin(theta/2)
            z = a[2]*sing(theta/2)
            
        (2) Quaternion to Axis-Angle:
            theta   = 2*arccos(w)
            a[0]    = x/sin(theta/2)
            a[1]    = y/sin(theta/2)
            a[2]    = z/sin(theta/2)
"""


import numpy as np

class Quaternion:
    """Quaternions for 3D rotations
    
        self.x = [w, x, y, z] : the Quaternion Quadruple
    """
    #static info:
    i = (0.,1.,0.,0.)
    j = (0.,0.,1.,0.)
    k = (0.,0.,0.,1.)
    def __init__(self, x):
        self.x = np.asarray(x, dtype=float)
        
    @classmethod
    def from_v_theta(cls, v, theta):
        """ Implement (1) above
        Construct quaternion (and normalize it) 
            v       : unit vector
            theta   : rotation angle 
        """
        theta = np.asarray(theta)
        v = np.asarray(v)
        
        s = np.sin(0.5 * theta)
        c = np.cos(0.5 * theta)
        vnrm = np.sqrt(np.dot(v,v))

        q = np.concatenate([[c], s * v / vnrm])
        return cls(q)

    def __repr__(self):
        return "Quaternion:\n" + self.x.__repr__()
        
    def __add__(self, other):
        return self.__class__( [(self[0]+other[0]),
                                (self[1]+other[1]),
                                (self[2]+other[2]),
                                (self[3]+other[3])])

    def __mul__(self, other):
        """ multiplication of two quaternions.
        see: Schneider, Eberly
        Geometric Tools for Computer Graphics p. 896"""
        if isinstance(other, Quaternion):
            prod = self.x[:, None] * other.x
    
            return self.__class__([(prod[0, 0] - prod[1, 1]
                                     - prod[2, 2] - prod[3, 3]),
                                    (prod[0, 1] + prod[1, 0]
                                     + prod[2, 3] - prod[3, 2]),
                                    (prod[0, 2] - prod[1, 3]
                                     + prod[2, 0] + prod[3, 1]),
                                    (prod[0, 3] + prod[1, 2]
                                     - prod[2, 1] + prod[3, 0])])
        if isinstance(other, float):
            Q = Quaternion([other, 0.,0.,0.])
            return self*Q
        if isinstance(other, np.ndarray):
            Q = Quaternion([0.,other[0],other[1],other[2]])
            return self*Q
            
                                 
    def conjugate(self):
        return self.__class__([ self.x[0],-self.x[1],-self.x[2],-self.x[3] ])

    def magnitude(self):
        return ((self*self.conjugate()).x).sum()
        
    def magnitude2(self):
        return np.sqrt(np.dot(self.x,self.x))
    
    def normalize(self):
        M = self.magnitude()
        Q = [el/M for el in self.x]
        return self.__class__(Q)
    
    def __div__(self, other):
        return NotImplemented
        
    def inverse(self):
        return self.conjugate()

    def as_v_theta(self):
        """Implement (2) above
        Starting wih a normalized Quaternion, Q.x = [w,x,y,z]
        Return the v, theta equivalent representation
        """
        # compute theta
        norm = np.sqrt(np.dot(self.x,self.x))
        theta = 2 * np.arccos(self.x[0] / norm)

        # compute the unit vector
        v = np.array(self.x[1:], order='F', copy=True)
        v /= np.sqrt(np.dot(v,v))

        return v, theta

    def as_rotation_matrix(self):
        """Return the rotation matrix of the (normalized) quaternion"""
        v, theta = self.as_v_theta()
        c = np.cos(theta)
        s = np.sin(theta)

        return np.array([[v[0] * v[0] * (1. - c) + c,
                          v[0] * v[1] * (1. - c) - v[2] * s,
                          v[0] * v[2] * (1. - c) + v[1] * s],
                         [v[1] * v[0] * (1. - c) + v[2] * s,
                          v[1] * v[1] * (1. - c) + c,
                          v[1] * v[2] * (1. - c) - v[0] * s],
                         [v[2] * v[0] * (1. - c) - v[1] * s,
                          v[2] * v[1] * (1. - c) + v[0] * s,
                          v[2] * v[2] * (1. - c) + c]])
                          
                          
