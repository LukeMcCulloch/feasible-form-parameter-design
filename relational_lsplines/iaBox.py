#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 19:48:12 2017

@author: luke
"""

import numpy as np
from itertools import islice, cycle
import copy
from extended_interval_arithmetic import ia
import extended_interval_arithmetic as iafuncs




## use the one in myspline.py
#class BoundsList(object):
#    """
#    myspline.py version Bounds used in THBspline Curves
#        DONE:  add method contains?
#        TODO: add disjunctive or
#        
#        TODO: switch to this one.
#    """
#    def __init__(self, bounds):
#        assert( np.asarray([isinstance(el, ia) for el in bounds]).any() )
#        self.bounds = [el for el in bounds]
#        self.isempty = np.asarray([el.isempty for el in bounds]).any()
#
#    def __getitem__(self, index):
#        return self.bounds[index]
#
#    def __repr__(self):
#        return 'BoundsList({})'.format(self.bounds)
#
#    def __str__(self):
#        return 'BoundsList({})'.format(self.bounds)
#
#    def __and__(self, other):
#        """Intersection of two boxes
#        """
#        new_dims = []
#        for ds, do in zip(self.bounds, other.bounds):
#            new_dims.append( ds & do )
#        return BoundsList(new_dims)
#
#    def __or__(self, other):
#        """Union of two boxes
#        """
#        new_dims = []
#        for ds, do in zip(self.bounds, other.bounds):
#            new_dims.append( ds | do )
#        return BoundsList(new_dims)
#
#    def add(self, interval):
#        self.bounds.append(interval)
#        return
#
#    def contains(self, other):
#        """is other in self?
#        """
#        ans = None
#        if isinstance(other, BoundsList):
#            check_dim = []
#            for ds, do in zip(self.bounds, other.bounds):
#                check_dim.append( ds.contains(do) )
#            iscontained = np.asarray([not el.isempty for el in check_dim]).any()
#            return iscontained
#        elif isinstance(other, float):
#            #print 'Hello!', self, other
#            check_dim = []
#            for ds in self.bounds:
#                check_dim.append( ds.contains(other) )
#            ans = [i for i, x in enumerate(check_dim) if x]
#            return ans #only indicies which contain the u location
#



class BoxList(object):
    """
        List of nD hyper_rectangles
    """
    def __init__(self, *boxes):
        self.boxes = boxes
    

    def __getitem__(self, index):
        return self.boxes[index]

    def __repr__(self):
        return 'BoxList({})'.format(self.boxes)

    def __str__(self):
        return 'BoxList({})'.format(self.boxes)
    
    def contains(self, *X):
        bc = self.get_box_containing(*X)
        if bc.isempty:
            return False
        else:
            return True
    
    def get_box_containing(self, *X):
        for box in self.boxes:
            if box.contains(*X):
                return box
        ndims = self.boxes[0].ndims
        empty_box = []
        for i in range(ndims):
            ebxi = ia(0.,0.)
            ebxi.isempty = True
            empty_box.append(ebxi)
        ans = Box(*empty_box)
        ans.isempty = True
        return ans
    
class Box(object):
    """
    nD hyper-rectangle
    
    bounds = [ia(xi,xs), ia(yi,ys), etc..]
    """
    def __init__(self, *X):
        self.ndims = len(X)
        self.indices = range(self.ndims)
        self._X = X
        self._bounds_u = X[0]
        self._bounds_v = X[1]
        self.isempty = False

    def __repr__(self):
        return 'Box({})'.format(self.X)

    def __str__(self):
        return 'Box({})'.format(self.X)
    
    def __getitem__(self, index):
        return self.X[index]
    
    def __call__(self):
        return self.X
    
    @property
    def bounds_u(self):
        return self._bounds_u
    
    @property
    def bounds_v(self):
        return self._bounds_v
    
    @property
    def X(self):
        return self._X
    @X.setter
    def X(self, X):
        for el in X:
            print el
        self.ndims = len(X)
        self.indices = range(self.ndims)
        self._X = X
        self._bounds_u = X[0]
        self._bounds_v = X[1]
        self.isempty = False
        return


    def __and__(self, other):
        """Intersection of two hyper-boxes
        """
        new_dims = []
        for ds, do in zip(self.X, other.X):
            new_dims.append( ds & do )
        return Box( *new_dims)

    def __or__(self, other):
        """Union of two hyper-boxes
        """
        new_dims = []
        for ds, do in zip(self.X, other.X):
            new_dims.append( ds | do )
        return Box(*new_dims)
    
    def contains(self, *X):
        """is the pt (X[0], X[1],...) in box?
        TODO:
        """
        check_dim = []
        for i,el in enumerate(X):
            check_dim.append( self[i].contains(el))
        iscontained = np.asarray([el for el in check_dim]).all()
        return iscontained#np.asarray([el for el in check_dim])



def gett(*x):
    return x
    
if __name__ == "__main__":
    b1 = Box(ia(0.,1.),ia(0.,1.))
    b2 = Box(ia(.3,.5),ia(.7,1.8))
    b3 = b1&b2
    b4 = b1|b2
    u=.8
    v=1.
    print b3.contains(u,v)
    #blist = BoxList(b3,b4)
    #b_new = blist.get_box_containing(.5,.5)