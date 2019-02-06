#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  5 12:40:31 2017

@author: luke
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from itertools import islice, cycle
import copy
from extended_interval_arithmetic import ia
from iaBox import BoxList, Box

from myspline import rbspline, BoundsList
import wavelet_matrices as wm
import curve as spline

np.set_printoptions(threshold= 200)

#def get_3D_instantiation_from_1part_shape
def instantiation3D_from_1part_shape(array,dim=3):
    array_shp = np.shape(array)
    return np.zeros((array_shp[0],array_shp[1],dim))



class level_properties(object):
    def __init__(self, nlevel,nu,nv):
        self.nlevel = nlevel
        self.nu = nu
        self.nv = nv
        
class levels(object):
    def __init__(self, *levels):
        self.levels = []
        for level in levels:
            self.levels.append(level)
            
class level_basis_data(object):
    """Class to hold 
    THBspline method compute_apb(self, u,v)
    data in nice and descriptive form
    
    This is form gathering active basis
    functions at all levels for a given (u,v)
    
    -which is still not yet specific enough
    """
    def __init__(self, surf, data):
        self.surf = surf
        self.level = surf.level
        self.uactive = data[0]
        self.vactive = data[1]
        self.uinactive = data[2]
        self.vinactive = data[3]

    def __repr__(self):
        return 'level{}basis_data()'.format(self.level)

    def __str__(self):
        return 'level{}basis_data()'.format(self.level)
    
    def print_(self):
        print 'uactive:'
        print self.uactive
        print 'vactive:'
        print self.vactive
        print 'uinactive:'
        print self.uinactive
        print 'vinactive:'
        print self.vinactive
        return

class Index(object):
    
    def __init__(self, surface):
        self.index = {}
        start_surf = surface.get_coarsest_level()
        self.index = self.index_hierarchical_surface(start_surf, 
                                                     level=0,
                                                     index = self.index)
        return
    
    def __call__(self, level):
        return self.index[level]
    
    def __getitem__(self, index):
        return self.index[index]

    def __repr__(self):
        return 'Index({})'.format(self.index)

    def __str__(self):
        return 'Index({})'.format(self.index)
    
    def index_hierarchical_surface(self, surface, level, index):
        if surface.child is None:
            index[level] = surface
            return index
        else:
            index[level] = surface
            return self.index_hierarchical_surface(surface.child, level+1)
    
    def addLevel(self, surface):
        top_surf = surface.get_finest_level()
        #bot_surf = surface.get_coarsest_level()
        #assert(surface.parent == bot_surf)
        #assert(top_surf == self.index[-1]) #dict keys, not list!!
        next = max(self.index.keys())+1
        """
        assert(top_surf == self.index[next-1])
        #"""
        self.index[next] = surface
        return next
    

def print_(mssg):
    print mssg

def noprint(mssg):
    return



