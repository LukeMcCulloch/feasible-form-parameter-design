#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 18:01:23 2018

@author: luke


Stuff that the THB sufaces and curves use to ... make themselves...
"""

import cPickle as pickle



level2n = { 0:4,
            1:5,
            2:7,
            3:11,
            4:19,
            5:35}
n2level = { 4:0,
            5:1,
            7:2,
            11:3,
            19:4,
            35:5}


class ProjectionGetter(object):
    """Placeholder to store the precomputed
    THB projection operators
    
    Parameters
    ----------
        
        rm_7x5  = rm_7x5 THB projection operator
        rm_11x7 = rm_11x7 THB projection operator
        rm_19x11 = rm_19x11 THB projection operator
        rm_35x19 = rm_35x19 THB projection operator
        rm_67x35 = rm_67x35 THB projection operator
        
    Returns
    ----------
        nothing
        
    Notes
    ----------
        This is designed to be fed to the
        
        function get_projection_matrices()
        when returned it will have all the matrices
    """
    def __init__(self):
        
        self.rm_7x5 = None
        self.rm_11x7 = None
        self.rm_19x11 = None
        self.rm_35x19 = None
        self.rm_67x35 = None
        
        self.map = {5:"rm_7x5",
                    7:"rm_11x7",
                    11:"rm_19x11",
                    19:"rm_35x19",
                    35:"rm_67x35"}
        
        self.inversemap = {"rm_7x5":5,
                           "rm_11x7":7,
                           "rm_19x11":11,
                           "rm_35x19":19,
                           "rm_67x35":35}
        
        self.exists = {5:False,
                       7:False,
                       11:False,
                       19:False,
                       35:False}
        
        self.check = {7:5,
                      11:7,
                      19:11,
                      35:19,
                      67:35}
        
        self.level2n = {0:4,
                        1:5,
                        2:7,
                        3:11,
                        4:19,
                        5:35}
        
        
        
    def __call__(self, level):
        n = self.level2n[level]
        exists = self.exists[n]
        name = self.map[n]
        if not exists:
            rm = self.get_projector(name)
            setattr(self, name, rm)
            self.exists[n] = True
            return rm
        else:
            return getattr(self, name)
             
    
    def get_projector(self, file_name):
        fileObject = open(file_name,'r')  
        rm = pickle.load(fileObject) 
        fileObject.close()
        return rm


