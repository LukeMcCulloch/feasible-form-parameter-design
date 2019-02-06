# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 17:00:58 2015

@author: lukemcculloch
"""

import numpy             as np       
import matplotlib.pyplot  as plt
import matplotlib.patches as patches
import matplotlib.path    as path
from mpl_toolkits.axes_grid.axislines import SubplotZero

# TLM B-Spline Curve class
import  curve                 as spline
from    initialValues import InitializeControlPoints, InitializeControlVertices

from  automatic_differentiation import ad
from  automatic_differentiation import IntervalAnalysis as IA
from  interval_arithmetic import ia
#from minDistance import plot_all, naive_minimizer
import copy
from   ADILS             import IntervalLagrangeSpline, Lagrangian
from   interval_analysis import IntervalAnalysis as IA
from   FormParameter     import FormParameterDict



class SimplePaver(object):
    """See Jaulin et. al. for subpaver development
    This sub paver is for use with the ad(ia) class combination
    Presently ad isn't really used, but this will pave the way
    so to speak
    """
    def __init__(self, solution_function, X, constraint, tol=.5):
        self.outside        = []
        self.inside         = []
        self.boundary       = []
        self.L              = [X]
        self._ntboxes       = len(self.L)
        self.n              = len(X)
        self.tol            = tol
        self.C              = constraint
        self.f              = solution_function 
        self.split_which    = 0
        self.maxit          = 3000
        return
    
    @property
    def ntboxes(self):
        return self._ntboxes
    
    @ntboxes.getter
    def ntboxes(self):
        return len(self.L) + \
                len(self.inside) + \
                len(self.outside) + \
                len(self.boundary)
                
    def eval_f(self, X):
        return self.f(X)
    
    def eval_system(self, X):#, C):
        return self.eval_f(X)# - C
        
    def split_box(self, X):
        l1 = copy.deepcopy(X)
        l2 = copy.deepcopy(X)
        xmid = X[self.split_which][0].value.getpoint(.5)
        a = copy.deepcopy(X[self.split_which][0])
        b = copy.deepcopy(X[self.split_which][0])
        a.value.sup = xmid
        b.value.inf = xmid
        a.real = (a.value.inf+a.value.sup)*.5
        b.real = (b.value.inf+b.value.sup)*.5
        l1[self.split_which][0] = a
        l2[self.split_which][0] = b
        return l1,l2
        
    def pick_direction_to_split(self, X):
        if X[0][0].value.width() < X[1][0].value.width():
            self.split_which = 1
        else:
            self.split_which = 0
        return 
        
    
    def compute_subpaving(self):
        inum = 0
        while len(self.L) >0 and inum < self.maxit:
            self.pave_one()
            inum += 1
        return
    
    def pave_one(self):
        C = self.C
        X = self.L.pop(-1)
        y = self.eval_system(X)#, C)
        test = y.value & C.value
        if y.value in C.value:
            self.boundary.append(X)
        elif y < C:
          self.inside.append(X)
        elif test.isempty:
            self.outside.append(X)
        elif y.value.width() < self.tol:
            self.boundary.append(X)
        else:
            self.pick_direction_to_split(X)
            Xmin,Xmax = self.split_box(X)
            self.L.append(Xmin)
            self.L.append(Xmax)
        return
        
    def plot(self):
        self.plot_subpaving()
        return

    def plot_subpaving(self):
        fig = plt.figure(1)
        ax = SubplotZero(fig, 111)
        fig.add_subplot(ax)
        for el in self.outside:
            self.plot_box(X = el, index = 0,
                          this_color = 'red', canvas = ax)
        for el in self.boundary:
            self.plot_box(X = el, index = 0,
                          this_color = 'yellow', canvas = ax)
        for el in self.inside:
            self.plot_box(X = el, index = 0,
                              this_color = 'green', canvas = ax)
        ax.set_xlim(-12., 12.)
        ax.set_ylim(-12., 12.)
        ax.axis('equal')
        plt.show()
        return
    
    def plot_box(self, X, index = 0, 
                 this_color = 'black',
                 canvas = None ):
        
        if canvas == None:
            fig = plt.figure(1)
            ax = fig.add_subplot(111)
        else:
            ax = canvas
        nrects = 1
        nverts = nrects*(1+3+1)
        verts = np.zeros((nverts, 2))
        codes = np.ones(nverts, int) * path.Path.LINETO
        
        ixpts = X[0]
        iypts = X[1]
        
        left    = []
        right   = []
        top     = []
        bottom  = []
        
        ix = ixpts[index]
        iy = iypts[index]
        xmin = ix.value.inf
        xmax = ix.value.sup
        ymin = iy.value.inf
        ymax = iy.value.sup
        
        left.append(xmin)
        right.append(xmax)
        bottom.append(ymin)
        top.append(ymax)  
        
        left    = np.asarray(left)
        right   = np.asarray(right)
        top     = np.asarray(top)
        bottom  = np.asarray(bottom)
        
        codes = np.ones(nverts, int) * path.Path.LINETO
        codes[0::5] = path.Path.MOVETO
        codes[4::5] = path.Path.CLOSEPOLY
        
        verts[0::5,0] = left
        verts[0::5,1] = bottom
        verts[1::5,0] = left
        verts[1::5,1] = top
        verts[2::5,0] = right
        verts[2::5,1] = top
        verts[3::5,0] = right
        verts[3::5,1] = bottom
        
        barpath = path.Path(verts, codes)
        patch = patches.PathPatch(barpath, 
                                  facecolor=this_color, 
                                  edgecolor='black', 
                                  alpha=.4)
        ax.add_patch(patch)
        if canvas == None:
            ax.set_xlim(-12., 12.)
            #ax.set_ylim(bottom.value.inf(), top.value.sup())
            ax.axis('equal')
            plt.show()
            return 
        else:
            return ax
    
            
    def print_sol(self):
        for el in self.boundary:
            print self.f(el) - self.C
        return
    
    def print_L(self):
        for el in self.L:
            print ''
            print el[0]
            print el[1]
        return


if __name__ == '__main__':
    
    def func1(X):
        
        x = X[0]
        y = X[1]
        return (x[0]**2 + y[0]**2)# * x[0].exp() 
        #return (x[0]**2 + x[1]**2 + y[0]**2 + y[1]**2)
    
    x = ad( ia(-10.,10.), N=2, dim = 0)
    y = ad( ia(-10.,10.), N=2, dim = 1)
    X = [[x],[y]]
    
    x = ad( ia(0.1,0.3), name = 'x', N=2, dim = 0)
    y = ad( ia(0.1,0.3), name = 'x', N=2, dim = 1)
    X1 = [[x],[y]]
    #C = ad(ia(.9,1.1), N=2, dim = -1)
    #C = ad(ia(4.7,5.3), N=2, dim = -1)
    C = ad( ia(4.0,5.9), N=2, dim = -1)
    #C = ia(ia(0.,2.)
    sp1 = SimplePaver(func1, X, C)
    self = sp1
    
    sp1.compute_subpaving()
    sp1.plot_subpaving()
    
    print self.ntboxes