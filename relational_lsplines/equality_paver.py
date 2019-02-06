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

from  automatic_differentiation import ad, VectorAD
from  automatic_differentiation import IntervalAnalysis as IA
from  interval_arithmetic import ia
#from minDistance import plot_all, naive_minimizer
import copy
from   ADILS             import IntervalLagrangeSpline, Lagrangian
from   FormParameter     import FormParameterDict



class UnconstrainedPaver(object):
    """See Jaulin et. al. for subpaver development
    This sub paver is for use with the ad(ia) class combination
    Presently ad isn't really used, but this will pave the way
    so to speak
    """
    def __init__(self, solution_function, X, tol=1.e-5):
        self.outside        = []
        self.inside         = []
        self.boundary       = []
        self.L              = [X]
        self._ntboxes       = len(self.L)
        self.n              = len(X[0])
        #self.N              = len(X)
        self.tol            = tol
        if tol<1.e-10: 
            print 'minimum tolerance overide'
            print 'minimum is 1.e-10'
            self.tol = 1.e-10
        #self.C              = constraint
        self.f              = solution_function 
        self.split_which    = 0
        self.maxit          = 300
        self.global_min     = 10000000.05
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
        
    def get_mid(self, X):
        mid = []
        for i in range(2):
            mid.append([])
            for j in range(self.n):
                mid[i].append(X[i][j].value.getpoint(.5))
        return mid
    
    def compute_subpaving(self):
        inum = 0
        while len(self.L) >0 and inum < self.maxit:
            self.pave_one()
            inum += 1
        return
    
    def pave_one(self):
        """Hansen page 133
        """
        tol_ = self.tol
        X = self.L.pop(-1)
        y = self.eval_system(X)#, C)
        y.n = len(y.hess)
        y = VectorAD.generate_numpy_repr(y)
        
        #1.) update the global min
        midx = self.get_mid(X)
        midf = self.f(midx)
        self.global_min = min(self.global_min, midf)
        
        
        gradient_consistent = self.check_gradient_consistency(y)
        if not gradient_consistent:
            print 'inconsistent gradient -> discard'
            self.outside.append(X)
            return
            
        hessian_consistent = self.check_gradient_consistency(y)
        if not hessian_consistent:
            print 'inconsistent hessian -> discard'
            self.outside.append(X)
            return
            
        """
        func = self.f
        Bspline_style = True
        #"""
        # try to reduce the box:
        converged = False
        X, nit, y, converged = IA.compute_interval_newton_basic(self.f,
                                                                     X,
                                                                     itmax=5,
                                                                     tol=tol_,
                                                                     Bspline_style = True)
        y = VectorAD.generate_numpy_repr(y)

        consistent_update = self.check_solution_consistency(X)
        small_enough = self.check_width(X)
        
        if consistent_update and ( not converged or not small_enough):
            print 'split the box'
            self.pick_direction_to_split(X)
            Xmin,Xmax = self.split_box(X)
            self.L.append(Xmin)
            self.L.append(Xmax)
            return
        #if  converged or small_enough:
        elif consistent_update and (converged or small_enough):            
            print 'accept solution'
            self.boundary.append(X)
            return
        elif not consistent_update:
            print 'inconsistent by Newtons method => Discard'
            self.outside.append(X)
            return
        else:
            print '-------------ERROR-------------'
            #self.outside.append(X)
            self.pick_direction_to_split(X)
            Xmin,Xmax = self.split_box(X)
            self.L.append(Xmin)
            self.L.append(Xmax)
            return
        
    def check_width(self, X):
        xw = []
        for i in range(self.n):
            xw.append(X[0][i].value.width() <= self.tol)
        for i in range(self.n):
            xw.append(X[0][i].value.width() <= self.tol)
        if False in xw:
            return False
        else:
            return True
        
    def check_solution_consistency(self, X):
        tarray = []
        for i in range(self.n):
            tarray.append(X[0][i].value.isempty)
        for i in range(self.n):
            tarray.append(X[1][i].value.isempty)
        if True in tarray:
            return False
        else:
            return True
        
    def check_gradient_consistency(self, f):
        tarray = np.asarray([f.igrad[:,0] <= 0.,0. <= f.igrad[:,1]])
        if False in tarray:
            return False
        else:
            return True
    
    def check_hessian_consistency(self, f):
        tarray = []
        for i in range(self.n*2):
            #tarray.append( f.ihess[i,i,0] >=  0.)
            tarray.append( [f.ihess[i,i,0] <=  0. , 0. <= f.ihess[i,i,0] ])
        tarray = np.asarray(tarray)
        if False in tarray:
            return False
        else:
            return True
        
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
        
    def plot_L(self):
        fig = plt.figure(1)
        ax = SubplotZero(fig, 111)
        fig.add_subplot(ax)
        for el in self.L:
            self.plot_box(X = el, index = 0,
                          this_color = 'blue', canvas = ax)
        ax.set_xlim(-15., 15.)
        ax.set_ylim(-15., 15.)
        ax.axis('equal')
        plt.show()
    
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

#
# Testing a function:
#
def func0(X): 
    x = X[0]
    y = X[1]
    return (x[0]**2 + y[0]**2)

def my_test1():
    x = ad( ia(-10.,10.), N=2, dim = 0)
    y = ad( ia(-10.,10.), N=2, dim = 1)
    X = [[x],[y]]
    sp1 = UnconstrainedPaver(func0, X)
    self = sp1
    #sp1.compute_subpaving()
    #sp1.plot_subpaving()
    #print self.ntboxes
    return X,sp1
#
# Hansen Test problems, G.O. version 1, page 141
#
def thcf(X):
    X1 = X[0][0]
    X2 = X[1][0]
    #return (X1*X1)*12. - (X1*X1*X1*X1)*6.3 + (X1*X1*X1*X1*X1*X1) + X2*(X2 - X1)*6.
    #return (X1**2)*12. - (X1*X1*X1**2)*6.3 + (X1*X1*X1*X1*X1**2) + X2*(X2 - X1)*6.
    #return (X1**2)*12. - (X1*X1*X1**2)*6.3 + (X1*X1*X1*X1*X1**2) + X2*(X2 - X1)*6.
    #return (X1**2)*12. - ((X1**2)**2)*6.3 + ((X1**2)**3)*1.0 + (X2**2 - X1*X2)*6.#X2*(X2 - X1)*6.#
    
    
    #return 12.*(X1**2) - 6.3*((X1**4)) + ((X1**2)**3) + 6.*(X2**2 - X1*X2)
    
    #return 12.*(X1**2) - ((X1**2)**2)*6.3 + ((X1**2)**3) + X2*(X2 - X1)*6.
    #return 12.*(X1**2) - ((X1**2)**2)*6.3 + ((X1**2)**3) + X2*(X2 - X1)*6.
    #return 12.*(X1*X1) - ((X1**2)**2)*6.3 + ((X1**2)**3) + X2*(X2 - X1)*6.
    #return (X1*X1)*12. - ((X1**2)**2)*6.3 + ((X1**2)**3) + X2*(X2 - X1)*6.
    #return (X1**2)*12. - (X1**4)*6.3 + (X1**6) + (X2**2 - X1*X2)*6.
    #return (X1**2)*12. - (X1**4)*6.3 + (X1**6) + X2*(X2 - X1)*6.
    #return 12.*(X1**2) - 6.3*(X1**4) + (X1**2)**3 + 6.*(X2**2 - X1*X2)
    #return 12.*(X1**2) - 6.3*((X1**4)) + (X1**2)**3 + 6.*X2*(X2 - X1)
    return 12.*(X1**2) - 6.3*(X1**4) + (X1**6) + 6.*X2*(X2 - X1)
    
    ## toy problem:
    #return X1**2 + X2**2 +12.*(X1**2) + (X1**2) - 6.*X2*(X2 - X1) #+ 4.*X2
    
def test_problem1():
    X1 = ad(ia(-1100.,1400.), name = 'X1', N=2, dim = 0)
    X2 = ad(ia(-1400.,1200.), name = 'X2', N=2, dim = 1)
    #X1 = ad(ia(-2.,4.), name = 'X1', N=2, dim = 0)
    #X2 = ad(ia(-4.,2.), name = 'X2', N=2, dim = 1)
    
    #X1 = ad(ia(-11.,14.), name = 'X1', N=2, dim = 0)
    #X2 = ad(ia(-14.,12.), name = 'X2', N=2, dim = 1)
    #X1 = ad(ia(-.11,.14), name = 'X1', N=2, dim = 0)
    #X2 = ad(ia(-.14,.12), name = 'X2', N=2, dim = 1)
    #C = ad( ia(-20.,20.), name = 'h0', N=2, dim = -1)
    X = [[X1],[X2]]
    """
    func = thcf
    tol = 1.e-10
    """
    sp = UnconstrainedPaver(thcf,X,tol = 1.e-10)
    #
    #sp.compute_subpaving()
    """ 
    X1 = ad(ia(-4.1,4.), name = 'X1', N=2, dim = 0)
    X2 = ad(ia(-4.,4.), name = 'X2', N=2, dim = 1)
    X = [[X1],[X2]]
    func = thcf
    tol = 1.e-20
    Bspline_style = True
    Xn, nit, y, converged  = IA.compute_interval_newton_basic(thcf,
                                                                      X,itmax=100,
                                                                      tol = 1.e-25)
    #"""
    return X, sp
    
if __name__ == '__main__':
    X,sp = test_problem1()
    x1 = X[0][0]
    x2 = X[1][0]
    self = sp
    """
    self.compute_subpaving()
    for i, el in enumerate(self.boundary):
        F = self.f(el)
        F.name = 'Minimum {}'.format(i)
        F.print_components()
        
        print 'at min location:'
        print el
    
    X1 = ad(ia(-4.,2.), name = 'X1', N=2, dim = 0)
    X2 = ad(ia(-2.,4.), name = 'X2', N=2, dim = 1)
    #"""
    """
    #X1 = ad(ia(-11.,14.), name = 'X1', N=2, dim = 0)
    #X2 = ad(ia(-14.,12.), name = 'X2', N=2, dim = 1)
    #X1 = ad(ia(-.11,.14), name = 'X1', N=2, dim = 0)
    #X2 = ad(ia(-.14,.12), name = 'X2', N=2, dim = 1)
    #C = ad( ia(-20.,20.), name = 'h0', N=2, dim = -1)
    func = thcf
    tol = 1.e-10
    Bspline_style = True
    V = [[X1],[X2]]
    #"""