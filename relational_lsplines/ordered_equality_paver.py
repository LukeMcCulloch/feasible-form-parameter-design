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


class box(object):
    def __init__(self, X, f=None):
        self.X = X
        self.f = f
        return
    
    def __repr__(self):
        return '{}: X={} f={}'.format(self.__class__.__name__,
                                  self.X,
                                  self.f)
    
    def __call__(self):
        return self.X
        
    def getKey(self):
        return self.f.value.inf
        
    def __cmp__(self, other):
        return self.getKey()
    #    if hasattr(other, 'f'):
    #        #return self.__cmp__(self.f.value.inf,other.f.value.inf)
    #        return self.f.value.__cmp__(other.f.value)

def getKey(item):
    return item.f.value.inf

class UnconstrainedPaver(object):
    """See Jaulin et. al. for subpaver development
    This sub paver is for use with the ad(ia) class combination
    """
    def __init__(self, solution_function, X, tol=1.e-5):
        self.outside        = []
        self.inside         = []
        self.boundary       = []
        self.L              = [box(X)]
        self._ntboxes       = len(self.L)
        self.n              = len(X[0])
        self.tol            = tol
        if tol<1.e-10: 
            print 'minimum tolerance overide'
            print 'minimum is 1.e-10'
            self.tol = 1.e-10
        #self.C              = constraint
        self.f              = solution_function 
        self.split_which    = 0
        self.maxit          = 10
        self.global_min     = 1.e10
        self.nc             = 0
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
        X = X()
        return self.f(X)
    
    def eval_system(self, X):#, C):
        return self.eval_f(X)# - C
        
    def split_box(self, X):
        X = X()
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
        l1 = box(l1)
        l2 = box(l2)
        #compute f on the box:
        l1.f = self.eval_system(l1)
        l2.f = self.eval_system(l2)
        return l1,l2
        
    def pick_direction_to_split(self, X):
        X = X()
        if X[0][0].value.width() < X[1][0].value.width():
            self.split_which = 1
        else:
            self.split_which = 0
        return 
    
    def get_mid(self, X):
        X = X()
        mid = []
        for i in range(2):
            mid.append([])
            for j in range(self.n):
                mid[i].append(X[i][j].value.getpoint(.5))
        return mid
    
    def compute_subpaving(self, which=-1,
                          plot=True,
                          maxit=None):
        """ which = 0  : breadth first
            which = -1 : depth first
        """
        print 'which = {}'.format(which)
        if maxit is None:
            maxit=self.maxit
        inum = 0
        #compute f for the initial box(es):
        #-ini f for consistency of opertion inside-
        for el in self.L:
            el.f = self.eval_system(el)
        #now pave it!:
        while len(self.L) >0 and inum < maxit:
            self.pave_one(which)
            self.L.sort(key=getKey, reverse=True)
            inum += 1
        if plot: self.plot()
        return
    
    def pave_one(self, which=-1):
        self.nc += 1
        X = self.L.pop(which)
        y = X.f     
        y.n = len(y.hess)
        y = VectorAD.generate_numpy_repr(y)
        
        #(3.)
        if y.value.inf > self.global_min:
            print 'not the global minimum'
            self.outside.append(X)
            return
        
        
        gradient_consistent = self.check_gradient_consistency(y)
        if not gradient_consistent:
            print 'inconsistent gradient -> discard'
            self.outside.append(X)
            return
            
        hessian_consistent = self.check_hessian_consistency(y)
        if not hessian_consistent:
            print 'inconsistent hessian -> discard'
            self.outside.append(X)
            return
        
        #upper bound on global min
        self.global_min = min(self.global_min, y.value.sup)
        
        """
        func = self.f
        Bspline_style = True
        #"""
        # try to reduce the box:
        converged = False
        gg = copy.deepcopy(X)
        # X = copy.deepcopy(gg)
        ilist, nit, y, converged, small_enough = IA.compute_interval_newton_basic_split2(self, 
                                                                                           func=self.f,
                                                                                            X=copy.deepcopy(X()),
                                                                                            itmax=1,
                                                                                            tol=self.tol,
                                                                                            Bspline_style = True)
        #X, nit, y, converged, small_enough, ilist=Xi, count, F, converged, small_enough, ilist
        assert(len(ilist)>0)
        for el in ilist:
            bel = box(el)
            bel.f = self.eval_system(bel)
            self.L.append(bel)
        
        X = self.L.pop(which)
        #X = box(X)
        y = VectorAD.generate_numpy_repr(y)
        
        
        self.global_min = min(self.global_min, y.value.sup)
        if y.value.inf > self.global_min:
            print 'not the global minimum'
            self.outside.append(X)
            return


        consistent_update = self.check_solution_consistency(X)
        consistent_update = consistent_update and self.check_gradient_consistency(y)
        consistent_update = consistent_update and self.check_hessian_consistency(y)
        
        small_enough = self.check_width(X)
        

            #if  converged or small_enough:
        #"""
        #if consistent_update and ( not converged or not small_enough):
        if consistent_update and ( not small_enough):
            print 'split the box'
            self.pick_direction_to_split(X)
            Xmin,Xmax = self.split_box(X)
            self.L.append(Xmin)
            self.L.append(Xmax)
            return
        #if  converged or small_enough:
        #elif consistent_update and (converged or small_enough): #needs to be converge not small for interval contraints
        elif consistent_update and ( small_enough):  
        #elif consistent_update and ( converged or small_enough):     
            print 'accept solution'
            print X
            self.boundary.append(X)
            return
        elif not consistent_update:
            print 'inconsistent by Newtons method => Discard'
            self.outside.append(X)
            return
        #"""
        else:
            print '-------------ERROR-------------'
            #self.outside.append(X)
            self.pick_direction_to_split(X)
            Xmin,Xmax = self.split_box(X)
            self.L.append(Xmin)
            self.L.append(Xmax)
            return
        
    def check_width(self, X):
        X = X()
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
        X=X()
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
            tarray.append( f.ihess[i,i,1] >=  0.)
            #tarray.append( [f.ihess[i,i,0] <=  0. , 0. <= f.ihess[i,i,1] ])
        tarray = np.asarray(tarray)
        if False in tarray:
            return False
        else:
            return True
    
    def check_boundary(self):
        self.gb = []
        self.global_min = ia.infinity
        self.global_solution = None
        #ming = self.global_min
        for X in self.boundary:
            X = X()
            y = self.eval_system(X)#, C)        
            y.n = len(y.hess)
            print y
            y = VectorAD.generate_numpy_repr(y)
            self.global_min = min(self.global_min, y.value.sup)
            if y.value.sup == self.global_min:
                self.global_solution = box(X)
            if y.value.inf > self.global_min:
                pass
            else:
                self.gb.append(box(X))
        return self.global_solution
    
    def try_min(self):
        temp = []
        for X in self.gb:
            X=X()
            converged = False
            X, nit, y, converged, small_enough = IA.compute_interval_newton_basic(self.f,
                                                                                     X,
                                                                                     itmax=10,
                                                                                     tol=self.tol,
                                                                                     Bspline_style = True)
            temp.append(box(X))
        self.gb = temp
        return
        
    def plot(self):
        self.plot_subpaving()
        return

    def plot_subpaving(self,plot_L = True):
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
        if plot_L:
            for el in self.L:
                self.plot_box(X = el, index = 0,
                              this_color = 'blue', canvas = ax)
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
        X=X()
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
            print '\n',el()
            el.f = self.f(el())
            print el.f# - self.C
        return
    
    def print_L(self):
        for el in self.L:
            print ''
            print el()[0]
            print el()[1]
        return
    
    def print_Lf(self):
        for el in self.L:
            print ''
            print el.f
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
    return 12.*(X1**2) - 6.3*(X1**4) + (X1**6) + 6.*(X2**2 - X1*X2)
    #return 12.*(X1**2) - 6.3*(X1**4) + (X1**6) + 6.*X2*(X2 - X1)
    
    ## toy problem:
    #return X1**2 + X2**2 +12.*(X1**2) + (X1**2) - 6.*X2*(X2 - X1) #+ 4.*X2

def schwefel1(X, n=50):
    x1 = X[0][0]
    x2 = X[1][0]
    f=0.
    for i in range(50):
        f += (x1+x2)**2
    return f
    
def schwefel2(X, n=50):
    """pb17 page 143
    """
    x1 = X[0][0]
    x2 = X[1][0]
    return (1.5 - x1 +x1*x2)**2 + \
           (2.25 - x1 + x1*x2**2)**2 + \
           (2.625 - x1 + x1*x2**3)**2

def schwefel3(X, n=50):
    """pb18 page 143
    """
    x1 = X[0][0]
    #x2 = X[1][0]
    f = 0.
    for i in range(len(X)):
        f += (x1 - X[i][0]**2)**2
        f += (X[i][0]-1.)**2
    return f
    
def schwefel32(X, n=50):
    """pb28 page 145
    """
    x1 = X[0][0]
    #x2 = X[1][0]
    f = 0.
    for i in range(len(X)):
        f += (x1 - X[i][0]**2)**2
        f += (1.-X[i][0])**2
    return f

def booth(X, n=2):
    """prob 20 page 143
    """
    x1 = X[0][0]
    x2 = X[1][0]
    return (x1 + 2.*x2 - 7.)**2 +\
            (2.*x1 + x2 -5.)**2
            
def rosenbrock(X):
    """prob 29 page 145
    """
    x1 = X[0][0]
    x2 = X[1][0]
    return 100.*(x2-x1**2)**2 + (x1-1.)**2
    
    
def test_problem2():
    X1 = ad(ia(-2.,4.), name = 'X1', N=2, dim = 0)
    X2 = ad(ia(-4.,2.), name = 'X2', N=2, dim = 1)
    X = [[X1],[X2]]
    sp = UnconstrainedPaver(thcf,X,tol = 1.e-6)    
    return
def test_problem1():
    #X1 = ad(ia(-1100.,1400.), name = 'X1', N=2, dim = 0)
    #X2 = ad(ia(-1400.,1200.), name = 'X2', N=2, dim = 1)
    #X1 = ad(ia(-2.,4.), name = 'X1', N=2, dim = 0)
    #X2 = ad(ia(-2.,4.), name = 'X2', N=2, dim = 1)
    X1 = ad(ia(-2.,4.), name = 'X1', N=2, dim = 0)
    X2 = ad(ia(-4.,2.), name = 'X2', N=2, dim = 1)
    #X1 = ad(ia(-4.,4.), name = 'X1', N=2, dim = 0)
    #X2 = ad(ia(-4.,4.), name = 'X2', N=2, dim = 1)
    #X1 = ad(ia(2.9,3.1), name = 'X1', N=2, dim = 0)
    #X2 = ad(ia(0.4,.6), name = 'X2', N=2, dim = 1)
    
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
    #1.0
    sp = UnconstrainedPaver(thcf,X,tol = 1.e-6)
    
    #16.
    sp = UnconstrainedPaver(schwefel1,X,tol = 1.e-6)
    
    #17.
    #sp = UnconstrainedPaver(schwefel2,X,tol = 1.e-6)

    #18.
    sp = UnconstrainedPaver(schwefel3,X,tol = 1.e-3)
    
    #20.
    sp = UnconstrainedPaver(booth,X,tol = 1.e-3)
    
    #28
    #sp = UnconstrainedPaver(schwefel32,X,tol = 1.e-3)
    
    #29
    sp = UnconstrainedPaver(rosenbrock,X,tol = 1.e-2)
    #
    #sp.compute_subpaving()
    """ 
    X1 = ad(ia(-4.1,4.), name = 'X1', N=2, dim = 0)
    X2 = ad(ia(-4.,4.), name = 'X2', N=2, dim = 1)
    X = [[X1],[X2]]
    func = thcf
    tol = 1.e-6
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
    self.tol=1.e-3
    which=-1
    #self.compute_subpaving(which = -1)
    self.compute_subpaving(which=-1,
                           plot=True,
                           maxit= 500)
    def do_one():
        self.compute_subpaving(which=-1,
                           plot=True,
                           maxit= 1)
        return
    def do_width_first():
        self.compute_subpaving(which=0,
                           plot=True,
                           maxit= 1)
        return
    #self.compute_subpaving(which = 0)
    #self.plot_subpaving()
    #for el in self.L:
    #    print el.f
        
    self.print_sol()
    """
    for i, el in enumerate(self.boundary):
        F = self.f(el())
        F.name = 'Minimum {}'.format(i)
        F.print_components()
        
        print 'at min location:'
        print el()
    #"""
    """
    self.maxit=1000
    self.compute_subpaving()
    for i, el in enumerate(self.boundary):
        F = self.f(el)
        F.name = 'Minimum {}'.format(i)
        F.value
        
        print 'at min location:'
        print el
        print '\n'
    """
    """
    self.compute_subpaving()
    for i, el in enumerate(self.boundary):
        F = self.f(el)
        F.name = 'Minimum {}'.format(i)
        F.value
        
        print 'at min location:'
        print el
        print '\n'
    """
    """
    X1 = ad(ia(-2.,4.), name = 'X1', N=2, dim = 0)
    X2 = ad(ia(-2.,4.), name = 'X2', N=2, dim = 1)
    
    #X1 = ad(ia(-11.,14.), name = 'X1', N=2, dim = 0)
    #X2 = ad(ia(-14.,12.), name = 'X2', N=2, dim = 1)
    #X1 = ad(ia(-.11,.14), name = 'X1', N=2, dim = 0)
    #X2 = ad(ia(-.14,.12), name = 'X2', N=2, dim = 1)
    #C = ad( ia(-20.,20.), name = 'h0', N=2, dim = -1)
    func = thcf
    #func = schwefel1
    #func = schwefel2
    tol = 1.e-10
    Bspline_style = True
    #V = [[X1],[X2]]
    X = [[X1],[X2]]
    if Bspline_style:
        x = X[0]
        y = X[1]
        m = len(x)
        x_update = []  
        for el in x:
            x_update.append(el)
        for el in y:
            x_update.append(el)
        Xi = np.asarray(x_update)
    else:
        Xi = copy.deepcopy(X)
        
    Xi = copy.deepcopy(Xi)
    n = len(X)
    saveXi = copy.deepcopy(Xi)
    xiold = [copy.deepcopy(Xi)]
    i=0
    a=IA.compute_interval_newton_basic(thcf, copy.deepcopy(X), itmax = 10, Bspline_style=True)
    
    #"""
    """
    
    #X1 = ad(ia(2.99,3.01), name = 'X1', N=2, dim = 0)
    #X2 = ad(ia(.49,.51), name = 'X2', N=2, dim = 1)
    
    X1 = ad(ia(-11.,14.), name = 'X1', N=2, dim = 0)
    X2 = ad(ia(-14.,12.), name = 'X2', N=2, dim = 1)
    #X1 = ad(ia(-.11,.14), name = 'X1', N=2, dim = 0)
    #X2 = ad(ia(-.14,.12), name = 'X2', N=2, dim = 1)
    #C = ad( ia(-20.,20.), name = 'h0', N=2, dim = -1)
    #func = thcf
    #func = schwefel1
    func = schwefel2
    tol = 1.e-10
    Bspline_style = True
    #V = [[X1],[X2]]
    X = [[X1],[X2]]
    if Bspline_style:
        x = X[0]
        y = X[1]
        m = len(x)
        x_update = []  
        for el in x:
            x_update.append(el)
        for el in y:
            x_update.append(el)
        Xi = np.asarray(x_update)
    else:
        Xi = copy.deepcopy(X)
        
    Xi = copy.deepcopy(Xi)
    n = len(X)
    saveXi = copy.deepcopy(Xi)
    xiold = [copy.deepcopy(Xi)]
    i=0
    
    b=IA.compute_interval_newton_basic(func, 
                                        copy.deepcopy(X), 
                                         itmax = 25, 
                                         tol = 1.e-5, 
                                         Bspline_style=True)
    #"""