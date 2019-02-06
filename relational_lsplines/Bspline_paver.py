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
from    initialValues import InitializeControlPoints, InitializeControlVertices, interval_bounds, lagrangian_bounds

from  automatic_differentiation import ad
from  automatic_differentiation import IntervalAnalysis as IA
from  interval_arithmetic import ia
#from minDistance import plot_all, naive_minimizer
import copy
from   ADILS             import IntervalLagrangeSpline, Lagrangian
from   interval_analysis import IntervalAnalysis as IA
from   FormParameter     import FormParameterDict


def check_good_pts(X):
    xpts = X[0]
    ypts = X[1]
    bad = True
    n = len(xpts)
    for i in range(n):
        if (1. in xpts[i].grad[0,i]):
            bad = False
        else:
            print 'bad X Vertices!!'
            bad = True
            break
        if (1. in ypts[i].grad[0,i+n]):
            bad = False
        else:
            print 'bad Y Vertices!!'
            bad = True

    return bad
    

class CurvePaver(object):
    def __init__(self, X, tol=.5):
        self.outside                = []
        self.inside                 = []
        self.boundary               = []
        self.bad                    = []
        self.L                      = [X]
        self._ntboxes               = len(self.L)
        self.dim                    = X.curve.dim
        self.n                      = X.curve.n
        self.tol                    = tol
        self.large                  = 1.e10
        self.small                  = 1.e-6
        self.split_which            = 0
        self.maxit                  = 25
        self.upper_bound            = self.large
        self.max_curve_degeneracy   = 4
        assert(self.small < self.tol)
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

    def exclude_by_upper_bound(self, f):
        """ If True discard
        """
        ans = False
        if self.upper_bound < f.value.inf: ans = True
        return ans
        
    def check_gradient_consistency(self, Lspline):
        """ If False discard
        """
        print 'Now check Gradient for overall consistency:'
        tarray = Lspline.ia.check_for_stationary_pt(Lspline.f.ider[Lspline.mask])
        if False in tarray:
            print 'infeasible gradient-> discard!'
            return False
        else:
            return True
        
    def check_constraint_consistency(self, Lspline):
        """ If False discard
        """
        print 'Checking Individual Constraint Consistency'
        all_individually_consistent = True
        for i in range(len(Lspline.Lagrangian.equality)):
            if not Lspline.monotonicity_constraint(i).contains(0.):
                all_individually_consistent = False
                break
        if not all_individually_consistent:
            print 'Lspline failed individual consistency checks -> discard!'
            return False
        else:
            return True
            
    def degenerate_vertex_width(self):
        """ If True discard
        """
        ok_count = 0    
        for vertex in Lspline.curve.interval_xpts[1:-1]:
            if vertex.value.width() < self.small:
                ok_count +=1
        for vertex in Lspline.curve.interval_ypts[1:-1]:
            if vertex.value.width() < self.small:
                ok_count +=1
        if ok_count > self.max_curve_degeneracy:
            print 'Something IS Wrong:'
            print 'the control points have collapsed'
            return True
        else:
            return False
            
    def list_exclude_by_constraint(self, i):
        for i in range(len(self.L)):
            if not self.L[i].get_constraint(i).contains(0.):
                self.L.pop(i)
        return
        
    def do_splitting(self, Lspline):
        self.pick_direction_to_split([Lspline.curve.interval_xpts,
                                      Lspline.curve.interval_ypts])
        Xmin,Xmax = self.split_box(Lspline)
        if self.check_constraint_consistency(Xmin):
            self.L.append(Xmin)
        if self.check_constraint_consistency(Xmax):
            self.L.append(Xmax)
        return
            
    def split_box(self, Lspline):
        X = [Lspline.curve.interval_xpts,
             Lspline.curve.interval_ypts]
        self.saveL = copy.deepcopy(Lspline)
        l1 = copy.deepcopy(X)
        l2 = copy.deepcopy(X)
        xmid = X[self.split_which][self.split_which_index].getpoint(.45)
        a = copy.deepcopy(X[self.split_which][self.split_which_index])
        b = copy.deepcopy(X[self.split_which][self.split_which_index])
        a.value.sup = xmid.value
        b.value.inf = xmid.value
        l1[self.split_which][self.split_which_index] = a
        l2[self.split_which][self.split_which_index] = b
        dummy_a = copy.deepcopy(Lspline)
        dummy_a.curve.interval_xpts = l1[0]
        dummy_a.curve.interval_ypts = l1[1]
        
        dummy_b = copy.deepcopy(Lspline)
        dummy_b.curve.interval_xpts = l2[0]
        dummy_b.curve.interval_ypts = l2[1]
        
        print 'dim = ',self.split_which
        print 'index = ',self.split_which_index
        return dummy_a,dummy_b
        
    def pick_direction_to_split(self, X):
        xpts = X[0]
        ypts = X[1]
        max_xlen = 0.
        max_ylen = 0.
        xi = 0
        yi = 0
        for i in range(1,self.n-1):
            test_x = max(max_xlen, xpts[i].width())
            test_y = max(max_ylen, ypts[i].width())
            if test_x > max_xlen:
                xi = i
                max_xlen = test_x
            if test_y > max_ylen:
                yi = i
                max_ylen = test_y
        if xpts[xi].width() < ypts[yi].width():
            self.split_which = 1
            self.split_which_index = yi
        else:
            self.split_which = 0
            self.split_which_index = xi
        return 
        
    
    def compute_subpaving(self, stack = -1):
        inum = 0
        while len(self.L) > 0 and inum < self.maxit:
            self.pave_one(stack)
            inum += 1
        return
        
    def compute_simple_paver(self, stack = -1):
        inum = 0
        while len(self.L) > 0 and inum < self.maxit:
            self.simple_pave_one(stack)
            inum+=1
        return
    
    def pave_one(self, which = -1):
        print '\nGetting new LSpline from the list'
        Lspline = self.L.pop(which)
        
        print 'check control point properties'
        ck_X = Lspline.get_ivertices()
        bad = check_good_pts(ck_X)
        if bad:
            print 'Something is Wrong:'
            print 'the control points have bad gradient informaton!!'
            print '->discard to bad list'
            self.bad.append(Lspline)
            return
        
        Lspline.compute_box_consistency(box_type = 'm')
        Lspline.contract_constraints()
        
        try:
            Lspline.interval_newton(itr=1) 
        except:
            if self.degenerate_vertex_width():
                print '->discard to inside list'
                self.inside.append(Lspline)
                return
            else:
                self.do_splitting(Lspline)
                return 

        check_no_sol = Lspline.prove_no_solution()
        if check_no_sol[0]:
            print 'no solution in this curve space -> discard'
            self.outside.append(Lspline)
            return 
            
        gradient_isnan = np.isnan(Lspline.f.ider).any()
        if gradient_isnan:
            print 'gradient is not finite.  Split'
            self.do_splitting(Lspline)
            return 
            
        gradient_consistent = self.check_gradient_consistency(Lspline)
        if not gradient_consistent:
            print 'inconsistent gradient -> discard'
            self.outside.append(Lspline)
            return
            
        constraint_consistent = self.check_constraint_consistency(Lspline)
        if not constraint_consistent:
            print 'at least on inconsistent constraint-> discard'
            self.outside.append(Lspline)
            return
        
        ## Hansen certain feasibility  (p 180) and ch 5.1
        ## we can use any certainly feasible pt
        ## to compute a least upperbound on f
        ## could try thin_f, since we already have it... 
            ##if we think its worth it to check it for consistency
        self.upper_bound = min(self.upper_bound, 
                               Lspline.f.value.sup)
        if self.exclude_by_upper_bound(Lspline.f):
            print 'curve is non optimal-> discard'
            self.outside.append(Lspline)
            return
  
        check_zero = Lspline.prove_simple_zero() 
        if check_zero[0]:
            print 'zero proven in this curve space'
            Lspline.interval_newton(itr=1)
            stuck  = Lspline.prove_stuck() 
            if stuck[0]:
                for el in Lspline.Lagrangian.equality:
                    Lspline.compute_box_consistency(which=[el], box_type = 'm')
                print 'stuck -> split'
                self.do_splitting(Lspline)
                return 
                
            elif not stuck[0]:
                print 'not stuck'
                Lspline.interval_newton(itr=5)
                #stuck  = Lspline.prove_stuck() 
                if Lspline.ia.check_for_convergence(Lspline.f.ider[Lspline.mask], self.tol):
                    print 'add solution to ->boundary'
                    self.boundary.append(Lspline)
                    return
                else:
                    print 'no convergence->return to top of the list'
                    self.L.append(Lspline)
                    return
        
        print 'end of checker -> split'
        self.do_splitting(Lspline)
        return
        
    def simple_pave_one(self, which = -1):
        Lspline = self.L.pop(which)
        constraint_consistent   = self.check_constraint_consistency(Lspline)
        if not constraint_consistent:
            self.outside.append(Lspline)
            return
        self.do_splitting(Lspline)
        return
        
    def plot(self, xi=0, yi=0):
        self.plot_subpaving(xi=0,yi=0)
        return

    def plot_subpaving(self, xi=0, yi=0):
        fig = plt.figure(1)
        ax = SubplotZero(fig, 111)
        fig.add_subplot(ax)
        for el in self.outside:
            self.plot_box(X = el, x_index = xi, y_index=yi,
                          this_color = 'red', canvas = ax)
        for el in self.boundary:
            self.plot_box(X = el, x_index = xi, y_index = yi,
                          this_color = 'yellow', canvas = ax)
        for el in self.inside:
            self.plot_box(X = el, x_index = xi, y_index = yi,
                              this_color = 'green', canvas = ax)
        ax.set_xlim(-12., 12.)
        ax.set_ylim(-12., 12.)
        ax.axis('equal')
        plt.show()
        return
    
    def plot_box(self, X, 
                 x_index = 0, y_index =0,
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
        
        ix = ixpts[x_index]
        iy = iypts[y_index]
        xmin = ix.min
        xmax = ix.max
        ymin = iy.min
        ymax = iy.max
        
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
            #ax.set_ylim(bottom.min(), top.max())
            ax.axis('equal')
            plt.show()
            return 
        else:
            return ax
            
    def print_sol(self):
        for el in self.boundary:
            print el.f
        return
    
    def print_L(self):
        for el in self.L:
            print ''
            print el
        return
    
    def print_f(self):
        for i,el in enumerate(self.L):
            print el.f
            print ''
        return
        
    def print_gradf(self):
        for i,el in enumerate(self.L):
            print el.f.ider
            print ''
        return
            
    def print_vertices(self, ini=0,fini=None):
        
        for i,el in enumerate(self.L[ini:fini]):
            print 'self.L[{}]'.format(i)
            print 'Gradient:'
            print el.f.ider
            print 'Vertices:'
            print el.curve.interval_xpts
            print ''
            print el.curve.interval_ypts
            print '---------------------'
        return

#
# Hansen Test problems, G.O. version 1, page 141
#
def thcf(X):
    X1 = X[0][0]
    X2 = X[1][0]
    return 12.*(X1**2) - 6.3*(X1**4) + (X1**6) + 6.*(X2**2 - X1*X2)
    #return 12.*(X1**2) - 6.3*(X1**4) + (X1**6) + 6.*X2*(X2 - X1)
    
def test_problem1():
    X1 = ad(ia(-12.,10.), name = 'X1', N=2, dim = 0)
    X2 = ad(ia(-11.,10.3), name = 'X2', N=2, dim = 1)
    C = ad( ia(-20.,20.), name = 'h0', N=2, dim = -1)
    X = [[X1],[X2]]
    """
    func = thcf
    tol = 1.e-10
    """
    sp = CurvePaver(thcf,X,C,tol = 0.5)
    #Xn, nit = IntervalAnalysis.compute_interval_newton_basic(thcf,X,itmax=100)
    sp.compute_subpaving()
    return
    
if __name__ == '__main__':
    option  = 'ini_test'
    test    = 'basic'
    
    if option =='ini_test':
        print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 12.
        xe = 12.
        ye = 0.
        alphab = 0.#-5.
        alphae = 0.#-15.
        Cab_given = 0.
        Cae_given = 0.
        xc = 4.
        yc = 4.
        curve_area = 72.#84.#
        slope = 'down'
        
        ae = alphae
        ab = alphab
        
        ini_v = InitializeControlVertices(alphae=ae, alphab=ab,
                                          Cab_given=0.,Cae_given=0.,
                                                              nCV = 7)
        curve = spline.Bspline(ini_v.vertices, k, nump)  
        wi = 3.5 #2.5#
        w  = .5 #.65
        ep = 1.e-10
        sp = 1.e-4
        x = ini_v.vertices[:,0]
        y = ini_v.vertices[:,1]
        interval_data={}
        #interval_data, small = interval_bounds(curve)
        # Max March 2015 ultimate test -Mod used August 2015
        """
        interval_data['xmin'] = [ 0.0,   0.01, 0.1, 0.2, 0.000000001, 0.000000001, 11.999999999 ] 
        interval_data['xmax'] = [ 0.0000000001,   11.9, 11.9, 11.9, 11.99, 11.99, 12.0000000001 ]
        interval_data['ymin'] = [11.9999999999,   0.01, 0.01, 0.0, 0.,0.,0.]# -0.000000001, -0.000000001, -0.000000001  ] #
        interval_data['ymax'] = [12.0000000001,  12.0, 12., 12., 12., 11., 0.00000000001]
        #"""
        """
        interval_data['xmin'] = [ 0.0,   0.01, 0.1, 0.2, 0.3, 0.4, 11.999999999 ] 
        interval_data['xmax'] = [ 0.0000000001,   11.6, 11.7, 11.8, 11.9, 11.99, 12.0000000001 ]
        interval_data['ymin'] = [11.9999999999,   0.01, 0.01, 0.01, 0.0, 0.0, -0.000000001]# -0.000000001, -0.000000001, -0.000000001  ] #
        interval_data['ymax'] = [12.0000000001,  12.0000000001, 12.0000000001, 12., 12., 11., 0.00000000001]
        #"""
        #wide and even March 2015 sucess
        """
        interval_data['xmin'] = [x[0],    x[1]-w-w,  x[2]-wi, x[3]-wi, x[4]-wi, x[5]-wi, x[6]-ep] 
        interval_data['xmax'] = [x[0]+ep, x[1]+wi, x[2]+wi, x[3]+wi, x[4]+wi, x[5]+w, x[6]]
        interval_data['ymin'] = [y[0]-ep, y[1]-ep, y[2]-wi, y[3]-wi, y[4],    y[5],    y[6]] 
        interval_data['ymax'] = [y[0],    y[1],    y[2],    y[3]+wi, y[4]+wi, y[5]+ep, y[6]+ep]
        #"""
        #consistent and even August 17 2015
        #"""
        interval_data['xmin'] = [x[0]-ep, x[1]-w,  x[2]-wi, x[3]-wi, x[4]-wi, x[5]-wi, x[6]-ep] 
        interval_data['xmax'] = [x[0]+ep, x[1]+wi, x[2]+wi, x[3]+wi, x[4]+wi, x[5]+w, x[6]]
        interval_data['ymin'] = [y[0]-ep, y[1]-ep, y[2]-wi, y[3]-wi, y[6]-ep,    y[6]-ep,    y[6]-ep] 
        interval_data['ymax'] = [y[0]+ep, y[0]+ep,    y[0]+ep,    y[3]+wi, y[4]+wi, y[5]+ep, y[6]+ep]
        #"""
        #interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        #"""
        interval_data['lmin'] = [-500.,-100000.,-100000.,-200000.,-200000.,-200000.]
        interval_data['lmax'] = [500.,100000.,100000.,200000.,200000.,200000.]
        #"""
        
        
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = curve_area)#value = AF.fuzzyNumber(60.,70.,80.) )#
        FPD.add_AngleConstraint(kind='equality', location = 0., value = alphab)#AF.fuzzyNumber(-5.,-2.,0.))#
        FPD.add_AngleConstraint(kind='equality', location = 1., value = alphae)
        FPD.add_CurvatureConstraint(kind='equality', location = 0., value = Cab_given)
        FPD.add_CurvatureConstraint(kind='equality', location = 1., value = Cae_given)
        #"""
        FPD.add_E1(kind='LS', weight = 1.0)
        FPD.add_E2(kind='LS', weight = 0.5)
        FPD.add_E3(kind='LS', weight = 0.5)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        
        L = Lagrangian(FPD)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        
#        Lspline.check_IA_correctness = True
#        self = Lspline
#        self.contract_constraints()
        
        Lspline.contract_constraints()
        sp1 = CurvePaver(Lspline)
        self = sp1
        print sp1.ntboxes
        #self.compute_subpaving()