# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 16:04:29 2016

@author: lukemcculloch
"""
import numpy as np
#import scipy as sp
#import scipy.misc
from itertools import islice, cycle

import matplotlib.pyplot as plt
import copy

import curve as spline
linear_vertices = spline.linear_vertices
import KV as kv
import knot_operations as ko
import wavelet_matrices as wm
import HierarchicalBspline as hb
import thbasis as thb

from extended_interval_arithmetic import ia


np.set_printoptions(threshold= 25)


class Box(object):
    """
        TODO:  add method contains?
    """
    def __init__(self, dims):
        assert( np.asarray([isinstance(dim, ia) for dim in dims]).any() )
        self.dims = [dim for dim in dims]
        self.isempty = np.asarray([dim.isempty for dim in dims]).any()
    
    def __and__(self, other):
        """Intersection of two boxes
        """
        new_dims = []
        for ds, do in zip(self.dims, other.dims):
            new_dims.append( ds & do )
        return Box(new_dims)

    def __or__(self, other):
        """Union of two boxes
        """
        new_dims = []
        for ds, do in zip(self.dims, other.dims):
            new_dims.append( ds | do )
        return Box(new_dims)
    
    def __in__(self, other):
        """is other in self?
        """
        check_dim = []
        for ds, do in zip(self.dims, other.dims):
            check_dim.append( ds in do )
        iscontained = np.asarray([el.isempty for el in check_dim]).any()
        return iscontained
        
        
class KDnode(object):
    def __init__(self, box, level, position=None, 
                 axis=0, dim=2,
                 a=None, b=None):
        self.box = box      # Kiss: lower left and upper right corners
                            # here it will be to ia numbers
                            # x = ia(left,right)
                            # y = ia(bottom,top)
        self.level      = level
        self.axis       = axis    #[0,1,...dim-1][axis] 
                                  # = active axis of next splitting
        #self.axislist  = axislist  #axislist = [0,1,...dim-1]
        self.dim        = dim
        dimensions      = cycle(range(dim))
        self.axis_cycle = islice(dimensions, axis, None)
        self.position   = position
        #child boxes: (binary subdivision is extensable...easiest?)
        self.a = a #left KDnode
        self.b = b #right KDnode
    
    
    def cycle(self):
        """Cycle around the dim of the tree
            use itertools cycle for more elegance (DONE)
            https://stackoverflow.com/questions/10986970/
            python-how-to-toggle-between-two-values
        """
        val = next(self.axis_cycle)
        return val
    
    def compute_axis(self):
        return self.cycle()
    
    def compute_position(self):
        return 0.5
    
    def insertbox(self, B, Q, L):
        """Box, Qnode, Level
            if B is Q.box: 
                Q.level = L
            else: 
                'visit all nodes with root {}.format(Q)
                if the level of that node is less than L, 
                increase it to L'
        """
        if B is Q.box:
            Q.level = L
        else:
            intersect = B & Q.box 
            for child in [Q.a, Q.b]:
                if child is not None:
                    if not intersect.isempty:
                        self.insertbox(intersect, child, L)
                    else:
                        """KISS
                            -compute splitting value <position>
                            -compute axis <new_axis>
                        """
                        # compute axis
                        new_axis = self.compute_axis()
                        #compute splitting position
                        pos = self.compute_position()
                        """JUTTLER
                           -position = 0.5 for dyadic boxing
                        """
                        #1.) create the box childbox of child
                        childbox = child.box
                        intersect_child = B & childbox
                        if not intersect_child.isempty:
                            # create the node child(node), setting
                            # it's box to childbox,
                            # it's level to Q.level
                            # and children to null
                            childnode = KDnode(childbox, Q.L, axis=new_axis )
                            self.insertbox(intersect_child, childnode, L)
                            
                            """
                            After each box insertion we perform 
                                a cleaning step, 
                            -visiting all sub–trees;
                            -deleting those where all nodes have the same level. 
                            This 
                            -reduces the depth of the tree 
                                to a minimal value 
                            and 
                            -optimizes the performance of all algorithms.
                            """
                        """END-JUTTLER"""
                            
                        
                        """KISS"""
                        """
                        Q.position = pos
                        Q.axis = new_axis
                        #create nodes Q.a, Q.b
                        self.setchild(B & Q.a, Q, 0)
                        self.setchild(B & Q.b, Q, 1)
                        #"""
                        """END-KISS"""
                        break
        return
                        
    def setchild(self, C, Q, side):
        """CurrentNode, QparentNode, SideQbox
            C : current node of the KD-tree
            Q : parent (node) of C
            compute box using Q.box, Q.axis, Q.position, side
        """
        box = copy.copy(Q.box)
        axis = Q.axis
        box[side]
        box = Box([])
        C.box = box
        C.level = Q.level
        return
    
    """
        Queries Defined:
            Kiss dissertation page 36-37
            Juttler paper, page 12-14
    """
    def isactive(self, bounds):
        """Query1:
        all the basis functions of this level 
        whose support is completely 
        contained in the box b 
        are active
        """
        return
    
    def ispassive(self, bounds):
        """Query@:
        all the basis functions of this level 
        whose support is (completely) 
        contained in the box b 
        are passive
        """
        return
    
    def iscontained(self, bounds):
        """Query3:
        returns the highest level l with the property that
        space Omega^{l} contains the bounds box B
        
        B & Omega^{l} not Empty 
        """
        return
    
    def highestlevel(self, bounds):
        """Query4: 
        returns the maximum value L 
        for which the following
        condition is satisfied:
        
        B & Omega^{L} not Empty , L = 0,...,N-1
        """
        return
    
    
    
class THBTree(object):
    """
        Kiss: 30  ~d=2 for THB spline surfaces~
        Juttler: 10
    """
    def __init__(self):
        return



class rbspline(spline.Bspline):
    
    def __init__(self, vertices, 
                 k, nump, 
                 t=None, trange=None, 
                 bounds = ia(0.,1.)):
        
        super(rbspline, self).__init__(vertices, k, nump, t=None, trange=None )
        
        if t is None:
            self.t_copy         = kv.make_knot_vector(k,len(vertices)-1,self.trange,type='open')
        else:
            self.t_copy         = t 
        self.t                  = copy.deepcopy(self.t_copy)
        self.bounds                 = bounds 
        self.upper_boundary_knots   = None      # knots which bound the coarser level of detail for C2 continuity
        self.lower_boundary_knots   = None      # knots which bound the finer  level of detail for C2 continuity
        self.level                  = 0         # 0 is the coarsest
        self.interior_knots         = []        # ?
        self.parent                 = None      # chg to quad tree for surf
        self.children               = []        # chg to quad tree for surf
        self.active_indicies        =  self.active(bounds)#set(self.offset_reference.keys())  # portions of the spline space that have associated vertices
        self.active_prange          = bounds
    
    @property
    def vertices(self):
        return self._vertices
    @vertices.setter
    def vertices(self, *args):
        vertices =  args[0]
        self._vertices = vertices
        if np.any(self.vertices_copy !=vertices):
            if len(vertices) != self.n:
                #print 'hierarchical edit: vertices changed'
                self.n = len(vertices)
                self.t = kv.make_knot_vector(self.k, self.n-1, np.array([0.,1.]), 'open')
                self.initialize_curve_basis()
                self.initialize_basis_matrices()
                self.compute_basis()
            for obsvr in self.hierarchical_obsvrs:
                obsvr.propagate_update()
            self.compute_curve()
            self.vertices_copy = copy.deepcopy(vertices)
        return
        
    def active_knots(self, bounds=None,spans=None):
        """ tells which knots (basis) and (indirectly control points)
        are active on this level of the hierarchy
        """
        if bounds is None:
            bi = self.bounds[0]
            be = self.bounds[1]
        else:
            bi = bounds[0]
            be = bounds[1]
        if spans is None:
            si = self.FindSpan(bi)
            se = self.FindSpan(be)
        else:
            si = spans[0]
            se = spans[1]
        k = self.k
        n = self.n
        t = self.t
        p =self.p
        return [si-p,se]
    
    def active_basis_at_u(self, u):
        return self.active_Bik_interval(u=u)
    
    
    def get_range_of_basis(self, basis_index):
        a = self.t[basis_index]
        b = self.t[basis_index+self.k]
        return ia(a,b)
    
    
    def enclosed_basis_in_range(self, bounds=None):
        if bounds is None: bounds = self.bounds
        assert(isinstance(bounds, ia))
        active_knots = self.active_knots(bounds)
        active_ranges = []
        contained_basis = []
        indices = []
        for basis_index in range(active_knots[0],active_knots[1]+1):
            indices.append(basis_index)
            active_ranges.append(
                    self.get_range_of_basis(basis_index))
        
        for basis_index, active_range in zip(indices,active_ranges):
            if active_range in bounds:
                contained_basis.append(basis_index)
        return contained_basis
        
        
    
    def dyadic_refinement(self):
        """
            -Builds the basis space 1 level finer
            -By computing knots at each midpoint of the knot interval
            -Puts a new curve on this basis space
        """
        t           = self.t
        k           = self.k
        n           = self.n
        Nseg        = n-k+1     # number of k order bezier segments in the curve
        #nknots = np.zeros((Nseg,1),float)
        rspline            = copy.deepcopy(self)
        rspline.base_curve = False
        for i in range(Nseg):
            #nknots[i] = (t[k+i]+t[k+i+1])/2.
            iknot = (t[k+i-1]+t[k+i])/2.
            tnew, vertices = ko.knot_insertion(rspline.k,
                                               rspline.vertices,
                                               rspline.t,
                                               iknot)
            #rspline.t           = tnew                             
            rspline.vertices    = vertices
            rspline.t           = tnew
        self.children = rspline
        rspline.parent = self
        return rspline
        
    
    def bounds_refinement(self, bounds, knots=None):
        """
            old, interesting, not used.
        """
        nbpoints     = 4
        if knots:# is not None:
            nipoints     = len(knots)
        else:
            nipoints = 0
        k           = self.k
        p           = self.p
        n           = self.n
        bi          = bounds[0]
        be          = bounds[1]
        si          = self.FindSpan(bi)
        se          = self.FindSpan(be)
        relknots    = range(si,se+2)
        Nseg        = n-k+1
        rspline                     = copy.deepcopy(self)
        rspline.level               = self.level + 1
        self.upper_boundary_knots   = [bi,bi,be,be]
        for i in range(nbpoints):
            #nknots[i] = (t[k+i]+t[k+i+1])/2.
            iknot = self.upper_boundary_knots[i] #(t[k+i-1]+t[k+i])/2.
            tnew, vertices, ir = ko.knot_insertion_hierarchical(
                                                   rspline.k,
                                                   rspline.vertices,
                                                   rspline.t,
                                                   iknot)
            #print 'vertices {}:{} will change'.format(ir-p,ir)
            #print vertices
            #rspline.t           = tnew
            rspline.vertices    = vertices
            rspline.t           = tnew
            rspline.lower_boundary_knots = [bi,bi,be,be]
        
        rspline.bounds = bounds
        rspline.active_prange = bounds
        rspline.active_offsets = []
        if knots is not None:
            for knot in knots:
                tnew, vertices, ir = ko.knot_insertion_hierarchical(
                                                   rspline.k,
                                                   rspline.vertices,
                                                   rspline.t,
                                                   knot)
                
                ith_kvtx = ir-rspline.k+2
                #rspline.offset_reference[ith_kvtx]    = [ith_kvtx,knot]
                rspline.active_offsets.append(ith_kvtx)
                rspline.vertices                            = vertices
                rspline.t                                   = tnew
                rspline.interior_knots.append(knot)
                
        #rspline.offsets = np.zeros(
        #                (rspline.n,rspline.dim),float)
        self.hierarchical_obsvrs.append(rspline)
        rspline.parent = self
        #rspline.active_offsets = np.asarray(rspline.active_offsets)
        rspline.active_indicies = rspline.active_knots(rspline.bounds)#set(rspline.offset_reference.keys())
        self.children.append(rspline)
        return rspline
    
    
    def isactive(self, pt):
        if pt in self.active_intervals:
            return True
        else:
            return False
    
    def get_active_curve(self, pt):
        if self.isactive(pt):
            return self
        else:
            for ch in self.children:
                return ch.get_active_curve(pt)
    
    def print_knots(self):
        for i, el in enumerate(self.t):
            print i, el
            
    def print_properties(self):
        print '# knots = {}'.format(len(self.t))
        print '# vertices = {}'.format(self.n)
        
    def bases(self, ends):
        """return a range 
        from end[0] to end[1], inclusive
        """
        #range(active_parent_basis[0],active_parent_basis[1]+1)
        return range(ends[0],ends[1]+1)
    
    def truncate_basis(self, bounds, u):
        """first show how the basic calc
        takes place

        Q1:  How do you know which basis
        funcs at the parent (V0) level
        to alter via
        parent = parent - parent*(sum (child))
        ans:  all of them which are active
        ----
        Q2:  is that formula identically 0??
        ans:  no:  the inner basis
        does not sum to 1 until
        you are ~really~ on the interior!
        """
        inner = self.children.enclosed_basis_in_range(bounds)
        outer = self.enclosed_basis_in_range(bounds)

        oBasis = self.TLMBasisFuns(u)
        iBasis = self.children.TLMBasisFuns(u)
        retBasis = np.zeros_like(iBasis)

        #truncate:
        #        active_parent_basis = self.active(u)
        #        active_child_basis = self.children.active(u)
        apb = self.bases(self.active(u))
        acb = self.bases(self.children.active(u))
        if u in bounds:
            
            """
            outer   : the indices of the parent basis functions 
                        totally enclosed in the refined space
            inner   : the indices of the child basis functions
                        totally enclosed in the refined space
            apb     : the indices of the active parent basis functions
                        at the current parametric location, u
            acb     : the indices of the active child basis functions
                        at the current parametric location, u
                        
            W*Sum(oBasis[apb] - oBasis[outer] ) + Sum(iBasis[inner])
            """
            w = 1. - sum(iBasis[inner])
            
            deleted = sum(oBasis[outer])
            oBasis[outer] = 0.
            
            
            e = sum(oBasis[apb])
            
            
            W = 0.
            
            if e !=0.:
                W = w/e
                #W = W/float(len(inner))
            #print 'u,W=',u,W
            #for i in inner:
            #    oBasis[apb] = oBasis[apb] - W*iBasis[i]
            oBasis[apb] = W*oBasis[apb]
                #for key in apb:
                #    oBasis[key] = oBasis[key] - oBasis[key]*iBasis[el]
            #for i in outer:
            #    oBasis[i] = 0.
            if False:
                print u, sum(oBasis[apb])  +sum(iBasis[inner]) 
        return oBasis
    
    

    def sum_thbasis(self, bounds, u):
        inner = self.children.enclosed_basis_in_range(bounds)
        outer = self.enclosed_basis_in_range(bounds)

        oBasis = self.TLMBasisFuns(u)
        iBasis = self.children.TLMBasisFuns(u)
        
        
        return 


    def plot_truncated_basis(self, bounds,
                             interval = ia(0.,1.),
                             nump=30, offset=0., scaler=None):
        p = self.p
        a = interval[0]
        b = interval[1]
        #bi = self.active(interval)
        s = np.linspace(a,b,nump,endpoint=True)
        #numb = bi[1] - bi[0] + 1
        basis = np.zeros((nump,self.n))
        for i,u in enumerate(s):
            #span = self.FindSpan(u)
            #self.BasisFuns(span,u,basis[i,span-p-bi[0]:span+1])
            basis[i] = self.truncate_basis(bounds,u)[:,0]
            basis[i][basis[i] < 0.] = 0.
        #plt.plot(s,)
        plt.plot(s,basis[:]-offset)
        return #basis
    
    
    def plot_a_basis(self, bounds, blist,
                             interval = ia(0.,1.),
                             nump=30, offset=0., scaler=None):
        p = self.p
        a = interval[0]
        b = interval[1]
        s = np.linspace(a,b,nump,endpoint=True)
        basis = np.zeros((nump,self.n))
        for i,u in enumerate(s):
            basis[i][blist] = self.TLMBasisFuns(u)[blist]

        plt.plot(s,basis[:]-offset)
        return 




class grid(object):
    def __init__(self,dim=2,udegree=3,vdegree=3):
        self.dim    = dim
        self.up      = udegree
        self.uk      = udegree+1
        self.vp      = vdegree
        self.vk      = vdegree+1
        self.ustart = [ 0 for el in range(self.uk)]
        self.uend = [ 1 for el in range(self.uk)]
        self.vstart = [ 0 for el in range(self.vk)]
        self.vend = [ 1 for el in range(self.vk)]
        return

if __name__ == '__main__':
    k=4
    nump=100
    start = [0.,12.]
    end = [12.,0.]
    num = 11 #makes c0c's 5 basis function be enclosed in ia(.2,.8)
    vv = spline.curvalinear_vertices(start,end,num)

    
    
    #c0 = hb.LevelSpline(vv,k,nump)rbspline
    c0 = rbspline(vv,k,nump)
    c0.print_properties()
    
    c0i = c0.bounds_refinement(bounds=ia(.75,.8),
                               knots=[.77,.775,.78])
    c0.print_properties()
    
    c0c = copy.copy(c0)
    c0c.hierarchical_obsvrs = []
    c1 = c0c.dyadic_refinement()
    c0.print_properties()
    
    #c1.t[4] = .54
    c1i = rbspline(c1.vertices,
                         k,
                         nump,
                         t=c1.t)
    
    c1i = c1i.bounds_refinement(bounds=ia(.75,.8),
                               knots=[.77,.775,.78])
#    
#    c1.plot_basis_interval([0.,1.],100,.0)
#    c0i.plot_basis_interval([.6,.9],10000,1.)
#    c0.plot_basis_interval([0.,1.],100,2.)
    
    
    self = c0c
    print c1.enclosed_basis_in_range(ia(.4,.72))
    
    
    
    print c1 is c0c.children
    print c1.parent is c0c
    
    """******************************************************
    """
    """bounds = ia(.25,.75)#"""
    #bounds = ia(.33,.67)
    bounds = ia(.2,.8)
    bounds = ia(.125,.875)
    
    bounds2 = ia(.3,.7)
    print c0c.enclosed_basis_in_range(bounds)
    print c1.enclosed_basis_in_range(bounds)
    
    def examinbasis(u):
        b0 = c0c.TLMBasisFuns(u)
        b1 = c1.TLMBasisFuns(u)
        for e1, e2 in zip(b0,b1):
            print e1, e2
        return
        
    
    b0 = c0c.TLMBasisFuns(bounds[0]+.01)
    b1 = c1.TLMBasisFuns(bounds[0]+.01)
    
    
    
    c1.plot_basis_interval([0.,1.],100,.0)
    c0c.plot_basis_interval([0.,1.],100,2.)
    
    c2 = c1.dyadic_refinement()
    
    inner = self.children.enclosed_basis_in_range(bounds)
    outer = self.enclosed_basis_in_range(bounds)
    u =.35
    u = 0.555555555556
    
    plt.close()
    
    
    stplt = 0.
    
    #total top
    c2.plot_basis_interval([0.,1.],250,stplt)
    
    mplt = stplt+1.
    #total middle
    c1.plot_basis_interval([0.,1.],250,mplt)
    
    mmplt = mplt+1.
    #middle:
    c2.plot_a_basis(bounds2,
                    blist=c2.enclosed_basis_in_range(bounds2),
                    nump=250,offset=mmplt)
    
    bbplt = mmplt+1.
    #truncated middle:
    c1.plot_truncated_basis(bounds2, nump=250, offset=bbplt)
    
    bbbplt = bbplt+1.
    #truncated bottom:
    c0c.plot_truncated_basis(bounds, nump=250, offset=bbbplt)
    
    cplt = bbbplt+1.
    #bottom:
    c0c.plot_truncated_basis(ia(.5,.51),offset=cplt,nump=1000) #none truncated
    
    bb = c0c.truncate_basis(bounds, u)
    bbb = c1.TLMBasisFuns(u)
    print sum(bbb[inner]) + sum(bb)
    
    
    u = 0.555555555556