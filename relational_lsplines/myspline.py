# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 16:04:29 2016

@author: lukemcculloch

key to a THB curve:
    -


A THB surface...

    - has dyadic refinement at all levels, 
        so all the possible basis at each level are full and symmetric
    - queries:
        set which basis are active at a given level
        
evaluation...

    Key idea:  evaluation always takes place
        as the tensor product of the
            p+1 u Basis, active at u
            p+1 v Basis, active at v
            p+1 x p+1 u Cntrl Pts(u,v)
            
        The truncated hierarchical basis just picks them out!
        
        basis spans figure out what gets truncated
        
        evaluation takes place at a single u,v pt
        as always.
        
    
THB basis plotting:
    
    rbspline_object.plot_thbasis() #plots the multilevel truncated basis
    
    so, for the examples below __main__:  c0c.plot_thbasis()
    
THB curve plotting:
    
    rbspline_object.plotcurve() # plots the multilevel curve, evaluated 
                                # in one of at least two possible THB ways
    
    so, for the examples below __main__:  c0c.plotcurve()

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
import extended_interval_arithmetic as iafuncs

import THButilities

np.set_printoptions(threshold= 40)


class Box(object):
    """
        TODO:  add method contains?
    """
    def __init__(self, dims):
        assert( np.asarray([isinstance(dim, ia) for dim in dims]).any() )
        self.dims = [dim for dim in dims]
        self.isempty = np.asarray([dim.isempty for dim in dims]).any()
        

    def __repr__(self):
        return 'Box({})'.format(self.dims)

    def __str__(self):
        return 'Box({})'.format(self.dims)
    
    def __getitem__(self, index):
        return self.dims[index]

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
#
#    def __in__(self, other):
#        """is other in self?
#        """
#        check_dim = []
#        for ds, do in zip(self.dims, other.dims):
#            check_dim.append( ds in do )
#        iscontained = np.asarray([el.isempty for el in check_dim]).any()
#        return iscontained
    def contains(self, u,v):
        """is other in self?
        """
        check_dim = []
        check_dim.append( self[0].contains(u))
        check_dim.append( self[1].contains(v))
        iscontained = np.asarray([el.isempty for el in check_dim]).any()
        return iscontained

class BoundsList(object):
    """
    Bounds used in THBspline Curves 
        
        DONE: add method contains?
        TODO: add disjunctive _or_
    """
    def __init__(self, bounds):
        # stupid
        #type issue: ia /= ia
        # occuring
        #assert( np.asarray([isinstance(el, ia) for el in bounds]).any() ),'type {}'.format(type(el))
        pbounds = [el for el in bounds]
        self.isempty = np.asarray([el.isempty for el in bounds]).any()
        #if len(self.bounds)>1:
        self._bounds = self.revise_make_simply_connected(pbounds)
        

    def __getitem__(self, index):
        return self.bounds[index]

    def __repr__(self):
        return 'BoundsList({})'.format(self.bounds)

    def __str__(self):
        return 'BoundsList({})'.format(self.bounds)
    
    @property
    def bounds(self):
        return self._bounds
    @bounds.setter
    def bounds(self, *args):
        bounds =  args[0]
        self._bounds = self.revise_make_simply_connected(bounds)
        return
        

    def __and__(self, other):
        """Intersection of two boxes
        """
        new_dims = []
        for ds, do in zip(self.bounds, other.bounds):
            new_dims.append( ds & do )
        return BoundsList(new_dims)
    
    
    
    def make_disjoint(self):
        """DEV: never use
        """
        e1 = self.bounds[0]
        dl = [e1]
        for el in self.bounds[1:]:
            nbl1 = [ bi.exclusive_disjunctive_elimination(el) for bi in dl]
            nbl2 = [ el.exclusive_disjunctive_elimination(bi) for bi in dl]
            nbl3 = [ bi&el for bi in dl ]
            dl = nbl1+nbl3 +nbl2
        dl = [el for el in dl if not el.isempty]
        return
    
    
    def revise_make_simply_connected(self, boundslist):
        save = copy.copy(boundslist)
        L1 = len(save)
        new = self.make_simply_connected(boundslist)
        if len(new)<L1:
            return self.revise_make_simply_connected(new)
        else:
            return new
        
    def make_simply_connected(self, boundslist=None):
        """given a bounds list, make sure there are
        no overlaps
        
        find some efficiencies here, later
        
        dev
        ----------
        
        overlapping bounds:
        boundslist = BoundsList([ia(0.03125, 0.09375),
                                ia(0.0625, 0.125),
                                ia(0.09375, 0.15625),
                                ])
    
        good bounds:
        boundslist = BoundsList([ia(.1,.28),
                                ia(.3,.5),
                                ia(.51,.72)])
    
        """
        if boundslist is None:
            boundslist = self.bounds
            
        def doit(startlist, searchlist):
            nbl = startlist
            nextbl = []
            for el in searchlist:
                dl = []
                while nbl:
                    bi = nbl.pop()
                    thing = el.smart_or(bi)
                    if isinstance(thing, tuple):
                        thing = list(thing)
                    else:
                        thing = [thing]
                    for new_el in thing:
                        if new_el not in dl:
                            dl.append(new_el)
                nextbl +=dl
            nextbl += startlist
            return nextbl
            
        startlist = [boundslist[0]]
        searchlist = boundslist[1:]
        nbl = doit(startlist, searchlist)
        
        
        for el in boundslist:
            nbl = doit([el],nbl)
        for el in boundslist:
            nbl = doit(nbl,[el])
        
        
        nbl.sort()
        return nbl
        
    def __or__(self, other):
        """Union of two boxes
        """
        new_dims = []
        for ds, do in zip(self.bounds, other.bounds):
            new_dims.append( ds | do )
        return BoundsList(new_dims)

    def add(self, interval):
        self.bounds.append(interval)
        return

    def contains(self, other):
        """is other in self?
        """
        ans = None
        if isinstance(other, BoundsList):
            check_dim = []
            for ds, do in zip(self.bounds, other.bounds):
                check_dim.append( ds.contains(do) )
            iscontained = np.asarray([not el.isempty for el in check_dim]).any()
            return iscontained
        elif isinstance(other, float):
            #print 'Hello!', self, other
            check_dim = []
            for ds in self.bounds:
                check_dim.append( ds.contains(other) )
            ans = [i for i, x in enumerate(check_dim) if x]
            return ans #only indicies which contain the u location

    def exlude_from_bounds(self, other):
        """exclude boundslist 2 from boundslist1
        """
        return iafuncs.interval_list_diff(self.bounds,other.bounds)

class KDnode(object):
    """
        Kiss: 30  ~d=2 for THB spline surfaces~
        Juttler: 10
    """
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
        (Kiss one by one approach, page 39)
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



class THBhierarchy(object):
    """
        Kiss: 38
        
        idea: take the basis functions out of the curve...
        
        but get them from the curve originally...
        
        
        a dyadically refined surface need not go around 
        generating a new longitudinal curve 
        for every dyadicly refined transverse curve
    """
    def __init__(self, root):
        self.root = root
        return

    def __getitem__(self, index):
        return

    def setlist(self, Q, X):
        """
            inputs:
                Q : root of the KD tree
                X : X[L] list or array of active functions on level L
        """
        for L in X:
            #TODO: create the set I for all basis of level L
            #       acting on Omega^{0}
            self.setbox(I, X[L])
        return

    def setbox(self, T, XL):
        """Kiss page 39-40
            inputs:
                I : box in the knot index space
                XL : list of active functions on level L

            L is a global variable??
        #TODO: create the axis aligned bounding box B
        #covering all cells of
        #level L which belong
        #to the supporst with indicies in I
        """

        return


class THBCurveLists(object):
    def __init__(self, tcurvelist, lcurvelist):
        self.tcurvelist = tcurvelist
        self.lcurvelist = lcurvelist
    
    def generate_arrays(self):
        return
    
    def access(self):
        return

class LevelData(object):
    def __init__(self, level, nu, nv,
                 ubounds, vbounds):
        self.level = level
        self.nu = nu
        self.nv = nv
        self.ubounds = ubounds
        self.vbounds = vbounds
        self.child = None
        self.parent = None

#class THBsurface(object):
#    """NOT USED!
#    
#        Note that with surfaces the parent child
#        hierarchy is no longer sufficient
#        Because under dyadic refinement,
#        new curves are created in the gaps 
#        between each directions curves at a higher level.
#        
#        There is not clear representation of which child went
#        with which parent level by level.
#        
#        basis support must be divorced from that old concept
#        
#        -----------!OR:
#        
#        use the master curves to derive the active basis 
#        in parent from child!?
#        
#        parents point all children at the level up master curves.
#        
#    """
#    """
#    parsimonious representation (1 curve to stand in for 
#                                 entire dyadic level
#                                 possible due to 
#                                 compound bounds):
#        
#        self.SCl['u'] = {   0: <class curve.Bspline>, 
#                            1: <class curve.Bspline>
#                            }
#        
#        self.SCl['v'] = {   0: <class curve.Bspline>
#                            }
#        
#    """
#    """
#    TODO:
#        1.)  fix initializer so that transverse hull curves
#        become lofted 
#        2.)  make sure all levels of the surface represent the
#        (numericaly exact) same geometry
#    """
#    def __init__(self, tcurvelist=None, lcurvelist=None):
#        self.tcurvelist = tcurvelist
#        self.lcurvelist = lcurvelist 
#        #self.lu = len(tcurvelist)
#        #self.lv = len(lcurvelist)
#        self.tlevels = {} #
#        self.ulevels = {} #
#        self.vlevels = {} #
#        self.uHLevels = {}#
#        self.vHLevels = {}#
#        self.Mlevels = {} #
#        self.level_vertices = {}
#        self.realbox = {} #
#        self.ubasis_boundsboxs = {}
#        self.vbasis_boundsboxs = {}
#        self.ubounds = BoundsList([ia(0.,1.)])
#        self.vbounds = BoundsList([ia(0.,1.)])
#        self.level_data_list = []
#        
#        
#    def sound_hierarchy(self):
#        """Depth should not matter?
#        """
#        depth_t = 0
#        depth_l = 0
#        for el in self.tcurvelist:
#            depth_t = max(depth_t, el.sound_hierarchy())
#        for el in self.lcurvelist:
#            depth_l = max(depth_l, el.sound_hierarchy())
#        #        if depth_t < depth_l:
#        #            diff = depth_l-depth_t
#        #            for el in range(diff):
#        #                
#        return
#    
#    def makebounds(self, ia_interval):
#        return BoundsList([ia_interval])
#        
#    def set_coefficients(self, mesh):
#        for i,el in enumerate(mesh):
#            self.c = el
#        return
#    
#    
#    
#    def eval1thb(self, u,v, ans = 0., level=0):
#        """
#            for the curve version see 
#            def rthbeval(self, u):
#        
#         unrefined 
#         sapce
#         over here
#        .       .    .
#        |       |    |
#        .       .    .   refined space lurks below over here
#        |       |    |
#        Bu_{i}  i1   Bu_{i+1}  i2
#        
#        if i1 is in between two knots points,
#        and one is refined and the other is not
#        then we are ~not enclosed~ on the refined basis
#        take the non-refined answer 
#        for
#            S = A U R
#            we would be in A, not R
#            we would be in the unrefined space
#        """
#        #udepth = self.nHLevels['u'], u,v, level):
#        """
#        
#        >>> gg.surface_point(.5,.5,0)
#        array([[ 6.,  0.]])
#        
#        >>> gg.HL['v'][0]
#        <class curve.Bspline>
#        >>> gg.HL['v'][0].eval(.5)
#        
#        array([[ 0.,  6.]])
#        >>> gg.HL['u'][0].eval(.5)
#        array([[ 6.,  0.]])
#        
#        """
#        """
#        for now, use level calls to get what you would
#        get if you had been in a self hosting tree
#        like THBcurve below, 
#        I mean -rbspline- below
#        """
#        ucurve,vcurve = self.get_level_master_curves(level)
#        p = ucurve.p
#        q = vcurve.p
#        uspan = ucurve.FindSpan(u)
#        vspan = vcurve.FindSpan(v)
#        #uBasisFuns = ucurve.TLMBasisFuns(u) #Nu  TODO:  change to truncated basis.....
#        #vBasisFuns = vcurve.TLMBasisFuns(v) #Nv
#        
#        uBasisFuns,ssf_u = ucurve.truncate_basis2(u) #Nu  TODO:  change to truncated basis.....
#        vBasisFuns,ssf_v = vcurve.truncate_basis2(v) #Nv
#        
#        
#        P = self.get_level_vertices(level)
#        
#        uind   = uspan-p
#        #ushape = np.shape(uBasisFuns)
#        #vshape = np.shape(vBasisFuns)
#        #Sarray = np.zeros((ushape[0],ushape[1],3))
#        S = 0.
#        for l in range(q+1):
#            temp = 0.
#            vind = vspan-q+l
#            #Sarray[l] = np.dot( uBasisFuns.T, P)
#            for k in range(0,p+1):
#                #print 'P[{},{}] = {}'.format(vind,uind+k,P[vind,uind+k])
#                temp += uBasisFuns[k]*P[uind+k][vind]#[uind+k][vind] #P[vind][uind+k] #
#            S += vBasisFuns[l]*temp
#        #S = np.dot(vBasisFuns.T, Sarray)
#        return S
#    
#    def evalthb(self, u,v):
#        S = 0.
#        for el in gg.level_vertices:
#            S = S + self.eval1(u,v,ans=0.,level=el)
#        return S
#    
#    
#    
#    def eval1(self, u,v, ans = 0., level=0):
#        """
#         unrefined 
#         sapce
#         over here
#        .       .    .
#        |       |    |
#        .       .    .   refined space lurks below over here
#        |       |    |
#        Bu_{i}  i1   Bu_{i+1}  i2
#        
#        if i1 is in between two knots points,
#        and one is refined and the other is not
#        then we are ~not enclosed~ on the refined basis
#        take the non-refined answer 
#        for
#            S = A U R
#            we would be in A, not R
#            we would be in the unrefined space
#        """
#        #udepth = self.nHLevels['u'], u,v, level):
#        """
#        
#        >>> gg.surface_point(.5,.5,0)
#        array([[ 6.,  0.]])
#        
#        >>> gg.HL['v'][0]
#        <class curve.Bspline>
#        >>> gg.HL['v'][0].eval(.5)
#        
#        array([[ 0.,  6.]])
#        >>> gg.HL['u'][0].eval(.5)
#        array([[ 6.,  0.]])
#        
#        """
#        """
#        for now, use level calls to get what you would
#        get if you had been in a self hosting tree
#        like THBcurve below, 
#        I mean -rbspline- below
#        """
#        ucurve,vcurve = self.get_level_master_curves(level)
#        p = ucurve.p
#        q = vcurve.p
#        uspan = ucurve.FindSpan(u)
#        vspan = vcurve.FindSpan(v)
#        #uBasisFuns = ucurve.TLMBasisFuns(u) #Nu  TODO:  change to truncated basis.....
#        #vBasisFuns = vcurve.TLMBasisFuns(v) #Nv
#        
#        uBasisFuns,ssf_u = ucurve.truncate_basis2(u) #Nu  TODO:  change to truncated basis.....
#        vBasisFuns,ssf_v = vcurve.truncate_basis2(v) #Nv
#        
#        
#        P = self.get_level_vertices(level)
#        
#        uind   = uspan-p
#        #ushape = np.shape(uBasisFuns)
#        #vshape = np.shape(vBasisFuns)
#        #Sarray = np.zeros((ushape[0],ushape[1],3))
#        S = 0.
#        for l in range(q+1):
#            temp = 0.
#            vind = vspan-q+l
#            #Sarray[l] = np.dot( uBasisFuns.T, P)
#            for k in range(0,p+1):
#                #print 'P[{},{}] = {}'.format(vind,uind+k,P[vind,uind+k])
#                temp += uBasisFuns[k]*P[uind+k][vind]#[uind+k][vind] #P[vind][uind+k] #
#            S += vBasisFuns[l]*temp
#        #S = np.dot(vBasisFuns.T, Sarray)
#        return S
#    
#    def eval(self, u,v):
#        S = 0.
#        for el in gg.level_vertices:
#            S = S + self.eval1(u,v,ans=0.,level=el)
#        return S
#        
#    def get_3D_instantiation_from_1part_shape(self, array,dim=3):
#        array_shp = np.shape(array)
#        return np.zeros((array_shp[0],array_shp[1],dim))
#    
#    def initialize_grid(self, 
#                        ku=4, kv=4, 
#                        nx=4, ny=4, nz=4, 
#                        nump=30, box=None):
#        """
#            TODO: fix this to take transverses initially
#            and store 
#            -those transverses,
#            -longitudinal loft curves
#            -and transverse surface curves to match
#            
#        """
#        if box is None:
#            box = Box([ia(0.,12.),
#                       ia(0.,12.), 
#                       ia(0.,12.)])
#    
#        self.ku = ku
#        self.kv = kv
#        self.realbox[0] = box
#        
#        self.ubasis_boundsboxs[0] = Box([ia(0.,1.),ia(0.,1.)])
#        self.vbasis_boundsboxs[0] = Box([ia(0.,1.),ia(0.,1.)])
#        
#        
#        self.lu = nx
#        self.lv = ny
#    
#        self.x = np.linspace(box[0][0],box[0][1],nx,endpoint=True)
#        self.y = np.linspace(box[1][0],box[1][1],ny,endpoint=True)
#        self.z = np.linspace(box[2][0],box[2][1],nz,endpoint=True)
#        self.mesh = np.outer(self.x,self.y)
#        
#        mx,my = np.meshgrid(self.x,self.y)
#        mz = np.zeros_like(mx)
#        self.mx = mx
#        self.my = my
#        self.mz = mz
#        
#        """Set the surface vertices:
#            
#        TLM, for some reason I commented these:
#        #array_shp = np.shape(self.mx)
#        #array_3d = np.zeros((array_shp[0],array_shp[1],3))  #array_3d[0,0][x1,x2,x3]
#        
#        uncommenting below:
#        """
#        array_shp = np.shape(self.mx)
#        array_3d = np.zeros((array_shp[0],array_shp[1],3))  #array_3d[0,0][x1,x2,x3]
#        array_3d = self.get_3D_instantiation_from_1part_shape(self.mx)
#        for i in range(array_shp[0]):
#            for j in range(array_shp[1]):
#                array_3d[i,j] = [self.mx[i,j],self.my[i,j],self.mz[i,j]]
#                
#        self.level_vertices[0] = array_3d
#        """Done
#        """
#        self.child = None
#        
#        ptsu = np.asarray([mx[0],
#                           np.zeros_like(my.T[0])]).T
#    
#        ptsv = np.asarray([np.zeros_like(mx[0]), 
#                           my.T[0] ]).T
#    
#        self.uHLevels[0] = rbspline(ptsu, ku, nump)
#        self.vHLevels[0] = rbspline(ptsv, kv, nump)
#        
#        
#        self.HL = {'u':self.uHLevels,
#                    'v':self.vHLevels}
#        
#        
#        
#        
#        #
#        #******************************** start new loft procedure
#        #
#        self.tlevels[0] = []
#        self.ulevels[0] = []
#        self.vlevels[0] = []
#        #store real boat refined transverse curves:
#        for i in range(self.lu):
#            pts = np.asarray([mx[i],my[i],mz[i]]).T
#            self.tlevels[0].append(rbspline(pts, ku,nump))
#        
#        #loft them longitudinally. 
#        for i in range(self.lv):
#            pts = np.asarray([mx[:,i],my[:,i],mz[:,i]]).T
#            self.vlevels[0].append(
#                    rbspline(pts, ku,nump,interpolate=True))
#            
#        # now set vertices.... for fast surface eval later.
#        nv_pts = self.vlevels[0][0].n
#        nu_pts = len(self.vlevels[0])
#        vertices = np.zeros((nu_pts,nv_pts,3),float) #possibly backwards
#        for row_i, curve in enumerate(self.vlevels[0]):
#            vertices[row_i,:] =  curve.vertices
#            
#        #now finalize the transverse net.
#        l2_ucurves = []
#        for j_col in range(nv_pts):
#            new_curve = rbspline(vertices[:,j_col],
#                                 k,
#                                 nump)
#            #new_curve.bounds = umc.bounds
#            #new_curve.parent = umc
#            l2_ucurves.append(new_curve)
#        #store real transverses:
#        self.ulevels[0] = l2_ucurves  
#        
#        #
#        #******************************** end new loft procedure
#        #
#            
#        self.Mlevels[0] = np.ones((self.lu,self.lv), int)
#        self.unlevels = 0
#        self.vnlevels = 0
#        self.nHLevels = {'u':self.unlevels,
#                        'v':self.vnlevels}
#        self.level = 0
#        
#        self.level_data_list.append(
#                                    LevelData(level=0,
#                                      nu = 4,
#                                      nv = 4,
#                                      ubounds = Box([ia(0.,1.),ia(0.,1.)]),
#                                      vbounds = Box([ia(0.,1.),ia(0.,1.)])  
#                                      ) )
#        return
#    
#    #    def refine_grid(self, dir_, 
#    #                    nL_refine_depth=1, 
#    #                    ubounds=None,
#    #                    vbounds=None):
#    def refine_grid(self, 
#                    nL_refine_depth=1, 
#                    ubounds=None,
#                    vbounds=None):
#        """
#            dir_            : character, u or v 
#            self.nHLevels   : {'u': 1, 'v': 1} 
#                                (current max depth in each dir_)
#            bounds          : BoundsList object
#            nL_refine_depth : number of levels (depth) 
#                                to add at once
#                                
#            Only Symmetric refinement makes sense
#            for a tensor product surface
#            of level l
#        """
#        """
#            TODO: modify dyadic refinement
#            to return knots and control points.
#            assemble all the new control points
#            so as to catch the new vertices...
#            
#            or just add surface knot insertion
#            and especially
#            surface knot refinement
#            
#            to the class
#        """
#        """----------------------------------------------
#        """
#        dir_ = 'u'
#        #first master curve in u at this level
#        old_level = self.nHLevels[dir_]
#        new_level = old_level +1
#        
#        old_ucurves,old_vcurves = self.get_level_all_curves(old_level)
#        
#        sc = self.HL[dir_]  #sc_u, sc_v = self.get_level_base_curves(old_level)
#        sc[new_level] = sc[old_level].dyadic_refinement()
#        if bounds is not None:
#            sc[new_level].bounds = sc[old_level].bounds
#        else:
#            sc[new_level].bounds = bounds
#        self.nHLevels[dir_] +=1
#        
#        #now the real curves:
#        sc[new_level].bounds
#        self.ulevels[new_level] = []
#        for i in range(self.lu):
#            lup_u_thbspline = self.ulevels[old_level][i].dyadic_refinement()
#            lup_u_thbspline.bounds = sc[new_level].bounds
#            self.ulevels[new_level].append(lup_u_thbspline)
#        
#        """----------------------------------------------
#        """
#        dir_ = 'v'
#        #first master curve in v at this level
#        sc = self.HL[dir_]
#        sc[new_level] = sc[old_level].dyadic_refinement()
#        if bounds is not None:
#            sc[new_level].bounds = sc[old_level].bounds
#        else:
#            sc[new_level].bounds = bounds
#        self.nHLevels[dir_] +=1
#        
#        #now the real curves:
#        sc[new_level].bounds
#        self.vlevels[new_level] = []
#        for i in range(self.lv):
#            lup_v_thbspline = self.vlevels[old_level][i].dyadic_refinement()
#            lup_v_thbspline.bounds = sc[new_level].bounds
#            self.vlevels[new_level].append(lup_v_thbspline)
#        
#        #TODO:  do this the right way
#        # interleaving does not work because the number 
#        # of points inserted varies by the curve.k and curve.n!
#        #
#        #now interleve the new control vertices as u and v curves:
#        #        for i in range(lup_u_thbspline.n):
#        #            # dyadic knot refinement inserts uncertain number
#        #            # of control points
#        #            pass
#        #new_pts = []
#        #for curve in ulevels[old_level]:
#            
#        # self.child = new pts array?
#        
#        """----------------------------------------------
#        Store level data
#        #"""
#        nu = lup_u_thbspline.n
#        nv = lup_v_thbspline.n
#        gg.level_data_list.append(
#                            LevelData(level=new_level,
#                                      nu=nu,
#                                      nv=nv,
#                                      ubounds = Box([ia(0.,1.),ia(0.,1.)]),
#                                      vbounds = Box([ia(0.,1.),ia(0.,1.)])  
#                                      ) )
#        return
#    
#    
#    
#    
#    def dyadic_loft_refinement(self, 
#                               level=None,
#                               ubounds=BoundsList([ia(0.0, 1.0)]),
#                                vbounds=BoundsList([ia(0.0, 1.0)]) ):
#        """
#        input: 
#            level : highest level to be refined
#                    this level already exists!
#            
#        use this!  --or something very like it, I think!
#        """
#        
#        """TODO:  make interior real curve
#        aware of their THB basis!
#        
#        mapping problem:  dyadic refinement in u
#        gives new curves in v that don't belong to 
#        previous curves in v
#        
#        give them a u location based on the knot loc
#        in the u space, the knot that they -represent-?
#        """
#        if level is None:
#            old_level = self.nHLevels['u']
#            assert(old_level == self.nHLevels['v']),'u and v refinement depth !=='
#            level = old_level
#            
#            
#        self.master_curve_refinement(level,ubounds,vbounds) #refine masters
#        umc, vmc = self.get_level_master_curves(level+1) #local refined masters to be emulated
#        
#        
#        #refine real transverses:
#        ucurves, vcurves = self.get_level_all_curves(level)
#        l2_ucurves = []
#        for curve in ucurves:
#            new_curve = curve.dyadic_refinement()
#            new_curve.bounds = umc.bounds
#            new_curve.parent = umc
#            l2_ucurves.append(new_curve) #refine in u
#            
#        #store real boat refined transverse curves:
#        self.tlevels[level+1] = l2_ucurves     #boat transverses!  not surface control vertices!
#        
#        
#        #loft them longitudinally.  Refine the loft
#        loft_set = self.loft_through_transverse(l2_ucurves)
#        l2_vcurves = []
#        for curve in loft_set:
#            new_curve = curve.dyadic_refinement()
#            new_curve.bounds = vmc.bounds
#            new_curve.parent = vmc
#            l2_vcurves.append(new_curve)
#        self.vlevels[level+1] = l2_vcurves
#        
#        # now set vertices.... for fast surface eval later.
#        nv_pts = l2_vcurves[0].n
#        nu_pts = len(l2_vcurves)
#        vertices = np.zeros((nu_pts,nv_pts,3),float) #possibly backwards
#        for row_i, curve in enumerate(l2_vcurves):
#            vertices[row_i,:] =  curve.vertices
#        
#        #now finalize the transverse net.
#        l2_ucurves = []
#        for j_col in range(nv_pts):
#            new_curve = rbspline(vertices[:,j_col],
#                                 k,
#                                 nump)
#            new_curve.bounds = umc.bounds
#            new_curve.parent = umc
#            l2_ucurves.append(new_curve)
#        #store real transverses:
#        self.ulevels[level+1] = l2_ucurves  
#        
#        self.level_vertices[level+1] = vertices
#        return
#    
#    
#    def master_curve_refinement(self, 
#                               level,
#                               ubounds=BoundsList([ia(0.0, 1.0)]),
#                                vbounds=BoundsList([ia(0.0, 1.0)]) ):
#        umaster_curve, vmaster_curve = self.get_level_base_curves(level)
#        
#        
#        
#        unew_curve = umaster_curve.dyadic_refinement()
#        unew_curve.bounds = ubounds
#        unew_curve.parent = umaster_curve
#        umaster_curve.children = unew_curve
#        self.nHLevels['u'] +=1
#        
#        
#        vnew_curve = vmaster_curve.dyadic_refinement()
#        vnew_curve.bounds = vbounds
#        vnew_curve.parent = vmaster_curve
#        vmaster_curve.children = vnew_curve
#        self.nHLevels['v'] +=1
#        
#        self.set_level_master_curves(unew_curve, vnew_curve, 
#                                     ubounds,vbounds,
#                                     level+1)
#        return
#    
#    
#    def loft_through_transverse(self, tcurve_set):#, unew_curve):
#        """
#        use this!  --or something like it
#        
#        vertices[:,j] = the transvere pts in column j
#        
#        lvertices[i,:] = the longitudinal pts in row i
#        """
#        nu_pts = tcurve_set[0].n
#        nv_pts = len(tcurve_set)
#        vertices = np.zeros((nu_pts,nv_pts,3),float) 
#        
#        #vertices[:,j] = the transvere pts in column j:
#        loft_set = []
#        for col_j, curve in enumerate(tcurve_set):
#            vertices[:,col_j] = curve.vertices 
#
#        
#        #lvertices[i,:] = the longitudinal pts in row i:
#        k = curve.k
#        nump = curve.nump
#        for row_i in range(nu_pts):
#            loft_set.append(rbspline(vertices[row_i,:],
#                                     k,
#                                     nump,
#                                     interpolate=True))
#        return loft_set
#    
#    
#    
#    #    def refine_grid_u(self, 
#    #                    nL_refine_depth=None, 
#    #                    bounds=1):
#    #        self.refine_grid('u',nL_refine_depth,bounds)
#    #        return
#    #    
#    #    def refine_grid_v(self, 
#    #                    nL_refine_depth=None, 
#    #                    bounds=1):
#    #        self.refine_grid('v',nL_refine_depth,bounds)
#    #        return
#    def active(self, span,p):
#        return range(span-p,span+1)
#    
#    
#        
#    def return_transverse_vertices(self, level, index):
#        return self.level_vertices[level][index,:]
#        
#    def return_longitudinal_vertices(self, level, index):
#        return self.level_vertices[level][:,index]
#        
#    
#    def get_level_all_curves(self, level):
#        return self.ulevels[level], self.vlevels[level]
#    
#    def get_level_base_curves(self, level):
#        return self.get_level_master_curves(level)
#    
#    def get_level_master_curves(self, level):
#        return self.HL['u'][level], self.HL['v'][level]
#    
#    def set_level_master_curves(self, 
#                                ucurve, 
#                                vcurve,
#                                ubounds, 
#                                vbounds,
#                                level):
#        ucurve.bounds = ubounds
#        vcurve.bounds = vbounds
#        self.HL['u'][level] = ucurve
#        self.HL['v'][level] = vcurve
#        return
#        
#    
#    def get_level_vertices(self, level):
#        """Adding a level of indirection here
#        because I am not sure how I want these
#        stored and accessed.
#        """
#        #get the vertices
#        vertices = self.level_vertices[level]
#        return vertices
#    
#    def surface_point(self, u,v, level):
#        """
#        
#        >>> gg.surface_point(.5,.5,0)
#        array([[ 6.,  0.]])
#        
#        >>> gg.HL['v'][0]
#        <class curve.Bspline>
#        >>> gg.HL['v'][0].eval(.5)
#        
#        array([[ 0.,  6.]])
#        >>> gg.HL['u'][0].eval(.5)
#        array([[ 6.,  0.]])
#        
#        """
#        ucurve,vcurve = self.get_level_base_curves(level)
#        p = ucurve.p
#        q = vcurve.p
#        uspan = ucurve.FindSpan(u)
#        vspan = vcurve.FindSpan(v)
#        uBasisFuns = ucurve.TLMBasisFuns(u) #Nu  TODO:  change to truncated basis.....
#        vBasisFuns = vcurve.TLMBasisFuns(v) #Nv
#        
#        
#        P = self.get_level_vertices(level)
#        
#        uind = uspan-p
#        ushape = np.shape(uBasisFuns)
#        vshape = np.shape(vBasisFuns)
#        Sarray = np.zeros((ushape[0],ushape[1],3))
#        S = 0.
#        for l in range(q+1):
#            temp = 0.
#            vind = vspan-q+l
#            #Sarray[l] = np.dot( uBasisFuns.T, P)
#            for k in range(0,p+1):
#                temp += uBasisFuns[k]*P[uind+k][vind]
#            S += vBasisFuns[l]*temp
#        #S = np.dot(vBasisFuns.T, Sarray)
#        return S
#    
#    def DumbKnotInsertion(self, curve, u):
#        new_curve = curve.dyadic_refinement()
#        return new_curve.vertices, new_curve.t
#    
#    def SurfaceKnotIns(self, u, level=None):
#        """
#        input:
#            UP : u direction knots
#            VP : v direction knots
#            Pw : control vertices of the surface
#            
#            r  : number of times to insert the knot
#            
#            k  : order of the direction u or v, depending 
#            which direction we are inserting in this time?
#            
#        output:
#            UQ : u direction knots
#            VQ : v direction knots 
#            Qw : control vertices of the surface
#            
#        """
#        ucurve,vcurve = self.get_level_base_curves(level)
#        if level is None:
#            bot = u
#            level
#        
#        p = self.ku-1
#        #ku = self.ku
#        UP = ucurve.t
#        np = ucurve.n
#        
#        q = self.kv-1
#        #kv = self.kv
#        VP = vcurve.t
#        mp = vcurve.n
#        
#        alpha = []
#        
#        dir_ = 'u'
#        k = ucurve.k
#        if dir_ == 'u':
#            #copy u knot vector as 
#            #UQ = np.copy(UP)
#            #UQ = np.empty_like(UP)
#            #UQ[:] = UP
#            #UQ = UP.copy()
#            UQ = []
#            VQ = []
#            for j in range(1,r+1):
#                L = k-p+j
#        return
#        
#    
#    def fill_level(self):
#        """fill in the sparse representation
#        with dyadic points
#        """
#        return
#    
#    
#    def surface_knot_insertion(self):
#        #p = 
#        return
#    
#    
#    
#    def plotpts(self, up=None, vp=None):
#        if up is None and vp is None:
#            up = self.ulevels[0]
#            vp = self.vlevels[0]
#        for el in up:
#            el.plotcurve()
#        for el in vp:
#            el.plotcurve()
#        return
#    
#    def plot_all_levels(self):
#        
#        for el in self.ulevels:
#            umc = self.HL['u'][el]
#            for uel in self.ulevels[el]:
#                assert(uel.bounds == umc.bounds),'umc bounds != uel bounds!'
#                uel.plotcurve()
#        for el in self.vlevels:
#            vmc = self.HL['v'][el]
#            for vel in self.vlevels[el]:
#                assert(vel.bounds == vmc.bounds),'vmc bounds != vel bounds!'
#                vel.plotcurve()
#        return
#    
#    def plot_THBcurves(self, level=0):
#        """TODO:  make interior real curve
#        aware of their THB basis!
#        
#        mapping problem:  dyadic refinement in u
#        gives new curves in v that don't belong to 
#        previous curves in v
#        
#        give them a u location based on the knot loc
#        in the u space, the knot that they -represent-?
#        """
#        el = level
#        umc = self.HL['u'][el]
#        for uel in self.ulevels[el]:
#            assert(uel.bounds == umc.bounds),'umc bounds != uel bounds!'
#            uel.plotcurve()
#        el = level
#        vmc = self.HL['v'][el]
#        for vel in self.vlevels[el]:
#            assert(vel.bounds == vmc.bounds),'vmc bounds != vel bounds!'
#            vel.plotcurve()
#        return
#    
#    def refine(self, box):
#        return


## End of THBsurface
##***************************************
##





## 
##***************************************
## Start of THBspline curve:
    

class rbspline(spline.Bspline):
    """THBspline Curves
    recursive (refined) basis spline
    this curve object is its own hierarchy
    my thb surfaces are set up the same way in thbsurface
    
    -ammended to actually use the inactive L+1 basis
    to compute the active L basis directly!
    -as opposed to subtracting the contribution of
    the active L+1 basis from the L basis to get
    the active L basis!
    (both are found in the LR spline literature
    only the addative (project inactive, no need to subtract)
    method is seen in the THB literature)
    (there is a paper out there that highlights
    efficiencies of each in differing circumstances)
    """
    def __init__(self, 
                 vertices,
                 k, 
                 nump,
                 t=None, 
                 trange=None,
                 bounds = BoundsList([ia(0.,1.)]),
                 level = 0,
                 interpolate=False,
                 rm = None,
                 rmi=None):

        super(rbspline, self).__init__(vertices, 
             k, nump, t=None, trange=None )

        if t is None:
            self.t_copy         = kv.make_knot_vector(k,len(vertices)-1,self.trange,type='open')
        else:
            self.t_copy         = t
        self.t                  = copy.deepcopy(self.t_copy)
        self.bounds                 = bounds
        self.upper_boundary_knots   = None      # knots which bound the coarser level of detail for C2 continuity
        self.lower_boundary_knots   = None      # knots which bound the finer  level of detail for C2 continuity
        self.level                  = level         # 0 is the coarsest
        level = THButilities.n2level[len(vertices)]
        self.level                  = level 
        self.interior_knots         = []        # ?
        self.parent                 = None      # chg to quad tree for surf
        self.children               = None        # chg to quad tree for surf
        self.active_indicies        = self.active_(bounds)#set(self.offset_reference.keys())  # portions of the spline space that have associated vertices
        self.active_prange          = bounds
        self._verbose               = False
        self.checksums_             = False
        self.tol                    = 1.e-5
        self.wgen                   = wm.FindP #stupid name
        self.FindP                  = wm.FindP
        self.FindQ                  = wm.FindQ
        self.Pget                   = THButilities.ProjectionGetter
        #if rm is None:
        #    rm = self.wgen(self.k-1,self.level+1)
        """# if you set this here, you must reset it if dyadically refining self.
        # not worrying about it for the moment
        # perhaps it is better to be lazy anyway
        #"""
        self.rm                     = rm
        self.rmi                    = rmi
        self.PQinverse              = None
        
        if interpolate:
            # Intrinsic Method Calls : Fully Initiazlize the B-Spline
            self.ivertices      = np.copy(self.vertices) #copy old vertices for plotting
            self.GlobalCurveInterp() #find new vertices and knot vector to interpolate the old vertices
        self.interpolate = interpolate

    def __call__(self, u):
        return self.eval(u)
        #return self.rthbeval(u)

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
    
    @property
    def level(self):
        return self._level
    @level.setter
    def level(self, l):
        self._level = l
        #
        # self.P is not used!
        #
        #assert(self.P is not None),'ERROR: THBspline curve projection (to fine) matrix should be given'
        #if self.P is None:
        #    self.P = wm.FindP(self.p,1)
        #
        ##self.P = wm.FindP(self.p,self.n-self.p) #NOV 2017 - do not do this!
        #self.Q = wm.FindQ(self.p,self.n-4,'min')
        return
    
    
    
    #there is no point with this until the curve already exists
    # because it is not designed to always know its parents
    # especially at startup!
    def compute_curve_THB(self):
        if self.verbose:
            print 'recomputing THB-bspline curve points'
        self.allCurveArray_THB()
        return
    
    def allCurveArray_THB(self, init=True):
        if init:
            self.children = None
            self.parent = None
            self.bounds = BoundsList([ia(0.,1.)])
        k=self.k
        #p=curve.p
        #U=curve.t
        P=self.project_vertices_to_this_level()
        """
            Return the array of all curve values
            for the standard parameterization"""
        self.allr[:,:,:]=0.0
        for d in range(self.dim): # loop over all spatial dimensions
            for j in range(k): #loop over all existing derivatives
                for i in range(self.n): # loop over all points
                    #implicit loop over all curve parameterization points:
                    self.allr[:,j,d]+=self.allbasis[:,j,i]*P[i,d]

        self.r = self.allr[:,0,:]
        self.dr = self.allr[:,1,:]
        
        return
#    
        
    def GlobalCurveInterp(self):
        """P&T algorithm page 369
            Global Curve Interpolation
            through n+1 points

            Inputs:
                n = self.n-1
                Q = self.vertices
                r = self.dim # defined on page 364
                p = self.p
                m = n+p+1
                U = self.t
            Output:
                P = the new vertices to replace Q
                in self.vertices
                U = the new knot vector
            
            Other definitions:
                uk_bar : parameter values assigned to Q_k 
                            - the parameter position 
                              at which the curve passes 
                              through the point to interpolated
                U      : the new knot vector
                
        """
        
        n = self.n-1
        Q = self.vertices
        r = self.dim
        p = self.p
        m = n+p+1
        U = self.t #comment out to JH mod form

        
        #Find the Chord Length: (9.5)
        d=0.
        store_sections = []
        for k in range(1,n+1):
            drule=0. 
            for kk in range(r):
                drule += (Q[k,kk]-Q[k-1,kk])**2
            stored = np.sqrt(drule)
            store_sections.append(stored)
        d = sum(store_sections)
            
        #Define the uk_bar: (9.6)
        uk_bar=np.zeros((n+1),float)
        uk_bar[n]=1.0
        for k in range(1,n):
            uk_bar[k]=uk_bar[k-1]+store_sections[k-1]/d
            
        #Compute the new Knot Vector U: (9.8)
        #U = np.zeros((m+1),float) #JH mod
        for j in range(1,n-p+1):
            U[j+p]=0.
            for i in range(j,j+p):
                U[j+p]+=uk_bar[i]
            U[j+p] = U[j+p]/p
        #U[m-p:m+2]=1.  #JH mod

        #initilize the A array (coefficient matrix):
        A = np.zeros((n+1,n+1),float)
        for i in range(n+1):
            span = self.FindSpan(uk_bar[i])
            self.BasisFuns(span,uk_bar[i],A[i,span-p:span+1])
        
        #invert the sytem of equations:
        P = np.linalg.solve(A,Q[:,:])

        #update the curve parameters
        self.vertices=P
        self.t=U
        return

    def active_(self, bounds=None):
        """dictionary of bounds index to start, stop
        inclusive list of all possibly active, 
        possibly refined, basis functions
        
        parameters
        ----------
        bounds = BoundsList (a list of ia bounds)
        """
        if bounds is None: 
            bounds = self.bounds
            #assert(len(bounds.bounds)==1)
            #bounds = bounds.bounds[0]
        actv_indicies = {}
        for i, el in enumerate(bounds):
            actv_indicies[i] = self.active(el)
        return actv_indicies
    
    def return_active_basis_and_indices(self,u):
        indices_list =  self.bases(self.active(u))
        return self.TLMBasisFuns(u)[indices_list], indices_list


    #    def active_basis(self, bounds=None,spans=None):
    #        """active_knots 
    #        is a misnomer.
    #        """
    #        return self.active_knots(self, bounds=None,spans=None)
    
    def get_active_basis_indices(self, bounds):
        if isinstance(bounds,float):
            return self.active_Bik_interval(u=bounds)
        else:
            bi = bounds[0]
            be = bounds[1]
            si = self.FindSpan(bi)
            se = self.FindSpan(be)
            return [si-self.p,se]
    
    def active_basis(self, bounds=None):
        return self.active_knots(bounds)
    def active_knots(self, bounds=None):#,spans=None):
        """ tells which BASIS (basis not Knots!!!) 
        and (indirectly control points)
        are active on this level of the hierarchy
        within supplied bounds
        
        mea maxima  - what a stupid 
        thing to totally mis-represent.
        
        using the basic
        find span
        [span-p:span]
        technique
        
        
        parameters
        ----------
            bounds
            spans
            
            
        notes
        ----------
        This is not the same as enclosed_basis_in_range
        -this gives any index of a basis which has inflence in 
        the bounds input here.
            
        same as curve.active, except this one assumes bounds is an interval
            -and can never be a scalar
            -use active_basis_at_u for that
        """
        return self.get_active_basis_indices(bounds)
            
    def active_basis_given_cv_index(self,i):
        """Given a basis func, aka C.V. index
        return the knots over which it is active
        
        ... wait, knots or vertices?
        what are you trying to do??
        
        this was wrong: 
            old return  self.t[i:i+self.p+3]
            
        Identical functions:
            -
        
        dev
        ----------
        try this:
            self.bases(self.active_knots(.4))
        vs this:
            self.active_knots_given_cv_index(self.FindSpan(.4))
            
            
            self.bases(self.active_knots(.4))
            Out[114]: [3, 4, 5, 6]
            
            self.active_knots_given_cv_index(self.FindSpan(.4))
            Out[115]: array([0.   , 0.125, 0.25 , 0.375])
            
        """
        self_bounds = self.get_range_of_basis(i)
        return self.enclosed_basis_in_range(self_bounds)

        
 
    
    def active_projective_indices_given_cv_index(self, i): 
        """given the ith basis index on this level
        return the basis indices at the next level
        which are encompassed by this basis function's 
        interval of activeness (interval where a 
        standard B-spline version is non-zero on this level) 
        """
        self_bounds = self.get_range_of_basis(i)
        #if self_bounds is None:
        #    self_bounds = ia(0.0,0.0)
        #    self_bounds.isempty = True
        if self.children is not None:
            ans = self.children.enclosed_basis_in_range(self_bounds)
            return ans 
        else:
            print 'WARNING: projecting bounds to a level that does not exist'
            return []
    
    def active_projective_bounds_given_cv_index(self, i): 
        """given the ith basis index on this level
        return the enclosure-bounds encompassed by the basis
        which replace it on the next level
        """
        indices = self.active_projective_indices_given_cv_index(i)
        if indices:
            a = self.children.get_range_of_basis(indices[0])
            b = self.children.get_range_of_basis(indices[-1])
            return  ia(a[0],b[0])
        else:
            self_bounds = ia(0.0,0.0)
            self_bounds.isempty = True
            return self_bounds
    

    def get_bounds_given_index(self,i):
        """
        return active bounds this level
        given basis index i
        
        TODO: cleanup.  This function
        is the same as:
        ----------
        * active_bounds_given_cv_index
        * get_range_of_basis
        """
        return self.get_range_of_basis(i)
    
    def active_bounds_given_cv_index(self,i):
        """
        return active bounds this level
        given basis index i
        
        TODO: cleanup.  This function
        is the same as:
        ----------
        * get_bounds_given_index
        * get_range_of_basis
        """
        bds = self.t[i:self.k+i+1] #another way to write it
        #note that python will exclude the k+i+1 'th index
        # due to the way that the vector[a:b] works
        # ie vector[a:b] returns [vector[a], vector[a+1]... vector[b-1]]
        return ia(bds[0],bds[-1])
    

    def active_basis_at_u(self, u):
        return self.active_Bik_interval(u=u)
    
    
    def given_ParentIndex_get_ChildIndex(self, parent_i):
        """
        current issue: later in the curve, 
        it starts to just give the end node..
        
        think about it some more 
            
            #child_check_bounds = parent.active_projective_bounds_given_cv_index(index)
            #child_selectable_basis = new_curve.active_knots(bounds = child_check_bounds)
        """
        assert(self.parent is not None),'Error, parent does not exist'
        if parent_i == 0:
            return 0
        
        parent = self.parent
        
        ga = parent.greville_abscissa(parent_i)
        child_real_basis = self.active_basis_at_u(ga)
        #child_real_basis[-1]+=1
        child_real_basis = range(*child_real_basis)
        cindex = child_real_basis[1]
        return cindex
    
    
    def given_Index_get_ChildIndex(self, i):
        """
        really an alias for 
        given_ParentIndex_get_ChildIndex
        
        
        """
        if i == 0:
            return 0
        child = self.children
        
        ga = self.greville_abscissa(i)
        child_real_basis = child.active_basis_at_u(ga)
        #child_real_basis[-1]+=1
        child_real_basis = range(*child_real_basis)
        cindex = child_real_basis[1]
        return cindex


    def get_range_of_basis(self, basis_index):
        """
        
        
        TODO: cleanup.  This function
        is the same as:
        ----------
        * get_bounds_given_index
        * active_bounds_given_cv_index
        
        
        dev
        ----------
        
        c0c.get_range_of_basis(0)
        >>> ia(0.0, 0.125)
            
        c0c.get_range_of_basis(10)
        >>> ia(0.875, 1.0)
        """
        a = self.t[basis_index]
        b = self.t[basis_index+self.k]
        return ia(a,b)

    
    def enclosed(self, bounds):
        return self.enclosed_basis_in_range(bounds)
    def enclosed_basis_in_range(self, bounds=None):
        """return the basis enclosed within this interval
        
        parameters
        ----------
        bounds:             a single 'ia' interval
        
        
        returns
        ----------
        contained_basis:    a list of indices
        
        note
        ----------
        Not the same as active_basis
        active_basis returns ALL basis which act on the interval
        
        This returen ONLY the basis totally enclosed by the interval
            
        dev
        ----------
        self.get_range_of_basis(2)
        >>> ia(0.0, 0.375)
        
        self.enclosed_basis_in_range(ia(.0,.375))
        >>> [0, 1, 2]
        
        self.enclosed_basis_in_range(ia(.0,.374))
        >>> [0, 1]
        
        rtbspline.active_knots : same as curve.active
        self.active(.5)
        >>> [4, 7]
        """
        if bounds is None: 
            bounds = self.bounds
            assert(len(bounds.bounds)==1)
            bounds = bounds.bounds[0]
        assert(isinstance(bounds, ia)), 'bounds is {}'.format(type(bounds) )
        #active_knots = self.active_knots(bounds) #this means active basis indices, not knots
        active_knots = self.active_basis(bounds) 
        active_ranges = []
        contained_basis = []
        indices = []
        #
        for basis_index in range(active_knots[0],active_knots[1]+1):
            indices.append(basis_index)
            active_ranges.append(
                    self.get_range_of_basis(basis_index))

        for basis_index, active_range in zip(indices,active_ranges):
            if active_range in bounds:
                contained_basis.append(basis_index)
        return contained_basis
    
    def enclosed_basis_curve(self):
        """Get the set of all
        possibly active basis on a given curve
        
        fixed to actually return what it says now!
        """
        eb = []
        for bds in self.bounds:
            eb += self.enclosed_basis_in_range(bds)
        list(set(eb)).sort()
        return eb



    def dyadic_refinement_detailed(self):
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
        rspline.level = self.level+1
        rspline.base_curve = False
        locations = []
        for i in range(Nseg):
            #nknots[i] = (t[k+i]+t[k+i+1])/2.
            iknot = (t[k+i-1]+t[k+i])/2.
            tnew, vertices, loc = ko.knot_insertion_hierarchical(rspline.k,
                                               rspline.vertices,
                                               rspline.t,
                                               iknot)
            locations.append((loc-p,loc))
            #print 'vertices {}:{} will change'.format(loc-p,loc)
            #rspline.t           = tnew
            rspline.vertices    = vertices
            rspline.t           = tnew
        self.children = rspline
        #self.children.append(rspline)
        rspline.parent = self
        return rspline, locations



    def dyadic_refinement(self,
                          make_projectors=False):
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
        rspline.level = self.level+1
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
        #if make_projectors:
        #    rspline.rm = rspline.FindP(rspline.p,rspline.level+1)
        self.children = rspline
        #self.children.append(rspline)
        rspline.parent = self
        return rspline


    def isactive(self, pt):
        """ 'Fixed' March 9, 2018
        active_bounds_this_pt = self.bounds.contains(pt)
        """
        #if pt in self.active_intervals:
        active_bounds_this_pt = self.bounds.contains(pt)
        if active_bounds_this_pt:
            return True
        else:
            return False

    def get_active_curve(self, pt):
        if self.isactive(pt):
            return self
        else:
            self.children.get_active_curve(pt)

    def print_knots(self):
        for i, el in enumerate(self.t):
            print i, el

    def print_properties(self):
        print '# knots = {}'.format(len(self.t))
        print '# vertices = {}'.format(self.n)

    def bases(self, ends):
        """return a range
        from end[0] to end[1], inclusive
        
        usage:  pass in 
        ends = self.bases(self.active(u))
        returns the bases active 
        at u
        """
        #range(active_parent_basis[0],active_parent_basis[1]+1)
        return range(ends[0],ends[1]+1)
    
    def clip_basis_outside_refinement(self, basis, interval):
        save = np.zeros_like(basis)
        save[interval] = copy.copy(basis[interval])
        basis = basis*0.
        basis[interval] = save[interval]
        return basis
            
    
    def get_one_level(self, u):
        """
        
        about the methods used
        --------------------
        self.bounds.contains(u) : returns the bounds box containing u
        self.bases() : 
        self.active() : elemental Bspline property: 
                            uses findspan to return the potentially infuential
                            range (non-zero) basis funcs
            NOTE: this actually lives way back in the curve.py class file!
        """
        active_parent_bounds_index = self.bounds.contains(u)
        len_apbi = len(active_parent_bounds_index)
        assert(len(active_parent_bounds_index)<2),'error: parent bounds active in two spans at once'
        if len_apbi==1:
            #it would be weird if bounds overlap.  Take the 1st
            parent_bounds = self.bounds[active_parent_bounds_index[0]]
        else:
            parent_bounds = ia(0.,0.)#self.bounds[0]#dummy_bounds#
        
        active_parents = self.enclosed_basis_in_range(
                                        parent_bounds
                                        ) #basis within Omega at this level
        apb = self.bases(self.active(u)) #active basis at u at this level
        return active_parents, apb
    
    def get_all_one_level(self):
        """Get the set of all
        possibly active basis on a given curve
        
        this is actually the same as self.enclosed_basis_curve
        
        this function does not care what
        particular basis may be non zero
        at a particular parametric location.
        
        use when projecting the whole
        truncated space
        """
        tot_within_bounds = []
        for interval_bound in self.bounds:
            active = self.enclosed_basis_in_range(
                                            interval_bound)
            tot_within_bounds += active
            
        list(set(tot_within_bounds)).sort()
        return  tot_within_bounds
    
    
    def truncate_this_level(self, oBasis, iBasis, apb):
        if self.rm is None:
            try:
                self.rm = self.Pget(self.level)
            except:
                self.rm = self.FindP(self.p,self.level+1)
            #self.rm = self.wgen(self.p,self.level+1)#
        oBasis[apb] = np.matmul(self.rm.T,iBasis)[apb]
        #oBasis[apb] = oBasis[apb] - np.matmul(self.rm.T,iBasis)[apb]
        return oBasis
    
    
    def do_recursion(self, u, basis_levels):
        dummy_bounds = ia(0.,0.)
        
        oBasis = basis_levels[-1]
        
        active_parents, apb = self.parent.get_one_level(u)
        oBasis = self.clip_basis_outside_refinement(oBasis, 
                                                        active_parents)
        
        iBasis = self.TLMBasisFuns(u)
        active_parents, apb = self.get_one_level(u)
        iBasis = self.clip_basis_outside_refinement(iBasis, 
                                                active_parents)
        
        child = self
        
        active_child_bounds_index = child.bounds.contains(u)
        len_acbi = len(active_child_bounds_index)
        assert(len(active_child_bounds_index)<2),'error: child bounds active in two spans at once'
        if len_acbi==1:
            child_bounds = child.bounds[active_child_bounds_index[0]]
            inner = child.enclosed_basis_in_range(child_bounds)
            basis_levels.append(iBasis)
            basis_levels = child.do_recursion(u, basis_levels)#.TLMBasisFuns(u)
            outer = self.parent.enclosed_basis_in_range(child_bounds) #to be eliminated
            
        else:
            child_bounds = dummy_bounds
            inner = []
            basis_levels.append(iBasis)
            basis_levels = child.do_recursion(u, basis_levels)##.TLMBasisFuns(u)
            outer = [] #none to be eliminated
            # if none are encl
        #iBasis = basis_levels[-1]
        
        if self.rm is None:
            try:
                self.rm = self.Pget(self.level)
            except:
                self.rm = self.FindP(self.p,self.level+1)
            #self.rm = self.wgen(self.p,self.level+1)#
        
        #iBasis[inner] = 0.
        oBasis[apb] = np.matmul(self.rm.T,iBasis)[apb]
        oBasis[outer] = 0.
        
        if self.children is not None:
            return self.children.do_recursion(u, basis_levels)
        else:
            return basis_levels
        
        

    def rthb_basis(self, u):
        """
        This one does it!!!
        Done:
            -build level from
            sum of all levels below it
            -not just one level below it!

        b = c0c.rthb_basis(.7)
        sum(b[0]) + sum(b[1]) + sum(b[2])

        notes:  this method always sounds to the bottom
        of the hierarchy and accrues all the way up!
        
        Issue for THBsurfaces:
            trying to roll through here is problamatic when 
            parents and children don't know each other
            directly!
        solution:
            always just use the -master_curves- 
            (+ local bounding boxes)
            to construct the basis functions.
        one pot. issue:
            ok to parent up the hierarchy
            but try to child down it
            and you will stay with mastercurves
            -this may be ok for basis
            as long as bounds are correct!
        """
        dummy_bounds = ia(0.,0.)
        basis_levels, inactive_basis_levels, lcurve = self.recursive_to_bottom(u)
        level_mats = []
        
        #active_mats = []
        #TODO in THBsurface: make lcurve.parent make sense.
        #note that it is only needed for basis, -whew-
        while lcurve.parent is not None: #TODO: make this work for surface curves!
            #child = lcurve
            lcurve = lcurve.parent
            iBasis =  basis_levels[-1]
            active_parents, apb = lcurve.get_one_level(u)
            oBasis = lcurve.TLMBasisFuns(u)
            child = lcurve.children

            active_child_bounds_index = child.bounds.contains(u)
            len_acbi = len(active_child_bounds_index)
            assert(len(active_child_bounds_index)<2),'error: child bounds active in two spans at once'
            if len_acbi==1:
                child_bounds = child.bounds[active_child_bounds_index[0]]
                inner = child.enclosed_basis_in_range(child_bounds)
                outer = child.parent.enclosed_basis_in_range(child_bounds) #to be eliminated

            else:
                child_bounds = dummy_bounds
                inner = []
                outer = [] #none to be eliminated
                # if none are encl
                #apb=[]


            if lcurve.rm is None:
                lcurve.rm  = lcurve.wgen(lcurve.p,lcurve.level+1)#
            
            #playing with wavelet subdivision
            """
            if lcurve.rmi is None:
                lcurve.rmi = lcurve.FindQ(lcurve.p,lcurve.level+1,'min')
                Qshape = np.shape(lcurve.rmi)
                Pshape = np.shape(lcurve.rm)
                ncol = Qshape[1] + Pshape[1] 
                #assert(ncol == Qshape[1]+Pshape[1]),'ERROR transform matrix must be SQUARE"
                PQM1 = np.zeros((Pshape[0],ncol))
                PQM1[:,:Pshape[1]] = lcurve.rm
                PQM1[:,Pshape[1]:] = lcurve.rmi
                #same as the slightly more intuitive:
                #PQM1[:,-Qshape[1]:] = lcurve.rmi
                lcurve.PQM1 = PQM1
                lcurve.AB1 = np.linalg.inv(PQM1)
                lcurve.wavelet_phi = np.zeros((lcurve.n+Qshape[1],1))
            #"""
                
                

            level_mats.append(lcurve.rm)
            #oBasis[apb] = oBasis[apb] -np.matmul(lcurve.rm[inner].T,iBasis[inner])[apb]
            #oBasis[apb] = np.matmul(lcurve.rm[inner].T,iBasis[inner])[apb]
#            cl = len(basis_levels)
#            matrix_reducer = np.matmul(level_mats[-1].T,iBasis) #usesactive only
#            for i in range(cl-1):
#                cBasis = basis_levels[i]
#                local_reducer = np.matmul(level_mats[i].T,cBasis)
#                for j in range(i+1,cl):
#                    #transBasis = basis_levels[j]
#                    local_reducer = np.matmul(level_mats[j].T,local_reducer)
#                matrix_reducer = matrix_reducer + local_reducer
#            oBasis[apb] = oBasis[apb] - matrix_reducer[apb] #mustsubtract since we used active
            #oBasis = matrix_reducer

            iBasis = inactive_basis_levels[-1]
            oBasis = np.matmul(level_mats[-1].T,iBasis)
            
            
            #playing with wavelet subdivision
            """
            lcurve.wavelet_phi[:] = 0.
            lcurve.wavelet_phi[0:lcurve.n] = oBasis
            c1_cktwoscale = np.matmul(iBasis.T,lcurve.PQM1)
            lcurve.wavelet_phi[lcurve.n:] = c1_cktwoscale.T[lcurve.n:]
            cb2_from_c1b = np.matmul(lcurve.wavelet_phi.T,lcurve.AB1)
            #"""
            """
            print 'level 2 basis funs from level 1 basis funs:'
            print cb2_from_c1b.T
            print 'difference between computed and actual:'
            print cb2_from_c1b.T - iBasis
            #"""
            #done

            #oBasis[apb] = oBasis[apb] - np.matmul(matrix_reducer.T,iBasis)[apb]
            #oBasis[apb] = oBasis[apb] -np.matmul(lcurve.rm[inner].T,iBasis[inner])[apb]
            #oBasis[apb] = oBasis[apb] -np.matmul(lcurve.rm[inner].T,iBasis[inner])[apb]
            oBasis[outer] = 0.

            active_basis = lcurve.clip_basis_outside_refinement(oBasis,
                                                    active_parents)

            inactive_basis = oBasis - active_basis
            basis_levels.append(active_basis)
            inactive_basis_levels.append(inactive_basis)
        return basis_levels
    
    def sum_basis(self, u):
        basis_list = self.rthb_basis(u)
        ans = 0.
        for b in basis_list:
            ans = ans + sum(b)
        return ans
    
    def rthbeval_old(self, u):
        """
        The correct ans comes from this guy!
        using the rthb_basis
        recursive (in form if not in code) thbasis
        
        notes:  this method always sounds to the bottom
        of the hierarchy and accrues all the way up!
        
        for the surfaces version see 
        def eval1thb
        """
        b = self.rthb_basis(u)
        ans=0.
        lcurve = self.get_bottom_curve()
        #while lcurve.parent is not None:
        for lel in b:
            
            verts = lcurve.vertices
            ans = ans + np.dot(lel.T,verts)
            lcurve = lcurve.parent
        return ans
    
    
    def build_restriction_operator(self):
        """Wavelets for computer graphics, 1995!
        """
        if self.rm is None: #control vertex projection (prolongation)
            try:
                self.rm = self.Pget(self.level)
            except:
                self.rm = self.FindP(self.p,self.level+1)
        if self.rmi is None: # 
            self.rmi = self.FindQ(self.p,self.level+1)
        Qshape = np.shape(self.rmi)
        Pshape = np.shape(self.rm)
        ncol = Qshape[1] + Pshape[1] 
        PQM1 = np.zeros((Pshape[0],ncol))
        PQM1[:,:Pshape[1]] = self.rm
        #print 'broadcasting'
        #print 'shape rmi = ',Qshape
        #print 'shape rm = ',Pshape
        PQM1[:,Pshape[1]:] = self.rmi
        self.PQM1 = PQM1
        self.PQinverse = np.linalg.inv(PQM1)
        #todo, get scipy linalg.LU working!
        A = self.PQinverse[:self.n,]    #scaling functions (cntrl pts)
        B = self.PQinverse[self.n:]     #wavelets
        return  A,B
        
    
    def restrict(self, vector=None,
                thiscurve=None,
                use_active=True):
        """ 1 Level vector restriction
        Project data from a Fine level-> Coarse Level
        
        parameters
        ----------
            vector      : vector of information about
                          the current control point vector at this level 
                          (e.g. the points themselves, or gradient data)
                          
            input_curve : the curve at the same level in the 
                          hierarchy as the vector of data
            
        assumptions:
        ----------
            vector is at the same level as self.rmi
            
        details
        ----------
            use_active
            
            -truncate this level
                
            -the return value does not know about 
            what is inactive at this level,
                          
        dev
        ----------
            c0c.vertices - c1.restrict()
            Out[28]: 
            array([[ 3.46944695e-17,  1.77635684e-15],
                   [-2.77555756e-17, -1.77635684e-15],
                   [ 3.33066907e-16,  3.55271368e-15],
                   [-2.22044605e-16, -2.66453526e-15],
                   [ 2.22044605e-16,  8.88178420e-16],
                   [ 1.33226763e-15,  4.44089210e-16],
                   [-8.88178420e-16,  2.22044605e-16],
                   [ 0.00000000e+00,  3.88578059e-16],
                   [ 1.77635684e-15, -6.93889390e-17],
                   [-3.55271368e-15,  5.55111512e-17],
                   [ 0.00000000e+00, -7.97972799e-17]])
    
            vb = np.matmul(c1.rmi,c1.vertices)
        """
        #----------------------------------------------------------------------
        if vector is None:
            vector = self.vertices
        #----------------------------------------------------------------------
        if thiscurve is None:
            this=self
            #this = thiscurve.get_finest_curve()
        else:
            this = thiscurve
        #----------------------------------------------------------------------
        P = this.vertices.copy()
        #----------------------------------------------------------------------
        A,B = this.parent.build_restriction_operator()
        #
        #**********************************************************************
        if not use_active:
            #if this.parent.PQinverse is None:
            # B: wavelet coefficients, aka detail coefficients
            return np.matmul(A,P)#np.matmul(vector.T,this.parent.PQinverse.T)
        else:
            #print 'restricting to coarse'
            active, inactive = this.get_all_active_inactice_indices()
            #P[inactive,:] = 0.
            #project to coarse level
            P =  np.matmul(A,P)
            current_active = this.parent.enclosed_basis_curve()
            #set priors to zero? - this was done above?
            P[current_active] = 0.
            
            active, inactive = this.parent.get_all_active_inactice_indices()
            lP = this.parent.vertices.copy()
            lP[inactive] = 0.
            
            P += lP
            return P
    
    def prolong(self, vector=None,
                thiscurve=None,
                use_active=True):
        """1 Level vector prolongation
        Project data from a coarse level -> fine level
        
        parameters
        ----------
            vector      : vector of information about
                          the current control point vector at this level 
                          (e.g. the points themselves, or gradient data)
                          
            input_curve : the curve at the same level in the 
                          hierarchy as the vector of data
                          
        assumptions:
        ----------
            vector is at the same level as self.rm 
            
        details
        ----------
            use_active:
            -truncate this level
            -the return value does not know about 
            what is inactive at this level
            -vector (incoming)
            will be prolonged THB stlye together
            with this levels vertices
            (maybe we could make a basis option too?)
        """
        #----------------------------------------------------------------------
        if vector is None:
            vector = self.vertices
        #----------------------------------------------------------------------
        if thiscurve is None:
            this=self
            #this = self.get_coarsest_curve
        else:
            this = thiscurve
        #----------------------------------------------------------------------
        if this.rm is None:
            this.rm = this.FindP(this.p,this.level+1)
        #**********************************************************************
        if not use_active:
            return  np.matmul(this.rm,vector)
        else:
            if this.rm is None:
                this.rm = this.FindP(this.p,this.level+1)
                
            #Step1 original level, #find active #zero out of bounds and #project up::
            P = vector.copy()  
            active, inactive = this.get_all_active_inactice_indices()
            #P[inactive,:] = 0. 
            dv = np.matmul(this.rm,P) #project to prolonged
            
            # prolong level, #find active #zero out of bounds 
            if this.children is not None:
                active, inactive = this.children.get_all_active_inactice_indices()
                lP = this.children.vertices.copy()
                lP[inactive,:] = 0.
                current_active = this.children.enclosed_basis_curve()
                
                #truncate the interior of the prolonged
                # vertices from the original level
                # prior level
                dv[current_active] = 0.
                
                P = dv+lP
            
            
            return  P
            
    def prolong_vertices_to_this_level(self,
                                       other_curve_at_some_level=None,
                                       vector=None):
        """
        equivalent to project_vertices_all_active
        -but pushes all vertices to the level of the original caller
        -conducting THB truncation along the way
        """
        if other_curve_at_some_level is None:
            other_curve_at_some_level = self.get_coarsest_curve()
            vector = other_curve_at_some_level.vertices
        
        if other_curve_at_some_level is self:
            return vector
        
        
        vector = other_curve_at_some_level.prolong(vector=vector)
        child = other_curve_at_some_level.children
        if child is not self:
            return self.prolong_vertices_to_this_level(child,
                                                          vector = vector)
            
        return vector
    
    def restrict_vertices_to_this_level(self, 
                                        other_curve_at_some_level=None,
                                        vector=None):
        if other_curve_at_some_level is None:
            other_curve_at_some_level = self.get_finest_curve()
            vector = other_curve_at_some_level.vertices
            
            #active, inactive = other_curve_at_some_level.get_all_active_inactice_indices()
            #vector[inactive] = 0.
            
        if other_curve_at_some_level is self:
            return vector
        vector = other_curve_at_some_level.restrict(vector=vector)
        parent = other_curve_at_some_level.parent
        if parent is not self:
            return self.restrict_vertices_to_this_level(parent,
                                                        vector=vector)
        return vector
    
    def project_vertices_to_this_level(self):
        vr = self.restrict_vertices_to_this_level()
        vp = self.prolong_vertices_to_this_level()
        active, inactive = self.get_all_active_inactice_indices()
        current_active = self.parent.enclosed_basis_curve()
        vr[current_active] = 0.
        vp[inactive] = 0.
        return vr+vp
    
    
    def project_vertices_all_active(self):
        """take the easiest to understand THB basis
        for curves and make it work independent of u
        
        project all active __
        to the highest level of detail
        evaluate the surace there using standard 
        B-spline surface represetation
        
        inverse of the other evaluation methods seen in this class
        
        this one evaluates from coarse to fine
        project then truncate
        
        
        The point of this:
            ----------
            -seperate basis evaluation from curve property evaluation
            -carry out all curve property evaluation on the finest level
            -make possible that standard path for doing
            variational design of THB surfaces.
        """
        #get the coarsest level:
        this = self.get_coarsest_curve()
        dummy_bounds = ia(0.,0.)
        count = 0
        P = this.vertices.copy() 
        while this is not None:
            
            if this.parent is None:
                P = this.vertices.copy()  
                active, inactive = this.get_all_active_inactice_indices()
                #P[inactive,:] = 0. #technically correct to truncate this level now!
                # it is correct to 0 out the inactive basis
                # in general, but since this is the highest level anyway
                # this does nothing. -- assuming nested spaces!
                #
                # Why does is the highest level not affected
                # by this you might ask?
                # A:  Because the highest level has no inactive 
                # indices according to it's own bounds!
                #
                # also:
                # P must be "inner-truncated" after it is prolongated!
                # that is the "real THB" truncation mechanism.
                # because the fine level knows how that should go.
                # not the current level!
                # it happens in code below.
            else:
                lP = this.vertices.copy()  
                
                #make the prolongation/projection-to-fine operator
                # if not found 
                if this.parent.rm is None:
                    this.parent.rm  = this.wgen(this.parent.p,
                                         this.parent.level+1)
                
                #active and inactive at 'this' level
                active, inactive = this.get_all_active_inactice_indices()
                #
                #*********************************
                #  deactivate inactive vertices 
                # (truncation of simple non-active vertices at this level)
                lP[inactive,:] = 0.
                #
                #*********************************
                # prepare to project the prior coarse level up to 
                # the current level
                dv = np.zeros_like(this.vertices)
                #
                #*********************************
                # prolongation:
                # project prior (lower) level up to current:
                dv[:] = np.matmul(this.parent.rm, P[:]) 
                #
                #*********************************
                #  deactivate inactive vertices (truncation)
                # of level before truncation
                #*********************************
                # bunch of deleted code in here!
                #*********************************
                # much simpler not to check exact location
                
                #this is actually the same as 'active', above
                #current_active = this.enclosed_basis_curve()
                
                #
                #*********************************
                # THB Truncation of the vector 
                # projected here from 
                # the prior level:
                # all those active at this level
                # entailed deactivating at the prior level
                #
                #dv[current_active,:] = 0.
                dv[active,:] = 0.
                #
                #*********************************
                # 
                lP += dv.copy()
                #
                P = lP.copy() #update P to this level, after projection and truncation
                    
            this = this.children
            count +=1
        return P
    
    
    def project_vertices_sensible(self,u):
        """This one is clear.
        but still reliant on u
        
        sensible - Project to make the easiest to understand THB basis
        for curves
        
        project all active ___
        to the highest level of detail
        evaluate the surace there using standard 
        B-spline surface represetation
        
        inverse of the other evaluation methods seen in this class
        
        this one evaluates from coarse to fine
        project then truncate
        
        
        The point of this:
            ----------
            -seperate basis evaluation from curve property evaluation
            -carry out all curve property evaluation on the finest level
            -make possible that standard path for doing
            variational design of THB surfaces.
        """
        #get the coarsest level:
        this = self.get_coarsest_curve()
        dummy_bounds = ia(0.,0.)
        count = 0
        P = this.vertices.copy() 
        while this is not None:
            
            if this.parent is None:
                P = this.vertices.copy()  
                active, inactive = this.get_active_inactice_indices(u)
                P[inactive,:] = 0. #technically correct to truncate this level now!
                # it is correct to 0 out the inactive basis
                # in general, but since this is the highest level anyway
                # this does nothing. -- assuming nested spaces!
                #
                # also:
                # P must be "inner-truncated" after it is prolongated!
                # that is the "real THB" truncation mechanism.
                # because the fine level knows how that should go.
                # not the current level!
                # it happens in code below.
            else:
                lP = this.vertices.copy()  
                
                #make the prolongation/projection-to-fine operator
                # if not found 
                if this.parent.rm is None:
                    this.parent.rm  = this.wgen(this.parent.p,
                                         this.parent.level+1)
                
                #active and inactive at 'this' level
                active, inactive = this.get_active_inactice_indices(u)
                #
                #*********************************
                #  deactivate inactive vertices 
                # (truncation of simple non-active vertices at this level)
                lP[inactive,:] = 0.
                #
                #*********************************
                # prepare to project the prior coarse level up to 
                # the current level
                dv = np.zeros_like(this.vertices)
                #
                #*********************************
                # prolongation:
                # project prior (lower) level up to current:
                dv[:] = np.matmul(this.parent.rm, P[:]) 
                #
                #*********************************
                #  deactivate inactive vertices (truncation)
                # of level before
                active_bounds_index = this.bounds.contains(u)
                len_acbi = len(active_bounds_index)
                if len_acbi==1:
                    bounds = this.bounds[
                            active_bounds_index[0]]
                    #current level basis enclosed in these bounds:
                    current_active = this.enclosed_basis_in_range(bounds)
                else:
                    mssg = 'ERROR: active bounds in two places at once'
                    assert(len_acbi<2),mssg
                    bounds = dummy_bounds
                    current_active = [] 
                
                #
                #*********************************
                # THB Truncation of the prior level:
                # 'current_active' at this level
                # entailed deactivating at the prior level
                #
                dv[current_active,:] = 0.
                #
                #*********************************
                # 
                lP += dv.copy()
                #
                P = lP.copy() #update P to this level, after projection and truncation
                    
            this = this.children
            count +=1
        return P
    
    
    def eval_new2(self,u):
        return self.eval(u)
    
    def eval(self,u):
        """Use the best documented
        version of vertex projection 
        to evaluate with the basis on 
        the highest level.
        
        dev
        ----------
        self = c0c
        for u in np.linspace(0.,1.,30):
            print self(u)-self.eval_new(u)
            print self.CurvePoint(u)-self.eval_new(u)
            print self.CurvePoint(u)-self(u)
            print '--- --- --- ---'
        """
        this = self.get_bottom_curve()
        #
        #P = self.project_vertices_sensible(u) #better
        #
        P = self.project_vertices_all_active() #much simpler!
        # stationary!
        #
        basis = this.TLMBasisFuns(u)
        return np.dot(basis.T,P)
    
    
    
    def project_vertices(self, u):
        """project all active ___
        to the highest level of detail
        evaluate the surface there using standard 
        B-spline surface represetation
        
        inverse of the other evaluation methods seen in this class
        
        this one evaluates from coarse to fine
        project then truncate
        
        
        The point of this:
            ----------
            -seperate basis evaluation from curve property evaluation
            -carry out all curve property evaluation on the finest level
            -make possible that standard path for doing
            variational design of THB surfaces.
        """
        #get the coarsest level:
        this = self.get_coarsest_curve()
        dummy_bounds = ia(0.,0.)
        count = 0
        P = this.vertices.copy() 
        while this is not None:
            
            if this.parent is None:
                P = this.vertices.copy()  
            else:
                lP = this.vertices.copy()  
                
                #make the prolongation/projection-to-fine operator
                # if not found 
                if this.parent.rm is None:
                    this.parent.rm  = this.wgen(this.parent.p,
                                         this.parent.level+1)
                
                #active and inactive at 'this' level
                active, inactive = this.get_active_inactice_indices(u)
                #
                #*********************************
                #  deactivate inactive vertices (truncation)
                lP[inactive,:] = 0.
                #
                #*********************************
                # prepare to project to the next finer level:
                final_dv = np.zeros_like(this.vertices)
                #
                #*********************************
                # prolongation:
                # project prior (lower) level up to current:
                final_dv[:] = np.matmul(this.parent.rm, P[:]) 
                #
                #*********************************
                #  deactivate inactive vertices (truncation)
                # of level before
                active_bounds_index = this.bounds.contains(u)
                len_acbi = len(active_bounds_index)
                if len_acbi==1:
                    bounds = this.bounds[
                            active_bounds_index[0]]
                    #coarse level basis enclosed in bounds:
                    outer = this.enclosed_basis_in_range(bounds)
                else:
                    mssg = 'ERROR: active bounds in two places at once'
                    assert(len_acbi<2),mssg
                    bounds = dummy_bounds
                    inner = []
                    outer = [] 
                    
                final_dv[outer,:] = 0.
                #
                #*********************************
                # 
                lP += final_dv.copy()
                #
                P = lP.copy() #update P to this level, after projection and truncation
                    
            this = this.children
            count +=1
        return P


    def eval_new(self,u):
        """
        
        dev
        ----------
        self = c0c
        for u in np.linspace(0.,1.,30):
            print self(u)-self.eval_new(u)
            print self.CurvePoint(u)-self.eval_new(u)
            print self.CurvePoint(u)-self(u)
            print '--- --- --- ---'
        """
        this = self.get_bottom_curve()
        P = self.project_vertices(u)
        basis = this.TLMBasisFuns(u)
        return np.dot(basis.T,P)
    
    
    
    def rthbeval(self, u):
        """
        The correct ans comes from this guy!
        using the rthb_basis
        recursive (in form if not in code) thbasis
        
        notes:  
        ----------
            this method always sounds to the bottom
            of the hierarchy and accrues all the way up!
            
            for the surfaces version see 
            def eval1thb
            
        about the methods used
        --------------------
            get_active_inactice_indices: return a list of active 
                                        and a list of inactive basis functions
                                        at u
            wgen : generate the coarse to fine projection matrices
                    which go from the current level to one dyadic level finner
                    effected here: dv[i] = np.matmul(this_curve.rm,dv[i])
            
            
            
        
        """
        b = self.rthb_basis(u)
        ans=0.
        lcurve = self.get_bottom_curve() #finest mesh!
        basis = lcurve.TLMBasisFuns(u)
        #c1_cktwoscale = np.matmul(lcurve.PQM1.T,basis) 
        
        #while lcurve.parent is not None:
        dv = {}
        
        for i,lel in enumerate(b):
            #dv[i] = np.multiply(lcurve.vertices, lel)
            dv[i] = lcurve.vertices.copy() #np.multiply(lcurve.vertices, lel)
            active, inactive = lcurve.get_active_inactice_indices(u)
            dv[i][inactive] = 0. #at a given level, discard inactive on that level
            if i ==0:
                pass
            else:
                #target_basis = lcurve.child.TLMbasis(u)
                #wishful_basis = lcurve.TLMBasisFuns(u)
                #this_basis = lel
                
                
                #c1_cktwoscale = np.matmul(lcurve.PQM1.T,target_basis) 
                #ck_twoscale = np.matmul(lcurve.PQM1.T,target_basis) 
                
                #c1_cktwoscale[0:lcurve.n] = this_basis
                #ck_twoscale[0:lcurve.n] = wishful_basis
                
                #this_scale = np.matmul(lcurve.AB1.T,c1_cktwoscale)  #c1_cktwoscale.copy()
                #ck_scale = np.matmul(lcurve.AB1.T,ck_twoscale) 
                
                this_curve = lcurve
                #active, inactive = this_curve.get_active_inactice_indices(u)
                #dv[i][inactive] = 0. #backwards from how they say it, 
                #because we are on the level of the basis dv[i] itself
                
                while this_curve.children is not None: #scale lower levels up to this level
                    
                    if this_curve.children.rm is None:
                        this_curve.children.rm  = this_curve.wgen(lcurve.children.p,
                                                                  lcurve.children.level+1)#
                    #this_scale = np.matmul(this_curve.parent.PQM1.T,basis) 
                    dv[i] = np.matmul(this_curve.rm,dv[i]) #coarse to fine, one level 
                    active, inactive = this_curve.children.get_active_inactice_indices(u)
                    dv[i][active] = 0. #after projecting to fine, discard active, 
                    #dv[i][inactive] = 0.
                    this_curve = this_curve.children
            #verts = lcurve.vertices
            #ans = ans + np.dot(lel.T,verts)
            
            lcurve = lcurve.parent
        
        #for key in dv:
        #    ans += np.dot(b[key].T,dv[key])
        ans = 0.
        P = 0.
        for key in dv:
            P += dv[key]
#        for i,lel in enumerate(b):
#            #ans += dv[i]
#            
#            #ans += sum(dv[i]) #use this one to complete dot products if using basis above
#            
#            ans += np.dot(basis.T,dv[i]) #use this one if only using vertices above

        ans = np.dot(basis.T,P)
        return ans
        

    def get_finest_curve(self):
        """get finest level curve
        more recognizable alias
        """
        return self.get_bottom_curve()
    
    def get_bottom_curve(self):
        """get finest level
        """
        if self.children is None:
            return self
        else:
            return self.children.get_bottom_curve()
    
    def get_coarsest_curve(self):
        """get coarsest level
        """
        if self.parent is None:
            return self
        else:
            return self.parent.get_coarsest_curve()
    
    #"""
    def recursive_to_bottom_old(self, u, basis_levels=None):
        """u=.7
        """
        if basis_levels is None:
            basis_levels = []
        if self.children is None:
            basis_levels.append(self.bottom_level(u))
            return basis_levels, self
        else:
            return self.children.recursive_to_bottom(u, basis_levels)
    #"""
    def recursive_to_bottom(self, 
                            u, 
                            basis_levels=None,
                            inactive_bl=None):
        """u=.7
        """
        if basis_levels is None:
            basis_levels = []
            inactive_bl = []
        if self.children is None:
            activeBasis, inactive_Basis = self.bottom_level(u)
            basis_levels.append(activeBasis)
            inactive_bl.append(inactive_Basis)
            return basis_levels, inactive_bl, self
        else:
            return self.children.recursive_to_bottom(u,
                                                     basis_levels,
                                                     inactive_bl)
    
    #"""
    def bottom_level_old(self, u):
        #assert(self.children is None)
        oBasis = self.TLMBasisFuns(u)
        active_parents, apb = self.get_one_level(u)
        oBasis = self.clip_basis_outside_refinement(oBasis, 
                                                active_parents)
        return oBasis
    #"""
        
    def bottom_level(self, u):
        #assert(self.children is None)
        oBasis = self.TLMBasisFuns(u)
        active_parents, apb = self.get_one_level(u)
        activeBasis = self.clip_basis_outside_refinement(oBasis,
                                                active_parents)
        inactive_Basis = oBasis - activeBasis
        return activeBasis, inactive_Basis
            
    
    def get_active_inactice_indices(self, u):
        total_indices = range(self.n)
        active, non_zero = self.get_one_level(u)
        ti = set(total_indices)
        ai = set(active)
        ii = ti - ai
        return list(ai), list(ii)
        
    def get_all_active_inactice_indices(self):
        total = range(self.n)
        total_in_space = self.get_all_one_level() #same as self.enclosed_basis_curve()
        ti = set(total_in_space)
        deli = list(set(total) - ti) #delete_i'th ones (set to zero in projector)
        return total_in_space, deli
        
        
    
    def truncate_basis2(self, u):
        """
        TODO: double check logic when 
        u is outside bounds.
        """
        dummy_bounds = ia(0.,0.)
        active_parent_bounds_index = self.bounds.contains(u)
        len_apbi = len(active_parent_bounds_index)
        assert(len(active_parent_bounds_index)<2),'error: parent bounds active in two spans at once'
        if len_apbi==1:
            parent_bounds = self.bounds[active_parent_bounds_index[0]]
        else:
            parent_bounds = self.bounds[0]#dummy_bounds#
        
        active_parents = self.enclosed_basis_in_range(
                                        parent_bounds
                                        ) #generally non-zero
        """
        active parents are self indices 
        where self basis functions are 
        within bounds.
        """
        apb = self.bases(self.active(u)) #actually active
        """
        -active
        Returns the indices of the 
        active non-zero basis functions
        this is the index range of the control vertices
        which have influence here
        
        -bases makes a list of them.
        """
        
        
        oBasis = self.TLMBasisFuns(u)
        
        if self.children is not None:
            child = self.children
            
            active_child_bounds_index = child.bounds.contains(u)
            len_acbi = len(active_child_bounds_index)
            assert(len(active_child_bounds_index)<2),'error: child bounds active in two spans at once'
            if len_acbi==1:
                child_bounds = child.bounds[active_child_bounds_index[0]]

                inner = child.enclosed_basis_in_range(child_bounds)
                iBasis = child.TLMBasisFuns(u)
                ssf = 1.
                #iBasis, ssf = child.truncate_basis2(u)#.TLMBasisFuns(u)
                outer = self.enclosed_basis_in_range(child_bounds) #to be eliminated
                
            else:
                child_bounds = dummy_bounds
                inner = []
                iBasis = child.TLMBasisFuns(u)
                ssf = 1.
                #iBasis, ssf = child.truncate_basis2(u)##.TLMBasisFuns(u)
                outer = [] #to be eliminated
                # if none are enclosed, then none are deleted.
            
            
            
            savebasis = copy.copy(iBasis)
            iBasis = self.clip_basis_outside_refinement(iBasis, 
                                                        inner)
            assert(np.any(savebasis == iBasis)),'iBasis failed at u = '.format(u)
            ssf = sum(iBasis)
            
            #savebasis = copy.deepcopy(oBasis)
            #oBasis = self.clip_basis_outside_refinement(oBasis, 
            #                                            active_parents) #active_parents) #apb
            #assert(np.any(savebasis == oBasis)),'oBasis failed at u = {}'.format(u)
            
            
            w = 1. - ssf#sum(iBasis[inner]) 
            level = self.level
            if self.rm is None:
                try:
                    self.rm = self.Pget(self.level)
                except:
                    self.rm = self.FindP(self.p,self.level+1)
                #self.rm = self.wgen(self.p,self.level+1)#go from L+1 -> to L # wgen == FindP
                
                #self.wQ = FindQ(self.p,self.level+1,'min')
            #if len(outer)>0:
            #    deleted = oBasis[outer].copy() #drop this down into the refined basis??
            oBasis[outer] = 0.

            e = sum(oBasis[apb])
            W = 0.

            if e !=0.:
                W = w/e
            
            #oBasis[apb] = W*oBasis[apb]
            
            #if np.linalg.norm(iBasis)>0.:
            #    oBasis[apb] = np.matmul(self.rm.T,iBasis)[apb]
                #oBasis = np.matmul(self.rm.T,iBasis)
                #oBasis[apb] = np.matmul(iBasis[apb], np.matmul(self.rm.T,iBasis)[apb])
            #if ssf<=0.:
            #    oBasis[apb] = oBasis[apb] - np.matmul(self.rm.T,iBasis)[apb]
            oBasis[apb] =  oBasis[apb] - np.matmul(self.rm.T,iBasis)[apb]
            #oBasis[apb] =  np.matmul(self.rm.T,iBasis)[apb]
            oBasis[oBasis < 0.] = 0.
            """
            This goes wrong for levels > 1 b/c 
            the truncation is only looking at the partial results,
            1 level down.
            """
            #oBasis = oBasis - np.matmul(self.rm.T,iBasis)
            #oBasis[apb] = np.matmul(self.rm.T,iBasis)[apb]
            
            #oBasis[apb] =  np.matmul(self.rm,oBasis)[apb]
            #oBasis[apb] =  oBasis[apb] - np.matmul(self.rm.T,iBasis)[apb]
            
            
#            needed! - otherwise it breaks the next level!
            savebasis = copy.copy(oBasis)
            oBasis = self.clip_basis_outside_refinement(oBasis, 
                                                        active_parents) #active_parents) #apb
            assert(np.any(savebasis == oBasis)),'oBasis failed at u = {}'.format(u)
#            
            #but this is late!  it cuts after scaling!!
            
            if self._verbose:
                print ''
                print 'level = ',self.level
                print 'potantially active this level = ',active_parents
                print 'active children = ',inner
                print 'apb = ',apb
                
            
        
            if self.checksums_:
                s1 = sum(oBasis[apb])  +sum(iBasis[inner])#sum(oBasis[apb]) + sum(iBasis[inner])
                print 's1 = ',s1
                s2 = sum(oBasis)  + sum(iBasis)
                print 's2 = ', s2
                print 'checksum: '+\
                   'u={},s1={}, s2={}, diff={}'.format(u,s1,s2,abs(s1-s2))
                
                assert(abs(float(s2-s1)) < self.tol),'failed checksum '+\
                   'u={},s1={}, s2={}, diff={}'.format(u,s1,s2,abs(s1-s2))
            
#            
#            savebasis = copy.copy(oBasis)
#            oBasis = self.clip_basis_outside_refinement(oBasis, 
#                                                        active_parents) #active_parents) #apb
#            assert(np.any(savebasis == oBasis)),'oBasis failed at u = {}'.format(u)
#            
            return oBasis, sum(oBasis)+ssf
        
        else:
            savebasis = copy.copy(oBasis)
            oBasis = self.clip_basis_outside_refinement(oBasis, 
                                                        active_parents) #active_parents) #apb
            assert(np.any(savebasis == oBasis)),'oBasis failed at u = {}'.format(u)
            
            if self._verbose:
                print ''
                print 'level = ',self.level
                print 'potentially active bottom this level = ',active_parents
            
            return oBasis, sum(oBasis)
    
    
    
    def CustomBasis(self,u,t):
        """
        This is not the way
        """
        i = self.FindSpan(u)
        k=self.k
        p=self.p
        U=t#self.t
        Nall = np.zeros((self.n,1),float)
        N = Nall[i-p:i+1]
        """
            i = span index
            u = point on the parameterization
            k = order of the basis function (passed for ease of code update)
            U = knot vector
            N = EmptyArray of Basis functions from [0...p], i.e. size k, at u.
            
        Output:
            N = Array of all Basis Function Values at u.  """
        #print N
        left=np.zeros((k),float)
        right=np.zeros((k),float)
        #p=k-1
        N[0]=1.0
        for j in range(1,k): #use k instead of p to match loop syntax for python vs C++
            left[j]  =u-U[i+1-j]
            right[j] =U[i+j]-u
            saved=0.0
            for r in range(0,j):
                temp  = N[r]/(right[r+1]+left[j-r])
                N[r]  = saved+right[r+1]*temp
                saved = left[j-r]*temp
            
            N[j]=saved
        return Nall
        

    def truncate_basis(self, u):
        """bounds is an interval!
        """
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

        """
        refined space: child.bounds

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

        """
           Missing proper return of the inner functions
        """
        """
            To try idea 2, 
            the basis evaluator must return
            the refined basis scaled by the deleted coarse basis
        """
        #inner = []
        #iBasis = []
        dummy_bounds = ia(0.,0.)
        if self.children is not None:
            child = self.children

            #""" Parent used??????
            active_parent_bounds_index = self.bounds.contains(u)
            len_apbi = len(active_parent_bounds_index)
            assert(len(active_parent_bounds_index)<2),'error: parent bounds active in two spans at once'
            if len_apbi==1:
                parent_bounds = self.bounds[active_parent_bounds_index[0]]
            else:
                parent_bounds = self.bounds[0]#dummy_bounds
                #assert(False),"parent_bounds = dummy_bounds! u = {}".format(u)
            #"""


            active_child_bounds_index = child.bounds.contains(u)
            len_acbi = len(active_child_bounds_index)
            assert(len(active_child_bounds_index)<2),'error: child bounds active in two spans at once'
            if len_acbi==1:
                child_bounds = child.bounds[active_child_bounds_index[0]]


                inner = child.enclosed_basis_in_range(child_bounds)
                iBasis = child.truncate_basis(u)#.TLMBasisFuns(u)
                #this is correct
                #because outer looks for the outer basis function(s) eclosed
                # in the bounds of the inner function's active range.
                outer = self.enclosed_basis_in_range(child_bounds) #to be eliminated
                oBasis = self.TLMBasisFuns(u)
                
            else:
                child_bounds = dummy_bounds
                inner = []
                iBasis = child.truncate_basis(u)##.TLMBasisFuns(u)

                #this is correct
                #because outer looks for the outer basis function(s) eclosed
                # in the bounds of the inner function's active range.
                outer = [] #to be eliminated
                oBasis = self.TLMBasisFuns(u) 
                # if none are enclosed, then none are deleted.


            #allowably active as defined by self.bounds
            active_parents = self.enclosed_basis_in_range(parent_bounds)


            #truncate:
                
            #active child bases used/useful?
            #acb = self.bases(child.active(u))

            #actully active parent basis 
            apb = self.bases(self.active(u)) #(not??) used when children are active
            
            #naive construction:
            if child_bounds.contains(u):
                
                
                

                """zero oBasis not enclosed in active_parents
                as these are not allowed past wherever the range is
                """
                """
                save = np.zeros_like(oBasis)
                save[active_parents] = copy.copy(oBasis[active_parents])
                oBasis = oBasis*0.
                oBasis[active_parents] = save[active_parents]
                #"""
                oBasis = self.clip_basis_outside_refinement(oBasis, 
                                                            active_parents)
                """zero iBasis components which are not in active_child:
                    yes - those must always be zero
                    inner is like active_parents
                    
                    inner is different than outer
                    outer is outer parents restricted to
                    just those parents enclosed by 
                    this limited interval: the inner bounds!
                """     
                """
                save = np.zeros_like(iBasis)
                save[inner] = copy.copy(iBasis[inner])
                iBasis = iBasis*0.
                iBasis[inner] = save[inner]
                #"""
                iBasis = self.clip_basis_outside_refinement(iBasis, 
                                                            inner)

                w = 1. - sum(iBasis[inner]) #says only these basis are live.
                                            # This is actually right!

                #deleted = sum(oBasis[outer])
                deleted = oBasis[outer].copy()
                if len(outer)>0:
                    oBasis[outer] = 0.

                e = sum(oBasis[apb])
                W = 0.

                if e !=0.:
                    W = w/e
                oBasis[apb] = W*oBasis[apb]
                if self._verbose:
                    print 'hello',u, sum(oBasis[apb])  +sum(iBasis[inner])
                
                """zero oBasis not enclosed in active_parents
                as these are never allowed
                """
                """
                save = np.zeros_like(oBasis)
                save[active_parents] = copy.copy(oBasis[active_parents])
                oBasis = oBasis*0.
                oBasis[active_parents] = save[active_parents]
                #"""
                
        
                if self.checksums_:
                    s1 = sum(oBasis[apb])  +sum(iBasis[inner])#sum(oBasis[apb]) + sum(iBasis[inner])
                    print 's1 = ',s1
                    s2 = sum(oBasis)  + sum(iBasis)
                    print 's2 = ', s2
                    print 'checksum: '+\
                       'u={},s1={}, s2={}, diff={}'.format(u,s1,s2,abs(s1-s2))
                    
                    assert(abs(float(s2-s1)) < self.tol),'failed checksum '+\
                       'u={},s1={}, s2={}, diff={}'.format(u,s1,s2,abs(s1-s2))
                
                

            else:  #child_bounds does not contain u
                

                """zero oBasis not enclosed in active_parents
                as these are never allowed
                
                Needed here to trucate outside bounds!
                and get the right CurvePoint at k=3
                """
                #iBasis = self.clip_basis_outside_refinement(iBasis, 
                #                                            inner)
                oBasis = self.clip_basis_outside_refinement(oBasis, 
                                                            active_parents)
                
                #"""
                #save = np.zeros_like(oBasis)
                #save[active_parents] = copy.copy(oBasis[active_parents])
                #oBasis = oBasis*0.
                #oBasis[active_parents] = save[active_parents]
                
                #"""
                #so no trucation is going to take place.
                #if self.bounds.contains(u):
                #   pass
                #else:
                #   pass # oBasis =  oBasis[" not inner"] incorrect b/c inner is for deleting oBasis[inner]
        elif self.children is None:
            #critical for getting the upper basis right:
            active_parent_bounds_index = self.bounds.contains(u) #which interval bounds box on this level?
            len_apbi = len(active_parent_bounds_index)
            assert(len(active_parent_bounds_index)<2),'error: parent bounds active in two spans at once'
            if len_apbi==1:
                parent_bounds = self.bounds[active_parent_bounds_index[0]]
            else:
                parent_bounds = dummy_bounds

            active_parents = self.enclosed_basis_in_range(parent_bounds) #ia(low,high)
            oBasis = self.TLMBasisFuns(u)
            #active_parents = self.enclosed_basis_in_range(parent_bounds)

            
            """zero oBasis components which are not in active_parents:
            """
            """
            save = np.zeros_like(oBasis)
            save[active_parents] = oBasis[active_parents]
            oBasis = oBasis*0.
            oBasis[active_parents] = save[active_parents]
            #"""
            oBasis = self.clip_basis_outside_refinement(oBasis, 
                                                        active_parents)
            #return oBasis

#            if self.bounds.contains(u):
#                return oBasis
#            else:
#                return oBasis*0.
        else:
            print 'error, how many children??'
            print self.children
        return oBasis


    def eval_thb_basis(self, u, basis_dict, L):
        """bounds is an interval
        #"""
        #        whichbounds = 0
        #        for i, el in bounds:
        #            if u in el:
        #                if whichbounds is None:
        #                    print 'foundsbounds {}, level {}'.format(i, L)
        #                    whichbounds = i
        #                else:
        #                    print 'error, overlapping bounds??'

        basis_dict[L] = self.truncate_basis(u)
        return  basis_dict


    def recursive_thb_basis(self, u, basis_dict, L):

        self.eval_thb_basis(u, basis_dict,L)
        if self.children is not None:
            self.children.recursive_thb_basis(u, basis_dict, L+1)
        else:
            print 'TODO! fix basis eval!'
        return
    

    def get_thb_basis(self, u):
        """bounds is an interval
        """
        if self.parent is not None: print 'warning, assembling thb basis starting above Omega^{0}'
        #ltot = self.sound_hierarchy()
        self.basis_dict = {}
        self.recursive_thb_basis(u, self.basis_dict, 0)
        return  self.basis_dict



    def eval_deprecated(self, u, ans=0.):
        span = self.FindSpan(u)
        #b = self.truncate_basis(u)
        b,ssf = self.truncate_basis2(u)
        ans = ans + np.dot(b[span-self.p:span+1].T,self.vertices[span-self.p:span+1,:])
        #ans = ans + np.dot(b.T,self.vertices)
        if self.children is not None:
            ans = self.children.eval(u, ans)
        return ans
    
    def eval1(self,u):
        b,ssf = self.truncate_basis2(u)
        return b
    
    def eval2(self, u, S=0.):
        """
        TODO:  check this THBsurface style
        """
        nh = self.sound_hierarchy()
        child = self
        S = 0.
        for el in range(nh):
            basis = child.eval1(u)
            S = S + np.dot(basis.T,child.vertices)
            if child.children is not None:
                child = child.children
            else:
                assert(nh==el),'nh != el!'
            
        #        basis = self.eval1(u)
        #        S += np.dot(basis.T,self.vertices)
        #        
        #        if self.children is not None:
        #            S =  self.children.eval2(u, S)
        #else:
            #return S
        return S
        
        
        
    
    def evs(self, u):
        ans = 0.
        child = self
        while child is not None:
            b,ssf = child.truncate_basis2(u)
            tmp = np.dot(b.T,child.vertices)
            #print child.level, tmp, ssf
            ans = ans + tmp
            child = child.children
        return ans
    
    def eval_basis(self, u):
        depth = self.sound_hierarchy()
        ans = 0.
        child = self
        ds = []
        vts = []
        for i in range(depth+1):
            b,ssf = child.truncate_basis2(u)
            tmp = np.dot(b.T,child.vertices)
            #print child.level, tmp, ssf
            ds.append(b)
            vts.append(child.vertices)
            child = child.children
        return ds, vts
    
    def check_ds(self, ds):
        s =0.
        for el in ds:
            s += sum(el)
        print s
        return 
    
    def ds_dot_vts(self, ds, vts):
        ans = 0.
        for el, vt in zip(ds,vts):
            ans += np.dot(el.T,vt)
        return ans
    
    




    def sound_hierarchy(self, level=-1):
        level += 1
        if self.children is not None:
            level = self.children.sound_hierarchy(level)
        return level




    #    def sum_thbasis(self, bounds, u):
    #        """bounds is an interval
    #        """
    #        inner = self.children.enclosed_basis_in_range(bounds)
    #        outer = self.enclosed_basis_in_range(bounds)
    #
    #        oBasis = self.TLMBasisFuns(u)
    #        iBasis = self.children.TLMBasisFuns(u)
    #
    #
    #        return


    def plot_thbasis(self,
                     interval = ia(0.,1.),
                     nump=1000, offset=0., scaler=None):
        """bounds is (currently) an interval
        
        dev
        ----------
        
         interval = ia(0.,1.)
         nump=1000
         offset=0.
         scaler=None
        """
        #p = self.p
        a = interval[0]
        b = interval[1]
        depth = self.sound_hierarchy()
        #bi = self.active(interval)
        s = np.linspace(a,b,nump,endpoint=True)
        #numb = bi[1] - bi[0] + 1
        basisl = []
        lcurve = self.get_bottom_curve()
        for i in range(depth+1):
            basisl.append(np.zeros((nump,lcurve.n)))
            lcurve = lcurve.parent
            
        for i,u in enumerate(s):
            #span = self.FindSpan(u)
            #self.BasisFuns(span,u,basis[i,span-p-bi[0]:span+1])
            basis_list= self.rthb_basis(u)
            for j in range(depth+1):
                basisl[j][i] = basis_list[j].T
            #basis[i][basis[i] < 0.] = 0. #crutch - don't rely on this!
        #plt.plot(s,)
        for j in range(depth+1):
            plt.plot(s,basisl[j][:]-offset)
            offset += 1.
        return #basis


    def plot_truncated_basis(self,
                             interval = ia(0.,1.),
                             nump=30, offset=0., scaler=None):
        """bounds is (currently) an interval
        """
        #p = self.p
        a = interval[0]
        b = interval[1]
        #bi = self.active(interval)
        s = np.linspace(a,b,nump,endpoint=True)
        #numb = bi[1] - bi[0] + 1
        basis = np.zeros((nump,self.n))
        for i,u in enumerate(s):
            #span = self.FindSpan(u)
            #self.BasisFuns(span,u,basis[i,span-p-bi[0]:span+1])
            bstuff,ssf = self.truncate_basis2(u)
            basis[i] = bstuff[:,0]
            #basis[i][basis[i] < 0.] = 0. #crutch - don't rely on this!
        #plt.plot(s,)
        plt.plot(s,basis[:]-offset)
        return #basis


    def plot_a_basis(self, blist,
                         interval = ia(0.,1.),
                         nump=30, offset=0., scaler=None,
                         bounds = None):
        """bounds is (currently) an interval
        """
        if bounds is None:
            bounds = self.bounds
        p = self.p
        a = interval[0]
        b = interval[1]
        s = np.linspace(a,b,nump,endpoint=True)
        basis = np.zeros((nump,self.n))
        for i,u in enumerate(s):
            basis[i][blist] = self.TLMBasisFuns(u)[blist,0]

        plt.plot(s,basis[:]-offset)
        return
    
    def plot3DmultiList(self, 
                        transverseList, 
                        longitudinalList):
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        
        self.biggest = 0.
        for curve in transverseList:
            ax = curve.plotcurvehierarchy(canvas = ax)
        for curve in longitudinalList:
            ax = curve.plotcurvehierarchy(canvas = ax)
        
        
        biggest = self.biggest
        axmax = max(self.extremes1d(0))
        aymax = max(self.extremes1d(1))
        if curve.dim ==3:
            azmax = max(self.extremes1d(2))
            biggest = max(axmax,aymax,azmax,biggest)
        elif curve.dim==2:
            biggest = max(axmax,aymax,biggest)
            
        #
        self.biggest = biggest
        #
        ax.set_xlim3d(-biggest*.0,1.*biggest)
        ax.set_ylim3d(-biggest*.5,.5*biggest)
        if curve.dim ==3:
            ax.set_zlim3d(-biggest*.5,.5*biggest)
        plt.show()
        return
    
    def plotcurve3d(self, 
                  interval = ia(0.,1.),
                  nump = 30,
                  canvas = None,
                  biggest=None):
        """
        plot a rbspline
        (THBspline curve)
        """
        colors = {0:'purple',
              1:'green',
              2:'blue',
              3:'white'}
        
        a = interval[0]
        b = interval[1]
        level = self.level
        cc = colors[level%3]
        s = np.linspace(a,b,nump,endpoint=True)
        cvec = np.zeros((nump,self.dim))
        for i,u in enumerate(s):
            cvec[i] = self.rthbeval(u) #eval on each level and sum up
            #cvec[i] = self.eval_new(u) #project fine eval there
            #cvec[i] = self.eval(u)
        
        if canvas is None:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            canvas = ax
            #canvas.show()
        
        canvas.plot(cvec[:,0], 
                    cvec[:,2],
                    cvec[:,1],
                    linestyle = "-",
                    color = cc)
        
        def pchl(chld, ci):
            clr = colors[ci]
            for bounds in chld.bounds:
                acv = chld.vertices[chld.enclosed_basis_in_range(bounds)]
                canvas.plot(acv[:,0],
                             acv[:,2], 
                             acv[:,1], 
                             color = clr,
                             marker = 's', 
                             alpha = .4)
            if chld.children is not None:
                ci = pchl(chld.children, ci)
            return (ci+1)%3
        
        ci = 0
        dm = pchl(self, ci)
        #
        if biggest is None: 
            try:
                biggest = self.biggest
            except:
                biggest = 0.
        axmax = max(self.extremes1d(0))
        aymax = max(self.extremes1d(1))
        if self.dim ==3:
            azmax = max(self.extremes1d(2))
            biggest = max(axmax,aymax,azmax,biggest)
        else:
            biggest = max(axmax,aymax,biggest)
        #
        self.biggest = biggest
        #
        canvas.set_xlim3d(-biggest*.0,1.*biggest)
        canvas.set_ylim3d(-biggest*.5,.5*biggest)
        if self.dim==3:
            canvas.set_zlim3d(-biggest*.5,.5*biggest)
        return canvas
    
    
    def plotcurvehierarchy(self, canvas=None,
                           biggest=None):
        curve = self.get_coarsest_curve()
        level = curve.level
        
        ax = canvas
        if biggest is None: 
            try:
                biggest = self.biggest
            except:
                biggest = 0.
                
                
        ax = curve.plotonecurve(level,
                                biggest=biggest,
                                canvas=ax)
        while curve.children is not None:
            curve = curve.children
            level = curve.level
            ax = curve.plotonecurve(level,
                                    biggest=biggest,
                                    canvas=ax)
            
        return ax
    
    def plotonecurve(self, level, 
                     biggest=None,
                     greville_abscissa=None,
                     isu=False,isv=False, 
                     nump=30, dim=3,
                     colorize_levels=True,
                     color_override=None,
                     plot_vertices=True,
                     canvas=None):
        """
        greville_abscissa : the parametric location in the perpindicular 
        parametric direction at which this curve 'lives'
        
        
        greville_abscissa = ga
        
        level = surf.level
        """
        if colorize_levels:
            colors = {0:'purple',
                  1:'green',
                  2:'blue',
                  3:'red',
                  4:'black'}
        else:
            colors = {0:'.65',
                      1:'.55',
                      2:'.28',
                      3:'.1'}
        cc = colors[level%3]
        if color_override is not None:
            cc = color_override
            
            
        if canvas is None:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            canvas = ax
            #canvas.show()
        
            
        curve = self
        
        if self.children is not None:
            ilist = self.get_active_disjoint_bounds()
        else:
            ilist = self.bounds.bounds
        
        if len(ilist)>0:
            for intval in ilist:
                a = intval[0]
                b = intval[1]
                s = np.linspace(a,b,nump,endpoint=True)
                cvec = np.zeros((nump,curve.dim))
                for i,u in enumerate(s):
                    #this check must be done outside 
                    # to ensure a consistent plotting interval 
                    # once we get to this point
                    #ck = self.child.bounds.get_box_containing()
                    #if ck.isempty:
                    cvec[i] = curve.CurvePoint(u) #assumes that this is a THB curve
                if curve.dim==2:
                    canvas.plot(cvec[:,0], 
                                cvec[:,1],
                                linestyle = "-",
                                color = cc)
                elif curve.dim ==3:
                    canvas.plot(cvec[:,2], 
                                cvec[:,0],
                                cvec[:,1],
                                linestyle = "-",
                                color = cc)
                if plot_vertices:
                    #vi = curve.bases(curve.active_knots(intval))
                    vi = curve.bases(curve.active_basis(intval))
                    cpts = curve.vertices[vi]
                    if curve.dim==2:
                        canvas.plot(cpts[:,0], 
                                    cpts[:,1],
                                    marker = "o",  
                                    alpha = .5,
                                    color=cc)
                    elif curve.dim ==3:
                        canvas.plot(cpts[:,2], 
                                    cpts[:,0],
                                    cpts[:,1],
                                    marker = "o",  
                                    alpha = .5,
                                    color=cc)
        #
        if biggest is None: biggest=0.
        axmax = max(self.extremes1d(0))
        aymax = max(self.extremes1d(1))
        if curve.dim==2:
            biggest = max(axmax,aymax,biggest)
        if curve.dim==3:
            azmax = max(self.extremes1d(2))
            biggest = max(axmax,aymax,azmax,biggest)
        #
        self.biggest = biggest
        #
        canvas.set_xlim3d(-biggest*.0,1.*biggest)
        canvas.set_ylim3d(-biggest*.5,.5*biggest)
        if curve.dim==3:
            canvas.set_zlim3d(-biggest*.5,.5*biggest)
        canvas.biggest = biggest
        return canvas
    
    
    
    
    
    def get_active_disjoint_bounds(self):
        """returns the disjoint set
        of active intervals 
        at this level, exluding the next level finer
        
        assumptions
        ----------
        -nested spaces!
        
        
        
        returns:
        ----------
         bounds list which excludes contributions from the next
         level down
         
        replaces
        ----------
         use of exclusive_disjunctive_elimination
         with exlude_from_bounds - uses ia list method
                             iafuncs.interval_list_diff
         instead!
         
        notes
        ----------
         it is not recursive.  By enforcing nested lists of
         bounds, we can be sure that by excluding the 
         bounds one level finer, we also catch all levels
         finer than that!
        
        """
        
        return self.bounds.exlude_from_bounds(self.children.bounds)
    
    
    def plotcurve(self, 
                  interval = ia(0.,1.),
                  nump = 30,
                  canvas = None):
        """
        plot a rbspline
        (THBspline curve)
        """
        if self.dim==3: return self.plotcurve3d(interval, nump, canvas)
        colors = {0:'blue',
                  1:'red',
                  2:'green'}
        a = interval[0]
        b = interval[1]
        s = np.linspace(a,b,nump,endpoint=True)
        cvec = np.zeros((nump,self.dim))
        for i,u in enumerate(s):
            cvec[i] = self.rthbeval(u) #eval on each level and sum up
            #cvec[i] = self.eval_new(u) #project fine eval there
            #cvec[i] = self.eval(u)
        
        if canvas is None:
            canvas = plt
        if self.dim==2:
            canvas.plot(cvec[:,0], cvec[:,1])
        
        def pchl(chld, ci):
            clr = colors[ci]
            for bounds in chld.bounds:
                acv = chld.vertices[chld.enclosed_basis_in_range(bounds)]
                canvas.plot(acv[:,0],
                             acv[:,1], 
                             color = clr,
                             marker = 's', 
                             alpha = .4)
            if chld.children is not None:
                ci = pchl(chld.children, (ci+1)%3)
            return ci%3
        
        ci = 0
        dm = pchl(self, ci)
        return
    
    def plotcurve_no_vertices(self, 
                  interval = ia(0.,1.),
                  nump = 30):
        """could get to be widly overcalled for some 
        curves if used in surface plotting??
        """
        return
    
    
    def plotcurve_onelevel(self, 
                  interval = ia(0.,1.),
                  nump = 30,
                  color_ci = 3,
                  canvas=None):
        if canvas is None:
            canvas = plt
        
        colors = {0:'blue',
                  1:'red',
                  2:'green',
                  3:'black'}
        a = interval[0]
        b = interval[1]
        s = np.linspace(a,b,nump,endpoint=True)
        cvec = np.zeros((nump,self.dim))
        for i,u in enumerate(s):
            cvec[i] = self.rthbeval(u)

        
        acv = self.vertices[self.enclosed_basis_in_range(
                                            self.bounds[0])]
        clr = colors[color_ci]
        
        if self.dim ==2:
            canvas.plot(cvec[:,0], cvec[:,1],color= clr)
        elif self.dim == 3:
            canvas.plot(cvec[:,0], 
                        cvec[:,1],
                        cvec[:,2],
                        color= clr)
        
        canvas.plot(acv[:,0],
                     acv[:,1], 
                     color = clr,
                     marker = 's', 
                     alpha = .4)
        
        return
        
        
    
    def plotvertices(self):
        plt.plot(self.vertices[:,0],
                     self.vertices[:,1], 
                     color = 'blue',
                     marker = 's', 
                     alpha = .4)



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

def make_knot_vertices(num, k=4, dim=2):
    """utility to rturn a set of starting points
        start   : starting point 
        end     : ending point
        num     : numberof vertices in the real curve
            (not the number of ~vertices~ here
            since this is not a real curve we are
            making!)
    """
    vertices = kv.make_knot_vector(k,num-1,np.array([0.,1.]) ,type='open')
    array = np.zeros((np.shape(vertices)[0],dim),float)
    array[:,0] = vertices.copy()
    array[:,1] = vertices.copy()
        
    return array
    
    
 
class knot_spline(rbspline):
    """
        nonuniform subdivision for THB refinement 
        subdivides the knots
        Let's make that happen with minimal effort
    """
    
    def knot_support(self, knot_index):
        #if knot_index < self.p+1:
        #    num_zeros = self.p+1 - knot_index
        
        #range(knot_index,knot_index+p+1) !!??
        return range(knot_index,knot_index+self.p+2) 
    
    
    def knot_support_guess(self, knot_index):
        #if knot_index < self.p+1:
        #    num_zeros = self.p+1 - knot_index
        
        #range(knot_index,knot_index+p+1) !!??
        return range(knot_index,knot_index+self.p+1) 
    
    
    def knot_refinement(self, knot_index):
        """
            Use Knot Insertion on Knots!
            It's KNOTS!!
        inputs:
            knot index : == control point index
        """
        k = self.k
        Nseg        = self.n-self.k+1 # number of k order 
                            #bezier segments in the curve
        original_knots = self.t.copy()
        tau = self.t[self.knot_support(knot_index)]
        new_knots = []
        for i in range(Nseg):            
            iknot = (tau[k+i-1]+tau[k+i])/2.
            new_knots.append( iknot)#(tau[i-1]+tau[i])/2. )
        a1=thb.subdivide_basis_not_SABIN(self.p,
                                      tau,
                                      tp=new_knots)
        a2=thb.subdivide_basis1(self.p, 
                             self.vertices,
                             tau,
                             tp=new_knots)
        a3=thb.subdivide_basis(self.p, 
                            self.vertices,
                            tau,
                            tp=new_knots)
        return 
    
    def knot_refinement_not_fixed(self, knot_index):
        """
            Use Knot Insertion on Knots!
            It's KNOTS!!
        inputs:
            knot index : == control point index
        """
        k = self.k
        Nseg        = self.n-self.k+1 # number of k order 
                            #bezier segments in the curve
        original_knots = self.t.copy()
        tau = self.t[self.knot_support(knot_index)]
        #tau = self.t[self.knot_support_guess(knot_index)]
        #new_knots = []
        #for i in range(1,len(tau)):
        #    new_knots.append( (tau[i-1]+tau[i])/2. )
        
        for i in range(Nseg):            
            iknot = (tau[k+i-1]+tau[k+i])/2.
            new_knots.append( iknot)
            
        a1=thb.subdivide_basis_not_SABIN(self.p,
                                      tau,
                                      tp=new_knots)
        return 


if __name__ == '__main__':
    k=4
    nump=100
    start = [0.,12.]
    end = [12.,0.]
    num = 4
    #num = 11#7# #makes c0c's 5 basis function be enclosed in ia(.2,.8)
    vv = spline.curvalinear_vertices(start,end,num)

    
    #knot_curve = spline.linear_vertices([0.,0.],[1.,1.],4)
    #knot_curve = make_knot_vertices(num=4,k=4)
    #kcv = knot_spline(knot_curve,k,nump)
    knot_index = 0
    #self = kcv
    #c0 = hb.LevelSpline(vv,k,nump)rbspline
    c0 = rbspline(vv,k,nump)
    #c0.print_properties()

    #c0i = c0.bounds_refinement(bounds=ia(.75,.8),
    #                           knots=[.77,.775,.78])
    #c0.print_properties()

    c0c = copy.copy(c0)
    c0c.hierarchical_obsvrs = []
    c0c = c0.dyadic_refinement()
    c0c = c0c.dyadic_refinement()
    c0c = c0c.dyadic_refinement()
    #c0c = c0c.dyadic_refinement()
    c0c.parent = None #TLM crucial to make it less deep a tree!
    c1 = c0c.dyadic_refinement()
    #c0.print_properties()

    #c1.t[4] = .54
    c1i = rbspline(c1.vertices,
                         k,
                         nump,
                         t=c1.t)

    #c1i = c1i.bounds_refinement(bounds=ia(.75,.8),
    #                           knots=[.77,.775,.78])
#
#    c1.plot_basis_interval([0.,1.],100,.0)
#    c0i.plot_basis_interval([.6,.9],10000,1.)
#    c0.plot_basis_interval([0.,1.],100,2.)


    self = c0c
    #print c1.enclosed_basis_in_range(ia(.4,.72))



    #print c1 is c0c.children
    #print c1.parent is c0c

    """******************************************************
    """
    """bounds = ia(.25,.75)#"""
    #bounds = ia(.33,.67)
    bounds = ia(.2,.8)
    bounds1 = ia(.1,.875)
    bounds2 = ia(.3125,.71875)

    #print c0c.enclosed_basis_in_range(bounds)
    #print c1.enclosed_basis_in_range(bounds)

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
    
    
    #c2.bounds = BoundsList([ia(.1,.28),ia(.3,.5),ia(.51,.72)])
    #c1.bounds = BoundsList([bounds1])
    #
    #
    #
    if True:
    
        #inner = []
        #for child in self.children:
        #    inner += child.enclosed_basis_in_range(bounds)
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
    
    
        #middle (c2). only active part:
        mmplt = mplt+1.
        #c2.bounds = BoundsList([bounds2])
        c2.bounds = BoundsList([ia(.1,.28),ia(.3,.5),ia(.51,.72)])
        #c2.plot_a_basis(blist=c2.enclosed_basis_in_range(bounds2),
        #                nump=250,
        #                offset=mmplt)
        #"""
        # rethink plotting?
        #c2.plot_truncated_basis(ia(0.,1.),
        #                        nump=250,
        #                        offset=mmplt)
        #"""
    
    
        bbplt = mmplt+1.
        #c1.bounds = BoundsList([bounds]) #blows up right now
        #c1.bounds = BoundsList([bounds2])
        #c1.bounds = BoundsList([ia(.2,.8)])
        #c1.bounds = BoundsList([ia(.0,.35), ia(.65,1.)])
    
        #truncated middle:
        #c1.bounds = BoundsList([ia(.15,.4), ia(.6,.85)])
        #c1.bounds = BoundsList([ia(.0,.48), ia(.52,1.)])
        #c1.bounds = BoundsList([ia(.2,.8)])
        c1.bounds = BoundsList([bounds1])
        #c1.plot_truncated_basis(nump=250, offset=bbplt)
    
        bbbplt = bbplt+1.
        #truncated bottom:
        #c0c.bounds = BoundsList([ia(.2,.35), ia(.65,.8)])
        #c0c.bounds = BoundsList([ia(0.,1.)])
        #c0c.bounds = BoundsList([ia(0.,.2), ia(.8,1.0)])
        #c0c.plot_truncated_basis(nump=250, offset=bbbplt)
    
        c0c.plot_thbasis(offset = mmplt)
    
        cplt = bbbplt+1.
        #bottom:
        #c0c.plot_truncated_basis(ia(.5,.51),offset=cplt,nump=1000) #none truncated
        #I decided truncation had to be ~automagic~...
        #breaking that beautiful functionality
        #
        c0c.bounds = BoundsList([ia(0.0, 1.0)])
        c0c.plot_a_basis(blist=c0c.enclosed_basis_in_range(ia(0.,1.)),
                        nump=250,
                        offset=cplt)
    
        
        if True:
            tol = 1.e-1
            u = .21
            self = c1
            callme = np.linspace(0.,1.,100,endpoint=True)
            for u in callme:
                b = c0c.truncate_basis(u)
                #bbb = c1.TLMBasisFuns(u) #inner
                bb = c1.truncate_basis(u)
                bbb = c2.truncate_basis(u)
                #print sum(bbb[inner]) + sum(bb)
                
                
                b,a = c0c.truncate_basis2(u)
                bb,aa = c1.truncate_basis2(u)
                bbb,aaa = c2.truncate_basis2(u)
                sum_ = sum(b) + sum(bb) + sum(bbb)
                #print u, 'sum = ',sum_
                if c0c.checksums_:
                    assert(abs(1.-sum_)<tol),'sum is not 1! but {}, u = {}'.format(sum_, u)
                    u = 0.323232323232
            """
            bounds = ia(0.125, 0.875)
            u=.2
            a = c0c.truncate_basis(bounds, u)
            b = c1.truncate_basis(bounds, u)
            c = c2.truncate_basis(bounds, u)
        
            gg = c0c.get_thb_basis(bounds, u)
        
            sum_ = 0.
            for el in range(len(gg)):
                sum_ += sum(gg[el])
            print sum_
            #"""
        
            la = BoundsList([bounds,bounds2])
            lb = BoundsList([ia(.5,1.4),ia(1.,2.)])
            
            u=0.30303030303 #worst case!  
            #
            #  during function:
            #               c0c.truncate_basis
            # there is a 
            #  ~small?~ screwup in iBasis[inner]   ??
            #
            #  due to function:
            #               enclosed_basis_in_range
            
            #u=.4
            u=.777
            ans = c0c.eval(u)
            anso = c0c.CurvePoint(u)
            ansoo = c1.CurvePoint(u)
            ansooo = c2.CurvePoint(u)
            #discrepency!
            
            #for u in callme:
            #    print u, np.linalg.norm(c0c.eval(u)-c0c.CurvePoint(u))
                
             
            #""" ~always usually certainly good~
            #del(c0c.children.children)
            #del(c0c.children)
            cd0 = copy.deepcopy(c0c)
            cd1 = cd0.dyadic_refinement()
            cd1.children = None
            cd1.bounds = BoundsList([ia(0.125, 0.875)])#BoundsList([ia(.3,.777)])#
            cd0.children = cd1
            #cd1.bounds = BoundsList([ia(.3,.777)])#
    #        cd0.plot_truncated_basis(nump=250, offset=-1.)
    #        cd1.plot_truncated_basis(nump=250, offset=1.)
    #        
    #        cd0.plot_a_basis(
    #                        blist=cd0.enclosed_basis_in_range(ia(0.,1.)),
    #                        nump=250,
    #                        offset = -2.)
    #        cd1.plot_a_basis(
    #                        blist=cd1.enclosed_basis_in_range(ia(0.,1.)),
    #                        nump=250,
    #                        offset = 2.)
            #cd1.bounds = BoundsList([ia(.25,.777)])
        if False:
            #"""
            u=.777
            ans = cd0.eval(u)
            anso = cd0.CurvePoint(u)
            ansoo = cd0.children.CurvePoint(u)
            
            tol = 8.e-6
            for u in callme:
                dd = cd0.eval(u)-cd0.CurvePoint(u)
                nmdd = np.linalg.norm(dd)
                #print u, nmdd
                assert(abs(nmdd)<tol),'difference found at u = {}'.format(u)
             
            #"""
            #"""
            u = .4
            self = c0c
            
            tol = 8.e-1
            for u in callme:
                dd = c0c.eval(u)-c0c.CurvePoint(u)
                nmdd = np.linalg.norm(dd)
                #print u, nmdd
                assert(abs(nmdd)<tol),'difference found at u = {}'.format(u)
            #"""
            
            
            for u in callme:
                print '\nu=',u
                ds,vts = c0c.eval_basis(u)
                c0c.check_ds(ds)
                ans = c0c.ds_dot_vts(ds,vts)
                assert(np.linalg.norm(ans - c0c.CurvePoint(u))<c0c.tol)
                #print 'ans = ',ans, c0c.CurvePoint(u)
            
            u=.777
            print '\nu=',u
            ds,vts = c0c.eval_basis(u)
            c0c.check_ds(ds)
            ans = c0c.ds_dot_vts(ds,vts)
            #print 'ans = ',ans, c0c.CurvePoint(u)
            assert(np.linalg.norm(ans - c0c.CurvePoint(u))<c0c.tol)
                
            """
             Kiss dissertation:  p 58 convertison to 
             tensor product patches...
            """
            u=.33
            a1,b1 = c1.return_active_basis_and_indices(u)
            a2,b2 = c2.return_active_basis_and_indices(u)
            
        ##plt.close()
#        gg = THBsurface()
#        gg.initialize_grid()
#        gg.dyadic_loft_refinement(0) 
#        gg.dyadic_loft_refinement(1)
        #gg.dyadic_loft_refinement(2) 
        
#        gg.dyadic_loft_refinement(2,
#                                  ubounds = BoundsList([ia(.2,.8)]),
#                                  vbounds = BoundsList([ia(.25,.75)]) 
#                                  )
#        gg.refine_grid(nL_refine_depth=1,
#               ubounds = BoundsList([ia(.2,.8)]),
#               vbounds = BoundsList([ia(.25,.75)]) 
#               )
        #gg.plotpts()
        
        
        #c2.bounds = BoundsList([ia(.1,.28),ia(.3,.5),ia(.51,.72)])
        #dplt = cplt+1.
        #c2.plot_truncated_basis(nump=250,
        #                offset=dplt)
        
        test_interpolated = False
        if test_interpolated:
            cic = rbspline(c0c.vertices,4,150,interpolate=True)
            c0c.plotcurve_detailed()
            cic.plotcurve()
        
        
        #        gg = THBsurface()
        #        gg.initialize_grid()
        #        self = gg
        u = 1.
        v = 1.
        ubounds=BoundsList([ia(0.0, 1.0)])
        vbounds=BoundsList([ia(0.0, 1.0)])
        level = 0
        
        
        print 'checking .777 '
        u = .777
        b1,s1 = c1.truncate_basis2(.777)
        b0,s0 = c0c.truncate_basis2(.777)
        print sum(b1) + sum(b0)
        
        
        
        cd0 = copy.deepcopy(c0c)
        cd1 = cd0.dyadic_refinement()
        cd1.children = None
        cd1.bounds = BoundsList([ia(0.125, 0.875)])#BoundsList([ia(.3,.777)])#
        cd0.children = cd1
        
        
        self = c0c
        u=.33