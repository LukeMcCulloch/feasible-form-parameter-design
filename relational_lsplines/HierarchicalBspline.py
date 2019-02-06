# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 19:02:54 2016

@author: lukemcculloch

Key thought:  

Linear decomposition holds for 
    fairness norms, area, moments
    into reference + offset
    so I can "tack on" an area bounded by the 
    new part of the curve only and it makes sense
    do specify local constraints (area, Xc, Yc)
    by update of control point offsets only.
For Tangent and curvature, locality holds
    and these can be computed using...
    and updated via...
"""
import numpy as np
import matplotlib.pyplot as plt
import copy
import curve as spline
import KV as kv
import knot_operations as ko

from extended_interval_arithmetic import ia

class LevelSpline(spline.Bspline):
    
    def __init__(self, vertices, k, nump, 
                 t=None, trange=None, bounds = ia(0.,1.) ):
        super(LevelSpline, self).__init__(vertices, k, nump, t=None, trange=None )
        
        if t is None:
            self.t_copy         = kv.make_knot_vector(k,len(vertices)-1,self.trange,type='open')
        else:
            self.t_copy         = t 
        self.t                  = copy.deepcopy(self.t_copy)
        self.bounds                 = bounds    # active parametric range
        self.upper_boundary_knots   = None      # knots which bound the coarser level of detail for C2 continuity
        self.lower_boundary_knots   = None      # knots which bound the finer detail level for C2 continuity
        self.level                  = 0         # 0 is the coarsest
        self.interior_knots         = []        # ?
        self.parent                 = None      # need to re-write as quad tree...?
        self.children               = []        # need to re-write as quad tree...?
        self.vertex_map             = [ [a,b] for a, b in zip(range(self.n),self.t[self.k-2:self.n+2]) ] # knots where basis live?
        self.vertex_map[1][1]       = (self.t[self.k-1]+ self.t[self.k])/2.
        self.vertex_map[-2][1]       = (self.t[self.n]+ self.t[self.n-1])/2.
        self.offset_reference       = {}        # Forsey Offset style...? 
        for el in self.vertex_map:
            self.offset_reference[el[0]] = [el[0],el[1]]
        self.active_indicies        =  set(self.offset_reference.keys())  # portions of the spline space that have associated vertices
        self.offsets                = None      # offset from lower level representation to this representation
        self.base_curve             = True      # TODO: if this is used anywhere it needs to be false for non base levels
        return
    
    def compute_Greville_Abscissae(self):
        """re-added in Nov 2016
        where did this go?
        """
        u = self.t
        p = self.p
        self.ga_v = np.zeros(self.n,float)
        
        self.ga_v1 = np.zeros(self.n,float)
        for i in range(self.n):
            self.ga_v1[i] = np.sum(u[i+1:i+p+1]) #sum stops 1 short just like loop
        self.ga_v1 = self.ga_v1/self.p
            
        #        self.ga_v2 = np.zeros(self.n,float)
        #        for i in range(self.n):
        #            for k in range(1,self.p+1):
        #                self.ga_v2[i] += u[i+k]
        #            self.ga_v2[i] = self.ga_v2[i]/self.p
        #        
        #        self.ga_v = (self.ga_v1+self.ga_v2) /2.
        return
        
    @property
    def offsets(self):
        return self._offsets
    
    @offsets.setter
    def offsets(self, *args):
        offsets =  args[0]
        self._offsets = offsets
        return 
    
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
        
    def compute_vertex_map(self):
        self.vertex_map             = [ [a,b] for a, b in zip(range(self.n),self.t[self.k-2:self.n+2]) ]
        self.vertex_map[1][1]       = (self.t[self.k-1]+ self.t[self.k])/2.
        self.vertex_map[-2][1]       = (self.t[self.n]+ self.t[self.n-1])/2.
        return
        
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
        return rspline
        
    def active_knots(self, bounds=None,spans=None):
        """ tells which knots and (indirectly control points)
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
    
    def get_basis_index(self, knot_index):
        """ See, e.g. Birks knot insertion notes
        """
        return knot_index - self.k + 2
        
    def shift_reference(self, starting_index, plus_shift):
        new_dict = copy.deepcopy(self.offset_reference)
        for i in range(starting_index,self.n):
            del new_dict[i]
            new_dict[i+plus_shift]  = self.offset_reference[i]
        self.offset_reference       = new_dict
        self.active_indicies        = set(self.offset_reference.keys())
        return
        
    def bounds_refinement(self, bounds, knots=None):
        """
        -input a vector of knots and the bounding interval
        -output a copy of the old curve, but with these knots inserted
        
        The point of maximum influence of a vertex
        is necessarily a place where the first derivative 
        of it's basis function
        is zero.  Because that's the peak of the basis function!
        
        offset_reference :
            -The ith vertex goes with the ith basis function.  That simple.
            -but the rth - k + 2 th knot goes with those!+
            (TLM 2017 : ith + k - 1 basis?)
        
        assumes : 
            -bounds are listed [lower , upper]
            -all knots are listed from low to high
            -inserted knots must always increase, or extra machinary needed.
        """      
        nbpoints     = 4
        if knots:# is not None:
            nipoints     = len(knots)
        else:
            nipoints = 0
        k           = self.k
        p           = self.p
        n           = self.n
        bi = bounds[0]
        be = bounds[1]
        si = self.FindSpan(bi)
        se = self.FindSpan(be)
        """ Span location is in front of the ith span
        """
        
        if si != se: 
            print 'interior point already extant - TODO?'
            print 'algorithm good???'
        relknots = range(si,se+2)
        #Nknots      = n+k
        #Nintknots   = n-k
        Nseg        = n-k+1
        #nknots = np.zeros((Nseg,1),float)
        
        #new level of the hierarchy
        rspline                     = copy.deepcopy(self)
        rspline.offset_reference    = {}
        rspline.level               = self.level + 1
        self.upper_boundary_knots = [bi,bi,be,be]
        for i in range(nbpoints):
            #nknots[i] = (t[k+i]+t[k+i+1])/2.
            iknot = self.upper_boundary_knots[i] #(t[k+i-1]+t[k+i])/2.
            tnew, vertices, ir = ko.knot_insertion_hierarchical(
                                                   rspline.k,
                                                   rspline.vertices,
                                                   rspline.t,
                                                   iknot)
            print 'vertices {}:{} will change'.format(ir-p,ir)
            print vertices
            #rspline.t           = tnew
            rspline.vertices    = vertices
            rspline.t           = tnew
            rspline.lower_boundary_knots = [bi,bi,be,be]
        
        
        rspline.bounds = bounds
        rspline.active_offsets = []
        if knots is not None:
            for knot in knots:
                tnew, vertices, ir = ko.knot_insertion_hierarchical(
                                                   rspline.k,
                                                   rspline.vertices,
                                                   rspline.t,
                                                   knot)
                
                ith_kvtx = ir-rspline.k+2
                rspline.offset_reference[ith_kvtx]    = [ith_kvtx,knot]
                rspline.active_offsets.append(ith_kvtx)
                rspline.vertices                            = vertices
                rspline.t                                   = tnew
                rspline.interior_knots.append(knot)
                
        rspline.offsets = np.zeros(
                    (rspline.n,rspline.dim),float)
        rspline.compute_vertex_map()
        self.hierarchical_obsvrs.append(rspline) #causes an issue if one curve is used as a factory
        rspline.parent = self
        rspline.active_offsets = np.asarray(rspline.active_offsets)
        rspline.active_indicies = set(rspline.offset_reference.keys())
        self.children.append(rspline)
        print 'shift reference'
        self.shift_reference(starting_index=si-p+1,
                             plus_shift = nipoints + nbpoints )
        return rspline
    
    def update_interior_offset(self, index, offsets):
        if index in self.offset_reference:
            self.offsets[index] = offsets
        else:
            print 'error, index is outside refinement!'
        return
    
    def update_all_interior_offsets(self, offsets):
        if self.base_curve:
           self.update_all_reference_values(offsets) 
        s = self.active_offsets[0]
        e = self.active_offsets[-1]
        self.offsets[s:e+1] = offsets
        return
    
    def update_interior_reference_value(self, index, value):
        print 'should only be used by the base curve'
        assert(self.base_curve)
        self.vertices[index] = value
        return
    
    def update_all_reference_values(self, value):
        print 'should only be used by the base curve'
        assert(self.base_curve)
        self.vertices = value
        return
        
    
        
    def propagate_update(self):
        print 'Updating HBspline!'
        k           = self.k
        p           = self.p
        n           = self.n
        dcurve = copy.deepcopy(self.parent)
        dcurve.hierarchical_obsvrs = []
        dcurve.children = None
        # boundary vertices are tied to the parent
        for iknot in self.lower_boundary_knots:
            tnew, vertices, ir = ko.knot_insertion_hierarchical(
                                                   dcurve.k,
                                                   dcurve.vertices,
                                                   dcurve.t,
                                                   iknot)
            dcurve.vertices    = vertices
            dcurve.t           = tnew
        
        # interior vertices are tied to the parent by the offset
        if self.interior_knots is not None:
            for iknot in self.interior_knots:
                tnew, vertices, ir = ko.knot_insertion_hierarchical(
                                                       dcurve.k,
                                                       dcurve.vertices,
                                                       dcurve.t,
                                                       iknot)
                dcurve.vertices    = vertices
                dcurve.t           = tnew
            
        #now update!
        self.vertices[:] = self.offsets[:] + dcurve.vertices[:]
        
        print 'updating'
        self.allCurveArray()
        return
        
    #def get_level_from_vertex(self, index):
    #    #fi = self.bounds[0]
    #    #fe = self.bounds[1]
    #    level = None
    #
    #    return level
        
    def plotLevel(self,bounds = None, color='black'):
        if bounds is None:
            fi = self.bounds[0]
            fe = self.bounds[1]
        else:
            fi = bounds[0]
            fe = bounds[1]
        lt = np.where(self.s > fi)
        gt = np.where(self.s < fe)
        interior = np.intersect1d(lt,gt)
        s = [fi] 
        r = [self.CurvePoint(fi)]
        for el in interior:
            s.append(self.s[el]) 
            r.append(self.r[el])
        s.append(fe)
        r.append(self.CurvePoint(fe))
        r=np.asarray(r)
        plt.plot(r[:,0],r[:,1])
        plt.show()
        return
    
    def print_vertex_map(self):
        for el in self.vertex_map:
            print el
        return
    def print_offset_map(self):
        for el in self.offset_reference:
            print el, self.offset_reference[el]
            
    def print_properties(self):
        print '# knots = {}'.format(len(self.t))
        print '# vertices = {}'.format(self.n)
            
def translate(*args):
    curve = args[0]
    if len(args)==3:
        translate = [args[1],args[2]]
    else:
        translate = args[1]
    v1 = copy.deepcopy(curve.vertices)
    v1 = v1 + translate
    curve.vertices = v1
    return curve

class BsplineHierarchy(object):
    """A hierarchical B-spline as a list 
        from low detail 
        to high detail
    """
    def __init__(self, splines):
        self.SplineHierarchy = splines
        return
        
if __name__ == '__main__':
    basic = False
    interactive = True
    if False:
        k=4
        nump=30
        start = [0.,12.]
        end = [12.,0.]
        num = 10
        vv = spline.curvalinear_vertices(start,end,num)
        
        vb = spline.Bspline(vv,k,nump)
        h0spline = LevelSpline(vv,k,nump)
        
        h1spline = h0spline.bounds_refinement([.75,.8])
        
        h0spline.plotcurve_detailed()
        h1spline.plotcurve_detailed()
        self = h1spline
        a=h0spline
        b=h1spline
        print self.active_knots()
        print self.active_Bik_interval(u=.75)
        print self.active_Bik_interval(u=.8)
        
        av = copy.deepcopy(a.vertices)
        av[3] = [5.,5.]
        av[8] = [4.,4.]
        h0spline.vertices = av
        h0spline.plotcurve_detailed()
        h1spline.plotcurve_detailed()
        
    if basic:
        k=4
        nump=100
        start = [0.,12.]
        end = [12.,0.]
        num = 10
        vv = spline.curvalinear_vertices(start,end,num)
        
        vb = spline.Bspline(vv,k,nump)
        h0spline = LevelSpline(vv,k,nump)
        
        h1spline = h0spline.bounds_refinement(bounds=[.75,.8],
                                              knots=[.77,.775,.78])
        
        #h0spline.plotcurve_detailed()
        #h1spline.plotcurve_detailed()
        self = h1spline
        a=h0spline
        b=h1spline
        print self.active_knots()
        print self.active_Bik_interval(u=.75)
        print self.active_Bik_interval(u=.8)
        
        for key in h1spline.offset_reference:
            #index = h1spline.get_basis_index(key) #junk
            #h1spline.offsets[key] = [5.,.5]
            h1spline.update_interior_offset(key, [5.,.5])
        
        
        av = copy.deepcopy(a.vertices)
        av[3] = [5.,5.]
        av[8] = [4.,4.]
        h0spline.vertices = av
        h0spline.plotcurve_detailed()
        h1spline.plotcurve_detailed()
        
        if False:
            for i in range(nump):
                print i, h0spline.s[i], h0spline.r[i] - h1spline.r[i]
        
        plt.close()        
        h0spline.plot(vertices = True)
        h1spline.plot(vertices = False)
        
        a = h0spline
        b = h1spline
    if interactive:
        k=4
        nump=100
        start = [0.,12.]
        end = [12.,0.]
        num = 10
        vv = spline.curvalinear_vertices(start,end,num)
        
        vb = spline.Bspline(vv,k,nump)
        h0spline = LevelSpline(vv,k,nump)
        
        h1spline = h0spline.bounds_refinement(bounds=[.75,.8],
                                              knots=[.77,.775,.78])
                                              
        bsh = BsplineHierarchy([h0spline,h1spline])