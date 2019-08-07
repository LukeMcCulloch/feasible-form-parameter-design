# 20130628
# TLM resturctured codes
# P&T basis
#TODO: curve plotter that normalized by height and length separately.
#
# array stored
# Automatically Differentiated
# Lagrange Optimized
# Basis - Splines

# July 2nd 2013, 4 am - placement of live basis values into full size vector
# July 8, 2013 implemented P&T derivatives through ADLspline methods

# August:  Added AD differential geometry via AD_FrenetSerret
#           added support for binary VTK files via EVTK

import numpy as np
# imports for plotting:
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from time import clock, time
from KV import make_knot_vector
import routines as routines  #crossproduct, make_knot_vector

from poly_integration import ibasis_quad2,   \
                                Mbasis_quad, \
                                dbasis_quad, \
                                VolKernelQuatdrature
from automatic_differentiation import ad

from nbspbasis import *         #LB methods

from itertools import product, combinations #quick box plotting...

import copy
from knot_operations import * #knot_removal, knot_insertion
from routines import set_curve_knots

from minDistance import naive_minimizer
import minDistance as mdtfinder


from Quaternion import Quaternion
from frames import Frame, DeepVector
import transformations as tfm
#/*!------------------------------------------------------------------------------------------------------------------------!*/#

#
#**********************************************************************
# An experiment from 2015
class FastDraggableText:
    lock = None  # only one can be animated at a time
    def __init__(self, tbox):
        self.tbox = tbox
        self.press = None
        self.background = None
        ##
        self.savend     = [] ##TLM
        self.savestart  = [] ##TLM

    def connect(self):
        """connect to all the events we need
        """
        self.cidpress = \
            self.tbox.figure.canvas.mpl_connect('button_press_event', 
                                                        self.on_press)
        self.cidrelease = \
            self.tbox.figure.canvas.mpl_connect('button_release_event', 
                                                        self.on_release)
        self.cidmotion = \
            self.tbox.figure.canvas.mpl_connect('motion_notify_event', 
                                                        self.on_motion)

    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.tbox.axes: return
        if FastDraggableText.lock is not None: return
        contains, attrd = self.tbox.contains(event)
        if not contains: return
        print 'event contains', self.tbox.xy
        x0, y0 = self.tbox.xy

        self.press = x0, y0, event.xdata, event.ydata
        FastDraggableText.lock = self

        # draw everything but the selected rectangle and store the pixel buffer
        canvas = self.tbox.figure.canvas
        axes = self.tbox.axes
        self.tbox.set_animated(True)
        canvas.draw()
        self.background = canvas.copy_from_bbox(self.tbox.axes.bbox)

        # now redraw just the rectangle
        axes.draw_artist(self.tbox)

        # and blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_motion(self, event):
        'on motion we will move the tbox if the mouse is over us'
        if FastDraggableText.lock is not self:
            return
        if event.inaxes != self.tbox.axes: return
        x0, y0, xpress, ypress = self.press
        
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        self.savend = [x0+dx,y0+dy]#TLM
        self.tbox.set_x(x0+dx)
        self.tbox.set_y(y0+dy)
        #self.tbox.xy = (dx,dy)

        canvas = self.tbox.figure.canvas
        axes = self.tbox.axes
        # restore the background region
        canvas.restore_region(self.background)
        canvas.restore_region(self.background)

        # redraw just the current rectangle
        axes.draw_artist(self.tbox)
        

        # blit just the redrawn area
        canvas.blit(axes.bbox)

    def on_release(self, event):
        'on release we reset the press data'
        
        if FastDraggableText.lock is not self:
            return

        self.press = None
        FastDraggableText.lock = None

        # turn off the tbox animation property and reset the background
        self.tbox.set_animated(False)
        self.background = None

        # redraw the full figure
        self.tbox.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.tbox.figure.canvas.mpl_disconnect(self.cidpress)
        self.tbox.figure.canvas.mpl_disconnect(self.cidrelease)
        self.tbox.figure.canvas.mpl_disconnect(self.cidmotion)
        
        

#
#**********************************************************************
#    Helper Functions (where should )
def move_vertices(move_vector,vertices):
    length,dim = np.shape(vertices)
    for d in range(dim):
        vertices[:,d] += move_vector[d]
    return vertices
def move_copy_vertices(move_vector,vertices):
    mv = copy.copy(vertices)
    return move_vertices(move_vector,mv)

#
#**********************************************************************
# 
def reverse_spline(curve):
    return Bspline(curve.vertices[::-1], curve.k, curve.nump)
#
#**********************************************************************
# 
def linear_vertices(start, end, num):
    """utility to rturn a set of starting points
        start   : starting point[x,y] 
        end     : ending point[x,y]
        num     : numberof vertices
    """
    start = np.asarray(start)
    end = np.asarray(end)
    D = end - start
    dim = len(D)
    vertices = []
    for i in range(dim):
        vi = np.linspace(start[i],end[i],num)
        vertices.append(vi)
    vertices = np.asarray(vertices)
    vertices = np.transpose(vertices)
    return vertices
    
#
#**********************************************************************
# 
def curvalinear_vertices(start,end,num):
    pi = np.pi
    s = np.linspace(0.,1.,num)
    vertices = []
    for i in s:
        x = np.cos(i*pi)
        y = np.cos(i*pi)
        x1 = start[0]*(1.-i)*x
        y1 = start[1]*(1.-i)*y
        x2 = -end[0]*i*x
        y2 = -end[1]*i*y
        X = x1*(1.-i) + x2*i
        Y = y1*(1.-i) + y2*i
        #XY = start*(1.-i) + end*i
        Y = y1*(1.-i) + y2*i
        vertices.append([X,Y])
    return np.asarray(vertices)
#
#**********************************************************************
# 
class Bspline(object):
    """
    B-Spline Curve

    1.) P&T formulation

        local variables for faster computation
        e.g.
        def someMethod(self):
            U=self.t, etc.
            do sometghing with U, etc...

    Attributes:
        vertices    :   vertices
        order       :   k
        degree      :   p
        knots       :   t

    #Basic NoteS:
            1.) A B-spline curve must have
                a greater number of vertices, N
                than its order, k
                if N==k, then we have a bezier curve.
            2.) The total number of basis functions == number of vertices
                but only some basis functions are necessarily active
        
    #Note: pp. 88 in P&T:
            A B-Spline is at least p-m times continuously differentiable at
            a knot of multiplicity m
            For instance, a O[3] curve with 1 interior knot:
                            p=2 m=1 so C=1
                            tan is 1st derivative
                            FONC(tan) is second derivative - discontinuous at the knot.

    Nov 2017: let's put B-splines in context.
    In a global coordinate context, that is
    """
    def __init__(self, vertices, k, nump=30, t=None, 
                 trange=None, FormParameters=None, 
                 IneqConstraints=None, 
                 surfFormParameters=None):
        """
            Assume knot vector ranges from 0 to 1. 
        """
        self.__version__        = '$Revision: 2.0'.split()[1]
        self.__usage__          = 'usage: Bspline(vertices, order)'
        self.verbose            = False
            
        self.vertices_copy      = np.asarray(vertices)  
        self.vertices           = np.asarray(vertices)  
        self.hierarchical_obsvrs= []
        self.dim                = len(vertices[0])
        if trange is None:
            self.trange         = np.array([0.,1.])  
        else:
            self.trange=trange
        self._t_observers       = []
        if t is None:
            self.t_copy         = make_knot_vector(k,len(vertices)-1,self.trange,type='open')
        else:
            self.t_copy         = t 
        self.t                  = copy.deepcopy(self.t_copy)
        self.n                  = len(vertices)
        self.p                  = k-1
        self.nump               = nump 
        ##
        ##**************************************************
        ## Curve Coordinates
        Origin = np.asarray([0.,0.,0.])
        x, y, z = np.identity(3)
        self.ccs = Frame( np.asarray([x,y,z]) ) #curve coordinate system
        self.origin = DeepVector(v = Origin,
                                B=self.ccs)
        self.rots = Quaternion.from_v_theta(x, 0.) #null rotation vector
        ##
        ##**************************************************
        ##
            
        
        self._k_observers       = []
        self._k_simple_observers= []
        self.k_copy             = k
        self.k                  = k
                     
        self.annotations        = []
        
        self.initialize_curve_basis()
        
        self.area               = None
        self.Mx                 = None
        self.My                 = None
        self.Xc                 = None
        self.Yc                 = None
        self.curvature          = None
        
        # Fairness Functional: Integrals of the Matrix of products of Basis Function Derivatives
        self.initialize_basis_matrices()

        # Use Automatic Differentiation (AD)
        # via the following handles:
        
        self.xpts = None
        self.ypts = None
        self.zpts = None
        
        self.lpts = None #AD lagrange multipliers
        self.epts = None #AD penalty lagrange multipliers - not used
        self.svls = None #AD penalty distance             - not used
        

        # Intrinsic Method Calls : Fully Initiazlize the B-Spline
        self.finish_init()


        #experimental: add constraints (NOT USED - did not go this way):
        self.FormParam              = FormParameters
        self.IneqConstraints        = IneqConstraints
        self.surfFormParam          = surfFormParameters
        
        self.data = {} # and data for optimizations
        

    def __str__(self):
        return 'curve.Bspline k = {}, n = {}'.format(self.k, self.n)
    def __repr__(self):
        return '<class curve.Bspline>'
        #return 'curve.py spline class: \nMethods for Fairness: \n.basis_matrix computes the \nmatrix of integrals of basis product derivatives \nneeded to compute the energy fairness measures \nMomentMatrices : Analogous Matrix for Area and Moments'

    def __call__(self, s=None):
        if s is None:
            self.plotcurve()
            return self.vertices, self.k, self.t
        return self.CurvePoint(s)

    @property
    def k(self):
        return self._k
    @property
    def t(self):
        return self._t
    @property
    def vertices(self):
        return self._vertices
    
    #@property
    #def annotations(self):
    #    return self._annotations
    
    @vertices.setter
    def vertices(self, *args):
        vertices =  args[0]
        self._vertices = vertices
        #if np.any(not self.vertices_copy == vertices): #PEP 8 Jan 17 2018 need to fix this
        if np.any(self.vertices_copy !=vertices):
            if not len(vertices) == self.n:
                #print 'len changed'
                self.n = len(vertices)
                self.t = make_knot_vector(self.k, self.n-1, np.array([0.,1.]), 'open')
                self.initialize_curve_basis()
                self.initialize_basis_matrices()
                self.compute_basis()
            self.compute_curve()
            self.vertices_copy = copy.deepcopy(vertices)
        return
    
    """
    @property
    def frame_vertices(self):
        return self._frame_vertices
    
    @frame_vertices.setter
    def frame_vertices(self, *args):
        frame_vertices =  args[0]
        self._frame_vertices = frame_vertices
        return
    #"""
        
    
    @k.setter
    def k(self, *args):
        k = args[0]
        self._k = k
        if self.k_copy != k:
            #print 'updating k'
            self.k_copy = k
            self.p = k-1
            self.t = make_knot_vector(k, self.n-1, np.array([0.,1.]), 'open')
            self.initialize_curve_basis()
            self.compute_basis()
            self.compute_curve()
            #self.finish_init()
            self.k_copy = k        
        return
        
    @t.setter
    def t(self, *args):        
        t = args[0]
        self._t = t
        #if np.any(self.t_copy != t):
        if not np.any(self.t_copy == t): #PEP 8 change Jan 17 2018
            
            self.initialize_curve_basis()
            self.initialize_basis_matrices()
            self.compute_basis()
            self.compute_curve()
        return
    
    
    def notify_observers(self, *args, **kwargs):
        callback_list = args[0]
        for callback in callback_list:
            print kwargs
            setattr( self, args[1], callback( *args[2]) )
        return
        
    def set_simple_observers(self, *args, **kwargs):
        callback_list = args[0]
        #for callback in callback_list:
        #    print kwargs
        #    setattr(self, args[1], args[2])
        for callback in callback_list:
            print kwargs
            setattr( self, args[1], callback( *args[2]) )
        return
    
    def notify_functions(self, *args, **kwargs):
        callback_list = args[0]
        for callback in callback_list:
            callback()
        return
    
    def bind_to(self, which_list, callback):
        which_list.append(callback)
        return
    
    def info(self, object, spacing=10, collapse=1):
        """Print methods and doc strings.
        
        Takes module, class, list, dictionary, or string."""
        methodList = [method for method in dir(object) if callable(getattr(object, method))]
        processFunc = collapse and (lambda s: " ".join(s.split())) or (lambda s: s)
        print "\n".join(["%s %s" %
                      (method.ljust(spacing),
                       processFunc(str(getattr(object, method).__doc__)))
                     for method in methodList])
    
    def initialize_curve_basis(self):
        self.s              = np.linspace(self.t[0],self.t[-1],self.nump,endpoint=True)  #param space
        self.r              = np.zeros((self.nump,self.dim),float)                              # store B-spline
        self.dr             = np.zeros((self.nump,self.dim),float)                              # store DB-spline
        self.span           = np.zeros((self.nump),int)                                # store non zero span of basis funcs
        self.basis          = np.zeros((self.nump,self.n),float)                        # store all basis func 
        self.dbasis         = np.zeros((self.nump,self.n),float) #1st derivative - old way check - not used any more
        self.allbasis       = np.zeros((self.nump,self.k,self.n),float)             # store all dbasis function values at every param loc. 
        self.allr           = np.zeros((self.nump,self.k,self.dim),float)           # store all curve & curve derivative values at every param
        return
    
    def initialize_basis_matrices(self):
        self.testM1         = np.zeros( (self.n,self.n), float) # write these the hard way - to test M1 computations once and for all
        self.M1             = np.zeros( (self.n,self.n), float) # these are actually used for the integrals - not M itself
        self.M2             = np.zeros( (self.n,self.n), float) # these are actually used for the integrals - not M itself
        self.M3             = np.zeros( (self.n,self.n), float) # these are actually used for the integrals - not M itself
        self.M              = np.zeros( (self.n,self.n,self.nump), float) # Store Non-Differentiated Products of Basis Functions Derivatives.  -ALL!
        self.E1             = None    # fairness functions
        self.E2             = None    # fairness functions
        self.E3             = None    # fairness functions
        self.AL             = None    # arc length
        self.AM1            = np.zeros( (self.n,self.n), float)         # area computations: storing Ibasis products
        self.AM2            = np.zeros( (self.n,self.n,self.n), float)  # moment computations: storing Ibasis products
        return
        
    def finish_init(self):
        # Intrinsic Method Calls : Fully Initiazlize the B-Spline
        self.SpanArray()     # span loc of every pt on the parameterization
        self.BasisArray()    # basis functions only               - store in curve.basis[]
        self.DBasisArray()   #basis functions and all derivatives - store in curve.allbasis[]
        self.allCurveArray() # compute the curve and all existing derivatives at every parameter value
        return
    def compute_basis(self):
        if self.verbose:
            print 'recomputing curve basis'
        self.SpanArray()     # span loc of every pt on the parameterization
        self.BasisArray()    # basis functions only               - store in curve.basis[]
        self.DBasisArray() 
        return
    def compute_curve(self):
        if self.verbose:
            print 'recomputing curve points'
        self.allCurveArray()
        return
    
    
    ##
    ##*****************************************  
    ## 3D transformations
    def rotate_surface(self, axis, angle):
        self.rots = Quaternion.from_v_theta(axis,
                                            angle)
        return
    def translate_surface(self, translation):
        return
    def reflect_surface(self,axis):
        return
    ##
    ##*****************************************  
    ##

    def FindSpan(self,u):
        """FindSpan: a Piegel and Tiller Algorithm
            Determines the knot Span Index
            Needs testing/editing for multiplicative knots.
            B/c it gets stuck when there are multiples
                of the same knot value on the interior
                of the knot vector.
                
                tests to see if the index changes
                 if not, and if the number of iterations exceeds 50, 
                 use whatever index it has gotten to (stuck on)
            Input:
                    n = number of vertices (contrast with P&T, who use # vertices -1)
                    p = degree of curve
                    u = knot position
                    U = knot vector

            Output:
                    mid = the knot span index
        """
        n=self.n
        p=self.p
        U=self.t
        if u==U[n+1]:
            return n-1 #special case n = m-p-1 #TLM mod 6-10-2013

        low  = p
        high = n+1

        mid = (low + high)/2
        count = 0
        max_count = 50
        while ((u < U[mid] or u >= U[mid+1]) and count<max_count):
            if u < U[mid]:
                high = mid
            else:
                low = mid
            mid = (low+high)/2
            count +=1
        if count>10:
            print 'warning, took {} iters to find span, max = {}'.format(count, max_count)
        return mid


    def FindNearestPoint(self, point):
        """
            Input:
                point = [x,y,z]
            Output:
                parameter location s, on the curve
        """
        return mdtfinder.naive_minimizer(self, point)


    def FindParam(self,x,index, tol=0.00000001):
        return self.FindPoint(x,index)    
    
    def FindPoint(self, 
                  x, 
                  index, 
                  tol=0.00000001,
                  maxiter=100):
        """Given an 'x' or 'y' coordinate,
            this algorithm will try to find its'
            coresponding point
            on the curve parameterization
            -------------------------------------------------
            
            Parameters
            ----------
                x       = 1-D point to find on the curve
                index   = 0 for x, 1 for y
                tol     = search tolerance
                
            Output
            ----------
                mid     = s location of x on the curve
                q       = [x,y] values of the curve at mid
                found   = success as True or False"""

        #binary search:
        
        high = self.trange[-1]
        low  = self.trange[0]
        mid  = (low + high)/2. #initial guess
        
        q = self.CurvePoint(mid) #vector valued curve value at parameterization mid
            
        count=0
        while((abs(q[index]-x))>tol and count<maxiter):
            if x<q[index]:
                high = mid
            else:
                low = mid
            mid = (low+high)/2.
            q = self.CurvePoint(mid)
            count += 1
        if count >= maxiter:
            print 'error, point not found'
            found = False
            return mid, q, found
        else:
            found = True
            q = self.CurvePoint(mid)
            return mid, q, found  
        
    
    def FindPointOmnidirectional(self, 
                  findx, 
                  index, 
                  tol=0.00000001,
                  sloc=None):
        """Given an 'x' or 'y' coordinate,
            this algorithm will try to find its'
            coresponding point
            on the curve parameterization
            -------------------------------------------------
            input:
                x       = 1-D point to find on the curve
                index   = {0 to find y using x},
            Output:
                mid     = s location of x on the curve
                q       = [x,y] values of the curve at mid
                
                
                
                
        TODO, TLM OCT 30 2017
            make this work for arbitrary locations
            so you can split the bulb curve.
            
            not used.
                    
                    """

        #binary search:
        
        high = self.trange[-1]
        low  = self.trange[0]
        mid  = (low + high)/2. #initial guess
        
        if sloc is None:
            q = self.CurvePoint(mid) #vector valued curve value at parameterization mid
        else:
            q = self.CurvePoint(sloc)
            
        count=0
        while((abs(q[index]-findx))>tol and count<100):
            if findx<q[index]:
                high = mid
            else:
                low = mid
            mid = (low+high)/2.
            q = self.CurvePoint(mid)
            count += 1
        if count >= 100:
            print 'error, point not found'
            found = False
            return mid, q, found
        else:
            found = True
            q = self.CurvePoint(mid)
            return mid, q, found

    
    

    
    def plot_basis_interval_scaled(self,interval, 
                                   nump=30, offset=0.,
                                   scaler=None):
        p = self.p
        a = interval[0]
        b = interval[1]
        bi = self.active(interval)
        s = np.linspace(a,b,nump,endpoint=True)
        numb = bi[1] - bi[0] + 1
        basis = np.zeros((nump,numb))
        for i,u in enumerate(s):
            span = self.FindSpan(u)
            self.BasisFuns(span,u,basis[i,span-p-bi[0]:span+1])
        #plt.plot(s,)
        plt.plot(s,basis[:]-offset)
        return basis
        
        
    def plot_basis_interval(self,interval, nump=30, offset=0.):
        p = self.p
        a = interval[0]
        b = interval[1]
        bi = self.active(interval)
        s = np.linspace(a,b,nump,endpoint=True)
        numb = bi[1] - bi[0] + 1
        basis = np.zeros((nump,numb))
        for i,u in enumerate(s):
            span = self.FindSpan(u)
            self.BasisFuns(span,u,basis[i,span-p-bi[0]:span+1])
        #plt.plot(s,)
        plt.plot(s,basis[:]-offset)
        return basis
    
    def plot_suppl2_basis_interval(self,interval, nump=30, offset=0.):
        
        curve1 = self
        curve2 = self.dyadic_refinement()
        
        p = curve1.p
        a = interval[0]
        b = interval[1]
        b1 = curve1.active(interval)
        b2 = curve2.active(interval)
        
        tot = int(np.round(nump/(b-a)))
        ma = int(np.round(a*tot) )
        mb = nump
        mc = int( np.round( (1.-b)*tot) )
        
        s1_a = np.linspace(0.,a,ma,endpoint=True)
        s1_b = np.linspace(a,b,mb,endpoint=True)
        s1_c = np.linspace(b,1.,mc,endpoint=True)
        s1 = np.linspace(0.,1.,tot-1,endpoint=True)
        s2 = np.linspace(a,b,nump,endpoint=True)
        
        numb = b2[1] - b2[0] + 1
        basis = np.zeros((nump,numb))
        for i,u in enumerate(s2):
            span = curve1.FindSpan(u)
            curve2.BasisFuns(span,u,basis[i,span-p-bi[0]:span+1])
        #plt.plot(s,)
        
        
        
        
        
        plt.plot(s,basis[:]-offset)
        return basis
    
    
    def eval_BasisPt(self, i, pt):
        N=np.zeros((self.k),float)
        return self.BasisFuns(i,pt,N)
        

    def BasisFuns(self,i,u,N):
        k=self.k
        U=self.t
        """
        P&T module to compute the
        Cox De-Boor Basis Functions
        Efficiently at a point

        Inputs:
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
        return N #array of all non zero basis functions at u
    
    def active_Bik_interval(self,i=None,u=None):
        """Returns the active set of basis functions
        
            Note: this is built to use 
            known span if coming from some algorithm
        
        parameters
        ----------
        i = ith span - location at which the point in question lives
        u = parametric location
        
        returns
        ----------
        [low,high] an closed integer interval of the basis functions
        """
        if u is None:
            assert(i is not None)
            return [i-self.p,i]
        else:
            i = self.FindSpan(u)
            return [i-self.p,i]
            
    def active(self, bounds):
        """Returns the indices of the 
        active non-zero basis functions
        this is the index range of the control vertices
        which have influence here
        """
        if isinstance(bounds,float):
            return self.active_Bik_interval(u=bounds)
        else:
            bi = bounds[0]
            be = bounds[1]
            si = self.FindSpan(bi)
            se = self.FindSpan(be)
            return [si-self.p,se]
            
    #    def active_knots(self,i):
    #        """Given a basis func, aka C.V. index
    #        return the knots over which it is active
    #        """
    #        print 'WARNING: questionable active knots'
    #        print 'where did you get this formula?'
    #        return self.t[i:i+self.p+3]
    
    
            
            
    def TLMBasisFuns(self,u):
        i = self.FindSpan(u)
        k=self.k
        p=self.p
        U=self.t
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
        return Nall #array of all non zero basis functions at u
    
    
            
            
    def BasisLocal(self,u):
        i = self.FindSpan(u)
        k=self.k
        p=self.p
        U=self.t
        N = np.zeros((self.k,1),float)
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
        return N #array of all non zero basis functions at u


    def TLMDBasisFuns(self,u):
        i = self.FindSpan(u)
        #k=self.k
        p=self.p
        n=self.k #for n<=p on page 72.  so not n=self.n.
        U=self.t
        
        DERS = np.zeros((self.k,self.n,1),float) 
        ders = DERS[:,i-p:i+1]
        #N = Nall[i-p:i+1]
        """
        P&T function to compute all the
            Non Vanishing Derivatives of 
            Cox De-Boor Basis Functions
            Efficiently at a point

        Inputs:
            i = span index
            u = point on the parameterization
            k = order of the basis function (passed for ease of code update)
            U = knot vector
            ders = EmptyArray of Basis functions from [0...p], i.e. size k, at u.
            
        Output:
            ders = Array of all Basis Function Values at u.   """

        left=np.zeros((self.k),float)
        right=np.zeros((self.k),float)
        
        ndu = np.zeros((self.k,self.k),float) #store basis funcs and knot differences
        a   = np.zeros((2,self.k),float) #store two most recently computed fows of ak,j ak-1,j
        
        ndu[0,0]=1.0
        for j in range(1,p+1):
            left[j]  = u-U[i+1-j]
            right[j] = U[i+j]-u
            saved = 0.0
            for r in range(j):
                                                    #upper tri
                ndu[j,r]    = right[r+1]+left[j-r]
                temp        = ndu[r,j-1]/ndu[j,r]
                                                    #lower tri
                ndu[r,j]    = saved+right[r+1]*temp
                saved = left[j-r]*temp
            ndu[j,j] = saved
            
        for j in range(p+1): #LOAD THE BASIS FUNCTIONS
            ders[0,j]=ndu[j,p]
            
        for r in range(p+1):
            s1=0
            s2=1
            a[0,0]=1.0
            for k in range(1,n):
                d=0.0
                rk = r-k
                pk = p-k
                if r>=k:
                    a[s2,0] = a[s1,0]/ndu[pk+1,rk]
                    d       = a[s2,0]*ndu[rk,pk]
                if rk>=-1:
                    j1=1
                else:
                    j1=-rk
                if ((r-1)<=pk):
                    j2=k-1
                else:
                    j2=p-r
                for j in range(j1,j2+1):
                    a[s2,j]=(a[s1,j]-a[s1,j-1])/ndu[pk+1,rk+j]
                    d += a[s2,j]*ndu[rk+j,pk]
                if r<=pk:
                    a[s2,k]=-a[s1,k-1]/ndu[pk+1,r]
                    d+=a[s2,k]*ndu[r,pk]
                ders[k,r] = d
                j=s1
                s1=s2
                s2=j
        r=p
        for k in range(1,n):
            for j in range(p+1):
                ders[k,j] *= r
            r*=(p-k)
        return DERS#ders


    def DersBasisFunc(self,i,u,ders):
        """
        P&T function to compute all the
            Non Vanishing Derivatives of 
            Cox De-Boor Basis Functions
            Efficiently at a point

        Inputs:
            i = span index
            u = point on the parameterization
            k = order of the basis function (passed for ease of code update)
            U = knot vector
            ders = EmptyArray of Basis functions from [0...p], i.e. size k, at u.
            
        Output:
            ders = Array of all Basis Function Values at u.   """
        
        #k=self.k
        p=self.p
        n=self.k #for n<=p on page 72.  so not n=self.n.
        U=self.t

        left=np.zeros((self.k),float)
        right=np.zeros((self.k),float)
        
        ndu = np.zeros((self.k,self.k),float) #store basis funcs and knot differences
        a   = np.zeros((2,self.k),float) #store two most recently computed fows of ak,j ak-1,j
        
        ndu[0,0]=1.0
        for j in range(1,p+1):
            left[j]  = u-U[i+1-j]
            right[j] = U[i+j]-u
            saved = 0.0
            for r in range(j):
                                                    #upper tri
                ndu[j,r]    = right[r+1]+left[j-r]
                temp        = ndu[r,j-1]/ndu[j,r]
                                                    #lower tri
                ndu[r,j]    = saved+right[r+1]*temp
                saved = left[j-r]*temp
            ndu[j,j] = saved
            
        for j in range(p+1): #LOAD THE BASIS FUNCTIONS
            ders[0,j]=ndu[j,p]
            
        for r in range(p+1):
            s1=0
            s2=1
            a[0,0]=1.0
            for k in range(1,n):
                d=0.0
                rk = r-k
                pk = p-k
                if r>=k:
                    a[s2,0] = a[s1,0]/ndu[pk+1,rk]
                    d       = a[s2,0]*ndu[rk,pk]
                if rk>=-1:
                    j1=1
                else:
                    j1=-rk
                if ((r-1)<=pk):
                    j2=k-1
                else:
                    j2=p-r
                for j in range(j1,j2+1):
                    a[s2,j]=(a[s1,j]-a[s1,j-1])/ndu[pk+1,rk+j]
                    d += a[s2,j]*ndu[rk+j,pk]
                if r<=pk:
                    a[s2,k]=-a[s1,k-1]/ndu[pk+1,r]
                    d+=a[s2,k]*ndu[r,pk]
                ders[k,r] = d
                j=s1
                s1=s2
                s2=j
        r=p
        for k in range(1,n):
            for j in range(p+1):
                ders[k,j] *= r
            r*=(p-k)
        return #ders

    def DersOneBasisFunc(curve,i,u,ders):
        """
        Compute derivatives of basis function Nip
            Input:
                p = order of the curve
                m = len(knot vector)-1
                U = knot vector
                i = span index at u - i.e. there are n of them, k=order nonzero at any point.
                u = parameter space position 
                n = number of derivatives
            Output:
                ders[k,k] = the ith set of basis function derivative values at u"""
        p=curve.p
        #m=len(curve.t)-1
        U=curve.t
        n=curve.k

        N = np.zeros((curve.k,curve.k),float)
        ND = np.zeros((curve.k),float)
        
        if (u<U[i] or u >= U[i+p+1]):  #local property
            for k in range(n):
                ders[k]=0.0
            return
        for j in range(p+1): #initialize zeroth degree functions
            if (u>=U[i+j] and u<U[i+j+1]):
                N[j,0]=1.0
            else:
                N[j,0]=0.0
        for k in range(1,p+1): # create full tringular table
            if N[0,k-1]==0.0:
                saved = 0.0
            else:
                saved = ((u-U[i])*N[0,k-1])/(U[i+k]-U[i])
            for j in range(p-k+1):
                Uleft = U[i+j+1]
                Uright = U[i+j+k+1]
                if N[j+1,k-1]==0.0:
                    N[j,k]=saved
                    saved=0.0
                else:
                    temp = N[j+1,k-1]/(Uright-Uleft)
                    N[j,k] = saved + (Uright-u)*temp
                    saved = (u-Uleft)*temp
        print N
        
        ders[0] = N[0,p]
        for k in range(1,n): # compute the derivatives
            for j in range(k+1):
                ND[j]=N[j,p-k]    
            for jj in range(1,k+1):
                if ND[0]==0.:
                    saved=0.0
                else:
                    saved = ND[0]/(U[i+p-k+jj]-U[i])
                    
                for j in range(k-jj+1):
                    Uleft = U[i+j+1]
                    #print 'j = {}, jj = {}, k= {}'.format(j,jj,k)
                    #print i+j+p+jj+1
                    #print U[i+j+p+jj+1]
                    Uright = U[i+j+p+jj+1]
                    if (ND[j+1]==0.0):
                        ND[j]=(p-k+jj)*saved
                        saved = 0.0
                    else:
                        temp = ND[j+1]/(Uright-Uleft)
                        ND[j] =(p-k+jj)*(saved-temp)
                        saved = temp
                    
            ders[k] = ND[0] # the kth derivative
            
        return

    

    def SpanArray(curve):
        """
            Return the array of all span indices
            for the standard parameterization"""
        for i,u in enumerate(curve.s):
            curve.span[i]=curve.FindSpan(u)
        return

    
    def BasisArray(curve):
        p=curve.p
        """
            Return the array of all basis functions arrays
            for the standard parameterization"""
        for i,u in enumerate(curve.s):
            curve.BasisFuns(curve.span[i],u,curve.basis[i,curve.span[i]-p:curve.span[i]+1])

        return



    def DBasisArray(curve):
        p=curve.p
        """
            Return the array of all 1st derivative basis functions arrays
            for the standard parameterization"""
        for i,u in enumerate(curve.s):
            curve.DersBasisFunc(curve.span[i],u,curve.allbasis[i,:,curve.span[i]-p:curve.span[i]+1])
        return


    #def CurvePoint(self,u,N):  # P&T algorithm
    def CurvePoint(self,u):     # python algorithm
        """
            P&T algorithm to find a point on a curve
            Inputs:
                n = number of vertices (P&T use # vertices -1)
                i = span index
                u = point on the parameterization
                k = order of the basis function (passed for ease of code update)
                U = knot vector
                N = EmptyArray of Basis functions from [0...k], i.e. size k, at u. (deprecated - not used this way anymore)
            
            Output:
                dot(N,P) = = 2D or 3D point on a B-Spline curve 
        """
        #n=self.n
        #k=self.k
        p=self.p
        #U=self.t
        P=self.vertices

        N=np.zeros((self.k),float)
        
        span = self.FindSpan(u)
        self.BasisFuns(span,u,N)
        #print N
        #C=0.0
        return np.dot(N,P[span-p:span+1,:])

    #def CurvePoint2(self,u,N):  # algorithm places min basis values in correct location within the "N" array
    def CurvePoint2(self,u):     # python algorithm
        """
            P&T algorithm to find a point on a curve
            Inputs:
                n = number of vertices (P&T use # vertices -1)
                i = span index
                u = point on the parameterization
                k = order of the basis function (passed for ease of code update)
                U = knot vector
                N = zeros array of Basis functions of size n. (deprecated)
            
            Output:
                dot(N,P) = 2D or 3D point on a B-Spline curve 
        """
        
        #n=self.n
        #k=self.k
        p=self.p
        #U=self.t
        P=self.vertices

        N=np.zeros((self.n),float)
        
        span = self.FindSpan(u)
        self.BasisFuns(span,u,N[span-p:span+1])
        #print N

        return np.dot(N,P)  # same as curvePoint
    
    def CurveDerivsAlg1(curve,u,CK):
        """P&T algorithm to find a point on a curve
            and all non-zero derivatives of the curve
            at that point

            Inputs:
                n  = number of vertices (P&T use # vertices -1) (this is standard from page 93 forward.)
                p  = order of the curve
                U  = knot vector
                P  = curve vertices
                u  = point on the parameterization
                d  = k - the number of nonvanishing derivatives
                CK = zeros array of Basis functions of size n.
            
            Output:
                CK = array of 2D or 3D points on a B-Spline curve or its derivatives"""
        #n = curve.n
        #p = curve.p
        #U = curve.t
        P = curve.vertices
        d = curve.k #a 4th order curve has no more than 3 non-zero derivatives.
        
        #du = min(d,p)
        #for k in range(p+1,d):
        #    CK[k]=0.0
        nders=np.zeros((curve.k,curve.n),float)
        span=curve.FindSpan(u)
        curve.DersBasisFunc(span,u,nders)
        for k in range(d):
            #CK[k]=0.0
            #for j in range(k):
            #CK[k]=np.dot(CK[k,:],P[span-p:span+1,:])
            #print CK[k]
            #print P
            CK[k]=np.dot(nders[k],P)
        
        return



    def allCurveArray(curve):
        #n=curve.n
        k=curve.k
        #p=curve.p
        #U=curve.t
        P=curve.vertices
        """
            Return the array of all curve values
            for the standard parameterization"""
        curve.allr[:,:,:]=0.0
        for d in range(curve.dim): # loop over all spatial dimensions
            for j in range(k): #loop over all existing derivatives
                for i in range(curve.n): # loop over all points
                    #implicit loop over all curve parameterization points:
                    curve.allr[:,j,d]+=curve.allbasis[:,j,i]*P[i,d]

        curve.r = curve.allr[:,0,:]
        curve.dr = curve.allr[:,1,:]
        
        return

    def knot_insertion(self, knot):
        t, vertices         = knot_insertion(self.k,
                                             self.vertices,
                                             self.t,
                                             knot)
        #self.n              = len(vertices)  
        self.vertices       = vertices
        self.t              = t    
        #self.initialize_curve_basis()
        #self.initialize_basis_matrices()
        #self.finish_init()
        return
    
    def knot_removal(self, index,
                     tol=.00001):
        """
        Remove a knot and move the control points
        to make the curve stay the same
        
        TODO: this is unsafe and unfinished
        it assumes the knot at index
        will be removed successfully once.
        """
        #if self.dim ==2:
        print 'Knot Removal Risky - Not finished code (but almost)'
        #        t,vertices = knot_removal(self.k-1,
        #                                  self.vertices,
        #                                  self.t,
        #                                  self.t[index],
        #                                  index)
        knot = self.t[index]
        #t,vertices = unsafe_knot_removal(self,knot,index)
        vertices = RemoveCurveKnot(self,
                                   u=knot,
                                   r=index,
                                   TOL=tol)
        
        #now we have the vertices
        #
        #
        #***********************************************
        #drop the removed vertex
        pws = list(np.shape(vertices))
        pws[0] = pws[0]-1
        Pw_final = np.zeros(pws,float)
        final_index = 0
        for i,vtx in enumerate(self.vertices):
            if i is not index-self.p+1:
                Pw_final[final_index] = vertices[i]
                final_index +=1
        #
        #***********************************************
        # drop the removed knot
        ti = self.t[0:index]
        te = self.t[index+1:]
        tv = np.zeros((len(ti)+len(te)),float)
        tv[:index] = ti
        tv[index:] = te
        self.vertices = Pw_final
        self.t = tv
        return
        #        else:
        #            self.knot_removal(index)
        #            return 
        
    def knot_removal3D(self,index):
        knot = self.t[index]
        knots, vertices = RemoveCurveKnot(self,
                                       u=knot,
                                       r=index)
        ti = self.t[0:index]
        te = self.t[index+1:]
        tv = np.zeros((len(ti)+len(te)),float)
        tv[:index] = ti
        tv[index:] = te
        self.vertices = vertices
        self.t = tv
        return
    
    def plot_knots(self, y=0):
        plt.plot(self.t,y*np.ones(len(self.t)),marker = 'o')
        plt.axis('equal')        
        return

    def CurveSplit(self, pt=.5):
        """ Routine to split a curve at a specified point
            into 2 curves ..
            
            ?
            of order equal 
            to the original curve order.? 
            - no - maybe not since the curves have different knots

            pt = paramatric position == a knot position...

            Note:
                takes either standard or interpolated splines
                
            Note 2: must improve by:
                    -any knot pt to be chosen!
            
            Note 3: make another routine to:
                    -curves come out with identical knot vectors
                      to the self curve! -this is crucial...
        """

        start = self.t[0]
        end   = self.t[-1]
        cc=copy.deepcopy(self)

        # 0.) start composing the new knot vectors:
        t1=list(cc.t[0:cc.k])+list(cc.t[(len(cc.t)-cc.k):len(cc.t)])
        t2=list(cc.t[0:cc.k])+list(cc.t[(len(cc.t)-cc.k):len(cc.t)])

        
        # 1.) determine interior knot set
        interiorKnots = sorted(interiorKnotsSet(cc.t))

        
        # 2.) determine the knot multiplicity at pt
        mult = multiplicity(cc.t)


        # 3.) make a list of all interior knots on each segment
        c1l = []
        c2l = []
        match_knotvector = []
        for knot in interiorKnots:
            if knot<pt:
                for i in range(mult[knot]):
                    c1l.append(knot)
            if knot>pt:
                for i in range(mult[knot]):
                    c2l.append(knot)
            elif knot==pt:  # In a sense, this is not an interior knot for these vertices! - be careful.
                for i in range(mult[knot]):
                    match_knotvector.append(knot) 
                


        # 4.)  break the curve:
        ##       insert knots at the point "pt" (k-mult[pt]) times.
        Flag = 0
        try:
            numberTimes = cc.k-mult[pt]  # 1 or moreknots exist at point = pt
            #Flag = 1
        except:
            numberTimes = cc.k
        for i in range(numberTimes):
            match_knotvector.append(pt)
            #cc.t, cc.vertices = knot_insertion(cc.k,cc.vertices,cc.t,pt)
            cc.knot_insertion(pt)
            #cc.n+=1
            
            
        # 5.) get the new knot vector multiplicity
        multcc = multiplicity(cc.t)  

        #match_knotvector = []
        #for knot in interiorKnots:
        #    if knot==pt:  # In a sense, this is not an interior knot for these vertices! - be careful.
        #        match_knotvector.append(knot)


        # Generalize the code below:
        """
        if Flag==1:  #then 1 or more knots exist at the split point
            v1=cc.vertices[0:cc.k,:]
            v2=cc.vertices[(cc.k):(2*cc.k),:]
            c1 = Bspline(v1, cc.k)
            c2 = Bspline(v2, cc.k)
            for knot in match_knotvector:
                c1.t,c1.vertices = knot_insertion(c1.k,c1.vertices,c1.t,pt)
                c2.t,c2.vertices = knot_insertion(c2.k,c2.vertices,c2.t,pt)
        #"""
        
            
        #if Flag==0:

        v1=cc.vertices[0:(cc.k+len(c1l)),:]
        v2=cc.vertices[(cc.k+len(c1l)):(2*cc.k+len(c1l)+len(c2l)),:]

        #length of the new knot vectors:
        len_t1 = cc.k+len(v1)
        len_t2 = cc.k+len(v2)

        dt = pt-start
        for i,knot in enumerate(c1l):
            t1.insert(cc.k+i,(knot-start)/dt)

        dt = end-pt
        for i,knot in enumerate(c2l):
            t2.insert(cc.k+i,(knot-pt)/dt)

        t1=np.asarray(t1)
        t2=np.asarray(t2)
        c1 = Bspline(v1, cc.k, 50, t1 )
        c2 = Bspline(v2, cc.k, 50, t2 )
        
        #if Flag==1:  #then 1 or more knots exist at the split point
        #    for knot in match_knotvector:
        #        c1.t,c1.vertices = knot_insertion(c1.k,c1.vertices,c1.t,pt)
        #        c2.t,c2.vertices = knot_insertion(c2.k,c2.vertices,c2.t,pt)
        #c1.n = len(c1.vertices) 
        #c2.n = len(c2.vertices) 
                
        return c1,c2
    
    def curve_join(curve1, curve2):
        """
            Assumes coincident end-start point:
            
            Joins the end of curve1
            to the start of curve2
            parameterwise
            
            knot_insertion(k,points,t,ta)
            k=order
            points = control vertices
            t = knot vector
            ta = where to insert the knot on s
        """
        
        for i in range(curve1.k-1):
            curve1.t, curve1.vertices = knot_insertion(curve1.k,curve1.vertices,curve1.t,1.)
            curve1 = Bspline(curve1.vertices, curve1.k, curve1.nump)
            ## let curve2 knot itself!
            dummy, curve2.vertices = knot_insertion(curve2.k,curve2.vertices,curve2.t,0.)
            curve2 = Bspline(curve2.vertices, curve2.k, curve2.nump)
        new_vertices = np.zeros((curve1.n+curve2.n, max(curve1.dim,curve2.dim)),float)
        ## error: the curves already fit:  the knot vectors have to be scaled and the curves combined
        
        
        ## the only thing to be done here is to simplify the curve as much as possible!
        return


    def basis_matrix(curve):
        """
        Method to Compute the M matrix of
        integrated products of
        basis function derivatives

        Attributes:

            N   = order of accuracy of the integration
            i   = the ith vertex
            k   = order of the curve            
            t   = knot vector
            d   = derivative order

        Note:
            This method is integral over the curve
            and so returns one value for each M1,M2,M3
            for each curve
            
        Physically, 
            basis function derivatives let 
            you contruct derivatives of the 
            curve (or higher dimenionsal analog)
            as functions of precomputed basis 
            derivatives
            -derivatives of basis funtions 
            are basis funtions of one less degree
            -this makes for nice 
            differentiation-by-projection
            operators, concerning derivatives of 
            functions represented on such 
            spaces.
            
            This is very different in character from 
            optimizing the design variables.
            Do not get them confused!
        """
        ## change to pass pre-evaluated dNik
        curve.d1    = 1
        curve.d2    = 2
        curve.d3    = 3
        #Calculate the order of the multiplied derivative basis functions:
        curve.N1    = (curve.p-curve.d1)+(curve.p-curve.d1)+1  #Order of basis product = e.g. (degree=3-1-1=1)+(3-1-1=1)+1=3
        curve.N2    = (curve.p-curve.d2)+(curve.p-curve.d2)+1  # Only degree exponents add!
        curve.N3    = (curve.p-curve.d3)+(curve.p-curve.d3)+1  # use order - therefore, after adding degrees, add one extra to get the order


        ## TLM mod - needs verification 4-17-2013 - ammended to zero, 7-7-2013
        v1=0
        v2=0
        v3=0
        

        curve.args1 = curve.k,curve.t,curve.d1  # note: passing order to args is inconsistent with implementation... right?
        curve.args2 = curve.k,curve.t,curve.d2
        curve.args3 = curve.k,curve.t,curve.d3

        #"""
        if (curve.k-curve.d3>v3):
            for i in range(curve.n):
                for j in range(curve.n):
                    curve.M1[i,j] = dbasis_quad(curve,curve.N1, i, j, *curve.args1) #integrate
                    curve.M2[i,j] = dbasis_quad(curve,curve.N2, i, j, *curve.args2) # each basis product
                    curve.M3[i,j] = dbasis_quad(curve,curve.N3, i, j, *curve.args3) # of derivatives
        elif (curve.k-curve.d2>v2):
            for i in range(curve.n):
                for j in range(curve.n):
                    curve.M1[i,j] = dbasis_quad(curve,curve.N1, i, j, *curve.args1) # of basis func.
                    curve.M2[i,j] = dbasis_quad(curve,curve.N2, i, j, *curve.args2) # individualy
        ##"""
        elif (curve.k-curve.d1>v1):
            for i in range(curve.n):
                for j in range(curve.n):
                    curve.M1[i,j] = dbasis_quad(curve,curve.N1, i, j, *curve.args1)
            
        return

    
    def greville_abscissa(self,i):
        """
            takes   : the ith knot
            returns : the knot average 
                      corresponding control point c_{i}
                      of the
                      - ith knot
                      - jth curve
        """
        d = self.p
        ga = 0.
        #for t in range(i,i+d): #farin definition
        for t in range(i+1,i+d+1):
            ga += self.t[t]
        return ga/d
    

    def MomentMatrices(curve):
        """
            Compute the AM1 - Area
            and
            Compute the AM2 - Moment
            matrices

            These are integrated products of Derivatives of B-Spline Basis functions.
        """
        if curve.verbose:
            print 'computing AM1 and AM2'
        #Compute the Area Matrices
        #of the BSpline as per usual:
        N1 = curve.k+curve.k-1
        #t1 = make_knot_vector((N1),curve.n-1,curve.trange,type='open')
        #curve.args1 = curve.k,t1,1
        curve.args1 = curve.k,curve.t,1
        
        #N2 = curve.k+curve.k+curve.p
        #t2 = make_knot_vector((N2),curve.n-1,curve.trange,type='open')
        #curve.args2 = curve.k,t2,1 
        
        curve.AM1[:,:] = 0.0
        for i in range(curve.n):            # loop over  basis values
            for j in range(curve.n):        # loop over dbasis values
                    curve.AM1[i,j] = ibasis_quad2(curve, N1, i, j, *curve.args1) #integration

        curve.AM1=curve.AM1-np.transpose(curve.AM1) #Stephan Harries Formulation, page 158


        # Compute the momment matrices
       
        N2 = curve.k+curve.k+curve.k-1 #written out explicitly
        t2 = make_knot_vector((N2),curve.n-1,curve.trange,type='open')
        #curve.args2 = curve.k,t2,1 
        curve.args2 = curve.k,curve.t,1
        
        for i in range(curve.n):            # loop over  basis values
            for j in range(curve.n):        # loop over dbasis values
                for l in range(curve.n):
                    curve.AM2[l,i,j] = Mbasis_quad(curve, N2, i, j, l, *curve.args2)
                    #note: take advantage of the way Numpy calculates 3D arrays to avoid a loop structure below! [l_first!,then i,j]
                    # note: j is the derivative index...
   



    def fairness(curve):
        """
            Method to compute the non-weighted fairness 
            functionals of a B-spline curve.

        """        
        curve.E1 = np.dot(np.dot(curve.vertices[:,0],curve.M1),curve.vertices[:,0])+np.dot(np.dot(curve.vertices[:,1],curve.M1),curve.vertices[:,1])
        curve.E2 = np.dot(np.dot(curve.vertices[:,0],curve.M2),curve.vertices[:,0])+np.dot(np.dot(curve.vertices[:,1],curve.M2),curve.vertices[:,1])
        curve.E3 = np.dot(np.dot(curve.vertices[:,0],curve.M3),curve.vertices[:,0])+np.dot(np.dot(curve.vertices[:,1],curve.M3),curve.vertices[:,1])

        return


    def precompute_arc_length_basis_matrix_integrals(self):
        self.pts_M_pts()
        return
    def pts_M_pts(curve):
        """
        Method to compute the "M" (non integrated)
        Matrix of products
        of basis function derivatives
        at a point s on curve parameterized
        by the knot vector, t.
        
        These are used for the arc length calculation I believe (Spet. 2014 - its been a while)
            
        Input:  B-spline curve
            Output: nxn matrix

        Note:
            This method is integral over the curve
            and so returns one value for each M1,M2,M3
            for each curve
        """
        
        #curve.args1 = curve.k,curve.t,curve.d1
        for ts in range(1,len(curve.s)-1,1):
            for i in range(curve.n):
                for j in range(curve.n):
                    curve.M[i,j,ts] = curve.allbasis[ts][1,i] * curve.allbasis[ts][1,j] #1st derivatives of the basi at param ts
        return
    
    def lift_1D(self,vertex,index):
        """increase the dimension of 
        the curve by 1 with the new dimension
        located at index position index
        """
        vv = [0,1,2]
        vv.pop(index)
        olv = self.vertices
        nlv = np.zeros((self.n,self.dim+1),float)
        nlv[:,vv] = olv[:]
        nlv[:,index] = vertex[index]
        self.vertices = nlv
        return
    
    
        
    def compute_arclength(self):
        xpts = self.vertices[:,0]
        ypts = self.vertices[:,1]
        presum = 0.0
        ssum   = 0.0
        for ts in range(1,len(self.s)-1,1):
            M = self.M[:,:,ts] #depends on self.pts_M_pts()
            presum   =  np.dot((np.dot(xpts,M)),xpts) # x part
            presum   =  presum + np.dot((np.dot(ypts,M)),ypts) # y part
            #print presum.value
            #presum   =  presum.sqrt()
            presum = np.sqrt(presum)
                        #presum is the scalar valued function we must integrate over the local parameter value:
            a        = self.s[ts-1]   #min 
            b        = self.s[ts]     #max
            ssum     = presum*(b-a) + ssum   # very simple integration
            presum   =  0.0
        self.AL=ssum
        return
        
    def compute_arclengthApprox(self):
        vertices = self.vertices
        cheap_sum = 0.
        for i in range(1,self.n):
            test = np.linalg.norm(vertices[i]-vertices[i-1])
            cheap_sum = cheap_sum + test
        self.ALapprox = cheap_sum
        #self.ALapprox = np.linalg.norm(vertices)
        return

    def compute_tangent(self, s, vertices=None):
        return self.curvetan(s,vertices)
    
    def curvetan(curve,s, vertices=None):
        if vertices is None:
            xpts = curve.vertices[:,0]
            ypts = curve.vertices[:,1]
        else:
            xpts = vertices[:,0]
            ypts = vertices[:,1]
        p=curve.p
    
        # Get the Dbasis functions at s:
        localBasis = np.zeros((curve.n,curve.n),float)
        span = curve.FindSpan(s)
        curve.DersBasisFunc(span,s,localBasis[:,span-p:span+1])
        #basis=[0.1,0.2,0.4,0.2,0.1] #test basis
        
        #Using standard form:
        xpts = curve.vertices[:,0]
        ypts = curve.vertices[:,1]

        qx = np.dot(xpts,localBasis[1])
        qy = np.dot(ypts,localBasis[1])


        #test = (qy/qx).atan()
        #test = np.arctan(qy/qx)
        test = np.arctan2(qy,qx)

        return test
    
    def curvetanZ(curve,s, vertices=None):
        if vertices is None:
            xpts = curve.vertices[:,0]
            zpts = curve.vertices[:,2]
        else:
            xpts = vertices[:,0]
            zpts = vertices[:,2]
        p=curve.p
    
        # Get the Dbasis functions at s:
        localBasis = np.zeros((curve.n,curve.n),float)
        span = curve.FindSpan(s)
        curve.DersBasisFunc(span,s,localBasis[:,span-p:span+1])

        qx = np.dot(xpts,localBasis[1])
        qz = np.dot(zpts,localBasis[1])


        #test = (qy/qx).atan()
        test = np.arctan(qz/qx)

        return test

    def compute_curvature(self, s, vertices=None):
        if vertices is None:
            xpts = self.vertices[:,0]
            ypts = self.vertices[:,1]
        else:
            xpts = vertices[:,0]
            ypts = vertices[:,1]
        localBasis = np.zeros((self.n,self.n),float)
        span = self.FindSpan(s)
        self.DersBasisFunc(span,s,localBasis[:,span-self.p:span+1])
        localBasis = localBasis
        qxdot   = np.dot(xpts,localBasis[1])
        qxddot  = np.dot(xpts,localBasis[2])
        qydot   = np.dot(ypts,localBasis[1])
        qyddot  = np.dot(ypts,localBasis[2])
        store = (qxdot*qyddot - qydot*qxddot)
        curvature = store/(np.sqrt(qxdot*qxdot + qydot*qydot)**3.)
        return curvature

    def curvearea(curve):
        """
            Automatic Differentiation
            of the Area
            Under a B-Spline Curve,
            in x-y space
            
            Use Harries formulation p.158 
        """
        xpts = curve.vertices[:,0]
        ypts = curve.vertices[:,1]
        if abs(np.linalg.norm(curve.AM1)) < .0000001:
            curve.MomentMatrices()
        #temp1   = np.dot(np.dot(xpts,np.transpose(curve.AM1)),ypts)
        temp1   = np.dot(ypts,np.dot(curve.AM1,xpts)) #cleaner
        store2  = (xpts[-1]*ypts[-1] - xpts[0]*ypts[0])*0.5
        curve.area = temp1*.5 + store2
        return
    
    def compute_area(self, x_axis = 0., vertices=None):
        if vertices is None:
            xpts = self.vertices[:,0]
            #ypts = self.vertices[:,1]
            ypts = self.offset_vertices(
                            axis = x_axis, which = 1)[:,1].T
        else:
            xpts = vertices[:,0]
            ypts = vertices[:,1]
        if abs(np.linalg.norm(self.AM1)) < .0000001:
            self.MomentMatrices()
        temp1   = np.dot(np.dot(xpts,np.transpose(self.AM1)),ypts)
        store2  = (xpts[-1]*ypts[-1] - xpts[0]*ypts[0])
        self.area = (temp1 + store2)*0.5
        return
    
    def compute_area_to_y(self, x_axis = 0.):
        xpts = self.offset_vertices(axis = x_axis, which = 0)[:,0].T
        #xpts = self.vertices[:,0]
        ypts = self.vertices[:,1]
        if abs(np.linalg.norm(self.AM1)) < .0000001:
            self.MomentMatrices()
        #temp1   = np.dot(np.dot(xpts,(-self.AM1).T),ypts)
        temp1   = -np.dot(ypts,np.dot(self.AM1,xpts))
        store2  = -(xpts[0]*ypts[0] - xpts[-1]*ypts[-1])
        self.area_to_y = (temp1 + store2)*0.5
        return
    
    def compute_moments(curve, vertices=None):
        if vertices is None:
            xpts = curve.vertices[:,0]
            ypts = curve.vertices[:,1]
        else:
            xpts = vertices[:,0]
            ypts = vertices[:,1]
        n = curve.n

        temp_array0=[]
        temp_array1=[]
        for l in range(n):
            temp_array0.append(list())
            temp_array1.append(list())
            for i in range(n):
                temp_array0[l].append( np.dot(xpts,np.transpose(curve.AM2[l,i,:])) )
                temp_array1[l].append( np.dot(ypts,np.transpose(curve.AM2[l,i,:])) )
    
        temp_array0=np.asarray(temp_array0)
        temp_array1=np.asarray(temp_array1)
        
        temp3  = np.dot(ypts,np.dot(ypts,temp_array0))
        temp4  = np.dot(ypts,np.dot(xpts,temp_array1))
        store3 = (xpts[-1]*ypts[-1]*ypts[-1])-(xpts[0]*ypts[0]*ypts[0])
        curve.Mx  = (store3+((temp3-temp4)*2.))*(1./6)
        
        temp5  = np.dot(xpts,np.dot(ypts,temp_array0))
        temp6  = np.dot(xpts,np.dot(xpts,temp_array1))
        store4 = (ypts[-1]*xpts[-1]*xpts[-1]) - (ypts[0]*xpts[0]*xpts[0])
        curve.My   = (store4 + (temp5-temp6))*(1./3.)
        return
        
    def computeCentroid(curve):
        curve.computeXc()
        curve.computeYc()
        return

    def computeXc(curve):
        curve.Xc = curve.My/curve.area
        return 

    def computeYc(curve):
        curve.Yc = curve.Mx/curve.area
        return
        
    def compute_moments_to_y(curve):
        xpts = curve.vertices[:,0]
        ypts = curve.vertices[:,1]
        n           = curve.n

        temp_array0=[]
        temp_array1=[]
        for l in range(n):
            temp_array0.append(list())
            temp_array1.append(list())
            for i in range(n):
                temp_array0[l].append( np.dot(xpts,np.transpose(curve.AM2[l,i,:])) )
                temp_array1[l].append( np.dot(ypts,np.transpose(curve.AM2[l,i,:])) )
    
        temp_array0=np.asarray(temp_array0)
        temp_array1=np.asarray(temp_array1)
        
        temp3  = np.dot(ypts,np.dot(ypts,temp_array0))
        temp4  = np.dot(ypts,np.dot(xpts,temp_array1))
        store3 = -( (xpts[-1]*ypts[-1]*ypts[-1])-(xpts[0]*ypts[0]*ypts[0]) )
        curve.My_to_y  = -(store3+((temp3-temp4)*2.))*(1./6)
        
        temp5  = np.dot(xpts,np.dot(ypts,temp_array0))
        temp6  = np.dot(xpts,np.dot(xpts,temp_array1))
        store4 = -( (ypts[-1]*xpts[-1]*xpts[-1]) - (ypts[0]*xpts[0]*xpts[0]) )
        curve.Mx_to_y   = -(store4 + (temp5-temp6))*(1./3.)
        return
        
    def computeXc_to_y(curve):
        curve.Xc_to_y = curve.My_to_y/curve.area_to_y
        return 

    def computeYc_to_y(curve):
        curve.Yc_to_y = curve.Mx_to_y/curve.area_to_y
        return
    
    
    def extremes1d(self, index):
        dmin = self.vertices[0,index]
        dmax = self.vertices[0,index]
        for pt in self.vertices:
            dmin = min(dmin,pt[index])
            dmax = max(dmax,pt[index])
        return (dmin,dmax)
    
    def smart_extremes(self):
        extremes = []
        for di in range(self.dim):
            extremes.append(self.extremes1d(index=di) )
        return extremes
        
    def extreme_C0(self):
        xmin = self.vertices[0,0]
        xmax = self.vertices[0,0]
        ymin = self.vertices[0,1]
        ymax = self.vertices[0,1]
        for el in self.vertices:
            xmin = min(xmin, el[0])
            xmax = max(xmax, el[0])
            ymin = min(ymin, el[1])
            ymax = max(ymax, el[1])
        return xmin,xmax,ymin,ymax

    def offset_vertices(self, axis = 0., which = 0):
        verts = copy.deepcopy(self.vertices)
        verts[:,which] = (verts[:,which] - axis)
        return verts
    
    def offset_vertices_for_area(self, axis = 0., which = 0):
        verts = copy.deepcopy(self.vertices)
        verts[:,which] = (axis - verts[:,which])
        return verts
    
    def nothelp__(self, key, help = None):
        """experimental - I forgot what this was about!
        """
        return self.__getattribute__(key, help=True)
        
    def help(self):
        """Recursive help documentation
        print 'to get help for the Bspline class,'
        print 'type '
        print '>>>help(Bspline)'
        """
        print 'to get help for the Bspline class,'
        print 'I bet you already knew this:'
        print 'type '
        print '>>>help(Bspline)'
        return

    def plot(self, color_='black', vertices = False):
        self.plotcurve(color_,vertices)
        return
    
    #"""   
    def plotcurve(self, color_='black',vertices=False):
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #fig.add_subplot(ax)
        if self.dim ==2:
            plt.plot(self.r[:,0],self.r[:,1], color=color_)
            if vertices:
                plt.plot(self.vertices[:,0],
                         self.vertices[:,1],
                         alpha=.4,
                         color = 'blue',
                         marker = 's', )
            plt.axis('equal')
            #plt.show()
        #self.ax = ax
        if self.dim ==3:
            self.plot3D()
        return
    #"""
    """
    def plotcurve(curve, color_='black'):
        if curve.dim ==2:
            plt.plot(curve.r[:,0],curve.r[:,1], color=color_)
            plt.axis('equal')
            plt.show()
        if curve.dim ==3:
            curve.plot3D()
        return
    #"""
    def on_pick(self, event):
        self.event = event
        artist = event.artist
        print 'on pick'
        print event.x, event.y
        self.store_event_pick = event
        return
    
    def get_ind_under_point(self, event):
        """get the index of the vertex under point if within epsilon tolerance
        """  
        print 'store event in gettin ind'
        self.store_event = event
        #xy = np.asarray(self.controlpoly.xy)
        #xyt = self.controlpoly.get_transform().transform(xy)
        #xt, yt = xyt[:, 0], xyt[:, 1]
        #print xt, yt
        #xy = np.asarray(self.controlpoly.xy)
        xyt = ax.get_transform().transform(self.vertices)
        xt, yt = xyt[:, 0], xyt[:, 1]
        print xt
        print yt
        #inv = ax.transData.inverted()   
        #xt,yt = self.vertices[:, 0], self.vertices[:, 1]
        #xydata = inv.transform((event.x,event.y))
        #print xydata
        #ex = xydata[0]
        #ey = xydata[1]
        print event.x, event.y
        d = np.sqrt((xt-event.x)**2 + (yt-event.y)**2)
        indseq = np.nonzero(np.equal(d, np.amin(d)))[0]
        ind = indseq[0]
        print ind
        self._ind = ind
        #if d[ind]>=5:
        #    ind = None
        return ind
        
    def button_press_callback(self, event):
        if self.verbose:print 'button_press_callback'
        print ' a mouse button is pressed'
        #if not self.showverts: return
        #if event.inaxes==None: return
        #if event.button != 1: return
        self._ind = self.get_ind_under_point(event)
        print 'got ind = ',self._ind
        #return #do not return
        
    def button_release_callback(self, event):
        print 'mouse button is released'
        #if not self.showverts: return
        if event.button != 1: return
        self._ind = None
        
    def motion_notify_callback(self, event):
        if self.verbose:print 'motion_notify_callback'
        #if not self.showverts: return
        if self._ind is None: return
        if event.inaxes is None: return
        if event.button != 1: return
        x,y = event.xdata, event.ydata
        print x,y
        self.vertices[self._ind,0] = x
        self.vertices[self._ind,1] = y
        self.compute_curve_plot()
        return 
        
    def compute_curve_plot(self):
        fig.canvas.draw()
        self.allCurveArray()
        
        line = plt.plot(self.r[:,0],
                self.r[:,1],
                color = 'black',
                picker=True)
        verts = plt.plot(self.vertices[:,0],
                 self.vertices[:,1], 
                 color = 'blue',
                 marker = 's', 
                 alpha = .4)
        fig.canvas.blit(ax.bbox)

        return
        
    def plotcurve_detailed(self, 
                           curvature = 'no', 
                           scale=1.0, 
                           color_='black',
                           picker=5, 
                           canvas = None,
                           normalize=False,
                           scale_curve=None,
                           fill=False,
                           fillcolor='green'):
        
        
        if canvas is None:
            fig, ax = plt.subplots()
        else:
            ax = canvas
            

        if self.dim ==2:
            
            if scale_curve is None:
                scale_curve = self
            if normalize:
                xbds = scale_curve.extremes1d(0)
                ybds = scale_curve.extremes1d(1)
                line = ax.plot(self.r[:,0]/xbds[1],
                        self.r[:,1]/ybds[1],
                        color = color_,
                        picker=False)
                        
                verts = ax.plot(self.vertices[:,0]/xbds[1],
                         self.vertices[:,1]/ybds[1], 
                         color = 'blue',
                         marker = 's', 
                         alpha = .4)
                if curvature == 'yes':
                    scale_ = scale
                    self.plotCurvature(scale = scale_)
                plt.axis('equal')
                
            else:
                xbds = 1.
                ybds = 1.
                line = ax.plot(self.r[:,0],
                        self.r[:,1],
                        color = color_,
                        picker=False)
                        
                verts = ax.plot(self.vertices[:,0],
                         self.vertices[:,1], 
                         color = 'blue',
                         marker = 's', 
                         alpha = .4)
                
                if curvature == 'yes':
                    scale_ = scale
                    self.plotCurvature(scale = scale_)
                plt.axis('equal')
            if fill:
                xmin, xmax,ymin, ymax = self.extreme_C0()
                ax.fill_between(self.r[:,0],
                                 ymin,
                                 self.r[:,1],
                                 facecolor = fillcolor, 
                                 alpha=.1,
                                 label='area')
                
            
        
        elif self.dim ==3:
            pass
            """
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.autoscale_view(tight=None, scalex=1., scaley=1., scalez=True)
            ax.plot(self.vertices[:,0],
                    self.vertices[:,1],
                    self.vertices[:,2],  
                    linestyle = "--")
            ax.plot(self.r[:,0], self.r[:,1], self.r[:,2])
            #if curvature == 'yes':
            #    scale_ = scale
            #    self.plotCurvature(scale = scale_)
            plt.axis('equal')
            plt.show()
            #"""
        return ax
        
    def plot_area_to_x(self, color='green'):
        xmin, xmax,ymin, ymax = self.extreme_C0()
        plt.fill_between(self.r[:,0],ymin,self.r[:,1],
                            facecolor = color, alpha=.1,
                            label='area')
        return
    
    def plot_area_to_x_shiftaxis(self, 
                                 color='green', 
                                 shift_axis = 0.):
        xmin, xmax,ymin, ymax = self.extreme_C0()
        plt.fill_between(self.r[:,0],shift_axis,self.r[:,1],
                            facecolor = color, alpha=.1,
                            label='area')
        return
        
    def plot_annotation(self, text = 'Hello World!', loc = (.5,.5)):
        ax = self.ax        
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        bbox_props = dict(boxstyle="rarrow,pad=0.3", fc="cyan", ec="b", lw=2)
        tbox = ax.text(loc[0],loc[1],text,ha="center", va="center",size=15)
        tbox.xy = loc
        tr = FastDraggableText(tbox)
        tr.connect()
        plt.show()
        return
    
    def plot3D(self, canvas = None, color='black'):
        """
            standardized a bit and commented out the experimental
            picker...
        """
        if canvas is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            ax = canvas
        dummy = self.vertices[-1]-self.vertices[0]
        l=dummy[0]
        w=dummy[1]
        h=dummy[2]
        ##ax = fig.add_subplot(111, projection='3d')
        #ax = fig.gca(projection='3d')
        #ax.autoscale_view(tight=None, scalex=1., scaley=1., scalez=True)
        try:
            ax.plot(self.vertices[:,0],
                    self.vertices[:,1],
                    self.vertices[:,2],   
                    marker="o",
                    linestyle = "-")
                    #picker=False)
            ax.plot(self.r[:,0], 
                    self.r[:,1], 
                    self.r[:,2])
        except:
            ax.plot(self.vertices[:,0],
                    self.vertices[:,1],  
                    marker="o",
                    linestyle = "-")
            ax.plot(self.r[:,0], self.r[:,1])        
        #ax.set_xlim3d(-l, l)
        #ax.set_ylim3d(-w,w)
        #ax.set_zlim3d(-h,h) 
        """
        fig.canvas.mpl_connect('pick_event', 
                                   self.on_pick)
        fig.canvas.mpl_connect('button_press_event', 
                               self.button_press_callback)
        fig.canvas.mpl_connect('button_release_event', 
                               self.button_release_callback)
        fig.canvas.mpl_connect('motion_notify_event', 
                               self.motion_notify_callback)
        #"""
        return ax
    
    def plot3D_hull_system(self, canvas = None, 
                           color='.01',
                           znotation=True,
                           plot_vertices=False):
        """
            standardized a bit and commented out the experimental
            picker...
        """
        if canvas is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            ax = canvas
        if plot_vertices:
            ax.plot(self.vertices[:,2],
                    self.vertices[:,0],
                    self.vertices[:,1],   
                    marker="o",
                    linestyle = "",
                    color = color)
        ax.plot(self.r[:,2], 
                self.r[:,0], 
                self.r[:,1],
                color = color)
        return ax
        
    def plot_subplot_curve(curve, canvas = None, 
                           color = None, Lspline = None):
        """Implement subplotting for the ADILS module
        """
        if canvas is None: return 'Error, subplot requires canvas'
        elif color is None:color = 'black'
        ax  = canvas
        ax.plot(curve.r[:,0],curve.r[:,1], color, marker = '')
        ax.plot(curve.vertices[:,0],curve.vertices[:,1], 
                color = 'blue', marker = 's', alpha = .4)
        return ax


    def plotAllDers(curve):
        k=curve.k
        plt.plot(curve.allr[:,0,0], curve.allr[:,0,1])
        plt.plot(curve.r[:,0], curve.r[:,1],linestyle = "--")
        #plt.plot(curve.r[:,0],curve.dr[:,1]/curve.dr[:,0])
        for j in range(1,k-1):
            plt.plot(curve.allr[:,0,0], curve.allr[:,j,1]/curve.allr[:,j,0], 
                     linestyle = "--")
            #plt.plot(curve.allr[:,j,0], curve.allr[:,j,1])
        plt.axis('equal')
        plt.show()
        return


    def plotbasis(curve,which=0):
        print 'plotting the {} level basis Bik'.format(which)
        plt.plot(curve.s,curve.allbasis[:,which])
        #plt.axis('equal')
        plt.ylabel('Value')
        plt.xlabel('Parameter Space')
        plt.show()

    def plot3d_wIvertices(curve):
        dummy = curve.vertices[-1]-curve.vertices[0]
        #curveSize = np.linalg.norm(dummy)
        l=dummy[0]
        w=dummy[1]
        h=dummy[2]
        #fig = plt.figure(figsize=(9.5,5.0))
        fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #ax = fig.gca(projection='3d')
        ax = fig.add_subplot(111, projection='3d')
        ax.autoscale_view(tight=None, scalex=1., scaley=1., scalez=True)
        try:
            ax.plot(curve.vertices[:,0],curve.vertices[:,1],curve.vertices[:,2],  linestyle = "--")
            ax.plot(curve.r[:,0], curve.r[:,1], curve.r[:,2])
        except:
            ax.plot(curve.vertices[:,0],curve.vertices[:,1],  marker="x",linestyle = "--")
            ax.plot(curve.r[:,0], curve.r[:,1])
        #ax.set_xlim3d(-50, 50)
        #ax.set_ylim3d(-50,50)
        #ax.set_zlim3d(-50,50)
        #ax.set_xlim3d(-curveSize, curveSize)
        #ax.set_ylim3d(-curveSize,curveSize)
        #ax.set_zlim3d(-curveSize,curveSize)
        
        #        ax.set_xlim3d(-l, l)
        #        ax.set_ylim3d(-w,w)
        #        ax.set_zlim3d(-h,h) 
        
        axmax = max(curve.extremes1d(0))
        aymax = max(curve.extremes1d(1))
        azmax = max(curve.extremes1d(2))
        biggest = max(axmax,aymax,azmax)
        #
        ax.set_xlim3d(-biggest*.0,1.*biggest)
        ax.set_ylim3d(-biggest*.5,.5*biggest)
        ax.set_zlim3d(-biggest*.5,.5*biggest)
        
        plt.show()
        return


    def plotcurveList(self, curveList):
        for curve in curveList:
            plt.plot(curve.r[:,0],curve.r[:,1])
            plt.axis('equal')
        plt.show()
        return



    def plot3DmultiList(curve, 
                        transverseList, 
                        longitudinalList, 
                        unlimited = True,
                        minx=0.,miny=0.,minz=0.,
                        limx=None,limy=None,limz=None,
                        view='z-vertical'):
        """
        """
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        
        if view == 'z-vertical':
            for item in transverseList:
                try:
                    ax.plot(item.ivertices[:,0],item.ivertices[:,1],item.ivertices[:,2],marker = "o",  linestyle = "none", color='k')
                except:
                    pass
                ax.plot(item.vertices[:,0],
                        item.vertices[:,1],
                        item.vertices[:,2],
                        marker = "4",  linestyle = "--", color='g')
                ax.plot(item.r[:,0], 
                        item.r[:,1], 
                        item.r[:,2], label='parametric curve', color='b')
            for item in longitudinalList:
                try:
                    ax.plot(item.ivertices[:,0],
                            item.ivertices[:,1],
                            item.ivertices[:,2],marker = "o",  linestyle = "none", color='g')
                except:
                    pass
                ax.plot(item.vertices[:,0],
                        item.vertices[:,1],
                        item.vertices[:,2],marker = "o",  linestyle = "--", color='g')
                ax.plot(item.r[:,0], item.r[:,1], item.r[:,2], label='parametric curve', color='r')
        #
        #**********************************************************************
        #
        elif view == 'y-vertical':
            for item in transverseList:
                try:
                    ax.plot(item.ivertices[:,0],item.ivertices[:,1],item.ivertices[:,2],marker = "o",  linestyle = "none", color='k')
                except:
                    pass
                ax.plot(item.vertices[:,0],
                        item.vertices[:,2],
                        item.vertices[:,1],
                        marker = "4",  linestyle = "--", color='g')
                ax.plot(item.r[:,0], 
                        item.r[:,2], 
                        item.r[:,1], label='parametric curve', color='b')
            for item in longitudinalList:
                try:
                    ax.plot(item.ivertices[:,0],
                            item.ivertices[:,2],
                            item.ivertices[:,1],marker = "o",  linestyle = "none", color='g')
                except:
                    pass
                ax.plot(item.vertices[:,0],
                        item.vertices[:,2],
                        item.vertices[:,1],marker = "o",  linestyle = "--", color='g')
                ax.plot(item.r[:,0], 
                        item.r[:,2], 
                        item.r[:,1], label='parametric curve', color='r')
                #ax.plot(item.r[:,0], item.r[:,1], item.r[:,2], label='parametric curve', color='r')
  
        #
        #**********************************************************************
        #
        elif view == 'x-vertical':
            for item in transverseList:
                try:
                    ax.plot(item.ivertices[:,0],item.ivertices[:,1],item.ivertices[:,2],marker = "o",  linestyle = "none", color='k')
                except:
                    pass
                ax.plot(item.vertices[:,1],
                        item.vertices[:,2],
                        item.vertices[:,0],
                      marker = "4",  linestyle = "--", color='g')
                ax.plot(item.r[:,1], 
                        item.r[:,2], 
                        item.r[:,0], label='parametric curve', color='b')
            for item in longitudinalList:
                try:
                    ax.plot(item.ivertices[:,1],
                            item.ivertices[:,2],
                            item.ivertices[:,0],marker = "o",  linestyle = "none", color='g')
                except:
                    pass
                ax.plot(item.vertices[:,1],
                        item.vertices[:,2],
                        item.vertices[:,0],marker = "o",  linestyle = "--", color='g')
                ax.plot(item.r[:,1], 
                        item.r[:,2], 
                        item.r[:,0], label='parametric curve', color='r')
                
                #ax.plot(item.r[:,0], item.r[:,1], item.r[:,2], label='parametric curve', color='r')
  

        #
        #**********************************************************************
        #
        #if unlimited:
        #todo abs/norm
        if limx is not None:
            ax.set_xlim3d(minx,limx)
            ax.set_ylim3d(miny,limy)
            ax.set_zlim3d(minz,limz)
        else:
            axmax = max(curve.extremes1d(0))
            aymax = max(curve.extremes1d(1))
            azmax = max(curve.extremes1d(2))
            biggest = max(axmax,aymax,azmax)
            #
            ax.set_xlim3d(-biggest*.0,1.*biggest)
            ax.set_ylim3d(-biggest*.5,.5*biggest)
            ax.set_zlim3d(-biggest*.5,.5*biggest)
        plt.show()

    def compute_curve_of_curvature(self, crange=None, scale=1.0):
        if crange is None:
            crange = self.nump
            s = self.s
        else:
            s = np.linspace(start=0.,stop=1.,num=crange)
        ortho_vector = np.asarray([0.,0.,1.])
        cvp = np.zeros((crange,2),float)
        #x=[]
        #y=[]
        for i in range(self.nump):
            cuv = self.compute_curvature(s[i])
            curvature = cuv
            slope     = self.compute_tangent(s[i])
            tanpoint = [np.cos(slope),np.sin(slope),0.]
            vector = crossproduct(ortho_vector, tanpoint)
            tanpoint = [self.r[i][0]+vector[0],self.r[i][1]+vector[1]]

            curvature_point = [self.r[i][0]+scale*curvature*vector[0],
                               self.r[i][1]+scale*curvature*vector[1]]
            #x.append(curvature_point[0])
            #y.append(curvature_point[1])   
            cvp[i] = curvature_point

        #x = np.asarray(x)
        #y = np.asarray(y)
        return cvp
        
    def plotCurvature(self, scale=1.0, color_ = 'grey', alpha=0.):
        ortho_vector = np.asarray([0.,0.,1.])
        plt.plot(self.r[:,0],self.r[:,1])
        alpha_ = alpha
        x = []
        y = []
        for i in range(self.nump):
            cuv = self.compute_curvature(self.s[i])
            curvature = cuv
            slope     = self.compute_tangent(self.s[i])
            tanpoint = [np.cos(slope),np.sin(slope),0.]
            vector = crossproduct(ortho_vector, tanpoint)
            tanpoint = [self.r[i][0]+vector[0],self.r[i][1]+vector[1]]

            curvature_point = [self.r[i][0]+scale*curvature*vector[0],
                               self.r[i][1]+scale*curvature*vector[1]]
            x.append(curvature_point[0])
            y.append(curvature_point[1])
            #plt.axis('equal')
            plt.grid()
            plt.plot([self.r[i][0],curvature_point[0]],
                     [self.r[i][1],curvature_point[1]], 
                    color = color_)
        
        x = np.asarray(x)
        y = np.asarray(y)
        plt.plot(x, y, color = color_, alpha = .4)
        plt.fill_between(self.r[:,0], 
                            y, 
                            self.r[:,1], 
                            facecolor = color_, 
                            alpha=alpha_,
                            label='curvature')
        #plt.xticks([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.])
        #plt.yticks([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.])
        #plt.show()

        return
        
    def plotCurvature_expensive(self, scale=1.0, color_ = 'grey', 
                                alpha=0.,nump=1,flip=False):
        """for planar curves
        
            debug:
                i=15
                alpha=0.
                scale = 3.5
                nump=5
        """
        num = self.nump*nump
        s = np.linspace(self.t[0],self.t[-1],num)
        r = np.zeros((num,self.dim),float) 
        for i,el in enumerate(s):
            r[i] = self.CurvePoint(el)
        ortho_vector = np.asarray([0.,0.,1.])
        plt.plot(r[:,0],r[:,1])
        alpha_ = alpha
        x = []
        y = []
        if not flip:
            for i in range(num):
                curvature = self.compute_curvature(s[i])
                slope     = self.compute_tangent(s[i])
                tanpoint = [np.cos(slope),np.sin(slope),0.]
                vector = crossproduct(ortho_vector, tanpoint)
                #tanpoint = [r[i][0]+vector[0],r[i][1]+vector[1]]
    
                curvature_point = [r[i][0]+scale*curvature*vector[0],
                                   r[i][1]+scale*curvature*vector[1]]
                x.append(curvature_point[0])
                y.append(curvature_point[1])
                
                #plt.grid()
                plt.plot([r[i][0],curvature_point[0]],
                         [r[i][1],curvature_point[1]], 
                        color = color_)
        elif  flip:
            for i in range(num):
                curvature = self.compute_curvature(s[i])
                slope     = self.compute_tangent(s[i])
                tanpoint = [np.cos(slope),np.sin(slope),0.]
                vector = crossproduct(ortho_vector, tanpoint)
                #tanpoint = [r[i][0]+vector[0],r[i][1]+vector[1]]
    
                curvature_point = [r[i][0]-scale*curvature*vector[0],
                                   r[i][1]-scale*curvature*vector[1]]
                x.append(curvature_point[0])
                y.append(curvature_point[1])
                
                #plt.grid()
                plt.plot([r[i][0],curvature_point[0]],
                         [r[i][1],curvature_point[1]], 
                        color = color_)
        
        x = np.asarray(x)
        y = np.asarray(y)
        plt.plot(x, y, color = color_, alpha = .4)
        plt.axis('equal')
        if alpha == 0.:
            pass
        else:
            plt.fill_between(r[:,0], 
                                y, 
                                r[:,1], 
                                facecolor = color_, 
                                alpha=alpha_,
                                label='curvature')
        #plt.xticks([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.])
        #plt.yticks([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.])
        #plt.show()

        return
    
    def plotCurvature_spines(self, scale=1.0, color_ = 'grey',factor=1):
        ortho_vector = np.asarray([0.,0.,1.])
        #plt.plot(self.r[:,0],self.r[:,1])
        x = []
        y = []
        for i in range(self.nump):
            cuv = self.compute_curvature(self.s[i])
            curvature = cuv
            slope     = self.compute_tangent(self.s[i])
            tanpoint = [np.cos(slope),np.sin(slope),0.]
            vector = crossproduct(ortho_vector, tanpoint)
            tanpoint = [self.r[i][0]+vector[0],self.r[i][1]+vector[1]]

            curvature_point = [self.r[i][0]+scale*curvature*vector[0],
                               self.r[i][1]+scale*curvature*vector[1]]
            x.append(curvature_point[0])
            y.append(curvature_point[1])
            
            plt.plot([self.r[i][0],curvature_point[0]],
                     [self.r[i][1],curvature_point[1]], 
                    color = color_)
        plt.grid()
        return
    
    def plotCurvature_nospines(self, 
                               scale=1.0, 
                               color_ = 'grey',
                               factor=2,
                               alpha=.0):
        hs = Bspline(self.vertices, self.k, self.nump*factor)
        ortho_vector = np.asarray([0.,0.,1.])
        #plt.plot(hs.r[:,0],hs.r[:,1])
        x = []
        y = []
        for i in range(hs.nump):
            cuv = hs.compute_curvature(hs.s[i])
            curvature = cuv
            slope     = hs.compute_tangent(hs.s[i])
            tanpoint = [np.cos(slope),np.sin(slope),0.]
            vector = crossproduct(ortho_vector, tanpoint)
            tanpoint = vector[:2] + hs.r[i]

            curvature_point = [hs.r[i][0]+scale*curvature*vector[0],
                               hs.r[i][1]+scale*curvature*vector[1]]
            x.append(curvature_point[0])
            y.append(curvature_point[1])
            #plt.grid()
            #plt.plot([hs.r[i][0],curvature_point[0]],
            #         [hs.r[i][1],curvature_point[1]], 
            #        color = color_)
        alpha_=alpha
        x = np.asarray(x)
        y = np.asarray(y)
        plt.plot(x, y, color = color_, alpha = .4)
        """# issue - fill between works between vertical points!
        # the area shaded below is not correct!
        plt.fill_between(hs.r[:,0], 
                            y, 
                            hs.r[:,1], 
                            facecolor = color_, 
                            alpha=alpha_,
                            label='curvature')
        #"""
        # the fill between below get the right outline at outer edge
        # at the expense of hacing the wrong inner border
        plt.fill_between(x, 
                            hs.r[:,0], 
                            y,
                            facecolor = color_, 
                            alpha=alpha_,
                            label='curvature')
        return

            
#/*!------------------------------------------------------------------------------------------------------------------------!*/#
        
class interpolatedBspline(Bspline):

    
    def __init__(self, vertices, k, nump=50, t=None):
        """
            Assume knot vector ranges from 0 to 1. 
        """
        self.verbose        = False
        self.vertices_copy  = np.asarray(vertices)  
        self.vertices       = np.asarray(vertices)  # vertices
        self.dim            = len(self.vertices[0])
        self.n              = len(self.vertices)    # number of vertices
        self.trange         = np.array([0.,1.])     # parametric space
        self.p              = k-1                   # degree
        #self.matchParam     = np.copy(self.vertices)
        self.nump           = nump
        
        self._t_observers   = []
        if t is None:
            self.t_copy         = make_knot_vector(k,len(vertices)-1,self.trange,type='open')
        else:
            self.t_copy         = t 
        self.t                  = copy.deepcopy(self.t_copy)
        self._k_observers   = []
        self._k_simple_observers = []
        self.k_copy         = k
        self.k              = k                     # order
        
        
        ##
        ##**************************************************
        ## Curve Coordinates
        Origin = np.asarray([0.,0.,0.])
        x, y, z = np.identity(3)
        self.ccs = Frame( np.asarray([x,y,z]) ) #curve coordinate system
        self.origin = DeepVector(v = Origin,
                                B=self.ccs)
        self.rots = Quaternion.from_v_theta(x, 0.) #null rotation vector
        ##
        ##**************************************************
        ##
        
        self.s              = np.linspace(self.t[0],self.t[-1],self.nump,endpoint=True)  #param space
        self.r              = np.zeros((self.nump,self.dim),float)                              # store B-spline
        self.dr             = np.zeros((self.nump,self.dim),float)                              # store DB-spline
        self.span           = np.zeros((self.nump),int)                                # store non zero span of basis funcs
        self.basis          = np.zeros((self.nump,self.n),float)                        # store all basis func 
        self.dbasis         = np.zeros((self.nump,self.n),float) #1st derivative - old way check
        self.allbasis       = np.zeros((self.nump,self.k,self.n),float)                        # store all dbasis func 
        self.allr           = np.zeros((self.nump,self.k,self.dim),float)

        self.area           = None
        self.Mx             = None
        self.My             = None
        self.curvature      = None

        # Fairness Functional: Integrals of the Matrix of products of Basis Function Derivatives
        self.testM1         = np.zeros( (self.n,self.n), float) # write these the hard way - to test M1 computations once and for all
        self.M1             = np.zeros( (self.n,self.n), float) # these are actually used for the integrals - not M itself
        self.M2             = np.zeros( (self.n,self.n), float) # these are actually used for the integrals - not M itself
        self.M3             = np.zeros( (self.n,self.n), float) # these are actually used for the integrals - not M itself
        self.M              = np.zeros( (self.n,self.n,self.nump), float) # Store Non-Differentiated Products of Basis Functions Derivatives.  -ALL!
        self.E1             = None    # fairness functions
        self.E2             = None    # fairness functions
        self.E3             = None    # fairness functions
        self.AL             = None    # arc length

        # Further Curve Properties:
        self.AM1            = np.zeros( (self.n,self.n), float)         # area computations: storing Ibasis products
        self.AM2            = np.zeros( (self.n,self.n,self.n), float)  # moment computations: storing Ibasis products

        # Use Automatic Differentiation
        self.xpts = None
        self.ypts = None
        self.zpts = None
        #self.lpts = None
        
        self.data = {}

        # Intrinsic Method Calls : Fully Initiazlize the B-Spline
        self.ivertices      = np.copy(self.vertices) #copy old vertices for plotting
        self.GlobalCurveInterp() #find new vertices and knot vector to interpolate the old vertices
        self.SpanArray()
        self.BasisArray()
        self.DBasisArray() #basis functions and all derivatives
        self.allCurveArray() # compute the curve and all existing derivatives at every parameter value


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
        self.ukbar = uk_bar
        return

    def GlobalCurveInterp_gen1_THB_compatible_ver0(self,vertices):
        """P&T algorithm page 369
            Global Curve Interpolation
            through n+1 points
            
            
            version 0: equally spaced ukbar 
            eually spaced real knots

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
                
        
        notes
        ----------
            See P&T, page 364-365, especially 9.7, choosing knots
            keeping in mind 9.3-9.5
            -this may give singular matrices when used in the wild.
                
        """
        
        n = self.n-1
        Q = vertices
        r = self.dim
        p = self.p
        m = n+p+1
        U = self.t #comment out to JH mod form

        
        #Find the Chord Length: (9.5)
        d=0.
        store_sections = []
        for k in range(1,n+1):
            drule=0. 
            #loop over dimensions
            for kk in range(r):
                drule += (Q[k,kk]-Q[k-1,kk])**2  #version 1: centripetal (not used)
            stored = np.sqrt(drule)
            store_sections.append(stored)
        d = sum(store_sections)
            
        #Define the uk_bar: (9.6)
        uk_bar=np.zeros((n+1),float)
        uk_bar[n]=1.0
        for k in range(1,n):
            #uk_bar[k]= uk_bar[k-1]+ ( store_sections[k-1]/d )
            #************************ 
            # version 0:
            # equally spaced ukbar (9.3)
            uk_bar[k]= k/float(n)
            #*************************
        #Compute the new Knot Vector U: (9.8)
        #U = np.zeros((m+1),float) #JH mod
        for j in range(1,n-p+1):
            U[j+p]=0.
            for i in range(j,j+p):
                #U[j+p]+=uk_bar[i]
                #************************ 
                #THB compatible knots (warning, possibly singular)
                U[j+p]+=j/(n-p+1.) 
                #especially if used with chord or centripetal parameterization
                #************************ 
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
        self.ukbar = uk_bar
        return

    def GlobalCurveInterp_gen1_THB_compatible_ver1(self,vertices):
        """P&T algorithm page 369
            Global Curve Interpolation
            through n+1 points
            
            
            version 1: Centripetal ukbar 
            eually spaced real knots

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
                
        
        notes
        ----------
            See P&T, page 364-365, especially 9.7, choosing knots
            keeping in mind 9.3-9.5
            -this may give singular matrices when used in the wild.
                
        """
        
        n = self.n-1
        Q = vertices
        r = self.dim
        p = self.p
        m = n+p+1
        U = self.t #comment out to JH mod form

        
        #Find the Chord Length: (9.5)
        d=0.
        store_sections = []
        for k in range(1,n+1):
            drule=0. 
            #loop over dimensions
            for kk in range(r):
                drule += (Q[k,kk]-Q[k-1,kk])**2  #version 1: centripetal
            stored = np.sqrt(drule)
            store_sections.append(stored)
        d = sum(store_sections)
            
        #Define the uk_bar: (9.6)
        uk_bar=np.zeros((n+1),float)
        uk_bar[n]=1.0
        for k in range(1,n):
            #************************ 
            # version 1:  
            #  using centripetal ukbar spacing
            uk_bar[k]= uk_bar[k-1]+ ( store_sections[k-1]/d )
            #************************ 
            #equally spaced ukbar (9.3)
            #uk_bar[k]= k/float(n)
            
        #Compute the new Knot Vector U: (9.8)
        #U = np.zeros((m+1),float) #JH mod
        for j in range(1,n-p+1):
            U[j+p]=0.
            for i in range(j,j+p):
                #U[j+p]+=uk_bar[i]
                #************************ 
                #THB compatible knots (warning, possibly singular)
                U[j+p]+=j/(n-p+1.) 
                #especially if used with chord or centripetal parameterization
                #************************ 
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
        self.ukbar = uk_bar
        return

    def GlobalCurveInterp_gen1_THB_compatible_ver2(self,vertices):
        """P&T algorithm page 369
            Global Curve Interpolation
            through n+1 points
            
            version 2: chord length ukbar parameterization
            eually spaced real knots

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
                
        
        notes
        ----------
            See P&T, page 364-365, especially 9.7, choosing knots
            keeping in mind 9.3-9.5
            -this may give singular matrices when used in the wild.
                
        """
        
        n = self.n-1
        Q = vertices
        r = self.dim
        p = self.p
        m = n+p+1
        U = self.t #comment out to JH mod form

        
        #Find the Chord Length: (9.5)
        d=0.
        store_sections = []
        for k in range(1,n+1):
            drule=0. 
            #loop over dimensions
            for kk in range(r):
                drule += abs((Q[k,kk]-Q[k-1,kk])) #version 2: chord length
            store_sections.append(drule)
        d = sum(store_sections)
            
        #Define the uk_bar: (9.6)
        uk_bar=np.zeros((n+1),float)
        uk_bar[n]=1.0
        for k in range(1,n):
            #************************ 
            # version 2:  
            #  using chord length ukbar spacing
            uk_bar[k]= uk_bar[k-1]+ ( store_sections[k-1]/d )
            #************************ 
            #equally spaced ukbar (9.3)
            #uk_bar[k]= k/float(n)
            
        #Compute the new Knot Vector U: (9.8)
        #U = np.zeros((m+1),float) #JH mod
        for j in range(1,n-p+1):
            U[j+p]=0.
            for i in range(j,j+p):
                #U[j+p]+=uk_bar[i]
                #************************ 
                #THB compatible knots (warning, possibly singular)
                U[j+p]+=j/(n-p+1.) 
                #especially if used with chord or centripetal parameterization
                #************************ 
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
        self.ukbar = uk_bar
        return


    #def plot3d_wIvertices(curve, l=50,w=50,h=50):
    def plot3d_wIvertices(curve):
        dummy = curve.vertices[-1]-curve.vertices[0]
        curveSize = np.linalg.norm(dummy)
        fig = plt.figure()
        #fig = plt.figure(figsize=(9.5,5.0))
        #ax = fig.add_subplot(111, projection='3d')
        #ax = fig.gca(projection='3d')
        ax = fig.add_subplot(111, projection='3d')
        #ax.autoscale_view(tight=None, scalex=1., scaley=1., scalez=True)
        try:
            ax.plot(curve.ivertices[:,0],curve.ivertices[:,1],curve.ivertices[:,2],marker = "o",  linestyle = "none")
            ax.plot(curve.vertices[:,0],curve.vertices[:,1],curve.vertices[:,2],  linestyle = "--")
            ax.plot(curve.r[:,0], curve.r[:,1], curve.r[:,2])
        except:
            ax.plot(curve.ivertices[:,0],curve.ivertices[:,1],marker = "o",  linestyle = "none")
            ax.plot(curve.vertices[:,0],curve.vertices[:,1],  marker="x",linestyle = "--")
            ax.plot(curve.r[:,0], curve.r[:,1])
        ax.set_xlim3d(-curveSize, curveSize)
        ax.set_ylim3d(-curveSize,curveSize)
        ax.set_zlim3d(-curveSize,curveSize)
        plt.show()

class curveSet(Bspline):
    """
        Set of Boundary Curves
        
        assumes standard curve type (interpolates the first and last point)
    """
    def __init__(self, curveiteratable, check_order = False):
        self.reverse_area = False
        if type(curveiteratable) == type([]):
            self.curveset = self.orderCurves(curveiteratable, flag=check_order)
            self.computeMoments()
        elif type(curveiteratable) == type({}):
            dum = []
            for key in curveiteratable:
                print 'curveiteratable = ',curveiteratable[key]
                dum.append(curveiteratable[key])
            self.curveset = dum
        else:
            print 'Error, curve set does not match available types'
        return
        
    def reverse_curve(self, curve):
        new_vertices = []
        for vertex in reversed(curve.vertices):
            new_vertices.append(vertex)
        new_vertices = np.asarray(new_vertices)
        new_curve = Bspline(new_vertices, curve.k, curve.nump)
        return new_curve   
        
    def orderCurves(self, curvelist, flag=False):
        #routines.crossproduct
        tol = .1
        
        origin  = curvelist[0].vertices[0]
        
        # check order of 1st two curves:
        print 'checking 1st curve orientation'
        orient = np.linalg.norm(curvelist[0].vertices[-1] - curvelist[1].vertices[0])
        #print orient
        if abs(orient) > tol:  
            print 'reverse second curve'
            curvelist[1] = self.reverse_curve(curvelist[1])
            orient = np.linalg.norm(curvelist[0].vertices[-1] - curvelist[1].vertices[0])
            if abs(orient) <tol:
                print 'aligned 1st and 2nd curve nose to tale with success!'
            else:
                print 'error, 1st and 2nd curves do not meet at the end of the first curve.'
        else:
            print '1st and 2nd curves area aligned'
        
        v1 = curvelist[0].vertices[-1] - origin
        v2 = curvelist[1].vertices[-1] - curvelist[0].vertices[0]
        self.poly_handedness = routines.crossproduct(v2,v1)[2]
        
        if self.poly_handedness < 0.:
            self.reverse_area = True
            
        print 'polyhandedness = ', self.poly_handedness         
        if flag == True:
            print 'checking poly orientation'

            for i in range(1,len(curvelist)):
                start   = curvelist[i].vertices[0] - curvelist[0].vertices[0]
                end     = curvelist[i].vertices[-1] - curvelist[0].vertices[0]
                v2 = start - origin
                v3 = end - origin
                sign = routines.crossproduct(v3,v2)
                print 'sign = ',sign
                if (sign[2] < 0. and self.poly_handedness > 0.):
                    print 'flip curve ', i
                    curvelist[i] = self.reverse_curve(curvelist[i])
                elif (sign[2] > 0. and self.poly_handedness < 0.):
                    print 'flip curve ', i
                    curvelist[i] = self.reverse_curve(curvelist[i])

            curveset = curvelist
        else:
            curveset = curvelist
            
            
        return curveset
    
    def computeFairness(self):
        for curve in self.curveset:
            curve.basis_matrix()
        for curve in self.curveset:
            curve.fairness
        return
      
    def computeMoments(self):
        for curve in self.curveset:
            curve.MomentMatrices()
        return
        
        
    def compute_area(self):
        new_v = copy.deepcopy(self.vertices)
        self.area = 0.
        self.computeMoments()
        #print 'printing curveset:',self.curveset
        for curve in self.curveset:
            print curve.vertices
            xpts = curve.vertices[:,0]
            ypts = curve.vertices[:,1]
    
            #Area calculations:
            curve_area       = -0.5*np.dot(np.dot(xpts,np.transpose(curve.AM1)),ypts)
            self.area += curve_area
        #if self.reverse_area == True:
        #    self.area *= -1.
        return
        
    def plot_set(self):
        for curve in self.curveset:
            #s = np.ma.masked_array(curve.s, np.idff(curve.s[0]>=0)
            x = curve.r[0]
            y = curve.r[1]
            plt.plot(x,y)
            
        return
        

#-----------------------------------------------------------------------------------------
class BsplineSurface(Bspline):

    def __init__(self, tcurveNet, lcurveNet, surfNet=None, 
                 kind='Loft',surface=None, surfacept=None, 
                 dim=3,
                 xax = (1.,0.,0.), 
                 yax = (0.,1.,0.), 
                 zax = (0.,0.,1.),
                 current_rot = None):
        """
        init assigns the attributes
        of the B-Spline Surface

        - Quick and Dirty way is relient on the curves which make up the surface

        Loops:
             innermost  : dimension
                middle  : longitudinal
                 outer  : transverse
        ie:
            [u,v][dim]
            
        """
        self.dim    = dim # 3D vertices are expected
        
        self.utrange = tcurveNet[0].trange
        self.uorder  = tcurveNet[0].k
        self.uknots  = tcurveNet[0].t
        self.un      = tcurveNet[0].n
        self.u       = tcurveNet[0].s
        self.unump   = tcurveNet[0].nump
        self.ubasis  = tcurveNet[0].basis

        self.vtrange = lcurveNet[0].trange
        self.vorder  = lcurveNet[0].k
        self.vknots  = lcurveNet[0].t
        self.vn      = lcurveNet[0].n
        self.v       = lcurveNet[0].s
        self.vnump   = lcurveNet[0].nump
        self.vbasis  = lcurveNet[0].basis
        
        self.tcurveNet  = tcurveNet # transverse curve network (only used for sizing the surface - I think...)
        self.lcurveNet  = lcurveNet # longitudinal curve network
        self.surfNet    = []        # final surface curve network - [transverse[longitudinal[x,y,z]]]
        self.surface    = np.zeros((self.unump,self.vnump,3),float) #storage of points on the surface
        self.surfacept  = surfacept # surface at an arbitrary point

        # cheap (from a coding perspective) way of getting the needed funcs
        # to evaluate surface properties at arbitrary u,v:
        self.uBasisFuns     = tcurveNet[0].BasisFuns 
        self.uDersBasisFunc = tcurveNet[0].DersBasisFunc
        self.uFindSpan      = tcurveNet[0].FindSpan
        
        self.vBasisFuns     = lcurveNet[0].BasisFuns
        self.vDersBasisFunc = lcurveNet[0].DersBasisFunc
        self.vFindSpan      = lcurveNet[0].FindSpan
        
        self.K=np.zeros((len(self.u),len(self.v)),float)    #gaussian curvature
        self.meanCurvature=np.zeros((len(self.u),len(self.v)),float)   #average curvature
        
        
        ##
        ##**************************************************
        ## Curve Coordinates
        Origin = np.asarray([0.,0.,0.])
        x, y, z = np.identity(3)
        self.ccs = Frame( np.asarray([x,y,z]) ) #curve coordinate system
        self.origin = DeepVector(v = Origin,
                                B=self.ccs)
        
        ## axes of rotation
        self._ax_yaw        = yax           # self._ax_LR = yax #(0, -1, 0) #Left Right
        self._ax_alt_yaw_   = (0., 0., 1.)     # self._ax_LR_alt = (0, 0, 1)
        self._ax_pitch      = xax           # self._ax_UD = xax #(1,0,0)
        self._ax_roll       = zax
        
        ## initial conditions:
        self.rots = Quaternion.from_v_theta(z, 0.) #null rotation vector
        # define the current rotation axis and angle
        if current_rot is None:
            self._current_rot = Quaternion.from_v_theta((0., 0., 1.), 0.)
        else:
            self._current_rot = current_rot
        self._current_trans_x = 0.
        self._current_trans_y = 0.
        self._current_trans_z = 0.
        ##
        ##**************************************************
        ##

        # volume kernal:
        self.VolumeKernel   = np.zeros((self.un,self.vn,self.un,self.vn,self.un,self.vn),float)
        self.Vol_kernel1 =  np.zeros((self.un,self.vn,self.un,self.vn),float)
        self.Vol_kernel2 =  np.zeros((self.un,self.vn),float)
        #self.Vol_kernelz =  np.zeros((self.un,self.vn,self.un,self.vn,self.un,self.vn),float)
        self.vol = 0.0
        #self.makeSurfNet()
        
        self.FormParmeters = {}
        
        # AD Surface:
        self.surfNetAD = None        
        
        if kind=='Loft':
            self.LoftSurface()  # replace with __call__  ?
        if kind=='bilinear':
            self.BidirectionSurface()

        self.closed = False
        
    
    
    ##
    ##*****************************************  
    ## 3D transformations
    
    def rotate_surface(self, axis, angle):
        """
        usage:
            to compose a rotar about the z axis, e.g.:
            self.rotate_surface( (0,0,1) , 0.5*np.pi )
        """
        self.rots = Quaternion.from_v_theta(axis,
                                            angle)
        
        self._current_rot = self._current_rot * self.rots
        
        Rs = (self._current_rot * self.rots).as_rotation_matrix()
        
        return self.update_loft(Rs = Rs)
    
    def translate_surface(self, dx=0., dy=0., dz=0.):
        self._current_trans_x = dx
        self._current_trans_y = dy
        self._current_trans_z = dz
        Ts = np.asarray([self._current_trans_x,
                         self._current_trans_y,
                         self._current_trans_z])
        return self.update_loft(Ts = Ts[:])
    
    def reflect_surface(self, point,dual_axis):
        """currently only reflection about
        planes dual to the given system
        are supported 
        
        inputs:
            dual_axis : the axis of the reflection plane
            
        """
        #v0 = np.zeros((4,1),float)
        #v1 = np.zeros((4,1),float)
        v0 = np.zeros(4,float)
        v1 = np.zeros(4,float)
        v0[3] = v1[3] = 1.0
        v0[0:3] = point[:]
        v1[0:3] = dual_axis[:]
        Rs = tfm.reflection_matrix(v0, v1)
        return self.update_loft(Rs = Rs[:3,:3])
    
    def translate_tcurves(self,Ts):
        nnet = []
        for curve in self.tcurveNet:
            vertices = curve.vertices[:,0] + Ts[0]
            vertices = curve.vertices[:,1] + Ts[1]
            vertices = curve.vertices[:,2] + Ts[2]
            curve.vertices = vertices
            nnet.append(curve)
        tcurveNet = nnet
        return tcurveNet
    
    def translate_lcurves(self,Ts):
        nnet = []
        for curve in self.lcurveNet:
            vertices = curve.vertices[:,0] + Ts[0]
            vertices = curve.vertices[:,1] + Ts[1]
            vertices = curve.vertices[:,2] + Ts[2]
            curve.vertices = vertices
            nnet.append(curve)
        lcurveNet = nnet
        return lcurveNet
    
    def rotate_tcurves(self,Rs):
        nnet = []
        for curve in self.tcurveNet:
            vertices = np.dot(curve.vertices,Rs)
            curve.vertices = vertices
            nnet.append(curve)
        tcurveNet = nnet
        return tcurveNet
    
    def rotate_lcurves(self,Rs):
        nnet = []
        for curve in self.lcurveNet:
            vertices = np.dot(curve.vertices[:],Rs)
            curve.vertices = vertices
            nnet.append(curve)
        lcurveNet = nnet
        return lcurveNet
    
    def update_loft(self,Rs=None,Ts=None):
        """TODO:
            make this more efficient
            by similar methods 
            as was done/started for B-spline curves
            i.e. don't instantiate all over again!
        """
        if Rs is not None:
            tcurveNet = self.rotate_tcurves(Rs)
            lcurveNet = self.rotate_lcurves(Rs)
        if Ts is not None:
            tcurveNet = self.translate_tcurves(Rs)
            lcurveNet = self.translate_lcurves(Rs)
        return BsplineSurface(tcurveNet,
                              lcurveNet,
                              kind='Loft')#,
                              #current_rot = self._current_rot)
    
    ##
    ##*****************************************  
    ## project-cuts
    def station_cuts(self, n):
        return
    
    def waterplane_cuts(self, n):
        return
    
    def profile_cuts(self, n):
        return
    ##
    ##*****************************************  
    ## compute
    def define_surfNet(self, surfNet):
        self.surfNet=surfNet
        return
        
    def define_surfNetAD(self, surfNetAD):
        self.surfNetAD = surfNetAD
        return


    def compute_surface(self):
        """assumes surfNet is already in place but shall change"""
        self.surface    = np.zeros((self.unump,self.vnump,3),float) 
        for w1,u in enumerate(self.u):
            for w2,v in enumerate(self.v):
                for i in range(self.un):
                    for j in range(self.vn):
                        self.surface[w1,w2,:] += self.surfNet[i][j]*self.ubasis[w1,i]*self.vbasis[w2,j]
        return
        
#    def compute_surface_general(self):
#        self.surface    = np.zeros((self.unump,self.vnump,3),float) 
#        for w1,u in enumerate(self.u):
#            for w2,v in enumerate(self.v):
#                active_vertices = get_active_vertices(u,v)
#                active_knots = get_active_knots(u,v)
#                self.surface[w1,w2,:] = 
#        return

    def LoftSurface(self):
        """Form a B-Spline Surface wich
            interpolates a transverse curve network
            via the vertices of the longitudinal curve network which interpolate
            some the control vertices of the transverse curves
            -a Harries surface"""

        # establish the surface vertex network:
        """
        if self.surfnet = None:
            for el in self.lcurveNet:
                self.surfNet.append(el.vertices)
            self.surfNet = np.asarray(self.surfNet)
        else:
            self.surfNet = None
            for el in self.lcurveNet:
                self.surfNet.append(el.vertices)
            self.surfNet = np.asarray(self.surfNet)
        #"""
        self.surfNet = []
        for el in self.lcurveNet:
            self.surfNet.append(el.vertices)
        self.surfNet = np.asarray(self.surfNet)
        
        # calculate the points on the surface:
        for w1,u in enumerate(self.u):
            for w2,v in enumerate(self.v):
                for i in range(self.un):
                    for j in range(self.vn):
                        self.surface[w1,w2,:] += self.surfNet[i][j]*self.ubasis[w1,i]*self.vbasis[w2,j]
                        
        return

    def BidirectionSurface(self):
        """
            Form a surface from a bi-directional curve network
            requirements:
            
        """

        tsteps = self.un  #tcurveNet[0].n
        lsteps = self.vn  #lcurveNet[0].n

        
        
        for i in range(tsteps):
            self.surfNet.append([])
            for j in range(lsteps):
                if i==0:
                    l=list(self.lcurveNet[0].vertices[j])
                    self.surfNet[i].append(list(l))
                elif j==0:
                    l=list(self.tcurveNet[0].vertices[i])
                    self.surfNet[i].append(list(l))
                elif i==lsteps-1:
                    l=list(self.lcurveNet[1].vertices[j])
                    self.surfNet[i].append(list(l))
                elif j==tsteps-1:
                    l=list(self.tcurveNet[1].vertices[i])
                    self.surfNet[i].append(list(l))
        
                
                else:
                    #fix with sequence of "linear curve generator"
                    l=[ (((tsteps-i-1.)/(tsteps-1.))*self.tcurveNet[0].vertices[i])/2.\
                        +((1.-((tsteps-i-1.)/(tsteps-1.)))*self.tcurveNet[1].vertices[i])/2.\
                        +(((lsteps-j-1.)/(lsteps-1.))*self.lcurveNet[0].vertices[j])/2.\
                        +((1.-((lsteps-j-1.)/(lsteps-1.)))*self.lcurveNet[1].vertices[j])/2. ]
                    self.surfNet[i].append(list(l[0]))
        
            self.surfNet[i] = np.asarray(self.surfNet[i])
            #print i,j,self.surfNet
                
        self.surfNet = np.asarray(self.surfNet)

        #print self.sufnet        
        for w1,u in enumerate(self.u):
            for w2,v in enumerate(self.v):
                for i in range(self.un):
                    for j in range(self.vn):
                        #pass
                        #print i,j
                        self.surface[w1,w2,:] += 0.+\
                            self.surfNet[i][j]*\
                            self.ubasis[w1,i]*\
                            self.vbasis[w2,j]
        
        return

    

    def calc_K(self):
        if self.surfNetAD is None:
            adsurf = surfAD(self)
            self.define_surfNetAD(adsurf)
        else:
            pass
        for w1,u in enumerate(self.u):
            for w2,v in enumerate(self.v):
                a = surf_gaussian_curvature_quick(self,u,v)
                #print surf_gaussian_curvature(self,u,v)
                #print a
                #print a.value
                #print a.der
                self.K[w1,w2]=a.value
        return
    

    def calc_meanCurvature(self):
        for w1,u in enumerate(self.u):
            for w2,v in enumerate(self.v):
                a = surf_mean_curvature(self,u,v)
                self.meanCurvature[w1,w2]=a.value
        return

    def eval_K(self, u,v):
        #a = AD object
        a = surf_gaussian_curvature(self,u,v)
        return a.value

    def eval_meanCurvature(self, u,v):
        a = surf_mean_curvature(self,u,v)
        return a.value

    


    def evalSurface(self, u,v):
        """
            Evaluate a suface at an arbitrary u,v location on the surface
                -this demonstrates why the knot vector must be
                constant for all u
                and seperately
                constant for all v """

        if (u<=1. and u>=0. and v<=1. and v>=0.):
            tcurveNet = self.tcurveNet
            lcurveNet = self.lcurveNet

            
            # compute basis functions in u:
            uN      = np.zeros((tcurveNet[0].n),float) 
            uspan   = tcurveNet[0].FindSpan(u)
            ubasis  = self.uBasisFuns(uspan,u,uN[uspan-tcurveNet[0].p:uspan+1])
            
            # compute basis functions in v:
            vN      = np.zeros((lcurveNet[0].n),float)
            vspan   = lcurveNet[0].FindSpan(v)
            vbasis  = self.vBasisFuns(vspan,v,vN[vspan-lcurveNet[0].p:vspan+1])

            # compute the surface at u,v:
            surf_point = 0.
            for i in range(self.tcurveNet[0].n):
                for j in range(self.lcurveNet[0].n):
                    surf_point += self.surfNet[i][j]*uN[i]*vN[j]
        else:
            surf_point = np.zeros((3),float)
        return surf_point
    
    def FindBoundaryLoc(self, pt):
        """return the location 
        on the u boundary and
        on th v boundaries 
        where the pt can be coordinated from on the interior.
        """
        return


    def findPtloc(self, pt, s, uv_search_Index, r3_index):
        """Function to find the v parameterization on a surface, given:
                -pt location in real space which we want to locate the parameterization of
                -s location in parameter space

            inputs:
                surface
                pt x or y or z location
                u or v location
                uv_search_Index = [0,1]   maps => [u,v]
                    0 to search for u
                    1 to search for v
                r3_index = [0,1,2] maps => [x,y,z]
                    0 to search using x
                    1 to search using y
                    2 to search using z

            outputs:
                mid     = u or v location
                q       = [x,y,z] surf evaluation
                found   = True or False
                


            An example:

                say you want to find the v location
                where:
                x=.45
                u=.4

                then:
                pt=.45
                s=.4
                uv_search_Index=1
                r3_index=0
        """

        span = [[self.utrange[0],self.utrange[-1]],[self.vtrange[0],self.vtrange[-1]]]

        # algorithm: bianary search
        # range in either u or v:
        high = span[uv_search_Index][-1]
        low  = span[uv_search_Index][0]
        mid  = (low + high)/2. #the initial guess

        search = [s]
        search.insert(uv_search_Index,mid)
        
        q=self.evalSurface(search[0],search[1]) #[x,y,z] location

        count=0
        while((abs(q[r3_index]-pt))>.00000001 and count<101):
            if pt<q[r3_index]:
                high = mid
            else:
                low = mid
            mid = (low+high)/2.
            search = [s]
            search.insert(uv_search_Index,mid)
            q=self.evalSurface(search[0],search[1]) 
            count += 1
        if count > 100:
            print 'error, point not found'
            found = False
            #return mid, q, found
        else:
            found = True
            #return mid, q, found    

        #mid = param loc
        #q = point in R2 or R3
        #found = True or False
        # count = number of iterations required
        return mid, q, found, count

    #    def findPt(self, pt):
    #        px = pt[0]
    #        py = pt[1]
    #        pz = pt[2]
    #        span = [[self.utrange[0],self.utrange[-1]],[self.vtrange[0],self.vtrange[-1]]]
    #        
    #        for uv_search_Index in [0,1]:
    #            high = span[uv_search_Index][-1]
    #            low  = span[uv_search_Index][0]
    #            mid  = (low + high)/2.
    #            
    #            search = [s]
    #            search.insert(uv_search_Index,mid)            
    #            q=self.evalSurface(search[0],search[1])
            
            
            

    def ComputeVolumeKernel(self):
        """Function to return the 6 dimensional
            matrix of integrated basis products
            intrinsic to fast volume calculation
            of the area enclosed by b-spline surfaces

            inputs:
                -surface
                
            outputs:
                -integrated volume kernal matrix
        """
        #order = (self.uorder-1+self.vorder-1)+\
        #        (self.uorder-2+self.vorder-1+self.uorder-1+self.vorder-2)+1
        
        order_u = (self.uorder-1)+\
                    (self.uorder-2+self.uorder-1)+1
        order_v = (self.vorder-1)+\
                    (self.vorder-1+self.vorder-2)+1
        #self.args = self.uorder,self.vorder,self.t
        order = [order_u,order_v]

        for iz in range(self.un):
            #print 'iz :::::::::::: >', iz
            for jz in range(self.vn):
                #print 'jz --->', jz
                for iy in range(self.un):
                    #print iy#'iy --',iy
                    for jy in range(self.vn):
                        #print 'jy', jy
                        for ix in range(self.un):
                            #print ix
                            for jx in range(self.vn):
                                self.VolumeKernel[ix,jx,iy,jy,iz,jz] = VolKernelQuatdrature(self,order,ix,jx,iy,jy,iz,jz)

        #self.VolumeKernel=self.VolumeKernel-np.transpose(self.VolumeKernel) #Stephan Harries Formulation
        
        return


#    def ComputeVolumeKernal(self):
#        """Function to return the 6 dimensional
#            matrix of integrated basis products
#            intrinsic to fast volume calculation
#            of the area enclosed by b-spline surfaces
#
#            inputs:
#                -surface
#                
#            outputs:
#                -integrated volume kernal matrix
#        """
#        order = (self.uorder-1+self.vorder-1)+1+\
#                (self.uorder-2+self.vorder-1+self.uorder-1+self.vorder-2)+1
#        
#
#        #self.args = self.uorder,self.vorder,self.t
#        
#
#        for iz in range(self.un):
#            print 'iz :::::::::::: >', iz
#            for iy in range(self.un):
#                print iy#'iy --',iy
#                for ix in range(self.un):
#                    print ix
#                    self.VolumeKernel[ix,:,iy,:,iz,:] = VolKernelQuatdrature(self,order,ix,jx,iy,jy,iz,jz)
#
#
#
#        for jz in range(self.vn):
#            print 'jz --->', jz
#            for jy in range(self.vn):
#            #print 'jy', jy
#                for jx in range(self.vn):
#                    self.VolumeKernel[:,jx,:,jy,:,jz] = VolKernelQuatdrature(self,order,ix,jx,iy,jy,iz,jz)
#
#        #self.VolumeKernel=self.VolumeKernel-np.transpose(self.VolumeKernel) #Stephan Harries Formulation
#        
#        return

    def CalculateSurfaceVolume(self):
        """Function to compute the volume enclosed by the surface(s)"""

        
        # Compute the momment matrices
        #up = self.uorder-1
        #vp = self.vorder-1
        #N2 = up + vp + ( (up-1+vp) + (up+vp-1) + (0) )
        
        #t2 = make_knot_vector((N2),curve.n-1,curve.trange,type='open')
        #self.args2 = self.k,t2,1 
        
        ## surf1.surfNet[jx,:][:,0] longi x
        ## surf1.surfNet[jy,:][:,1]  longi y
        ## surf1.surfNet[jz,:][:,2]  longi z
        
        ## surf1.surfNet[:,ix][:,0] tran x
        ## surf1.surfNet[:,iy][:,1] tran y
        ## surf1.surfNet[:,iz][:,2] tran z
        
        #for iz in range(self.un):
        self.vol=0.
#        for jz in range(self.vn):
#            #print 'iz :::::::::::: >', iz
#            #for jz in range(self.vn):
#            for ix in range(self.un):
#                #print 'jz --->', jz
#                #for iy in range(self.un):
#                for jx in range(self.vn):
#                    #print iy#'iy --',iy
#                    #for jy in range(self.vn):
#                    for iy in range(self.un):
#                        #print 'jy', jy
#                        #for ix in range(self.un):
#                        for jy in range(self.vn):
#                            #print ix
#                            #for jx in range(self.vn):
#                            for iz in range(self.un):
#                                self.vol += \
#                                np.dot( 
#                                np.dot( self.surfNet[:,ix][:,0] ,
#                                np.dot( 
#                                np.dot( self.surfNet[:,iy][:,1] ,
#                                np.dot( 
#                                np.dot( self.surfNet[:,iz][:,2] ,\
#                                self.VolumeKernel[:,:,:,:,:,:]),
#                                self.surfNet[jy,:][:,1])),
#                                self.surfNet[jx,:][:,0])),
#                                self.surfNet[jz,:][:,2] )

        for jz in range(self.vn):
            for ix in range(self.un):
                for jx in range(self.vn):
                    for iy in range(self.un):
                        for jy in range(self.vn):
                            for iz in range(self.un):
                                self.vol += \
                                self.surfNet[ix][jx,0] * \
                                self.surfNet[iy][jy,1] * \
                                self.surfNet[iz][jz,2] * \
                                self.VolumeKernel[ix,jx,iy,jy,iz,jz] #* \
                                #self.surfNet[iy,:][jy,1] * \
                                #self.surfNet[ix,:][jx,0] * \
                                #self.surfNet[iz,:][jz,2] 
                                
#        self.vol = np.dot( self.surfNet[jz][:,2]   ,
#                          np.dot( self.surfNet[:,ix][:,0] ,
#                                 np.dot( self.surfNet[jx][:,0]   ,
#                                        np.dot( self.surfNet[:,iy][:,1] ,
#                                               np.dot( self.surfNet[jy][:,1]   , 
#                                                      np.dot( self.surfNet[:,iz][:,2] ,self.VolumeKernel[:,:,:,:,:,:]))))))
#                                                      
#        self.vol = np.dot( self.surfNet[jz][:,2]   ,
#                          np.dot( self.surfNet[:,ix][:,0] ,
#                                 np.dot( self.surfNet[jx][:,0]   ,
#                                        np.dot( self.surfNet[:,iy][:,1] ,
#                                               np.dot( self.surfNet[jy][:,1]   , 
#                                                      np.dot( self.surfNet[:,:][:,2] ,self.VolumeKernel[:,:,:,:,:,:]))))))
#
#        self.vol = np.dot( np.dot( self.surfNet[:,ix][:,0],np.dot( np.dot( self.surfNet[:,iy][:,1],
#                                  np.dot(  np.dot( self.surfNet[:,iz][:,2],self.VolumeKernel[:,:,:,:,:,:]),
#                                         self.surfNet[jz][:,2])) ,self.surfNet[jy][:,1] )),self.surfNet[jx][:,0])
#        
#        self.vol = np.dot( np.dot( self.surfNet[:,:][:,0],np.dot( np.dot( self.surfNet[:,:][:,1],
#                                  np.dot(  np.dot( self.surfNet[:,:][:,2],self.VolumeKernel[:,:,:,:,:,:]),
#                                         self.surfNet[:][:,2])) ,self.surfNet[:][:,1] )),self.surfNet[:][:,0])

#        for iz in range(self.un):
#            for jz in range(self.vn):
#                for iy in range(self.un):
#                    for jy in range(self.vn):
#                        self.Vol_kernel1[iy,jy,iz,jz]=np.dot( self.surfNet[:,iz][:,0] ,np.dot(self.VolumeKernel[:,:,iy,jy,iz,jz],self.surfNet[jy,:][:,0]))

        return


    def closeHull(self):
        
        tbcrv=[tCurveNet[0],tCurveNet[-1]]
        lbcrv=[lCurveNet[0],lCurveNet[-1]]

        #Parameterize the space:

        #Distance from aft t curve to fwd t curve:
        xtd = []
        ytd = []
        ztd = []
        for i in range(tCurveNet[0].n):
            xtd.append(tCurveNet[-1].vertices[i,0]-tCurveNet[0].vertices[i,0])
            ytd.append(tCurveNet[-1].vertices[i,1]-tCurveNet[0].vertices[i,1])
            ztd.append(tCurveNet[-1].vertices[i,2]-tCurveNet[0].vertices[i,2])

        #Step size in each direction, at each point on the param:
        dx = [] #dx[0] is the step size of the keel points, dx[-1] step size of wL
        dy = []
        dz = []
        loc_n = lCurveNet[0].n
        for i in range(loc_n):
            dx.append(xtd[i]/(loc_n - 2))
            dy.append(ytd[i]/(loc_n - 2))
            dz.append(ztd[i]/(loc_n - 2))
        

        #make new transverse nets along the parameterization
        ticrv = [np.copy(tCurveNet[0])]
        


        #Insert knots at the corner to sperate the pieces:

        
        
        return




    def plotSurface(self, 
                    innerBox=None, 
                    OuterBox=None, 
                    limx=None,limy=None,limz=None):
        """
            Plot a surface
        """
        
        tCurveNet=self.tcurveNet
        lCurveNet=self.lcurveNet
        surfnet=self.surfNet
        fig = plt.figure(figsize=(9.5,5.0))
        #---- First subplot
        rect = fig.add_subplot(1, 1, 1).get_position()
        ax = Axes3D(fig, rect)
        #ax.plot_wireframe(self.surface[:,:,0], self.surface[:,:,1], self.surface[:,:,2], rstride=10, cstride=10, color = 'r')
        ax.plot_surface(self.surface[:,:,0], 
                        self.surface[:,:,1], 
                        self.surface[:,:,2],  
                        rstride=20, 
                        cstride=5, 
                        color='r',
                        alpha=.3)
        #ax.plot_surface(self.surface[:,:,0], -self.surface[:,:,1], self.surface[:,:,2],  rstride=20, cstride=5, color='r')

        #"""
        #"""
        #ax.plot(tnet[:,0],tnet[:,1],0.)
        for item in tCurveNet:  
            ax.plot(item.r[:,0], 
                    item.r[:,1], 
                    item.r[:,2], label='parametric curve', color = 'b')
            ax.plot(item.vertices[:,0],
                    item.vertices[:,1],
                    item.vertices[:,2],
                    marker = "o",  linestyle = "--", color = '0.75')
        for item in lCurveNet:  
            ax.plot(item.r[:,0], 
                    item.r[:,1], 
                    item.r[:,2], label='parametric curve', color = 'g')
            ax.plot(item.vertices[:,0],
                    item.vertices[:,1],
                    item.vertices[:,2],
                    marker = "x",  linestyle = "--", color = '0.75')
        for item in surfnet:
            ax.plot(item[:,0],
                    item[:,1],
                    item[:,2],
                    marker = "x",  linestyle = "--", color='g')
        if innerBox:
            r=innerBox[0]
            t=innerBox[1]
            u=innerBox[2]
            for s, e in combinations(np.array(list(product(r,t,u))), 2):
                if np.sum(np.abs(s-e)) == r[1]-r[0] or np.sum(np.abs(s-e)) == t[1]-t[0] or np.sum(np.abs(s-e)) == u[1]-u[0]:
                    ax.plot3D(*zip(s,e), color="b")
        if OuterBox:
            r=OuterBox[0]
            t=OuterBox[1]
            u=OuterBox[2]
            for s, e in combinations(np.array(list(product(r,t,u))), 2):
                if np.sum(np.abs(s-e)) == r[1]-r[0] or np.sum(np.abs(s-e)) == t[1]-t[0] or np.sum(np.abs(s-e)) == u[1]-u[0]:
                    ax.plot3D(*zip(s,e), color="b")
            
        #"""
#        if limx:
#            ax.set_xlim3d(0., limx)
#            ax.set_ylim3d(0.,limy)
#            ax.set_zlim3d(0.,limz)
        axmax = max(curve.extremes1d(0))
        aymax = max(curve.extremes1d(1))
        azmax = max(curve.extremes1d(2))
        biggest = max(axmax,aymax,azmax)
        #
        ax.set_xlim3d(-biggest*.0,1.*biggest)
        ax.set_ylim3d(-biggest*.5,.5*biggest)
        ax.set_zlim3d(-biggest*.5,.5*biggest)
        plt.show()
        return


    def plotMultiSurface(self, surfList=[], limx=None,limy=None,limz=None):
        """
            Plot surfaces together
        """
        if surfList==[]:
            tCurveNet=self.tcurveNet
            lCurveNet=self.lcurveNet
            surfnet=self.surfNet
            fig = plt.figure(figsize=(9.5,5.0))
            #---- First subplot
            rect = fig.add_subplot(1, 1, 1).get_position()
            ax = Axes3D(fig, rect)
            #ax.plot_wireframe(self.surface[:,:,0], self.surface[:,:,1], self.surface[:,:,2], rstride=10, cstride=10, color = 'r')
            ax.plot_surface(self.surface[:,:,0], self.surface[:,:,1], self.surface[:,:,2],  rstride=20, cstride=5, color='r')
            ax.plot_surface(self.surface[:,:,0], -self.surface[:,:,1], self.surface[:,:,2],  rstride=20, cstride=5, color='r')

            #"""
            #"""
            #ax.plot(tnet[:,0],tnet[:,1],0.)
            for item in tCurveNet:  
                ax.plot(item.r[:,0], item.r[:,1], item.r[:,2], label='parametric curve', color = 'b')
                ax.plot(item.vertices[:,0],item.vertices[:,1],item.vertices[:,2],marker = "o",  linestyle = "--", color = '0.75')
            for item in lCurveNet:  
                ax.plot(item.r[:,0], item.r[:,1], item.r[:,2], label='parametric curve', color = 'g')
                ax.plot(item.vertices[:,0],item.vertices[:,1],item.vertices[:,2],marker = "x",  linestyle = "--", color = '0.75')
            for item in surfnet:
                ax.plot(item[:,0],item[:,1],item[:,2],marker = "x",  linestyle = "--", color='g')
            
            #"""
            if limx:
                ax.set_xlim3d(-limx, limx)
                ax.set_ylim3d(-limy,limy)
                ax.set_zlim3d(-limz,limz)
            #plt.show()
        else:
            fig = plt.figure(figsize=(9.5,5.0))
            rect = fig.add_subplot(1, 1, 1).get_position()
            ax = Axes3D(fig, rect)
            for surf in surfList:
                tCurveNet=surf.tcurveNet
                lCurveNet=surf.lcurveNet
                surfnet=surf.surfNet
                
                #---- First subplot
                
                #ax.plot_wireframe(surface.surface[:,:,0], surface.surface[:,:,1], surface.surface[:,:,2], rstride=10, cstride=10, color = 'r')
                ax.plot_surface(surf.surface[:,:,0], surf.surface[:,:,1], surf.surface[:,:,2],  rstride=20, cstride=5, color='r')
                #ax.plot_surface(surf.surface[:,:,0], -surf.surface[:,:,1], surf.surface[:,:,2],  rstride=20, cstride=5, color='r')

                #"""
                #"""
                #ax.plot(tnet[:,0],tnet[:,1],0.)
                for item in tCurveNet:  
                    ax.plot(item.r[:,0], item.r[:,1], item.r[:,2], label='parametric curve', color = 'b')
                    ax.plot(item.vertices[:,0],item.vertices[:,1],item.vertices[:,2],marker = "o",  linestyle = "--", color = '0.75')
                for item in lCurveNet:  
                    ax.plot(item.r[:,0], item.r[:,1], item.r[:,2], label='parametric curve', color = 'g')
                    ax.plot(item.vertices[:,0],item.vertices[:,1],item.vertices[:,2],marker = "x",  linestyle = "--", color = '0.75')
                for item in surfnet:
                    ax.plot(item[:,0],item[:,1],item[:,2],marker = "x",  linestyle = "--", color='g')
            if limx:
                ax.set_xlim3d(-limx, limx)
                ax.set_ylim3d(-limy,limy)
                ax.set_zlim3d(-limz,limz)
        plt.show()
        return       

#-----------------------------------------------------------------------------------------
class TsplineSurface(BsplineSurface):
    """
        Efficient data structures for T-spline modeling
        
        Ache, C. and Berkhan, V.
        Leibniz University
        
    """

    def __init__(self, tcurveNet, lcurveNet, surfNet=None, kind='Loft',surface=None, surfacept=None, dim=3):
        """
        init assigns the attributes
        of the B-Spline Surface

        - Quick and Dirty way is relient on the curves which make up the surface

        Loops:
             innermost  : dimension
                middle  : longitudinal
                 outer  : transverse
        ie:
            [u,v][dim]
        """
        self.dim    = dim # 3D vertices are expected
        
        self.utrange = tcurveNet[0].trange
        self.uorder  = tcurveNet[0].k
        self.uknots  = []#tcurveNet[0].t
        self.un      = tcurveNet[0].n
        self.u       = tcurveNet[0].s
        self.unump   = tcurveNet[0].nump
        self.ubasis  = []#tcurveNet[0].basis

        self.vtrange = lcurveNet[0].trange
        self.vorder  = lcurveNet[0].k
        self.vknots  = []#lcurveNet[0].t
        self.vn      = lcurveNet[0].n
        self.v       = lcurveNet[0].s
        self.vnump   = lcurveNet[0].nump
        self.vbasis  = []#lcurveNet[0].basis
        
        self.tcurveNet  = tcurveNet # transverse curve network (only used for sizing the surface - I think...)
        self.lcurveNet  = lcurveNet # longitudinal curve network
        self.surfNet    = []        # final surface curve network - [transverse[longitudinal[x,y,z]]]
        self.knotNet    = []        # suface knot system
        self.surface    = np.zeros((self.unump,self.vnump,3),float) #storage of points on the surface
        self.surfacept  = surfacept # surface at an arbitrary point

        # cheap (from a coding perspective) way of getting the needed funcs
        # to evaluate surface properties at arbitrary u,v:
        self.uBasisFuns     = tcurveNet[0].BasisFuns 
        self.uDersBasisFunc = tcurveNet[0].DersBasisFunc
        self.uFindSpan      = tcurveNet[0].FindSpan
        
        self.vBasisFuns     = lcurveNet[0].BasisFuns
        self.vDersBasisFunc = lcurveNet[0].DersBasisFunc
        self.vFindSpan      = lcurveNet[0].FindSpan
        
        self.K=np.zeros((len(self.u),len(self.v)),float)    #gaussian curvature
        self.meanCurvature=np.zeros((len(self.u),len(self.v)),float)   #average curvature

        # volume kernal:
        self.VolumeKernel   = np.zeros((self.un,self.vn,self.un,self.vn,self.un,self.vn),float)
        self.Vol_kernel1 =  np.zeros((self.un,self.vn,self.un,self.vn),float)
        self.Vol_kernel2 =  np.zeros((self.un,self.vn),float)
        #self.Vol_kernelz =  np.zeros((self.un,self.vn,self.un,self.vn,self.un,self.vn),float)
        self.vol = 0.0
        #self.makeSurfNet()
        if kind=='Loft':
            self.LoftSurface()  # replace with __call__  ?
        elif kind=='bilinear':
            self.BidirectionSurface()
        else:
            pass
            #self.compute_tsurface()
        self.closed = False
        
        
    def define_knots(self):
        for curve in self.tcurveNet:
            self.uknots.append(curve.t)
        self.uknots = np.asarray(self.uknots)
        for curve in self.lcurveNet:
            self.vknots.append(curve.t)
        self.vknots = np.asarray(self.vknots)
        return
        
        
    def define_basis(self):
        # assumes all curves have same nump as each surface curve
        for curve in self.tcurveNet:
            self.ubasis.append(curve.basis)
        self.ubasis = np.asarray(self.ubasis)
        for curve in self.lcurveNet:
            self.vbasis.append(curve.basis)
        self.vbasis = np.asarray(self.vbasis)
        return
    
    
    def define_surfNet(self, surfNet):
        self.surfNet=surfNet
        return


    def compute_surface(self):
        """assumes surfNet is already in place but shall change"""
        self.surface    = np.zeros((self.unump,self.vnump,3),float) 
        for w1,u in enumerate(self.u):
            for w2,v in enumerate(self.v):
                for i in range(self.un):
                    for j in range(self.vn):
                        self.surface[w1,w2,:] += self.surfNet[i][j]*self.ubasis[w1,i]*self.vbasis[w2,j]
        return
        
        
    def compute_surfNet(self):
        # establish the surface vertex network:
        for el in self.lcurveNet:
            self.surfNet.append(el.vertices)
        self.surfNet = np.asarray(self.surfNet)
        return
        
        
    def compute_tsurface(self):
        """assumes surfNet is already in place but shall change"""
        self.compute_surfNet()
        self.surface    = np.zeros((self.unump,self.vnump,3),float) 
        W1 = len(self.tcurveNet)
        W2 = len(self.lcurveNet)
        
        for w1, tcurve in zip(range(W1),self.tcurveNet):
            un = tcurve.n
            ut = tcurve.t
            ubasis = tcurve.basis
            print 'tcuve stuff = ', un
            print ubasis
            for w2, lcurve in zip(range(W2),self.lcurveNet):
                vn = lcurve.n
                vt = lcurve.t
                vbasis = lcurve.basis
                print 'lcuve stuff = ', vn
                print vbasis
                for i in range(un):
                    for j in range(vn):
                        print 'w1,w2,i,j = ', w1,w2,i,j
                        print 'self.surfNet[i][j] =',self.surfNet[i][j]
                        print 'ubasis[w1,i] =',ubasis[w1,i]
                        print 'vbasis[w2,j] =',vbasis[w2,j]
                        self.surface[w1,w2,:] += self.surfNet[i][j]*ubasis[w1,i]*vbasis[w2,j]
        return

    def LoftSurface(self):
        """Form a B-Spline Surface wich
            interpolates a portion (or all) of the transverse curve network
            via the vertices of the longitudinal curve network which interpolate
            some (or all) of the transverse vertices network
            -a Harries surface"""

        # establish the surface vertex network:
        for el in self.lcurveNet:
            self.surfNet.append(el.vertices)
            self.knotNet.append(el.t)
        self.surfNet = np.asarray(self.surfNet)
        self.knotNet = np.asarray(self.knotNet)

        # calculate the points on the surface:
        for w1,u in enumerate(self.u):
            for w2,v in enumerate(self.v):
                for i in range(self.un):
                    for j in range(self.vn):
                        self.surface[w1,w2,:] += self.surfNet[i][j]*self.ubasis[w1,i]*self.vbasis[w2,j]
                        
        return
    
    def plotSurface(self, limx=None,limy=None,limz=None):
        """
            Plot a surface
        """
        tCurveNet=self.tcurveNet
        lCurveNet=self.lcurveNet
        surfnet=self.surfNet
        fig = plt.figure(figsize=(9.5,5.0))
        #---- First subplot
        rect = fig.add_subplot(1, 1, 1).get_position()
        ax = Axes3D(fig, rect)
        ax.plot_surface(self.surface[:,:,0], self.surface[:,:,1], self.surface[:,:,2],  rstride=20, cstride=5, color='r')
        
        #ax.plot_surface(self.surface[:,:,0], -self.surface[:,:,1], self.surface[:,:,2],  rstride=20, cstride=5, color='r')

        #"""
        #"""
        #ax.plot(tnet[:,0],tnet[:,1],0.)
        for item in tCurveNet:  
            ax.plot(item.r[:,0], item.r[:,1], item.r[:,2], label='parametric curve', color = 'b')
            ax.plot(item.vertices[:,0],item.vertices[:,1],item.vertices[:,2],marker = "o",  linestyle = "--", color = '0.75')
        for item in lCurveNet:  
            ax.plot(item.r[:,0], item.r[:,1], item.r[:,2], label='parametric curve', color = 'g')
            ax.plot(item.vertices[:,0],item.vertices[:,1],item.vertices[:,2],marker = "x",  linestyle = "--", color = '0.75')
        #for item in surfnet:
        #    ax.plot(item[:,0],item[:,1],item[:,2],marker = "x",  linestyle = "--", color='g')

        if limx:
            ax.set_xlim3d(0., limx)
        if limy:
            ax.set_ylim3d(0.,limy)
        if limz:
            ax.set_zlim3d(0.,limz)
        plt.show()
        return
#-----------------------------------------------------------------------------------------
def main():
    
    return 
    
    

def make_longitudinal_curve_net(tCurveNet):
    lCurveNet = []
    for j in range(tCurveNet[0].n):
        lCurveNet.append(construct_jth_curve_from_net(tCurveNet, 0, len(tCurveNet), j))
    return lCurveNet
    
def testTspline():
    xmin = 0.
    xmax = 50.
    ymin = 0.
    ymax = 5.
    zmax = 0.
    tc1 = np.asarray([[0.,0.,0.],[10.,0.,0.],[20.,0.,0.],[30.,0.,0.],[40.,0.,0.],[50.,0.,0.]])
    tc2 = np.asarray([[0.,1.,0.],[20.,1.,0.],[30.,1.,0.],[40.,1.,0.],[50.,1.,0.]])
    tc3 = np.asarray([[0.,2.,0.],[20.,2.,0.],[30.,2.,0.],[40.,2.,0.],[50.,2.,0.]])
    tc4 = np.asarray([[0.,3.,0.],[10.,3.,0.],[20.,3.,0.],[40.,3.,0.],[50.,3.,0.]])
    tc5 = np.asarray([[0.,4.,0.],[10.,4.,0.],[20.,4.,0.],[30.,4.,0.],[40.,4.,0.],[50.,4.,0.]]) 
    tc6 = np.asarray([[0.,5.,0.],[10.,5.,0.],[20.,5.,0.],[50.,5.,0.]]) 
    tnet = np.asarray([tc1,tc2,tc3,tc4,tc5, tc6])  

    lc1 = tnet[0]
    lc2 = tnet[1] 
    lc3 = tnet[2] 
    lc4 = tnet[3]
    lc5 = tnet[4]
    lc6 = tnet[5]
    lnet = np.asarray([lc1,lc2,lc3,lc4,lc5,lc6])
    
    k=4
    nump = 50
    
    tCurveNet = []
    for array in tnet:
        tCurveNet.append(Bspline(array,k,nump))

    lCurveNet = []
    for array in lnet:
        lCurveNet.append(Bspline(array,k,nump))       
       
    surf = TsplineSurface(tCurveNet, lCurveNet, kind = 'TSpline')
    #surf.plotSurface(None,None,xmax,ymax,10.)
    surf.plotSurface(xmax,ymax,10.)
    
    return surf
    

    
if __name__ == '__main__':
    
    
    maxsurf = False
    k=4
    nump=300
    vertices=np.asarray([   [-4.,0.], [-2.,0.],[0.,0.],
                                [2.,0.],
                                [3.0,0.],
                                [6.0,5.],
                                [9.0,12.],
                                [11.,12.],
                                [12.,12.],[14.,12.],[16.,12.]])
    v0 = Bspline(vertices, k, nump)
    curve = copy.deepcopy(v0)

    vertices=np.array([[0.,0.],[1.,1.],[2.,2.],[3.,3.],[4.,4.]])
    a = Bspline(vertices,4,50)
    b1,b2=a.CurveSplit(.5)
    curve1 = b1
    curve2 = b2
    ta = 0.
    points = curve2.vertices
    k = curve2.k
    t = curve2.t
    #self = curve1
    
    vertices=np.asarray([   [-4.,0.,0.], [-2.,0.,0.],[0.,0.,0.],
                                [2.,0.,0.],
                                [3.0,0.,0.],
                                [6.0,5.,0.],
                                [9.0,12.,0.],
                                [11.,12.,0.],
                                [12.,12.,0.],[14.,12.,0.],[16.,12.,0.]])
    v0 = Bspline(vertices, k, nump)
    #v0.plot3DmultiList([v0],[])
    self = v0
    index = 7
    curve = self
    v0.knot_insertion(.099)
    #v0.plot3DmultiList([v0],[])
    #v0.knot_removal(4)
       

    #plt.close()
    def read_maxsurf_curves(filename = 'simple_surf.txt'):
        maxsurf_file  = open(filename)
        lines =  maxsurf_file.readlines()
        maxsurf_file.close()
        data = []
        curves = {}
        for line in lines:
            data.append(line.split())
        for el in data:
            curve_name = str(el[0] + el[1])
            if curve_name in curves:
                pass
            else:
                curves[curve_name] = []
            curves[curve_name].append([ float(a) for a in el[3:6]])
        for curve in curves:
            curves[curve] = np.asarray(curves[curve])
        return curves
            
    def make_maxsurf_curves(sets_of_curve_points, k=4, nump = 30):
        curve_list = []
        for curve in sets_of_curve_points:
            curve_list.append(Bspline(sets_of_curve_points[curve],k, nump))
        return curve_list
        
    def swap_indices(curve,i,j):
        if isinstance(curve, Bspline):
            new_vertices = copy.deepcopy(curve.vertices)
            new_vertices[:,i] = curve.vertices[:,j]
            new_vertices[:,j] = curve.vertices[:,i]
        elif isinstance(curve, np.ndarray):
            new_vertices = copy.deepcopy(curve)
            new_vertices[:,i] = curve[:,j]
            new_vertices[:,j] = curve[:,i]
        return new_vertices
        
    def reverse_curve(curve):
        if isinstance(curve, Bspline):
            new_vertices = []
            for vertex in reversed(curve.vertices):
                new_vertices.append(vertex)
            new_vertices = np.asarray(new_vertices)
            new_curve = Bspline(new_vertices, curve.k, curve.nump)
        return new_curve
    
    def flip_sign(vertices,i):
        vertices[:,i] = -vertices[:,i]
        return vertices
        
    """      
    ## curve 1-------------
    curves = read_maxsurf_curves()
    for curve in curves:
        curves[curve] = swap_indices(curves[curve],1,2)
    curves = make_maxsurf_curves(curves)
    #curves[1] = reverse_curve(curves[1])
    curve_set = curveSet(curves, check_order = True)
    
    for curve in curves:
        curve.MomentMatrices()
        curve.curvearea()
        print curve.area
        
    curve_set.compute_area()
    print 'curve_set area = ',curve_set.area
    
    #-----Curve 2
    curves = read_maxsurf_curves('simple_surf_4.txt')
    for curve in curves:
        curves[curve] = swap_indices(curves[curve],0,1)
        curves[curve] = swap_indices(curves[curve],1,2)
        #curves[curve] = flip_sign(curves[curve],0)
    curves = make_maxsurf_curves(curves)
    #curves[1] = reverse_curve(curves[1])
    curve_set = curveSet(curves, check_order = False)
    #curves[2] = reverse_curve(curves[2])
    for curve in curves:
        curve.MomentMatrices()
        curve.curvearea()
        print curve.area
        
    curve_set.compute_area()
    print 'curve_set area = ',curve_set.area
    #"""
    k=4
    nump=30
    l0=linear_vertices((0.,12.),(12.,0.),7)
    l0 = Bspline(l0,k=4,nump=100)
    l1=linear_vertices((0.,12.),(12.,0.),8)
    l1 = Bspline(l1,k=4,nump=100)
    
    l0t = copy.deepcopy(l0)
    l2 = linear_vertices((0.,12.),(12.,0.),14)
    l2 = Bspline(l2,k=4,nump=100)
    
    l0t.knot_insertion(.125)
    l0t.knot_insertion(.375)
    l0t.knot_insertion(.625)
    l0t.knot_insertion(.875)
    
    l0.plot_knots(y=-.1)
    l0t.plot_knots(y=0.)
    l0t.plotbasis(which=0)
    
    a=0.25
    b=0.50
    c=0.75
    
    start = [0.,12.]
    end = [12.,0.]
    num = 10
    vv = curvalinear_vertices(start,end,num)
    
    vb = Bspline(vv,k,nump)
    
    l0.verbose=True
    plt.close()
    l0.plotcurve_detailed()
    self = l0
    
    
    if False: #check ability to accept a knot vector on start
        k=4
        vv = curvalinear_vertices(start,end,5)
        vb = Bspline(vv,k,nump)
        to=copy.copy(vb.t)
        to[4]=.51
        vc = Bspline(vv,k,nump,t=to)
    
    import wavelet_matrices as wm
    help(wm)
    
    c1 = curve1
    P = wm.FindP(3,1)
    Q = wm.FindQ(3,1)
    
    PQM1 = np.zeros((c1.n+1,c1.n+1))
    PQM1[:,0:4] = P
    PQM1[:,4] = Q.T
    
    AB1 = np.linalg.inv(PQM1)
    
    #page 84:
    PHI2 = np.matmul(PQM1,AB1)
    
    #page 82, the big mystery:
    #Coarse = np.zeros()
    
    #page 83 synthesis:
    P2 = wm.FindP(3,2)
    Q2 = wm.FindQ(3,2)
    #TODO:... figure it out!
    
    """
    The B-spline filter bank
    is described on page 95
    """
    
    k=3
    nump=30
    qc =linear_vertices((0.,12.),(12.,0.),4)
    qc = Bspline(qc,k=3,nump=30)
    P2 = wm.FindP(2,1)
    Q2 = wm.FindQ(2,1)
    PQM_2 = np.zeros((qc.n,qc.n))
    PQM_2[:,0:3] = P2
    PQM_2[:,3] = Q2.T
    
    """
    Possible explanation:
        you need to start at level 2
        and work down to level 1
    """
    
    ii = interpolatedBspline(c1.vertices,k=4,nump=30)



    k=4
    nump=300
    vertices=np.asarray([[1.,0.,0.],
                         [1.,0.,3.], #[2.0,6.],
                         [1.,6.0,12.],
                         [1.,12.,12.]])
    v0 = Bspline(vertices, k, nump)
    plt.close()
    """
    v0.plotcurve_detailed()
    v0.t
    #"""
    v0.knot_insertion(.5)
    #v0.knot_removal(4)
    
    
    self = v0
    curve=self
    
    index=4
    u=.5
    r=4
    s=1
    num=1
    t=1
    #v0.plot3D()
    #plt.close()
    
    vertices = np.asarray([[  8.06747597,  13.35873976],
                           [  8.41295643,  13.66065505],
                           [ 10.78800443,  15.0847025 ],
                           [ 21.93061342,  15.0847025 ],
                           [ 49.50732787,  15.0847025 ],
                           [ 77.56538942,  15.0847025 ],
                           [ 93.64522185,  15.0847025 ],
                           [ 93.97358995,   8.0270766 ],
                           [ 95.71354028,   4.34564096],
                           [ 96.48102433,   3.30656573]])
    knots = np.asarray([ 0.        ,  0.        ,  
                        0.        ,  0.        ,  0.03933712,
                        0.21569948,  0.39206184,  0.56842421,  
                        0.74478657,  0.92114893,  1.        , 
                        1.        ,  1.        ,  1.        ])
                        
    CPK_cdr_b = Bspline(vertices,k=4,nump=30,t=knots)
    
    
    fullverts = np.asarray([[   0.        ,    0.        ],
                           [   3.55524179,    3.61325094],
                           [   6.91945396,    9.37268781],
                           [   6.89367322,   15.0847025 ],
                           [  21.93061342,   15.0847025 ],
                           [  49.50732787,   15.0847025 ],
                           [  77.56538942,   15.0847025 ],
                           [  93.64522185,   15.0847025 ],
                           [  94.04778224,    6.43245933],
                           [  97.26878801,    1.60429372],
                           [ 100.53889507,    0.        ]])
    fullknts = np.asarray([ 0.   ,  0.   ,  0.   ,  0.   ,  0.125,  0.25 ,  0.375,  0.5  ,
        0.625,  0.75 ,  0.875,  1.   ,  1.   ,  1.   ,  1.   ])
                        
    CProfile2D = Bspline(fullverts,k=4,nump=30,t=fullknts)
    
    
    iknts = interiorKnotsSet(CProfile2D.t)
    
    CPK_cdr_b.knot_insertion(.125)
    
    #homotopy knot transformation
    # knot_homotopy
    # rename to knot_deformation
    mknot = knot_deformation(CPK_cdr_b.t,
                          CProfile2D.t,10)
    
    """
        knot = self.t[index]
        #t,vertices = unsafe_knot_removal(self,knot,index)
        vertices = RemoveCurveKnot(self,
                                   u=knot,
                                   r=index,
                                   TOL=tol)
    """
    index = 4
    curve = copy.copy(CPK_cdr_b)
    knot = curve.t[index]
    u = knot
    r=index  #index
    s=1
    num=1
    TOL = 1.0000001
    
    
    
    
    vertices=np.asarray([   [-4.,0.,0.], [-2.,0.,0.],[0.,0.,0.],
                                [2.,0.,0.],
                                [3.0,0.,0.],
                                [6.0,5.,0.],
                                [9.0,12.,0.],
                                [11.,12.,0.],
                                [12.,12.,0.],[14.,12.,0.],[16.,12.,0.]])
                        
    a = Bspline(vertices,4,30)
    b = interpolatedBspline(vertices,4,30)
    c = copy.deepcopy(b)
    c.GlobalCurveInterp_gen1_THB_compatible_ver0(vertices)
    d = copy.deepcopy(b)
    d.GlobalCurveInterp_gen1_THB_compatible_ver1(vertices)
    e = copy.deepcopy(b)
    e.GlobalCurveInterp_gen1_THB_compatible_ver2(vertices)
    
    b.plot3d_wIvertices()
    """
    c.plot3d_wIvertices()
    d.plot3d_wIvertices()
    e.plot3d_wIvertices()
    
    print sum(np.linalg.norm(el) for el in c.smart_extremes())
    print sum(np.linalg.norm(el) for el in d.smart_extremes())
    print sum(np.linalg.norm(el) for el in e.smart_extremes())
    
    cksum = {0:sum(np.linalg.norm(el) for el in c.smart_extremes()),
             1:sum(np.linalg.norm(el) for el in d.smart_extremes()),
             2:sum(np.linalg.norm(el) for el in e.smart_extremes())}
    
    key = min(cksum, key=cksum.get)
    key = min(cksum, key=lambda k: cksum[k]) #same thing
    #"""