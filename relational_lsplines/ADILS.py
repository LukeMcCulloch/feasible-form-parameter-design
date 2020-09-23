"""
    Generate B-spline curve
    via
    form parmeter design 
    with Interval Analysis
    
    TODO: 
        Dimension Agnostic - 2D only at the moment
        
    Guiding facts:
        
        
"""
import numpy as np
import copy
import matplotlib.pyplot  as plt
import matplotlib.patches as patches
import matplotlib.path    as path


from mpl_toolkits.axes_grid.axislines import SubplotZero
from matplotlib.transforms import BlendedGenericTransform
from functools import partial

import  curve         as     spline
from  automatic_differentiation import ad
from  automatic_differentiation import IntervalAnalysis as IA
from  interval_arithmetic import ia #FUN todo: try 
from  extended_interval_arithmetic import ia as eia
# making an extended arithmetic solver! 
# right now this must use non-extended interval arithmetic because
# the inversion of the interval Hessian is 
# only coded for non-singular interval matrices
#
#from    fuzzylogic        import very, extremely, somewhat, slightly, up, down, tri, utri, trap, not_this
#from    interval_analysis import IntervalAnalysis as IA
from    FormParameter     import FormParameterDict, generalized_aattractor
from    initialValues     import InitializeControlPoints, InitializeControlVertices, \
                                  interval_bounds, lagrangian_bounds
from    plots             import Plotter

#np.set_printoptions(precision=2)
#*****************************************
# THB spline methods:
from myspline import rbspline, BoundsList



###############################################################################
class LazyFPDConstraintInterface(object):
    """A simplified Lazy Constraint Construction Class
    
    Goals of this class:
        *I'd like to make the lagrange curve solver easier to use
            (especially programatically :- by automation)
            for the more esoteric contraints or 
            delay setup for any constraints
        *make restart easier for knot vector deformations
        
    This is a bit of functional lazy code to head in that direction
        
        *use case 1: development of the split curve CPKeel
        in the bare hull maker, I want to specify some FPD constraints
        but I don't really want to pick the curve I want to optimize
        just yet.  It's too messy to read that way.
    """
    def __init__(self):
        self.lazyconstraints = []   #this will hold the rules 
                                    #ready for when FPD (curve) is availible
                                    #(FPD memoizes basis data by the way)
        self.lazyspecial = [] #Lagrangian.sp = {this stuff}
    
    def Constraint_Area(self, value, location=None, 
                           kind = 'equality',weight=1.0):
        if location is not None:
            print 'note: location is not used with the area constraint'
        
        def set_constraint(FPD):
            """transform the form parameter class
            and return it
            """
            FPD.add_AreaConstraint(value = value,
                                    kind=kind,
                                    weight=weight)
            return FPD
        
        self.lazyconstraints.append( set_constraint )
        
        return #set_constraint   # do not return 
                                 # handle it internally for ease of programming 
                                 # later
                                 
                                 
    
    def Constraint_Xc(self, value, location=None, 
                           kind = 'equality',weight=1.0):
        if location is not None:
            print 'note: location is not used with the area constraint'
        
        def set_constraint(FPD):
            """transform the form parameter class
            and return it
            """
            FPD.add_XcConstraint(value = value,
                                 kind=kind,
                                 weight=weight)
            return FPD
        
        self.lazyconstraints.append( set_constraint )
        return #set_constraint
    
    
    
    
    def Constraint_Tangent(self, value, location, 
                           kind = 'equality',weight=1.0):
        def set_constraint(FPD):
            """transform the form parameter class
            and return it
            """
            FPD.add_AngleConstraint(value = value,
                                    location=location,
                                    kind=kind,
                                    weight=weight)
            return FPD
        
        
        self.lazyconstraints.append( set_constraint )
        return 
    
    
    
    
    def Constraint_Curvature(self, value, location, 
                           kind = 'equality',weight=1.0):
        def set_constraint(FPD):
            """transform the form parameter class
            and return it
            """
            FPD.add_CurvatureConstraint(value = value,
                                        location=location,
                                        kind=kind,
                                        weight=weight)
            return FPD
        
        
        self.lazyconstraints.append( set_constraint )
        return 
    
    
    def Constraint_ArcLength(self, value, 
                           kind = 'equality',weight=1.0):
        def set_constraint(FPD):
            """transform the form parameter class
            and return it
            """
            FPD.add_ArcLength(value = value,
                              kind=kind,
                              weight=weight)
            return FPD
        
        
        self.lazyconstraints.append( set_constraint )
        return 
    def Constraint_ArcLengthApprox(self, value, 
                           kind = 'equality',weight=1.0):
        def set_constraint(FPD):
            """transform the form parameter class
            and return it
            """
            FPD.add_ArcLengthApprox(value = value,
                                    kind=kind,
                                    weight=weight)
            return FPD
        
        
        self.lazyconstraints.append( set_constraint )
        return 
    
    
    

    def Constraint_Vertex(self, index_tuple, value, location, 
                           kind = 'equality', weight=1.0, doall=False):
        """index_tuple : (index_of_vertex, dimension_to_constrain)
        """
        def set_constraint(FPD):
            """transform the form parameter class
            and return it
            """
            ind = index_tuple[0]
            dim = index_tuple[-1]
            if dim==0 or doall:
                FPD.add_xVertexConstraint(value = value,
                                    location=location,
                                    kind=kind,
                                    weight=weight,
                                    index = ind)
            if dim==1 or doall:
                FPD.add_yVertexConstraint(value = value,
                                    location=location,
                                    kind=kind,
                                    weight=weight,
                                    index = ind)
            if dim==2 or doall:
                FPD.add_zVertexConstraint(value = value,
                                    location=location,
                                    kind=kind,
                                    weight=weight,
                                    index = ind)
            return FPD
        
        self.lazyconstraints.append( set_constraint )
        return 
    
    
    def Constraint_Fix(self, index_tuple, value, location=None, 
                           kind = 'equality', weight=1.0, doall=False):
        """index_tuple : (index_of_vertex, dimension_to_constrain)
        """
        def set_constraint(FPD):
            """transform the form parameter class
            and return it
            """
            ind = index_tuple[0]
            dim = index_tuple[-1]
            if dim==0 or doall:
                FPD.add_xFixity(index=ind,
                                value=value)
            if dim==1 or doall:
                FPD.add_yFixity(index=ind,
                                value=value)
            if dim==2 or doall:
                FPD.add_zFixity(index=ind,
                                value=value)
            return FPD
        
        self.lazyconstraints.append( set_constraint )
        return 
    
    
    def Constraint_Special(self, gfunc_closure):
        """index_tuple : (index_of_vertex, dimension_to_constrain)
        """
        def set_constraint(FPD):
            """transform the form parameter class
            and return it
            """
            return gfunc_closure(FPD.curve)
        
        self.lazyspecial.append( set_constraint )
        return 
    
    


###############################################################################
class CurveDesigner(object):
    """Simplified Interface to Lspline Curve Solver
    
    NOT COMPLETE, NOT USED
    
    Use Lagrangian directly to add constraints
    which do not conform to these inputs.
    
    Parameters and constraints 
    the B-spline curve must satisfy
    ----------
        xb           : 1st  control point, index 0
        yb           : 1st  control point, index 1
        xe           : last control point, index 0
        ye           : last control point, index 1
        alphab       : tangent at the curve start
        alphae       : tangent at the curve end
        Cab          : curvature at the curve start
        Cae          : curvature at the curve start
        area         : area between the curve and the x-axis (index 0)
        Xc           : 'x' location of the center of area of the curve
        k            : curve Order (degree+1)
        nump         : number of evaluation stations when drawing or otherwise 
                       approximating the curve 
        num_vertices : number of control vertices for the curve
        
    Returns
    ----------
    Notes
    ----------
        *Note finished!
        *Possibly not to be finished until after the PhD?
    """
    def __init__(self,  curve=None,
                        xb=None,
                        xe=None,
                        yb=None,
                        ye=None,
                        alphab=None,
                        alphae=None,
                        Cab=None,
                        Cae=None,
                        area = None,
                        Xc = None,
                        Yc = None,
                        k=4,
                        nump=30,
                        num_vertices=None):
        if curve is None:
            assert(xb is not None),'ERROR: curve or endpoints required, xb missing'
            assert(yb is not None),'ERROR: curve or endpoints required, yb missing'
            assert(xe is not None),'ERROR: curve or endpoints required, xe missing'
            assert(yb is not None),'ERROR: curve or endpoints required, ye missing'
            self.xb     = xb
            self.yb     = yb
            self.xe     = xe
            self.ye     = ye
            vertices    = spline.linear_vertices([xb,yb],
                                                 [xe,ye],
                                                 num=num_vertices)
            self.curve = spline.Bspline(vertices,k,nump)
        else:
            assert(xb is None),'ERROR: curve given, xb overspecified'
            assert(yb is None),'ERROR: curve given, yb overspecified'
            assert(xe is None),'ERROR: curve given, xe overspecified'
            assert(yb is None),'ERROR: curve given, ye overspecified'
            self.xb     = curve.vertices[0,0]
            self.yb     = curve.vertices[0,1]
            self.xe     = curve.vertices[-1,0]
            self.ye     = curve.vertices[-1,1]
            self.curve  = curve
            k = curve.k
            nump = curve.nump
        
        
        self.__version__= '$Revision: 0.1'.split()[1]
        self.__usage__  = 'usage: CurveDesigner(curve data...)'
        self.verbose    = False
        self.k          = k
        self.nump       = nump
        self.alphab     = alphab
        self.alphae     = alphae
        self.Cab        = Cab
        self.Cae        = Cae
        self.area       = area
        self.Xc         = Xc
        self.Yc         = Yc
        
    def __call__(self):
        return
        


###############################################################################
class FPDIdeal(object):
    """FPD data, rule, and curve class
        input: 
            -form parameters for a curve
            -rules defining its validity
    NOT COMPLETE, a version of this is really used in 
    GeometryGenerator.py
    """
    def __init__(self, xb=0.,yb=12.,xe=12.,ye=0.,
                 alphab=None,alphae=None,
                 Cab_given=None,Cae_given=None,
                 area=None,xc=None,yc=None,
                 ymax=None,
                 nCV = 7, 
                 stationloc=None, base=None,
                 nump=30, k=4,
                 tol = 1.e-7):
        self.xb         = xb
        self.yb         = yb
        self.xe         = xe
        self.ye         = ye
        self.alphab     = alphab
        self.alphae     = alphae
        self.Cab_given  = Cab_given
        self.Cae_given  = Cae_given
        self.area       = area
        self.ymax       = ymax
        self.xc         = xc
        self.yc         = yc
        self.nCV        = nCV
        self.nump       = nump
        self.stationloc = stationloc,
        self.base       = base
        self.k          = k
        self.tol        = tol
        self.base       = base
        
        if self.xc  is 'default':
            self.xc = self.xe/2.
        
    def make_curve(self):
        self.vertices = InitializeControlVertices(self.xb,
                                               self.yb,
                                               self.xe,
                                               self.ye,
                                               self.alphab,
                                               self.alphae,
                                               self.Cab_given,
                                               self.Cae_given,
                                               self.area,
                                               self.xc,
                                               self.yc,
                                               self.nCV)
        
        self.curve = spline.Bspline(self.vertices.vertices,
                                    self.k,self.nump)
        return #self.curve
    
    def solve_curve_ini(self):
        """
            STAGE 1
        """
        interval_data, small = interval_bounds(self.curve)
        FPD = self.setup_basic()
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                                    interval_data, small, 1.e4)
        self.Lspline = IntervalLagrangeSpline(self.curve,
                                         L, data = interval_data)
        
        self.Lspline.curve.pts_M_pts()
        self.Lspline.curve.compute_arclength()
        self.Lspline.optimize(stop=30)
        return
    
    def solve_curve_fini(self,Xc=None):
        """
            STAGE 2
        """
        interval_data, small = interval_bounds(self.Lspline.curve)
        
        if Xc is None:
            FPD = self.setup_stage2()
        else:
            FPD = self.setup_stage2(Xc)
        
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, 
                                    interval_data, small, 1.e4)
        self.Lspline = IntervalLagrangeSpline(self.curve,
                                         L, data = interval_data)
        
        self.Lspline.curve.pts_M_pts()
        self.Lspline.curve.compute_arclength()
        
        self.Lspline.optimize(stop=30)
        return
    
        
    def setup_basic(self):
        FPD = FormParameterDict(self.curve)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', 
                                weight = .5)
        FPD.add_AreaConstraint(kind='equality', 
                               value = self.area)
        
        return FPD
        
    def setup_stage2(self, Xc = None):
        interval_data, small = interval_bounds(
                self.Lspline.curve)
        
        E1_0 = copy.deepcopy(self.Lspline.curve.E1)
        E2_0 = copy.deepcopy(self.Lspline.curve.E2)
        E3_0 = copy.deepcopy(self.Lspline.curve.E3)
        S_0  = copy.deepcopy(self.Lspline.curve.AL)
        
        FPD = FormParameterDict(self.curve)
        if Xc is None and self.xc is None:
            pass
        else:
            if Xc is None: Xc = self.xc
            FPD.add_XcConstraint(kind = 'equality', 
                                 value = Xc)
        FPD.add_E1(kind='LS', weight = .5/E1_0.value)
        FPD.add_E2(kind='LS', weight = .5/E2_0.value)
        FPD.add_E3(kind='LS', weight = .5/E3_0.value)
        FPD.add_ArcLengthApprox(kind='LS', 
                                weight = .5/S_0)
        FPD.add_AreaConstraint(kind='equality', 
                               value = self.area)
        
        """
        TODO: match form curves
        """
        FPD.add_yPointConstraint(kind= 'max', 
                                 value=self.ymax, 
                                 location=0.5 )
        
        if (np.sqrt((self.ymax-self.yb)**2)<self.tol):
            
            dx = (self.xe-self.xb)/2.
            
            FPD.add_yVertexConstraint(kind = 'equality', 
                                      value=self.ymax, 
                                      index =1)
            FPD.add_yVertexConstraint(kind = 'equality', 
                                      value=self.ymax, 
                                      index =2)
            FPD.add_xVertexConstraint(kind = 'equality', 
                                      value=dx, 
                                      index =2)
        
        
            
        
        return FPD
    
    def check_area_to_box_relation(self):
        dx = self.xe-self.xb
        dy1 = self.ye - self.yb
        dy2 = self.ymax - self.yb
        dy = max(dy1,dy2)
        Amax = dx*dy
        assert Amax>self.area, \
             "Area {} > {}feasibility fail".format(Amax,
                                            self.area)
        return  
    


###############################################################################
#class dummy(object):
#    def __init__(self):
#        return
###############################################################################
#method_dict = {'equality':self.compute,
#                    'max':self.compute_max,
#                    'min':self.compute_min,
#                    'LS':self.compute_LS}
###############################################################################
class SolverPostProcessor(object):
    """Post Process an Lspline.Lagrangian
    """
    def __init__(self,Lspline):
        self.Lspline = Lspline
        self.method_dict = {'equality':Lspline.Lagrangian.equality,
                           'LS':Lspline.Lagrangian.obj,
                           'max':Lspline.Lagrangian.maximum,
                           'min':Lspline.Lagrangian.minimum}
        self.badparameters = {}
        self.tol = Lspline.tol
        self.parse_results()
        
    def find_normalization(self):
        L = self.Lspline.Lagrangian
        val = 0.
        for key in L.obj:
            fp = L.obj[key]
            ttype, desired, actual, weight = self.getdata(fp)
            #print ttype, desired, actual, weight
            diff = self.difference(desired,actual)
            val = max(diff,val)
        return val
        
    def parse_results(self):
        self.badparameters['equality'] = self.parsetype('equality')
        self.badparameters['max'] = self.parsetype('max')
        self.badparameters['min'] = self.parsetype('min')
        
        self.normalization_constant = self.find_normalization
        return
    
    def parsetype(self, typenm):
        """
        Notes
        ----------
            the computed value will represent different things
            depending on the constraint.
            For, instance, in the case of VertexConstraints, 
            computed_value it is the difference
            between desired and actual value.
            
            -this is presently a bug in the 
            checking.
            -workaround is to include a 
            convergence check first.
            -if that is ok, then things 
            are ok!
            -if not, then we need a local check 
            of the gradient to be 
            `most general` in checking for 
            bad constraints.
            -but even this is complicated by
            slack variables?* (havent gone back to
            think about this bit)
        """
        sub_dict_to_analyze = self.method_dict[typenm]
        results = {}
        for key in sub_dict_to_analyze:
            
            fp = sub_dict_to_analyze[key]
            
            ttype, desired, actual, weight = self.getdata(fp)
            
            diff = self.difference(desired, actual)
            
            ck_failtrue = self.ck(diff)
            
            #see note above for why this became so baroque
            # break out by constraint type, or just 
            # check gradient first in the future.
            li = sub_dict_to_analyze[key].index
            if self.Lspline.conv > self.Lspline.tol:
                if ck_failtrue and self.Lspline.mask[li]:
                    if np.linalg.norm(self.Lspline.f.grad[0,li])>0.:
                        #print 'detailing li = ',li
                        results[key] = {'type':ttype,
                                        'value':actual,
                                        'desired_value':'pass_value',
                                        'pass_value':desired,
                                        'weight':weight}
        return results
    
    def getdata(self, fp):
        """Getting the constarint data:
        
        returns:
            -------
            -type: geometric type of constaint (area, tangent, etc)
            -pass_value: Given during setup, this is the form parameter
            -computed_value: value updates when the curve updates
            -weight: not used unless this is an objective (or penalty function)
        """
        return fp.type, fp.pass_value, fp.computed_value, fp.weight
    
    def difference(self, expected, actual):
        try:
            return np.linalg.norm(actual-expected)
        except:
            return np.linalg.norm(actual.value-expected)
            
    
    def ck(self, diff):
        if diff>self.tol:
            return True
        else:
            return False
        
        
###############################################################################
class FritzJohn(object):
    def __init__(self):
        return
###############################################################################
class Lagrangian(object):
    def __init__(self, FPD, fritzjohn_conditions=False,tol=1./(10.**7),
                 parent=None,child=None):
        self.curve      = FPD.curve
        self.FPD        = FPD #now store this to facilitate multigrid contraints using a native FPD prolongation method
        self.obj        = {}
        self.equality   = {}
        self.minimum    = {}
        self.maximum    = {}
        self.fix        = {}    #fixity, 1 DOF per 'constraint'
        self.sp         = {}    #implement function constraints, forces, cusps, 
                                #(I thought this had cusps... where did that go?)
                                #are they special though?
        self.method_dict = {'equality':self.equality,
                            'max':self.maximum,
                            'min':self.minimum,
                            'fix':self.fix,
                            'LS':self.obj}
        #****************************************
        # simple hierarchical
        # limitations with projecting constraints...
        self.parent = parent
        self.child = child
        #****************************************
        self.fjc = fritzjohn_conditions
        if self.fjc:
            self.fj  = FritzJohn()
        self.parse_form_parameters(FPD)
        self.parse_fixed_parameters(FPD)
        self.N = self.compute_dim_space()    
        self.tol = tol
        
    
    def reset(self,
              save_obj=True,
              save_equality=True,
              save_min=True,
              save_max=True,
              save_fix=True,
              save_sp=True):
        if not save_obj:
            self.obj = {}
        if not save_equality:
            self.equality = {}
        if not save_min:
            self.minimum = {}
        if not save_max:
            self.maximum = {}
        if not save_fix:
            self.fix = {}
        if not save_sp:
            self.sp = {}
        
        self.method_dict = {'equality':self.equality,
                            'max':self.maximum,
                            'min':self.minimum,
                            'LS':self.obj}
        return
    
    def switch_all_of_one_Lagrangian_kind(self,ini_kind='equality',
                                          fini_kind = 'LS'):
        for key in self.method_dict[ini_kind]:
            self.modifyLagrangian(key = key,
                                  kind = ini_kind,
                                  newkind = fini_kind
                                  )
        self.method_dict[ini_kind] = {}
        setattr(self, ini_kind, {})
        self.N = self.compute_dim_space()  
        return
    
    def restore_Lagrangian(self, FPD = None):
        if FPD is None:
            FPD = self.FPD
        
        self.obj        = {}
        self.equality   = {}
        self.minimum    = {}
        self.maximum    = {}
        self.method_dict = {'equality':self.equality,
                            'max':self.maximum,
                            'min':self.minimum,
                            'fix':self.fix,
                            'LS':self.obj}
        
        self.parse_form_parameters(FPD)
        self.parse_fixed_parameters(FPD)
        self.N = self.compute_dim_space() 
        return
    
    def modifyLagrangian(self, key, kind, newkind=None, 
                         new_value=None, new_weight=None):
        """Hooks to modify the Lagrangian
        
        notes
        ----------
        use this to switch from 
        LS interpolation to equality when solving for 
        longitudinals.
        -nice way to initiate a solution when the 
        system is overly constrained at the coarse level
        
        """
        
        if newkind is not None:
            self.method_dict[newkind][key] = self.method_dict[kind][key]
            #del(self.method_dict[kind][key])
            kind = newkind
            
        if new_value is not None:
            self.method_dict[kind][key].pass_value = new_value
            self.method_dict[kind][key].given = new_value
        
        if new_weight is not None:
            self.method_dict[kind][key].weight = new_weight
        return
    
    
    
    def new_prolongation_level(self, refined_curve, FPD=None):
        """Once a THB spline has been refined (prolongated)
        use this to refine the Lagrangian while keeping constraints
        intact
        
        Adds a level to the Lagrangian Class
        
        notes
        ---------
        formerly called prolongate
        
        This adds another fine level to the Lagrangian level system
        it does not prolongate the solution!
        
        may need to rename to 'add_prolongation_level'
        or similar.
        (same for Lagrangian.prolongation and FPD.prolongation)
        """
        if FPD is None:
            FPD = self.FPD.new_prolongation_level(refined_curve)
        
        self.child =  Lagrangian(FPD,parent=self)
        return self.child
    
        
    def parse_form_parameters(self, FPD):
        for el in FPD.FormParam:
            self.method_dict[FPD.FormParam[el].kind][el] = FPD.FormParam[el]
        return
        
    
    def parse_fixed_parameters(self, FPD):
        for el in FPD.FixParam:
            self.method_dict[FPD.FixParam[el].kind][el] = FPD.FixParam[el]
        return
    
    
    def print_objective_functions(self):
        """Post Processing
        """
        for el in self.obj:
            name = self.obj[el].type
            value = self.obj[el].computed_value.value
            dvalue = self.obj[el].pass_value
            weight = self.obj[el].weight
            print '\n {} objective => {}'.format(el, name)
            print 'value = {}'.format(value)
            print 'desired value = {}'.format(dvalue)
            print 'weight = {}'.format(weight)
        return
            
    def print_form_parameters(self, kind='equality'):
        """Post Processing
        """
        if kind is 'all':
            self.print_form_parameters(kind='eq')
            self.print_form_parameters(kind='ls')
            self.print_form_parameters(kind='maxmin')
            return
        elif 'eq' in kind:
            print '--------------------------------------------------'
            print '\nEquality Constraints:'
            for el in self.equality:
                ttype = self.equality[el].kind
                expected = self.equality[el].pass_value
                actual = self.equality[el].computed_value.value
                weight = self.equality[el].weight
                print '\n {} constraint => {}'.format(el,self.equality[el].type)
                if 'angle' in ttype:
                    print 'expected = {}, found = {}'.format(np.rad2deg(expected), actual)
                    diff = np.linalg.norm(np.rad2deg(expected) -  actual)
                else:
                    print 'expected = {}, found = {}'.format(expected, actual)   
                    diff = np.linalg.norm(expected -  actual) 
                print 'difference = {}'.format(diff)
                print 'weight = {}'.format(weight)
                
        elif 'ls' in kind:
            print '--------------------------------------------------'
            print '\nLeast Square Objectives:'
            for el in self.obj:
                desired = self.obj[el].pass_value
                actual = self.obj[el].computed_value.value
                weight = self.obj[el].weight
                print '\n {} objective => {}'.format(el,self.obj[el].type)
                print 'desired = {}, found = {}'.format(desired, actual) 
                print 'weight = {}'.format(weight)
        else:
            print '--------------------------------------------------'
            print '\nMax Constraints:'
            for el in self.maximum:
                expected = self.maximum[el].pass_value
                actual = self.maximum[el].computed_value.value
                weight = self.maximum[el].weight
                print '\n {} constraint => {}'.format(el,self.maximum[el].type)
                print 'expected = {}, found = {}'.format(expected, actual) 
                print 'weight = {}'.format(weight)
            print '--------------------------------------------------'
            print '\nMin Constraints:'
            for el in self.minimum:
                expected = self.minimum[el].pass_value
                actual = self.minimum[el].computed_value.value
                weight = self.minimum[el].weight
                print '\n {} constraint => {}'.format(el,self.minimum[el].type)
                print 'expected = {}, found = {}'.format(expected, actual) 
                print 'weight = {}'.format(weight)
        return
    
    
    
            
    def get_anomolous_form_parameters(self):
        """Post Processing
        """
        aeq = []
        amax = []
        amin = []
        print '\nChecking Equality Constraints'
        for el in self.equality:
            ttype = self.equality[el].kind
            expected = self.equality[el].pass_value
            actual = self.equality[el].computed_value.value
            weight = self.equality[el].weight
            if 'angle' in ttype:
                diff = np.linalg.norm(np.rad2deg(expected) -  actual)
            else: 
                diff = np.linalg.norm(expected -  actual) 
            if diff>self.tol:
                aeq.append(el)
        print '\nChecking Max Constraints:'
        for el in self.maximum:
            expected = self.maximum[el].pass_value
            actual = self.maximum[el].computed_value.value
            weight = self.maximum[el].weight
            if expected<actual:
                diff = np.linalg.norm(expected -  actual) 
                if diff>self.tol:
                    amax.append(el)
        print '\nChecking Min Constraints:'
        for el in self.minimum:
            expected = self.minimum[el].pass_value
            actual = self.minimum[el].computed_value.value
            weight = self.minimum[el].weight
            print '\n {} constraint => {}'.format(el,self.minimum[el].type)
            print 'expected = {}, found = {}'.format(expected, actual) 
            print 'weight = {}'.format(weight)
            if expected>actual:
                diff = np.linalg.norm(expected -  actual) 
                if diff>self.tol:
                    amin.append(el)
        return aeq,amax,amin
                
    
    def compute_dim_space(self):
        """Compute the dimensions
        of the Lagrangian optimization space
        """
        dim     = self.curve.dim
        n       = dim*self.curve.n
        ilen    = len(self.minimum) + len(self.maximum)
        self.ilen = ilen
        nfp     = len(self.equality) + 2*ilen
        Ndim       = n + nfp
        if self.fjc:
            Ndim += 1
        return Ndim
#    
#    def prolong_space(self, curve):
#        """
#        4,5,7,11,19,35,...
#        """
#        dim     = curve.dim
#        n       = dim*curve.n
#        
#        ilen    = len(minimum) + len(maximum)
#        ilen = ilen
#        nfp     = len(equality) + 2*ilen
#        Ndim       = n + nfp
#        if fjc:
#            Ndim += 1
#        return Ndim
    
    
    def setup_derivatives(self):
        """Instantiates the independent variables of the optimization as
        scalarAD points
        """
        dim  = self.curve.dim
        Ndim = self.N
        ilen = self.ilen
        
        count = 0
        pts = []
        for j in range(dim):
            p = []
            for i, pt in enumerate(self.curve.vertices):
                value = self.curve.vertices[i,j]
                p.append( ad(value, dim=count, N=Ndim, of_scalars = True) )
                count += 1 #could count with  count = i+self.curve.n*j instead - don't get lazy!
            pts.append(p)
            
        for el in self.equality:
            self.equality[el].index = count
            self.equality[el].Lagrange = ad(0., dim=count, N=Ndim, of_scalars = True)
            count +=1
            
        for el in self.minimum:
            self.minimum[el].index       = count
            self.minimum[el].slack_index = count+ilen
            self.minimum[el].Lagrange    = ad(1.0, dim=count, N=Ndim, of_scalars = True )
            self.minimum[el].slack       = ad(1.0, dim=count+ilen, N=Ndim, of_scalars=True )
            count +=1
            
        for el in self.maximum:
            self.maximum[el].index = count
            self.maximum[el].slack_index = count+ilen
            self.maximum[el].Lagrange = ad(1.0, dim=count, N=Ndim, of_scalars = True )
            self.maximum[el].slack = ad(1.0, dim=count+ilen, N=Ndim, of_scalars = True )
            count +=1
        
        if self.fjc:
            self.fj.weight = 1.0
            self.fj.Lagrange = ad(1.0, dim=count, N=Ndim, of_scalars = True )
            count += 1
            
        ptslist = []
        for td in range(dim):
            ptslist.append(pts[td])
        return ptslist#np.asarray(pts[0]), np.asarray(pts[1])
        
    def setup_simple_interval(self, imin = 5.0, imax = 5.0 ):
        curve = self.curve
        interval_x = []
        interval_y = []
    
        count = 0
        for i,vertex in enumerate(curve.xpts):
            a = copy.deepcopy(vertex)
            c = copy.deepcopy(vertex)
            a.value = a.value - imin    
            c.value = c.value + imax
            if a.value < 0.:
                a.value = 0.
            interval_x.append( ad( ia(a.value,c.value), N=self.N, dim = count, name = 'x{}'.format(i) ) )
            count += 1
        
        for i,vertex in enumerate(curve.ypts):
            a = copy.deepcopy(vertex)
            c = copy.deepcopy(vertex)
            a.value = a.value - imin
            c.value = c.value + imax
            if a.value < 0.:
                a.value = 0.
            interval_y.append( ad( ia(a.value,c.value), N=self.N, dim = count, name = 'y{}'.format(i)  ) )
            count += 1
        
             
        for el in self.equality:
            a = copy.deepcopy(self.equality[el].Lagrange)
            c = copy.deepcopy(self.equality[el].Lagrange)
            a.value = a.value - imin
            c.value = c.value + imax
            if a.value < 0.:
                a.value = 0.
            self.equality[el].interval_Lagrange = ad( ia(a.value, c.value), N=self.N, dim = count, name = 'lambda{}'.format(el)  )
            count += 1
            
        #TODO:  add interval inequalities
    
        #        ilen = len(IneqConstraints)
        #        for key in IneqConstraints:
        #            IneqConstraints[key]['interval_Lagrange'] = scalarAD(1.0,np.matrix(Imatrix[count]),np.matrix(np.zeros((N,N),float)))
        #            IneqConstraints[key]['interval_slack']    = scalarAD(1.0,np.matrix(Imatrix[count+ilen]),np.matrix(np.zeros((N,N),float))) 
        #            count +=1
        
        if self.fjc:
            self.fj.weight = 1.0
            self.fj.interval_Lagrange = ad( ia(0.,1.), N=self.N, dim=count, name='lambda_fj')
            count += 1
            
        ixpts = np.asarray(interval_x)
        iypts = np.asarray(interval_y)
    
        return ixpts, iypts
        
    def interval_Lagrange_ranged(self, curve,
                     ixmin, ixmax, 
                     iymin, iymax, 
                     ilmin=None, ilmax=None,
                     fjcmin = None, fjcmax = None):

        
        interval_x = []
        interval_y = []
        
        count = 0
        for i,vertex in enumerate(curve.xpts):
            a = copy.deepcopy(vertex)
            c = copy.deepcopy(vertex)
            a.value = ixmin[i]
            c.value = ixmax[i]
            interval_x.append(ad(ia(a.value, c.value), N=self.N, dim=count ) )
            count +=1 
        
        for i,vertex in enumerate(curve.ypts):
            a = copy.deepcopy(vertex)
            c = copy.deepcopy(vertex)
            a.value = iymin[i]
            c.value = iymax[i]
            interval_y.append(ad(ia(a.value, c.value), N=self.N, dim=count ) )
            count +=1 
        
        if ilmin is not None:
            for i, el in enumerate(self.equality):
                a = copy.deepcopy(self.equality[el].Lagrange)
                c = copy.deepcopy(self.equality[el].Lagrange)
                a.value = ilmin[i]
                c.value = ilmax[i]
                self.equality[el].interval_Lagrange = ad(ia(a.value, c.value), 
                                                             N=self.N, dim=count )
                count +=1 
        
        if self.fjc:
            print 'self.fjc = {}'.format(self.fjc)
            self.fj.weight = 1.0
            self.fj.interval_Lagrange = ad( ia(fjcmin,fjcmax), 
                                           N=self.N, dim=count, name='lambda_fj')
            count += 1
        ixpts = np.asarray(interval_x)
        iypts = np.asarray(interval_y)
    
        return ixpts, iypts
                
    def equality_compare(self):
        for el in self.Lagrangian.equality:
            print '-----',el,'-----'
            print self.equality[el].Lagrange
            print self.equality[el].interval_Lagrange
            print 'contains? = ',\
                    self.equality[el].interval_Lagrange.contains(
                            self.equality[el].Lagrange)
        return    
#   
###############################################################################
#
#   Consider writing a meta-method to combine 
#   FPD, Lagrangian, and IntervalLagrangeSpline
#   into one caller
#
# then just prolongate once to refine
#
#
###############################################################################
#
class IntervalLagrangeSpline(object):
    #def __init__(self, curve, Lagrangian, Objhash, Weight=None, points=None, s=None, data={}):
    def __init__(self, curve, Lagrangian, fritz_john_conditions = False,
                 Weight=None, points=None, s=None, data={},
                 parent=None,child=None,
                 THBactiveset=None,
                 newfixity=None):
        self.__version__    = '$Revision: 0.1'.split()[1]
        self.__usage__      = 'usage: IntervalLagrangeSpline(curve, Lagrangian, Weight, points, s, data)'
        self.verbose        = False
        self.curve          = curve
        self.Lagrangian     = Lagrangian
        
        condition='kkt'
        if fritz_john_conditions: condition='fj'
        self.conditions     = {'kkt':0,'fj':1}
        
        self.condition      = self.conditions[condition]
        self.f              = None
        self.fjc            = fritz_john_conditions
        #self.Objhash        = Objhash
        self.nitr           = 0
        self.Weight         = Weight 
        self.points         = points
        self.s              = s
        self.data           = data
        self.sor            = 1.0   # Successive Relaxation Coefficient - not needed for conergence
        self.conv = None
        self.setup()
        self.create_mask() 
        self.THBactiveset   = THBactiveset
        self.set_fixity()#THBactiveset)
        self.ia = IA(self.mask) #to be singleton? - pass mask to routines...
        self.instate_AD_vertices() #self.curve.xpts, self.curve.ypts
        self.check_IA_correctness = False
        self.type_preference = {'area':0,'angle':0,'curvature':0,'Xc':1,'Yc':1}
        self.tpl = ['area','angle','curvature','moment']
        #
        self.restore = {}
        #
        #**************************************************
        # hierarchy
        self.parent = parent
        self.child = child
        #
        #**************************************************
        #
        if 'xmin' in self.data:
            self.instate_interval_vertices()
        else:
            self.instate_interval_vertices(i_min = 0.05, i_max = 0.05)
        return
        
    def __call__(self, vertices_, Lagrangian_, stop_ = 15):
        self.optimize(vertices=vertices_, 
                      stop = stop_, 
                      Lagrangian = Lagrangian_)
        return
    
    def setup(self, conditions='fj'):
        self.Inumber=15
        self.tol=1./(10.**7)
        self.map_variables()
        self.IniVertices    = copy.copy(self.curve.vertices)
        self.nfv            = self.curve.n*self.curve.dim 
        self.nec            = len(self.Lagrangian.equality)
        self.nic            = len(self.Lagrangian.minimum) + len(self.Lagrangian.maximum)
        if self.condition == self.conditions['kkt']:
            self.n              = self.nfv + self.nec + 2*self.nic
        elif self.condition == self.conditions['fj']:
            self.n              = self.nfv + self.nec + 2*self.nic + 1
        #*******************************************
        # THB optimization of curves
        #self.maxgrads = np.zeros((self.n,1),float)
        self.maxgrads = np.zeros((self.curve.n*self.curve.dim,1),float)
        self.initialgrads = np.zeros_like(self.maxgrads)
        # store the maximum gradient for each component
        self.refine_by_these_grad_indices = None
        return
    
    def switch_all_of_one_Lagrangian_kind(self,ini_kind='equality',
                                          fini_kind = 'LS'):
        #
        self.Lagrangian.switch_all_of_one_Lagrangian_kind(
                                            ini_kind,fini_kind)
        self.setup()
        self.create_mask() 
        self.set_fixity()#self.THBactiveset)
        self.ia = IA(self.mask) 
        self.instate_AD_vertices() 
        return
    
    def restore_Lagrangian(self):
        """
        dev
        ----------
        from  automatic_differentiation import IntervalAnalysis as IA
        """
        
        refine_by_these_grad_indices = self.refine_by_these_grad_indices
        self.Lagrangian.restore_Lagrangian()
        self.setup()
        self.create_mask() 
        self.set_fixity()#self.THBactiveset)
        self.ia = IA(self.mask) 
        self.instate_AD_vertices() 
        self.refine_by_these_grad_indices = refine_by_these_grad_indices
        return
    
    def normalize(self, norms=None, scale=.5):
        from utility_optimization import package_norms
        L = self.Lagrangian
        if norms is None:
            norms = package_norms(L, recalc=False)
        E1_0    = norms[0]
        E2_0    = norms[1]
        E3_0    = norms[2]
        S_0     = norms[3]
        for el in L.obj:
            fp = L.obj[el]
            if fp.type == 'E1':
                fp.weight = scale/E1_0
            elif fp.type == 'E2':
                fp.weight = scale/E2_0
            elif fp.type == 'E3':
                fp.weight = scale/E3_0
            elif fp.type == 'AL':
                fp.weight = scale/S_0
        return
    
    def refine(self):
        """Use self.new_prolongation_level to build up 
        a Lspline solver for a new dyadily refined curve
        
        notes
        ----------
        it DOES matter which level you call refine with
        
        because self.refine_by_these_grad_indices
        knows only about solving it has done on its level
        
        
        
        
        """
        print 'refine called...'
        #if self.refine_by_these_grad_indices is None:
        #   refine_by_these_grad_indices = []
        #
        #*********************************************************************
        # get bounds
        boundslist = []
        if self.refine_by_these_grad_indices:
            for i in self.refine_by_these_grad_indices:
                cvi = self.get_vertex_index_by_any_solver_component_index(i)
            
                if cvi is not None and self.mask[i]:
                    #issue : curve does not exist yet:
                    self.curve.dyadic_refinement()
                    boundslist.append(
                        self.curve.active_projective_bounds_given_cv_index(cvi) 
                                )
                    del(self.curve.children)
                    self.curve.children = None
            #
            #******************************************************************
            # refine
            lenbds = len(boundslist)
            if lenbds>0:
                if self.child is  None:
                    #formerly called self.prolongate
                    bb = BoundsList(boundslist)
                    bb = BoundsList( bb.make_simply_connected() )
                    print 'refining Lspline'
                    self.new_prolongation_level(
                            boundslist=bb)#,
                            #THBactiveset = self.refine_by_these_grad_indices
                            #        )
                    #bb=self.child.curve.bounds.bounds
                    #self.child.bounds = BoundsList(bb)
                    
                else:
                    self.update_child() #update fixity
        return
    def refine_by_force(self, boundslst = None):
        """Use self.new_prolongation_level to build up 
        a Lspline solver for a new dyadily refined curve
        
        notes
        ----------
        it DOES matter which level you call refine with
        
        because self.refine_by_these_grad_indices
        knows only about solving it has done on its level
        
        
        
        
        """
        if boundslst is None:
            boundslst = [eia(0.,1.)]
        print 'forced refine called...',boundslst
        #if self.refine_by_these_grad_indices is None:
        #   refine_by_these_grad_indices = []
        
        if self.child is  None:
            bb = BoundsList(boundslst)
            #bb.make_simply_connected()
            bb = BoundsList( bb.make_simply_connected() )
            print 'refining Lspline'
            self.new_prolongation_level(
                    boundslist=bb)#,
                    #THBactiveset = self.refine_by_these_grad_indices
                    #        )
            #print 'child bounds = ',self.child.curve.bounds
        else:
            self.update_child() #update fixity
        return
    
    def change_bounds(self,refine_by_these_grad_indices=None):
        #originalbounds = self.curve.bounds
        print 'before change bounds, bounds = '
        print self.curve.bounds
        
        if refine_by_these_grad_indices is None:
            refine_by_these_grad_indices = self.refine_by_these_grad_indices
        boundslist = []
        for i in refine_by_these_grad_indices:
            cvi = self.get_vertex_index_by_any_solver_component_index(i)
        
            if cvi is not None and self.mask[i]:
                boundslist.append(
                        #self.curve.active_projective_bounds_given_cv_index(cvi) #change bounds does not look at the next level
                        self.curve.active_bounds_given_cv_index(cvi)  #change bounds looks at this level
                            )
        if boundslist:
            #self.curve.bounds = BoundsList(boundslist)
            #
            for bds in boundslist:
                self.curve.bounds.bounds.append(bds)
            #
            self.curve.bounds = BoundsList(
                    self.curve.bounds.make_simply_connected() )
        print 'after change bounds, bounds = '
        print self.curve.bounds
        return
    
#    def get_refinement_indices(self, refine_by_these_grad_indices=None):
#        """
#        use 
#        self.refine_by_these_grad_indices = np.argpartition(self.f.grad,-4)[0,-4:]
#        
#        recursively
#        to always get the 3 or fewer real 
#        independent control vetices for 
#        refinement
#        """
#        if refine_by_these_grad_indices is None:
#            refine_by_these_grad_indices = self.refine_by_these_grad_indices
#            boundslist = []
#        boundslist = []
#        for i in refine_by_these_grad_indices:
#            cvi = self.get_vertex_index_by_any_solver_component_index(i)
#        
#            if cvi is not None and self.mask[i]:
#                boundslist.append(self.curve.active_projective_bounds_given_cv_index(cvi) 
#                            )
#        #
#        return
    
    def update_child(self, THBactiveset=None):
        """update fixities
        -hybrid thb solver
        """
        if THBactiveset is None:
            THBactiveset=self.refine_by_these_grad_indices
        assert(self.child is not None),'ERROR: called update child - but Lspline has not been refined'
        #if self.child is not None:
        print 'setting fixties of child curve'
        print 'hybrid solve - fixity where THB truncates'
        self.child.set_fixity()
                #THBactiveset=THBactiveset)
        return
    
    def new_prolongation_level(self, 
                               boundslist,
                               THBactiveset=None,
                               newfixity=None):
        """adds a NEW prolongation LEVEL to Lspline
        IntervalLagrangeSpline Class
        
        -different than merely projecting vectors up to a finer,
        pre-existing level
        
        notes
        ---------
        This adds another fine level to the Lspline level system
        it does not prolongate the solution!
        
        may need to rename to 'add_prolongation_level'
        or similar.
        (same for Lagrangian.prolongation and FPD.prolongation)
        
        -Set to work from any level - it will sound to finest
        and then refine from there
        """
        #if THBactiveset is None:
        #    THBactiveset = self.refine_by_these_grad_indices
        
        thisCurve = self.curve.get_finest_curve()
        newcurve = thisCurve.dyadic_refinement()
        
        #print 'ini bounds = ',boundslist
        #boundslist.make_simply_connected()
        boundslist = BoundsList( 
                boundslist.make_simply_connected() )
        newcurve.bounds = boundslist
        #print 'ini child bounds = ',newcurve.bounds
        
        thisLspline = self.get_finest_Lspline()
        assert(thisLspline.curve is thisCurve),'Lagrangian Levels do not match THBcurve levels!'
        LgChild = thisLspline.Lagrangian.new_prolongation_level(
                                            refined_curve=newcurve)
        #thisLg.child is LChild 
        
        #print  'newcurve.bounds = ',newcurve.bounds
        Lspline = IntervalLagrangeSpline(newcurve, 
                                         LgChild, 
                                         data = {},
                                         parent=self)#,
                                         #THBactiveset=THBactiveset)
        #print 'Lspline.curve.bounds = ',Lspline.curve.bounds
        #quick fix: (NO!)
        #Lspline.mask = Lspline.masker(Lspline.Lagrangian) #fails!
        Lspline.error_code = 'None'
        self.child = Lspline
            
        return
    
    def get_finest_Lspline(self):
        """get finest level
        """
        if self.child is None:
            return self
        else:
            return self.child.get_finest_Lspline()
    
    def get_coarsest_Lspline(self):
        """get coarsest level
        """
        if self.parent is None:
            return self
        else:
            return self.parent.get_coarsest_Lspline()
        
        
    def create_mask(self):
        self.mask = self.masker(self.Lagrangian)
        return
    
    
    def set_fixity(self, thbactive=None, new_fix=None):
        """set fixed any inactive thb vertices
        
        BUG here:  self.mask[:] = False
        is setting the constraints to null
        and so on
        
        notes
        ----------
        -this is for the hybrid THB solver
        that uses fixities to avoid moving THB vertices
        on levels where they are not active.
        
        must be called before self.instate_AD_vertices
        or the vertices will/or/could be wrong in the solver!
        
        avoid this by using a true projective update of the 
        vertices throughout the hierarchy.
        """
        self.compute_fixity_mask_and_vertices(self.Lagrangian)
        if thbactive is not None:
            for el in thbactive:
                assert(el<self.curve.n*self.curve.dim),'ERROR, {} is not a vertex index'.format(el)
            #self.mask[:] = False
            #self.mask[thbactive] = True
        if new_fix is not None:
            for el in new_fix:
                self.mask[el] = False
        return
    
    def get_ivertices(self):
        return [self.curve.interval_xpts,self.curve.interval_ypts]
    
    def get_live_vertex_indicies(self):
        return np.asarray(range(self.n))[self.mask]
    
    def masker(self, Lagrangian=None):
        if Lagrangian is None:
            Lagrangian = self.Lagrangian
            
        curve  = self.curve
        #nfv    = self.nfv#curve.n*curve.dim 
        #nec    = self.nec#len(Lagrangian.equality)
        #nic    = self.nic#len(Lagrangian.minimum) + len(Lagrangian.maximum)

        """    
            Check AD implementation for a better way... 
            issue:  e.g. xpts dot basis .... must make sense
        """
        mask=[]
        
        for i in range(self.n):#nfv+nec+2*nic):
#            if (i==0 or i==(curve.n-1) or i==(curve.n) or i==(2*curve.n-1)):
#                mask.append(0)
#            else:
#                mask.append(1)
                
            #if (i+1 % curve.n == 0 or \
            #    i+1 %curve.n+1 and i<=curve.n*curve.dim:
            if curve.dim == 3:
                if (i==0              or \
                    i==(curve.n-1)    or \
                    i==(curve.n)      or \
                    i==(2*curve.n-1)  or \
                    i==(2*curve.n)    or \
                    i==(3*curve.n-1)):
                    mask.append(0)
                else:
                    mask.append(1)
            elif curve.dim == 2:
                if (i==0             or \
                    i==(curve.n-1)   or \
                    i==(curve.n)     or \
                    i==(2*curve.n-1)):
                    mask.append(0)
                else:
                    mask.append(1)
        
        #for j in range(curve.dim):
        #    for i in range(self.n):
        
        # institute curvedegeneracy for knuckles and cusps.
        # see Marcus Bole's ship design optimization work
        if 'knuckle' in self.curve.data:
            for knuckle_dict in self.curve.data['knuckle']:
                #knuckle_vertex = knuckle_dict['freevertex']
                for move_vertex in knuckle_dict['coincident']:
                    mask[move_vertex]=0
                    mask[move_vertex+curve.n]=0
                    
        #
        #*************************************
        # do fixity masks
        for fx in Lagrangian.fix:
            ff = Lagrangian.fix[fx]
            raw_index = ff.index_of_vertex
            if ff.type is 'xfix':
                mask[raw_index] = 0
            elif ff.type is 'yfix':
                print raw_index
                mask[curve.n+raw_index] = 0
            elif ff.type is 'zfix':
                mask[2*curve.n+raw_index] = 0
            
            
        mask = np.array(mask, dtype=bool)
        return mask
    
    def compute_fixity_mask_and_vertices(self, Lagrangian):
        #
        #*************************************
        # do fixity values
        for fx in Lagrangian.fix:
            ff = Lagrangian.fix[fx]
            raw_index = ff.index_of_vertex
            value = ff.given
            #print 'compute_fixity_mask_and_vertices got fix value = ',value,' type = ',ff.type
            if ff.type is 'xfix':
                self.curve.vertices[raw_index,0] = value
            elif ff.type is 'yfix':
                self.curve.vertices[raw_index,1] = value
            elif ff.type is 'zfix':
                self.curve.vertices[raw_index,2] = value
        return
        
    
    def map_variables(self):
        """Map a dim x n array of vertices
        to a n*dim vector of vertices
        
        
        """
        self.vertex_map = {}
        for j in range(self.curve.dim):
            for i in range(self.curve.n):
                self.vertex_map[i+j*self.curve.n] = (j,i)
        return
    
    
    def map_index_to_vertex(self, i):
        return self.vertex_map[i] #i%self.n
    
    
    def get_vertex_index_by_any_solver_component_index(self,i):
        """if you find that a gradient component is above 
        threshold for refinement,
        this tells you which vertex/basis to refine.
        def
        Then use 
        
        rbspline.active_projective_bounds_given_cv_index(i=vtuple[1])
        
        to get the proper bounds to tell the 
        THB curve evaluator to evaluate this one as refined.
        
        """
        #mssg = 'ERROR: attempted to select lagrangian constraint instead of vertex index'
        #assert(i<self.curve.n*self.curve.dim),mssg
        if i<self.curve.n*self.curve.dim:
            vtuple = self.map_index_to_vertex(i)
            return vtuple[1]
        else:
            return None
         
    
    def get_vertex_by_any_index(self,i):
        rows = range(self.curve.dim)
        vtuple = self.map_index_to_vertex(i)
        
        X = []
        for d in rows:
            X.append(self.vertices[d][vtuple[1]] )
        return np.asarray(X)
        
    
    def instate_AD_vertices(self):
        """ Define the AD conrol point and constraint parmeters"""
        ptslist = self.Lagrangian.setup_derivatives()
        if self.curve.dim == 2:
            self.curve.xpts, self.curve.ypts = ptslist
            self.curve.ADvertices = ptslist
        elif self.curve.dim ==3:
            self.curve.xpts, self.curve.ypts,self.curve.zpts = ptslist
            self.curve.ADvertices = ptslist
        else:
            self.curve.ADvertices = ptslist
            print 'Work in Progress, FAILURE, '
            print 'this program only works in 2D or 3D, '
            print 'get started fixing that here!'
        return
        
    
    def instate_interval_vertices(self, i_min = None, i_max = None):
        """
            instantiate a set of interval design variables
            using the fuzzyNumber class
        """
        if i_min is not None and i_max is not None:

            self.curve.interval_xpts, self.curve.interval_ypts = self.Lagrangian.setup_simple_interval(imin=i_min, imax=i_max)
        else:
            if self.fjc:
                self.curve.interval_xpts, self.curve.interval_ypts = self.Lagrangian.interval_Lagrange_ranged( self.curve,
                                                                                                  self.data['xmin'],self.data['xmax'],
                                                                                                  self.data['ymin'], self.data['ymax'],
                                                                                                  self.data['lmin'], self.data['lmax'],
                                                                                                  self.data['fjcmin'], self.data['fjcmax'])
            else:
                self.curve.interval_xpts, self.curve.interval_ypts = self.Lagrangian.interval_Lagrange_ranged( self.curve,
                                                                                              self.data['xmin'],self.data['xmax'],
                                                                                              self.data['ymin'], self.data['ymax'],
                                                                                              self.data['lmin'], self.data['lmax'])
        return
    
        
    def extreme_C0(self):
        """extremes
        """
        if self.curve.dim == 2:
            print '2D extemes'
            xmin = self.curve.interval_xpts[0].value.inf
            xmax = self.curve.interval_xpts[0].value.sup
            ymin = self.curve.interval_ypts[0].value.inf
            ymax = self.curve.interval_ypts[0].value.sup
            for el in self.curve.interval_xpts[1:]:
                xmin = min(xmin, el.value.inf)
                xmax = max(xmax, el.value.sup)
            for el in self.curve.interval_ypts[1:]: ##added 10-5-2015 - probable bug fix
                ymin = min(ymin, el.value.inf)
                ymax = max(ymax, el.value.sup)
            return xmin,xmax,ymin,ymax
        elif self.curve.dim == 3:
            print '3D extremes'
            xmin = self.curve.interval_xpts[0].value.inf
            xmax = self.curve.interval_xpts[0].value.sup
            ymin = self.curve.interval_ypts[0].value.inf
            ymax = self.curve.interval_ypts[0].value.sup
            zmin = self.curve.interval_zpts[0].value.inf
            zmax = self.curve.interval_zpts[0].value.sup
            for el in self.curve.interval_xpts[1:]:
                xmin = min(xmin, el.value.inf)
                xmax = max(xmax, el.value.sup)
            for el in self.curve.interval_ypts[1:]: ##added 10-5-2015 - probable bug fix
                ymin = min(ymin, el.value.inf)
                ymax = max(ymax, el.value.sup)
            for el in self.curve.interval_zpts[1:]: ##added 10-5-2015 - probable bug fix
                zmin = min(zmin, el.value.inf)
                zmax = max(zmax, el.value.sup)
            return xmin,xmax,ymin,ymax,zmin,zmax
        else:
            print 'interval extremes not presently supported for {(N)D, N>3}'
            return 
            
        
    def contract_constraints(self, vertices=None, L=None, update=True, key=None):
        """
            Hansen hull and box conistency, page 212
        """
        if  vertices is None:
            vertices = [self.curve.interval_xpts,
                        self.curve.interval_ypts]
        if  L is None:
            L = self.Lagrangian
        
        if key is None:
            for key in L.equality:
                equality_constraint = L.equality[key]
                if equality_constraint.has_contractor:
                    print equality_constraint.type
                    constraint_value = equality_constraint.pass_value
                    vertices = equality_constraint.contractor(vertices, constraint_value)
        else:
            equality_constraint = L.equality[key]
            if equality_constraint.has_contractor:
                print equality_constraint.type
                constraint_value = equality_constraint.pass_value
                vertices = equality_constraint.contractor(vertices, constraint_value)
        if update == True:
            self.curve.interval_xpts = vertices[0]
            self.curve.interval_ypts = vertices[1]
            return
        else:
            return vertices
    
    def get_sliced_vertices(self, index, slice_low, slice_high, xy='x', vertices=None):
        if  vertices is None:
            vertices = [self.curve.interval_xpts,
                        self.curve.interval_ypts]
        xy = {'x':0, 'y':1}[xy]
        v = copy.deepcopy(vertices)
        v[xy][index].value.inf = slice_low
        v[xy][index].value.sup = slice_high
        return v
            
    def get_objective(self, key, vertices=None, L=None):
        """
            key = the key of the objective
            
            returns : the obj f
            
        """
        if  vertices is None:vertices = [self.curve.interval_xpts,
                                         self.curve.interval_ypts]
        if  L is None:L = self.Lagrangian
        pass_value  = L.obj[key].pass_value
        fp          = L.obj[key].given
        kind        = L.obj[key].kind
        return L.obj[key].computeUpdate(kind, vertices, pass_value) -fp
        
    def get_constraint(self, key, vertices=None, L=None):
        """
            key = the key of the constraint
            
            returns : the constraint h where
                      h = (actual_value - form_parameter)
            
            a consistent contraint includs 0.
        """
        if  vertices is None:vertices = [self.curve.interval_xpts,
                                         self.curve.interval_ypts]
        if  L is None:L = self.Lagrangian
        pass_value  = L.equality[key].pass_value
        fp          = L.equality[key].given
        kind        = L.equality[key].kind
        return L.equality[key].computeUpdate(kind, vertices, pass_value) -fp
        

        
    def get_partial_constraint(self, key, vertices = None, L=None):
        if  L is None:
            L = self.Lagrangian
        if  vertices is None:
            vertices = [self.curve.interval_xpts,
                        self.curve.interval_ypts]
        #pass_value  = L.equality[key].pass_value
        fp          = L.equality[key].given
        kind        = L.equality[key].kind
        f = partial(L.equality[key].computeUpdate, kind )
        def lf(vertices, fmp):
            return f(vertices)-fmp
        return partial(lf,fmp=fp)
    
    def monotonicity_constraint(self, key, vertices=None, L = None):
        """
            key = the key of the constraint
        """
        if  vertices is None:
            vertices = [self.curve.interval_xpts,
                        self.curve.interval_ypts]
        if  L is None:
            L = self.Lagrangian
        
        f1 = self.get_partial_constraint(key)
        nv_min, nv_max = self.ia.monotonic_contraction(f1(vertices),
                                                       vertices, 
                                                       self.vertex_map)
        f_standard = f1(vertices)
        fmin = f1(nv_min)
        self.cmin = self.curve.area
        fmax = f1(nv_max)
        self.cmax = self.curve.area
        
        d1 = copy.deepcopy(f_standard)
        d1.value.inf = fmin.value.inf
        d1.value.sup = fmax.value.sup
        """
        for i in range(self.n):
            d1.grad[0,i].inf = min(fmin.grad[0,i].inf,
                                    fmin.grad[0,i].sup,
                                    fmax.grad[0,i].inf,
                                    fmax.grad[0,i].sup)
            d1.grad[0,i].sup = max(fmin.grad[0,i].inf,
                                    fmin.grad[0,i].sup,
                                    fmax.grad[0,i].inf,
                                    fmax.grad[0,i].sup)
            for j in range(self.n):
                d1.hess[i,j].inf = min(fmin.hess[i,j].inf,
                                        fmin.hess[i,j].sup,
                                        fmax.hess[i,j].inf, 
                                        fmax.hess[i,j].sup)
                d1.hess[i,j].sup = max(fmin.hess[i,j].inf,
                                        fmin.hess[i,j].sup,
                                        fmax.hess[i,j].inf, 
                                        fmax.hess[i,j].sup)
        #"""
#                
        return d1#ia(fmin.inf, fmax.sup)
        
        
        
    def box_consistency(self, key, 
                        vertex_tpl,
                        vertices=None, 
                        L = None, max_it = 5): #dino bc
        """
            key = the key of the constraint

            current structure is inefficient:
                    you grab fp, etc, every i,j pair
                    but it's constant up to key: 
                        therefore you're grabbing the same value each time.
        """
        if  vertices is None:
            vertices = [self.curve.interval_xpts,
                        self.curve.interval_ypts]
        if  L is None:
            L = self.Lagrangian
            
        pass_value  = L.equality[key].pass_value
        fp          = L.equality[key].given
        kind        = L.equality[key].kind
        def func(V,C):
            return L.equality[key].computeUpdate(kind, V, pass_value) - C
            
        return self.ia.compute_scalar_newton(func, fp, 
                                             vertices, vertex_tpl, max_it=5)
                                             
    def monotonic_box_consistency(self, key, 
                                  vertex_tpl,
                                  vertices=None, 
                                  L = None, max_it = 5): #dino bc
        """
            key = the key of the constraint
        """
        if  vertices is None:
            vertices = [self.curve.interval_xpts,
                        self.curve.interval_ypts]
        if  L is None:
            L = self.Lagrangian
            
        pass_value  = L.equality[key].pass_value
        fp          = L.equality[key].given
        kind        = L.equality[key].kind
        
        func = partial(self.monotonicity_constraint,key)
            
        return self.ia.montonic_compute_scalar_newton(func, 
                                                 vertices, 
                                                 vertex_tpl, 
                                                 max_it=5)
        
    def compute_box_consistency(self, 
                                which = None, 
                                xyz = [0,1], 
                                verts=[1,2,3,4,5],
                                max_it = 5, 
                                box_type = 'natural'):
        """
            TODO:  usemask to iterate
            over the live vertices
        """
        if which is None:
            keys = self.Lagrangian.equality.keys()
        else:
            keys = which
            
        if box_type == 'm':
            for key in keys:
                try:
                    print 'monotonic box consistency : ', key, self.Lagrangian.equality[key].type
                    for i in xyz:
                        for j in verts: # TODO: use MASK!
                            vertices = self.monotonic_box_consistency(key, (i,j),
                                                            [self.curve.interval_xpts,
                                                             self.curve.interval_ypts],
                                                             self.Lagrangian
                                                             )
                        self.curve.interval_xpts = vertices[0]
                        self.curve.interval_ypts = vertices[1]
                except:
                    pass
        elif box_type == 'natural':
            for key in keys:
                try:
                    print 'standard box consistency : ', key, self.Lagrangian.equality[key].type
                    for i in xyz:
                        for j in verts: # TODO: use MASK!
                            vertices = self.box_consistency(key, (i,j),
                                                            [self.curve.interval_xpts,
                                                             self.curve.interval_ypts],
                                                             self.Lagrangian
                                                             )
                        self.curve.interval_xpts = vertices[0]
                        self.curve.interval_ypts = vertices[1]
                except:
                    pass
        return
    

        
    
    def optimize(self, vertices=None, stop=None, Lagrangian=None):
        self.error_code = None
        if stop is None:
            Inumber = self.Inumber
        else:
            Inumber   =  stop
        if vertices   is None: vertices = self.curve.ADvertices
        if Lagrangian is None : Lagrangian = self.Lagrangian
        count           = 0
        self.conv       = 1.0
        self.weight_f   = 1.
        self.delta_z    = 0.
        self.history_Xc = [0.]
        while (self.conv>self.tol and count<Inumber):
            self.count = count
            self.f      = self.compute_lagrangian(vertices, Lagrangian)
            self.delta_z     = self.invertReduced(self.f)
            if self.delta_z is 'NAN': #need to clean up this bail out!
                self.error_code = 'NAN'
                print 'NAN' #this nan is defined as a string in function def invertReduced()
                return 'NAN'
            vertices    = self.upDateVerticesReduced(self.delta_z, 
                                                     vertices, 
                                                     Lagrangian)
            self.conv   = self.ReducedConvergenceCriteria(self.f)
            count += 1
            if self.verbose: print ':',count
        
        if count == Inumber:
            print 'WARNING: Lspline max iteration stop'
            self.error_code = 'max_iter'
        if np.isnan(self.conv): #func ReducedConvergenceCriteria is math -> returns true nan if system is singular
            self.nitr = count
            print 'ERROR, NAN encountered in solution returning to original Coordinates'
            self.curve.vertices = self.IniVertices
            self.error_code = 'NAN'
            return 'NAN'
        else:
            self.nitr = count
            self.vertices = vertices
            self.update_Bspline(vertices)
            
        self.postprocess = SolverPostProcessor(self)
        return #vertices
    
    
    def optimize_thb(self, vertices=None, 
                     stop=None, Lagrangian=None,
                     relax_grad_bounds = False,
                     restart=False,
                     ceiling = 1.e-10,
                     floor = None,
                     min_iter_Dbds=3):
        """Recipe for a THB solver
        
        Goal: eval as THB (just evaluate vertices projected to the finest level) during solve
        -some of the curve should be able to remain low-res
        -use gradient to mark for refinement.. (some facility already exists)
        (i.e. for real use at the refined level)
        
        There is no need to project the gradient back 
        when you may simply deform the fine level vertices and project those back.
        -except maybe choosing which fine level vertices to refine...?
        
        
        dev
        ----------
        self = Lspline.get_finest_Lspline()
        Lagrangian = self.Lagrangian   
        vertices = self.curve.ADvertices
        
        relax_grad_bounds = False
        restart=False
        ceiling = 1.e-2
        floor = None
        do_thb=True
        """
        #assert(self.child is None),'Error: THB optimization must be called on finest level'
            
        self.error_code = None
        if stop is None:
            Inumber = self.Inumber
        else:
            Inumber   =  stop
        if vertices   is None: vertices = self.curve.ADvertices
        if Lagrangian is None : Lagrangian = self.Lagrangian
        if floor is None: floor = self.tol
        count           = 0
        self.conv       = 1.0
        self.weight_f   = 1.
        self.delta_z    = 0.
        self.history_Xc = [0.]
        
        badgrad = 0 #if things get worse over time, 
        #break outside the solver to smooth
        csldbdx = 0
        thbvertices = None
        while (self.conv>self.tol and count<Inumber):
            if self.child is not None:
                this = self.get_finest_Lspline()
                return this.optimize_thb(vertices,stop,Lagrangian,
                                         relax_grad_bounds,
                                         restart,
                                         ceiling,
                                         floor)
            #else we are on the right level.  Keep going.
            thbvertices = self.curve.project_vertices_all_active()
            vertices = self.push_vertices_to_AD_values(thbvertices)
            """ #now what!:
            
            -build a mask for the various levels?
            -perturb vertices at this lowest level based on the
            THB bounds.... something like that.
            
            case1:
                #THB total projection:
                thb_vertices = Lspline.curve.project_vertices_sensible(.5)
                b1 = Lspline.curve.children.TLMBasisFuns(.5)
            
            case2:
                #Boier-Martin
                thb_basis = Lspline.curve.rthb_basis(.5)
                cvl0 = Lspline.curve.vertices
                cvl1 = Lspline.curve.children.vertices
                #(the question with this one has always been
                # how does one extend it to surfaces?)
            
            case3!: 
                #THB total projection
                thb_vertices = Lspline.curve.project_vertices_all_active()
                b1 = Lspline.curve.children.TLMBasisFuns(.5)
                # 
                # this is the really one that competes with case 2
                # (and is clearly compatible with THBsurfaces)
                
            #"""
            """
            for case 1, thb_vertices are UNAFFECTED by choice of u
            """
            self.count = count
            self.f      = self.compute_lagrangian(vertices, Lagrangian)
            #**************************************************************
            # THB gradient - initial 
            #**************************************************************
            # save the initial gradient for global comparison
            if not restart and count == 0:
                self.initialgrads = self.f.grad[0].T
            #**************************************************************
            #**************************************************************
            self.delta_z     = self.invertReduced(self.f)
            if self.delta_z is 'NAN': #need to clean up this bail out!
                self.error_code = 'NAN'
                print 'NAN' #this nan is defined as a string in function def invertReduced()
                return 'NAN'
            vertices    = self.upDateVerticesReduced(self.delta_z, 
                                                     vertices, 
                                                     Lagrangian)
            self.conv   = self.ReducedConvergenceCriteria(self.f)
            count += 1
            csldbdx +=1
            if self.verbose: print ':',count
            #**************************************************************
            # THB gradient : monitoring changes
            #**************************************************************
            #
            #
            #elementwise maximum of the vertex mapped elements of the graadient:
            #self.tempgrads = np.asarray(
            #        np.maximum(self.maxgrads,abs(self.f.grad[0].T)  ) )
            #of the vertex mapped elements of the graadient:
            self.tempgrads = np.asarray(
                    np.maximum(self.maxgrads,abs(self.f.grad[0,:self.curve.n*self.curve.dim].T)  ) )
            #
            #elementwise > less than
            self.tempgreater = np.greater(abs(self.f.grad[0,:self.curve.n*self.curve.dim].T),self.maxgrads)
            #***
            #sort gradient indices in order of largest
            #refine_by_these_grad_indices = np.argsort(abs(self.f.grad[0,:self.curve.n*self.curve.dim]))[::-1]
            #now threshold them
            #
            dobreak=False
            refine_by_these_grad_indices = np.flatnonzero(abs(self.tempgrads)>ceiling)

            if len(refine_by_these_grad_indices)>0:
                # sort indices by magnitude of gradient
                self.refine_by_these_grad_indices = np.argsort(abs(self.f.grad[0,refine_by_these_grad_indices]))[0,::-1]
            
                self.refine_by_these_grad_indices = np.asarray(self.refine_by_these_grad_indices)
                if np.shape(self.refine_by_these_grad_indices)[1]>8:
                    self.refine_by_these_grad_indices = list(np.asarray(self.refine_by_these_grad_indices[0,:8]))
                    #self.refine()
                    #break
                else:
                    self.refine_by_these_grad_indices = list(np.asarray(self.refine_by_these_grad_indices[0,:]))
                    boundslist = []
                    for i in self.refine_by_these_grad_indices:
                        
                        cvi = self.get_vertex_index_by_any_solver_component_index(i)
        
                        if  cvi is not None and self.mask[i]:
                            boundslist.append(i)
                    self.refine_by_these_grad_indices = boundslist
                
                if len(self.refine_by_these_grad_indices)>0:
                    if self.parent is not None and self.child is None:
#                        print 'Change Bounds'
#                        self.change_bounds(self.refine_by_these_grad_indices)
#                        print self.curve.bounds
                        #if self.curve.bounds[0].inf == 0. and self.curve.bounds[1].sup == 1.:
                        #    break
                        if csldbdx>=min_iter_Dbds:
                            csldbdx = 0
                            #print 'Change Bounds'
                            #self.change_bounds(self.refine_by_these_grad_indices)
                            #print self.curve.bounds
                            #if self.curve.bounds.contains(ia(0.,1.)):
                            if self.curve.bounds[0].contains(ia(0.001,.999)):
                                pass
                            else:
                                print 'Change Bounds'
                                self.change_bounds(self.refine_by_these_grad_indices)
                                print self.curve.bounds
#                                print 'Break'
#                                break
#                            else:
#                                print 'Change Bounds'
#                                self.change_bounds(self.refine_by_these_grad_indices)
#                                print self.curve.bounds
#                                #    self.refine()
#                                #    break
                    #else:
                    #    Inumber = count
                    #    print 'gradient-interrupt.  count = ',count
                    #    self.optimize_iterations = count
            #**************************************************************
            #
            
            #**************************************************************
            # project updated vertices back to their levels
            #**************************************************************
            self.update_thb_all_levels(vertices)
            #**************************************************************
            
        if count == Inumber:
            print 'Max Iteration Reached'
            self.error_code = 'max_iter'
        if np.isnan(self.conv): #func ReducedConvergenceCriteria is math -> returns true nan if system is singular
            self.nitr = count
            print 'ERROR, NAN encountered in solution returning to original Coordinates'
            self.curve.vertices = self.IniVertices
            self.error_code = 'NAN'
            return 'NAN'
        else:
            self.nitr = count
            self.vertices = vertices #update AD list of vertices
            #**************************************************************
            # project updated vertices back to their levels
            # instead of the usual update
            #**************************************************************
            self.update_thb_all_levels(vertices)
            #**************************************************************
            
        #if self.conv>self.tol:
        self.postprocess = SolverPostProcessor(self)
        return #vertices
    
    def interval_newton(self, itr = 5):
        self.interval_optimize(stop = itr)
        return
    def interval_optimize(self, vertices=None, 
                          stop=None, L = None, solver = 'N',
                          test_check = None, point = None, 
                          check_IA_correctness=False, sor=None):
        
        if stop is None:
            Inumber = self.Inumber
        else:
            Inumber = stop
        if  vertices is None: vertices = [self.curve.interval_xpts,
                                          self.curve.interval_ypts]
        if  L is None:L = self.Lagrangian
        
        if sor is not None:
            self.sor = sor
        self.check_IA_correctness = check_IA_correctness  
        count       = 0
        self.conv   = 0
        #self.ia = IA(self.mask) #attaching interval class - maybe make this a singleton? - move to mask instantiation point
        ## we coud choose our algorithms at runtime...
        # anyway, 
        while (self.conv<(self.n-4) and count<Inumber):
            self.count = count
            if point is None:
                point = .5
            self.v_old = copy.deepcopy(vertices)
            
            #for point in [0.25,.5,.75]:
            # Hansen starting interval, etc.. page? 
            thin_vertices, thin_Lagrange = self.ia.get_thin_expansion_point(point, vertices, L)#self.Lagrangian)
            thin_vertices = np.asarray(thin_vertices)
            
            # Normalize Fritz John Conditions
            #thin_Lagrange = self.norm_lagrangian_multipliers(thin_Lagrange)
            #L = self.norm_lagrangian_multipliers(L)
            
            try:
                # traditional f:
                self.thin_f = self.compute_interval_lagrangian(thin_vertices, thin_Lagrange) #L) #
                #self.f      = self.compute_interval_lagrangian(vertices, L)
                
                #monotonic f:
                #self.thin_f = self.compute_monotonic_interval_lagrangian(thin_vertices, thin_Lagrange) #L)##no point?
                self.f      = self.compute_monotonic_interval_lagrangian(vertices, L)
            except:
                print 'error in functional computation'
                #self.thin_f = 
            
            #inersion:
            if solver == 'N':
                interval_z  = self.ia.compute_interval_newton( self.f, 
                                                              self.thin_f, 
                                                              L, 
                                                              vertices,
                                                              thin_vertices,
                                                              thin_Lagrange,
                                                              fjc = self.fjc)
                # using 77 test newton2
#                interval_z  = self.ia.compute_interval_newton3( self.f, 
#                                                              self.thin_f, 
#                                                              L, 
#                                                              vertices,
#                                                              thin_vertices,
#                                                              thin_Lagrange)
                #print interval_z
                self.xnx = interval_z
                vertices, self.Lagrangian = self.IntervalupDateVertices(interval_z, vertices, L)
            elif solver =='K':
                interval_z  = self.ia.interval_krawczyk_inversion(self.f, L, vertices)
                vertices, self.Lagrangian = self.IntervalupDateVertices(interval_z, vertices, L)
            else:
                interval_z  = self.ia.interval_krawczyk_inversion(self.f, L, vertices)
                vertices, self.Lagrangian = self.IntervalupDateVertices(interval_z, vertices, L)
                interval_z  = self.ia.compute_interval_newton(self.f, self.thin_f, L, vertices)
                vertices, self.Lagrangian = self.IntervalupDateVertices(interval_z, vertices, L)
            
            count += 1
            print ':',count
        if np.isnan(self.conv):
            self.nitr = count
            self.curve.vertices = self.IniVertices
            print 'Warning, NAN encountered in solution returning to original Coordinates'
        else:
            self.nitr = count
            """ssor could be detrimental:  extra interval operations
            can be good or bad!  This needs more thought, or to be checked out.
            """
            for i in range(self.curve.dim):
                for j in range(self.curve.n):
                    vertices[i][j] = vertices[i][j]*self.sor +self.v_old[i][j]*(1.-self.sor)
            self.interval_vertices = vertices
            self.update_Bspline(vertices)
        return
        
    
    
        
    
            
    def upDateVerticesReduced(self, delta_z, vertices, Lagrangian):
        """
            Function to update the "Reduced Matrix" vertices and Lagrange multipliers
        """

        count   = 0  #denotes which lagrangian variable are we updating
        #for j,verts in zip(range(self.curve.dim),[self.curve.xpts,self.curve.ypts]):
        for j in range(self.curve.dim):
            for i in range(self.curve.n):
                if self.mask[i+j*self.curve.n]==True:
                    vertices[j][i].value -= float(delta_z[self.mask][count])
                    count +=1
        
                
        for el in Lagrangian.equality:
            ei = Lagrangian.equality[el].index
            Lagrangian.equality[el].Lagrange.value -= float(delta_z[ei])
            count +=1
        for el in Lagrangian.minimum:
            if self.mask[count]==True:
                li = Lagrangian.minimum[el].index
                si = Lagrangian.minimum[el].slack_index
                Lagrangian.minimum[el].Lagrange.value -= float(delta_z[li])
                Lagrangian.minimum[el].slack.value -= float(delta_z[si])
            count +=1
        for el in Lagrangian.maximum:
            if self.mask[count]==True:
                li = Lagrangian.maximum[el].index
                si = Lagrangian.maximum[el].slack_index
                Lagrangian.maximum[el].Lagrange.value -= float(delta_z[li])
                Lagrangian.maximum[el].slack.value -= float(delta_z[si])
            count +=1
        if self.fjc:
            self.Lagrangian.fj.Lagrange.value -= float(delta_z[-1]) #ei almost looks good :(
            #self.Lagrangian.fj.Lagrange.value = abs(self.Lagrangian.fj.Lagrange.value)
            count +=1
            
        if 'knuckle' in self.curve.data: #then the requisite points better be masked out or its all junk!
            for knuckle_dict in self.curve.data['knuckle']:
                knuckle_vertex = knuckle_dict['freevertex']
                for move_vertex in knuckle_dict['coincident']:
                    #self.curve.vertices[move_vertex,0] = self.curve.vertices[knuckle_vertex,0]
                    self.curve.xpts[move_vertex].value = self.curve.xpts[knuckle_vertex].value
                    #self.curve.vertices[move_vertex,1] = self.curve.vertices[knuckle_vertex,1]
                    self.curve.ypts[move_vertex].value = self.curve.ypts[knuckle_vertex].value
        return vertices
        
        
    def OLDIntervalupDateVertices(self, z_list, vertices, Lagrange):
        """DONE: check for zero width updates 
        - do not kill the interval for the next iteration.
        """
        self.conv = 0
        for i in range(1,self.curve.n-1):
                       
            iwidthx = z_list[i-1].value.width()
            if abs(iwidthx) > 0.0:
                vertices[0][i] = z_list[i-1] 
            iwidthy = z_list[i+self.curve.n-3].value.width()
            if abs(iwidthy) > 0.0:
                vertices[1][i] = z_list[i+self.curve.n-3]
               
            if abs(iwidthx) < self.tol:
                self.conv +=1
            if abs(iwidthy) < self.tol:
                self.conv +=1
                
        for el, j in zip(Lagrange.equality, range(self.curve.dim*self.curve.n, 
                         self.curve.dim*self.curve.n + self.nec) 
                        ):
            iwidth = z_list[j-4].value.width()
            if abs(iwidth) > 0.0:
                Lagrange.equality[el].interval_Lagrange = z_list[j-4]
            if abs(iwidth) < self.tol:
                self.conv +=1
        return vertices, Lagrange
    def IntervalupDateVertices(self, z_list, vertices, Lagrange):
        """DONE: check for zero width updates 
        - do not kill the interval for the next iteration.
        """
        self.conv = 0
        for i in range(1,self.curve.n-1):
                       
            iwidthx = z_list[i-1].value.width()
            test = z_list[i-1].value & vertices[0][i].value
            if test in vertices[0][i].value:
                vertices[0][i].value = test
            iwidthy = z_list[i+self.curve.n-3].value.width()
            test = z_list[i+self.curve.n-3].value & vertices[1][i].value
            if test in vertices[1][i].value:
                vertices[1][i].value = test
               
            if abs(iwidthx) < self.tol:
                self.conv +=1
            if abs(iwidthy) < self.tol:
                self.conv +=1
                
        for el, j in zip(Lagrange.equality, range(self.curve.dim*self.curve.n, 
                         self.curve.dim*self.curve.n + self.nec) 
                        ):
            iwidth = z_list[j-4].value.width()
            test = Lagrange.equality[el].interval_Lagrange.value & z_list[j-4].value
            if test in Lagrange.equality[el].interval_Lagrange.value:
                Lagrange.equality[el].interval_Lagrange.value = test
            if abs(iwidth) < self.tol:
                self.conv +=1
                
        if self.fjc:
            iwidth = z_list[-1].value.width()
            test = self.Lagrangian.fj.interval_Lagrange.value & z_list[-1].value
            if test in self.Lagrangian.fj.interval_Lagrange.value:
                self.Lagrangian.fj.interval_Lagrange.value = test
            if abs(iwidth) < self.tol:
                self.conv +=1
        return vertices, Lagrange
        
        
        
    def update_Bspline(self, vertices):
        if isinstance(vertices[0][0].value, ia):
            self.curve.interval_xpts = vertices[0]
            self.curve.interval_ypts = vertices[1]
        elif isinstance(vertices[0][0].value, float):
            new_vertices = np.zeros((self.curve.n, self.curve.dim))
            for j in range(self.curve.dim):
                for i in range(self.curve.n):
                    new_vertices[i,j] = vertices[j][i].value
            self.curve.vertices = new_vertices
        return
    
    def push_vertices_to_AD_values(self, vertices):
        for di in range(self.curve.dim):
            for thbi in range(self.curve.n):
                self.curve.ADvertices[di][thbi].value = vertices[thbi][di]
        return self.curve.ADvertices
    
    def push_AD_values_to_vertices(self, ADvertices):
        new_vertices = np.zeros((self.curve.n, self.curve.dim))
        for di in range(self.curve.dim):
            for thbi in range(self.curve.n):
                new_vertices[thbi][di] = self.curve.ADvertices[di][thbi].value
        return new_vertices
    
    def update_thb_vertices_this_level(self, vertices):
        self.curve.vertices = vertices
        
    def update_thb_all_levels(self, ADvertices):
        """for level in THBLspline
        THB-project the vertices to that level
        and update them there with update_Bspline.
        """
        #make a vertices vector that the projector expects
        vertices = self.push_AD_values_to_vertices(ADvertices)
        self.update_thb_vertices_this_level(vertices)
        this = self
        while this.parent is not None:
            vector = this.curve.restrict(vertices)
            this.parent.update_thb_vertices_this_level(vector)
            this = this.parent
        return
        
        
                
    def is_active(self, sm):
        """
            is sm is negative,
            then the constraint
            is violated
            
            -> must use strict inequality
            -> s=0 is active!
        """
        if sm >0.0: #self.tol:# 
            active = False
        else:
            active = True
        return active   
        
    def compute_interval_lagrangian(self, vertices=None, L=None, Lagrange_vector =None):
        f               = 0.
        self.store_AF = None
        self.store_AF = []
        
        if  vertices is None: vertices = [self.curve.interval_xpts,
                                          self.curve.interval_ypts]
        if  L is None:L = self.Lagrangian
        
        if self.check_IA_correctness:
            print'Interval inf < interval sup '
            print 'check Not implemented at this time'
            #            for el in L.obj:
            #                pass_value  = L.obj[el].pass_value
            #                new         = L.obj[el].computeUpdate('LS', vertices, pass_value)
            #                L.obj[el].computed_value = new
            #                if isinstance(new, ad):
            #                    self.test_hessian_consistency(new.hess)
            #                    self.test_gradient_consistency(new.grad)
            #                f           = (L.obj[el].weight*new) + f
            #                
            #            for el in L.equality:
            #                pass_value  = L.equality[el].pass_value
            #                fp          = L.equality[el].given
            #                Lambda           = L.equality[el].interval_Lagrange
            #                new = L.equality[el].computeUpdate('equality', vertices, pass_value)
            #                L.equality[el].computed_value = new
            #                self.test_hessian_consistency(new.ihess)
            #                self.test_gradient_consistency(new.igrad)
            #                f           = Lambda*(new-fp)*(-1.) +f
                            
        
        else:
            
            for el in L.obj:
                pass_value  = L.obj[el].pass_value
                new = L.obj[el].computeUpdate('LS', vertices, pass_value) 
                L.obj[el].computed_value = new
                f           = f + (L.obj[el].weight*new)
                self.store_AF.append({'el':el,
                                      'obj':new,
                                      'f':(L.obj[el].weight)*(new) })
            if self.fjc:
                f = f*L.fj.interval_Lagrange
                
            for el in L.equality:
                pass_value  = L.equality[el].pass_value
                fp          = L.equality[el].given
                Lambda      = L.equality[el].interval_Lagrange
                new = L.equality[el].computeUpdate('equality', vertices, pass_value) 
                L.equality[el].computed_value = new
                f           = f + Lambda*(new - fp)*(-1.) #was -
                self.store_AF.append({'el':el,
                                      'c':new,
                                      'f':Lambda*(new-fp), 
                                      'l':Lambda})            
        return f
        
    def compute_monotonic_interval_lagrangian(self, vertices=None, L=None, Lagrange_vector =None):
        f=0.
        self.store_AF = None
        self.store_AF = []
        
        if  vertices is None: vertices = [self.curve.interval_xpts,
                                          self.curve.interval_ypts]
        if  L is None:L = self.Lagrangian
        
        for el in L.obj:
            pass_value  = L.obj[el].pass_value
            new = L.obj[el].computeUpdate('LS', vertices, pass_value) 
            L.obj[el].computed_value = new
            f           = f + (L.obj[el].weight*new)
            self.store_AF.append({'el':el,
                                  'obj':new,
                                  'f':(L.obj[el].weight)*(new) })
        if self.fjc:
                f = f*L.fj.Lagrange
                
        for el in L.equality:
            Lambda      = L.equality[el].interval_Lagrange
            value       = self.monotonicity_constraint(el, vertices)
            L.equality[el].computed_value = value#L.equality[el].computeUpdate('equality', vertices, pass_value) 
            f           = f - Lambda*value
            self.store_AF.append({'el':el,
                                  'c':value,
                                  'f':Lambda*value, 
                                  'l':Lambda})
            
        return f
        
    def compute_lagrangian(self, vertices=None, L=None):
        """Compute and Sum every component of the 
        Lagrangian
        
        Parameters
        ----------
            vertices (optional) : vertices to be used for computation
                                    by default these will be 
                                    the Auto Diff vertices 
                                    of the curve
                                type=list of lists of 
                                each separate spactial dimension of
                                the vertices
                                each component of each list,
                                that is, each component of each 
                                control vertex, is of type=AD
            
        Returns
        ----------
            f : the result of the Lagrangian evaluation
                type = AD (auto diff) object (1 object only)
                
                
        Notes
        ----------
            the computed value will represent different things
            depending on the constraint.
            For, instance, in the case of VertexConstraints, 
            computed_value it is the difference
            between desired and actual value.
        """
        if vertices is None: vertices = self.curve.ADvertices
        if L is None : L = self.Lagrangian
        
        self.store_AD = None
        self.store_AD = []
        
        f               = 0.
        
        
        for el in L.obj:
            pass_value  = L.obj[el].pass_value
            new         = L.obj[el].computeUpdate('LS', vertices, pass_value)
            L.obj[el].computed_value = new
            f           = (L.obj[el].weight*new) + f
            self.store_AD.append({'el':el,
                                  'obj':new,
                                  'f':(L.obj[el].weight)*(new) })
        if self.fjc:
            f = f*L.fj.Lagrange
        
        
        for el in L.equality:
            pass_value  = L.equality[el].pass_value
            fp          = L.equality[el].given
            Lambda      = L.equality[el].Lagrange
            new         = L.equality[el].computeUpdate('equality', vertices, pass_value)
            L.equality[el].computed_value = new
            f           = Lambda * (new - fp)*(-1.) + f
            #f           = Lambda * (new - fp) + f
            self.store_AD.append({'el':el,
                                  'c':new,
                                  'f':Lambda * (new - fp)*(-1.) , 
                                  'l':Lambda})



        for el in L.minimum:
            li = L.minimum[el].index
            si = L.minimum[el].slack_index
            pass_value  = L.minimum[el].pass_value
            g = L.minimum[el].computeUpdate('min', vertices, pass_value)
            L.minimum[el].computed_value = g
            Lambda           = L.minimum[el].Lagrange
            slack       = L.minimum[el].slack 
            slack.value = -np.sqrt(g.value)
            # -np.sqrt(abs(g.value))
            #abs(g.value)  #-g.value
            L.minimum[el].slack.value = slack.value
            sm          = g.value#max(g.value,0.)
            active      = self.is_active(sm)
            self.mask[li]  = active
            if (active):
                if self.verbose:
                    print 'min active, sm = {}'.format( sm )
                    print vertices
                    print 'slack  is zero for an active constraint'
                slack.value = 0.
                L.minimum[el].slack.value = 0. #slack.value #0.
                G  =  Lambda*( g-(slack**2))
                #f = G*(-1.)  + f # + g#.5*g**2
                f = G + f #new
                self.mask[li]               = True
                self.mask[si]               = True
                L.minimum[el].active        = True
                L.minimum[el].active_slack  = True
            else:
                if self.verbose:
                    print 'min not active, sm = {}'.format( sm )
                #L.value = 0.
                #L.minimum[el].Lagrange.value = 0.
                G  =  Lambda*( g-(slack**2))
                f = G  + f # + g
                self.mask[li]               = False
                self.mask[si]               = False
                L.minimum[el].active        = False
                L.minimum[el].active_slack  = False
                if 'barrier' in self.curve.data:
                    print 'todo: real barriers'
                    if self.curve.data['barrier'] == True:
                        f += L*1./g 
                        
        for el in L.maximum:
            li = L.maximum[el].index
            si = L.maximum[el].slack_index
            pass_value  = L.maximum[el].pass_value
            g           = L.maximum[el].computeUpdate('max', vertices, pass_value)#*(-1.)
            L.maximum[el].computed_value = g
            Lambda      = L.maximum[el].Lagrange
            slack       = L.maximum[el].slack 
            slack.value = -np.sqrt(g.value)
            #slack.value = -np.sqrt(abs(g.value))
            #abs(g.value)  #-g.value
            L.maximum[el].slack.value = slack.value
            sm          = g.value#max(g.value,0.)
            active      = self.is_active(sm)
            self.mask[li]  = active
            if (active):
                if self.verbose:
                    print 'max active, sm = {}'.format( sm )
                    print 'slack  is zero for an active constraint'
                    print vertices
                slack.value = 0. #TLM Feb 2018 (seriously)
                L.maximum[el].slack.value = 0. #slack.value #0.
                G  =  Lambda*( g-(slack**2))#*(-1.) #reversed from floor
                f = G*(-1.)  + f #+ g#.5*g**2
                self.mask[li]               = True
                self.mask[si]               = True
                L.maximum[el].active        = True
                L.maximum[el].active_slack  = True
            else:
                if self.verbose:
                    print 'max not active, sm = {}'.format( sm )
                #L.value = 0.
                #L.maximum[el].Lagrange.value = 0.
                G  =  Lambda*( g-(slack**2))
                f = G*(-1.)  + f  #+ g
                self.mask[li]               = False
                self.mask[si]               = False
                L.maximum[el].active        = False
                L.maximum[el].active_slack  = False
                if 'barrier' in self.curve.data:
                    print 'todo: real barriers'
                    if self.curve.data['barrier'] == True:
                        f += L*1./g
        ##
        ## LagrangeSpline Forces:
        ##
        #scale = f.value
        #scale = 1.0
        """
        if 'force' in self.curve.data:
            for vertex, weightx, weighty in zip(self.curve.data['force']['vertices'],self.curve.data['force']['weights'][0],self.curve.data['force']['weights'][1]):
                
                xattract = scale*weightx*((self.curve.xpts[vertex]-1000.)**2)
                yattract = scale*weighty*((self.curve.ypts[vertex]-1000.)**2)
                print 'force adding ',xattract.value, yattract.value
                f += xattract
                f += yattract
        #"""
        ## LagrangeSurface Curve Forces:
        ##
        """
        if 'force' in self.curve.data:
            if self.curve.data['force']==True:
                f += 15000.*np.sqrt((self.curve.xpts[1]-1000.)**2)#*(-1.)
                f += 5000.*np.sqrt((self.curve.xpts[2]-1000.)**2)#*(-1.)
                f += 5000.*np.sqrt((self.curve.ypts[2])**2)#*(-1.)
        #"""
        
        #if self.fjc:
        #    f = f*L.fj.Lagrange
        
        for el in self.Lagrangian.sp:
            f += self.Lagrangian.sp[el](vertices)
            
        return f
        
        

    def invertReduced(self, f):
        mask = self.mask
        """
           Code for using the Reduced Equations
               Alter AD methods to have end vertices as simple floats to eliminate this step
        """
        #apply a "mask" to remove things like the x and y end vertices from the system of equations.
        FONCf=f.grad[:,mask]
        SOSCf=f.hess[:,mask]
        SOSCf=SOSCf[mask,:]
    
        # set the system of equations for inversion (matrix type requirement - eliminate with np.array instead of np.matrix):
        FONCflist=[]
        for i in range(len(np.transpose(FONCf))):
            FONCflist.append(FONCf[0,i])
        FONCflist=np.asarray(FONCflist)
        
        # Invert the system of equations:
        delta_z = np.zeros(self.n,float)
        try:
            delta_z[self.mask] = np.linalg.solve(SOSCf,FONCflist)
        except:
            print 'Error: could not invert linear system of equations'
            delta_z = 'NAN'
        return delta_z
        
    def mag_lagrangian_multipliers(self, L):
        """
            something like a magnitude goes here
        """
        mag = 0.
        for el in L.equality:
            mag += (L.equality[el].interval_Lagrange**2)#.sqrt()
        for el in L.minimum:
            mag += L.minimum[el].interval_Lagrange
        for el in L.maximum:
            mag += L.maximum[el].interval_Lagrange
        #for F.John Conitions
        #for el in L.obj:
        #    mag += L.obj[el].weight
        if self.fjc:
            mag += L.fj.interval_Lagrange#weight
        return mag
        
    def norm_lagrangian_multipliers(self, L=None):
        if L is None:L = self.Lagrangian
        nmag = self.mag_lagrangian_multipliers(L)#.max.value
        for el in L.equality:
            #L.equality[el].interval_Lagrange = L.equality[el].interval_Lagrange*(nmag**-1.)
            num = L.equality[el].interval_Lagrange/np.sqrt(nmag.value.sup)#*(nmag**-1.)
            L.equality[el].interval_Lagrange.value.inf = num.value.inf
            L.equality[el].interval_Lagrange.value.sup = num.value.sup
        for el in L.minimum:
            #L.minimum[el].interval_Lagrange = L.minimum[el].interval_Lagrange*(nmag**-1.)
            num = L.minimum[el].interval_Lagrange/nmag.value.sup
            L.minimum[el].interval_Lagrange.value.inf = num.value.inf
            L.minimum[el].interval_Lagrange.value.sup = num.value.sup
        for el in L.maximum:
            #L.maximum[el].interval_Lagrange = L.maximum[el].interval_Lagrange*(nmag**-1.)
            num = L.maximum[el].interval_Lagrange/nmag.value.sup
            L.maximum[el].interval_Lagrange.value.inf = num.value.inf
            L.maximum[el].interval_Lagrange.value.sup = num.value.sup
        if self.fjc:
            num = L.fj.interval_Lagrange/nmag.value.sup
            L.fj.interval_Lagrange.value.inf = num.value.inf
            L.fj.interval_Lagrange.value.sup = num.value.sup
        return L
    
    def ReducedConvergenceCriteria(self, f):
        mask = self.mask
        max1 = np.amax(np.transpose(f.grad[:,mask]))
        max2 = np.amin(np.transpose(f.grad[:,mask]))    
        conv = np.sqrt(max1*max1+max2*max2)
        return conv
        
#    def set_point_vertices(self, pt=None):
#        self.interval_vertices = [self.curve.interval_xpts,
#                                  self.curve.interval_ypts]
#        if isinstance(pt,float):
#            pts = []
#            for xyz in 
    
    def prove_simple_zero(self):
        which_verts = None
        zero_exists = True
        vi          = self.get_live_vertex_indicies()
        which_verts = {'sol exists':[],
                       'sol indeterminate':[]}
        try:
            for i,el in zip(vi, self.xnx):
                if el.prove_zero:
                    which_verts['sol exists'].append(i)
                elif not el.prove_zero:
                    which_verts['sol indeterminate'].append(i)
                    zero_exists = False
        except:
            self.interval_optimize(stop = 1)
            zero_exists, which_verts = self.prove_simple_zero()
        return zero_exists, which_verts
    
    def prove_no_solution(self):
        which_verts = None
        no_solution = False
        vi          = self.get_live_vertex_indicies()
        which_verts = {'nonexistence':[],
                       'sol indeterminate':[]}
        try:
            for i,el in zip(vi, self.xnx):
                if el.prove_nosol:
                    which_verts['nonexistence'].append(i)
                    no_solution = True
                elif not el.prove_nosol:
                    which_verts['sol indeterminate'].append(i)
                    
        except:
            self.interval_optimize(stop = 1)
            no_solution, which_verts = self.prove_no_solution()
        return no_solution, which_verts
    
    def prove_stuck(self):
        which_verts = None
        stuck       = True
        vi          = self.get_live_vertex_indicies()
        which_verts = {'changed':[],
                       'unchanged':[]}
        try:
            for i,el in zip(vi, self.xnx):
                if el.prove_stuck:
                    which_verts['unchanged'].append(i)
                elif not el.prove_stuck:
                    which_verts['changed'].append(i)
                    stuck = False
        except:
            self.interval_optimize(stop = 1)
            stuck, which_verts = self.prove_stuck()
        return stuck, which_verts
    
    def test_hessian_consistency(self, hess):
        it = range(len(hess))
        for i in it:
            for j in it:
                bigI = np.argmax(hess[i,j])
                if bigI != 1:
                    if ( hess[i,j].inf==hess[i,j].sup ):
                        pass#print 'ok'
                    else:
                        print i,j, 'max = ', np.argmax(hess[i,j]), hess[i,j]
        return
        
    def test_gradient_consistency(self, grad):
        it = range(len(grad))
        for i in it:
            bigI = np.argmax(grad[i])
            if bigI !=2:
                if grad[i,0] ==grad[i,1]==grad[i,2]:
                    pass #print 'ok'
                else:
                    print i, 'max = ',bigI, grad[i]
        return
        
    def check_vertices_width(self, vertices = None):
        if  vertices is None:
            vertices = [self.curve.interval_xpts,
                        self.curve.interval_ypts]
                                
        
        return
    ## TODO:
    #def break_plot(self):
    #    plt.close()
    #    return
    #    
    #@break_plot
    def get_interval_xpts(self):
        return self.curve.interval_xpts
    def get_interval_ypts(self):
        return self.curve.interval_ypts
    
    def display(self, ideal = 'scalar', mytitle = 'default title', canvas = None, 
                        x_axis_name = 'y', y_axis_name = 'z', closeit_=True):
        if canvas is not None:
            print 'error, ADILS display needs a canvas'
            fig = canvas
        if ideal != 'scalar': return self.display_interval()
        x_axis_name_ = x_axis_name
        y_axis_name_ = y_axis_name
        self.nice_plot = Plotter(title = mytitle, 
                                 x_axis_name = x_axis_name_,
                                 y_axis_name = y_axis_name_)
        xmin,xmax,ymin,ymax = self.extreme_C0()
        dx_ = xmax - xmin
        dy_ = ymax - ymin
        self.nice_plot.set_range(xi = xmin, xe = xmax, dx = dx_,
                                 yi = ymin, ye = ymax, dy = dy_)
        self.nice_plot.set_ticks()
        self.nice_plot.plot_this(plotfunc = self.curve.plot_subplot_curve )
        for el, hue in zip(self.Lagrangian.equality, self.nice_plot.tableau20):
#            if self.Lagrangian.equality[el].type == 'curvature' and 'ini' in mytitle:
#                self.nice_plot.plot_this(
#                    plotfunc = self.Lagrangian.equality[el].object.plot,
#                    object_ = self.Lagrangian.equality[el],
#                    Lspline_ = self,
#                    const_num = el,
#                    docurvature=False)#, 
#                    #color_='green')
#            else:
            self.nice_plot.plot_this(
                plotfunc = self.Lagrangian.equality[el].object.plot,
                object_ = self.Lagrangian.equality[el],
                Lspline_ = self,
                const_num = el)#, 
                #color_='green')
        self.nice_plot.legend_loc[1] -=1.
        self.nice_plot.plot(Lspline = self)
        #self.nice_plot.save_pdf(filename=mytitle, ftype='.png')
        self.nice_plot.save_pdf(filename=mytitle, ftype='.pdf', closeit=closeit_)
        return
    
    def display2(self, ideal = 'scalar', mytitle = 'default title', canvas = None, 
                        x_axis_name = 'y', y_axis_name = 'z'):
        if canvas is not None:
            print 'error, ADILS display needs a canvas'
            fig = canvas
        if ideal != 'scalar': return self.display_interval()
        x_axis_name_ = x_axis_name
        y_axis_name_ = y_axis_name
        self.nice_plot = Plotter(title = mytitle, 
                                 x_axis_name = x_axis_name_,
                                 y_axis_name = y_axis_name_)
        xmin,xmax,ymin,ymax = self.extreme_C0()
        dx_ = xmax - xmin
        dy_ = ymax - ymin
        self.nice_plot.set_range(xi = xmin, xe = xmax, dx = dx_,
                                 yi = ymin, ye = ymax, dy = dy_)
        self.nice_plot.set_ticks()
        self.nice_plot.plot_this(plotfunc = self.curve.plot_subplot_curve )
        for el, hue in zip(self.Lagrangian.equality, self.nice_plot.tableau20):
#            if self.Lagrangian.equality[el].type == 'curvature' and 'ini' in mytitle:
#                self.nice_plot.plot_this(
#                    plotfunc = self.Lagrangian.equality[el].object.plot,
#                    object_ = self.Lagrangian.equality[el],
#                    Lspline_ = self,
#                    const_num = el,
#                    docurvature=False)#, 
#                    #color_='green')
#            else:
            self.nice_plot.plot_this(
                plotfunc = self.Lagrangian.equality[el].object.plot,
                object_ = self.Lagrangian.equality[el],
                Lspline_ = self,
                const_num = el)#, 
                #color_='green')
        self.nice_plot.legend_loc[1] -=1.
        self.nice_plot.plot(Lspline = self)
        #self.nice_plot.save_pdf(filename=mytitle, ftype='.png')
        self.nice_plot.save_pdf(filename=mytitle, ftype='.pdf')
        return
    
    def display_interval(self):
        print 'interval publication plotting function is TBD'
        return
        
    def plotcurve(self):
        self.curve.plot()
        return
    def plot_constraints(self, canvas = None):
        if canvas is None:
            title = 'default constraint plot'
            fig = plt.figure(1)
            ax = fig.add_subplot(111)
            fig.add_subplot(ax)
        else:
            ax = canvas
            
        for el in [1,2]:#self.Lagrangian.equality:
            constraint = self.Lagrangian.equality[el].computeUpdate
            axi = constraint.plot(canvas = ax, Lspline = self)
            fig.add_subplot(axi)
            
        if canvas is None:
            plt.title(title)
            plt.show()
            return 
        else:
            return ax
        
    def plot(self):
        self.plotILspline()
        return
    def plotILspline(self, 
                     title              ='Lspline', 
                     annotations        = None, 
                     orderedtextlist    = None, 
                     plot_cv_interval   = True,
                     canvas             = None, 
                     interior_color     = None):
        
        if canvas is None:
            fig = plt.figure(1)
            ax = SubplotZero(fig, 111)
            fig.add_subplot(ax)
            xmax = self.curve.interval_xpts[-1].value.sup
            xmin = self.curve.interval_xpts[0].value.inf
            ymax = self.curve.interval_ypts[-1].value.sup
            ymin = self.curve.interval_ypts[0].value.inf
            plt.xticks([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.])
            plt.yticks([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.])
            plt.xticks(np.linspace(xmin,xmax,10.))
            plt.yticks(np.linspace(ymin,ymax,10))
            
            plt.grid()
            for direction in ["xzero", "yzero"]:
                ax.axis[direction].set_axisline_style("-|>")
                ax.axis[direction].set_visible(True)
            for direction in ["left", "right", "bottom", "top"]:
                ax.axis[direction].set_visible(False)
        else:
            #fig = figure
            ax  = canvas
            #fig.add_subplot(ax)
            
        if interior_color is None:
            interior_color = 'green'
        else:
            pass
            
        Vmin,Vmax   = self.getIntervalVertices()
        k           = self.curve.k
        nump        = 30
        minCurve    = spline.Bspline(Vmin, k, nump)
        maxCurve    = spline.Bspline(Vmax, k, nump)
        start       = minCurve.vertices[0]
        end         = maxCurve.vertices[-1]
        
        
        ax.text(0, 1.05, 'y', transform=BlendedGenericTransform(ax.transData, ax.transAxes), ha='center')
        ax.text(1.05, 0, 'x', transform=BlendedGenericTransform(ax.transAxes, ax.transData), va='center')
    

    
        ax.plot(minCurve.r[:,0],minCurve.r[:,1], linewidth=1.2, color="black")
        ax.plot(maxCurve.r[:,0],maxCurve.r[:,1], linewidth=1.2, color="black")
        #ax.fill_between(minCurve.r[:,0], minCurve.r[:,1], maxCurve.r[:,1], facecolor='blue', alpha=0.5)
        #ax.fill_between(minCurve.r[:,0], maxCurve.r[:,1], minCurve.r[:,1], facecolor='blue', alpha=0.5)
        if annotations is not None:
            for label in annotations:
                ax.plot(label[1],label[2],label[3])
                plt.annotate(label[0],xy=(label[1],label[2]),xytext=(10,-10),
                                 textcoords='offset points')
        if orderedtextlist is not None:
            for otxt in orderedtextlist:
                ax.text(otxt[0],otxt[1],otxt[2])
                

        if plot_cv_interval == True:
            ax = self.plot_interval_control_vertices(canvas = ax, this_color = interior_color)
        
        ax.plot(minCurve.vertices[:,0], minCurve.vertices[:,1], marker="x",linestyle = "--", color = "red")
        ax.plot(maxCurve.vertices[:,0], maxCurve.vertices[:,1], marker="x",linestyle = "--", color = "blue")
        ax.axis('equal')
        plt.xlim(start[0]-.1,end[0]+.1)
        plt.ylim(end[1]-.1,start[1]+.1)
        #Rotated_Plot = ndimage.rotate(plt, 90)
        #plt.imshow(Rotated_Plot, cmap=plt.cm.gray)        
        
        
        if canvas is None:
            plt.title(title)
            plt.show()
            return 
        else:
            return ax
    
    def getIntervalVertices(self):
        """Return two arrays of the interval curve min and max control points
        To be used in plotting only
        """
        Vmin = []
        Vmax = []
        for x,y in zip(self.curve.interval_xpts, self.curve.interval_ypts):
            Vmin.append([x.value.inf,y.value.inf])
            Vmax.append([x.value.sup,y.value.sup])
        Vmin = np.asarray(Vmin)
        Vmax = np.asarray(Vmax)
        return  Vmin,Vmax
            
    def plot_interval_control_vertices(self, title=None, annotations = None, orderedtextlist = None, canvas = None, this_color = None):
        """
            A somewhat opaque function to plot numerous boxes in matplotlib.
            see:
                http://matplotlib.org/users/path_tutorial.html
            for the write up.
        """
        if canvas is None:
            fig = plt.figure(1)
            ax = fig.add_subplot(111)
        else:
            ax = canvas
            
        if this_color is None:
            this_color = 'green'
        else:
            pass
        
        nrects = self.curve.n
        nverts = nrects*(1+3+1)
        verts = np.zeros((nverts, 2))
        codes = np.ones(nverts, int) * path.Path.LINETO
        
        ixpts = self.curve.interval_xpts
        iypts = self.curve.interval_ypts
        
        left    = []
        right   = []
        top     = []
        bottom  = []
        for ix,iy in zip(ixpts,iypts):
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
        patch = patches.PathPatch(barpath, facecolor=this_color, edgecolor='black', alpha=.5)
        ax.add_patch(patch)

        if canvas is None:
            ax.set_xlim(left[0], right[-1])
            ax.set_ylim(bottom.min(), top.max())
            plt.show()
            return 
        else:
            return ax
            
    def plot_interval_by_index(self, index, title=None, annotations = None, orderedtextlist = None, canvas = None, this_color = None):
        """
            A somewhat opaque function to plot numerous boxes in matplotlib.
            see:
                http://matplotlib.org/users/path_tutorial.html
            for the write up.
        """
        if canvas is None:
            fig = plt.figure(1)
            ax = fig.add_subplot(111)
        else:
            ax = canvas
            
        if this_color is None:
            this_color = 'green'
        else:
            pass
        
        nrects = self.curve.n
        nverts = nrects*(1+3+1)
        verts = np.zeros((nverts, 2))
        codes = np.ones(nverts, int) * path.Path.LINETO
        
        ixpts = self.curve.interval_xpts
        iypts = self.curve.interval_ypts
        
        left    = []
        right   = []
        top     = []
        bottom  = []
        #for ix,iy in zip(ixpts,iypts):
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
        patch = patches.PathPatch(barpath, facecolor=this_color, edgecolor='black', alpha=.5)
        ax.add_patch(patch)

        if canvas is None:
            ax.set_xlim(0., 20.)
            #ax.set_ylim(bottom.min(), top.max())
            ax.axis('equal')
            plt.show()
            return 
        else:
            return ax
            
    def plot_all_intervals(self):
        for i in range(self.curve.n):
            self.plot_interval_by_index(i)
        return

    def print_lagrange(self):
        for el in self.Lagrangian.equality:
            print self.Lagrangian.equality[el].Lagrange
            print self.Lagrangian.equality[el].interval_Lagrange
            print 'Lagrange in Interva_Lagrange?'
            print self.Lagrangian.equality[el].Lagrange.value in self.Lagrangian.equality[el].interval_Lagrange.value
        return
    def translate(self, *args):
        if len(args)==2:
            translate = [args[0],args[1]]
        else:
            translate = args[0]
        v1 = copy.deepcopy(self.curve.vertices)
        v1 = v1 + translate
        self.curve.vertices = v1
        return
        
    def rotate(self, theta):
        """cc-wise 2D rotation
        about the origin
        """
        theta = np.deg2rad(theta)
        s = np.sin(theta)
        c = np.cos(theta)
        v = copy.deepcopy(self.curve.vertices)
        rotmat = np.asarray([[c,-s],
                             [s,c]])
        rv = []
        for el in v:
            rv.append(np.dot(rotmat,el.T))
        rv= np.asarray(rv)
        self.curve.vertices = rv                                                
        return
    
    def flip(self, axis):
        v = copy.deepcopy(self.curve.vertices)
        for i in range(self.curve.n):
            v[i,axis] = -v[i,axis]
        self.curve.vertices = v
        return    

###############################################################################
#
#******************************************************************************
#
def test(*args):
    for el in args:
        for pt in el:
            print pt
    return
    

    
#
#******************************************************************************
#
if __name__ == '__main__':
    #from utility_optimization import package_norms
    option = 'ini_test'
    
    if option =='ini_test':
        interval    = 1
        scalar      = 0
        adaptive    = 0
        ws=10.
        w1=1.0
        w2=1.#0.00001
        w3=0.001
        objWt = {'ws':ws,
                 'w1':w1,
                 'w2':w2,
                 'w3':w3,
                 'Inumber':15}
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
        curve_area = 77.#72.#ia(68.,76.)#84.#
        """ For the ia area, the curve ends up
        with ad(ia[65.839695389,78.160508845]) 5 iterations total
        an enclosure of the constraint
        ad(ia[65.8438631525,78.156338144]) 7 iterations total
        """
        slope = 'down'
        
        ae = alphae
        ab = alphab
        
        ini_v = InitializeControlVertices(alphae=ae, alphab=ab,Cab_given=0.,Cae_given=0.)
        #ini_v.vertices = ini_v.vertices+1
        curve = spline.Bspline(ini_v.vertices, k, nump)  
        #curve.verbose = True
            
        wi = 3.5 #2.5#
        w  = .5 #.65
        ep = 1.e-10
        sp = 1.e-4
        #enforce functional 'good sense' at the interval level:
        x = ini_v.vertices[:,0]
        y = ini_v.vertices[:,1]

        interval_data = {}
        # from the paver:
        """
        interval_data['xmin'] = [ 0.0,   0.01, 0.1, 0.2, 0.3, 0.4, 11.999999999 ] 
        interval_data['xmax'] = [ 0.0000000001,   11.3, 11.4, 11.5, 11.9, 11.99, 12.0000000001 ]
        interval_data['ymin'] = [11.9999999999,   0.01, 0.01, 0.0, 0.,0.,0.]# -0.000000001, -0.000000001, -0.000000001  ] #
        interval_data['ymax'] = [12.0000000001,  12.0, 12., 12., 12., 11., 0.00000000001]
        #"""
        # Max March 2015 ultimate test -Mod used August 2015
        #"""
        interval_data['xmin'] = [ 0.0,   0.01, 0.1, 0.2, 0.000000001, 0.000000001, 11.999999999 ] 
        interval_data['xmax'] = [ 0.0000000001,   12., 12., 12., 12., 12., 12.0000000001 ]
        interval_data['ymin'] = [11.9999999999,   0.01, 0.01, 0.0, 0.,0.,0.]# -0.000000001, -0.000000001, -0.000000001  ] #
        interval_data['ymax'] = [12.0000000001,  12.0, 12., 12., 12., 11., 0.00000000001]
        #"""
        # conistent Max March 2015 ultimate test -Mod used August 2015
        """
        interval_data['xmin'] = [ -0.0000000001,   0.0000000001,0.0000000001,0.0000000001, 0.0000000001,0.0000000001, 11.999999999 ] 
        interval_data['xmax'] = [ .0000000001,   12.0000000001,12.0000000001,12.0000000001,12.0000000001,12.0000000001, 12.0000000001 ]
        interval_data['ymin'] = [11.9999999999,   0.0000000001,0.0000000001,0.0000000001,-0.0000000001,-0.0000000001,-0.0000000001 ] 
        interval_data['ymax'] = [12.0000000001,  12.0000000001,12.0000000001,12.0000000001,12.0000000001,12.0000000001, .00000000001]
        #"""
        """ #TBD, this way:
        interval_data = interval_bounds(curve)
        """
        #consistent and even August 17 2015
        """
        interval_data['xmin'] = [x[0]-ep, x[1]-w,  x[2]-wi, x[3]-wi, x[4]-wi, x[5]-wi, x[6]-ep] 
        interval_data['xmax'] = [x[0]+ep, x[1]+wi, x[2]+wi, x[3]+wi, x[4]+wi, x[5]+w, x[6]]
        interval_data['ymin'] = [y[0]-ep, y[1]-ep, y[2]-wi, y[3]-wi, y[6]-ep,    y[6]-ep,    y[6]-ep] 
        interval_data['ymax'] = [y[0]+ep, y[0]+ep,    y[0]+ep,    y[3]+wi, y[4]+wi, y[5]+ep, y[6]+ep]
        #"""
        #October 2015 -...
        #consistent and not cheating in y  ---  Standard one it seems, in March 2017
        """
        interval_data['xmin'] = [x[0]-ep, x[1]-w,  x[2]-wi, x[3]-wi, x[4]-wi, x[5]-wi, x[6]-ep] 
        interval_data['xmax'] = [x[0]+ep, x[1]+wi, x[2]+wi, x[3]+wi, x[4]+wi, x[5]+.52, x[6]]
        interval_data['ymin'] = [y[0]-ep, y[1]-w, y[2]-wi, y[3]-wi, y[6]-ep,    y[6]-ep,    y[6]-ep] 
        interval_data['ymax'] = [y[0]+ep, y[0]+ep,    y[0]+ep,    y[3]+wi, y[4]+wi, y[5]+w, y[6]+ep]
        #"""
        """
        interval_data['xmin'] = [x[0]-ep, x[1]-w,  x[2]-w, x[3]-wi, x[4]-wi, x[5]-wi,  x[6]-ep] 
        interval_data['xmax'] = [x[0]+ep, x[1]+wi, x[2]+wi, x[3]+wi, x[4]+w, x[5]+w,   x[6]]
        interval_data['ymin'] = [y[0]-ep, y[1]-wi, y[2]-wi, y[3]-wi, y[6]-w, y[6]-w,  y[6]-ep] 
        interval_data['ymax'] = [y[0]+ep, y[0]-ep, y[0]-ep, y[3]+wi, y[4]+wi, y[5]+wi,  y[6]+ep]
        #"""
        
        #wide and even March 2015 sucess
        """
        interval_data['xmin'] = [x[0],    x[1]-w-w,  x[2]-wi, x[3]-wi, x[4]-wi, x[5]-wi, x[6]-ep] 
        interval_data['xmax'] = [x[0]+ep, x[1]+wi, x[2]+wi, x[3]+wi, x[4]+wi, x[5]+w, x[6]]
        interval_data['ymin'] = [y[0]-ep, y[1]-ep, y[2]-wi, y[3]-wi, y[4],    y[5],    y[6]] 
        interval_data['ymax'] = [y[0],    y[1],    y[2],    y[3]+wi, y[4]+wi, y[5]+ep, y[6]+ep]
        #"""
        #modified wide and even March 2015 sucess
        """
        interval_data['xmin'] = [x[0],    .3,      x[2]-wi, x[3]-wi, x[4]-wi, x[5]-wi, x[6]-ep] 
        interval_data['xmax'] = [x[0]+ep, x[1]+wi, x[2]+wi, x[3]+wi, x[4]+wi, 11.8,    x[6]]
        interval_data['ymin'] = [y[0]-ep, y[1]-ep, y[2]-wi, y[3]-wi, 0.,  0.1,     y[6]] 
        interval_data['ymax'] = [y[0],    y[0]-sp, 12.,  y[3]+wi, y[4]+wi, y[5]+ep, y[6]+ep]
        #"""
        
        
        
        #working non 0 starting angles (use box consistency extensively):
        """ # also works with flat start and end
        interval_data['xmin'] = [x[0],    x[1]-w,  x[2]-w,  x[3]-wi, x[4]-wi, x[5]-wi, x[6]-ep] 
        interval_data['xmax'] = [x[0]+ep, x[1]+wi, x[2]+wi, x[3]+wi, x[4]+w, x[5]+w,  x[6]] #changed x[4] +wi to x[4] + w 
        interval_data['ymin'] = [y[0]-ep, 8.,       6.,     y[3]-wi,   0.,      0.,    y[6]] 
        interval_data['ymax'] = [y[0],    12.,      12.,    y[3]+wi,   6.,      4.,    y[6]+ep]
        #"""
        """ #fails with flat start and end
        interval_data['xmin'] = [x[0],    x[1]-w,  x[2]-w,  x[3]-wi, x[4]-wi, x[5]-wi,  x[6]-ep] 
        interval_data['xmax'] = [x[0]+ep, x[1]+wi, x[2]+wi, x[3]+wi, x[4]+w,  x[5]+w,   x[6]] #changed x[4] +wi to x[4] + w 
        interval_data['ymin'] = [y[0]-ep, 8.,       6.,     y[3]-wi, y[6]+sp, y[6]+sp,  y[6]] 
        interval_data['ymax'] = [y[0],    y[0]-sp, y[0]-sp, y[3]+wi,   6.,      4.,     y[6]+ep]
        #"""
        """ #fails with flat start and end
        interval_data['xmin'] = [x[0],        sp,      sp,       x[3]-wi,     6.,      8.,    x[6]-ep] 
        interval_data['xmax'] = [x[0]+ep,     4.,       6.,      x[3]+wi,    12.-sp,  12.-sp,    x[6]]
        interval_data['ymin'] = [y[0]-ep,     8.,       6.,      y[3]-wi,   y[6]+sp, y[6]+sp,    y[6]] 
        interval_data['ymax'] = [y[0],    y[0]-sp,   y[0]-sp,    y[3]+wi,     6.,      4.,    y[6]+ep]
        #"""
        """ #almost succeeds with flat start and end
        interval_data['xmin'] = [x[0],        sp,      sp,       x[3]-wi,     6.,      8.,     x[6]-ep] 
        interval_data['xmax'] = [x[0]+ep,     4.,       6.,      x[3]+wi,    12.-sp,  12.-sp,  x[6]]
        interval_data['ymin'] = [y[0]-ep,     8.,       6.,      y[3]-wi,    y[6],    y[6],    y[6]] 
        interval_data['ymax'] = [y[0],       y[0],   y[0],    y[3]+wi,       6.,      4.,      y[6]+ep]
        #"""
        
        #  ultimate test
        """
        interval_data['xmin'] = [ 0.0,   0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 11.999999999 ] 
        interval_data['xmax'] = [ 0.0000000001,   12., 12., 12., 12., 12., 12. ]
        interval_data['ymin'] = [11.9999999999,   0.0, 0.0, 0.0, 0.0, 0.0, 0.  ] 
        interval_data['ymax'] = [12.0,  12.0, 12., 12., 12., 12., 0.00000000001]
        #"""
        
       
        
        #wide lagrange start
        """
        #interval_data['lmin'] = [-500.,-500.,-500.,-500.,-500.]
        interval_data['lmin'] = [0.e-15,0.e-15,0.e-15,0.e-15,-5000.,-5000.,-5000.,-5000.,-5000.,-5000.,-5000.]
        #interval_data['lmin'] = [0.e-15,0.e-15,0.e-15,0.e-15,0.e-15]
        #interval_data['lmax'] = [500.,500.,500.,500.,500.]
        interval_data['lmax'] = [50000.,5000.,5000.,5000.,5000.,5000.,5000.,5000.,5000.,5000.,5000.,5000.]
        #"""
        
        """
        #interval_data['lmin'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ]
        interval_data['lmin'] = [-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.]
        interval_data['lmax'] = [ 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]
        #"""
        #"""
        interval_data['lmin'] = [-500.,-100000.,-100000.,-200000.,-200000.,-200000.]
        #interval_data['lmin'] = [0.e-15,0.e-15,0.e-15,0.e-15,0.e-15]
        interval_data['lmax'] = [500.,100000.,100000.,200000.,200000.,200000.]
        #FormParam = Formparmeter(curve, kind='equality', location=0., value = 0.)
        #FormParam.add_AngleConstraint()
        #"""
        
        #interval_data['fjcmin'] = -600.
        #interval_data['fjcmax'] = 500.
        interval_data['fjcmin'] = 0.
        interval_data['fjcmax'] = 1.
     
        
        #curve.data['barrier']=False
        
        #---------------------------------------    
        #FormParam.computeUpdate('equality', curve.xpts, curve.ypts)
        
        #t1 = FormParameter(kind='equality', location=0., value = 0.)
        
        large = 12.
        small = 0.
        
        FPD = FormParameterDict(curve) 
        
        #"""
        #FPD.add_AreaConstraint(kind='min', value = 60., weight = 1.)#value = AF.fuzzyNumber(60.,70.,80.))
        #FPD.add_AreaConstraint(kind='max', value = 85., weight = 1.)        
        #FPD.add_AreaConstraint(kind='LS',  value = 72., weight = .5) 
        
        #FPD.add_AreaConstraint(kind='min', value = 98., weight = 1.)#value = AF.fuzzyNumber(60.,70.,80.))
        #FPD.add_AreaConstraint(kind='max', value = 55., weight = 1.) 
        #"""
        FPD.add_AreaConstraint(kind='equality', value = curve_area)#value = AF.fuzzyNumber(60.,70.,80.) )#
        #"""
        #FPD.add_XcConstraint(kind='LS', value =12. , weight = 2000.) #3.35 12.

        #FPD.add_XcConstraint(kind='equality', value = 3.375, weight = 1. )
        #FPD.add_XcConstraint(kind='equality', value = 3.1, weight = 10. )#3.35 12.
        
        #FPD.add_YcConstraint(kind='equality', value = 4.85)
        
        #FPD.add_xPointConstraint(kind = 'equality', location = .5, value = 6.5)
        #FPD.add_yPointConstraint(kind = 'equality', location = .5, value = 5.5)
        
        #FPD.add_xVertexConstraint(kind = 'max', index = 3, value = 5.7, weight = 1000.)
        #FPD.add_yVertexConstraint(kind = 'min', index = 3, value = 5., weight = 1000.)
        #"""
        FPD.add_AngleConstraint(kind='equality', location = 0., value = alphab)#AF.fuzzyNumber(-5.,-2.,0.))#
        FPD.add_AngleConstraint(kind='equality', location = 1., value = alphae)
        FPD.add_CurvatureConstraint(kind='equality', location = 0., value = Cab_given)
        FPD.add_CurvatureConstraint(kind='equality', location = 1., value = Cae_given)
        #"""
        FPD.add_E1(kind='LS', weight = 1.)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
        
        
        #FPD.add_ArcLength(kind='LS', weight = 1.)
        #FPD.add_XcConstraint(kind='equality', value=3.7)
        

        

         
        if scalar:
            L = Lagrangian(FPD, fritzjohn_conditions=False)
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data, fritz_john_conditions = False)
            Lspline.curve.verbose = True
            self = Lspline
            vertices = [self.curve.xpts, self.curve.ypts]
            #self = Lspline
            #vertices = Lspline.optimize(vertices, stop = 35, Lagrangian = self.Lagrangian)
            Lspline.optimize(vertices, stop = 35, Lagrangian = self.Lagrangian)
            print Lspline.curve.area
            self.f      = self.compute_lagrangian(vertices, L)
            a,b = np.linalg.eig(self.f.hess)
            #vertices = Lspline.optimize(vertices)
            Lspline.optimize(vertices)
            #self.adaptive_optimiztion(vertices)
            

 
        if interval:
            L = Lagrangian(FPD, fritzjohn_conditions=False)
            Lspline = IntervalLagrangeSpline(curve, 
                                             L, 
                                             data = interval_data,
                                             fritz_john_conditions = False)
            Lspline.check_IA_correctness = False
            self = Lspline
            vertices = [self.curve.interval_xpts, self.curve.interval_ypts]
            
            #if False:
            go = False
            gok = False
            Test = False
            if go:
                self.contract_constraints([self.curve.interval_xpts,
                                           self.curve.interval_ypts], 
                                           self.Lagrangian)
            if gok:
                self.interval_optimize([self.curve.interval_xpts,self.curve.interval_ypts],solver='K',stop=2)
                self.interval_optimize([self.curve.interval_xpts,self.curve.interval_ypts],solver='both',stop=2)
            if go:
                self.compute_box_consistency(which=[1],xyz=[1],verts=[1],max_it=10)
                self.compute_box_consistency(which=[2],xyz=[1],verts=[5],max_it=10)
                self.compute_box_consistency(which=[3],xyz=[1],verts=[1,2],max_it=10)
                self.compute_box_consistency(which=[4],xyz=[1],verts=[4,5],max_it=10)
                #self.interval_optimize([self.curve.interval_xpts,self.curve.interval_ypts],15)
                self.compute_box_consistency(which=[1],xyz=[0],verts=[1],max_it=10)
                self.compute_box_consistency(which=[2],xyz=[0],verts=[5],max_it=10)
                self.compute_box_consistency(which=[3],xyz=[0],verts=[1,2],max_it=10)
                self.compute_box_consistency(which=[4],xyz=[0],verts=[4,5],max_it=10)
            if go:
                self.interval_optimize([self.curve.interval_xpts,self.curve.interval_ypts],3)
                    
            if False:
                self.interval_optimize([self.curve.interval_xpts,self.curve.interval_ypts],solver='K',stop=5)
                self.interval_optimize([self.curve.interval_xpts,self.curve.interval_ypts],solver='both',stop=2)
                #self.curve.interval_xpts = vertices[0]
                #self.curve.interval_ypts = vertices[1]
                #self.interval_optimize(vertices, 3)
            if Test:
                self.compute_box_consistency(which=[0])#, box_type = 'm')
                #vertices = [self.curve.interval_xpts,
                #            self.curve.interval_ypts]
                self.compute_box_consistency(which=[1],xyz=[1],verts=[1],max_it=10)
                self.compute_box_consistency(which=[2],xyz=[1],verts=[5],max_it=10)
                self.compute_box_consistency(which=[3],xyz=[1],verts=[1,2],max_it=10)
                self.compute_box_consistency(which=[4],xyz=[1],verts=[4,5],max_it=10)
                #self.interval_optimize([self.curve.interval_xpts,self.curve.interval_ypts],15)
                self.compute_box_consistency(which=[1],xyz=[0],verts=[1],max_it=10)
                self.compute_box_consistency(which=[2],xyz=[0],verts=[5],max_it=10)
                self.compute_box_consistency(which=[3],xyz=[0],verts=[1,2],max_it=10)
                self.compute_box_consistency(which=[4],xyz=[0],verts=[4,5],max_it=10)
                #            
                self.contract_constraints([self.curve.interval_xpts,
                                               self.curve.interval_ypts], 
                                               self.Lagrangian)
                self.interval_optimize(stop=5)
            
            if False:
                #self.interval_optimize([self.curve.interval_xpts,self.curve.interval_ypts],solver='K',stop=5)
                #self.interval_optimize([self.curve.interval_xpts,self.curve.interval_ypts],solver='both',stop=2) 
                self.compute_box_consistency(which=[0],box_type = 'm')
                #self.compute_box_consistency(which=[9],box_type = 'm')
                self.compute_box_consistency(which=[1],xyz=[1],verts=[1],max_it=10,box_type = 'm')
                self.compute_box_consistency(which=[2],xyz=[1],verts=[5],max_it=10,box_type = 'm')
                self.compute_box_consistency(which=[3],xyz=[1],verts=[1,2],max_it=10,box_type = 'm')
                self.compute_box_consistency(which=[4],xyz=[1],verts=[4,5],max_it=10,box_type = 'm')
                #self.interval_optimize([self.curve.interval_xpts,self.curve.interval_ypts],15)
                self.compute_box_consistency(which=[1],xyz=[0],verts=[1],max_it=10,box_type = 'm')
                self.compute_box_consistency(which=[2],xyz=[0],verts=[5],max_it=10,box_type = 'm')
                self.compute_box_consistency(which=[3],xyz=[0],verts=[1,2],max_it=10,box_type = 'm')
                self.compute_box_consistency(which=[4],xyz=[0],verts=[4,5],max_it=10,box_type = 'm')
                #self.contract_constraints([self.curve.interval_xpts,
                #                               self.curve.interval_ypts], 
                #                               self.Lagrangian)
                self.interval_optimize(stop=5)
            #self.plotILspline()
                
            #self.contract_constraints()
            #self.interval_optimize(stop=5)

            print ''
            print 'startng box'
            for a, b in zip(interval_data['xmin'], interval_data['xmax']):
                print a,b
            print ''
            print 'ending box:'
            def check_it():
                self.optimize([self.curve.xpts,self.curve.ypts])
                for a,b in zip(self.curve.xpts, self.curve.interval_xpts):
                    print a,b
                    print 'consistent? : ', b.value.contains(a.value)
                for a,b in zip(self.curve.ypts, self.curve.interval_ypts):
                    print a,b
                    print 'consistent? : ', b.value.contains(a.value)
                return
            #check_it()
            
            #"""
            print 'setting variables'
            self.localBasis = L.equality[3].computeUpdate.localBasis
            #self.localBasis = L.equality[4].computeUpdate.localBasis
            V = [self.curve.interval_xpts,
                 self.curve.interval_ypts]
            v = [self.curve.xpts,self.curve.ypts]
            constraint = 0.
            xpts = vertices[0]
            ypts = vertices[1]
            from utilities import vector_AND_
            #vertices = [xpts,ypts]
            
            #self.localBasis = self.Lagrangian.equality[4].computeUpdate.localBasis
            #"""
            point = .5
            L = self.Lagrangian
            
        def set_to_best():
            self.curve.interval_xpts[0].value.inf   = 0.-1.e-14
            self.curve.interval_xpts[0].value.sup   = 0.+1.e-14
            
            self.curve.interval_xpts[1].value.inf   = 1.-1.e-14
            self.curve.interval_xpts[1].value.sup   = 1.+1.e-14
            
            self.curve.interval_xpts[2].value.inf   = 3.-1.e-14
            self.curve.interval_xpts[2].value.sup   = 3.+1.e-14
            
            self.curve.interval_xpts[3].value.inf   = 6.-1.e-14
            self.curve.interval_xpts[3].value.sup   = 6.+1.e-14
            """
            self.curve.interval_xpts[3].value.inf   = 3.
            self.curve.interval_xpts[3].real.value  = 6.
            self.curve.interval_xpts[3].value.sup   = 9.
            """
            
            self.curve.interval_xpts[4].value.inf   = 9.-1.e-14
            self.curve.interval_xpts[4].value.sup   = 9.+1.e-14
            
            self.curve.interval_xpts[5].value.inf   = 11.-1.e-14
            self.curve.interval_xpts[5].value.sup   = 11.+1.e-14
            
            self.curve.interval_xpts[6].value.inf   = 12.-1.e-14
            self.curve.interval_xpts[6].value.sup   = 12.+1.e-14
            #__________________________________________________
            self.curve.interval_ypts[0].value.inf   = 12.-1.e-14
            self.curve.interval_ypts[0].value.sup   = 12.+1.e-14
            
            self.curve.interval_ypts[1].value.inf   = 12.-1.e-14
            self.curve.interval_ypts[1].value.sup   = 12.+1.e-14
            
            self.curve.interval_ypts[2].value.inf   = 12.-1.e-14
            self.curve.interval_ypts[2].value.sup   = 12.+1.e-14
            
            self.curve.interval_ypts[3].value.inf   = 6.-1.e-14
            self.curve.interval_ypts[3].value.sup   = 6.+1.e-14
            """
            self.curve.interval_ypts[3].value.inf   = 1.
            self.curve.interval_ypts[3].value.sup   = 11.
            """
            
            self.curve.interval_ypts[4].value.inf   = 0.-1.e-14
            self.curve.interval_ypts[4].value.sup   = 0.+1.e-14
            
            self.curve.interval_ypts[5].value.inf   = 0.-1.e-14
            self.curve.interval_ypts[5].value.sup   = 0.+1.e-14
            
            self.curve.interval_ypts[6].value.inf   = 0.-1.e-14
            self.curve.interval_ypts[6].value.sup   = 0.+1.e-14
            return
 