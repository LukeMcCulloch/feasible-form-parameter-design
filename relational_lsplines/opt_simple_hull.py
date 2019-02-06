# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 16:16:11 2016

@author: lukemcculloch



notes:
    simple_hull_rules class: HullDesignNode
        -contains the rules for the bare hull
        imported as follows,
        from simple_hull_rules import HullDesignNode as hullclp

    hull_from_simple_designspace class: hullmethods
         -takes a point in the design parameter 
        configuration space
        -and uses Nowacki Harries FPD methods 
        to generate a set of hull curves
        imported as follows,
        import hull_from_simple_designspace as hullmethods 

File:
    purpose:    This file generates a random set of mutually consistent
                form parameters for design of a bare hull-form
    contains:   a robust loop method for narrowing a design space
                of form parameters to a thin design vector of parameters
                to be transformed into a B-spline surface
                
                

Installation:
    cd to relational_lsplines directory with setup.py
    and enter at the prompt:
    
        pip install .
    
    or if updating with new code, enter:
        
        pip install . --update
"""
"""
Little irks:
    *Access for the bare hull and bulbous bow thin design spaces are 
    slightly inconsistent, though it's sort of arcane:
    
        bare_hull thinned design lives in two spots:
            
        *SD.design_space.rgp.env
        *SD.thin_hull_design.hullconstraints
        
        bbow thinned design lives in two spots:
        *SD.bbow.rgp.env
        *SD.thin_bbow_design.rgp.env

Bugs:
"""
import numpy as np
import matplotlib.pyplot as plt
import itertools
#import pickle as pickle
import cPickle as pickle
#from interval_arithmetic import ia
#
import FileTools
#
from extended_interval_arithmetic import ia #needed!
import sqKanren as lp
#
# HullDesignNode contains the rules for the bare hull.
# i.e. Relational Rules (old style):
from simple_hull_rules import HullDesignNode as hullclp #relational rules, no FPD!

#
# Below next,
# 'hullmethods' takes a point in the 
# design parameter configuration space
# and uses Nowacki Harries FPD methods to generate a set of hull curves
# i.e., FDP methods:
import hull_from_simple_designspace as hullmethods #FPD curve and bare hull generation only, no rules!
#
import BulbRules as bbr #new style relational rules
import GeometryGenerator as ggr #FPD
#
#import hull_use_HICLP
#
import sobol
import random 
import copy

from frames import Frame,DeepVector

import thbsurface as thbspline
#ProjectionGetter = thbspline.THButilities.ProjectionGetter


from ADILS import SolverPostProcessor
import utility_optimization as uopt

#give a repl to the user
# language help!
import repl
#
#rouding table output
from decimal import getcontext, Decimal

#tricks that were used to compose... cartesian products
#cartesian product from sets:
#    http://metapython.blogspot.com/2011/10/cartesian-product-for-sets.html
#instancemethod = lambda f: (lambda *a, **kw: f(*a, **kw))
#class mset(set):
#    __mul__ = __rmul__ = instancemethod(itertools.product)
#class strset(str):
#    __mul__ = __rmul__ = instancemethod(itertools.product)
    

#class DesignNode(object):
#    """NOT USED
#    wraps the States Class
#    logical resolution produces children 
#    e.g. with contracted states
#    
#    States      : cargo
#    parent      : parent
#    children    : any number of Design Nodes Derived From this one
#    """
#    def __init__(self,states=None, parent=None):
#        self._states = states 
#        self.parent = parent
#        self.children = []
#        
#    @property
#    def states(self):
#        return self._states
#    @states.setter
#    def states(self, states):
#        self.children += DesignNode(states, parent=self)
#        
#    #    @property
#    #    def children(self):
#    #        return self._children
#    #    @children.setter
#    #    def states(self, states):
#    #        self._children += DesignNode(states, parent=self)
#    #        
#    #    @property
#    #    def parent(self):
#    #        return self._parent
#    #    @parent.setter
#    #    def states(self, parent):
#    #        self._parent = parent
#    
#    
#
#
#
##class DesignTree(hullclp):
##    """This is like the States Class
##    Except contraction produce children with contracted states
##    """
##    def __init__(self, root):
##        self.root = Node(root)
##    
##    #    @property
##    #    def root(self):
##    #        return self._root
##    #    @root.setter
##    #    def root(self, state):
##    #        return
#
#

def close_3():
    """Closes 3 plots at a time - because I am lazy!
        Parameters
        ----------
            none
        Returns
        ----------
            nothing
            
        Notes
        ----------
    """
    plt.close()
    plt.close()
    plt.close()
    return

def get_projection_matrices(ProjectionGetterInstance): 
    """Closes 3 plots at a time - because I am lazy!
        Parameters
        ----------
            none
        Returns
        ----------
            nothing
            
        Notes
        ----------
    """
    #import cPickle as pickle
                    
    file_Name3 = "rm_11x7"
    file_Name4 = "rm_19x11"
    file_Name5 = "rm_35x19"
    file_Name6 = "rm_67x35"
    file_Name7 = "rm_7x5"
    
    ## 7x5
    fileObject = open(file_Name7,'r')  
    rm_7x5 = pickle.load(fileObject) 
    fileObject.close()
    
    ## 11x7
    fileObject = open(file_Name3,'r')  
    rm_11x7 = pickle.load(fileObject) 
    fileObject.close()
    
    ## 19x11
    fileObject = open(file_Name4,'r')  
    rm_19x11 = pickle.load(fileObject) 
    fileObject.close()

    ## 35x19
    fileObject = open(file_Name5,'r')  
    rm_35x19 = pickle.load(fileObject) 
    fileObject.close()
    
    ## 67x35
    fileObject = open(file_Name6,'r')  
    rm_67x35 = pickle.load(fileObject) 
    fileObject.close()
    
    
    ProjectionGetterInstance.rm_7x5 = rm_7x5
    ProjectionGetterInstance.rm_11x7 = rm_11x7
    ProjectionGetterInstance.rm_19x11 = rm_19x11
    ProjectionGetterInstance.rm_35x19 = rm_35x19
    ProjectionGetterInstance.rm_67x35 = rm_67x35
    return ProjectionGetterInstance

class TestParameters(object):
    """
        Waring, do not print designs_
        designs_ is the product space of
        the design variables.  This could be quite large!
        
        DONE: make variable based sizing of sobol sample space 
    """
    def __init__(self, spec, N=10):
        self.N = N
        self.testparam = self.generate_test_parameters(spec)
        self.designs_ = None
        self.keys = None # design parameter value names in Cartesian product space ordering
        self.generate_search_space()
        return
        
    def generate_test_parameters(self, spec):
        tp = {}
        for key in spec.__dict__:
            tp[key] = self.generate_quasirandom_sequence(
                    spec.__dict__[key])
        return tp
    
    def generate_quasirandom_sequence(self,
                                      allow_range):        
        s1 = sobol.sobolSeq([1,1],[1,1])
        inf = allow_range[0]
        sup = allow_range[1]
        r = sup-inf
        x = np.asarray(
                [ next(s1) for _ in range(self.N) ]
            )
        return r*x+inf
        
    def generate_search_space(self):
        """
            Builds the Cartesian Product Space of Design test Parameters
            and
            sets keys : the correct mapping from parmeter name
                    to design parameter value in the design accessor
        """
        lists = []
        keys = []
        for key in self.testparam:
            keys.append(key)
            lists.append(self.testparam[key])
        designs = itertools.product(*lists)
        self.designs_ = designs
        self.keys = keys
        return
    
    def get_design_(self):
        for el in self.designs_:
            yield el
    
    def next_design(self):
        return self.designs_.next()#flatten(self.design_.next())
    
            

class DesignSpecification(object):
    def __init__(self,
                 lwl = None,
                 draft = None,
                 bwl = None,
                 vol = None,
                 Cb = None,
                 Clcg = None,
                 LCG = None,
                 Cmidshp = None,
                 Cwp = None,
                 Ccp = None
                 ):
        self.lwl = lwl
        self.draft = draft
        self.bwl = bwl
        self.vol = vol
        self.Cb = Cb
        self.LCG = LCG
        self.Clcg = Clcg
        self.Cmidshp = Cmidshp
        self.Cwp = Cwp
        self.Ccp = Ccp
        
        
    def print_specification(self):
        print 'lwl ',self.lwl
        print 'draft ',self.draft
        print 'bwl ',self.bwl
        print 'vol ',self.vol
        print 'Cb ',self.Cb
        print 'LCG ',self.LCG
        print 'Clcg ',self.Clcg
        print 'Cmidshp ',self.Cmidshp
        print 'Cwp',self.Cwp
        print 'Ccp ',self.Ccp
        
    
    #def set_ranges_for_details(self):
    #    self.Acp = 

#class DesignParams(object):
#    def __init__(self,
#                 lwl,
#                 draft,
#                 bwl,
#                 vol,
#                 LCG,
#                 Clcg,
#                 Cmidshp,
#                 Cb,
#                 Cwp,
#                 Ccp):
#        self.lwl = lwl
#        self.draft = draft
#        self.bwl = bwl
#        self.vol = vol
#        self.Cb = Cb
#        self.Cwp = Cwp
#        self.LCG = LCG
#        self.Clcg = Clcg
#        self.Cmidshp = Cmidshp
#        self.Ccp = Ccp
#        
#    def print_params(self):
#        print 'lwl ',self.lwl
#        print 'draft ',self.draft
#        print 'bwl ',self.bwl
#        print 'vol ',self.vol
#        print 'Cb ',self.Cb
#        print 'Cwp',self.Cwp
#        print 'LCG ',self.LCG
#        print 'Cmidshp ',self.Cmidshp
#        print 'Ccp ',self.Ccp
    
        

#class DesignSpace(hull_use_HICLP.hullclp):
class DesignSpace(object):
    """DEPRECATED
    -should still work though.
    
    Parameters:
        spec : a Design_Specification object
    """
    SMALL = .1
    def __init__(self, spec, verbose=False):
        self._verbose = verbose
        self.spec = spec
        self.sobol_seq = sobol.sobolSeq([1,1],[1,1])#replaced with random.random for now
        #alternatvely, initialize a true IA hull
        self.hdp = hullclp() #hull_use_HICLP.hullclp()
        self.hdp.lwl = ia(spec.lwl[0],spec.lwl[1])
        self.hdp.draft = ia(spec.draft[0],spec.draft[1])
        self.hdp.bwl = ia(spec.bwl[0],spec.bwl[1])
        self.hdp.vol = ia(spec.vol[0],spec.vol[1])
        self.hdp.Cb = ia(spec.Cb[0],spec.Cb[1])
        self.hdp.Cwp = ia(spec.Cwp[0],spec.Cwp[1])
        self.hdp.Ccp = ia(spec.Ccp[0],spec.Ccp[1])
        self.hdp.LCG = ia(spec.LCG[0],spec.LCG[1])
        self.hdp.Clcg = ia(spec.Clcg[0],spec.Clcg[1])
        self.hdp.Cmidshp = ia(spec.Cmidshp[0],spec.Cmidshp[1])
        self.hdp.AC_revise()
        #"""
        #self.hdp.states = self.hdp.shape_bow_section_inizers()
        #self.hdp.states = self.hdp.shape_stern_section_inizers()
        #self.hdp.states = self.hdp.SAC_run_area_ini()
        #self.hdp.states = self.hdp.SAC_Xc_init()
        self.hdp.ini_areas()
        self.hdp.sac_XC_ini()
        self.hdp.states = self.hdp.LCG_constraints()
        #"""
        self.hdp.AC_revise()
        self.feasible_designs = []
        self.infeasible_designs = []
        #"""
        #    def get_design_parameter_list(self):
        #        
        #        return
        
    def factorial_search(self, test_parameters):#, designs = None):
        """
            test_parameters :: object of class TestParameters
            
            reference bckwds w/ str :: lwl <-> candidate[str(ds.hdp.lwl)]
        """
        
        """
            reference this way!
            print 'lwl = ', candidate[str(ds.hdp.lwl)]
            
            TODO: now copy the design space and reduce
            per this specific design
            then
        """
        SMALL = DesignSpace.SMALL
        tp = test_parameters
        keys = tp.keys
        
        #candidate = tp.next_design()
        for i, el in enumerate(tp.designs_):
            print 'design ',i
            candidate = dict(zip(keys, el)) #flatten(el)))
            design = copy.deepcopy(self.hdp)
            
            design.lwl = ia(candidate[str(design.lwl)]-SMALL,
                             candidate[str(design.lwl)]+SMALL)
            design.bwl = ia(candidate[str(design.bwl)]-SMALL,
                             candidate[str(design.bwl)]+SMALL)
            design.draft = ia(candidate[str(design.draft)]-SMALL,
                             candidate[str(design.draft)]+SMALL)
            design.vol = ia(candidate[str(design.vol)]-SMALL,
                             candidate[str(design.vol)]+SMALL)
            design.Cb = ia(candidate[str(design.Cb)]-SMALL,
                             candidate[str(design.Cb)]+SMALL)
            design.Cwp = ia(candidate[str(design.Cwp)]-SMALL,
                             candidate[str(design.Cwp)]+SMALL)
            design.Ccp = ia(candidate[str(design.Ccp)]-SMALL,
                             candidate[str(design.Ccp)]+SMALL)
            design.LCG = ia(candidate[str(design.LCG)]-SMALL,
                             candidate[str(design.LCG)]+SMALL)
            design.Cmidshp = ia(candidate[str(design.Cmidshp)]-SMALL,
                             candidate[str(design.Cmidshp)]+SMALL)
            if design.valid_design:
                self.feasible_designs.append(design)
            else:
                self.infeasible_designs.append(design)
        #print 'lwl = ', candidate[str(ds.hdp.lwl)] solve for the boat!
        return
        
    def save_space(self, space=None):
        if space is None:
            self.saved_space = copy.deepcopy(self.hdp)
        else:
            self.saved_space = space
        return
        
    def get_space(self):
        return self.saved_space
        
        
    def tree_search(self, 
                     Nd = 1,
                     initial_state = None,
                     sobol_seq=None):
        """
            sobol 'dive' into the design space...
            reference bckwds w/ str :: lwl <-> candidate[str(ds.hdp.lwl)]
        """
        """
            Each time we choose next parameter 
                to pick values from
            we move one level further down the tree 
        """
        """
            starts out as depth first search
        """
        self.saved_space = copy.deepcopy(self.hdp)
        #
        if sobol_seq is None:sobol_seq=self.sobol_seq
        if initial_state is None: initial_state = self.hdp#hullclp( ds.hdp.states[0])
        #
        SMALL = DesignSpace.SMALL
        #
        design = initial_state
        keys = design.keys
        k1 = set(design.Primary)
        k2 = set(design.Coefficients)
        k3 = set(design.Areas)
        #k4 = set(design.BowCurveCoeff) & set(design.SternCurveCoeff)
        ksac = set(design.list_SAC)
        kr = set(keys) - k1 - k2 - k3 - ksac
        
        
        design_sq = lp.States(self.hdp.states)
        
        #self.hdp.states = design_sq.states
        self.design_sq = design_sq
        design_sq = self.search_loop(design_sq,
                                     k2,
                                     sobol_seq)
        design_sq = self.search_loop(design_sq,
                                     k1,
                                     sobol_seq)
        design_sq = self.search_loop(design_sq,
                                     ksac,
                                     sobol_seq)
        design_sq = self.search_loop(design_sq,
                                     k3,
                                     sobol_seq)
        #        design_sq = self.search_loop(design_sq,
        #                                     k4,
        #                                     sobol_seq)
        design_sq = self.search_loop(design_sq,
                                     kr,
                                     sobol_seq)
        self.feasible_designs.append(design_sq.get_latest_child())
        return
        
    def search_loop(self, 
                    design_sq, 
                    searchlist,
                    generator):
        for key in searchlist:
            loop = True
            c = 0
            while loop is True:
                c += 1
                x = generator.next() 
                x = random.random()
                design_sq = self.set_this(key,x,design_sq)
                kid = design_sq.get_latest_child()
                if len(kid(key))==0:
                    s = kid.parent
                    if c>=2:
                        s=kid.parent.parent
                    self.hdp.states = s.states
                    self.infeasible_designs.append(s.children.pop(0))
                    tst = s(key)
                    try:
                        if tst.width() < SMALL and not tst.isempty:
                            loop=False
                    except:
                        loop=True
                else:
                    loop=False
        self.design_sq = design_sq
        return design_sq
            
        
    def set_this(self, key,x,design_sq):
        """TODO,
        if something goes wrong, back up one level
        and breadth first search
        """
        val = self.hdp(key)[0]
        if isinstance(val, ia):
            val = self.hdp(key)[0].getpoint(x)
            val = ia(val,val)
            self.hdp.__setattr__(key.name,val)
            lastnode = design_sq.get_latest_child()
            nextnode = lp.States(self.hdp.states, parent=lastnode )
            lastnode.children.append(nextnode)
            self.hdp.states = nextnode.states
            print 'done ', key,'=>', val
        else:
            print 'skipped ', key
        return design_sq
        
    
        
    def print_widths(self, design=None):
        if design is None:
            design = self.hdp
        for el in design.keys:
            if isinstance(design(el)[0], ia):
                print el, design(el)[0].width()
            else:
                print el, design(el)
        return
        
    def check_widths(self, design=None):
        ok = True
        if design is None:
            design = self.hdp
        for el in design.keys:
            if isinstance(design(el)[0], ia):
                if design(el)[0].width() > 0.:
                    ok=False
                    print el
        return ok

def print_fractions():            
    print '\nPrimary Coefficients'
    self.hdp.print_primary()
    print '\nSAC Coefficients'
    self.hdp.print_SAC_stats()
    print 'hull2 fwd Xc = ', hull2.Lsplines.SAC_entrance.curve.Xc
    print 'hull2 aft Xc = ', hull2.Lsplines.SAC_run.curve.Xc
    print '\n Entrance'
    print hull2.sac_co.Xe / hull2.sac_co.Le
    #print '\n mid'
    #print hull2.sac_co.Xm / hull2.sac_co.Lm
    print '\n Run'
    a = hull2.sac_co.Le + hull2.sac_co.Lm
    print (hull2.sac_co.Xr-a) / hull2.sac_co.Lr
    a= hull2.Lsplines.SAC_entrance.curve.Xc
    print '\nXc fwd = ',a, self.hdp(self.hdp.SAC_fwd_Xc)
    a = hull2.Lsplines.SAC_run.curve.Xc
    print 'Xc aft = ',a,  self.hdp(self.hdp.SAC_run_Xc)
    
    hull2.SAC.compute_area()
    a0 = hull2.SAC.area
    a1 = hull2.Lsplines.SAC_entrance.curve.area
    hull2.Lsplines.SAC_mid.curve.compute_area()
    a2 = hull2.Lsplines.SAC_mid.curve.area
    a3 = hull2.Lsplines.SAC_run.curve.area
    
    
    print 'SAC fwd area = ',a1,self.hdp(self.hdp.SAC_entrance_area)
    print 'SAC ,mid area = ',a2,self.hdp(self.hdp.SAC_mid_area)
    print 'SAC aft area = ',a3,self.hdp(self.hdp.SAC_run_area)
    print 'SAC area = ',a0,self.hdp(self.hdp.vol)
    print 'area from piece = ',a1+a2+a3
    return
    

def print_midsection():
    self.hdp.print_areas()
    self.hdp.print_SAC_stats()
    print 'Amsh = ',self.hdp(self.hdp.Amsh)
    print 'SAC_mid_len = ',self.hdp(self.hdp.SAC_mid_len)
    print 'lfsac = ',self.hdp(self.hdp.lfsac)
    print '\n Amsh*SAC_mid_len :'
    print self.hdp(self.hdp.Amsh)[0]*self.hdp(self.hdp.SAC_mid_len)[0]
    print '\nSAC_mid_area = ',self.hdp(self.hdp.SAC_mid_area)
    
    hull2.Lsplines.SAC_mid.curve.compute_area()
    print 'Lsplines SAC middle area = ',hull2.Lsplines.SAC_mid.curve.area
    print 'should be = ', self.hdp(self.hdp.Amsh)[0]*self.hdp(self.hdp.SAC_mid_len)[0]
    print 'should be = ', self.hdp(self.hdp.Amsh)[0]*self.hdp(self.hdp.lfsac)[0]
    print 'should be consistent with ', self.hdp(self.hdp.SAC_mid_area)
    return
    
    
def give_results(hull1):
    hull2 = copy.deepcopy(hull1)
    hull2.define_hull()
    hull2.definestations()#!!
    hullmethods.compute_hull_and_transverse_curves(hull2,
                                                   sacflat_=True)
    
    hullmethods.loft_hull(hull2)
    return hull2
    
    
def plot_results(hull2):
    hull2.Lsplines.DWL.curve.plot()
    hull2.Lsplines.CProfile.curve.plot()
    hull2.Lsplines.SAC.curve.plot()
    #hull2.Lsplines.SAC_entrance.curve.plot()
    #hull2.Lsplines.SAC_mid.curve.plot()
    #hull2.Lsplines.SAC_run.curve.plot()
    
    hullmethods.hull_gui(hull2)
    
    hull2.DWL.plot3DmultiList(hull2.tcurvenet,hull2.lcurvenet,
                              limx=100.,limy=100.,limz=100.)
    
    
    hull2.DWL.plot3DmultiList(hull2.lcurvenet,hull2.tcurvenet,
                              limx=100.,limy=100.,limz=100.)
    return
    

def interface_to_bbow(hull_states_object):
    """
    Add a bulbous bow to the new 
    randomly generated bare hull

    Parameters
    ----------
        hull_states_object : simple_hull_rules HullDesignNode imported as hullclp
            i.e. this is the old fashion hull rules container class with AC_revise
    Notes
    ----------
        hull_states_object = self.hdp
    """
    bare_hull   = hull_states_object
    ship_beam   = hull_states_object.bwl
    ship_depth  = hull_states_object.draft
    ship_Lpp    = hull_states_object.lwl  #TODO:  kraft wants Lpp!
    ship_Amshp  = hull_states_object.Amsh
    ship_Acp    = hull_states_object.Acp
    ship_Awp    = hull_states_object.Awp
    ship_Vol    = hull_states_object.vol
    #    BulbObject = bbr.BulbGenerator(ship_beam,
    #                                   ship_depth,
    #                                   ship_Lpp,
    #                                   ship_Amshp,
    #                                   ship_Acp,
    #                                   ship_Awp,
    #                                   ship_Vol,
    #                                   bare_hull)
    BulbObject = ggr.GeometryGenerator(ship_beam,
                                   ship_depth,
                                   ship_Lpp,
                                   ship_Amshp,
                                   ship_Acp,
                                   ship_Awp,
                                   ship_Vol,
                                   bare_hull)
    print BulbObject.rgp.env
    
#    BulbObject.A_mid = BulbObject.A_mid == ia(10.,20.) #issue: order of assignment "matters"
#    BulbObject.BulbBeam = BulbObject.BulbBeam == ia(5.,10.)
#    BulbObject.CBulbMid = BulbObject.CBulbMid == ia(.5,1.)
#    BulbObject.CBulbWtrPln = BulbObject.CBulbWtrPln == ia(.5,1.)
#    #BulbDepth = BulbDepth == ia(-10.1,10.)
#
#    
#    #bgp.add_one_rule(BulbDepth,'BulbDepth')
#    BulbObject.rgp.add_one_rule(BulbObject.A_mid,'A_mid')
#    BulbObject.rgp.add_one_rule(BulbObject.BulbBeam,'BulbBeam')
#    BulbObject.rgp.add_one_rule(BulbObject.CBulbMid,'CBulbMid')
#    BulbObject.rgp.add_one_rule(BulbObject.CBulbWtrPln,'CBulbWtrPln')
#    #BulbObject.bgp.add_one_rule(BulbObject.BulbDepth,'BulbDepth') 
    
    BulbObject.rgp.compute_fresh_rules_graph()
    BulbObject.tree_search()
    
    
    
    BulbObject.get_curves()
    
    BulbObject.inib_ideal.solve_curve_ini()
    BulbObject.inib_ideal.solve_curve_fini()
    
    BulbObject.mb_ideal.solve_curve_ini()
    BulbObject.mb_ideal.solve_curve_fini()
    #BulbObject.mb_ideal.Lspline.curve.plotcurve_detailed()
    BulbObject.mb_ideal.Lspline.curve.plotcurve()
    
    BulbObject.bwl_ideal.solve_curve_ini()
    BulbObject.bwl_ideal.solve_curve_fini()
    #BulbObject.bwl_ideal.Lspline.curve.plotcurve_detailed()
    BulbObject.bwl_ideal.Lspline.curve.plotcurve()
    
    
    BulbObject.fsb_ideal.solve_curve_ini()
    BulbObject.fsb_ideal.solve_curve_fini()
    #BulbObject.bwl_ideal.Lspline.curve.plotcurve_detailed()
    BulbObject.fsb_ideal.Lspline.curve.plotcurve()
    
    
    BulbObject.cp_ideal.solve_curve_ini()
    BulbObject.cp_ideal.solve_curve_fini()
    #BulbObject.cp_ideal.Lspline.curve.plotcurve_detailed()
    BulbObject.cp_ideal.Lspline.curve.plotcurve()
    
    BulbObject.specialized_bow_development_nose()
    BulbObject.make_bulb()
    
    print 'Xc = ',BulbObject.mb_ideal.Lspline.curve.Xc
    print 'Xc = ',BulbObject.bwl_ideal.Lspline.curve.Xc
    print 'Xc = ',BulbObject.cp_ideal.Lspline.curve.Xc
    return BulbObject

def interface_to_bbow_make_thin_bbow(hull_states_object,
                                     bbow_interface_section_area,
                                     bbow_length):
    """
    Develop a thin set of form parameters
    for a bulbous bow whihc is compatible with
    a bare hull

    Parameters
    ----------
        hull_states_object : simple_hull_rules HullDesignNode imported as hullclp
            i.e. this is the old fashion hull rules container class with AC_revise
    Notes
    ----------
        hull_states_object = self.hdp
    """
    bare_hull   = hull_states_object
    ship_beam   = hull_states_object.bwl
    ship_depth  = hull_states_object.draft
    ship_Lpp    = hull_states_object.lwl  #TODO:  kraft wants Lpp!
    ship_Amshp  = hull_states_object.Amsh
    ship_Acp    = hull_states_object.Acp
    ship_Awp    = hull_states_object.Awp
    ship_Vol    = hull_states_object.vol
    
    A_mid = bbow_interface_section_area
    
    BulbObject = ggr.GeometryGenerator(ship_beam,
                                   ship_depth,
                                   ship_Lpp,
                                   ship_Amshp,
                                   ship_Acp,
                                   ship_Awp,
                                   ship_Vol,
                                   bare_hull,
                                   given_A_mid=A_mid,
                                   given_length = bbow_length)
    print BulbObject.rgp.env
    
    
    BulbObject.rgp.compute_fresh_rules_graph() 
    
    
    
    
    BulbObject.tree_search() 
    
    return BulbObject

def interface_to_bbow_make_bare_hull(BulbObject, 
                                     skip_search=False):
    """Use FPD to develop a B-spline hull form
    from a thin set of form parameters (from relational rules)

    Parameters
    ----------
        BulbObject : thin set of form parameters describing the bbow
        
    
    Notes
    ----------
        a BulbObject has it's own ability to FPD itself into a hull
        matchin it's description.
    """
    if not skip_search:
        #skip this if making a new hull in simple_hull_rules_language!
        BulbObject.tree_search()  #moved to complex_hull_design_and_generation
    
    
    BulbObject.get_curves()
    
    BulbObject.inib_ideal.solve_curve_ini()
    BulbObject.inib_ideal.solve_curve_fini()
    
    BulbObject.mb_ideal.solve_curve_ini()
    BulbObject.mb_ideal.solve_curve_fini()
    
    BulbObject.bwl_ideal.solve_curve_ini()
    BulbObject.bwl_ideal.solve_curve_fini()
    
    
    BulbObject.fsb_ideal.solve_curve_ini()
    BulbObject.fsb_ideal.solve_curve_fini()
    
    
    BulbObject.cp_ideal.solve_curve_ini()
    BulbObject.cp_ideal.solve_curve_fini()
    
    BulbObject.specialized_bow_development_nose()
    BulbObject.make_bulb()
    
    print 'Xc = ',BulbObject.mb_ideal.Lspline.curve.Xc
    print 'Xc = ',BulbObject.bwl_ideal.Lspline.curve.Xc
    print 'Xc = ',BulbObject.cp_ideal.Lspline.curve.Xc
    return BulbObject


def api(design_spec):
    """
    NOT The API for _Everything_
     - instead this shows why we need a class 
    for all the scripting which 
    currently rides on top of _Everything_
    
    USE ShipDesigner instead of this!
    
    Parameters:
        design_spec = class DesignSpecification
    Notes:
        basic stuff that first goes in
        when starting of a boat design run
        
        TODO: move to some kind of class structure
        see the base class...
    """
    spec = design_spec
    self = DesignSpace(spec)
    self.tree_search()
    print '\n---------------------------------------------------------------\n'
    print 'initial design space:'
    self.print_widths()
    print '\n---------------------------------------------------------------\n'
    hull1 = hullmethods.Hull(States = self.hdp, 
                             which_state=0,
                             revise=False)
    hull2 = give_results(hull1)
    plot_results(hull2)
    bbow = interface_to_bbow(self.hdp)
    bbow.mb_ideal.Lspline.curve.plotcurve()
    bbow.bwl_ideal.Lspline.curve.plotcurve()
    bbow.cp_ideal.Lspline.curve.plotcurve()
    bbow.tcurvenet[0].plot3DmultiList(bbow.tcurvenet,
                  bbow.lcurvenet,
                  limx=10.,limy=10.,limz=10.)
    
    hull2.hullsurf.plotSurface(limx=10.,limy=10.,limz=10.)
    bbow.hullsurf.plotSurface(limx=10.,limy=10.,limz=10.)
    
    
    hull2.hullsurf.plotMultiSurface(
                    [hull2.hullsurf,bbow.hullsurf],
                    limx=100.,limy=100.,limz=100.)
    
    
    b2 = copy.deepcopy(bbow.hullsurf)
    #b2 = b2.rotate_surface( 
    #                (0.,0.,1.), np.pi/2.)
    b2 = b2.reflect_surface(np.asarray([0.,0.,1.]),
                            np.asarray([0.,0.,1.]))
    b2.compute_surface()
    
    b2.plotMultiSurface([bbow.hullsurf,b2],
                        limx=100.,limy=100.,limz=100.)
   
    return


class ShipDesigner(object):
    """
        Parameters
        ----------
        Returns
        ----------
        Developer Notes
        ----------
            this is a base class
            to facilitate the starting of a boat design run
            
            parameters:
                design_specification : DesignSpecification class instance
        
        Organization
        ----------
            1 bare hull design selection (relational rules)
            2 bare hull FPD (form parameter design)
            3 bbow design selection (relational rules)
            4 bbow FPD
        
        
        TODO
        ----------
            -if the user passes into a stage, have that stage 
            return the new object for the next stage,
            giving the user the option to control the sequence,
            instead of just call all the functions mechanically.
        
        
            -Done: How about generic hooks to store bare hull and bbow?
            -How about ensuring the bulb and hull 
            do not impenge on one another, ever?
            -How about (geo1+geo2=geo3) generic hooks to combine surfaces?? -probably not-
                -problem 1: infinite ways to combine
                -problem 2: ?do we need more problems that problem 1?
            -Avbstract, generalize and make nice:
                -step 1-1: make geometry Lagrangian hackable (a la already done in the GUI)
                -step 1-2: Rules system for Lagrangian curve design editing/hacking
                -step 2: rules for bare hull should be linguistic (a la geometry_generator)
    
    
        
        Consider
        ----------
            -The menu system is really a set of environments 
            -the current environment can be swapped for another via command
                *but only up to a branch
                *or back to the 'root'
            -this is much like the rules system!
                *it even uses the same representation (dict)
            -but there is no need for any type of unification here, right?
                (however exotic or pedestrian it might look)
    
    """
    def __init__(self, design_specification=None, 
                                         design_space=None, verbose=False):
        
        if design_specification is not None:
            self.design_space = DesignSpace(design_specification)
        elif design_space is not None:
            self.design_space = design_space
        self.verbose = verbose
        self.THB_projections = None
        self.bareblockextremes = None
        self.postprocessor = SolverPostProcessor
        #
        #*******************************************
        # interface 
        self.prompt = 'Random|ShipDesigner>'
        self.keep_going = False
        self.command_found = False
        self._current_menu = []
        self.ans_chr = ''
        self.ini = True
        self.blank = "''"
        self.sucess_finding_command_this_time = False
        self.availible_plot_commands = {}
        self.availible_print_commands = {}
        self.context_help = 'program-initialized'
        self.help_dictionary = {}
        self.hp = self.hyperparameter_basic() #default hyper parameters
        
    
    #
    #**************************************************************************
    # Menu Interface
    
    def ShipDesigner(self, start=False,verbose=True):
        """Interactive Program
        
        Menu System:
            Main Menyu
        """
        self.verbose = verbose
        self.keep_going = True
        self.ini = True
        self._current_menu.append(self.main_menu)
        self.lang_loop()
        return
        
#        def splash():
#            print 'interface not started'
#            print 'To start, use any of the following inputs at the prompt:'
#            print  ''
#            for key in doit:
#                print key
#            print '' 
#            print ' Or type '
#            print '>   quit'
#            print ' to quit.'
            
        
#        doit = ['y','t',True]
#        dont = ['quit']
#        if start:
#            self.verbose = verbose
#            self.keep_going = True
#            self._current_menu.append(self.main_menu)
#            self.lang_loop()
#        else:
#            print 'ShipDesigner interface,'
#            self.sinput = raw_input('start?>')
#            try:
#                if self.sinput.split()[0][0].lower() in doit:
#                    self.ShipDesigner(start=True)
#                elif self.sinput.split()[0][0].lower() in dont:
#                    self.end_program()
#                else:
#                    splash()
#                    #return ShipDesigner()
#            except:
#                print 'invalid input'
#                splash()
#                #splash()
#                #ShipDesigner(True)
#        return
        
    def exit_(self, cmnd):
        if self._current_menu:
            discard = self._current_menu.pop()
            #if self.verbose: print 'Leaving = ',discard, ' menu'
            if self._current_menu:
                self._current_menu.pop()(self.blank)
        
        print 'TOP MENU'
        self.main_menu(self.blank)
        return
    
    def input_get_command(self,cmnd):
        """append defaults to the options availible"""
        self._cur_options[self.blank]   = self.command_not_found
        #self._cur_options['end']       = self.exit_
        self._cur_options['back']       = self.exit_
        self._cur_options['help']       = self.help_menu
        self._cur_options['info']       = self.list_commands
        self._cur_options['menu']       = self.list_commands
        self._cur_options['quit']       = self.end_program
        self._cur_options['status']     = self.status
        
        #print 'testing ln 1024 we have cmnd = ',cmnd
        test = cmnd[0].split()[0].lower()
        
        if test in self._cur_options:
            #print 'found command ',test
            self.sucess_finding_command_this_time = True
            return test
        else:
            #print 'cound not find command ',test
            return self.blank
        
        
    
    def help_menu(self, ans_chr=None):
        print '--------------------------------------------------------------'
        print '\n HELP MENU TO BE CONSTRUCTED?\n'
        print '--------------------------------------------------------------'
        self.list_commands()
        
        usage = '   Go back to the last menu'
        self.help_dictionary['back']       = usage
        
        usage = '   Go to the help menu - (you are currently here)'
        self.help_dictionary['help']       = usage
        
        usage = '   Current Command Options'
        self.help_dictionary['info']       = usage
        
        usage = '   Current Command Options'
        self.help_dictionary['menu']       = usage
        
        usage = '   Quit the Program'
        self.help_dictionary['quit']       = usage
        
        usage = '   Current Status'
        self.help_dictionary['status']     = usage
        
        for key in self._cur_options:
            if key in self.help_dictionary:
                print ' command=> ','"',key,'"'+\
                '\n    usage: ',self.help_dictionary[key],'\n'
            else:
                print 'command: ','"',key,'"'+\
                '\n    usage:  N/A; Help information not implemented yet. \n'
        
        
        print '--------------------------------------------------------------'
        print '\ncontext sensitive help:\n'
        print '--------------------------------------------------------------'
        print self.context_help
        print '--------------------------------------------------------------'
        
        return
    
    
    def list_commands(self, ans_chr=None):
        print '--------------------------------------------------------------'
        print '\n The availible commands are:\n'
        print '--------------------------------------------------------------'
        for key in self._cur_options:
            if key is not self.blank:
                print key
        print '--------------------------------------------------------------'
        return
    
        
    def command_not_found(self, ans_chr=None):
        if not self.sucess_finding_command_this_time:
            print 'command not found'
        self.command_found = False
        return 
    
    
    def end_program(self, ans_chr=None):
        print 'Thank you For trying this beta release!'
        self.keep_going = False
        self._current_menu = []
        return
    
    
    def main_menu(self, ans_chr=None):
        """MAIN MENU"""
        self._current_menu.append(self.main_menu)
        self.prompt = 'R|SD>'#.menu.main>'
        
        messg = 'This is the main menu.  \n '+\
                '   Hull Design proceeds from here in a series of steps. \n '+\
                '   (1) Enter the "designspace" menu to use commands to generate a design space. \n '+\
                '   (2) Enter the "hull" menu to use commands to generate a single design. \n '+\
                '           This will be a set of self consistent form parameters. \n '+\
                '           which make up the constraints for the generation \n '+\
                '           of the global bare ship hull form curves. \n '+\
                '           that is: \n '+\
                '                   -SAC (Sectional Area Curve) \n '+\
                '                   -DWL (Design Water Line) \n '+\
                '                   -CPK (Center Plane Keel Profile \n '+\
                '   (3) Enter the "bulb" menu to do the same thing for\n '+\
                '           the design of a bulbous bow.\n '+\
                '   (4) Enter the "combine" menu\n '+\
                '           to ask ShipDesigner to combine the two surfaces\n '+\
                '           into one final Truncated Hierarchical B-spline Surface\n '+\
                ' \n Note:'+\
                '   *Help menus will be coming online gradually. \n'+\
                '         -kind of like the rest of this program!\n '+\
                '        *In a given menu, try typing "help" to print information \n'+\
                '         about what the menu is for.\n '
        self.context_help = messg
            
        self._cur_options = {'designspace':self.design_space_menu,
                             'hull':self.hull_menu,
                             'bulb':self.bow_menu,
                             'plot':self.plot_menu,
                             'combine':self.combine_menu}
        
        
        usage = '   Main Menu of the Program. -No need to type this one.\n'+\
                '   Start here and \n'+\
                '   (1) enter "designspace" to generate a designspace based on rules \n'+\
                '   (2) enter "hull" to generate a designspace based on rules \n'+\
                '   (3) enter "bulb" to generate a designspace based on rules \n'
                
        self.help_dictionary['main_menu']       = usage
        
        if self.ini:
            self.welcome()
            self.ini = False
            self.list_commands()
            return
            
        option = self.input_get_command(ans_chr)#, comnddict=self._cur_options)
        #if self.verbose:print 'found option = ',option
        
        self._cur_options[option](ans_chr)
        return
    
    
    
    
    def design_space_menu(self, ans_chr=None):
        """design space menu"""
        self.prompt = 'R|SD>'#:menu.design.space>'
        
        #if ans_chr[0] == 'designspace':
        #    self._current_menu.append(self.design_space_menu)
        #    return
        
        self._current_menu.append(self.design_space_menu)
        
        
        messg = 'This is the design space menu.  \n '+\
                '   This is step (1) of the design process. \n '+\
                '   *Simply enter "generate_design_space" \n '+\
                '     to generate a design space from default Form Parameter Rules. \n '+\
                '   *Presently default rules are coded in another file; sorry about that!\n '+\
                '    In fact, creation and usage of rules is found in  \n'+\
                '    file:  "simple_hull_rules_langage.py" \n'+\
                '    function:  "make_hull_rules_net" \n'+\
                '    \n'+\
                '   *Instantiate a design variable there with code like this: \n\n'+\
                '      Cb = lp.PStates(name="Cb")   \n\n'+\
                '   *Set a rule with code like this: \n\n'+\
                '      Cb = Cb == vol/(lwl*bwl*draft)  \n\n'+\
                '   (where the other symbols, vol, lwl,bwl,and draft are lp.PStates too.) \n'+\
                '   *Plug a rule into the relational interval database with: \n\n'+\
                '      hullrulesnet.set_rules( Cb )  \n\n'+\
                '   *Tell the rules to update the designspace with \n\n'+\
                '      hullrulesnet.rgp.compute_fresh_rules_graph()  \n\n'+\
                '   *That is it for manual rules creation! \n'+\
                '    \n '+\
                '   *But here in this menu we assume the rules over there '+\
                '    are just the ones we want. \n '+\
                ' \n'+\
                '    So enter "generate_design_space" to use them!'+\
                '    Then enter "back" to go back to the main menu.'
                
                
                
                
                
                
        
        self.context_help = messg
            
        self._cur_options = {'generate_design_space':self.generate_design_space}
            
        option = self.input_get_command(ans_chr)#, comnddict=options)
        #if self.verbose:print 'found option = ',option
        
        self._cur_options[option](ans_chr)
        return
        
    def add_rules(self, ans_chr=None):
        """rules input"""
        return
    
    
    def hull_menu(self, ans_chr=None):
        """hull menu"""
        self.prompt = 'R|SD>'#:menu.hull.design>'
        self._current_menu.append(self.hull_menu)
        #print 'menu:hulldesign'
        #if self.verbose: print 'current input: ',ans_chr
            
        self._cur_options = {'generate_form_parameters':self.bare_hull_random_design_selection_new,
                             'make_geometry':self.make_bare_hull}
            
        option = self.input_get_command(ans_chr)#, comnddict=options)
        #if self.verbose:print 'found option = ',option
        
        self._cur_options[option](ans_chr)
        return
    
    
    
    def bow_menu(self, ans_chr=None):
        """bow menu"""
        self.prompt = 'R|SD>'#:menu.bow.design>'
        self._current_menu.append(self.bow_menu)
        #print 'menu:hulldesign'
        #if self.verbose: print 'current input: ',ans_chr
            
        self._cur_options = {'generate_form_parameters':self.bulbous_bow_random_design_selection,
                             'make_geometry':self.make_bulbous_bow}
            
        option = self.input_get_command(ans_chr)#, comnddict=options)
        #if self.verbose:print 'found option = ',option
        
        self._cur_options[option](ans_chr)
        return
    
    def combine_menu(self, ans_chr=None):
        """bow menu"""
        self.prompt = 'R|SD>'#:menu.bow.design>'
        self._current_menu.append(self.combine_menu)
        #print 'menu:hulldesign'
        #if self.verbose: print 'current input: ',ans_chr
            
        self._cur_options = {'combine':self.make_THBspline_complex_hull}
            
        option = self.input_get_command(ans_chr)#, comnddict=options)
        #if self.verbose:print 'found option = ',option
        
        self._cur_options[option](ans_chr)
        return
    
    
    def plot_menu(self, ans_chr=None):
        """plotting menu"""
        self.prompt = 'R|SD>'#:menu.design.space>'
        
        self._current_menu.append(self.plot_menu)
        
            
        self._cur_options = {'plot_THB':self.plot_THB_ship_hull,
                             'plot_bare_hull':self.plot_bare_hull_results}
            
        option = self.input_get_command(ans_chr)
        
        self._cur_options[option](ans_chr)
        return
    
    
    def print_menu(self, ans_chr=None):
        """plotting menu"""
        self.prompt = 'R|SD>'#:menu.design.space>'
        
        self._current_menu.append(self.plot_menu)
        
            
        self._cur_options = {'bare_hull_design_space':self.plot_THB_ship_hull,
                             'bulb_design_space':self.plot_bare_hull_results}
            
        option = self.input_get_command(ans_chr)
        
        self._cur_options[option](ans_chr)
        return
    
    
    def status(self, ans_chr=None):
        """status function"""
        print '\n current menu:',self._current_menu[-1].__doc__
        
        print ' current options: '
        for op in self._cur_options:
            print '    ',op
        print ''
        print 'Type a command.'
        return 
    
    
    def lang_loop(self):
        """MAIN LOOP
        """
        while (self.keep_going):
            
            if not self._current_menu:
                self._current_menu.append(self.main_menu)
            if self.keep_going:
                self._current_menu.pop()('')
            
            
            while self.keep_going:
                self.sinput = raw_input(self.prompt)
                self.sucess_finding_command_this_time = False
                if not self.sinput:
                    print '>', self.sinput
                    self.keep_going = True
                else:
                    print '>', self.sinput
                    self.ans_chr = self.sinput.split()
                    if not self._current_menu:
                        if self.keep_going:
                            if self.verbose: print 'appending to current menu'
                            self._current_menu.append(self.main_menu)
                        else:
                            break
                    else:                        
                        if self.keep_going:
                            self._current_menu.pop()(self.ans_chr)
                        else:
                            break
                #if self.verbose: print 'current menu in lang loop is = ',self._current_menu
                    
                    #self._current_menu.pop()(self.ans_chr)
                    #self.keep_going = True
        print 'End Program'
        return

    def welcome(self):
        print '--------------------------------------------------------------\n'
        print '   Welcome to the (Random) Ship Hull Designer'
        print ''
        print '      (Currently only random mode is implemented ;) '
        print ''
        print '      But random mode is designed to approximate worst case'
        print '      ship design conditions. '
        print '      ...monkeys banging on a typewriter...'
        print '      This is intended as a proof of concept - that ShipDesigner'
        print '      could aid a machine learning program in designing ships'
        print '      by providing a background structure upon which to iterate'
        print '      and learn.'
        print ''
        print '   So we deomonstrate:'
        print '      Given a design space and nothing else, '
        print '      ...'
        print '      Generate a valid ship model!'
        print ''
        print ''
        print ' enter a command to get started:'
        print '--------------------------------------------------------------\n'
        return


    def design_commands(self, ncom):
        comnd = ['design space','hull','bulb','combine']
        return
    
    #
    #**************************************************************************
    # Hyperparameter routines
    def hyperparameter_basic(self):
        class Messenger:
            """python data 'pattern'
            http://python-3-patterns-idioms-test.readthedocs.io/en/latest/Messenger.html
            """
            def __init__(self, **kwargs):
                self.__dict__ = kwargs
                
        hp = Messenger(midship_design_vector = {1:'s',
                                                2:'s',
                                                3:'s',},
                       fwd_fairness_location_ratio      = .3,
                       fwd_transition_location_ratio    = .6,
                       aft_transition_location_ratio    = .4,
                       aft_fairness_location_ratio      = .65,
                       remove_fwdfairness = False,
                       heavy_config = False,
                       remove_aft_fos = False,
                       #
                       #***************************
                       # automation parameters:
                       #
                       ignore = [3,10],
                       dropout = [6,7],
                       use_in_linear_interp = [0,1,2,
                                               4,5,
                                               6,7,
                                               8,9,
                                               11,13],
                       use_in_longi_solver = [0,1,2,
                                               3,4,5,
                                               6,7,
                                               8,9,10,
                                               11,12,13],
                       aft_drop = .25,
                       bulb_volume_fraction = .35,
                       bulb_length_fraction = .67
                       )
        return hp
    
    
    def hyperparameter_search(self):
        """
        inputs
        ----------
            (when finished this will be) 
            ...designed to be used just after a feasible set of parametric
            constraints has been found using the interval constraint rules
            network (so we have 1 design)
            
            Basically freeze the design and allow random, or learned
            search over 
            the hyperparameters 
            of the form parameter design class
            
            using interval constraint spaces again, of course!
            
            To be re-written with the following in mind:
                
            Idea:  basically label all possible transverse curves with
            a unique index, 0-14 inclusive,
            
            Then input here can designate
            - which curves for linear interpolation
            - which curves for equality constraints in the longi solver
            - which curves for least squares approximation (LS) in the longi solver
            
            and, independent of the curve choice we have 
            
            - nose constraint
            - [0.,1.] placement of two fairness curves fore and aft (some may be turned off/switched to LS)
            
            This stuff constitutes a mixed integer 
            architecture specification set.
            
            Next we need an objective function
        """
        class hyptopt(object):
            def __init__(self):
                self.fwd_fairness = lp.PStates(name='fwd_fairness')
                self.fwd_transition = lp.PStates(name='fwd_transition')
                self.aft_fairness = lp.PStates(name='aft_fairness')
                self.aft_transition = lp.PStates(name='aft_transition')
                self.midship_fos_factor1 = ['l','s','h']
                self.midship_fos_factor2 = ['l','s','h']
                self.midship_fos_factor3 = ['l','s','h']
                self.extract_midship_curve_list = []# currently expressed as 
                                                    # FPD Hull attribute self.heavy_config
                self.extract_aft_curve_list = []    # currently expressed 
                                                    # as FPD Hull attribute self.remove_aft_fos
                self.extract_bow_curve_list = []    #nothing implemented on this front
            
            
            def _setup(self):
                
                #
                #*******************************************
                #
                return
        
        
        ff = lp.Variable('fwdfairness')
        ft = lp.Variable('fwdtransition')
        af = lp.Variable('aftfairness')
        at = lp.Variable('afttransition')
        mff1 = lp.Variable('mff1')
        mff2 = lp.Variable('mff2')
        mff3 = lp.Variable('mff3')
        
        quicklist = [ff,ft,af,at,mff1,mff2,mff3]
        
        s = lp.State(values={ff:None,
                          ft:None,
                          af:None,
                          at:None,
                          mff1:None,
                          mff2:None,
                          mff3:None})
        hyptuples = lp.States(s)
        
        
        hclass = hyptopt()
        
        
        alt = lp.State(values={ff:hclass.fwd_fairness,
                          ft:hclass.fwd_transition,
                          af:hclass.aft_fairness,
                          at:hclass.aft_transition,
                          mff1:None,
                          mff2:None,
                          mff3:None})
        #def with_vars(vars):
        #    ['l','s','h'] = vars
        #    return Goal.eq(from_list([x,2,z]),
        #                      from_list([1,y,3]))
        #goal = lp.Goal.bind(hclass.midship_fos_factor2,mff2)
        #alt = goal(alt)
        
        
        
        test = copy.deepcopy(hyptuples)
        for vr in quicklist:
            test = (test == ( vr, alt(vr) ) )
        
        
        hclass.ff = ff
        hclass.ft = ft
        hclass.af = af
        hclass.at = at
        hclass.mff1 = mff1
        hclass.mff2 = mff2
        hclass.mff3 = mff3
        
        return hyptuples, hclass
    #
    #**************************************************************************
    # Design routines
    
    def generate_design_space(self, ans_chr=None):
        #print '\n menu:designspace \n'
        if self.verbose: print 'Generating Fresh Design Space'
        import simple_hull_rules_language as shr
        make_hull_rules_net = shr.make_hull_rules_net
        combine_space_with_rules = shr.combine_space_with_rules
        HullGeometryGenerator = shr.HullGeometryGenerator
        DS = DesignSpecification(
                            lwl = (50.,130.),
                            draft = (15.,25.),
                            bwl = (25.,40.),
                            #vol=(18000.,20000.),
                            #vol=(22000.,28000.),
                            #vol=(12000.,25000.),
                            vol=(5000.,25000.),
                            LCG = (45.,68.),
                            Clcg = (.485,.515), #location Coeff for the LCG
                            Cb = (0.5,0.89),
                            Cwp = (.84,.94),#water plane coefficient
                            Cmidshp = (0.7,.96), #midship coefficient
                            Ccp = (.84,.94) #centerplane coefficient
                            )
        hullrulesnet = make_hull_rules_net(HullGeometryGenerator(rtype='gauss',
                                                             verbose=True))
        self.design_space = combine_space_with_rules(DS,
                                                hullrulesnet)
        return
    
    def get_bare_hull_extremes(self, dmy=None):
        """extremes of the symetric portion of the hull
        -to be done only after narrowing is complete
        
        self.lX => length extremes
        self.bx => half width extremes
        self.dx => draft extremes
        """
        message = 'ShipDesigner: get_bare_hull_extremes'
        if self.verbose: print message
        #
        ck = self.design_space.get_thick_list()
        if ck: print 'ERROR, warning: '+\
                        'hull extremes being taken from '+\
                        'a thick interval valued design space'
        lxb = 0.
        lxe = self.design_space.get_from_rgp('lwl')[1][0][1]
        byb = 0.
        bye = self.design_space.get_from_rgp('bwl')[1][0][1]/2.
        dzb = 0.
        dze = self.design_space.get_from_rgp('draft')[1][0][1]
        class blockextremes(object):
            """Tiny class holding data on the 
            extremes of the symetric half-bare-hull
            """
            def __init__(self, lb,le,bb,be,db,de):
                self.lb = lb
                self.le = le
                self.bb = bb
                self.be = be
                self.db = db
                self.de = de
            def __str__(self):
                print 'extreme lengths: ',self.lb, self.le
                print 'extreme widths: ',self.bb, self.be
                print 'extreme depths: ',self.db, self.de
                return 
        return blockextremes(lxb,lxe,byb,bye,dzb,dze)
        
    
    #
    #**************************************************************************
    # Design routines
    def bare_hull_random_design_selection_new(self,
                                              dmy=None,
                                              maxiter=None,
                                              hyperparameters = None,
                                              hand_tune=False):
        """generate a full set of thin hull constraints 
        for one single hull
        via random.random selection
        from the total design space designated via language rules
        using the underlying sqKanren logic and rules graph processing
        
        Parameters
        ----------
            design space
        Returns
        ----------
            - a single thin design from within the design space
            - this design is sutable for Form Paramter Design
        Notes
        ----------
        
        """
        self.design_space.tree_search(maxACiter=maxiter)
        #
        # nput of hyperparameters at hull instantiation:
        self.thin_hull_design = hullmethods.Hull(
                                    States = self.design_space.rgp.env, 
                                    which_state=0,
                                    revise=False,
                                    hyperparameter_input = hyperparameters)
        #but maybe easier to see if you use the other interfaces?
        #
        for name in self.design_space.__dict__:
            dname = getattr(self.design_space, name)
            if isinstance(dname, lp.PStates):
                if self.verbose: print dname
                mkname,env_value = self.design_space.get_from_rgp(dname)
                if self.verbose: print 'setting PStates',dname.name,' => mkVar ',mkname,' : ',env_value,'\n\n\n'
                setattr(self.thin_hull_design.hullconstraints, 
                        name,
                        mkname)
        self.bareblockextremes = self.get_bare_hull_extremes()
        print '-----------------------------------------------------'
        print 'Design Space selection Complete'
        return
    
    #
    #**************************************************************************
    # Design routines
    def bare_hull_random_design_selection(self,dmy=None,maxiter=None):
        """generate a full set of thin hull constraints 
        for one single hull
        via sobol random selection
        from the total design space designated via hard coded rules
        using the underlying sqKanren logic
        
        Parameters
        ----------
            design space
        Returns
        ----------
            - a single thin design from within the design space
            - this design is sutable for Form Paramter Design
        Notes
        ----------
        
        """
        message = '\nbare_hull_random_design_selection: \n'+\
        ' generate a full set of thin hull constraints \n'
        if self.verbose: print message
        
        self.design_space.tree_search(maxACiter=maxiter)
        self.design_space.print_widths()
        ok = False
        cnt= -1
        while True:
            cnt += 1
            print cnt
            ok = self.design_space.check_widths()
            if ok:
                hull1 = hullmethods.Hull(States = self.design_space.hdp, 
                                         which_state=0,
                                         revise=False)
                break
            else:
                #            print 'no good'
                #            break
                try:
                    self.design_space.tree_search()
                except:
                    print 'error, climb the tree'
                    kparent = self.design_sq.get_latest_child().parent
                    self.design_space.hdp.states = kparent.states
            if cnt>5:
                print 'tree count > 5'
                break
        if ok:
            print 'feasible design reached'
            self.thin_hull_design = hull1
        else:
            print 'ERROR: warning, possibly infeasible hullform'
            self.thin_hull_design = hullmethods.Hull(
                                        States = self.design_space.hdp, 
                                        which_state=0,
                                        revise=False)
        self.bareblockextremes = self.get_bare_hull_extremes()
        return
    
    #
    #**************************************************************************
    # Design routines
    def bare_hull_random_design_selection_simplified(self,dmy=None):
        """Like the above (for the old style rules), but without try except while loop
        """
        message = '\nBare_hull_random_design_selection_simplified: \n '+\
        ' generate a full set of thin hull constraints \n'
        if self.verbose: print message
        
        self.design_space.tree_search()
        self.design_space.print_widths()
        print '\n---------------------------------------------------------------\n'
        print 'initial design space:'
        self.design_space.print_widths()
        print '\n---------------------------------------------------------------\n'
        self.thin_hull_design = hullmethods.Hull(States = self.design_space.hdp, 
                                                 which_state=0,
                                                 revise=False)
        self.bareblockextremes = self.get_bare_hull_extremes()
        return
        
    #
    #**************************************************************************
    # Design routines
    def bulbous_bow_random_design_selection(self,dmy=None):
        """generate a full set of thin bulbous bow constraints 
        for one single bulbous bow
        via random selection
        from the total design space
        
        Parameters
        ----------
            design space
        Returns
        ----------
            - a single thin design from within the design space
            - this design is sutable for Form Paramter Design
        Notes
        ----------
        
        """
        message = '\nBulbous_bow_random_design_selection: \ngenerate a full set '+\
        'of thin bulbous bow constraints \n'
        if self.verbose: print message
        
        #from bare hull make random bbow thin design parameter set:
        #self.thin_bbow_design = interface_to_bbow_make_thin_bbow(self.design_space.hdp) 
        self.thin_bbow_design = \
                interface_to_bbow_make_thin_bbow(
                            self.thin_hull_design.hullconstraints,
                            bbow_interface_section_area = None,
                            bbow_length = None) 
        
        #self.thin_bbow_design = \
        #        interface_to_bbow_make_thin_bbow(
        #                    self.thin_hull_design.hullconstraints,
        #                    bbow_interface_section_area = self.hull.A_mid,
        #                    bbow_length = self.hull.BulbLength) 
        
        return
    
    
    #
    #**************************************************************************
    # Design routines
    def bare_hull_learning_design_selection(self,dmy=None):
        """TBD
        """
        return
    
    def bulbous_bow_learning_design_selection(self,dmy=None):
        """TBD
        """
        return
     
    
    #
    #**************************************************************************
    # Design routines
    def make_bare_hull(self, dmy=None, thin_hull_design=None,
                       hyperparameters = None):
        """Create a bare hull B-spline surface from
        a thin set of hull constraints
        
        Parameters
        ----------
        Returns
        ----------
        Notes
        ----------
            Use only after narrowing the 
            bare hull design space to
            a thin design... by either random selection
            or learned selection
            
        dev
        ----------
            import hull_from_simple_designspace as hullmethods 
        """
        message = '\n make_bare_hull: \nCreate a bare hull B-spline surface from\n'+\
        ' a thin set of hull constraints \n via FPD \n'
        #if self.verbose: print message
        print message
        
        if thin_hull_design is None: 
            thin_hull_design = self.thin_hull_design
        # --- found in function 'give_results():'
        hull2 = copy.deepcopy(thin_hull_design)
        #**************************************************
        # where to put the check...hello monads
        if hyperparameters is not None:
            hull2.change_hyperparameters(hyperparameters)
        #**************************************************
        hull2.define_hull()
        hull2.definestations()
        
        hullmethods.compute_hull_and_transverse_curves(hull2,
                                                       sacflat_=True)
        self.hull = hull2
        hullmethods.loft_hull(hull2,
                              refine=True,
                              extrema=self.bareblockextremes)
        # ---
        self.hull = hull2
        return
    
    
    
    #
    #**************************************************************************
    # Design routines
    def make_bulbous_bow(self,dmy=None):
        """FPD for bbow
        
            Use only after narrowing the 
            bulbous bow design space to generate
            a thin design... by either random selection
            or learned selection
        """
        message = '\nMake_bulbous_bow: \ngenerate a B-spline bulbous bow \n '+\
        'from a thin set of bbow constraints \n via FPD \n'
        if self.verbose: print message
        
        self.bbow = interface_to_bbow_make_bare_hull(self.thin_bbow_design,
                                                     skip_search=dmy)
        return
    
    
    
        
    
    #
    #**************************************************************************
    # Integrated Design routines
    def complex_hull_design_and_generation(self,
                                           dmy=None,
                                           maxiter=None,
                                           hyperparameters = None,
                                           hand_tune=False):
        
        ## SD.bare_hull_random_design_selection_new
        """generate a full set of thin hull constraints 
        for one single hull
        via random.random selection
        from the total design space designated via language rules
        using the underlying sqKanren logic and rules graph processing
        
        Parameters
        ----------
            design space
        Returns
        ----------
            - a single thin design from within the design space
            - this design is sutable for Form Paramter Design
        Notes
        ----------
        
        """
        message = '\ncomplex_hull_design_and_generation: \ngenerate a full set '+\
        'of bare hull bow constraints \n'+\
        'and then bbow constraints. \n'+\
        'Then use them in a consistent hull generation scheme \n'
        #if self.verbose: print message
        print message
        #
        
        
        ##
        ## Basic Hull Random Design selection
        ##
        self.design_space.tree_search(maxACiter=maxiter)
        
        
        #
        # nput of hyperparameters at hull instantiation:
        self.thin_hull_design = hullmethods.Hull(
                                    States = self.design_space.rgp.env, 
                                    which_state=0,
                                    revise=False,
                                    hyperparameter_input = hyperparameters)
        
        
        #but maybe easier to see if you use the other interfaces?
        #
        for name in self.design_space.__dict__:
            dname = getattr(self.design_space, name)
            if isinstance(dname, lp.PStates):
                if self.verbose: print dname
                mkname,env_value = self.design_space.get_from_rgp(dname)
                if self.verbose: print 'setting PStates',dname.name,' => mkVar ',mkname,' : ',env_value,'\n\n\n'
                setattr(self.thin_hull_design.hullconstraints, 
                        name,
                        mkname)
        self.bareblockextremes = self.get_bare_hull_extremes()
        print '-----------------------------------------------------'
        print 'Design Space selection Complete'
        
        
        
        ##
        ## Do Bulbous Bow:
        ##
        """generate a full set of thin bulbous bow constraints 
        for one single bulbous bow
        via random selection
        from the total design space
        
        Parameters
        ----------
            design space
        Returns
        ----------
            - a single thin design from within the design space
            - this design is sutable for Form Paramter Design
        Notes
        ----------
        
        """
        
        message = '\nBulbous_bow_random_design_selection: \ngenerate a full set '+\
        'of thin bulbous bow constraints \n'
        if self.verbose: print message
        
        #from bare hull make random bbow thin design parameter set:
        #self.thin_bbow_design = interface_to_bbow_make_thin_bbow(self.design_space.hdp) 
        self.thin_bbow_design = \
                interface_to_bbow_make_thin_bbow(
                            self.design_space.rgp.env,
                            bbow_interface_section_area = None,
                            bbow_length = None) 
        
        self.thin_bbow_design.tree_search() 
        
        
        
        
        
        
        ##
        ## bare_hull_random_design_selection_new:
        """generate a full set of thin hull constraints 
        for one single hull
        via random.random selection
        from the total design space designated via language rules
        using the underlying sqKanren logic and rules graph processing
        
        Parameters
        ----------
            design space
        Returns
        ----------
            - a single thin design from within the design space
            - this design is sutable for Form Paramter Design
        Notes
        ----------
        
        """
        self.thin_hull_design = hullmethods.Hull(
                                    States = self.design_space.rgp.env, 
                                    which_state=0,
                                    revise=False,
                                    hyperparameter_input = hyperparameters,
                                    bulbconstraints = self.thin_bbow_design.rgp)
        #but maybe easier to see if you use the other interfaces?
        #
        for name in self.design_space.__dict__:
            dname = getattr(self.design_space, name)
            if isinstance(dname, lp.PStates):
                if self.verbose: print dname
                mkname,env_value = self.design_space.get_from_rgp(dname)
                if self.verbose: print 'setting PStates',dname.name,' => mkVar ',mkname,' : ',env_value,'\n\n\n'
                setattr(self.thin_hull_design.hullconstraints, 
                        name,
                        mkname)
        self.bareblockextremes = self.get_bare_hull_extremes()
        print '-----------------------------------------------------'
        print 'Design Space selection Complete'
        
        

        return
    
    #
    #**************************************************************************
    # Design routines
    def make_THB_bare_hull(self,dmy=None):
        """
        Parameters
        ----------
        Returns
        ----------
        Notes
        ----------
        Issues
        ----------
            *presently if an inner level has the same 
            outer upper bounds in some DOF, as it's parent, 
            this causes an 
            error in plotting.
        
        
        Dev
        ----------
            import thbsurface as thbspline
            import relational_lsplines
            BoxList = relational_lsplines.iaBox.BoxList
            from myspline import rbspline
            
            from opt_simple_hull import get_projection_matrices
            
            box = thbspline.Box( ia(.0,.8125),ia(.0,1.) )
            
            box = thbspline.Box( ia(.0,1.),ia(.0,1.) )
            
            box1 = thbspline.Box( ia(.0,.8125), ia(.0,.34375) )
            box2 = thbspline.Box( ia(0.,.4375), ia(.34375,1.) )  
            
            FLIPED:
                
            box1 = thbspline.Box( ia(.0, .34375) , ia(.0, .8125) )
            box2 = thbspline.Box( ia(.34375, 1.) , ia(0., .4375) )  
            
            thbspline.BoxList( box1,box2 )
            
            
            changing from .8 to .5
            
            0.34375
            
            SD.THBhull.child.child.bounds[0][0].sup = .34375
            SD.THBhull.child.child.bounds[0][1].sup = .15
            
            
            SD.THBhull.child.child.bounds[1][0].inf = 0.34375
            SD.THBhull.child.child.bounds[1][0].sup = 1.
            
            SD.THBhull.child.child.bounds[1][1].inf = 0.
            SD.THBhull.child.child.bounds[1][1].sup = .1
            
            
            SD.THBhull.child.child.bounds
            Out[19]: BoxList((Box((ia(0.0, 0.34375), ia(0.0, 0.5))), Box((ia(0.34375, 1.0), ia(0.0, 0.3)))))
        
        
        ISSUE!
        -you must feed this animal something which has 11
        transverse curves!
        """
        #not strictly needed
        from myspline import rbspline
        art_tnet = []
        for i in range(self.hull.lcurvenet[0].n):
            art_tnet.append(rbspline(self.hull.hullsurf.surfNet[:,i],
                                     k=4,
                                     nump=30))
        #
        #
        #hull2 = self.hull
        if self.THB_projections is None:
            self.THB_projections = get_projection_matrices(
                            thbspline.THButilities.ProjectionGetter())
        
        t1 = thbspline.THBsurface(vertices = self.hull.hullsurf.surfNet,
                                  tcurvelist = art_tnet,
                                  #tcurvelist = self.hull.tcurvenet,
                                  ulevel = 2,
                                  lcurvelist = self.hull.lcurvenet,
                                  vlevel = 3,
                                  hullcurvelist=self.hull.tcurvenet, #TODO: be careful - I've added this Jan5, 12:48 am
                                  trm = self.THB_projections.rm_11x7,
                                  lrm = self.THB_projections.rm_19x11)
        
        t1.expand_initial_curve_list(override=True) #moved here from THBhull+THBbulb function 
        t1.set_top()
        
        
        
        #box1 = thbspline.Box( ia(.0, .0) , ia(.0, .0) )
        #box2 = thbspline.Box( ia(.0, 0.) , ia(0., .0) ) 
        
        #complications from not nesting the surface are "high"
        box1 = thbspline.Box( ia(.0, .34375) , ia(.0, .15) )
        box2 = thbspline.Box( ia(.34375, 1.) , ia(0., .1) )  
        
        
        t11 = t1.dyadic_loft_refinement(
                bounds = thbspline.BoxList( box1,box2),
                            trm = self.THB_projections.rm_19x11,
                            lrm = self.THB_projections.rm_35x19
                        )
                    
        #box1 = thbspline.Box( ia(.0, .34375) , ia(.0, .8125) )
        #box2 = thbspline.Box( ia(.34375, 1.) , ia(0., .4375) )  
        
        
        box1 = thbspline.Box( ia(.0, .34375) , ia(.0, .15) )
        box2 = thbspline.Box( ia(.34375, 1.) , ia(0., .1) )  
        
        
        t12 = t11.dyadic_loft_refinement(bounds = thbspline.BoxList(box1,box2),  
                                        trm = self.THB_projections.rm_35x19,
                                        lrm = self.THB_projections.rm_67x35
                                        ) 
                    
        t1.compute_surface()
        self.THBhull = t1
        return
    
    """
    OLD STUFF:
        
        
        t11 = t1.dyadic_loft_refinement(
                    bounds = thbspline.BoxList(  
                            thbspline.Box( ia(.0,.8125),      ##.5625
                                           ia(.0,.34375) ) ),  ##.3125
                            trm = self.THB_projections.rm_19x11,
                            lrm = self.THB_projections.rm_35x19
                            )
        
        t12 = t11.dyadic_loft_refinement(bounds = thbspline.BoxList(  
                            thbspline.Box( ia(.0,.8125),      #8125  #.4375 ## 7-> 7+3+3=13
                                           ia(.0,.34375) ) ),  #.40625  #.3125    ## 5 -> 5+3+3 = 11
                            trm = self.THB_projections.rm_35x19,
                            lrm = self.THB_projections.rm_67x35
                            ) 
        
    """
    
    """
    NEW STUFF:
        
        
        t11 = t1.dyadic_loft_refinement(
                    bounds = thbspline.BoxList(  
                            thbspline.Box( ia(.0,.8125),      ##.5625
                                           ia(.0,.34375) ),  ##.3125
                            thbspline.Box( ia(0.,.4375),      
                                           ia(.34375,1.) )  
                                ),
                            trm = self.THB_projections.rm_19x11,
                            lrm = self.THB_projections.rm_35x19
                            )
                    
        t12 = t11.dyadic_loft_refinement(bounds = thbspline.BoxList(
                            thbspline.Box( ia(.0,.8125),      #8125  #.4375 ## 7-> 7+3+3=13
                                               ia(.0,.34375) ),  #.40625  #.3125    ## 5 -> 5+3+3 = 11
                            thbspline.Box( ia(0.,.4375),      
                                           ia(.34375,1.) )  
                                ),  
                            trm = self.THB_projections.rm_35x19,
                            lrm = self.THB_projections.rm_67x35
                            ) 
                            
            
        box = thbspline.Box( ia(.0,.8125),ia(.0,1) )
        
        box1 = thbspline.Box( ia(.0,.8125), ia(.0,.34375) )
        box2 = thbspline.Box( ia(0.,.4375), ia(.34375,1.) )  
        
        thbspline.BoxList( box1,box2 )
        
        
    """
    
    """
    HOWTO CHANGE BOUNDS:
        
        SD.THBhull.child.child.set_bounds(thbspline.BoxList(
                            thbspline.Box( ia(.0,.8125),
                                           ia(.0,.34375) ),
                            thbspline.Box( ia(0.,.4375),      
                                           ia(.34375,1.) )  
                                ))
                            
                            
                            
    """
    
    
    #
    #**************************************************************************
    # Design routines
    
    def make_THB_bulb(self,dmy=None):
        """
        Parameters
        ----------
            dmy : interactive input (not yet used)
        Returns
        ----------
        Notes
        ----------
        """
        
        bbow = self.bbow
        if self.THB_projections is None:
            self.THB_projections = get_projection_matrices(
                        thbspline.THButilities.ProjectionGetter())
        
        t2 = thbspline.THBsurface(vertices = self.bbow.hullsurf.surfNet,
                                  tcurvelist = self.bbow.tcurvenet,
                                  ulevel = 2,
                                  lcurvelist = self.bbow.lcurvenet,
                                  vlevel = 1,
                                  trm = self.THB_projections.rm_11x7,
                                  lrm = self.THB_projections.rm_7x5)
        t3 = t2.reflect_surface((0.,0.,1.),(0.,0.,1.))
        t3 = t3.reflect_surface((1.,0.,0.),(1.,0.,0.))
        t3 = t3.rotate_surface((0.,0.,1.),np.pi/4.)
        t3 = t3.translate_surface( dz=-np.amin(t3.vertices[:,:,2]) )
        t3.set_top()
        
        self.THBbow = t3
        return
    
    def translate_THB_bulb(self, surface, translation_vector= None, ans_chr=None):
        """uses translate_surface from BsplineSurface class in curve.py 
        -which in turn uses the update_loft helper function
        
        -this is untested on higher levels.  
        -but is safe on the lowest level.
        
        """
        surf = copy.deepcopy(surface)
        try:
            if translation_vector is None:
                surf = surf.translate_surface( 
                        dz=-np.amin(surf.vertices[:,:,2]) )
            else:
                dx_ = translation_vector[0]
                dy_ = translation_vector[1]
                dz_ = translation_vector[2]
                surf = surf.translate_surface(dx=dx_, dy=dy_, dz=dz_)
        except:
            print 'translation failed'
            surf = copy.deepcopy(surface)
        return surf
    
    #
    #**************************************************************************
    # Design routines
    def make_THB_nose_flat(self,dmy=None):
        """
        Parameters
        ----------
            dmy : interactive input (not yet used)
            
        Returns
        ----------
            
        
        Action
        ----------
            gives a flat nose to the bare hull in THB form
            
            
        Dev
        ----------
            import thbsurface as thbspline
            BoxList = relational_lsplines.iaBox.BoxList
            
            
            box = thbspline.Box( ia(.0,.8125),ia(.0,1.) )
            
            box1 = thbspline.Box( ia(.0,.8125), ia(.0,.34375) )
            box2 = thbspline.Box( ia(0.,.4375), ia(.34375,1.) )  
            
            thbspline.BoxList( box1,box2 )
        
        """
        #        hull = self.THBhull
        #        tapb = hull.compute_total_active_by_level()
        #        
        #        
        #        t3 = copy.deepcopy(self.THBhull)
        #        
        #        
        #        tapb = hull.compute_total_active_by_level()
        #        
        #        t3_tapb = tapb[2]
        #        t3.vertices[t3_tapb.uactive[0],t3_tapb.vactive[0]]
        
        
        t0vertices = self.THBhull.child.child.tcurvelist[0].vertices
        
        self.THBhull.child.child.tcurvelist[1].vertices[1:,2] = t0vertices[1:,2]
        
        
        self.THBhull.child.child.set_vertices_from_tcurves()
        """ now set the longitudinals from the vertices:
        """
        self.THBhull.child.child.set_THB_curves_from_vertex_net(override_u=False,
                                                                override_v=True)
        
        # now
        # self.THBhull.child.child.vertices[:,0]
        # and
        # self.THBhull.child.child.vertices[:,1]
        # match
        return
    
    #
    #**************************************************************************
    # Design routines
    def add_THB_bulbous_bow_to_bare_hull(self,dmy=None):
        """
        Parameters
        ----------
            Fully generated THB hull and THB bow
        Returns
        ----------
            creates a THB surface incorporating both
            the original bare hull
            and the original bulbous bow
        Notes
        ----------
        (1) You must account for the fact that the transverse
                hull curves control vertices must incorporate some of what was the
                bbow longitudinal topmost curve control vertices.
        
        (2) Therefore, in order to match the shape properly,
                as we move progressively aft, the next transverse hull curve 
                collects progressively more and more on the topmost
                longitudinal of the bulbous bow.
                
        (3) The longitudinal CP _of _the _bbow runs aft to fwd
                but the transverse cp curve runs bottom to top,
                corresponding to running fwd to aft, along this portion.
                
            
        """
        hull = self.THBhull
        t3 = self.THBbow
        hull.set_top()
        
        tapb = hull.compute_total_active_by_level()
        
        t3_tapb = tapb[2]
        t3.vertices[t3_tapb.uactive[0],t3_tapb.vactive[0]]
        
        """bblcp1 will be a crucial component here...
        """
        bblcp1 = t3.lcurvelist[-1].vertices[::-1]  #flip top longitudinal bbow curve in z as per (3) above
        
        """
            *Match the fwd transverse control vertices of the bare hull, 
                (just the indicies that correspond with the bbow transverse vertices)
            *But you must account for the fact that the transverse
                hull curve control vertices must incorporate some of what was the
                bbow longitudinal topmost curve control vertices.
            *Align the vertices[n1:n1+n2] leftover vertices with the 
                topmost longitudinal control vertices in the bbow control net
        """
        #1.) fully forward
        n1 = t3.tcurvelist[-1].n
        n2 = t3.lcurvelist[-1].n-0
        hull.child.child.tcurvelist[0].vertices[0:n1] = t3.tcurvelist[-1].vertices
        hull.child.child.tcurvelist[0].vertices[n1:n1+n2] = bblcp1
        
        """
            *Match the next set of transverse control vertices of the bare hull, 
                (just the indicies that correspond with the bbow transverse vertices 
                at the same column number)
            *Align the vertices[n1:n1+n2] leftover vertices with the 
                topmost longitudinal control vertices in the bbow control net
        """
        #2.) nextmost
        n1 = t3.tcurvelist[-2].n
        n2 = t3.lcurvelist[-1].n-1
        hull.child.child.tcurvelist[1].vertices[0:n1] = t3.tcurvelist[-2].vertices
        hull.child.child.tcurvelist[1].vertices[n1:n1+n2] = bblcp1[1:]
        """match up with transverses.  Run along the top bbow 
        longitudinal where needed
        """
        #3.) nextmost
        n1 = t3.tcurvelist[-3].n
        n2 = t3.lcurvelist[-1].n-2
        hull.child.child.tcurvelist[2].vertices[0:n1] = t3.tcurvelist[-3].vertices
        hull.child.child.tcurvelist[2].vertices[n1:n1+n2] = bblcp1[2:]
        """match up with transverses.  Run along the top bbow 
        longitudinal where needed
        """
        #4.) nextmost
        n1 = t3.tcurvelist[-4].n
        n2 = t3.lcurvelist[-1].n-3
        hull.child.child.tcurvelist[3].vertices[0:n1] = t3.tcurvelist[-4].vertices
        hull.child.child.tcurvelist[3].vertices[n1:n1+n2] = bblcp1[3:]
        """match up with transverses.  Run along the top bbow 
        longitudinal where needed
        """
        try:
            #5.) nextmost
            n1 = t3.tcurvelist[-5].n
            n2 = t3.lcurvelist[-1].n-4
            hull.child.child.tcurvelist[4].vertices[0:n1] = t3.tcurvelist[-5].vertices
            hull.child.child.tcurvelist[4].vertices[n1:n1+n2] = bblcp1[4]
        except:
            print 'reduced vertices for bbow'
        #
        """Fix the hull vertices from the tcurves:
        """
        hull.child.child.set_vertices_from_tcurves()
        """ now set the longitudinals from the vertices:
        """
        hull.child.child.set_THB_curves_from_vertex_net(override_u=False,
                                                        override_v=True)
        """fixup the base curves?  --This didn't have to be done right now
        """
        #hull.expand_initial_curve_list(override=True)
        self.THBhull = hull
        return
    
    
    #
    #**************************************************************************
    # Design routines
    def make_THBspline_complex_hull(self,dmy=None):
        """wrapper to combine hulls as THB objects
        -make_THB_bare_hull()
        -make_THB_bulb()
        -add_THB_bulbous_bow_to_bare_hull()
        
        Parameters
        ----------
        Returns
        ----------
        Notes
        ----------
            blending curves, frenet frame matching
            boolean blends, bilinear blends, 
            (coons, hermite) =>  (quadratic, cubic) see Nowaki page 170
        """
        message = '\nMake_THBspline_complex_hull: \ngenerate the full hull \n '+\
        'from the bare hull and the bulbous bow \n Represent this as a THBspline \n'
        if self.verbose: print message
        
        if self.THB_projections is None:
            self.THB_projections = get_projection_matrices(
                                        thbspline.THButilities.ProjectionGetter())
        print '\n make THB bare hull'
        self.make_THB_bare_hull()
        print '\n nose to 90'
        self.make_THB_nose_flat()
        print '\n make THB bulbous bow'
        self.make_THB_bulb()
        self.add_THB_bulbous_bow_to_bare_hull()
        
        return
    #
    #**************************************************************************
    # THB with Lspline Optimization
    def THB_Lspline(self):
        """Main Idea:
            -Decide to fair some level of the THB curves 
            and enforce some constraints.
                (I recommned the level with highest fineness 
                OR 
                to optize a lower level then project the changes to the
                higher level vertices before optimizing that level
                for yet details which are even more fine grained)
            
            -Hold fixed any vertices not active at that level
            Optimize the 'whole curve' in this partially fixed condition 
            (the matrix will be small)
            -Fix large scale constraints at their starting values, or small
            deviations therof.
            -This is primarily for local constraints, like derivatives,
            or local interpolation.
            
        """
        return
    def THB_CenterFOS(self):
        """
        use the current intermediate representation
        to enforce flat of side.
        """
        box1 = thbspline.Box( ia(.0, 1.) , ia(.0, 1.) )
        box2 = thbspline.Box( ia(0., 1.) , ia(0., .1) )  
        bounds = thbspline.BoxList( box1,box2)
        
        self.THBhull.child.bounds = bounds
        return
    #
    #**************************************************************************
    # Plotting 
    def plot_hull_curves_of_form(self, dmy=None,  hull=None):
        """
            parameters:
                hull : hull form generated with relational_lsplines
        """
        if hull is None: hull = self.hull
        hull.Lsplines.DWL.curve.plot()
        hull.Lsplines.CProfile.curve.plot()
        hull.Lsplines.SAC.curve.plot()
        return
    
    def plot_bare_hull_results(self, ans_chr=None):
        """plot the transverse and longitudinal 
        curves which define the bare hull
        (the transverse curves are the true hull
        transverses.  Longitudinals 
        interpolate the transverse vertices and
        thus make up the hull surface control vertex net)
        
        import hull_from_simple_designspace as hullmethods
        """
        hull2 = self.hull
        self.plot_hull_curves_of_form(hull2)
        
        #hullmethods.hull_gui(hull2)            #Lspline-able curves
        hullmethods.hull_gui_no_Lsplines(hull2) #forego it to get 3D
        
        curve = hull2.DWL
        #xb,yb,zb = self.plot_get_extremes()
        biggest = hull2.LengthDWL
        xb = (-biggest*.0,1.*biggest)
        yb = (-biggest*.5,.5*biggest)
        zb = (-biggest*.5,.5*biggest)
        curve.plot3DmultiList(hull2.tcurvenet,
                              hull2.lcurvenet,
                              view='x-vertical',
                              minx=xb[0],
                              miny=yb[0],
                              minz=zb[0],
                              limx=xb[1],
                              limy=yb[1],
                              limz=zb[1])
        
        #        curve.plot3DmultiList(hull2.tcurvenet,
        #                              hull2.lcurvenet,
        #                              view='y-vertical')
        #        curve.plot3DmultiList(hull2.tcurvenet,
        #                              hull2.lcurvenet,
        #                              view='z-vertical')
        return
    
    
    
    
    def plot_bulbous_bow_form(self, ans_chr=None):
        bbow = self.bbow
        bbow.mb_ideal.Lspline.curve.plotcurve()
        bbow.bwl_ideal.Lspline.curve.plotcurve()
        bbow.cp_ideal.Lspline.curve.plotcurve()
        bbow.hullsurf.plotSurface(limx=10.,limy=10.,limz=10.)
        return
    
    def plot_bulbous_bow(self, ans_chr=None):
        bbow = self.bbow
        bbow.tcurvenet[0].plot3DmultiList(bbow.tcurvenet,
                          bbow.lcurvenet,
                          limx=10.,limy=10.,limz=10.)
        return
    
    def plot_get_extremes(self):
        hull = self.THBhull
        #
        biggest = max(hull.get_extremes())
        xb = (-biggest*.0,1.*biggest)
        yb = (-biggest*.5,.5*biggest)
        zb = (-biggest*.5,.5*biggest)
        return xb,yb,zb
        
    def plot_THB_ship_hull(self):
        """note that the transverses in this plot
        are not the hull form transverses
        but instead the u compliment to
        the longitudinal's v, so to speak.
        """
        hull = self.THBhull
        hull.compute_surface()
        #
        xb,yb,zb = self.plot_get_extremes()
        #
        hull.plotSurface()
        #
        curve = hull.hullcurvelist[0]
        xb,yb,zb = self.plot_get_extremes()
        curve.plot3DmultiList(hull.child.child.tcurvelist,
                              hull.child.child.lcurvelist,
                              view='x-vertical',
                              minx=xb[0],
                              miny=yb[0],
                              minz=zb[0],
                              limx=xb[1],
                              limy=yb[1],
                              limz=zb[1])
        
        curve.plot3DmultiList(hull.hullcurvelist,
                              hull.lcurvelist,
                              view='x-vertical',
                              minx=xb[0],
                              miny=yb[0],
                              minz=zb[0],
                              limx=xb[1],
                              limy=yb[1],
                              limz=zb[1])
        
        
        
        #        curve = hull.child.child.tcurvelist[0]
        #        curve.plot3DmultiList(hull.child.tcurvelist,
        #                              hull.child.child.lcurvelist,
        #                              minx=-100.,
        #                             miny=-100.,
        #                             minz=-100.,
        #                             limx=100.,
        #                             limy=100.,
        #                             limz=100.)
        return
    
    def plot_fine_hull_curves_3view(self,resolution_level='high'):
        """3 plots, testing the ability to show 
        plots amenable to viewing various aspects of the hull in 3d,
        since native matplotlib is 'hard to drive around in 3D'
        
        EDIT: the old name of this funciton was quite misleading.
        These are NOT THB curves.  These are just 
        the B-spline surface curves at the highest level of detail!
        
        inputs
        ----------
        resolution_level : high = finest level of curve detail
                            med = intermediate level of curve detail
                            low = lowest level of curve detail
                            
        
        
    
    
        
        note
        ----------
        
        -The transverses in this plot
        are not the hull form transverses
        but instead the u compliment to
        the longitudinal's v, so to speak.
        
        """
        resolution = {'low':0,
                      'med':1,
                      'high':2}
        level = resolution[resolution_level]
        
        hull = self.THBhull
        hull.compute_surface()
        
        curve = hull.hullcurvelist[0]
        
        this = hull.get_surface_by_level(level)
        
        xb,yb,zb = self.plot_get_extremes()
        curve.plot3DmultiList(this.tcurvelist,
                              this.lcurvelist,
                              view='x-vertical',
                              minx=xb[0],
                              miny=yb[0],
                              minz=zb[0],
                              limx=xb[1],
                              limy=yb[1],
                              limz=zb[1])
        #        curve.plot3DmultiList(this.tcurvelist,
        #                              this.lcurvelist,
        #                              view='y-vertical')
        #        curve.plot3DmultiList(this.tcurvelist,
        #                              this.lcurvelist,
        #                              view='z-vertical')
        return
    
    def get_template(self):
        # lines = testlines
        pre = '/home/luke/Documents/computational_naval_architecture/projects/relational_hull_design'
        directory = pre+'/relational_lsplines/relational_lsplines'
        filename =   'TemplateTable.tex'
        testlines = FileTools.GetLines(directory,filename)
        for line in testlines:
            print line
        return testlines
    
    def export_latex_results_table(self):
        """
        ship designer to latex
        
        """
        def roundthis(Q,p=3):
            if isinstance(Q,float):
                return round(Q,p)
            return Q.roundedreturn(p)
        #https://gist.github.com/jackiekazil/6201722
        #from decimal import getcontext, Decimal
        getcontext().prec = 3
        #
        self.hull.get_results()
        self.bbow.get_results()
        ## 
        ##################################################################
        ## map kanren to latex
        tl = self.design_space.get_node_list()
        k2l = {}
        #l2k = {}
        for el in tl:
            k2l[el] = el.name
            #l2k[el.name] = None
        l2k={'Acp'  : r'\acp',#done in code
             'Amsh' : r'\am', #done in code
             'Awp'  : r'\awp',#done in code
             'Cb'   : r'\cb',
             #'Ccp'  : r'\ccp',
             'Ccp'  : r'\CPK',
             'Cdl'  : r'\dlr', #Disp to Length Ratio, approx page 93, BasicRules.tex 
             'Clb'  : r'\Clb',
             'Clcg' : r'\Clcg',
             'Cmidshp'  : r'\cmdshp',
             'Cp'   : r'\cp',
             #'Cwp'  : r'\cwp',
             'Cwp'  : r'\CWL',
             'LCG'  : r'\lcb',
             'bwl'  : r'\beam',
             'disp' : r'\volDisp',
             'draft': r'\draft',
             'lfcp' : r'\lenFCPK',
             'lfsac': r'\lfsac',
             'lfwl' : r'\lenFWL', #\lenFSAC
             'lwl'  : r'\lwl',
             'vol'  : r'\volDisp'} #done in code
        il2k = {}
        for key in l2k:
            v = l2k[key]
            il2k[v] = key
            
        
        b2k = {'CBulbMid'       : r'\cbulbmdshp',
               'CBulbWtrPln'    : r'\cwpbulb',
               'CBulbCtrPln'    : r'\ccpbulb',
               'CBulbBlock'     : r'\cbbulb',
               'CBulbPrismatic' : r'\cpbulb',
               'Cbb'            : r'\Cbb',
               'Clpr'           : r'\Clpr',
               'Czb'            : r'\Czb',
               'Cabt'           : r'\Cabt',
               'Cabl'           : r'\Cabl',
               'Cvpr'           : r'\Cvpr',
               'A_lateral'      : r'\alateralbulb',
               'A_mid'          : r'\ambulb',
               'A_BBwl'         : r'\awpbulb',
               'A_lateral'      : r'\acpbulb',
               'BulbVolume'     : r'\volDispbulb',
               'BulbBeam'       : r'\beambulb',
               'BulbDepth'      : r'\draftbulb',
               'BulbLength'     : r'\lwlbulb'}
        ib2k = {}
        for key in b2k:
            v = b2k[key]
            ib2k[v] = key
        
        
        
        #for el in k2l:
        #    var, val = SD.design_space.rgp.get_name_and_value(el)
        #    print l2k[var.name]
            
        lines = self.get_template()
        #
        #  line = lines[25]
        #
        for i,line in enumerate(lines):
            #
            tokens = line.split()
            #
            if not tokens:
                continue #skips to top of loop
            ident = tokens[0]
            if '$' in ident:
                keyv = ident[1:-1]
                
                latexdict   = il2k
                latex_body  = self.hull
                latex_design  = self.design_space
                #if True:
                try:
                    var, val = latex_design.rgp.get_name_and_value(
                                                        latexdict[keyv])
                    actual = latex_body.results_dict[latexdict[keyv]]
                    print 'HEY ', actual
                    if actual is not None:
                        actual = latex_body.get_value(actual)
                        difference = 100.*(actual-val[0].midpoint())/actual
                        print 'HEY ', difference
                    else:
                        difference = None
                    d = dict(val = roundthis(val[0],3),
                             act = roundthis(actual,3),
                             dif = roundthis(difference,3) )
                    lines[i] = line.format(**d)
                    print line
                except:
                    print 'pass'
                
                
                
                latexdict   = ib2k
                latex_body  = self.bbow
                latex_design  = self.bbow
                #if True:
                try:
                    var, val = latex_design.rgp.get_name_and_value(
                                                        latexdict[keyv])
                    actual = latex_body.results_dict[latexdict[keyv]]
                    print 'HEY ', actual
                    if actual is not None:
                        #actual = latex_body.get_value(actual)
                        difference = 100.*(actual-val[0].midpoint())/actual
                        print 'HEY ', difference
                    else:
                        difference = None
                    d = dict(val = roundthis(val[0],3),
                             act = roundthis(actual,3),
                             dif = roundthis(difference,3) )
                    lines[i] = line.format(**d)
                    print line
                except:
                    pass
            else:
                print line
                
        
        the_filename = 'barehull_table.tex'
        with open(the_filename, 'w') as f:
            for line in lines:
                f.write(line)
        return
    
    #
    #**************************************************************************
    # Quick Demo of random design
    def demorun(self, verbosity=False,dmy=None):
        print '\n--------------------------------------------------------------'
        print 'Start Generation of Hull with Bulbous Bow '
        print 'representation type: THB-spline'
        print '--------------------------------------------------------------\n'
        self.verbose = verbosity
        self.bare_hull_random_design_selection_new()
        self.make_bare_hull()
        self.bulbous_bow_random_design_selection()
        self.make_bulbous_bow()
        #self.make_THBspline_complex_hull()
        print '\n--------------------------------------------------------------'
        print 'DONE'
        print '--------------------------------------------------------------\n'
        return
    
    
    #
    #**************************************************************************
    # Importing
    def Import(self, option=None, dmy=None):
        
        
        #
        #*********************************************************************
        # 
        def printoptions():
            """help menu for the import functions
            of the class ShipDesigner
            """
            message = '\nAt present bare hull design space import is '+\
                                            'implemented: option=designspace\n'
            message2 = '\n\n'
            print 'Welcome to the import menu'
            print ''
            print 'Thank you for using this hull generation program!'
            print ''
            print message
            print ''
            #print message2
            print ''
            print 'Availible options are:'
            for key in options:
                print ' ',key
            return
        
        
        def importdesignspace():
            """assumes the design space is extant
            (will only work if feasible with current design space)
            
            let user decide to search (it never hurts to search a thin designspace)
            """
            if dmy is not None:
                file_Name1 = dmy
            else:
                file_Name1 = "designspace"
                
            print 'design space import: ', file_Name1
            fileObject = open(file_Name1,'r') 
            ds=pickle.load(fileObject)  
            fileObject.close()
            
            
            print '-----------------------------------------------------'
            print 'starting import of design-space file name = ',file_Name1
            print 'building design database from file...'
            for name in ds:
                var = getattr(self.design_space,name)
                if self.verbose:
                    if isinstance(var, lp.Variable):
                        print 'var ', var,' is ',lp.Variable
                    elif isinstance(var, lp.PStates):
                        print 'var ', var,' is ',lp.PStates
                    else:
                        print 'var ', var,' is ',type(var)
                value = ds[name]
                mkvar, mkval = self.design_space.get_from_rgp(name)
                # set new value as this rule:
                var = var == value[0]
                self.design_space.rgp.add_one_rule(var,var.name)
                self.design_space.rgp.compute_fresh_rules_graph()
                print name,' => ',value[0]
            
            print '-----------------------------------------------------'
            print 'successful import finished'
            self.thin_hull_design = hullmethods.Hull(
                                        States = self.design_space.rgp.env, 
                                        which_state=0,
                                        revise=False)
            
            
            for name in self.design_space.__dict__:
                dname = getattr(self.design_space, name)
                if isinstance(dname, lp.PStates):
                    if self.verbose: print dname
                    mkname,env_value = self.design_space.get_from_rgp(dname)
                    if self.verbose: print 'setting PStates',dname.name,' => mkVar ',mkname,' : ',env_value,'\n\n\n'
                    setattr(self.thin_hull_design.hullconstraints, 
                            name,
                            mkname)
            self.bareblockextremes = self.get_bare_hull_extremes()
            print '-----------------------------------------------------'
            print 'Design Space selection Complete'
            return
        
        
        
        def importdesignspaceSimple():
            """assumes the design space is extant
            (will only work if feasible with current design space)
            
            let user decide to search (it never hurts to search a thin designspace)
            """
            if dmy is not None:
                file_Name1 = dmy
            else:
                file_Name1 = "designspace"
                
            print 'design space import: ', file_Name1
            fileObject = open(file_Name1,'r') 
            ds=pickle.load(fileObject)  
            fileObject.close()
            
            
            print '-----------------------------------------------------'
            print 'starting import of design-space file name = ',file_Name1
            print 'building design database from file...'
            for name in ds:
                var = getattr(self.design_space,name)
                if self.verbose:
                    if isinstance(var, lp.Variable):
                        print 'var ', var,' is ',lp.Variable
                    elif isinstance(var, lp.PStates):
                        print 'var ', var,' is ',lp.PStates
                    else:
                        print 'var ', var,' is ',type(var)
                value = ds[name]
                mkvar, mkval = self.design_space.get_from_rgp(name)
                # set new value as this rule:
                var = var == value[0]
                self.design_space.rgp.add_one_rule(var,var.name)
                self.design_space.rgp.compute_fresh_rules_graph()
                print name,' => ',value[0]
            
            print '-----------------------------------------------------'
            print 'successful import finished'
            self.thin_hull_design = hullmethods.Hull(
                                        States = self.design_space.rgp.env, 
                                        which_state=0,
                                        revise=False)
            
            
            for name in self.design_space.__dict__:
                dname = getattr(self.design_space, name)
                if isinstance(dname, lp.PStates):
                    if self.verbose: print dname
                    mkname,env_value = self.design_space.get_from_rgp(dname)
                    if self.verbose: print 'setting PStates',dname.name,' => mkVar ',mkname,' : ',env_value,'\n\n\n'
                    setattr(self.thin_hull_design.hullconstraints, 
                            name,
                            mkname)
            self.bareblockextremes = self.get_bare_hull_extremes()
            print '-----------------------------------------------------'
            print 'Design Space selection Complete'
            return
        
        #
        #*********************************************************************
        # 
        options = {'designspace':importdesignspace,
                   'help':printoptions,
                   'option':printoptions}
        
        if option is None:
            option = 'help'
        
        options[option.lower()]()
        return 
    
    #
    #**************************************************************************
    # Exporting
    
    def export(self,option=None,dmy=None):#,*kwds):
        """At present only bare hull curve exports are
        supported.
        """
        #key = kwds
        #
        #*********************************************************************
        # 
        def printoptions():
            """help menu for the export functions
            of the class ShipDesigner
            """
            #if key:
            #    print ' help : ',key
            #    print help(options[key[0]])
            #else:
            message = '\nAt present bare hull curve exports are '+\
                                            'implemented: option=rhino\n'
            message2 = '\nAnd THB surface exporting is implemented '+\
                        'as option=pickle\n'
            print 'Welcome to the export menu'
            print ''
            print 'Thank you for using this hull generation program!'
            print ''
            print message
            print ''
            print message2
            print ''
            print 'Availible options are:'
            for key in options:
                print ' ',key
            return
        
        #
        #*********************************************************************
        # 
        def exportdesignspace():
            """
            Parameters
            ----------
            Returns
            ----------
            Notes
            ----------
            
            dev
            ----------
            self = SD
            import cPickle as pickle
            
            """
            import cPickle as pickle#hack attempt
            quick_design = {}
            for var in self.design_space.rgp.vars:
                try:
                    attr = self.design_space.__getattribute__(var.name)
                    mk_name, mk_val = self.design_space.get_from_rgp(attr)
                    quick_design[mk_name.name] = mk_val
                except:
                    pass
            file_Name1 = "designspace"
            fileObject = open(file_Name1,'wb') 
            pickle.dump(quick_design,fileObject)  
            fileObject.close()
            return
        
        #
        #*********************************************************************
        # 
        def exportrhino():
            """export transverse and longitudinal curves
            to rhino
            """
            self.hull.export_rhino()
            self.bbow.export_curves(tcurvenet = self.THBbow.tcurvelist,
                                    lcurvenet = self.THBbow.lcurvelist)
            return
        
        def exportrhino_surface(filename=None):
            if filename is None:
                self.THBhull.IssueRhinoSurface(the_filename='SingleHull.txt')
            else:
                self.THBhull.IssueRhinoSurface(the_filename=filename)
                
            return
        
#        def exportLatex_table():
#            if dmy is not None:
#                file_Name1 = dmy
#            else:
#                if isinstance(dmy, dict):
#                    file_Name1 = dmy["fname"]
#                else:
#                    file_Name1 = "designspace"
#                    
#            flines = []
#            flines.append('\begin{small} \n')
#            flines.append('	\begin{table}  \n')
#            flines.append('		\begin{tabular}{llcccc} \n')
#            
#            return
        
        
        def pickle():
            """
            Parameters
            ----------
            Returns
            ----------
            Notes
            ----------
            
            dev
            ----------
            import cPickle as pickle
            """
            file_Name1 = "barehullsurf"
            fileObject = open(file_Name1,'wb') 
            pickle.dump(self.THBhull,fileObject)  
            fileObject.close()
            
            
            file_Name2 = "bbowsurf"
            fileObject = open(file_Name2,'wb') 
            pickle.dump(self.THBbow,fileObject)  
            fileObject.close()
            return
        
        
        #
        #*********************************************************************
        # 
        options = {'rhino':exportrhino,
                   'designspace':exportdesignspace,
                   #'designspaceSimple':importdesignspaceSimple,
                   'latex':self.export_latex_results_table,
                   'pickle':pickle,
                   'help':printoptions,
                   'option':printoptions}
        
        if option is None:
            option = 'help'
        
        options[option.lower()]()
        #
        #*********************************************************************
        #
        return
    
    
    
    #
    #**************************************************************************
    # Post Process
    def postprocess_bare_hull_stats(self, option=None,
            detailed=False,dmy=None):
        return self.print_bare_hull_stats(option,detailed)
    
    
    def print_bare_hull_stats(self, option=None,
            detailed=False,dmy=None):
        """
        
        Parameters
        ----------
            option : availible options can be queried by 
                    running this command with the option help, 
                    or with no command at all
        Returns
        ----------
            Nothing is returned.  This function prints a 
            variety of post processed hull form information
            
        Notes
        ----------
            currently supporting:
                Printing rules design parameters and corresponding 
                bare hull results 
        
        TODO: 
        ----------
            make the Database able to get from the Lsplines?
        """
        #
        #*********************************************************************
        # 
        def printoptions():
            print 'Welcome to the postprocessing menu'
            print ''
            print 'Thank you for using this hull generation program!'
            print ''
            print 'Availible options are:'
            for key in options:
                print ' ',key
            return
        #
        #*********************************************************************
        # 
        def printflats():
            print '--------------------------------------------------'
            print 'Postprocessing: Checking hull curve flat specifications '
            print '                 against the thin design'
            print '                 determined by the rules'
            print '--------------------------------------------------'
            print 'bare hull SAC'
            print '   flat portion:'
            print '      start : ',self.hull.stations.FOSAC[0]
            print '      stop  : ',self.hull.stations.FOSAC[1]
            print '     length : {}'.format(self.hull.stations.FOSAC[1] - 
                                             self.hull.stations.FOSAC[0])
            var,val = self.design_space.get_from_rgp('lfsac')
            print '   Rules Constraint:'
            print '      {} value: {}'.format(var,val)
            print '--------------------------------------------------'
            print 'bare hull WaterLine'
            print '   flat portion:'
            print '      start : ',self.hull.stations.FOWL[0]
            print '      stop  : ',self.hull.stations.FOWL[1]
            print '      diff  : {}'.format(self.hull.stations.FOWL[1] - 
                                             self.hull.stations.FOWL[0])
            var,val = self.design_space.get_from_rgp('lfwl')
            print '   Rules Constraint:'
            print '      {} value: {}'.format(var,val)
            print '--------------------------------------------------'
            print 'bare hull CenterPlane'
            print '   flat portion:'
            print '      start : ',self.hull.stations.FOCP[0]
            print '      stop  : ',self.hull.stations.FOCP[1]
            print '      diff  : {}'.format(self.hull.stations.FOCP[1] - 
                                             self.hull.stations.FOCP[0])
            var,val = self.design_space.get_from_rgp('lfcp')
            print '   Rules Constraint:'
            print '      {} value: {}'.format(var,val)
            print '--------------------------------------------------'
            return
        #
        #*********************************************************************
        # 
        def printSAC_to_hull():
            print '--------------------------------------------------'
            print 'Postprocessing: Checking SAC characteristics'
            print '                 against the centerplane'
            print '                 and waterplane design'
            print '                 in rule and form parameter terms'
            print '--------------------------------------------------'
            disp,disp_val = self.design_space.get_from_rgp('disp')
            lwl,lwl_val = self.design_space.get_from_rgp('lwl')
            draft,draft_val = self.design_space.get_from_rgp('draft')
            bwl,bwl_val = self.design_space.get_from_rgp('bwl')
            lfsac,lfsac_val = self.design_space.get_from_rgp('lfsac')
            Awp,Awp_val = self.design_space.get_from_rgp('Awp')
            Acp,Acp_val = self.design_space.get_from_rgp('Acp')
            Cp,Cp_val = self.design_space.get_from_rgp('Cp')
            Cb,Cb_val = self.design_space.get_from_rgp('Cb')
            Ccp,Ccp_val = self.design_space.get_from_rgp('Ccp')
            Cwp,Cwp_val = self.design_space.get_from_rgp('Cwp')
            return
        #
        #*********************************************************************
        # 
        def printSAC():
            """
            """
            print '--------------------------------------------------'
            print 'bare hull SAC'
            Lspline = self.hull.Lsplines.SAC
            Lspline.Lagrangian.print_form_parameters()
            Lspline.Lagrangian.print_objective_functions()
            if detailed:
                print 'Gradient of Active Vertices and Multipliers at convergence:'
                Lspline.f.grad.T[Lspline.mask.T]
            print '--------------------------------------------------'
            print 'postprocessing...'
            print 'bare hull SAC: summary of Lspline Optimization issues'
            ppf = self.postprocessor(Lspline)
            print ppf.badparameters
            #
            #******************************************************************
            # Database Vs Lspline
            print '--------------------------------------------------'
            print 'postprocessing...  '
            print '             Checking SAC specifications '
            print '                 against the thin design'
            print '                 determined by the rules \n'
            #
            #******************************************************************
            # area
            print '--------------------------------------------------'
            print 'bare hull SAC'
            var,val = self.design_space.get_from_rgp('vol')
            print '   Dispacement'
            print '       Rules Constraint:'
            print '      {} value: {}\n'.format(var,val)
            sac_area = self.hull.Lsplines.SAC.curve.area #many ways to get this!
            print '   Displacement from optimization'
            print '                       of SAC itself: {}\n'.format(sac_area)
            print ''
            print '   diff : {}'.format(val[0]-sac_area.value)
            print '--------------------------------------------------'
            #
            #******************************************************************
            # Xc
            print '--------------------------------------------------'
            print 'bare hull Sectional Area Curve (SAC)'
            var,val = self.design_space.get_from_rgp('LCG')
            print '   LCG '
            print '      Rules Constraint:'
            print '      {} value: {}\n'.format(var,val)
            sac_Xc = self.hull.Lsplines.SAC.curve.Xc #many ways to get this!
            print '   Xc from optimization'
            print '                       of SAC itself: {}\n'.format(sac_Xc)
            print ''
            print '   diff : {}'.format(val[0]-sac_Xc.value)
            print '--------------------------------------------------'
            #
            #******************************************************************
            # flats, see flats
            distance = self.hull.stations.FOSAC[1] - \
                                             self.hull.stations.FOSAC[0]
            print '--------------------------------------------------'
            print 'bare hull SAC'
            print '   flat portion'
            print '   from optimization:'
            print '      start : ',self.hull.stations.FOSAC[0]
            print '      stop  : ',self.hull.stations.FOSAC[1]
            print '      diff  : {}'.format(distance)
            var,val = self.design_space.get_from_rgp('lfsac')
            print '   Rules Constraint:'
            print '      {} value: {}'.format(var,val)
            print ''
            print '   diff : {}'.format(val[0]-distance)
            print '--------------------------------------------------'
            return
        #
        #*********************************************************************
        # 
        def printDWL():
            print '--------------------------------------------------'
            print 'bare hull WaterLine (DWL)'
            Lspline = self.hull.Lsplines.SAC
            Lspline.Lagrangian.print_form_parameters()
            Lspline.Lagrangian.print_objective_functions()
            if detailed:
                print 'Gradient of Active Vertices and Multipliers at convergence:'
                Lspline.f.grad.T[Lspline.mask.T]
            print '--------------------------------------------------'
            print 'postprocessing...'
            print 'bare hull DWL: summary of Lspline Optimization issues'
            ppf = self.postprocessor(Lspline)
            print ppf.badparameters
            #
            #******************************************************************
            # Database Vs Lspline
            print '--------------------------------------------------'
            print 'Postprocessing: '
            print '             Checking DWL specifications '
            print '                 against the thin design'
            print '                 determined by the rules'
            print '--------------------------------------------------'
            #
            #******************************************************************
            # area
            print '--------------------------------------------------'
            print 'bare hull DWL'
            var,val = self.design_space.get_from_rgp('Awp')
            print '   Waterplane Area (Awp)'
            print '       Rules Constraint:'
            print '       {} value: {}\n'.format(var,val)
            print '       Rules Constraint:'
            print '       Used as Half Waterplane area in design:'
            print '       {} half value: {}\n'.format(var,val[0]*.5)
            sac_area = self.hull.Lsplines.DWL.curve.area #many ways to get this!
            print '       Awp from optimization'
            print '                of the  DWL itself: {}'.format(sac_area)
            print ''
            print '   diff : {}'.format(val[0]*.5-sac_area.value)
            print '--------------------------------------------------'
            return
        #
        #*********************************************************************
        # 
        def printCPK():
            print '--------------------------------------------------'
            print 'bare hull CenterProfile Keel'
            Lspline = self.hull.Lsplines.SAC
            Lspline.Lagrangian.print_form_parameters()
            Lspline.Lagrangian.print_objective_functions()
            if detailed:
                print 'Gradient of Active Vertices and Multipliers at convergence:'
                Lspline.f.grad.T[Lspline.mask.T]
            print '--------------------------------------------------'
            print 'postprocessing...'
            print 'bare hull CProfile Keel: summary of Lspline Optimization issues'
            ppf = self.postprocessor(Lspline)
            print ppf.badparameters
            #
            #******************************************************************
            # Database Vs Lspline
            print '--------------------------------------------------'
            print 'Postprocessing: '
            print '             Checking CProfile specifications '
            print '                 against the thin design'
            print '                 determined by the rules'
            print '--------------------------------------------------'
            #
            #******************************************************************
            # area
            print '--------------------------------------------------'
            print 'bare hull CProfile Curve (CPK)'
            var,val = self.design_space.get_from_rgp('Acp')
            print '   Centerplane Area (Acp)'
            print '       Rules Constraint:'
            print '       {} value: {}\n'.format(var,val)
            sac_area = self.hull.Lsplines.CProfile.curve.area #many ways to get this!
            print '       Acp from optimization'
            print '                of the  CProfile itself: {}'.format(sac_area)
            print ''
            print '   diff : {}'.format(val[0]-sac_area.value)
            print '--------------------------------------------------'
            return
        #
        #*********************************************************************
        # 
        def printLongitudinals():
            lcurvenet = self.hull.Lsplines.lcurves
            print ''
            for i,Lspline in enumerate(lcurvenet):
                print '--------------------------------------------------'
                print '  Post Processing Interior Longitudinal Lspline ',i
                print '--------------------------------------------------'
                ppf = self.postprocessor(Lspline)
                print ppf.badparameters
                print '--------------------------------------------------'
                print 'method two'
                aeq, amax, amin = Lspline.Lagrangian.get_anomolous_form_parameters()
                print 'Bad Equality Constraints:'
                print aeq
                print 'Bad Max Constraints:'
                print amax
                print 'Bad Min Constraints:'
                print amin
                print ''
                
            return
        #
        #*********************************************************************
        # 
        options = {'flats':printflats,
                   'sac':printSAC,
                   'dwl':printDWL,
                   'cpk':printCPK,
                   'Longitudinals':printLongitudinals,
                   'help':printoptions,
                   'option':printoptions}
        if option is None:
            option = 'help'
        options[option]()
        return
    
            
    
    def store(self,dmy=None):
        """
        Parameters
        ----------
        Returns
        ----------
        Notes
        ----------
        """
        print "please use self.export('pickle') in the future"
        file_Name1 = "barehullsurf"
        fileObject = open(file_Name1,'wb') 
        pickle.dump(self.THBhull,fileObject)  
        fileObject.close()
        
        
        file_Name2 = "bbowsurf"
        fileObject = open(file_Name2,'wb') 
        pickle.dump(self.THBbow,fileObject)  
        fileObject.close()
        
        
        return

#
#******************************************************************************
# END

if __name__ == '__main__':
    SD = ShipDesigner()
    #    print 'start ShipDesigner?'
    #    sinput = raw_input('enter n to avoid it>')
    #    if sinput is 'n':
    #        pass
    #    else:
    #        SD.ShipDesigner()
    SD.ShipDesigner()
if False:
    import simple_hull_rules_language 
    try:
        from pycallgraph import PyCallGraph
        from pycallgraph.output import GraphvizOutput
        graphviz = GraphvizOutput()
        graphviz.output_file = 'TLMcode.png'
        pycall_availible = True
    except:
        print 'no module pycallgraph, please profile some other way!'
        pycall_availible = False
        
    
    
    # True  # False
    dofull      = True  #make bare hull and bbow
    pickleit    = True   #store it
    profileit   = False
    #check = False
    """
        Design spec:
        Ccp -> fraction of depth*lwl which makes CPkeel area
        Cwp -> fraction of bwl*lwl which makes DWL area
        
        When Cb is high and Cwp/Ccp are low, DWL/CPKeel end up wavey
    """
    if dofull:
        LITTLE = 1.e-2
        hull_template = hullclp()# hull_use_HICLP.hullclp()
        
        
        spec = DesignSpecification(
                                lwl = (100.,130.),
                                draft = (15.,25.),
                                bwl = (25.,40.),
                                #vol=(18000.,20000.),
                                #vol=(22000.,28000.),
                                vol=(12000.,25000.),
                                LCG = (49.,66.),
                                Clcg = (.485,.515),
                                Cb = (0.5,0.89),
                                Cwp = (.85,.99),#water plane coefficient
                                Cmidshp = (0.7,.95),
                                Ccp = (.85,.99) #centerplane coefficient
                                )
                                
                                #Cb = (0.7,0.95),
                                #Cwp = (.7,.95),#water plane coefficient
                                #Cmidshp = (0.9,0.99),
                                #Ccp = (.7,0.95))#centerplane coefficient
        """Construct spec sets
            to do clustering
            
            eg 
            -High WP area low draft High vol fwd Xc vessels
            vs
            -Low WP area deep draft High vol fwd Xc vessesl
            
            etc...
            
        """
        """Add rules to relate vol
            with 
            Cwp*depth
            Ccp*bwl
                
            fix longitudinal maker bug
            
            straighten CP and DWL (fairness func??)
            
            how about a trapazoidal area rule
            relating BFC/SFC to Midship Curves
            -to try and keep the DWL and CPKeel
            inside the bounds, non-doubly re-entrant
        
        """
        """
            self = hull2.sac_co
            
            idea: leave the fore and aft fairness 
                curves -floating-
            
                find the right spot on he SAC curve
                and stick them there at build time.
            
            issue: that might mess up the 
                DWL and CPkeel
                
            idea:
                design SAC, DWL, CPkeel using CLP
                -then infer what BFC,SFC are allowed
                
            current practice:
                -SAC, BFC,SFC are exact
                -imposing implicit area requrements
                on DWL and CPkeel
                -these need explicit relations
                -or else to be eliminated
        #"""
        tp = TestParameters(spec,N=1)
            
        ds = DesignSpace(spec)
        sobol_seq=None
        initial_state = None
        Nd = 1
        SMALL = DesignSpace.SMALL
        i=0
        w=.1
        
        self = ds
        
        self.tree_search()
        self.print_widths()
        #from hull_from_designspace import Hull as design_builder
        #"""
        ok = False
        cnt= -1
        while True:
            cnt += 1
            print cnt
            ok = self.check_widths()
            if ok:
                hull1 = hullmethods.Hull(States = self.hdp, 
                                         which_state=0,
                                         revise=False)
                break
            else:
                #            print 'no good'
                #            break
                try:
                    self.tree_search()
                except:
                    print 'error, climb the tree'
                    kparent = self.design_sq.get_latest_child().parent
                    self.hdp.states = kparent.states
            if cnt>5:
                print 'tree count > 5'
                break
        #"""
        if ok:
            hull2 = give_results(hull1)
            plot_results(hull2)
            
            canvas = hullmethods.hull_gui(hull2)
            
            #"""
            bbow = interface_to_bbow(self.hdp)
            
            bbow.mb_ideal.Lspline.curve.plotcurve()
            bbow.bwl_ideal.Lspline.curve.plotcurve()
            bbow.cp_ideal.Lspline.curve.plotcurve()
            #"""
            
            
            #hull2.make_bulb(bbow)
            #bbow.specialized_bow_development_nose()
            #bbow.make_bulb()
            #curve = bbow.tcurvenet[0]
            bbow.tcurvenet[0].plot3DmultiList(bbow.tcurvenet,
                          bbow.lcurvenet,
                          limx=10.,limy=10.,limz=10.)
            
            hull2.hullsurf.plotSurface(limx=10.,limy=10.,limz=10.)
            bbow.hullsurf.plotSurface(limx=10.,limy=10.,limz=10.)
            
            
            hull2.hullsurf.plotMultiSurface(
                            [hull2.hullsurf,bbow.hullsurf],
                            limx=100.,limy=100.,limz=100.)
            
            
            b2 = copy.deepcopy(bbow.hullsurf)
            #b2 = b2.rotate_surface( 
            #                (0.,0.,1.), np.pi/2.)
            b2 = b2.reflect_surface(np.asarray([0.,0.,1.]),
                                    np.asarray([0.,0.,1.]))
            b2.compute_surface()
            
            b2.plotMultiSurface([bbow.hullsurf,b2],
                                limx=100.,limy=100.,limz=100.)
            #            
            #            
            #            bbow.hullsurf.plotMultiSurface(
            #                            [hull2.hullsurf,bbow.hullsurf],
            #                            limx=100.,limy=100.,limz=100.)
            #            
            #            
            #            b2.plotMultiSurface(
            #                            [hull2.hullsurf,b2],
            #                            limx=100.,limy=100.,limz=100.)
            #            
            #            
            #            b3 = copy.deepcopy(bbow.hullsurf)
            #            b3 = bbow.hullsurf.rotate_surface( 
            #                            (0.,1.,0.), np.pi/2.)
            #            b3.compute_surface()
            #            
            #            b3.plotMultiSurface(
            #                            [hull2.hullsurf,b3],
            #                            limx=100.,limy=100.,limz=100.)
            
            
                
            ##
            ##**********************************************
            ## 
            if pickleit:
                    
                file_Name3 = "rm_11x7"
                file_Name4 = "rm_19x11"
                file_Name5 = "rm_35x19"
                file_Name6 = "rm_67x35"
                file_Name7 = "rm_7x5"
                
                ## 7x5
                fileObject = open(file_Name7,'r')  
                rm_7x5 = pickle.load(fileObject) 
                fileObject.close()
                
                ## 11x7
                fileObject = open(file_Name3,'r')  
                rm_11x7 = pickle.load(fileObject) 
                fileObject.close()
                
                ## 19x11
                fileObject = open(file_Name4,'r')  
                rm_19x11 = pickle.load(fileObject) 
                fileObject.close()
        
                ## 35x19
                fileObject = open(file_Name5,'r')  
                rm_35x19 = pickle.load(fileObject) 
                fileObject.close()
                
                ## 67x35
                fileObject = open(file_Name6,'r')  
                rm_67x35 = pickle.load(fileObject) 
                fileObject.close()
                
                t1 = thbspline.THBsurface(vertices = hull2.hullsurf.surfNet,
                                          tcurvelist = hull2.tcurvenet,
                                          ulevel = 2,
                                          lcurvelist = hull2.lcurvenet,
                                          vlevel = 3,
                                          trm = rm_11x7,
                                          lrm = rm_19x11)
            ##
            ##
                t2 = thbspline.THBsurface(vertices = bbow.hullsurf.surfNet,
                                          tcurvelist = bbow.tcurvenet,
                                          ulevel = 2,
                                          lcurvelist = bbow.lcurvenet,
                                          vlevel = 1,
                                          trm = rm_11x7,
                                          lrm = rm_7x5)
            ##
            ##
                file_Name1 = "barehullsurf"
                fileObject = open(file_Name1,'wb') 
                pickle.dump(t1,fileObject)  
                fileObject.close()
                
                file_Name2 = "bbowsurf"
                fileObject = open(file_Name2,'wb') 
                pickle.dump(t2,fileObject)  
                fileObject.close()
                
                ##
                ##**********************************************
                ## PROJECTION MATRICES (do once) 
                ##              - really, I mean it, just do this once ever
                if False:
                    file_Name3 = "rm_11x7"
                    fileObject = open(file_Name3,'wb') 
                    pickle.dump(t1.trm,fileObject)  
                    fileObject.close()
                    
                    file_Name4 = "rm_19x11"
                    fileObject = open(file_Name4,'wb') 
                    pickle.dump(t1.child.trm,fileObject)  
                    fileObject.close()
                    
                    file_Name5 = "rm_35x19"
                    fileObject = open(file_Name5,'wb') 
                    pickle.dump(t1.child.child.trm,fileObject)  
                    fileObject.close()
                    
                    file_Name6 = "rm_67x35"
                    fileObject = open(file_Name6,'wb') 
                    pickle.dump(t1.child.child.lrm,fileObject)  
                    fileObject.close()
                    
                    file_Name7 = "rm_7x5"
                    fileObject = open(file_Name7,'wb') 
                    pickle.dump(t2.lrm,fileObject)  
                    fileObject.close()
                
                
                #https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map
                #can't pickle objets with classmethods and such:
                #                file_Name3 = "bbowBsurf"
                #                fileObject = open(file_Name3,'wb') 
                #                pickle.dump(hull2.hullsurf,fileObject)  
                #                fileObject.close()
                
            #"""
                
            else: #don't pickle it
                ##
                ##**********************************************
                ## Don't pickle it
                
                t1 = thbspline.THBsurface(vertices = hull2.hullsurf.surfNet,
                                          hullcurvelist=hull2.tcurvenet,
                                          tcurvelist = hull2.tcurvenet,
                                          ulevel = 2,
                                          lcurvelist = hull2.lcurvenet,
                                          vlevel = 3,
                                          trm = rm_11x7,
                                          lrm = rm_19x11)
                
                
                
                t2 = thbspline.THBsurface(vertices = bbow.hullsurf.surfNet,
                                          hullcurvelist = bbow.tcurvenet,
                                          tcurvelist = bbow.tcurvenet,
                                          ulevel = 2,
                                          lcurvelist = bbow.lcurvenet,
                                          vlevel = 1,
                                          trm = rm_11x7,
                                          lrm = rm_7x5)
            ##
            ##**********************************************
            ##  end pickleit choices
            ##
            t1.set_top()
            t11 = t1.dyadic_loft_refinement(
                        bounds = thbspline.BoxList(  
                                thbspline.Box( ia(.0,.8125),      ##.5625
                                               ia(.0,.34375) ) ),  ##.3125
                                trm = rm_19x11,
                                lrm = rm_35x19
                                )
                        
            t12 = t11.dyadic_loft_refinement(bounds = thbspline.BoxList(  
                                thbspline.Box( ia(.0,.8125),      #8125  #.4375 ## 7-> 7+3+3=13
                                               ia(.0,.34375) ) ),  #.40625  #.3125    ## 5 -> 5+3+3 = 11
                                trm = rm_35x19,
                                lrm = rm_67x35
                                ) 
                        
            t1.compute_surface()
            t1.plotSurface(limx=100.,limy=100.,limz=100.)
            hull = t1 #for use in the end where hull and bbow are combined.
            
            
            
            
            t3 = t2.reflect_surface((0.,0.,1.),(0.,0.,1.))
            t3 = t3.reflect_surface((1.,0.,0.),(1.,0.,0.))
            t3 = t3.rotate_surface((0.,0.,1.),np.pi/4.)
            t3 = t3.translate_surface( dz=-np.amin(t3.vertices[:,:,2]) )
            t3.set_top()
            #            t22 = t2.dyadic_loft_refinement(bounds = thbspline.BoxList(  
            #                    thbspline.Box( ia(.125,.875), 
            #                                  ia(.125,.875) ) )    )
            #            
            #            t23 = t22.dyadic_loft_refinement(bounds = thbspline.BoxList(  
            #                                    thbspline.Box( ia(.38,.62), 
            #                                                      ia(.38,.62) ) )   ) 
            #            

            t3.compute_surface()
            t3.plotSurface(limx=10.,limy=10.,limz=10.)
            
            
            #"""
            
            #s0 = thbspline.THBsurface(verbose=True)
            #s0 = s0.dyadic_loft_refinement()
            #s0 = s0.dyadic_loft_refinement()
            #s0 = s0.dyadic_loft_refinement()
            #s0.set_top()
            #s1 = s0.dyadic_loft_refinement()
            
            #t1.eval(.3,.8)
            #t1.surface_point(.3,.8)
            
            #t2.eval(.3,.8)
            #t2.surface_point(.3,.8)
            
            """check the difference between B-spline and THBspline
            #suface evaluation for a given THBsurface called hull:
            #"""
            #"""
            if profileit:
                uopt.THBspline_to_Bspline_ckdiff_evaluator(t1)#,mod_level=2)
                uopt.THBspline_to_Bspline_ckdiff_evaluator(t3)
            #"""
                
        
        ##
        ##**********************************************
        ##  If we did not do a full bare-hull & bulb generation
        ## then get the geometry from pre-pickled files.
        else:
            
            def set_list(plist):
                for p in plist:
                    #x = self.sobol_seq.next()
                    x=.5
                    val = self.hdp(p)[0].getpoint(x)
                    print 'setting ', p,val
                    val = ia(val,val)
                    self.hdp.__setattr__(p.name,val)
                return
                
            def fa():
                ds = DesignSpace(spec)
                self = ds
                return self
            
            self = fa()
            set_list(self.hdp.hull_list)
            set_list(self.hdp.Coefficients)
            set_list(self.hdp.Areas)
        #"""
        print 'plot_results(hull2)'
        t1.set_top()
        av = t1.compute_true_nonzero_basis(.1,.1)
        tapb = t1.compute_total_active_by_level()
    ##
    ##*************************************************************************
    ##
    else: #we are starting 
        #from the existing pickled objects entirely
        #unpickle
        file_Name1 = "barehullsurf"
        file_Name2 = "bbowsurf"
        
        file_Name3 = "rm_11x7"
        file_Name4 = "rm_19x11"
        file_Name5 = "rm_35x19"
        file_Name6 = "rm_67x35"
        file_Name7 = "rm_7x5"
        
        ## BARE HULL
        fileObject = open(file_Name1,'r')  
        hull = pickle.load(fileObject) 
        fileObject.close()
        
        ## BBOW
        fileObject = open(file_Name2,'r')  
        t2 = pickle.load(fileObject)  
        fileObject.close()
        
        
        ## 7x5
        fileObject = open(file_Name7,'r')  
        rm_7x5 = pickle.load(fileObject) 
        fileObject.close()
        
        
        ## 11x7
        fileObject = open(file_Name3,'r')  
        rm_11x7 = pickle.load(fileObject) 
        fileObject.close()
        
        ## 19x11
        fileObject = open(file_Name4,'r')  
        rm_19x11 = pickle.load(fileObject) 
        fileObject.close()
        
        ## 35x19
        fileObject = open(file_Name5,'r')  
        rm_35x19 = pickle.load(fileObject) 
        fileObject.close()
        
        ## 67x35
        fileObject = open(file_Name6,'r')  
        rm_67x35 = pickle.load(fileObject) 
        fileObject.close()
        
        
        
        surfNet = hull.vertices
        tcurvenet = hull.tcurvelist
        lcurvenet = hull.lcurvelist
        handle = tcurvenet[0]
        
        
        ##
        ##*****************************************************
        ## profiling THBsurface 
        ##   -creation
        ##   -refinement
        ##   -evalutation
        #if profileit and pycall_availible:
        #    with PyCallGraph(output=graphviz):
            
        if True:
            #
            #**********************************************************
            # the BARE HULL
            hull  = thbspline.THBsurface(vertices = surfNet,
                                          tcurvelist = tcurvenet,
                                          ulevel = 2,
                                          lcurvelist = lcurvenet,
                                          vlevel = 3,
                                          hullcurvelist=hull.hullcurvelist,
                                          trm = rm_11x7,
                                          lrm = rm_19x11)
            hull.set_top()
            """
            import thbsurface 
            from myspline import rbspline
            #"""
            #
            #
            #"""######
            ##
            ##****************************************************************
            ## prepare to add bulbous bow to the BARE HULL:
            #11X19#
            t11 = hull.dyadic_loft_refinement(
                        bounds = thbspline.BoxList(  
                                thbspline.Box( ia(.0,.8125),      ##.5625
                                               ia(.0,.34375) ) ),  ##.3125
                                trm = rm_19x11,
                                lrm = rm_35x19
                                )
            #19X35#
            #"""
            #"""#####
            t12 = t11.dyadic_loft_refinement(bounds = thbspline.BoxList(  
                                thbspline.Box( ia(.0,.8125),      #8125  #.4375 ## 7-> 7+3+3=13
                                               ia(.0,.34375) ) ),  #.40625  #.3125    ## 5 -> 5+3+3 = 11
                                trm = rm_35x19,
                                lrm = rm_67x35
                                )
            #"""
            #35X35#
            #t13 = t12.dyadic_loft_refinement(bounds = thbspline.BoxList(  
            #                        thbspline.Box( ia(.0,.4), 
            #                                      ia(.0,.375) ) )    )
            #t14 = t13.dyadic_loft_refinement(bounds = thbspline.BoxList(  
            #                        thbspline.Box( ia(.0,.4), 
            #                                      ia(.0,.375) ) )    )
            ##
            ##****************************************************************
            ## 
            """hull test, DWL was not yet dyadic-compatible (FIXED)
            t11 = hull.dyadic_loft_refinement(bounds = thbspline.BoxList(  
                                    thbspline.Box( ia(.25,.75), 
                                                  ia(.25,.75) ) )    )
            t12 = t11.dyadic_loft_refinement(bounds = thbspline.BoxList(  
                                    thbspline.Box( ia(.3,.7), 
                                                      ia(.3,.7) ) )   ) 
            #"""
            """ nice looking set:
            t11 = hull.dyadic_loft_refinement(bounds = thbspline.BoxList(  
                                    thbspline.Box( ia(.1,.9), 
                                                  ia(.125,.875) ) )    )
            t12 = t11.dyadic_loft_refinement(bounds = thbspline.BoxList(  
                                    thbspline.Box( ia(.3125,.6875), 
                                                  ia(.3125,.6875) ) )    )
            #"""
            #"""######
            #
            #**********************************************************
            #
            # t2.rotate_surface((0.,0.,1.),np.pi/2.)
            #
            # Now the actual BULBOUS BOW
            t2 = thbspline.THBsurface(vertices = t2.vertices,
                                          tcurvelist = t2.tcurvelist,
                                          ulevel = 2,
                                          lcurvelist = t2.lcurvelist,
                                          vlevel = 1,
                                          trm = rm_11x7,
                                          lrm = rm_7x5)
            
            t3 = t2.reflect_surface((0.,0.,1.),(0.,0.,1.))
            t3 = t3.reflect_surface((1.,0.,0.),(1.,0.,0.))
            t3 = t3.rotate_surface((0.,0.,1.),np.pi/4.)
            t3 = t3.translate_surface( dz=-np.amin(t3.vertices[:,:,2]) )
            t3.set_top()
            
            
            
            #t2.compatibility()
            #t3.compatibility()
            #            t2.plotMultiSurface([t2,t3],
            #                                limx=100.,
            #                                limy=100.,
            #                                limz=100.)
                        #"""
            """#turn off refinement here
            t21 = t2.dyadic_loft_refinement(bounds = thbspline.BoxList(  
                                    thbspline.Box( ia(.25,.75), 
                                                  ia(.1,.9) ) )    )
            #"""
            """
            t22 = t21.dyadic_loft_refinement(bounds = thbspline.BoxList(  
                                    thbspline.Box( ia(.3125,.6875), 
                                                  ia(.125,.875) ) )    )
            #""" #does this show a bug?:
            #        t21 = t2.dyadic_loft_refinement()
            #        t22 = t21.dyadic_loft_refinement()
            
            #"""
            
            # BBOW refinement to show the levels:
            """ #
            t21 = t2.dyadic_loft_refinement(bounds = thbspline.BoxList(  
                    thbspline.Box( ia(.125,.875), 
                                  ia(.125,.875) ) )    )
            t22 = t21.dyadic_loft_refinement(bounds = thbspline.BoxList(  
                    thbspline.Box( ia(.3,.8), 
                                  ia(.3,.8) ) )    )
            
            #"""
            
            #bare hull THB surface evaluation:
            hull.compute_surface()
            hull.plotSurface(limx=100.,
                             limy=100.,
                             limz=100.,
                             minx=-100.,
                             miny=-100.,
                             minz=-100.)
            
            #"""
            #BBOW THB surface evaluation:
            #t2.compute_surface()
            #t2.plotSurface()
            #t2.plotSurface(limx=10.,limy=10.,limz=10.)
            #"""
            """
            final_dv = t2.store_final_dv
            intermediate_dv = t2.store_intermediate_dv
            nu = t2.store_nu
            nv = t2.store_nv
            lrm  = t2.store_lrm
            trm  = t2.store_trm
            lP = t2.store_lP
            P = t2.store_P
            this = t2.child
            #"""
            
            print '\n\nBare Hull Evaluation'
            print hull.eval(.3,.8)
            print hull.surface_point(.3,.8)
            print t11.surface_point(.3,.8)
            #print t12.surface_point(.3,.8)
            print hull.eval_old(.3,.8)
            
            
            """
            print '\n\nBulbBow Evaluation'
            print t2.eval(.3,.8)
            print t2.surface_point(.3,.8)
            print t21.surface_point(.3,.8)
            print t22.surface_point(.3,.8)
            #"""
            
            #        s0 = thbspline.THBsurface(verbose=False)
            #        s0 = s0.dyadic_loft_refinement()
            #        s0 = s0.dyadic_loft_refinement()
            #        #s0.vertices = hull.vertices
            #        #s0 = s0.dyadic_loft_refinement()
            #        s0.set_top()
            #        s1 = s0.dyadic_loft_refinement()
            #        s0.compute_surface()
            
            
            #s0.vertices[:] = hull.vertices[:]
            #s1.vertices[:] = t11.vertices[:]
            
            #        print '\n\ntrace element differences'
            #        for u,v in zip(np.linspace(0.,1.,10,endpoint=True),np.linspace(0.,1.,10,endpoint=True)):
            #            print '\n ------'
            #            print ' {},{}'.format(u,v)
            #            print hull.eval(u,v) - hull.eval_old(u,v)
            #            print hull.eval(u,v) - hull.surface_point(u,v)
            #            print hull.eval_old(u,v) - hull.surface_point(u,v)
            #            print 
            #            print hull.eval(u,v) - s0.eval(u,v) 
            #            print hull.surface_point(u,v) - s0.eval(u,v) 
            #            print ''
            
            
            """check the difference between B-spline and THBspline
            #suface evaluation for a given THBsurface called hull:
            #"""
            
            """
            hullcurves, ucurves, vcurves = hull.get_level_all_curves(hull.child.child.level)
            
            ucurves[0].plot3DmultiList(vcurves,hullcurves,
                   limx=100.,limy=100.,limz=100.)
            ucurves[0].plot3DmultiList(vcurves,ucurves,
                   limx=100.,limy=100.,limz=100.)
            #"""
            
            #"""
            t3.compute_surface()
            
            hull.compatibility()
            t3.compatibility()
            
            t2.compatibility()
            #t2.compute_surface()
            
            hull.plotMultiSurface([hull,t3],#,t2],
                                  limx=100.,limy=100.,limz=100.)
            #"""
            
        
        hull.set_top()
        av = hull.compute_true_nonzero_basis(.1,.1)
        tapb = hull.compute_total_active_by_level()
        
        t3_tapb = tapb[2]
        t3.vertices[t3_tapb.uactive[0],t3_tapb.vactive[0]]
        
        if profileit:
            uopt.THBspline_to_Bspline_ckdiff_evaluator(hull)#,mod_level=2)
            uopt.THBspline_to_Bspline_ckdiff_evaluator(t2)


    
    if True:
        #"""
        ##
        ##*********************************************************************
        ## ADD THE BULBOUS BOW:
        ## 
        ## incorporating the bow the implicit way (direct -> not blending)
        ## 
        ## principle method:  setup the transverse bulb incorporating curves
        ## then from that reconstruct the augmented hull vertices
        ## only at top level of refinement of course
        ##
        """
            curve 1 of 5:
                
            *) the longitudinal CP of the bbow runs aft to fwd
                    but the transverse cp curve runs bottom to top,
                    corresponding to running fwd to aft, along this portion.
                    
                    
            t3.lcurvelist[-1].vertices[0]
            Out[171]: array([  8.26802445e-16,   5.14141409e+00,   8.34651979e+00])
            
            t3.lcurvelist[-1].vertices[-1]
            Out[172]: array([  6.16632013e-16,   3.83448373e+00,   8.05712149e-01])
            
            so flip it in z with:
                
                bblcp1 = t3.lcurvelist[-1].vertices[::-1]
                
                starts at .8057
                runs to 8.34
                smaller z is forward
            
            *) the bbow transverse curves, after spinning and mirroring,
                    run
                    
                    curve.plot3DmultiList(t3.lcurvelist,
                                          [t3.tcurvelist[-1],t3.tcurvelist[-2]],
                                            minx=-100.
                                            ,miny=-100.,
                                            minz=-100.,
                                            limx=100.,
                                            limy=100.,
                                            limz=100.)
                    
                    so you want 
                        t3.tcurvelist[-1] 
                    together with 
                        bblcp1
                    to augment the CP transverse fwd curve of the bare hull
                
                    (note that yes, thankfully, tcurve still start keel, midship,
                    and run to waterline, breadth at their station.)
            
        """
        #for all new bulb maker curves:)
        hull.set_top()
        av = hull.compute_true_nonzero_basis(.1,.1)
        tapb = hull.compute_total_active_by_level()
        
        t3_tapb = tapb[2]
        t3.vertices[t3_tapb.uactive[0],t3_tapb.vactive[0]]
        
        L3_ini_vertices = copy.deepcopy(hull.child.child.vertices)
        
        bblcp1 = t3.lcurvelist[-1].vertices[::-1] #flip top long curve in z
        """
        """
        #1.)
        n1 = t3.tcurvelist[-1].n
        n2 = t3.lcurvelist[-1].n-0
        hull.child.child.tcurvelist[0].vertices[0:n1] = t3.tcurvelist[-1].vertices
        hull.child.child.tcurvelist[0].vertices[n1:n1+n2] = bblcp1
        
        """
            curve 2 of 5
            
            curve.plot3DmultiList([t3.tcurvelist[0],t3.tcurvelist[1]],
                                  t3.lcurvelist,
                                minx=-100.
                                ,miny=-100.,
                                minz=-100.,
                                limx=100.,
                                limy=100.,
                                limz=100.)
        """
        #2.)
        n1 = t3.tcurvelist[-2].n
        n2 = t3.lcurvelist[-1].n-1
        hull.child.child.tcurvelist[1].vertices[0:n1] = t3.tcurvelist[-2].vertices
        hull.child.child.tcurvelist[1].vertices[n1:n1+n2] = bblcp1[1:]
        """
        """
        #3.)
        n1 = t3.tcurvelist[-3].n
        n2 = t3.lcurvelist[-1].n-2
        hull.child.child.tcurvelist[2].vertices[0:n1] = t3.tcurvelist[-3].vertices
        hull.child.child.tcurvelist[2].vertices[n1:n1+n2] = bblcp1[2:]
        """
        """
        #4.)
        n1 = t3.tcurvelist[-4].n
        n2 = t3.lcurvelist[-1].n-3
        hull.child.child.tcurvelist[3].vertices[0:n1] = t3.tcurvelist[-4].vertices
        hull.child.child.tcurvelist[3].vertices[n1:n1+n2] = bblcp1[3:]
        """
        """
        #5.)
        n1 = t3.tcurvelist[-5].n
        n2 = t3.lcurvelist[-1].n-4
        hull.child.child.tcurvelist[4].vertices[0:n1] = t3.tcurvelist[-5].vertices
        hull.child.child.tcurvelist[4].vertices[n1:n1+n2] = bblcp1[4]
        """
        """
        ## now get vertices from the tcurve you just set:
        hull.child.child.set_vertices_from_tcurves()
        ## now set the longitudinals from the vertices (inefficient - sets tcurves too):
        hull.child.child.set_THB_curves_from_vertex_net()
        """
        for jfwd,jbckwd in zip(t3_tapb.vactive[0:5],
                               t3_tapb.vactive[0:5][::-1]):
            #print jfwd, jbckwd
            hull.child.child.vertices[t3_tapb.uactive[0:7],jfwd] = \
                            t3.vertices[:,jbckwd]
        ##       
        ## transverse boundary:
        for jfwd in t3_tapb.vactive[5:8]:
            #print jfwd
            hull.child.child.vertices[t3_tapb.uactive[0:7],jfwd] = \
                            t3.vertices[:,0]
        
        ##
        ## longitudinal boundary:
        for ifwd in t3_tapb.uactive[7:10]:
            hull.child.child.vertices[ifwd,t3_tapb.vactive[0:5]] = \
                t3.vertices[-1,:]
        
        #"""
        hull.expand_initial_curve_list(override=True)
        #hull.child.expand_initial_curve_list(override=True)
        #hull.child.child.expand_initial_curve_list(override=True)
        ##
        ##*********************************************************************
        ##
        #"""
        
        #j=0
        #for i in :
            
        def plotTHBhull():
            """
                plotTHBhull()
                
                hull.compute_surface()
            """
            hull.compute_surface()
            hull.plotSurface()
            curve = hull.tcurvelist[0]
            
            xb,yb,zb = self.plot_get_extremes()
            curve.plot3DmultiList(hull.child.child.lcurvelist,
                                  hull.child.child.tcurvelist,
                                  minx=xb[0],
                                  miny=yb[0],
                                  minz=zb[0],
                                  limx=xb[1],
                                  limy=yb[1],
                                  limz=zb[1])
            
#            curve.plot3DmultiList(hull.lcurvelist,
#                                  hull.hullcurvelist,
#                                  minx=-100.,
#                                 miny=-100.,
#                                 minz=-100.,
#                                 limx=100.,
#                                 limy=100.,
#                                 limz=100.)
            
            
            
#            curve = hull.child.child.tcurvelist[0]
#            curve.plot3DmultiList(hull.child.child.lcurvelist,
#                                  hull.tcurvelist,
#                                 minx=-100.,
#                                 miny=-100.,
#                                 minz=-100.,
#                                 limx=100.,
#                                 limy=100.,
#                                 limz=100.)
            
            return
            
        def reset():
            """
                reset()
            """
            for i in t3_tapb.uactive:
                for j in t3_tapb.vactive:
                    hull.child.child.vertices[i,j] = L3_ini_vertices[i,j]
            return
        
        hull.child.child.set_THB_curves_from_vertex_net()
        hull.tcurvelist[0].plot3DmultiList(hull.tcurvelist,
                                           hull.lcurvelist)
        """
        
        hull.tcurvelist[0].plot3DmultiList(hull.child.child.tcurvelist,[],
                                           minx=-100.,
                                           miny=-100.,
                                           minz=-100.,
                                           limx=100.,
                                           limy=100.,
                                           limz=100.)
        
        hull.tcurvelist[0].plot3DmultiList(hull.child.child.tcurvelist[0:4],[],
                                           minx=-100.,
                                           miny=-100.,
                                           minz=-100.,
                                           limx=100.,
                                           limy=100.,
                                           limz=100.)
        
        hull.tcurvelist[0].plot3DmultiList(hull.child.child.tcurvelist[2:5],[],
                                           minx=-100.,
                                           miny=-100.,
                                           minz=-100.,
                                           limx=100.,
                                           limy=100.,
                                           limz=100.)
        
        
        hull.tcurvelist[0].plot3DmultiList(hull.child.child.lcurvelist,[],
                                           minx=-100.,
                                           miny=-100.,
                                           minz=-100.,
                                           limx=100.,
                                           limy=100.,
                                           limz=100.)
        
        
        hull.tcurvelist[0].plot3DmultiList([hull.child.child.tcurvelist[0]],
                                           [hull.child.child.tcurvelist[1]],
                                           minx=-100.,
                                           miny=-100.,
                                           minz=-100.,
                                           limx=100.,
                                           limy=100.,
                                           limz=100.)
        
        
        
        hull.tcurvelist[0].plot3DmultiList(hull.child.child.lcurvelist[5:8],
                                           hull.child.child.lcurvelist[:5],
                                           minx=-100.,
                                           miny=-100.,
                                           minz=-100.,
                                           limx=100.,
                                           limy=100.,
                                           limz=100.)
        
        
        hull.tcurvelist[0].plot3DmultiList(hull.child.child.lcurvelist[5:8],
                                           [hull.child.child.tcurvelist[0]],
                                           minx=-100.,
                                           miny=-100.,
                                           minz=-100.,
                                           limx=100.,
                                           limy=100.,
                                           limz=100.)
        
        hull.tcurvelist[0].plot3DmultiList([hull.child.child.tcurvelist[0]],[],
                                            minx=-100.
                                            ,miny=-100.,
                                            minz=-100.,
                                            limx=100.,
                                            limy=100.,
                                            limz=100.)
        
        hull.tcurvelist[0].plot3DmultiList([hull.child.child.tcurvelist[0],
                                            hull.child.child.tcurvelist[4]],[],
                                            minx=-100.
                                            ,miny=-100.,
                                            minz=-100.,
                                            limx=100.,
                                            limy=100.,
                                            limz=100.)
        
        hull.tcurvelist[0].plot3DmultiList(hull.child.child.tcurvelist[0:2],[],
                                            minx=-100.
                                            ,miny=-100.,
                                            minz=-100.,
                                            limx=100.,
                                            limy=100.,
                                            limz=100.)
        
        hull.tcurvelist[0].plot3DmultiList(hull.child.child.tcurvelist[0:3],[],
                                            minx=-100.
                                            ,miny=-100.,
                                            minz=-100.,
                                            limx=100.,
                                            limy=100.,
                                            limz=100.)
        
        hull.tcurvelist[0].plot3DmultiList(hull.child.child.tcurvelist[0:4],[],
                                            minx=-100.
                                            ,miny=-100.,
                                            minz=-100.,
                                            limx=100.,
                                            limy=100.,
                                            limz=100.)
        
        hull.tcurvelist[0].plot3DmultiList(hull.child.child.tcurvelist[4:9],[],
                                            minx=-100.
                                            ,miny=-100.,
                                            minz=-100.,
                                            limx=100.,
                                            limy=100.,
                                            limz=100.)
        
        curve.plot3DmultiList([t3.lcurvelist[0],t3.lcurvelist[-1]],
                               [hull.child.child.tcurvelist[3]],
                                            minx=-100.
                                            ,miny=-100.,
                                            minz=-100.,
                                            limx=100.,
                                            limy=100.,
                                            limz=100.)
        
        
        
        curve.plot3DmultiList(hull.child.child.tcurvelist,
                                  hull.child.child.lcurvelist,
                                  minx=-100.,
                                 miny=-100.,
                                 minz=-100.,
                                 limx=100.,
                                 limy=100.,
                                 limz=100.)
        
        curve.plot3DmultiList(hull.child.tcurvelist,
                              hull.child.lcurvelist,
                              minx=-100.,
                             miny=-100.,
                             minz=-100.,
                             limx=100.,
                             limy=100.,
                             limz=100.)
    
        curve.plot3DmultiList(hull.tcurvelist,
                              hull.lcurvelist,
                              minx=-100.,
                             miny=-100.,
                             minz=-100.,
                             limx=100.,
                             limy=100.,
                             limz=100.)
        #"""
        
        #hull.numpu = 30
        #hull.numpv = 30
        plotTHBhull()
        
        
        if profileit:
            uopt.THBspline_to_Bspline_ckdiff_evaluator(hull,mod_level=2)
            uopt.THBspline_to_Bspline_ckdiff_evaluator(t3)