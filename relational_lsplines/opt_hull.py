# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 16:16:11 2016

@author: lukemcculloch

caresian product from sets:
    http://metapython.blogspot.com/2011/10/cartesian-product-for-sets.html

"""
"""
TODO: fix Xc and len SAC parts
to always work!
"""
import numpy as np
import itertools
#from interval_arithmetic import ia
#
from extended_interval_arithmetic import ia
import sqKanren as lp
#from design_space import HullSpace as hullclp
from design_tree import HullDesignNode as hullclp
from design_tree import DesignTree
#from hull_from_designspace import Hull as design_builder #builds a single design instance from  a thin HullCLP 
import hull_from_designspace as hullmethods
#
#import hull_use_HICLP
#
import sobol
import copy




instancemethod = lambda f: (lambda *a, **kw: f(*a, **kw))
class mset(set):
    __mul__ = __rmul__ = instancemethod(itertools.product)
class strset(str):
    __mul__ = __rmul__ = instancemethod(itertools.product)
    

class DesignNode(object):
    """wraps the States Class
    logical resolution produces children 
    e.g. with contracted states
    
    States      : cargo
    parent      : parent
    children    : any number of Design Nodes Derived From this one
    """
    def __init__(self,states=None, parent=None):
        self._states = states 
        self.parent = parent
        self.children = []
        
    @property
    def states(self):
        return self._states
    @states.setter
    def states(self, states):
        self.children += DesignNode(states, parent=self)
        
    #    @property
    #    def children(self):
    #        return self._children
    #    @children.setter
    #    def states(self, states):
    #        self._children += DesignNode(states, parent=self)
    #        
    #    @property
    #    def parent(self):
    #        return self._parent
    #    @parent.setter
    #    def states(self, parent):
    #        self._parent = parent
    
    



#class DesignTree(hullclp):
#    """This is like the States Class
#    Except contraction produce children with contracted states
#    """
#    def __init__(self, root):
#        self.root = Node(root)
#    
#    #    @property
#    #    def root(self):
#    #        return self._root
#    #    @root.setter
#    #    def root(self, state):
#    #        return





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
        self.Cmidshp = Cmidshp
        self.Cwp = Cwp
        self.Ccp = Ccp
        
    
    #def set_ranges_for_details(self):
    #    self.Acp = 

class DesignParams(object):
    def __init__(self,
                 lwl,
                 draft,
                 bwl,
                 vol,
                 LCG,
                 Cmidshp,
                 Cb,
                 Cwp,
                 Ccp):
        self.lwl = lwl
        self.draft = draft
        self.bwl = bwl
        self.vol = vol
        self.Cb = Cb
        self.Cwp = Cwp
        self.LCG = LCG
        self.Cmidshp = Cmidshp
        self.Ccp = Ccp
        
    def print_params(self):
        print 'lwl ',self.lwl
        print 'draft ',self.draft
        print 'bwl ',self.bwl
        print 'vol ',self.vol
        print 'Cb ',self.Cb
        print 'Cwp',self.Cwp
        print 'LCG ',self.LCG
        print 'Cmidshp ',self.Cmidshp
        print 'Ccp ',self.Ccp
    
        

#class DesignSpace(hull_use_HICLP.hullclp):
class DesignSpace(object):
    SMALL = .1
    def __init__(self, spec, verbose=False):
        self._verbose = verbose
        self.spec = spec
        self.sobol_seq = sobol.sobolSeq([1,1],[1,1])
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
        self.hdp.Cmidshp = ia(spec.Cmidshp[0],spec.Cmidshp[1])
        self.hdp.AC_revise()
        #"""
        self.hdp.states = self.hdp.shape_bow_section_inizers()
        self.hdp.states = self.hdp.shape_stern_section_inizers()
        self.hdp.states = self.hdp.SAC_run_area_ini()
        self.hdp.states = self.hdp.SAC_Xc_init()
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
        k4 = set(design.BowCurveCoeff) & set(design.SternCurveCoeff)
        ksac = set(design.list_SAC)
        kr = set(keys) - k1 - k2 - k3 - k4 -ksac
        
        
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
        design_sq = self.search_loop(design_sq,
                                     k4,
                                     sobol_seq)
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
    print 'Xc fwd = ',a, self.hdp(self.hdp.SAC_fwd_Xc)
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
    hullmethods.compute_hull(hull2)
    
    
    hull2.Lsplines.DWL.curve.plot()
    hull2.Lsplines.CProfile.curve.plot()
    hull2.Lsplines.SAC_entrance.curve.plot()
    hull2.Lsplines.SAC_mid.curve.plot()
    hull2.Lsplines.SAC_run.curve.plot()
    
    hullmethods.hull_gui(hull2)
    
    hull2.DWL.plot3DmultiList(hull2.tcurvenet,hull2.lcurvenet)
    return hull2
    
    
if __name__ == '__main__':
    """
        Design spec:
        Ccp -> fraction of depth*lwl which makes CPkeel area
        Cwp -> fraction of bwl*lwl which makes DWL area
        
        When Cb is high and Cwp/Ccp are low, DWL/CPKeel end up wavey
    """
    LITTLE = 1.e-2
    hull_template = hullclp()# hull_use_HICLP.hullclp()
    
    spec = DesignSpecification(
                            lwl = (110.,116.),
                            draft = (10.,20.),
                            bwl = (10.,20.),
                            #vol=(18000.,20000.),
                            vol=(25000.,35000.),
                            LCG = (54.,59.),
                            Cb = (0.6,0.85),
                            Cwp = (.6,.89),#water plane coefficient
                            Cmidshp = (0.89,0.99),
                            Ccp = (.6,0.99))#centerplane coefficient
                            
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
        
        