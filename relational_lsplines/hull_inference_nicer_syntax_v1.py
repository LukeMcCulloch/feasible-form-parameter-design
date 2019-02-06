# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 18:44:37 2016

@author: lukemcculloch

help(lp.run_br)
"""
"""
@author: lukemcculloch

probably best called:
    Constraint Imperative Programming
    
if a function result changes, which internals changed?  I cannot know, ergo it is inefficient.

I need a map of variables to varables that could change (done var to arc(s) where arcs naturally unify vars)

Propogation of Information:
    August 6, 2016
    TODO:
        when a var state changes naturlly an arc, 
        the var doesn't find out toi do on the changed set.
        I need a trigger from state change 
        to put it back on the set of vars that changed
        
    TODO:
        multiple state
        -During unification, if a var isn't consistent, make a new state?
        -go back and expand the constraint so the return multiple states?
"""

#import tokenize
import numpy as np
import copy
from automatic_differentiation import ad
from interval_arithmetic import ia

import uKanren as lp




class Hull_NOPE(lp.States):
    """Hull Object
    High level input system for consistent 
    hull design algorithm direction
    using constraint logic programming
    
    This class holds methods which are themselves 
    constrants built out of more primitive constraints
    fromthe uKanren module
    
    HLISFCHDADUCLP \n
    CLIPFFPHD
    TODO:find a recursive acronym for this
    
    Arcs        : Hull Parameter Relationships
    Variables   : High Level Hull Design Parameters
    
    Potentially improve/shorten by using memoization decorators
    see
    https://wiki.python.org/moin/PythonDecoratorLibrary
    
    ...Fix to work on lists of states?  (proper ukanren)
    ...Overload variable to have arithmetic opertions 
    which explicitly take a state?
    ...Genaricize functions and turn input into a parsable language?
    --------------------------------------------------
    """
    def __init__(self, state=None):
        super(Hull, self).__init__(self._setup_())
        #self.graph = {} #not used
        #compatability:
        self.state = state or self._setup_()

    
    
    def _setup_(self):
        """Set up the design variables
        with maps to their usage in arcs
        """
        self._draft   = lp.Variable('draft',
                                    arcs=[self.compute_Cb,
                                          self.compute_FOS,
                                          self.compute_Ccp,
                                          self.compute_FOK,
                                          self.compute_drafts])
        
        self._dsmax   = lp.Variable('dsmax',
                                    arcs=[self.compute_midship_coefficient,
                                          self.compute_drafts])
        
        self._vol     = lp.Variable('vol',
                                    arcs=[self.compute_Cb,
                                          self.compute_FOSAC,
                                          self.compute_FOS,
                                          self.compute_SAC_section_properties])
                                    
        self._Cb      = lp.Variable('Cb',
                                    arcs=[self.compute_Cb,
                                          self.compute_prismatic_coefficient])
        #
        #-------------------------SAC LCG
        #         
        self._LCG      = lp.Variable('LCG',
                                    arcs=[self.compute_LCG])   
        self._SAC_fwd_Xc = lp.Variable('SAC_fwd_Xc',
                                    arcs=[self.compute_LCG])   
        self._SAC_mid_Xc = lp.Variable('SAC_mid_Xc',
                                    arcs=[self.compute_LCG])   
        self._SAC_run_Xc = lp.Variable('SAC_run_Xc',
                                    arcs=[self.compute_LCG])
                                    
        #
        #-------------------------SAC LCG
        #
        self._lwl     = lp.Variable('lwl',
                                    arcs=[self.compute_Cb,
                                          self.compute_Cwp,
                                          self.compute_flat_relations,
                                          self.compute_Ccp,
                                          self.compute_LCG,
                                          self.compute_SAC_section_properties])
                                          
        self._bwl     = lp.Variable('bwl',
                                    arcs=[self.compute_Cb,
                                          self.compute_Cwp,
                                          self.compute_midship_coefficient,
                                          self.compute_FOWL,
                                          self.compute_FOS])
        #waterplane area
        self._Awp     = lp.Variable('Awp',
                                    arcs=[self.compute_Cwp,
                                          self.compute_FOWL]) 
                                    
        self._Cwp     = lp.Variable('Cwp',
                                    arcs=[self.compute_Cwp])
                                    
        self._Acp     = lp.Variable('Acp',
                                    arcs=[self.compute_Ccp,
                                          self.compute_FOK])
        self._Ccp     = lp.Variable('Ccp',
                                    arcs=[self.compute_Ccp])
        
        #Area midship - area of largest midship section
        self._Amsh    = lp.Variable('Amsh',
                                    arcs=[self.compute_midship_coefficient,
                                          self.compute_FOSAC,
                                          self.compute_FOS])
        
        #midship coeff
        self._Cmidshp    = lp.Variable('Cmidshp',
                                    arcs=[self.compute_midship_coefficient,
                                          self.compute_prismatic_coefficient]) 
        
        self._Cp      = lp.Variable('Cp',
                                    arcs=[self.compute_prismatic_coefficient]) 
        
        
        self._lfwl   = lp.Variable('lfwl',
                                    arcs=[self.compute_FOWL,
                                         self.compute_flat_relations])
        
        self._lfos   = lp.Variable('lfos',
                                    arcs=[self.compute_FOS,
                                         self.compute_flat_relations])
        
        self._Afos   = lp.Variable('Afos',
                                    arcs=[])
                                    
        self._lfsac  = lp.Variable('lfsac',
                                   arcs=[self.compute_FOSAC,
                                         self.compute_flat_relations,
                                         self.compute_SAC_section_properties])
                                    
        self._lfcp  = lp.Variable('lfcp',
                                   arcs=[self.compute_flat_relations,
                                          self.compute_FOK])
        ##
        ## Bow Fairness curve
        ##
        self._bbfc  = lp.Variable('bbfc',
                                   arcs=[self.compute_bow_fairness_section])
        self._dbfc  = lp.Variable('dbfc',
                                   arcs=[self.compute_bow_fairness_section])
        self._Abfc  = lp.Variable('Abfc',
                                   arcs=[self.compute_bow_fairness_section])
        self._Cbfc  = lp.Variable('Cbfc',
                                   arcs=[self.compute_bow_fairness_section])
        ##
        ## Stern Fairness curve
        ##
        """TBD"""
        self._bsfc  = lp.Variable('bsfc',
                                   arcs=[self.compute_stern_fairness_section])
        self._dsfc  = lp.Variable('dsfc',
                                   arcs=[self.compute_stern_fairness_section])
        self._Asfc  = lp.Variable('Asfc',
                                   arcs=[self.compute_stern_fairness_section])
        self._Csfc  = lp.Variable('Csfc',
                                   arcs=[self.compute_stern_fairness_section])
        ##
        ## Multi SAC
        ##
        self._SAC_entrance_len = lp.Variable('SAC_entrance_len',
                                             arcs = [self.compute_SAC_section_properties])
        self._SAC_mid_len = lp.Variable('SAC_mid_len',
                                             arcs = [self.compute_SAC_section_properties])
        self._SAC_run_len = lp.Variable('SAC_run_len',
                                             arcs = [self.compute_SAC_section_properties])
        self._SAC_entrance_area = lp.Variable('SAC_entrance_area',
                                             arcs = [self.compute_SAC_section_properties])
        self._SAC_mid_area = lp.Variable('SAC_mid_area',
                                             arcs = [self.compute_SAC_section_properties])
        self._SAC_run_area = lp.Variable('SAC_run_area',
                                             arcs = [self.compute_SAC_section_properties])
        ##
        ##
        ##
        s  = lp.State(values={self._draft    : None,
                              self._dsmax    : None,
                              self._vol      : None,
                              self._LCG      : None,
                              self._Cb       : None,
                              self._lwl      : None,
                              self._bwl      : None,
                              self._Awp      : None,
                              self._Amsh     : None,
                              self._Cwp      : None,
                              self._Cmidshp  : None,
                              self._Cp       : None,
                              self._lfos     : None,
                              self._lfwl     : None,
                              self._lfsac    : None,
                              self._Afos     : None,
                              self._bbfc     : None,
                              self._dbfc     : None,
                              self._Abfc     : None,
                              self._Cbfc     : None,
                              self._bsfc     : None,
                              self._dsfc     : None,
                              self._Asfc     : None,
                              self._Csfc     : None,
                              self._SAC_entrance_len    : None,
                              self._SAC_mid_len         : None,
                              self._SAC_run_len         : None,
                              self._SAC_entrance_area   : None,
                              self._SAC_mid_area        : None,
                              self._SAC_run_area        : None,
                              self._SAC_fwd_Xc          : None,
                              self._SAC_mid_Xc          : None,
                              self._SAC_run_Xc          : None})
        self._set_observers_()
        return s
    def _set_observers_(self):
        """Set up the the sets of changed variables and their arcs 
        """
        self.arc_set = set([]) #active arcs
        self.var_set = set([]) #active variables 
        return
    
    def get_arc_set(self):
        """Find the constraints that have inputs which changed
        """
        for var in self.var_set:
            self.arc_set = self.arc_set.union(set(var.arcs))
        return
    
    def AC_revise(self, print_=False, maxit = 10):
        """Arc Consistency for the Hull Parameters
        """
        self.get_arc_set()
        if print_:
            print 'arcs of interest:\n',self.arc_set
        self.var_set = set([])
        count=0
        while len(self.arc_set)>0 and count<maxit:
            for arc in self.arc_set:
                #arc = self.arc_set.pop()
                if print_: print '\ndoing ', arc
                arc()
            self.arc_set = set([])
            self.get_arc_set()
            self.var_set = set([])
            count += 1
            #print count
        return
    
    def equal(self, c1,c2):
        """lp.Goal.eq, 
        with 
        hull constraint graph propogation enabled
        
        How to handle when the both are lp.Variable ?
        answer, use the more general lp.eq(anything, anything)
        
        The purpose of this function is to 
        add self.c1 to the set of vars to 
        trigger Arc Consistency
        """
        #if isinstance(lp.Variable, c1):
        #    var = c1
        #    val
        #assert(not isinstance(c2, lp.Variable))
        ck1 = self.state.value_of(c1)
        goal = lp.Goal.eq(c1, c2)
        self.state = goal(self.state)[0]
        ck2 = self.state.value_of(c1)
        if not (ck1 == ck2):
            self.var_set = self.var_set.union(set([c1]))
        self.AC_revise()
        return
        
    @property
    def states(self):
        return self._states
    @states.setter
    def states(self, states):
        self._states = states
    
    def clean_state(self, state, special=None):
        """Launguage is not quite up to par yet
        making do, clean extranious vars thatcrop up
        during rule computations
        """
        snm = [el.name for el in special]
        dkeys = []
        for key in state.values:
            if not isinstance(key, lp.Variable):
                dkeys.append(key)
            elif key.name in snm:
                dkeys.append(key)
        for key in dkeys:
            del(state.values[key])
        return state
        
    def set_updates(self, state, statei,vars):
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst)) 
        return
    
    def get_val(self, var):
        return self.state.value_of(var)
    
    @property
    def lwl(self):
        """Length of the Waterline
        """
        return self._lwl
    @lwl.setter
    def lwl(self, lwl):
        self.equal(self.lwl,lwl)

    @property
    def bwl(self):
        """Beam of the Waterplane
        """
        return self._bwl
    @bwl.setter
    def bwl(self, bwl):
        self.equal(self.bwl,bwl)
        
    @property
    def draft(self):
        """Draft at deepest point
        """
        return self._draft
    @draft.setter
    def draft(self, draft):
        self.equal(self.draft,draft)
    
    @property
    def dsmax(self):
        """midship draft dso
        TODO: rename
        """
        return self._dsmax
    @dsmax.setter
    def dsmax(self, dsmax):
        self.equal(self.dsmax,dsmax)
        
        
    @property
    def vol(self):
        return self._vol
    @vol.setter
    def vol(self, vol):
        print 'setting vol'
        self.equal(self.vol,vol)
        
        
    @property
    def Cb(self):
        return self._Cb
    @Cb.setter
    def Cb(self, Cb):
        self.equal(self.Cb,Cb)
    
    #
    #-------------------------SAC LCG
    #
    @property
    def LCG(self):
        """Longitudinal Center of Gravity of the vessel
        """
        return self._LCG
    @LCG.setter
    def LCG(self, val):
        self.equal(self.LCG,val)
    
    @property
    def SAC_fwd_Xc(self):
        """Xc of SAC fwd section
        """
        return self._SAC_fwd_Xc
    @SAC_fwd_Xc.setter
    def SAC_fwd_Xc(self, val):
        self.equal(self.SAC_fwd_Xc,val)
    
    @property
    def SAC_mid_Xc(self):
        """Xc of SAC mid section
        """
        return self._SAC_mid_Xc
    @SAC_mid_Xc.setter
    def SAC_mid_Xc(self, val):
        self.equal(self.SAC_mid_Xc,val)
    
    @property
    def SAC_run_Xc(self):
        """Xc of SAC aft section
        """
        return self._SAC_run_Xc
    @SAC_run_Xc.setter
    def SAC_run_Xc(self, val):
        self.equal(self.SAC_run_Xc,val)
    #
    #-------------------------
    #
        
    @property
    def Awp(self):
        return self._Awp
    @Awp.setter
    def Awp(self, Awp):
        self.equal(self.Awp,Awp)
        
        
    @property
    def Cwp(self):
        return self._Cwp
    @Cwp.setter
    def Cwp(self, Cwp):
        self.equal(self.Cwp,Cwp)
        
    @property
    def Acp(self):
        """Area Center Plane"""
        return self._Acp
    @Acp.setter
    def Acp(self, Acp):
        self.equal(self.Acp,Acp)
    
    @property
    def Ccp(self):
        """Coefficient Center Plane"""
        return self._Ccp
    @Ccp.setter
    def Ccp(self, Ccp):
        self.equal(self.Ccp,Ccp)
        
    @property
    def Amsh(self):
        return self._Amsh
    @Amsh.setter
    def Amsh(self, Amsh):
        self.equal(self.Amsh,Amsh)

        
    @property
    def Cmidshp(self):
        return self._Cmidshp
    @Cmidshp.setter
    def Cmidshp(self, Cmidshp):
        self.equal(self.Cmidshp,Cmidshp)

        
    @property
    def Cp(self):
        return self._Cp
    @Cp.setter
    def Cp(self, Cp):
        self.equal(self.Cp,Cp)
    
    @property
    def lfwl(self):
        return self._lfwl
    @lfwl.setter
    def lfwl(self, lfwl):
        self.equal(self.lfwl,lfwl)
        
    @property
    def lfos(self):
        return self._lfos
    @lfos.setter
    def lfos(self, lfos):
        self.equal(self.lfos,lfos)
    
    @property
    def lfsac(self):
        return self._lfsac
    @lfsac.setter
    def lfsac(self, lfsac):
        self.equal(self.lfsac,lfsac)
        
    @property
    def lfcp(self):
        return self._lfcp
    @lfcp.setter
    def lfcp(self, lfcp):
        self.equal(self.lfcp,lfcp)
        
    @property
    def Afos(self):
        return self._Afos
    @Afos.setter
    def Afos(self, Afos):
        self.equal(self.Afos,Afos)
    ##
    ## Bow Fairness curve
    ##
    @property
    def bbfc(self):
        return self._bbfc
    @bbfc.setter
    def bbfc(self, bbfc):
        self.equal(self.bbfc,bbfc)
    @property
    def dbfc(self):
        return self._dbfc
    @dbfc.setter
    def dbfc(self, dbfc):
        self.equal(self.dbfc,dbfc)
    @property
    def Abfc(self):
        return self._Abfc
    @Abfc.setter
    def Abfc(self, Abfc):
        self.equal(self.Abfc,Abfc)
    @property
    def Cbfc(self):
        return self._Cbfc
    @Cbfc.setter
    def Cbfc(self, Cbfc):
        self.equal(self.Cbfc,Cbfc)
    ##
    ## Stern Fairness curve
    ##
    @property
    def bsfc(self):
        return self._bsfc
    @bsfc.setter
    def bsfc(self, bsfc):
        self.equal(self.bsfc,bsfc)
    @property
    def dsfc(self):
        return self._dsfc
    @dsfc.setter
    def dsfc(self, dsfc):
        self.equal(self.dsfc,dsfc)
    @property
    def Asfc(self):
        return self._Asfc
    @Asfc.setter
    def Asfc(self, Asfc):
        self.equal(self.Asfc,Asfc)
    @property
    def Csfc(self):
        return self._Csfc
    @Csfc.setter
    def Csfc(self, other):
        self.equal(self.Csfc,other)
    ##
    ## SAC in 3 parts
    ##
    ## SAC length
    @property
    def SAC_entrance_len(self):
        return self._SAC_entrance_len
    @SAC_entrance_len.setter
    def SAC_entrance_len(self, other):
        self.equal(self.SAC_entrance_len,other)
    @property
    def SAC_mid_len(self):
        return self._SAC_mid_len
    @SAC_mid_len.setter
    def SAC_mid_len(self, other):
        self.equal(self.SAC_mid_len,other)
    @property
    def SAC_run_len(self):
        return self._SAC_run_len
    @SAC_run_len.setter
    def SAC_run_len(self, other):
        self.equal(self.SAC_run_len,other)
    ## SAC area
    @property
    def SAC_entrance_area(self):
        return self._SAC_entrance_area
    @SAC_entrance_area.setter
    def SAC_entrance_area(self, other):
        self.equal(self.SAC_entrance_area,other)
    @property
    def SAC_mid_area(self):
        return self._SAC_mid_area
    @SAC_mid_area.setter
    def SAC_mid_area(self, other):
        self.equal(self.SAC_mid_area,other)
    @property
    def SAC_run_area(self):
        return self._SAC_run_area
    @SAC_run_area.setter
    def SAC_run_area(self, other):
        self.equal(self.SAC_run_area,other)
    ##
    ##
    ##
    def compute_FOSAC(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.flat_of_sac_constraint()
        return
        
    def compute_LCG(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.LCG_constraints()
        return
    
    def compute_drafts(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.relate_drafts()
        return
        
    def compute_FOWL(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.flat_of_wl_constraint()
        return
        
    def compute_FOK(self, **args):
        """Flat of Keel
        """
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.flat_of_keel_constraint()
        return
        
    def compute_FOS(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.flat_of_side_constraint()
        return
        
        
    def compute_Cb(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.block_coefficient()
        return
    
    def compute_Cwp(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.waterplane_coefficient()
        return
    
    def compute_Ccp(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.center_plane_coefficient()
        return
        
    def compute_midship_coefficient(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.midship_coefficient()
        return
        
    def compute_prismatic_coefficient(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.prismatic_coefficient()
        return
    
    def compute_flat_relations(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.relate_flats()
        return
    
    def compute_bow_fairness_section(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.bow_fairness_section_area()
        return
    
    def compute_stern_fairness_section(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.stern_fairness_section_area()
        return
    
    def compute_SAC_section_properties(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.SAC_three_parts()
        self.state = self.SAC_run_area_consistency()
        self.state = self.SAC_entrance_len_consistency()
        self.state = self.constrain_SAC_Xc()
        return
        

    def set_generic(self, func, **args):
        """not used: --not fully thought out--
            would be used as a decorator 
            for the -compute_some_hull_property- 
            functions above, to shorten their code
        """
        def wrap(**args):
            if args:
                for arg in args:
                    key = getattr(self, arg) 
                    self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
            self.state = func()
            return
        return wrap
    
    def var_from_string(self, key):
        if isinstance(key, str):
            return self.__dict__['_'+str(key)]
        if isinstance(key, lp.Variable):
            return self.__dict__['_'+str(key.name)]
    
    def block_coefficient(self):
        """
        lpp     : length between perpendiculars or wL length
        br      : some breadth (e.g. breadth at the waterline)
        draft   : draft
        vol     : volumetric displacement/water weight
        
                Vol/(Lwl*Bwl*Draft) = Cb
        """
        lpp = self.lwl
        br = self.bwl
        draft = self.draft
        vol = self.vol
        Cb = self.Cb
        s = self.state
        
        c1 = lp.Variable('c1')
        c2 = lp.Variable('c2')
    
        #state = copy.copy(s)
        s.values[c1]=None
        s.values[c2]=None
        
        #print s
        
        goal = lp.Goal.both(lp.Goal.both(lp.Goal.mulo(lpp,br,c1),
                                         lp.Goal.mulo(c1,draft,c2)),
                            lp.Goal.divo(vol,c2,Cb)
                            )

        vars = [lpp,br,draft,vol,Cb]
        statei = copy.copy(s)
        state = goal(goal(goal(s)[0])[0])[0] #forward chaining!  maybe decorate or otherwise move outside
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))
        
        del(state.values[c1])
        del(state.values[c2])
        state = self.clean_state(state,[c1,c2])
        return state 
    
    
    def waterplane_coefficient(self):
        """
        lwl     : waterline length
        bwl     : waterline breadth
        Awp     : wateplane area
        Cwp     : waterplane coefficient
        
                Awp/(lwl*bwl) = Cwp
        """
        lwl = self.lwl
        bwl = self.bwl
        Awp = self.Awp
        Cwp = self.Cwp
        s = self.state
        
        c1 = lp.Variable('c1')
        #state = copy.copy(s)
        s.values[c1]=None
        
        #wrap it in another function in uKanren
        goal = lp.Goal.both(lp.Goal.mulo(lwl,bwl,c1),
                            lp.Goal.divo(Awp,c1,Cwp)
                            )
        
        vars = [lwl,bwl,Awp,Cwp]
        statei = copy.copy(s)
        state = goal(goal(s)[0])[0]
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))
        
        del(state.values[c1])
        state = self.clean_state(state,[c1])
        return state
        
        
    def center_plane_coefficient(self):
        """
        lwl     : waterline length
        draft   : max depth
        Acp     : center plane area
        Ccp     : center plane coefficient
        
                Acp/(Lwl*draft) = Ccp
        """
        lwl     = self.lwl
        draft   = self.draft
        Acp     = self.Acp
        Ccp     = self.Ccp
        s       = self.state
        
        c1 = lp.Variable('c1')
        #state = copy.copy(s)
        s.values[c1]=None
        
        #wrap it in another function in uKanren
        goal = lp.Goal.both(lp.Goal.mulo(lwl,draft,c1),
                            lp.Goal.divo(Acp,c1,Ccp)
                            )
        
        vars = [lwl,draft,Acp,Ccp]
        statei = copy.copy(s)
        state = goal(goal(s)[0])[0]
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))
        
        del(state.values[c1])
        state = self.clean_state(state,[c1])
        return state

        
    
        
    def midship_coefficient(self):
        """
        bwl     : max wL breadth
        Amsh    : Midship area
        Dsmax   : Draft at station of max section area
        Cmidshp : Midship Coefficient
        
                Amsh/(Bwl*Dsmax) = Cmidshp
        """
        bwl = self.bwl
        dsmax = self.dsmax
        Amsh = self.Amsh
        Cmidshp = self.Cmidshp
        s = self.state
        
        c1 = lp.Variable('c1')
        #state = copy.copy(s)
        s.values[c1]=None
        goal = lp.Goal.both(lp.Goal.mulo(dsmax,bwl,c1),
                            lp.Goal.divo(Amsh,c1,Cmidshp)
                            )
                            
        vars = [dsmax,bwl,Amsh,Cmidshp]
        statei = copy.copy(s)
        
        state = goal(s)[0]
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))

        del(state.values[c1])
        state = self.clean_state(state,[c1])
        return state
        
        
    def prismatic_coefficient(self):
        """
        Cb / Cmidshp = Cp
        """
        Cb = self.Cb
        Cmidshp = self.Cmidshp
        Cp = self.Cp
        s = self.state
        
        goal = lp.Goal.divo(Cb,Cmidshp,Cp)
        vars = [Cb,Cmidshp,Cp]
        statei = copy.copy(s)
        state = goal(s)[0]
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))
        

        state = self.clean_state(state,[])
        return state
    
    def constrain_SAC_Xc(self):
        lwl     = self.lwl
        LCG     = self.LCG
        vol     = self.vol
    
        SAC_entrance_area=self.SAC_entrance_area
        SAC_mid_area=self.SAC_mid_area
        SAC_run_area=self.SAC_run_area
        
        SAC_entrance_len=self.SAC_entrance_len
        SAC_mid_len=self.SAC_mid_len
        SAC_run_len=self.SAC_run_len
        
        SAC_fwd_Xc = self.SAC_fwd_Xc
        SAC_mid_Xc = self.SAC_mid_Xc
        SAC_run_Xc = self.SAC_run_Xc
        
        
        s       = self.state
        vars = [lwl,
                LCG,
                vol,
                SAC_entrance_area,
                SAC_mid_area,
                SAC_run_area,
                SAC_entrance_len,
                SAC_mid_len,
                SAC_run_len,
                SAC_fwd_Xc,
                SAC_mid_Xc,
                SAC_run_Xc]
        statei = copy.copy(s)
        
        c1 = lp.Variable('c1')
        s.values[c1]=None
        c2 = lp.Variable('c2')
        s.values[c2]=None
        c3 = lp.Variable('c3')
        s.values[c3]=None
        
        sc1 = lp.Variable('sc1')
        s.values[sc1]=None
        sc2 = lp.Variable('sc2')
        s.values[sc2]=None
        sc3 = lp.Variable('sc3')
        s.values[sc3]=None
        sct = lp.Variable('sct')
        s.values[sct]=None
        
        a1 = lp.Variable('a1')
        s.values[a1]=None
        a2 = lp.Variable('a2')
        s.values[a2]=None
        
        
        len_entrance = self.get_val(self.SAC_entrance_len)
        len_mid = self.get_val(self.SAC_mid_len)
        #len_run = self.get_val(self.SAC_run_len)
        
        if isinstance(len_mid, ia) and isinstance(len_entrance, ia):
            g = lp.Goal.eq(SAC_mid_Xc,(.5*len_mid)+len_entrance)
            s = g(s)[0]
        
        # entrance Xc must be in the length of the entrance somewhere:
        g = lp.Goal.lto(SAC_fwd_Xc,SAC_entrance_len)
        s = g(s)[0]
        # run Xc must be in the length of the run somewhere:
        if isinstance(len_mid, ia) and isinstance(len_entrance, ia):
            xcrun_min = len_entrance+len_mid
            g = lp.Goal.gto(SAC_run_Xc,xcrun_min)
            s = g(s)[0]
        g = lp.Goal.lto(SAC_run_Xc,lwl)
        s = g(s)[0]
        
        # Sum_i ( A_i*Xc_i ) = A_tot * Xc_tot
        # SAC -area- == some real volume
        g = lp.Goal.mulo(SAC_entrance_area,SAC_fwd_Xc,sc1) #A_e*Xc_e = sc1
        s = g(s)[0]
        g = lp.Goal.mulo(SAC_mid_area,SAC_mid_Xc,sc2) #A_m*Xc_m = sc2
        s = g(s)[0]
        g = lp.Goal.mulo(SAC_run_area,SAC_run_Xc,sc3) #A_r*Xc_r = sc3
        s = g(s)[0]
        
        g = lp.Goal.mulo(vol,LCG,sct) #A_tot * Xc_tot = sct
        s = g(s)[0]
        
        g = lp.Goal.addo(sc1,sc2,c1)
        s = g(s)[0]
        g = lp.Goal.addo(c1,sc3,c2)
        s = g(s)[0]
        
        g = lp.Goal.eq(c2,sct)
        s = g(s)[0]
        
        self.set_updates(s,statei,vars)  
        state = self.clean_state(s,[c1,c2,c3,
                                    sc1,sc2,sc3,sct,
                                    a1,a2])
        return state
    
    def LCG_constraints(self):
        lwl     = self.lwl
        LCG     = self.LCG
        #lfsac   = self.lfsac
        s       = self.state

        vars = [lwl,LCG]#, lfsac]
        statei = copy.copy(s)
        #g = lp.Goal.lto(LCG,lfsac)
        #s = g(s)[0]
        g = lp.Goal.lto(LCG,lwl)
        state = g(s)[0]
        
        self.set_updates(state,statei,vars)  
        
        state = self.clean_state(state,[])
        return state
        
    def flat_of_sac_constraint(self):
        """Flat of SAC (FOSAC)
        This constraint says that
        The boxed volume availibe to the Flat of SAC curve <= vol
        
            lfsac : length of flat portion of the SAC curve
        
            relation:
                lfsac*bwl*dsmax  <= vol
            better relation:
                Amsh*lfsac <= vol
            one more:
                lfsac<lwl
        issue/DONE:
            if exists SAC_flat_min 
            then exists vol_min
            
        DONE: update to use _real_ < language!
        DONE: use Amsh instead of bwl*dsmax
        """
        lwl     = self.lwl
        lfsac   = self.lfsac
        Amsh    = self.Amsh
        SAC_mid_area = self.SAC_mid_area
        vol     = self.vol 
        SAC_mid_len = self.SAC_mid_len
        s       = self.state

        vars = [lfsac,Amsh,vol,SAC_mid_len,SAC_mid_area]
        
        statei = copy.copy(s)
        c1 = lp.Variable('c1')
        s.values[c1]=None
        c2 = lp.Variable('c2')
        s.values[c2]=None
        c3 = lp.Variable('c3')
        s.values[c3]=None
        
        # Amsh*lfsac <= vol
        g = lp.Goal.both(lp.Goal.divo(vol,Amsh,c1),
                         lp.Goal.lto(lfsac,c1) )
        s = g(s)[0]
        
        g = lp.Goal.mulo(Amsh,lfsac,c3)
        s = g(s)[0]
        g = lp.Goal.lto(c3,SAC_mid_area)
        s = g(s)[0]
        
        # lfsac<lwl
#        g = lp.Goal.both(lp.Goal.subo(lwl,lfsac,c2),
#                         lp.Goal.lto(lfsac,c2) )
        g = lp.Goal.lto(lfsac,lwl)  #the issue is that lwl is thin!!
        s = g(s)[0]
        g = lp.Goal.eq(SAC_mid_len,lfsac)
        state = g(s)[0]
        
        
        self.set_updates(state,statei,vars)  
        #del(state.values[c1])
        #del(state.values[c2])
        #del(state.values[c3])
        state = self.clean_state(state,[c1,c2,c3])
        return state
        
    def flat_of_wl_constraint(self):
        """Flat of WL (FOWL)
        This constraint says that
        The boxed area availibe to the Flat of WL curve <= Awp
        
            lfwl : length of flat portion of the waterline area curve
        
            relation:
                lfwl*bwl  <= Awp
        DONE
        """
        bwl     = self.bwl
        Awp     = self.Awp
        lfwl    = self.lfwl 
        s       = self.state
        vars = [lfwl,bwl,Awp]
        
        statei = copy.copy(s)
        c = lp.Variable('c')
        s.values[c]=None
        
        # lfwl <= Awp/bwl
        g = lp.Goal.both(
                lp.Goal.both(lp.Goal.divo(Awp,bwl,c),
                         lp.Goal.lto(lfwl,c)),
                lp.Goal.divo(Awp,bwl,c))
        
        # lfwl <= Awp/bwl
#        g = lp.Goal.both(lp.Goal.divo(Awp,bwl,c),
#                         lp.Goal.lto(lfwl,c))
        state = g(s)[0]
        
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))
        
        del(state.values[c])
        state = self.clean_state(state,[c])
        return state
        

    def flat_of_keel_constraint(self):
        """Flat of cL (FOCP)
        This constraint says that
        The boxed area availibe to the Flat of centerplane curve <= Acp
        
            lfcp : length of flat portion of the waterline area curve
        
            relation:
                lfcp*bwl  <= Acp
        DONE
        """
        draft   = self.draft
        Acp     = self.Acp
        lfcp    = self.lfcp
        s       = self.state
        vars = [lfcp,draft,Acp]
        
        statei = copy.copy(s)
        c = lp.Variable('c')
        s.values[c]=None
        
        # lfcp <= Acp/draft
        g = lp.Goal.both(
                lp.Goal.both(lp.Goal.divo(Acp,draft,c),
                         lp.Goal.lto(lfcp,c)),
                lp.Goal.divo(Acp,draft,c))
        
        # lfcp <= Acp/draft
        #        g = lp.Goal.both(lp.Goal.divo(Acp,draft,c),
        #                         lp.Goal.lto(lfcp,c))
        state = g(s)[0]
        
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))
        
        state = self.clean_state(state,[])
        return state
        
    
    def flat_of_side_constraint(self):
        """
            flat of side, FOS, is the distance:
                from : end of flat of sac
                  to : farthest edge of flat of area
        
        lfos*bwl*draft < ( vol - lfsac*Amsh )
        <--> ?same? as:
        lfos < ( lwl - lfsac )
        
        other relations, for use with AFOS:
            AFOS*bwl < Length(FOS)*draft*bwl
            AFOS < Length(FOS)*draft
            Length(FOS)*draft > AFOS        
        """
        lwl     = self.lwl
        bwl     = self.bwl
        draft   = self.draft
        lfos    = self.lfos
        vol     = self.vol
        lfsac   = self.lfsac
        Amsh    = self.Amsh
        s       = self.state
        
        vars = [bwl,draft,lfos,vol,lfsac,Amsh]
        
        statei = copy.copy(s)
        c1 = lp.Variable('c1')
        s.values[c1]=None
        c2 = lp.Variable('c2')
        s.values[c2]=None
        c3 = lp.Variable('c3')
        s.values[c3]=None
        c4 = lp.Variable('c4')
        s.values[c4]=None
        c5 = lp.Variable('c5')
        s.values[c5]=None
        
        """#ok but needed?"""
        # lfos*bwl*draft < ( vol - lfsac*Amsh )
#        g = lp.Goal.both(
#                lp.Goal.both(
#                    lp.Goal.both(lp.Goal.both(lp.Goal.mulo(bwl,draft,c1),
#                                              lp.Goal.mulo(lfsac,Amsh,c2)),
#                                 lp.Goal.both(lp.Goal.subo(vol,c2,c3),
#                                              lp.Goal.divo(c3,c1,c4))),
#                    lp.Goal.both(lp.Goal.lto(lfos,c4),
#                                 lp.Goal.both(lp.Goal.divo(c3,c1,c4),
#                                              lp.Goal.subo(vol,c2,c3)))
#                            ),
#                 lp.Goal.both(lp.Goal.mulo(lfsac,Amsh,c2),
#                              lp.Goal.mulo(bwl,draft,c1))
#                              )
#        s = g(s)[0]  #ok but needed?

 
        """#succinct but doesn't back prop"""
        """# lfos*bwl*draft < ( vol - lfsac*Amsh )"""
        g=lp.Goal.both(
                    lp.Goal.both(lp.Goal.both(lp.Goal.mulo(bwl,draft,c1),
                                              lp.Goal.mulo(lfsac,Amsh,c2)),
                                 lp.Goal.both(lp.Goal.subo(vol,c2,c3),
                                              lp.Goal.divo(c3,c1,c4))),
                    lp.Goal.lto(lfos,c4)
                    )
        s = g(s)[0] #less sprawling but doesn't back prop.
        

        
        #very important relations below 
        #(should be subsummed into relate_flats but relate flats is failing...):
        #or just re-write flats with these.
#        """# lfsac  < (lwl - lfos)"""
#        g = lp.Goal.both(lp.Goal.subo(lwl,lfos,c5),
#                         lp.Goal.lto(lfsac,c5) )
#        s = g(s)[0]
#        
        
#        """# lfos  < (lwl - lfsac)"""
        g = lp.Goal.both(lp.Goal.subo(lwl,lfsac,c5),
                         lp.Goal.lto(lfos,c5) )
        #s = g(s)[0]
                    
        state = g(s)[0]
        
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))
        
        del(state.values[c1])
        del(state.values[c2])
        del(state.values[c3])
        del(state.values[c4])
        del(state.values[c5])
        state = self.clean_state(state,[c1,c2,c3,c4,c5])
        return state
    
    def relate_drafts(self):
        draft   = self.draft
        dsmax     = self.dsmax
        s       = self.state
        vars = [draft,dsmax]
        statei = copy.copy(s)
        g = lp.Goal.lto(draft, dsmax)
        s = g(s)[0]
        self.set_updates(s,statei,vars)  
        state = self.clean_state(s,[])
        return state
        
    def relate_flats(self):
        """
            Lfsac + Lfos <= LfwL
            LfwL <= lwl
        """
        lwl     = self.lwl
        lfwl    = self.lfwl
        lfsac   = self.lfsac
        lfos    = self.lfos
        lfcp    = self.lfcp
        s       = self.state
        
        statei  = copy.copy(s)
        c = lp.Variable('c')
        s.values[c]=None
        c1 = lp.Variable('c1')
        s.values[c1]=None
        c2 = lp.Variable('c2')
        s.values[c2]=None
        
        vars = [lwl,lfwl,lfsac,lfos]
        
        g = lp.Goal.lto(lfos,lwl)
        s = g(s)[0]
        g = lp.Goal.lto(lfsac,lwl)
        s = g(s)[0]
        g = lp.Goal.lto(lfsac,lfwl)
        s = g(s)[0]
        g = lp.Goal.lto(lfcp,lfwl)
        s = g(s)[0]
        
        #TLM Dec 16, 2016
        na = ia(0.,.01)
        g = lp.Goal.gto(lfos,na)
        s = g(s)[0]
        #quick fix to make lfos >0.
        
        g = lp.Goal.both(
                lp.Goal.both(lp.Goal.addo(lfsac,lfos,c),
                         lp.Goal.lto(c,lfwl)),
                lp.Goal.addo(lfsac,lfos,c)) #propogate back down must be done locally
        s = g(s)[0]
        
#        g = lp.Goal.both(lp.Goal.addo(lfsac,lfos,c),
#                         lp.Goal.lto(c,lfwl))
#        s = g(s)[0]
                         
        """# lfsac  < (lwl - lfos)"""
        g = lp.Goal.both(lp.Goal.subo(lwl,lfos,c2),
                         lp.Goal.lto(lfsac,c2) ) #not this!
        s = g(s)[0]
        
        """# lfos  < (lwl - lfsac)"""
#        g = lp.Goal.both(lp.Goal.subo(lwl,lfsac,c2),
#                         lp.Goal.lto(lfos,c2) )
#        s = g(s)[0]
        
        
        """
        ------------------------------------
        TODO: make lfwl = lfsac + lfos
            ?at the moment it only evenly distributes?
        ------------------------------------
        """
        g = lp.Goal.lto(lfos,lfwl) 
        s = g(s)[0]
        
        state = g(s)[0]
        
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))
        
        del(state.values[c])
        del(state.values[c1])
        del(state.values[c2])
        state = self.clean_state(state,[c,c1,c2])
        return state
        
    
    def stern_profile(self):
        lwl     = self.lwl
        lfcp    = self.lfcp
        Acp    = self.Acp
        s       = self.state
        
        vars = [lwl,lfcp,Acp]
        statei = copy.copy(s)
        c1 = lp.Variable('c1')
        s.values[c1]=None
        c2 = lp.Variable('c2')
        s.values[c2]=None
        
        #g = lp.Goal.
        #chg_lst = [var for var in vars if not state(var) == statei(var)]
        #self.var_set = self.var_set.union(set(chg_lst))
        
        del(state.values[c1])
        del(state.values[c2])
#        del(state.values[c3])
#        del(state.values[c4])
#        del(state.values[c5])
        state = self.state
        state = self.clean_state(state,[c1,c2])
        return state
    
    def area_FOS_constraint(self):
        """ Area_FOS
        min(area(FOS)) * bwl  < vol_under(FOS)
        
        IMPORTANT:
            AFOS*bwl  is a min bound on vol
            
            Do not use it to restrict vol max!
        
        This constraint concerns the volume under the
        curve which defines the fwd end of the flat of side
        
        A similar curve could define the aft end
        of the flat of side as well
        """
        bwl     = self.bwl
        draft   = self.draft
        Afos    = self.Afos
        vol     = self.vol
        s       = self.state
        
        vars = [Afos,bwl,vol]
        statei = copy.copy(s)
        c1 = lp.Variable('c1')
        s.values[c1]=None
        c2 = lp.Variable('c2')
        s.values[c2]=None
        
        g = lp.Goal.both(lp.Goal.both(lp.Goal.mulo(bwl,draft,c1),
                                      lp.Goal.divo(vol,c1,c2)),
                         lp.Goal.lto(Afos,c2))
        state = g(s)[0]

        
        statei = copy.copy(s)
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))
        
        del(state.values[c1])
        del(state.values[c2])
        state = self.clean_state(state,[c1,c2])
        return state
    def SAC_mid_prop(self):
        lfsac = self.lfsac
        SAC_mid_len = self.SAC_mid_len
        SAC_mid_area = self.SAC_mid_area
        s       = self.state
        
        vars = [lfsac,SAC_mid_len,SAC_mid_area]
        statei = copy.copy(s)

        #g = lp.Goal.lto(SAC_mid_area,vol_flatsac)
        #s = g(s)[0]
        g = lp.Goal.lto(SAC_mid_len,lfsac)
        state = g(s)[0]
        
        statei = copy.copy(s)
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))
        
        
        state = self.clean_state(state,[])
        return state
    def SAC_area(self, section, sections=None):
        """quasi generic function
        to ensurethat the areas of the SAC sections
        always equal the total volume of the vessel
        """
        vol = self.vol
        SAC_entrance_area = self.SAC_entrance_area
        SAC_mid_area = self.SAC_mid_area
        SAC_run_area = self.SAC_run_area
        if sections is None:
            vars = [vol,SAC_entrance_area,
                SAC_mid_area,SAC_run_area]
                
            sections = [SAC_entrance_area,
                        SAC_mid_area,
                        SAC_run_area]
        else:
            vars = [vol]
            for el in sections:
                vars.append(el)
        index = sections.index(section)
        this_section = sections.pop(index)
        
        s       = copy.copy(self.state)
        
        statei = copy.copy(s)
        c1 = lp.Variable('c1')
        s.values[c1]=None
        c2 = lp.Variable('c2')
        s.values[c2]=None
        
        g= lp.Goal.addo(sections[0],sections[1],c1)
        s = g(s)[0]
        g = lp.Goal.subo(vol,c1,c2)
        s = g(s)[0]
        g = lp.Goal.eq(this_section,c2)
        state = g(s)[0]
        
        statei = copy.copy(s)
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))
        
        del(state.values[c1])
        del(state.values[c2])

        state = self.clean_state(state, [c1,c2])
        return state
    def SAC_len(self, section, sections=None):
        """quasi generic function
        to ensurethat the lengths of the SAC sections
        always equal th length of the lwl
        """
        lwl = self.lwl
        SAC_entrance_len = self.SAC_entrance_len
        SAC_mid_len = self.SAC_mid_len
        SAC_run_len = self.SAC_run_len
        
        if sections is None:
            vars = [lwl,SAC_entrance_len,
                SAC_mid_len,SAC_run_len]
            sections = [SAC_entrance_len,
                        SAC_mid_len,
                        SAC_run_len]
        else:
            vars = [lwl]
            for el in sections:
                vars.append(el)
                
        index = sections.index(section)
        this_section = sections.pop(index)
        
        s       = self.state
        
        statei = copy.copy(s)
        c1 = lp.Variable('c1')
        s.values[c1]=None
        c2 = lp.Variable('c2')
        s.values[c2]=None
        
        g= lp.Goal.addo(sections[0],sections[1],c1)
        s = g(s)[0]
        g = lp.Goal.subo(lwl,c1,c2)
        s = g(s)[0]
        g = lp.Goal.eq(this_section,c2)
        state = g(s)[0]
        
        statei = copy.copy(s)
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))
        
        del(state.values[c1])
        del(state.values[c2])
        state = self.clean_state(state, [c1,c2])
        return state
    
    def SAC_run_area_consistency(self):
        vol = self.vol
        SAC_entrance_area = self.SAC_entrance_area
        SAC_mid_area = self.SAC_mid_area
        SAC_run_area = self.SAC_run_area
        s       = self.state
        statei = copy.copy(s)
        
        vars = [vol,
                SAC_entrance_area,
                SAC_mid_area,
                SAC_run_area]
                
        c1 = lp.Variable('c1')
        s.values[c1]=None
        c2 = lp.Variable('c2')
        s.values[c2]=None
        
        g = lp.Goal.both(lp.Goal.addo(SAC_run_area,SAC_mid_area,c1),
                             lp.Goal.addo(SAC_entrance_area,c1,vol))
        s = g(s)[0]
        g = lp.Goal.both(lp.Goal.addo(SAC_entrance_area,SAC_mid_area,c2),
                             lp.Goal.addo(SAC_run_area,c2,vol))
        state = g(s)[0]
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))
        
        del(state.values[c1])
        del(state.values[c2])
        state = self.clean_state(state, [c1,c2])
        return state
    
    def SAC_entrance_len_consistency(self):
        lwl = self.lwl
        SAC_entrance_len = self.SAC_entrance_len
        SAC_mid_len = self.SAC_mid_len
        SAC_run_len = self.SAC_run_len
        s       = self.state
        statei = copy.copy(s)
        
        vars = [lwl,
                SAC_entrance_len,
                SAC_mid_len,
                SAC_run_len]
                
        c1 = lp.Variable('c1')
        s.values[c1]=None
        
        g = lp.Goal.both(lp.Goal.addo(SAC_entrance_len,SAC_mid_len,c1),
                             lp.Goal.addo(SAC_run_len,c1,lwl))
        s = g(s)[0]
        g = lp.Goal.both(lp.Goal.addo(SAC_run_len,SAC_mid_len,c1),
                             lp.Goal.addo(SAC_entrance_len,c1,lwl))
        state = g(s)[0]
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))
        
        del(state.values[c1])
        state = self.clean_state(state, [c1])
        return state

    def SAC_three_parts(self):
        lwl     = self.lwl
        lfsac   = self.lfsac
        vol     = self.vol
        #Cb      = self.Cb
        SAC_entrance_len = self.SAC_entrance_len
        SAC_mid_len = self.SAC_mid_len
        SAC_run_len = self.SAC_run_len
        SAC_entrance_area = self.SAC_entrance_area
        SAC_mid_area = self.SAC_mid_area
        SAC_run_area = self.SAC_run_area
        s       = self.state
        statei = copy.copy(s)
        
        vars = [lfsac,
                vol,
                lwl,
                SAC_entrance_len,
                SAC_mid_len,
                SAC_run_len,
                SAC_entrance_area,
                SAC_mid_area,
                SAC_run_area]
        statei = copy.copy(s)
        c1 = lp.Variable('c1')
        s.values[c1]=None
        c2 = lp.Variable('c2')
        s.values[c2]=None
        c3 = lp.Variable('c3')
        s.values[c3]=None
        
        s = self.SAC_mid_prop() #conistency with len!!
        #
        #-----------SAC Length Conistency
        #
        s = self.SAC_len(section=SAC_run_len,
                         sections = [SAC_entrance_len,
                                    SAC_mid_len,
                                    SAC_run_len])
        s = self.SAC_len(section=SAC_entrance_len,
                         sections = [SAC_entrance_len,
                                    SAC_mid_len,
                                    SAC_run_len])
        s = self.SAC_len(section=SAC_mid_len,
                         sections = [SAC_entrance_len,
                                    SAC_mid_len,
                                    SAC_run_len])
        
        #
        #-----------SAC volume Consistency
        #
                                    
        s = self.SAC_area(section = SAC_entrance_area,
                          sections = [SAC_entrance_area,
                                      SAC_mid_area,
                                      SAC_run_area])
        s = self.SAC_area(section = SAC_mid_area,
                          sections = [SAC_entrance_area,
                                      SAC_mid_area,
                                      SAC_run_area])
        s = self.SAC_area(section = SAC_run_area,
                          sections = [SAC_entrance_area,
                                      SAC_mid_area,
                                      SAC_run_area])
        #s = self.SAC_mid_prop()
        g = lp.Goal.both(lp.Goal.addo(SAC_entrance_len,SAC_mid_len,c1),
                        lp.Goal.addo(SAC_run_len,c1,lwl))
        s = g(s)[0]
        #g = lp.Goal.lto(c2,lwl)
        #s = g(s)[0]
        g = lp.Goal.both(lp.Goal.addo(SAC_run_len,SAC_mid_len,c1),
                        lp.Goal.addo(SAC_entrance_len,c1,lwl))
        s = g(s)[0]
        """#g = lp.Goal.eq(SAC_mid_area,Amsh*SAC_mid_len) #documentation of thinking only
        """
        
        #g = lp.Goal.eq(SAC_mid_len,lfsac)
        #s = g(s)[0]
        
        s.values[c1]=None
        s.values[c2]=None
        #g = lp.Goal.both(
        #        lp.Goal.both(lp.Goal.addo(SAC_run_area,SAC_mid_area,c1),
        #                     lp.Goal.addo(SAC_entrance_area,c1,vol)),
        #        lp.Goal.both(lp.Goal.addo(SAC_entrance_area,SAC_mid_area,c3),
        #                lp.Goal.addo(c3,SAC_run_area,vol,)))
        g = lp.Goal.both( lp.Goal.addo(SAC_run_area,SAC_mid_area,c1),
                             lp.Goal.addo(SAC_entrance_area,c1,vol))
        #assert(len(g(s))==1)
        s = g(s)[0]
        """TODO: fix consistency issue between SAC lengths and SAC areas!!!
        """
        
        s = self.SAC_entrance_len_consistency()
        #s = self.SAC_run_area_consistency() #conistency with len!!
        #s = self.SAC_mid_prop() #conistency with len!!
        
        g = lp.Goal.gto(SAC_entrance_area, ia(0.,0.))
        s = g(s)[0]
        g = lp.Goal.gto(SAC_run_area, ia(0.,0.))
        s = g(s)[0]
        
        g = lp.Goal.both(lp.Goal.addo(SAC_entrance_area,SAC_mid_area,c2),
                        lp.Goal.addo(SAC_run_area,c2,vol))
        state = g(s)[0]
        state = s
        
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))
        
        #del(state.values[c1])  ##BUG!! things are getting weird -this call fails!-
        del(state.values[c2])
        del(state.values[c3])
        state = self.clean_state(state, [c1,c2,c3])
        return state

    def bow_fairness_section_area(self):
        """
        TBD if this is chosen to be used.
        bfc   -Bow Fairness Curve-
        bhere   : breadth at the bow fairness curve bbfc
        Ahere   : area abfc
        Dhere   : Draft at station dbfc
        Cmidshp : Midship Coefficient Cbfc
        
                Abfc/(bbfc*dbfc) = Cbfc
        """
        bbfc   = self.bbfc
        draft  = self.draft
        bwl    = self.bwl
        dbfc   = self.dbfc
        Abfc   = self.Abfc
        Cbfc   = self.Cbfc
        s      = self.state
        vars   = [bbfc,dbfc,draft,Abfc,Cbfc]
        statei = copy.copy(s)
        c1     = lp.Variable('c1')
        s.values[c1]=None
        
        goal = lp.Goal.both(
                    lp.Goal.both(lp.Goal.mulo(bbfc,dbfc,c1),
                            lp.Goal.divo(Abfc,c1,Cbfc)),
                    lp.Goal.mulo(bbfc,dbfc,c1))  
        s = goal(s)[0]
        
        goal = lp.Goal.lto(Cbfc,ia(1.,1.))
        s = goal(s)[0]
        goal = lp.Goal.lto(bbfc,bwl)
        s = goal(s)[0]
        goal = lp.Goal.lto(dbfc,draft)
        state = goal(s)[0]
        
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))

        del(state.values[c1])
        state = self.clean_state(state,[c1])
        return state

    def stern_fairness_section_area(self):
        """
        TBD if this is chosen to be used.
        sfc   -Stern Fairness Curve-
        bhere   : breadth at the bow fairness curve bbfc
        Ahere   : area abfc
        Dhere   : Draft at station dbfc
        Cmidshp : Midship Coefficient Cbfc
        
                Abfc/(bbfc*dbfc) = Cbfc
        """
        bsfc   = self.bsfc
        draft  = self.draft
        bwl    = self.bwl
        dsfc   = self.dsfc
        Asfc   = self.Asfc
        Csfc   = self.Csfc
        s      = self.state
        vars   = [bsfc,dsfc,draft,Asfc,Csfc]
        statei = copy.copy(s)
        c1     = lp.Variable('c1')
        s.values[c1]=None
        
        goal = lp.Goal.both(
                    lp.Goal.both(lp.Goal.mulo(bsfc,dsfc,c1),
                            lp.Goal.divo(Asfc,c1,Csfc)),
                    lp.Goal.mulo(bsfc,dsfc,c1))  
        s = goal(s)[0]
        
        goal = lp.Goal.lto(Csfc,ia(1.,1.))
        s = goal(s)[0]
        goal = lp.Goal.lto(bsfc,bwl)
        s = goal(s)[0]
        goal = lp.Goal.lto(dsfc,draft)
        state = goal(s)[0]
        
        chg_lst = [var for var in vars if not state(var) == statei(var)]
        self.var_set = self.var_set.union(set(chg_lst))

        del(state.values[c1])
        state = self.clean_state(state,[c1])
        return state
    
    def print_state(self, var=None):
        if var is None:
            print '\n------Hull---- State-----\n'
            print '\n Principle Particulars:'
            print '{}={}'.format(self.lwl, self.state(self.lwl) )
            print '{}={}'.format(self.bwl, self.state(self.bwl) )
            print '{}={}'.format(self.draft, self.state(self.draft) )
            print '{}={}'.format(self.vol, self.state(self.vol) )
            print '{}={}'.format(self.Cp, self.state(self.Cp) )     
            print '\n Max Draft:'
            print '{}={}'.format(self.dsmax, self.state(self.dsmax) )
            print '\n Water Plane Design:'
            print '{}={}'.format(self.Awp, self.state(self.Awp) )
            print '{}={}'.format(self.Cwp, self.state(self.Cwp) )
            print '\n Keel Profile, (Center Plane) Design:'
            print '{}={}'.format(self.Acp, self.state(self.Acp) )
            print '{}={}'.format(self.Ccp, self.state(self.Ccp) )
            print '\n Midship Design:'
            print '{}={}'.format(self.Amsh, self.state(self.Amsh) )
            print '{}={}'.format(self.Cmidshp, self.state(self.Cmidshp) )
            print '{}={}'.format(self.Cb, self.state(self.Cb) )
            print '\n Design of Flats:'
            print '{}={}'.format(self.lfwl, self.state(self.lfwl) )
            print '{}={}'.format(self.lfos, self.state(self.lfos) )
            print '{}={}'.format(self.lfcp, self.state(self.lfcp) )
            print '\n bow fairness curve:'
            print '{}={}'.format(self.bbfc, self.state(self.bbfc) )
            print '{}={}'.format(self.dbfc, self.state(self.dbfc) )
            print '{}={}'.format(self.Abfc, self.state(self.Abfc) )
            print '{}={}'.format(self.Cbfc, self.state(self.Cbfc) )
            print '\n stern fairness curve:'
            print '{}={}'.format(self.bsfc, self.state(self.bsfc) )
            print '{}={}'.format(self.dsfc, self.state(self.dsfc) )
            print '{}={}'.format(self.Asfc, self.state(self.Asfc) )
            print '{}={}'.format(self.Csfc, self.state(self.Csfc) )
            print '\n Destailed SAC (3X3-parts):'
            print '{}={}'.format(self.SAC_entrance_len, 
                                self.state(self.SAC_entrance_len) )
            print '{}={}'.format(self.SAC_mid_len, 
                                self.state(self.SAC_mid_len) )
            print '{}={}'.format(self.SAC_run_len, 
                                self.state(self.SAC_run_len) )
            print '\n{}={}'.format(self.SAC_entrance_area, 
                                self.state(self.SAC_entrance_area) )
            print '{}={}'.format(self.SAC_mid_area, 
                                self.state(self.SAC_mid_area) )
            print '{}={}'.format(self.SAC_run_area, 
                                self.state(self.SAC_run_area) )
            print '\n{}={}'.format(self.SAC_fwd_Xc, 
                                self.state(self.SAC_fwd_Xc) )
            print '{}={}'.format(self.SAC_mid_Xc, 
                                self.state(self.SAC_mid_Xc) )
            print '{}={}'.format(self.SAC_run_Xc, 
                                self.state(self.SAC_run_Xc) )
            print '{}={}'.format(self.LCG, self.state(self.LCG) )
            print '{}={}'.format(self.lfsac, self.state(self.lfsac) )
            print '\n--------------------\n'
        else:
            if var in self.state.values:
                print '{}=>{}'.format(var, self.state(var))
            else:
                try:
                    key = self.__dict__['_'+str(var)]
                    print '{}=>{}'.format(key, self.state(key))
                except:
                    print 'Var not found'
        return
        
def two(func):
    def fwrap():
        return '2', func()
    return fwrap

@two
def one():
    return '1!'        
if __name__=='__main__':

    
    a = ia(0.,100.)
    b = ia(50.,150.)
    c = a&b
    
    
    
    self = Hull()
    
    #goal = lp.Goal.eq(self.lwl, ia(2000.,300000.))
    #self.state = goal(self.state)[0]
    self.vol = ia(2000.,300000.)
    self.lwl = ia(120.,120.)
    self.bwl = ia(20.,20.) 
    
#    self.lfwl = ia(8.,100.) #need to set eary to enforce <=  !!!
#    self.lfsac = ia(1.,15.)
#    self.lfos = ia(5.,9.)
#    #self.lfcp = ia(0.,1000.)
#
#    #self.Cb = ia(.5,.7) 
#    self.draft = ia(20.,20.) 
#    self.Cwp = ia(.2,.4) 
#    self.Cmidshp =  ia(.7,.99) 
#    self.dsmax =  ia(18.,20.)
#    
#    #self.Acp=ia(0.,1000.)
#    self.Ccp = ia(.2,.4)
#    print self.print_state()
#    
#    
#    self.compute_prismatic_coefficient()
#    print self.print_state()
#    
#    self.compute_Cb()
#    print self.state
#    
#    #self.compute_Cb(**{'Cb':ia(.5,.55)})
#    self.compute_Cb()
#    print self.state
#    
#    self.compute_Cwp()
#    print self.state
#    
#    self.compute_midship_coefficient()
#    print self.state
#    
#    self.compute_prismatic_coefficient()
#    print self.state
#    
#    
#
#    self.print_state()
#    print '\nNow sanity checks\n'
#    print self.var_set
#    print ''
#    self.print_state('Cb')
#    self.print_state('Cp')
#    print 'now change Cp and Cb'
#    self.Cp = ia(0.,1.)
#    self.Cb = ia(0.,1.)
#    print ''
#    self.print_state('Cb')
#    self.print_state('Cp')
#    print '\nvariables that have changed:\n',self.var_set
#    self.print_state()
#    #self.AC_revise(print_=True)
#    self.print_state()
#    
#    #c1 = lp.Variable('c1')
#    #self.state.assign(c1, ia(1.,2.))
#    
#
#    
#    #self.AC_revise()
#    md = self.state(self.vol).getpoint(.5)
#    
#    #self.AC_revise()
#    #self.arc_set = set([])
#    #self.vol = ia(md,md)
#    #self.vol = ia(md-.01,md+.01)
#    
#    self.Cb = ia(.49,.51)
#    self.print_state()
#    #self.AC_revise()
#    self.Cb.arcs[0]()
#    self.Cb.arcs[1]()
#    self.print_state()
#    
#    #g = lp.Goal.eq(self.vol,ia(md,md))
#    #s = g(self.state)[0]
#    #self.AC_revise(print_=True)
#    
#    
#    #self.lfsac = ia(20.,76.)
#    
#    #bow fairness curve sectional area info:
#    self.Cbfc = ia(.1,1.)
#    self.Abfc = self.state.value_of(self.Amsh)/(2.*2.5)
#    self.bbfc = self.state.value_of(self.bwl)/(2.*2.5)
#    self.print_state()
#    
#    #stern fairness curve sectional area info:
#    self.Csfc = ia(.1,1.)
#    self.Asfc = self.state.value_of(self.Amsh)/(2.*2.5)
#    self.bsfc = self.state.value_of(self.bwl)/(2.*2.5)
#    
#
#
#    #section SAC properties
#    # lengths:
#    self.state = self.SAC_mid_prop()
#    x=self.state(self.lfsac).getpoint(.9)
#    self.SAC_mid_len = ia(x,x)
#    re = (120.0 - x)/2.
#    #self.SAC_entrance_len = ia(re,re)
#    self.SAC_run_len = ia(re,re)    
#    #self.state = self.SAC_mid_prop()
#    # volumes:
#    v = self.state(self.vol)
#    sma = self.state(self.SAC_mid_area) #SAC area == Volume
#    re = (v-sma)/2.
#    self.SAC_entrance_area = re
#    #c1 = lp.Variable('c1')
#    #g = lp.Goal.addo(self.SAC_entrance_area,self.SAC_run_area,c1)
#    #s=g(self.state)[0]
#    #g = lp.Goal.addo(self.SAC_entrance_area,self.SAC_run_area,c1)
#    #s=g(self.state)[0]
#    #-------------------------#
#    self.print_state()
"""    
lwl=ia(120.0,120.0)
bwl=ia(20.0,20.0)
draft=ia(20.0,20.0)
vol=ia(23520.0,24480.0)
Cb=ia(0.49,0.509999999999)
dsmax=ia(18.0,20.0)
Awp=ia(480.0,960.0)
Cwp=ia(0.2,0.4)
Amsh=ia(252.0,396.0)
Cmidshp=ia(0.7,0.99)
Cp=ia(0.494949494949,0.728571428569)
lfsac=ia(0.0,100.0)
lfwl=ia(0.0,100.0)
lfos=ia(0.0,100.0)

self = Hull()
self.lwl=ia(120.0,120.0)
self.bwl=ia(20.0,20.0)
self.draft=ia(20.0,20.0)
self.vol=ia(23520.0,24480.0)
self.Cb=ia(0.49,0.509999999999)
self.dsmax=ia(18.0,20.0)
self.Awp=ia(480.0,960.0)
self.Cwp=ia(0.2,0.4)
self.Amsh=ia(252.0,396.0)
self.Cmidshp=ia(0.7,0.99)
self.Cp=ia(0.494949494949,0.728571428569)
self.lfsac=ia(0.0,300.0)
self.lfwl=ia(0.0,300.0)
self.lfos=ia(0.0,300.0)

lfsac=ia(0.0,97.1428571428)
lfwl=ia(0.0,48.0)
lfos=ia(0.0,61.2)
#"""