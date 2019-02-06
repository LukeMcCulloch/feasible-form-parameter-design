# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 22:41:18 2016

@author: lukemcculloch

CONTAINS:  Not Quite BASIC NAVAL ARCH RULES

This file is from before I simplified life by following the FPD hierarchy
here we are still able to design, e.g. transverse curves
such as the bow and stern fairness curves
Up at the same level as the hull form curves.
This is not recommended practice at this time.  -TLM, April 2017
"""

#import tokenize
#import numpy as np
import copy
#from automatic_differentiation import ad
#from interval_arithmetic import ia
from extended_interval_arithmetic import ia
#from hull_inference_ob_graph import Hull as hullclp

#import uKanren as lp #original
#import eKanren as lp#nicer! NOT backwards compatible
import sqKanren as lp

#meta?
#from functools import wraps #preserve metadata if I start decorating heavily...
#from inspect import signature #not availible

"""
NOTES:
    -breadths and beams are full beams.  
    -Areas and area ratios use breadth as such

TODO: 
    -Done: accessor on states which returns a map from states to valsues
    for that particular parameter.  Actually this lives in states.
    
"""

class C3PartCurve(object):
    def __init__(self, cls_st,
                 total_length, l1,l2,l3,
                 total_area, a1,a2,a3,
                 Xc, x1,x2,x3):
        self.cls_st         = cls_st
        self.total_length   = total_length
        self.l1             = l1
        self.l2             = l2
        self.l3             = l3
        self.dl1 = lp.Variable('dl1')
        self.dl2 = lp.Variable('dl2')
        self.total_area     = total_area
        self.a1             = a1
        self.a2             = a2
        self.a3             = a3
        self.da1 = lp.Variable('da1')
        self.da2 = lp.Variable('da2')
        self.Xc             = Xc
        self.x1             = x1
        self.x2             = x2
        self.x3             = x3
        self.dx1 = lp.Variable('dx1')
        self.dx2 = lp.Variable('dx2')
        
    
    def C3_curve(self, aggregate_quant, section, sections,
                 f1=ia(.3,.37),f2=ia(.3,.37), niter=3):
        """ generic 3 part curve function: called with any section
        to ensure that the __quantities__ of the 3 sections
        always equal the total quantity of the vessel
        """
        vars = [aggregate_quant]
        for el in sections:
            vars.append(el)
        index = sections.index(section)
        this_section = sections.pop(index)
        
        s = self.cls_st.states
        statesi = copy.copy(s)
        
        if aggregate_quant == self.total_length:
            #print 'found len = ', self.total_length.name
            c1 = self.dl1
            c2 = self.dl2
        elif aggregate_quant == self.total_area:
            #print 'found area = ', self.total_area.name
            c1 = self.da1
            c2 = self.da2
        elif aggregate_quant == self.Xc:
            #print 'found Xc = ', self.Xc.name
            c1 = self.dx1
            c2 = self.dx2
        else:
            print 'Error:3 part curve: Could Not Find aggreagate_quant!'
        
        dlist = []
        
        states = (self.cls_st + (sections[0],sections[1],c1))
        for i in range(niter):
            states = (states - (aggregate_quant,c1,c2))
            states = (states == (this_section,c2))
            states = (states + (sections[0],sections[1],c1))
        
        
        self.cls_st.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars,
                         dlist=dlist)
        states = self.cls_st.clean_states(states.states,dlist)
        return states

        
        
def list_setter(design,x,designlist):
    """
        if using coefficients list
        be sure and initialize to [0.,1.]
        to set sesible vars here!
    """
    #print designlist
    for el in designlist:
        mvar = design(el)[0]
        #print el, mvar
        if isinstance(mvar, ia):
            val = mvar.getpoint(x)
            #print 'int', val
            val = ia(val,val)
            design.__setattr__(el.name,val)
        else:
            print 'mavar is a variale'
    return
        
def set_lots():
    x=1.
    list_setter(self, x, self.Coefficients)
    list_setter(self, x, self.Primary)
    list_setter(self, x, self.BowCurveCoeff)
    list_setter(self, x, self.SternCurveCoeff)
    #list_setter(self, x, self.Areas)
    #list_setter(self, x, self.list_SAC)
    return
    


class DList(list):
    def __init__(self):
        pass
    def getpoint(self, pt):
        return self[0].getpoint(pt)


#class DesignTree(object):
#    def __init__(self, parent=None, children=None):
#        self.parent = parent
#        if children is None:
#            self.children = []
#        else:
#            self.children = children
            
#class DesignTree(list):
class DesignTree(lp.States):
    """Holds the tree of design states
    """
    #def __init__(self, states, parent=None, children=None):
    #list.__init__(self, states)
    #        self.parent = parent
    #        if children is None:
    #            self._children = []
    #        else:
    #            self.children = children
        
    def get_latest_child(self):
        lc = len(self.children)
        #nst = len(self)
        if lc == 0:
            return lp.States(self.states)#, nst
            #return lp.States(self[:])#, nst
        else:
            return self.get_child()
    
    def __call__(self, x):
        """we could inherit this call from States
        """
        d = []
        #for st in self: #inherit from list
        for st in self.states:
            d.append(st(x))
        return d
    
    def __len__(self):
        return len(self.states)
    
            
class CompoundCurve(object):
    def __init__(self):
        pass         
            
class HullDesignNode(lp.States):
    """
        Key Insight: This class should intialize the design space
        and pass that initialized state off to a design tree
            -desisions made in that tree act on [state] that is fed
                back into this class
            -But new states are sent from here to the design tree
            -Which can do what it wants with them
    """
    SMALL = 1.e-3
    def __init__(self, strd_coeff=True, parent=None):
        self.nstates = 0
        self._verbose = False
        self.valid_designs = True
        self._states = None
        #super(HullDesignNode, self).__init__( states = DesignTree([self._setup_()]))
        super(HullDesignNode, self).__init__( states = [self._setup_()] )
        self.nstates = len(self.states)
        self.keys = [el for el in self.states[0].values.keys() 
                    if isinstance(el, lp.Variable)]
        self.make_lists()
        self.parent = parent
        if strd_coeff:
            self.ini_coeff()
            self.ini_lengths()
            self.ini_areas()
            self.ini_locations()
            self.compute_stern_fairness_section_ini()
        self.SAC_rules = C3PartCurve(self,
                                     self.lwl, 
                                     self.SAC_entrance_len,
                                     self.SAC_mid_len,
                                     self.SAC_run_len,
                                     self.vol,
                                     self.SAC_entrance_area,
                                     self.SAC_mid_area,
                                     self.SAC_run_area,
                                     self.LCG,
                                     self.SAC_fwd_Xc,
                                     self.SAC_mid_Xc,
                                     self.SAC_run_Xc)
            

    def print_coeff(self):
        for key in self.Coefficients:
            print key, self(key)
        return
    def print_primary(self):
        for key in self.Primary:
            print key, self(key)
        return
    def print_coeff_bow_curve(self):
        for key in self.BowCurveCoeff:
            print key, self(key)
        return
    def print_coeff_stern_curve(self):
        for key in self.SternCurveCoeff:
            print key, self(key)
        return
        
    def print_SAC_stats(self):
        for key in self.list_SAC:
            print key, self(key)
        return
        
    def print_areas(self):
        for key in self.Areas:
            print key, self(key)
        return
        
    def make_lists(self):
        
        self.hull_list = [self.vol,
                          #self.Cb,
                          self.lwl,
                          self.bwl,
                          self.draft,
                          self.LCG]
                          
        self.Coefficients = [self.Cb,
                             self.Cmidshp,
                             self.Cp,
                             self.Cwp,
                             self.Ccp,
                             self.Cbfc,
                             self.Csfc]
        self.Areas = [self.Amsh,
                      self.Awp,
                      self.Acp]#,
                      #self.Abfc,
                      #self.Asfc,
                      #self.SAC_entrance_area,
                      #self.SAC_mid_area,
                      #self.SAC_run_area]
                      
        self.Primary = [self.Cb,
                        self.lwl,
                        self.vol,
                        self.draft,
                        self.bwl,
                        self.LCG]#,
                        #self.Cb,
                        #self.Cwp,
                        #self.Ccp,
                        #self.LCG,
                        #self.Cmidshp]
        
        self.BowCurveCoeff = [self.Abfc,
                              self.Cbfc,
                              self.bbfc,
                              self.dbfc]
        
        self.SternCurveCoeff = [self.Asfc,
                                self.Csfc,
                                self.bsfc,
                                self.dsfc]
        self.list_SAC = [self.SAC_entrance_area,
                         self.SAC_run_area,
                         self.SAC_mid_area,
                         self.SAC_entrance_len,
                         self.SAC_run_len,
                         self.SAC_mid_len,
                         self.SAC_fwd_Xc,
                         self.SAC_mid_Xc,
                         self.SAC_run_Xc]
        
        self.glue = [self.loc_bfc,self.loc_sfc]
        self.alists=self.Coefficients +\
                    self.Areas+\
                    self.Primary+\
                    self.BowCurveCoeff+\
                    self.SternCurveCoeff+\
                    self.list_SAC+\
                    self.glue
        return
        
    def _setup_(self):
        """Set up the design variables
        of the design space
        """
                                    
        #
        #------------------------- Hull Primary
        #
        self._lwl     = lp.Variable('lwl',
                                    arcs=[self.compute_Cb,
                                          self.compute_Cwp,
                                          self.compute_Ccp,
                                          self.compute_flat_relations,
                                          self.compute_SAC_section_properties,
                                          self.compute_LCG,
                                          self.ini_lengths])
                                          
        self._bwl     = lp.Variable('bwl',
                                    arcs=[self.compute_Cb,
                                          self.compute_Cwp,
                                          self.compute_midship_coefficient,
                                          self.compute_FOWL,
                                          self.compute_FOS,
                                          self.derive_bow_fairness_section_loc_rule_fromDWL,
                                          self.derive_stern_fairness_section_loc_rule_fromDWL])
        self._draft   = lp.Variable('draft',
                                    arcs=[self.compute_Cb,
                                          self.compute_Ccp,
                                          self.compute_drafts,
                                          self.compute_FOK,
                                          self.compute_FOS,
                                          self.compute_midship_coefficient,
                                          self.derive_bow_fairness_section_loc_rule_fromCPKeel,
                                          self.derive_stern_fairness_section_loc_rule_fromCPKeel])
        
        self._dsmax   = lp.Variable('dsmax',
                                    arcs=[self.compute_midship_coefficient,
                                          self.compute_drafts])
        
        self._vol     = lp.Variable('vol',
                                    arcs=[self.compute_Cb,
                                          self.compute_SAC_section_properties,
                                          self.compute_FOS,
                                          self.compute_FOSAC,
                                          self.ini_areas,
                                          self.ini_lengths,])
                                    
        self._Cb      = lp.Variable('Cb',
                                    arcs=[self.compute_Cb,
                                         self.compute_prismatic_coefficient])
        #
        #-------------------------SAC LCG
        #         
        self._LCG      = lp.Variable('LCG',
                                    arcs=[self.compute_LCG,
                                          self.compute_SAC_section_properties,
                                          self.sac_XC_ini])   
        self._SAC_fwd_Xc = lp.Variable('SAC_fwd_Xc',
                                    arcs=[self.compute_LCG,
                                          self.sac_XC_ini,
                                          self.ini_lengths,
                                          self.ini_areas])   
        self._SAC_mid_Xc = lp.Variable('SAC_mid_Xc',
                                    arcs=[self.compute_LCG,
                                          self.sac_XC_ini,
                                          self.ini_lengths,
                                          self.ini_areas])   
        self._SAC_run_Xc = lp.Variable('SAC_run_Xc',
                                    arcs=[self.compute_LCG,
                                          self.sac_XC_ini,
                                          self.ini_lengths,
                                          self.ini_areas])
        
        #
        #------------------------- waterplane area
        #
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
                                          self.compute_SAC_section_properties,
                                          self.compute_FOS,
                                          self.compute_FOSAC,
                                          self.derive_bow_fairness_section_loc_rule_fromSAC,
                                          self.derive_stern_fairness_section_loc_rule_fromSAC])
        
        #midship coefficient
        self._Cmidshp    = lp.Variable('Cmidshp',
                                    arcs=[self.compute_midship_coefficient,
                                          self.compute_prismatic_coefficient]) 
        #prismatic coefficient
        self._Cp      = lp.Variable('Cp',
                                    arcs=[self.compute_prismatic_coefficient]) 
        
        
        self._lfwl   = lp.Variable('lfwl',
                                    arcs=[self.compute_FOWL,
                                         self.compute_flat_relations])
        
        self._lfos   = lp.Variable('lfos',
                                    arcs=[self.compute_FOS,
                                         self.compute_flat_relations])
        
        #self._Afos   = lp.Variable('Afos',
        #                            arcs=[self.compute_FOS])
                                    
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
        self._loc_bfc = lp.Variable('loc_bfc',
                                    arcs=[self.derive_bow_fairness_section_loc_rule_fromSAC,
                                          self.derive_bow_fairness_section_loc_rule_fromDWL,
                                          self.derive_bow_fairness_section_loc_rule_fromCPKeel])
        self._bbfc  = lp.Variable('bbfc',
                                   arcs=[self.compute_bow_fairness_section_area,
                                         self.derive_bow_fairness_section_loc_rule_fromDWL])
        self._dbfc  = lp.Variable('dbfc',
                                   arcs=[self.compute_bow_fairness_section_area,
                                         self.derive_bow_fairness_section_loc_rule_fromCPKeel])
        self._Abfc  = lp.Variable('Abfc',
                                   arcs=[self.compute_bow_fairness_section_area,
                                         self.derive_bow_fairness_section_loc_rule_fromSAC])
        self._Cbfc  = lp.Variable('Cbfc',
                                   arcs=[self.compute_bow_fairness_section_area])
        ##
        ## Stern Fairness curve
        ##
        """TBD"""
        self._loc_sfc = lp.Variable('loc_sfc',
                                    arcs=[self.derive_stern_fairness_section_loc_rule_fromSAC,
                                          self.derive_stern_fairness_section_loc_rule_fromDWL,
                                          self.derive_stern_fairness_section_loc_rule_fromCPKeel])
        self._bsfc  = lp.Variable('bsfc',
                                   arcs=[self.compute_stern_fairness_section_area,
                                         self.derive_stern_fairness_section_loc_rule_fromDWL])
        self._dsfc  = lp.Variable('dsfc',
                                   arcs=[self.compute_stern_fairness_section_area,
                                         self.derive_stern_fairness_section_loc_rule_fromCPKeel])
        self._Asfc  = lp.Variable('Asfc',
                                   arcs=[self.compute_stern_fairness_section_area,
                                         self.derive_stern_fairness_section_loc_rule_fromSAC])
        self._Csfc  = lp.Variable('Csfc',
                                   arcs=[self.compute_stern_fairness_section_area])
        ##
        ## SAC length
        ##
        self._SAC_entrance_len = lp.Variable('SAC_entrance_len',
                                             arcs = [self.compute_SAC_section_properties,
                                                     self.ini_lengths,
                                                     self.ini_areas,
                                                     self.derive_bow_fairness_section_loc_rule_fromSAC,
                                                     self.derive_bow_fairness_section_loc_rule_fromDWL,
                                                     self.derive_bow_fairness_section_loc_rule_fromCPKeel])
        self._SAC_mid_len = lp.Variable('SAC_mid_len',
                                             arcs = [self.compute_SAC_section_properties,
                                                     self.ini_lengths,
                                                     self.ini_areas])
        self._SAC_run_len = lp.Variable('SAC_run_len',
                                             arcs = [self.compute_SAC_section_properties,
                                                     self.ini_lengths,
                                                     self.ini_areas,
                                                     self.derive_stern_fairness_section_loc_rule_fromSAC,
                                                     self.derive_stern_fairness_section_loc_rule_fromDWL,
                                                     self.derive_stern_fairness_section_loc_rule_fromCPKeel])
        ##
        ## SAC area
        ##
        self._SAC_entrance_area = lp.Variable('SAC_entrance_area',
                                             arcs = [self.compute_SAC_section_properties,
                                                     self.ini_areas,
                                                     self.ini_lengths])
        self._SAC_mid_area = lp.Variable('SAC_mid_area',
                                             arcs = [self.compute_SAC_section_properties,
                                                     self.ini_areas,
                                                     self.ini_lengths])
        self._SAC_run_area = lp.Variable('SAC_run_area',
                                             arcs = [self.compute_SAC_section_properties,
                                                     self.ini_areas,
                                                     self.ini_lengths])
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
                              #self._Afos    : None,
                              self._loc_bfc  : None,
                              self._bbfc     : None,
                              self._dbfc     : None,
                              self._Abfc     : None,
                              self._Cbfc     : None,
                              self._loc_sfc  : None,
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
        #self.arc_sets = [set([])] #active arcs
        #self.var_sets = [set([])] #active variables 
        self.arc_sets = set([]) #active arcs
        self.var_sets = set([]) #active variables 
        return
    
    def get_arc_sets(self):
        """Find the constraints that have inputs which changed
        """
        #lac = len(self.arc_sets)
        #lvs = len(self.var_sets)
        #for c, var_set in enumerate(self.var_sets):
            #if c>=lac:
            #    self.arc_sets.append(set())
        for var in self.var_sets:
                #self.arc_sets[c] = self.arc_sets[c].union(set(var.arcs))
            self.arc_sets = self.arc_sets.union(set(var.arcs))
#            else:
#                for var in var_set:
#                    self.arc_sets[c] = self.arc_sets[c].union(set(var.arcs))
        return
    
    def AC_revise(self, print_=False, maxit = 15):
        """Arc Consistency for the Hull Parameters
        """
        self.get_arc_sets()
        #if print_:
        #    print 'arcs of interest:\n',self.arc_sets
        #self.var_sets = [set([])]
        self.var_sets = set([])
        
        #for c, arc_set in enumerate(self.arc_sets):
        count=0
        while len(self.arc_sets)>0 and count<maxit:
            for arc in self.arc_sets:
                if print_ or self._verbose: print '\ndoing ', arc
                arc()
            self.arc_sets = set([])
            self.get_arc_sets()
            self.var_sets = set([])
            #            self.arc_sets[c] = set([])
            #            self.get_arc_sets()
            #            self.var_sets[c] = set([])
            count += 1
            print count
            
        return
    
    def equal(self, c1,c2):
        """Updated lp.Goal.eq, 
        with 
        hull constraint graph propogation enabled
        
        How to handle when the both are lp.Variable ?
        answer, use the more general lp.eq(anything, anything)
        
        The purpose of this function is to 
        add self.c1 to the set of vars to 
        trigger Arc Consistency
        
        
        var set to be handle by hull class
        """
        
        l1 = len(self.states)
        statesi = copy.copy(self.states)
        #
        assert(l1==self.nstates),'nstates not equal nstates?'
        #before = [st.value_of(c1) for st in self.states]
        # unify
        si = (self == (c1, c2))
        self.states = si.states
        #after = [st.value_of(c1) for st in self.states]
        self.set_updates(states=self.states, 
                         states_old = statesi,
                         vars=[c1])
        self.AC_revise()
        #self.delete_redundant_states()
        if self.nstates == 0:
            self.valid_designs = False
        return
    
    def clean_states(self, states, special=None):
        """TODO: Update for many states...
        Launguage is not quite up to par yet
        making do, clean extranious vars thatcrop up
        during rule computations
        """
        return states
        for state in states:
            snm = [el.name for el in special]
            dkeys = []
            for key in state.values:
                #if not isinstance(key, lp.Variable):
                #    dkeys.append(key)
                if key.name in snm:
                    dkeys.append(key)
            for key in dkeys:
                del(state.values[key])
        return states
        
    def set_updates(self, states, states_old, vars, dlist=None, tol=.001):
        """TODO: really make this about changes 
        from old state 
        to new states(s)
        
        
        THIS IS GOING TO GET OUT OF WACK!!
        TODO??
        """
        l1 = len(states) #new
        l2 = len(states_old) #old
        #iters = np.
        #if l1==l2:
        c = 0
        if dlist is None:
            for state, statei in zip(states, states_old):
                chg_lst = [var for var in vars if not state(var) == statei(var)]
                #chg_lst = [var for var in vars if not self.diff( ( state(var), statei(var),tol))]
                var_set = self.var_sets.union(set(chg_lst)) 
                self.var_sets = var_set
        else:
            nlist = [el.name for el in dlist]
            var_set = set()
            for state, statei in zip(states, states_old):
                chg_lst = [var for var in vars if (not state(var) == statei(var) and not (state(var).name  in nlist))]             
                #var_set = var_set.union(set(chg_lst)) 
                var_set = self.var_sets.union(set(chg_lst)) 
                self.var_sets = var_set
                #--
                #var_set = self.var_sets[c].union(set(chg_lst)) 
                #self.var_sets[c] = var_set
        #    c+=1
#        else:
#            for c, state in enumerate(states):
#                chg_lst = vars
#                if c>=l2:
#                    #self.var_sets = list(self.var_sets)
#                    var_set = set(chg_lst)
#                    self.var_sets.append(set(var_set))
#                    #self.var_sets = set(self.var_sets)
#                else:
#                    var_set = self.var_sets[c].union(set(chg_lst)) 
#                    self.var_sets[c] = var_set
        return
    
    def set_val(self, var, val):
        for state in self.states:
            state.values[var] = val
        return 
    
    def get_val(self, var):
        return self(var)
    
    
    def equality_check(self):
        n = self.nstates
        equal_states = {}
        for i in range(n):
            equal_states[i]=[]
            for j in range(n):
                if i != j:
                    if self.equality_check_one(self.states[i],
                                               self.states[j]):
                        equal_states[i].append(j)
        return equal_states
        
    def equality_check_one(self,st1,st2):
        equality = True
        for k1 in self.keys:
            
            if st1(k1) != st2(k1):
                equality = False
                if self._verbose:
                    print k1, st1(k1)
                    print k1, st2(k1)
        return equality
    
    def delete_redundant_states(self):
        eq = self.equality_check()
        save = [0]
        delete = []
        for key in eq:
            if key not in save and key not in delete:
                save.append(key)
            for el in eq[key]:
                if el in save:
                    pass
                else:
                   delete.append(el) 
        delete = list(set(delete))
        delete.reverse()
        
        if self._verbose:
            print delete
        for el in delete:
            if self._verbose:
                print 'deleting redundant state {}'.format(el)
            del self.states[el]
            self.nstates -= 1
        return
        
    @property
    def states(self):
        return self._states
    @states.setter
    def states(self, states):
        l1 = self.nstates
        self._states = states#DesignTree(states)
        l2 = len(states)
        if l1 != l2:
            self.nstates=l2
            self.delete_redundant_states()
        
    #    @property
    #    def children(self):
    #        if self._children is None:
    #            self._children = []
    #        return self._children
    #    @children.setter
    #    def children(self, states):
    #        states.parent = self._states
    #        self._children.append(states)
        
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
        
    #    @property
    #    def Afos(self):
    #        return self._Afos
    #    @Afos.setter
    #    def Afos(self, Afos):
    #        self.equal(self.Afos,Afos)
        
    ##
    ## Bow Fairness curve
    ##
    @property
    def loc_bfc(self):
        return self._loc_bfc
    @loc_bfc.setter
    def loc_bfc(self, loc_bfc):
        self.equal(self._loc_bfc,loc_bfc)
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
    def loc_sfc(self):
        return self._loc_sfc
    @loc_sfc.setter
    def loc_sfc(self, loc_sfc):
        self.equal(self._loc_sfc,loc_sfc)
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
        self.states = self.flat_of_sac_constraint()
        return
    def compute_LCG(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.LCG_constraints()
        return
    
    def compute_drafts(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.relate_drafts()
        return
        
    def compute_FOWL(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.flat_of_wl_constraint()
        return
        
    def compute_FOK(self, **args):
        """Flat of Keel
        """
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.flat_of_keel_constraint()
        return
        
    def compute_FOS(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.flat_of_side_constraint()
        return
        
        
    def compute_Cb(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.block_coefficient()
        return
    
    def compute_Cwp(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.waterplane_coefficient()
        return
    
    def compute_Ccp(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.center_plane_coefficient()
        return
        
    def compute_midship_coefficient(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.midship_coefficient()
        return
        
    def compute_prismatic_coefficient(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.prismatic_coefficient()
        return
    
    def compute_flat_relations(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.relate_flats()
        return
    
    def compute_bow_fairness_section_area(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.bow_fairness_section_area()
        return
        
    
    def derive_bow_fairness_section_loc_rule_fromSAC(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.bow_fairness_curve_location_from_SAC()
        return
        
    
    def derive_bow_fairness_section_loc_rule_fromDWL(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.bow_fairness_curve_location_from_DWL()
        return
        
    
    def derive_bow_fairness_section_loc_rule_fromCPKeel(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.bow_fairness_curve_location_from_CPKeel()
        return
    
    
    
    def compute_bow_fairness_section_ini(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.shape_bow_section_inizers()
        return
    
    def compute_stern_fairness_section_ini(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.shape_stern_section_inizers()
        return
        
    
    def compute_stern_fairness_section_area(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.stern_fairness_section_area()
        return
    
    def derive_stern_fairness_section_loc_rule_fromSAC(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.stern_fairness_curve_location_from_SAC()
        return
    
    def derive_stern_fairness_section_loc_rule_fromDWL(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.stern_fairness_curve_location_from_DWL()
        return
    
    def derive_stern_fairness_section_loc_rule_fromCPKeel(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.stern_fairness_curve_location_from_CPKeel()
        return
    
    def compute_SAC_section_properties(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        #self.states = self.SAC_mid_prop()
        #self.states = self.SAC_three_parts() #SAC_entrance_area
        self.states = self.SAC_rules.C3_curve(self.lwl,
                                              section = self.SAC_entrance_len,
                                              sections = [  self.SAC_entrance_len,
                                                            self.SAC_mid_len,
                                                            self.SAC_run_len])
        self.states = self.SAC_rules.C3_curve(self.vol,
                                              section = self.SAC_entrance_area,
                                              sections = [  self.SAC_entrance_area,
                                                            self.SAC_mid_area,
                                                            self.SAC_run_area])
    #        self.states = self.SAC_rules.C3_curve(self.lwl,
    #                                              section = self.SAC_mid_len,
    #                                              sections = [  self.SAC_entrance_len,
    #                                                            self.SAC_mid_len,
    #                                                            self.SAC_run_len])
    #        self.states = self.SAC_rules.C3_curve(self.vol,
    #                                              section = self.SAC_mid_area,
    #                                              sections = [  self.SAC_entrance_area,
    #                                                            self.SAC_mid_area,
    #                                                            self.SAC_run_area])
#        self.states = self.SAC_rules.C3_curve(self.lwl,
#                                              section = self.SAC_run_len,
#                                              sections = [  self.SAC_entrance_len,
#                                                            self.SAC_mid_len,
#                                                            self.SAC_run_len])
    #        self.states = self.SAC_rules.C3_curve(self.vol,
    #                                              section = self.SAC_run_area,
    #                                              sections = [  self.SAC_entrance_area,
    #                                                            self.SAC_mid_area,
    #                                                            self.SAC_run_area])
        self.states = self.SAC_run_area_consistency()
        #self.states = self.SAC_entrance_len_consistency()
        self.states = self.constrain_SAC_Xc() #SAC_run_area
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
    
    def ini_coeff(self):
        """
            make sure all coefficient are in [0.,1.].
        """
        for key in self.Coefficients:
            self.__setattr__(key.name,ia(0.,1.))
        self.AC_revise()
        return
        
    def ini_lengths(self):
        """
            make sure all areas are >= 0.
        """
        self.states = self.SAC_length_ini()
        return
        
    def ini_areas(self):
        """
            make sure all areas are >= 0.
            
            how abou limits on how small a SAC section can be?
        """
        self.states = self.SAC_run_area_ini()
        return
    
    def ini_locations(self):
        #self.states = self.bow_fairness_curve_location()
        self.derive_bow_fairness_section_loc_rule_fromSAC()
        self.derive_stern_fairness_section_loc_rule_fromSAC()
        self.derive_stern_fairness_section_loc_rule_fromCPKeel()
        self.derive_bow_fairness_section_loc_rule_fromDWL()
        self.derive_bow_fairness_section_loc_rule_fromCPKeel()
        self.derive_stern_fairness_section_loc_rule_fromDWL()
        #self.states = self.derive_stern_fairness_section_loc_rule_fromSAC()
        return
    
    def sac_XC_ini(self):
        self.states = self.SAC_Xc_init()
        return
        
    def block_coefficient(self):
        """
        lpp     : length between perpendiculars or wL length
        br      : some breadth (e.g. breadth at the waterline)
        draft   : draft
        vol     : volumetric displacement/water weight
        
                Vol/(Lwl*Bwl*Draft) = Cb
        """
        lwl = self.lwl
        br = self.bwl
        draft = self.draft
        vol = self.vol
        Cb = self.Cb
        s = self.states
        statesi = copy.copy(s)
        
        c1 = lp.Variable('c1')
        c2 = lp.Variable('c2')
        self.set_val(c1,None)
        self.set_val(c2,None)
        dlist = [c1,c2]

        states = (self * (lwl,br,c1) * (c1,draft,c2) / (vol,c2,Cb) ) 
        states = (states * (lwl,br,c1) * (c1,draft,c2) / (vol,c2,Cb) ) 
        states = (states * (lwl,br,c1) * (c1,draft,c2) / (vol,c2,Cb) ) 
        
        
        vars = [lwl,br,draft,vol,Cb]
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,dlist)
        return states
    
    
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
        s = self.states
        statesi = copy.copy(s)
        
        c1 = lp.Variable('c1')
        self.set_val(c1,None)
        
        states = (self * (lwl,bwl,c1) / (Awp,c1,Cwp) )
        states = (states * (lwl,bwl,c1) / (Awp,c1,Cwp) )
        
        vars = [lwl,bwl,Awp,Cwp]
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        
        states = self.clean_states(states.states,[c1])
        return states
        
        
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
        s       = self.states
        statesi = copy.copy(s)
        
        c1 = lp.Variable('c1')
        self.set_val(c1,None)
        dlist = [c1]
        
        #wrap it in another function in uKanren
        #        goal = lp.Goal.both(lp.Goal.mulo(lwl,draft,c1),
        #                            lp.Goal.divo(Acp,c1,Ccp))
                            
        states = (self * (lwl,draft,c1) / (Acp,c1,Ccp) )
        states = (states * (lwl,draft,c1) / (Acp,c1,Ccp) )
        
        vars = [lwl,draft,Acp,Ccp]
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,dlist)
        return states

        
    
        
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
        draft = self.draft
        Amsh = self.Amsh
        Cmidshp = self.Cmidshp
        s       = self.states
        statesi = copy.copy(s)
        
        c1 = lp.Variable('c1')
        self.set_val(c1,None)
        
        states = (self == (draft,dsmax))
        states = (self * (bwl,dsmax,c1) / (Amsh,c1,Cmidshp))
        states = (states * (bwl,dsmax,c1) / (Amsh,c1,Cmidshp))
                          
        vars = [dsmax,bwl,Amsh,Cmidshp,draft]
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[c1])
        return states
        
        
        
    def prismatic_coefficient(self):
        """
        Cb / Cmidshp = Cp
        """
        Cb      = self.Cb
        Cmidshp = self.Cmidshp
        Cp      = self.Cp
        
        s       = self.states
        statesi = copy.copy(s)
        
        #print 'cb = ',self(Cb)[0]
        #print 'cmidshp = ',self(Cmidshp)[0]
        #print 'cp = ',self(Cp)[0]
        
        """
        cb = ia(0.95, 0.95)
        cmidshp = ia(0.989583333336, 0.989583333337)
        cp = ia(0.96, 0.96)
        """
        
        states = (self / (Cb,Cmidshp,Cp) )
        states = (states * (Cmidshp,Cp,Cb) )
        
        vars = [Cb,Cmidshp,Cp]
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[])
        return states
        
        
    def relate_drafts(self):
        draft   = self.draft
        dsmax   = self.dsmax
        
        s       = self.states
        statesi = copy.copy(s)
        
        #states = (self <= (draft,dsmax))
        states = (self == (draft,dsmax))
        
        vars = [draft,dsmax]
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[])
        return states
    
    def LCG_constraints(self):
        lwl     = self.lwl
        LCG     = self.LCG
        
        s       = self.states
        statesi = copy.copy(s)
        states = (self <= (LCG,lwl))
        
        vars = [lwl,LCG]
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[])
        return states
        
        
    ##
    ## ------------------------------------------------- Curve Flat Interaction
    ##
    def relate_flats(self, niter=3):
        """
            Lfsac + Lfos <= LfwL
            LfwL <= lwl
        """
        lwl     = self.lwl
        lfwl    = self.lfwl
        lfsac   = self.lfsac
        lfos    = self.lfos
        lfcp    = self.lfcp
        
        s       = self.states
        statesi = copy.copy(s)
        
        c = lp.Variable('c')
        self.set_val(c,None)
        c1 = lp.Variable('c1')
        self.set_val(c1,None)
        c2 = lp.Variable('c2')
        self.set_val(c2,None)
        
        vars = [lwl,lfwl,lfsac,lfos,lfcp]

        
        states = (self <= (lfos, lwl))
        states = (states <= (lfwl, lwl)) #Feb 3 change!
        states = (states <= (lfcp, lwl)) #Feb 3 change!
        states = (states <= (lfsac, lwl)) #Feb 3 change!
        for i in range(niter):
            states = (states <= (lfsac, lwl))
            states = (states <= (lfsac, lfwl))
            states = (states <= (lfcp, lfwl))
        
        
            """# lfos  < (lwl - lfsac)"""# implemented in flat_of_side_constraint
            states = ( 
               (states - (lwl,lfsac,c1)) <= (lfos,c1) 
            )
            
            states = (states + (lfsac,lfos,c))
            states = (states <= (c,lfwl) )
            states = (states + (lfsac, lfos, c) )
            
            """# lfsac  < (lwl - lfos)"""
            states = (states - (lwl,lfos,c2) )
            states = (
               (states <= (lfsac,c2)) <= (lfos,lfwl)
             )
            states = (states <= (lfos, lwl))
            """
            ------------------------------------
            TODO: make lfwl = lfsac + lfos
                ?at the moment it only evenly distributes?
            ------------------------------------
            """

        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[c,c1,c2])
        return states
        
        
    ##
    ## ------------------------------------------------- SAC 
    ##
    
    
    
    def SAC_length_ini(self):
        lwl                 = self.lwl
        SAC_entrance_len    = self.SAC_entrance_len
        SAC_mid_len         = self.SAC_mid_len
        SAC_run_len         = self.SAC_run_len
        #loc_bfc             = self.loc_bfc
        #loc_sfc             = self.loc_sfc
        s       = self.states
        statesi = copy.copy(s)
        
        vars = [lwl,
                SAC_entrance_len,
                SAC_mid_len,
                SAC_run_len]#,
                #loc_bfc,
                #loc_sfc]
        
        #        dum = self(vol)
        #        if len(dum)>0 and isinstance(dum[0], ia):
        #            vmax = max([el.sup for el in dum])
        #            dv = ia(0.,vmax)
        #            states = (self <= (SAC_mid_area,dv))
        #        else:
        #            states = (self <= (SAC_run_area,SAC_mid_area))
        
        states = (self * (ia(.3,.7),lwl,SAC_entrance_len))
        #states = (states * (ia(.25,.4),lwl,SAC_mid_len))
        states = (states * (ia(.3,.7),lwl,SAC_run_len))
        
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[])
        return states
        
    
    
    
    def SAC_run_area_ini(self):
        vol                 = self.vol
        SAC_entrance_area   = self.SAC_entrance_area
        SAC_mid_area        = self.SAC_mid_area
        SAC_run_area        = self.SAC_run_area
        s       = self.states
        statesi = copy.copy(s)
        
        vars = [vol,
                SAC_entrance_area,
                SAC_mid_area,
                SAC_run_area]
        
        #        dum = self(vol)
        #        if len(dum)>0 and isinstance(dum[0], ia):
        #            vmax = max([el.sup for el in dum])
        #            dv = ia(0.,vmax)
        #            states = (self <= (SAC_mid_area,dv))
        #        else:
        #            states = (self <= (SAC_run_area,SAC_mid_area))
        
        states = (self <= (SAC_mid_area,vol))
        states = (states <= (SAC_run_area,SAC_mid_area))
        states = (states <= (SAC_entrance_area,SAC_mid_area) )
        
        states = (states * (ia(.3,.7),vol,SAC_entrance_area))
        states = (states * (ia(.3,.7),vol,SAC_run_area))
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[])
        return states
        
        
    def SAC_Xc_init(self):
        """
            SAC_mid_Xc <= lwl
            SAC_fwd_Xc <= SAC_mid_Xc
            SAC_mid_Xc <= SAC_run_Xc
            SAC_run_Xc <= lwl
        """
        lwl         = self.lwl
        SAC_fwd_Xc  = self.SAC_fwd_Xc
        SAC_mid_Xc  = self.SAC_mid_Xc
        SAC_run_Xc  = self.SAC_run_Xc
        s           = self.states
        statesi     = copy.copy(s)
        
        vars = [lwl,
                SAC_fwd_Xc,
                SAC_mid_Xc,
                SAC_run_Xc]
        
        dum = self(lwl)
        if len(dum)>0 and isinstance(dum[0], ia):
            lmax = max([el.sup for el in dum])
            dl = ia(0.,lmax)
            states = (self <= (SAC_mid_Xc,dl))
        else:
            states = (self <= (SAC_fwd_Xc,SAC_mid_Xc))
        states = (states <= (SAC_fwd_Xc,SAC_mid_Xc))
        states = (states <= (SAC_mid_Xc,SAC_run_Xc))
        states = (states <= (SAC_run_Xc,lwl))
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[])
        return states
    
    
    
    
    def SAC_run_area_consistency(self, niter=3):
        """TODO
        emulate sac_run_len_conistency
        """
        vol                 = self.vol
        SAC_entrance_area   = self.SAC_entrance_area
        SAC_mid_area        = self.SAC_mid_area
        SAC_run_area        = self.SAC_run_area
        Amsh                = self.Amsh
        SAC_mid_len         = self.SAC_mid_len
        lfsac           = self.lfsac
        s       = self.states
        statesi = copy.copy(s)
        
        vars = [vol,
                SAC_entrance_area,
                SAC_mid_area,
                SAC_run_area,
                Amsh,
                SAC_mid_len,
                lfsac]
                
        c1 = lp.Variable('c1')
        self.set_val(c1,None)
        c2 = lp.Variable('c2')
        self.set_val(c2,None)
        
        
        #        dummmy = self(vol)
        #        if len(dummmy)>0 and isinstance(dummmy[0], ia):
        #            dum = max([el.sup for el in dummmy])
        #            dmin = min([el.sup for el in dummmy])
        #            dum = dum*.5
        #            dmin=dmin*.2
        #            dum = ia(dmin,dum)
        #            states = (self <= (SAC_entrance_area,dum))
        #            states = (states <= (SAC_mid_area,dum))
        #            states = (states <= (SAC_run_area,dum))
        #            #-----
        #            states = (states <= (SAC_entrance_area,vol))
        #            states = (states <= (SAC_mid_area,vol))
        #            states = (states <= (SAC_run_area,vol))
        #        else:
        states = (self <= (SAC_entrance_area,vol))
        states = (states <= (SAC_mid_area,vol))
        states = (states <= (SAC_run_area,vol))
            
        states = ( (states + (SAC_run_area,SAC_mid_area,c1) ) + (SAC_entrance_area,c1,vol))
        for i in range(niter):
            states = (states + (SAC_entrance_area,SAC_mid_area,c2) )
            states = (states + (SAC_run_area,c2,vol) )
            states = (states + (SAC_run_area,SAC_mid_area,c1) )
            states = (states + (SAC_entrance_area,c1,vol))
            states = (states * (Amsh,SAC_mid_len,SAC_mid_area))    #TLM 1/24/2017 !!!!!!!!!!!!!!!!!!!!!!!
            states = (states == (SAC_mid_len,lfsac) )
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[c1,c2])
        return states
    
    
    def SAC_entrance_len_consistency(self, niter=2):
        """
        NOTE: SAC_mid_Xc is always assumed to be
        in the exact middle of the flat SAC portion!
        """
        lwl = self.lwl
        SAC_entrance_len = self.SAC_entrance_len
        SAC_mid_len = self.SAC_mid_len
        SAC_run_len = self.SAC_run_len
        s       = self.states
        statesi = copy.copy(s)
        
        vars = [lwl,
                SAC_entrance_len,
                SAC_mid_len,
                SAC_run_len]
        
        c1 = lp.Variable('c1')
        self.set_val(c1,None)
        c2 = lp.Variable('c2')
        self.set_val(c2,None)
        dlist = [c1,c2]
        #
        #********** Rules added 12-19-2016 to initilize designs
        #
        dummmy = self(lwl)
        if len(dummmy)>0 and isinstance(dummmy[0], ia):
            dum = max([el.sup for el in dummmy])
            dmin = min([el.sup for el in dummmy])
            dum = dum*.5
            dmin = dmin*.2
            dum = ia(dmin,dum)
            states = (self <= (SAC_entrance_len,dum))
            states = (states <= (SAC_mid_len,dum))
            states = (states <= (SAC_run_len,dum))
            #-----
            states = (states <= (SAC_entrance_len,lwl))
            states = (states <= (SAC_mid_len,lwl))
            states = (states <= (SAC_run_len,lwl))
        else:
            states = (self <= (SAC_entrance_len,lwl))
            states = (states <= (SAC_mid_len,lwl))
            states = (states <= (SAC_run_len,lwl))
            
        for i in range(niter):            
            #
            #*********
            #                             
            states = (states + (SAC_entrance_len,SAC_mid_len,c1) )
            states = (states + (SAC_run_len,c1,lwl) )
            states = (states + (SAC_run_len,SAC_mid_len,c2) )
            states = (states + (SAC_entrance_len,c2,lwl) )
            states = (states <= (SAC_entrance_len,lwl))
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars,
                         dlist=dlist)
        states = self.clean_states(states.states,dlist)
        return states
        
        
    def flat_of_sac_constraint(self, niter = 3):
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
        lwl             = self.lwl
        lfsac           = self.lfsac
        Amsh            = self.Amsh
        SAC_mid_area    = self.SAC_mid_area
        vol             = self.vol 
        SAC_mid_len     = self.SAC_mid_len
        
        s       = self.states
        statesi = copy.copy(s)

        vars = [lfsac,Amsh,vol,SAC_mid_len,SAC_mid_area]
        
        c1 = lp.Variable('c1')
        self.set_val(c1,None)
        #c2 = lp.Variable('c2')
        #self.set_val(c2,None)
        c3 = lp.Variable('c3')
        self.set_val(c3,None)
        
        states = (self / (vol,Amsh,c1) )
        for i in range(niter):
            states = (states <= (lfsac,c1))
            states = (states * (Amsh,lfsac,c3) )
            states = (states <= (c3,SAC_mid_area) )
            states = (states <= (lfsac,lwl) ) #the issue is that lwl is thin!!
            states = (states == (SAC_mid_len,lfsac) )
            states = (states / (vol,Amsh,c1) )
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[c1,c3])
        return states
        
        
    def SAC_mid_prop(self):
        lfsac           = self.lfsac
        SAC_mid_len     = self.SAC_mid_len
        SAC_mid_area    = self.SAC_mid_area
        s       = self.states
        statesi = copy.copy(s)
        
        
        states = (self == (SAC_mid_len,lfsac) )
        
        vars = [lfsac,SAC_mid_len,SAC_mid_area]
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[])
        return states
        
    
    #    def SAC_branch(self, section, sections=None):
    #        """mak a more generic function: called with any section
    #        to set up rules for the SAC
    #        
    #        Just a thought - not implementing anything here yet...
    #        """
    #        vol = self.vol
    #        lwl = self.lwl
    #        vars = [vol,lwl]
    #        
    #        for el in sections:
    #            vars.append(el)
    #        index = sections.index(section)
    #        this_section = sections.pop(index)
    #        ##
    #        ## Do something...
             ## ...work at the arsenal..
    #        ##
    #        s = self.states
    #        statesi = copy.copy(s)
    #        
    #        return
        
        
    def SAC_area(self, section, sections=None, niter=3):
        """quasi generic function: called with any section
        to ensure that the areas of the SAC sections
        always equal the total volume of the vessel
        """
        vol                 = self.vol
        if sections is None:
            SAC_entrance_area   = self.SAC_entrance_area
            SAC_mid_area        = self.SAC_mid_area
            SAC_run_area        = self.SAC_run_area
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
        
        s = self.states
        statesi = copy.copy(s)
        c1 = lp.Variable('c1')
        c2 = lp.Variable('c2')
        self.set_val(c1,None)
        self.set_val(c2,None)
        
        dlist = [c1,c2]
        
        states = (self + (sections[0],sections[1],c1))
        for i in range(niter):
            states = (states - (vol,c1,c2))
            states = (states == (this_section,c2))
            states = (states + (sections[0],sections[1],c1))
        
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars,
                         dlist=dlist)
        states = self.clean_states(states.states,dlist)
        return states
        
        
        
    def SAC_len(self, section, sections=None, niter=3):
        """ensure that the lengths of the SAC sections
        always equal th length of the lwl
        """
        lwl                 = self.lwl
        if sections is None:
            SAC_entrance_len    = self.SAC_entrance_len
            SAC_mid_len         = self.SAC_mid_len
            SAC_run_len         = self.SAC_run_len
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
        
        # setup dummy vars
        s = self.states
        statesi = copy.copy(s)
        c1 = lp.Variable('c1')
        c2 = lp.Variable('c2')
        self.set_val(c1,None)
        self.set_val(c2,None)
        dlist = [c1,c2]
        
        # Rules:
        states = (self + (sections[0],sections[1],c1) )
        for i in range(niter):
            states = (states - (lwl,c1,c2) )
            states = (states == (this_section,c2) )
            states = (states + (sections[0],sections[1],c1) )
        
        #cleanup
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars,
                         dlist=dlist)
        states = self.clean_states(states.states,dlist)
        return states
        
        

    def SAC_three_parts(self, niter=5):
        lwl     = self.lwl
        lfsac   = self.lfsac
        vol     = self.vol
        #Cb      = self.Cb
        SAC_entrance_len    = self.SAC_entrance_len
        SAC_mid_len         = self.SAC_mid_len
        SAC_run_len         = self.SAC_run_len
        SAC_entrance_area   = self.SAC_entrance_area
        SAC_mid_area        = self.SAC_mid_area
        SAC_run_area        = self.SAC_run_area
        
        
        vars = [lfsac,
                vol,
                lwl,
                SAC_entrance_len,
                SAC_mid_len,
                SAC_run_len,
                SAC_entrance_area,
                SAC_mid_area,
                SAC_run_area]
        
        #s = self.SAC_mid_prop() #conistency with len!!
        #
        #-----------SAC Length Conistency
        #
        self.states = self.SAC_len(section=SAC_run_len,
                         sections = [SAC_entrance_len,
                                    SAC_mid_len,
                                    SAC_run_len])
        self.states  = self.SAC_len(section=SAC_entrance_len,
                         sections = [SAC_entrance_len,
                                    SAC_mid_len,
                                    SAC_run_len])
        self.states  = self.SAC_len(section=SAC_mid_len,
                         sections = [SAC_entrance_len,
                                    SAC_mid_len,
                                    SAC_run_len])
        #
        #-----------SAC volume Consistency
        #
        self.states  = self.SAC_area(section = SAC_entrance_area,
                          sections = [SAC_entrance_area,
                                      SAC_mid_area,
                                      SAC_run_area])
        self.states  = self.SAC_area(section = SAC_mid_area,
                          sections = [SAC_entrance_area,
                                      SAC_mid_area,
                                      SAC_run_area])
        self.states  = self.SAC_area(section = SAC_run_area,
                          sections = [SAC_entrance_area,
                                      SAC_mid_area,
                                      SAC_run_area])
        
        self.states = self.SAC_entrance_len_consistency()
        #
        #-----------local rules
        # -set them after non-local rule comutation-
        #
        s = self.states
        statesi = copy.copy(s)
        #
        c1 = lp.Variable('c1')
        self.set_val(c1,None)
        
        c2 = lp.Variable('c2')
        self.set_val(c2,None)
        
        c3 = lp.Variable('c3')
        self.set_val(c3,None)
        
        c4 = lp.Variable('c4')
        self.set_val(c4,None)
        
        dlist = [c1,c2,c3,c4]
        
        #s = self.SAC_mid_prop()
        states = (self + (SAC_entrance_len,SAC_mid_len,c1) )
        
        for i in range(niter):
            states = (states + (SAC_run_len,c1,lwl) )
            
            #g = lp.Goal.lto(c2,lwl)
            #states = (states <= (c2,lwl))
            
            states = (states + (SAC_run_len, SAC_mid_len, c4) )
            states = (states + (SAC_entrance_len, c4, lwl) )
            """#g = lp.Goal.eq(SAC_mid_area,Amsh*SAC_mid_len) #documentation of thinking only
            """
            
            #g = lp.Goal.eq(SAC_mid_len,lfsac)
            #s = g(s)[0]
            
            #s.values[c1]=None
            #s.values[c2]=None
            #g = lp.Goal.both(
            #        lp.Goal.both(lp.Goal.addo(SAC_run_area,SAC_mid_area,c1),
            #                     lp.Goal.addo(SAC_entrance_area,c1,vol)),
            #        lp.Goal.both(lp.Goal.addo(SAC_entrance_area,SAC_mid_area,c3),
            #                lp.Goal.addo(c3,SAC_run_area,vol,)))
            #
            #g = lp.Goal.both( lp.Goal.addo(SAC_run_area,SAC_mid_area,c2),
            #                     lp.Goal.addo(SAC_entrance_area,c2,vol))
            #assert(len(g(s))==1)
            #s = g(s)[0]
            #
            states = (states + (SAC_run_area,SAC_mid_area,c2))
            states = (states + (SAC_entrance_area,c2,vol))
            """TODO: fix consistency issue between SAC lengths and SAC areas!!!
            """
            
            #self.states = self.SAC_entrance_len_consistency() #moved up
            #s = self.SAC_run_area_consistency() #conistency with len!!
            #s = self.SAC_mid_prop() #conistency with len!!
            
            #        g = lp.Goal.gto(SAC_entrance_area, ia(0.,0.))
            #        s = g(s)[0]
            #        g = lp.Goal.gto(SAC_run_area, ia(0.,0.))
            #        s = g(s)[0]
            
            states = (states >= (SAC_entrance_area,ia(0.,0.)) )
            states = (states >= (SAC_run_area, ia(0.,0.)) )
            
            states = (states + (SAC_entrance_area,SAC_mid_area,c3) )
            states = (states + (SAC_run_area,c3,vol) )
            
            states = (states + (SAC_entrance_len,SAC_mid_len,c1) )
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars,
                         dlist=dlist)
        states = self.clean_states(states.states,dlist)
        return states
        
        
    
    def constrain_SAC_Xc(self, niter=4):
        lwl     = self.lwl
        LCG     = self.LCG
        vol     = self.vol
    
        SAC_entrance_area   = self.SAC_entrance_area
        SAC_mid_area        = self.SAC_mid_area
        SAC_run_area        = self.SAC_run_area
        
        SAC_entrance_len    = self.SAC_entrance_len
        SAC_mid_len         = self.SAC_mid_len
        SAC_run_len         = self.SAC_run_len
        
        SAC_fwd_Xc = self.SAC_fwd_Xc
        SAC_mid_Xc = self.SAC_mid_Xc
        SAC_run_Xc = self.SAC_run_Xc
        
        
        s = self.states
        statesi = copy.copy(s)
        
        
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
        
        c1 = lp.Variable('c1')
        self.set_val(c1,None)
        c2 = lp.Variable('c2')
        self.set_val(c2,None)
        c3 = lp.Variable('c3')
        self.set_val(c3,None)
        
        sc1 = lp.Variable('sc1')
        self.set_val(sc1,None)
        sc2 = lp.Variable('sc2')
        self.set_val(sc2,None)
        sc3 = lp.Variable('sc3')
        self.set_val(sc3,None)
        sct = lp.Variable('sct')
        self.set_val(sct,None)
        
        a1 = lp.Variable('a1')
        self.set_val(a1,None)
        a2 = lp.Variable('a2')
        self.set_val(a2,None)
        
        
        b1 = lp.Variable('b1')
        self.set_val(b1,None)
        b2 = lp.Variable('b2')
        self.set_val(b2,None)
        b3 = lp.Variable('b3')
        self.set_val(b3,None)
        b4 = lp.Variable('b4')
        self.set_val(b4,None)
        b5 = lp.Variable('b5')
        self.set_val(b5,None)
        
        f1 = lp.Variable('f1')
        self.set_val(f1,None)
        
        dlist = [c1,c2,c3,
                 sc1,sc2,sc3,sct,
                 a1,a2,
                 b1,b2,b3,b4,b5,
                 f1]
        states = (self * (ia(.5,.5), SAC_mid_len, b1) )
        for i in range(niter):
            
            states = (states * (ia(.5,.5), SAC_mid_len, b1) )
            states = (states + (b1,SAC_entrance_len,b2) )
            states = (states == (SAC_mid_Xc,b2) )
            
            
            states = (states * (ia(.64,.68), SAC_entrance_len, SAC_fwd_Xc))
            
            # entrance Xc must be in the length of the entrance somewhere:
            states = (states <= (SAC_fwd_Xc,SAC_entrance_len) )        
        
            # run Xc must be in the length of the run somewhere:
            states = (states + (SAC_mid_len, SAC_entrance_len, b3) )
            states = (states >= (SAC_run_Xc,b3))
            
            
            # (lwl - SAC_mid_len - SAC_entrance_len) => b4
            states = (states - (lwl, b3, b4) )
            # (ia(.3,.37) * b4 ) => b5
            states = (states * (ia(.32,.36),b4,b5))
            states = (states + (b3,b5,SAC_run_Xc))
            
            
            
        
        
            #        # Sum_i ( A_i*Xc_i ) = A_tot * Xc_tot
            #        # SAC -area- == some real volume
            #"""
            states = (states * (vol,LCG,sct) ) #total vol * LCG = sct
            states = (states == (c2,sct)) #equiv of area barycenters
            #
            states = (states * (SAC_entrance_area,SAC_fwd_Xc,sc1))
            states = (states * (SAC_mid_area,SAC_mid_Xc,sc2) )
            states = (states * (SAC_run_area,SAC_run_Xc,sc3) )
            states = (states + (sc1,sc2,c1) )
            states = (states + (c1,sc3,c2) ) #sum of area time center of area
            #"""
            states = (states * (ia(.5,.5), SAC_mid_len, b1) )
            #states = (states * (ia(.1,.3), SAC_run_len, SAC_run_Xc))
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,dlist)
        return states
    
        
        
    
        
    ##
    ## ------------------------------------------------- Water Line
    ##
    def flat_of_wl_constraint(self, niter=2):
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
        s = self.states
        statesi = copy.copy(s)
        
        vars = [lfwl,bwl,Awp]
        
        
        c = lp.Variable('c')
        self.set_val(c,None)
        
        states = (self / (Awp,bwl,c) )
        for i in range(niter):
            states = (states <= (lfwl,c) )
            states = (states / (Awp,bwl,c) )
        
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[c])
        return states
    
    
    ##
    ## ------------------------------------------------- Keel Profile
    ##
    def flat_of_keel_constraint(self, niter=2):
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
        s = self.states
        statesi = copy.copy(s)
        vars = [lfcp,draft,Acp]
        
        c = lp.Variable('c')
        self.set_val(c,None)
        
        """ lfcp <= Acp/draft"""
        states = ((self / (Acp,draft,c)) <= (lfcp,c))
        for i in range(niter):
            states = ((states / (Acp,draft,c)) <= (lfcp,c))
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[c])
        return states
    
    
    ##
    ## ------------------------------------------------- Flat of Side
    ##
    
    
    def flat_of_side_constraint(self, niter=6):
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
        s       = self.states
        statesi = copy.copy(s)
        
        vars = [bwl,draft,lfos,vol,lfsac,Amsh,lwl]
  
  
        c1 = lp.Variable('c1')
        self.set_val(c1,None)
        c2 = lp.Variable('c2')
        self.set_val(c2,None)
        c3 = lp.Variable('c3')
        self.set_val(c3,None)
        c4 = lp.Variable('c4')
        self.set_val(c4,None)
        c5 = lp.Variable('c5')
        self.set_val(c5,None)
        
        states = ( (self * (bwl,draft,c1))  * (lfsac,Amsh,c2) )
        for i in range(niter):
            """#needed?"""
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
            states = ( ( (states - (vol,c2,c3))  / (c3,c1,c4) ) <= (lfos,c4))
    
            
            #very important relations below 
            #(should be subsummed into relate_flats but relate flats is failing...):
            #or just re-write flats with these.
            #        """# lfsac  < (lwl - lfos)"""
            #        g = lp.Goal.both(lp.Goal.subo(lwl,lfos,c5),
            #                         lp.Goal.lto(lfsac,c5) )
            #        s = g(s)[0]
            #        
            
            """# lfos  < (lwl - lfsac)"""
            #        g = lp.Goal.both(lp.Goal.subo(lwl,lfsac,c5),
            #                         lp.Goal.lto(lfos,c5) )
            states = ( (states * (bwl,draft,c1))  * (lfsac,Amsh,c2) )
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[c1,c2,c3,c4,c5])
        return states

        
    ##
    ## ------------------------------------------------- fore-n-aft fairness curves 
    ##
        
    def shape_bow_section_inizers(self):
        bbfc    = self.bbfc
        draft   = self.draft
        bwl     = self.bwl
        dbfc    = self.dbfc
        Abfc    = self.Abfc
        Cbfc    = self.Cbfc
        Amsh    = self.Amsh
        s       = self.states
        statesi = copy.copy(s)
        
        vars   = [bbfc,dbfc,draft,Abfc,Cbfc,Amsh]
  
        #c1 = lp.Variable('c1')
        #self.set_val(c1,None)
        #dlist = [c1]
        
        #bwlmax = self.states(bwl).sup
        #draftmax = self.states(draft).sup
        #states = (states== ( bbfc,ia(0.1,bwlmax) ) )
        #states = (states== ( dbfc,ia(.1,draftmax) ) )
        
        states = (self == ( Cbfc,ia(.4,.8) ) )   
        #states = (states * (Abfc,ia(.1,.5),Amsh) )    
        states = (states <= (bbfc,bwl))
        states = (states <= (dbfc,draft))
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[])
        return states

    def bow_fairness_section_area(self, niter=2):
        """
        bfc   -Bow Fairness Curve-
        bbfc    : breadth  @bfc
        Abfc    : area  @bfc
        dbfc    : Draft  @bfc
        Cbfc    : Area Coefficient @bfc
        
                Abfc/(bbfc*dbfc) = Cbfc
        """
        bbfc    = self.bbfc
        draft   = self.draft
        bwl     = self.bwl
        dbfc    = self.dbfc
        Abfc    = self.Abfc
        Cbfc    = self.Cbfc
        s = self.states
        statesi = copy.copy(s)
        
        vars   = [bbfc,draft,bwl,dbfc,Abfc,Cbfc]
        
        c1     = lp.Variable('c1')
        self.set_val(c1,None)
        dlist = [c1]
        
        
        states = (self <= (Cbfc, ia(1.,1.)))
        states = (states * (bbfc,dbfc,c1) / (Abfc,c1,Cbfc))
        for i in range(niter):
            states = (states <= (bbfc,bwl) )
            states = (states <= (dbfc,draft) )
            states = (states * (bbfc,dbfc,c1) )
            states = (states / (Abfc,c1,Cbfc))
            
        
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,dlist)
        return states

    def bow_fairness_curve_location_from_SAC(self, niter=2):
        """a somewhat linear relation between 
        bfc curve location and the ratio
        of Area of BFC to area of Midship
        
        
        Abfc    : area  @bfc
        Amsh    : area mid ship
        loc_bfc : position (station) of th bow fairness curve
        
            
            loc_bfc = (Abfc/Amsh)*SAC_entrance_len
        """
        Abfc    = self.Abfc
        loc_bfc = self.loc_bfc
        Amsh    = self.Amsh
        SAC_entrance_len = self.SAC_entrance_len
        #
        s = self.states
        statesi = copy.copy(s)
        #
        vars   = [Abfc,loc_bfc,Amsh,SAC_entrance_len]
        
        c1     = lp.Variable('c1')
        self.set_val(c1,None)
        dlist = [c1]
        
        states = (self <= (loc_bfc, SAC_entrance_len) )
        states = (states >= (Abfc, ia(0.1,0.1)) )
        states = (states / (Abfc,Amsh,c1))
        for i in range(niter):
            states = (states * (c1,SAC_entrance_len,loc_bfc))
            #states = (states == (c1,ia(0.2,.9)))
            states = (states / (Abfc,Amsh,c1))
            
        #states = self
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,dlist)
        return states

    def bow_fairness_curve_location_from_DWL(self, niter=2):
        """a somewhat linear relation between 
        bfc curve location and the ratio
        of Area of BFC to area of Midship
        
        
        Abfc    : area  @bfc
        Amsh    : area mid ship
        loc_bfc : position (station) of th bow fairness curve
        
            
            loc_bfc = (Abfc/Amsh)*SAC_entrance_len
        """
        bbfc    = self.bbfc
        loc_bfc = self.loc_bfc
        bwl    = self.bwl
        SAC_entrance_len = self.SAC_entrance_len
        #
        s = self.states
        statesi = copy.copy(s)
        #
        vars   = [bbfc,loc_bfc,bwl,SAC_entrance_len]
        
        c1     = lp.Variable('c1')
        self.set_val(c1,None)
        dlist = [c1]
        

        states = (self <= (loc_bfc, SAC_entrance_len) )
        states = (states / (bbfc,bwl,c1))
        for i in range(niter):
            states = (states * (c1,SAC_entrance_len,loc_bfc))
            #states = (states == (c1,ia(0.2,.9)))
            states = (states / (bbfc,bwl,c1))
            
        #states = self
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,dlist)
        return states

    def bow_fairness_curve_location_from_CPKeel(self, niter=2):
        """a somewhat linear relation between 
        bfc curve location and the ratio
        of Area of BFC to area of Midship
        
        
        Abfc    : area  @bfc
        Amsh    : area mid ship
        loc_bfc : position (station) of th bow fairness curve
        
            
            loc_bfc = (Abfc/Amsh)*SAC_entrance_len
        """
        dbfc    = self.dbfc
        loc_bfc = self.loc_bfc
        draft    = self.draft
        SAC_entrance_len = self.SAC_entrance_len
        #
        s = self.states
        statesi = copy.copy(s)
        #
        vars   = [dbfc,loc_bfc,draft,SAC_entrance_len]
        
        c1     = lp.Variable('c1')
        self.set_val(c1,None)
        dlist = [c1]
        

        states = (self <= (loc_bfc, SAC_entrance_len) )
        states = (states / (dbfc,draft,c1))
        for i in range(niter):
            states = (states * (c1,SAC_entrance_len,loc_bfc))
            #states = (states == (c1,ia(0.2,.9)))
            states = (states / (dbfc,draft,c1))
            
        #states = self
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,dlist)
        return states
    
    def shape_stern_section_inizers(self):
        bsfc   = self.bsfc
        draft  = self.draft
        bwl    = self.bwl
        dsfc   = self.dsfc
        Asfc   = self.Asfc
        Csfc   = self.Csfc
        Amsh   = self.Amsh
        s       = self.states
        statesi = copy.copy(s)
        
        vars   = [bsfc,dsfc,bwl,draft,Asfc,Csfc,Amsh]
        
        #c1 = lp.Variable('c1')
        #self.set_val(c1,None)
        #dlist = [c1]
        
        #bwlmax = self.state(bwl).sup
        #draftmax = self.state(draft).sup
        #states = (states == (bsfc,ia(0.1,bwlmax)))
        #states = (states == (dsfc,ia(0.1,draftmax)) )
        
        states = (self == (Csfc,ia(0.4,.8)))   
        #states = (states * (Asfc,ia(.1,.5),Amsh) ) 
        states = (states <= (bsfc,bwl))
        states = (states <= (dsfc,draft))
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[])
        return states


    def stern_fairness_section_area(self, niter=2):
        """
        TBD if this is chosen to be used.
        sfc   -Stern Fairness Curve-
        bsfc    : breadth at the stern fairness curve
        Asfc    : area of the stern fairness curve
        dsfc    : Draft at stern fairnness curve
        Cfsc    : stern fairness Coefficient Csfc
        
        Asfc/(bsfc*dsfc) = Csfc
        """
        draft  = self.draft
        bwl    = self.bwl
        bsfc   = self.bsfc
        dsfc   = self.dsfc
        Asfc   = self.Asfc
        Csfc   = self.Csfc
        
        s = self.states
        statesi = copy.copy(s)
        
        c1 = lp.Variable('c1')
        self.set_val(c1,None)
        
        
        states = (self <= (Csfc, ia(1.,1.)))
        states = (states * (bsfc,dsfc,c1) / (Asfc,c1,Csfc) )
        for i in range(niter):
            states = (states <= (bsfc, bwl))
            states = (states <= (dsfc, draft))
            states = (states * (bsfc,dsfc,c1) / (Asfc,c1,Csfc) )
        
        vars   = [bsfc,dsfc,draft,Asfc,Csfc]
        
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[c1])
        return states

    def stern_fairness_curve_location_from_SAC(self, niter=5):
        """a somewhat linear relation between 
        sfc curve location and the ratio
        of Area of SFC to area of Midship
        
        
        Asfc    : area  @bfc
        Amsh    : area mid ship
        loc_sfc : position (station) of th bow fairness curve
        
            
            loc_bfc = (1.-(Abfc/Amsh))*SAC_run_len + SAC_entrance_len + SAC_mid_len
        """
        Asfc                = self.Asfc
        loc_sfc             = self.loc_sfc
        Amsh                = self.Amsh
        SAC_run_len         = self.SAC_run_len
        SAC_mid_len         = self.SAC_mid_len
        SAC_entrance_len    = self.SAC_entrance_len
        lwl                 = self.lwl
        #
        s = self.states
        statesi = copy.copy(s)
        #
        vars   = [Asfc,loc_sfc,
                  Amsh,SAC_run_len, 
                  SAC_entrance_len,
                  SAC_mid_len,
                  lwl]
        
        c1     = lp.Variable('c1')
        self.set_val(c1,None)
        c2     = lp.Variable('c2')
        self.set_val(c2,None)
        c3     = lp.Variable('c3')
        self.set_val(c3,None)
        c4     = lp.Variable('c4')
        self.set_val(c4,None)
        dlist = [c1,c2,c3,c4]
        
        states = (self <= (loc_sfc, lwl) )
        states = (states >= (Asfc, ia(0.1,0.1)) )
        states = (states / (Asfc,Amsh,c1))
        for i in range(niter):
            states = (states - (ia(1.,1.),c1,c2)) #does not require extra work to get c2!
            states = (states * (c2,SAC_run_len,c3))
            states = (states + (SAC_entrance_len,SAC_mid_len,c4))
            states = (states + (c4,c3,loc_sfc))
            #states = (states == (c1,ia(0.2,.9)))
            states = (states / (Asfc,Amsh,c1))
            

        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,dlist)
        return states

    def stern_fairness_curve_location_from_DWL(self, niter=5):
        """a somewhat linear relation between 
        sfc curve location and the ratio
        of Area of SFC to area of Midship
        
        
        Asfc    : area  @bfc
        Amsh    : area mid ship
        loc_sfc : position (station) of th bow fairness curve
        
            
            loc_bfc = (1.-(Abfc/Amsh))*SAC_run_len + SAC_entrance_len + SAC_mid_len
            
        TODO:  DWL in three parts
            ->such that each part lines up with SAC in 3 parts?
        """
        bsfc                = self.bsfc         #was Asfc
        loc_sfc             = self.loc_sfc
        bwl                = self.bwl           #was Amsh
        SAC_run_len         = self.SAC_run_len
        SAC_mid_len         = self.SAC_mid_len
        SAC_entrance_len    = self.SAC_entrance_len
        #
        s = self.states
        statesi = copy.copy(s)
        #
        vars   = [bsfc,loc_sfc,bwl,SAC_run_len, SAC_entrance_len,SAC_mid_len]
        
        c1     = lp.Variable('c1')
        self.set_val(c1,None)
        c2     = lp.Variable('c2')
        self.set_val(c2,None)
        c3     = lp.Variable('c3')
        self.set_val(c3,None)
        c4     = lp.Variable('c4')
        self.set_val(c4,None)
        dlist = [c1,c2,c3,c4]
        
        states = (self / (bsfc,bwl,c1))
        for i in range(niter):
            states = (states - (ia(1.,1.),c1,c2)) #does not require extra work to get c2!
            states = (states * (c2,SAC_run_len,c3))
            states = (states + (SAC_entrance_len,SAC_mid_len,c4))
            states = (states + (c4,c3,loc_sfc))
            #states = (states == (c1,ia(0.2,.9)))
            states = (states / (bsfc,bwl,c1))
            

        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,dlist)
        return states

    def stern_fairness_curve_location_from_CPKeel(self, niter=5):
        """a somewhat linear relation between 
        sfc curve location and the ratio
        of Area of SFC to area of Midship
        
        
        Asfc    : area  @bfc
        Amsh    : area mid ship
        loc_sfc : position (station) of th bow fairness curve
        
            
            loc_bfc = (1.-(Abfc/Amsh))*SAC_run_len + SAC_entrance_len + SAC_mid_len
            
        TODO:  DWL in three parts
            ->such that each part lines up with SAC in 3 parts?
        """
        dsfc                = self.dsfc         #was Asfc
        loc_sfc             = self.loc_sfc
        draft               = self.draft           #was Amsh
        SAC_run_len         = self.SAC_run_len
        SAC_mid_len         = self.SAC_mid_len
        SAC_entrance_len    = self.SAC_entrance_len
        #
        s = self.states
        statesi = copy.copy(s)
        #
        vars   = [dsfc,loc_sfc,draft,SAC_run_len, SAC_entrance_len,SAC_mid_len]
        
        c1     = lp.Variable('c1')
        self.set_val(c1,None)
        c2     = lp.Variable('c2')
        self.set_val(c2,None)
        c3     = lp.Variable('c3')
        self.set_val(c3,None)
        c4     = lp.Variable('c4')
        self.set_val(c4,None)
        dlist = [c1,c2,c3,c4]
        
        states = (self / (dsfc,draft,c1))
        for i in range(niter):
            states = (states - (ia(1.,1.),c1,c2)) #does not require extra work to get c2!
            states = (states * (c2,SAC_run_len,c3))
            states = (states + (SAC_entrance_len,SAC_mid_len,c4))
            states = (states + (c4,c3,loc_sfc))
            #states = (states == (c1,ia(0.2,.9)))
            states = (states / (dsfc,draft,c1))
            

        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,dlist)
        return states
        
        
        
        
    def print_state(self):
        keys = self.keys
        for i,state in enumerate(self.states):
            print '\n state {}'.format(i)
            for key in keys:
                if isinstance(key, lp.Variable):
                    print '{}={}'.format(key,state(key))
        return
    
          
        
    def print_widths(self):
        keys = self.keys
        for i,state in enumerate(self.states):
            print '\n state {}'.format(i)
            for key in keys:
                if isinstance(key, lp.Variable):
                    print '{}={}'.format(key,state(key)[0].width())
        return
    
        
    
if __name__ == "__main__":
    self = HullDesignNode()
    #"""
    self.dsmax = ia(-5.,30.)
    #self.dsmax = ia(0.,30.)
    #self.Cb = ia(-1.,1.)
    self.Cb = ia(0.,1.)
    self.vol = ia(2000.,300000.)
    self.lwl = ia(100.,120.)
    self.bwl = ia(0.,20.) 
    #"""
    #self.Cb = ia(-1.,1.)
    #    self.vol = ia(2000.,300000.)
    #    self.lwl = ia(50.,120.)
    #    self.bwl = ia(5.,20.) 
    #self.draft = ia(5.5,20.)
    #self.AC_revise()
    #print '-----------------'
    #print self
    #
    self.Cp = ia(-1.,1.)
    #
    #self.Cp = ia(0.,1.)
    #self.AC_revise()
    #print '-----------------'
    #print self
    #self.draft = ia(-5.5,20.)
    #
    self.draft = ia(0.,20.)
    #
    #self.Awp = ia(-.1,1.)
    
    #"""
    #root = lp.States(self.states)
    root = DesignTree(self.states)
    print 'design root, lwl = ',root(self.lwl)
    self.states = root.states
    self.lwl = ia(110.,120.)
    b = lp.States(self.states, parent=root)
    root.children.append(b)
    print 'design a, lwl = ',root(self.lwl)
    print 'design b from a, lwl = ',b(self.lwl)
    #"""
    self.Amsh = ia(400.,400.)
    #self.Abfc = ia(100.,100.)
    self.print_coeff_bow_curve()
    print self(self.SAC_entrance_len)
    print self(self.loc_bfc)
    #self.loc_bfc = ia(8.25,8.25)
    print self(self.loc_bfc)
    
    #self._verbose = True
    #self.Asfc = ia(100.,100.)
    self.print_coeff_stern_curve()
    print self(self.SAC_run_len)
    print self(self.loc_sfc)
    self.AC_revise()
    #self.SAC_mid_len = ia(0.,100)
    print self(self.loc_sfc)
    #self.loc_sfc = ia(79.75,79.75)
    #self.loc_sfc = ia(80.,85.)
    #self.loc_sfc = ia(47.75,47.75)
    #self.loc_sfc = ia(57.75,57.75)
    #self.loc_sfc = ia(159.0,159.0)
    #    self.AC_revise()
    #    print self(self.loc_sfc)
    #    self.AC_revise()
    #    print self(self.loc_sfc)
    #    self.AC_revise()
    #    print self(self.loc_sfc)