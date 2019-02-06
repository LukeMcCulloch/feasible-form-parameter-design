#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 18:09:31 2017

@author: luke
"""

import copy
from extended_interval_arithmetic import ia
import sqKanren as lp

GraphScaleNode = lp.GraphScaleNode
RulesGraphProcessor = lp.RulesGraphProcessor

#==============================================================================
# class GraphScaleNode(object):
#     """
#         Takes a compiled rule, fundict
#         and its list of variables, varlist
#     """
#     def __init__(self, fundict, varlist, name=None):
#         self.name = name
#         self.fundict = fundict
#         self.vars = varlist
#         
# 
# class RulesGraphProcessor(object):
#     """
#        Let this class handle the state 
#        to ensure consistent assignment of vars 
#        and rules using them
#     """
#     def __init__(self, nodes=[], env=None,
#                  verbose=True):#funclist, varlist):
#         self.nodes = nodes
#         self.arc_sets = set([]) #active arcs
#         self.var_sets = set([]) #active variables
#         #for func in funclist:
#         #    self.add(func, varlist)
#         if env is None: env = lp.States(lp.State())
#         self.env = env
#         self._verbose = verbose
#         
# #    def add_one_rule(self, 
# #                     funclist, 
# #                     varlist, 
# #                     rulename=None,
# #                     state=None):
#     def add_one_rule(self, 
#                      rule_node,
#                      rulename=None):
#         
#         st = self.env
#         
#         st, vars_ = rule_node.construct(rule_node, st, {})
#         dflist = rule_node.compile(rule_node,vars_)
#         
#         self.env = st
#         
#         self.nodes.append(GraphScaleNode(dflist, 
#                                          vars_, 
#                                          rulename))
#         pass
#     
#     def compute_rules_graph(self):
#         for rule in self.nodes:
#             self.compute_one_rule(rule)
#             self.AC_revise()
#         return
#     
#     def compute_one_rule(self, rule, state=None):
#         """
#             DONE TODO: Issue with dummy vars_
#             Need to scrub nonexistant var 
#             or else it will trigger 
#             in
#             set_updates  !!
#             --FIXED below--
#         """
#         if state is None:
#             state = self.env
#         fundict = rule.fundict #actually a dict
#         vars_ = rule.vars
#         save_state = copy.deepcopy(state)
#         new_state = self.env
#         for el in fundict:
#             new_state = self.compute_one_relation(fundict[el], 
#                                                   state)
#         self.env = new_state
#         if 'unknown_var' in vars_: del(vars_['unknown_var'])
#         self.set_updates(new_state.states, save_state.states, vars_)
#         return
#         
#     def compute_one_relation(self, rule, state):
#         new_state = rule(state)
#         return new_state
#     
#         
#     def set_updates(self, states, states_old, vars_, tol=.001):
#         """TODO: make sure you compare the same state from
#         old to new.  That way you really do track real
#         changes.
#         
#         -Put a unique identifier on each state 
#         
#         inputs:  States.states, states_old.states
#         (dont be fooled!)
#         """
#         for state, statei in zip(states, states_old):
#             chg_lst = [vars_[el] for el in vars_ if (not state(vars_[el]) == statei(vars_[el]) )]
#             #chg_lst = [var for var in vars if not self.diff( ( state(var), statei(var),tol))]
#             var_set = self.var_sets.union(set(chg_lst)) 
#             self.var_sets = var_set
#         return 
#         
#     
#     
#     
#     def get_arc_sets(self):
#         """Find the constraints that have inputs which changed
#         """
#         #lac = len(self.arc_sets)
#         #lvs = len(self.var_sets)
#         #for c, var_set in enumerate(self.var_sets):
#             #if c>=lac:
#             #    self.arc_sets.append(set())
#         for var in self.var_sets:
#                 #self.arc_sets[c] = self.arc_sets[c].union(set(var.arcs))
#             self.arc_sets = self.arc_sets.union(set(var.arcs))
# #            else:
# #                for var in var_set:
# #                    self.arc_sets[c] = self.arc_sets[c].union(set(var.arcs))
#         return
#     
#     
#     
#     def AC_revise(self, print_=False, maxit = 8):
#         """Arc Consistency for the Hull Parameters
#         """
#         self.get_arc_sets()
#         #if print_:
#         #    print 'arcs of interest:\n',self.arc_sets
#         #self.var_sets = [set([])]
#         self.var_sets = set([])
#         
#         #for c, arc_set in enumerate(self.arc_sets):
#         count=0
#         while len(self.arc_sets)>0 and count<maxit:
#             for arc in self.arc_sets:
#                 if print_ or self._verbose: print '\ndoing ', arc
#                 self.env = arc(self.env)
#             self.arc_sets = set([])
#             self.get_arc_sets()
#             self.var_sets = set([])
#             #            self.arc_sets[c] = set([])
#             #            self.get_arc_sets()
#             #            self.var_sets[c] = set([])
#             count += 1
#             print count
#             
#         return
# 
#==============================================================================


if __name__ == "__main__":
    print '\n\n\n Do Not Use\n\n\n'
    print 'BulbRules2 '
    print 'testing '
    
    A_BT = lp.PStates(name='A_BT')
    A_BL = lp.PStates(name='A_BL')
    Bb = lp.PStates(name='Bb')
    c = lp.PStates(name='c')
    d = lp.PStates(name='d')
    
    
    C_BB = lp.PStates(name='C_BB')
    A_Bb = lp.PStates(name='A_Bb')
    
    
    A_Bmid = lp.PStates(name='A_Bmid')
    A_Llateral = lp.PStates(name='A_Llateral')
    A_Lflat = lp.PStates(name='A_Lflat')
    
    A_mid = lp.PStates(name='A_mid')
    A_lateral = lp.PStates(name='A_lateral')
    A_flat = lp.PStates(name='A_flat')
    
    BulbBeam = lp.PStates(name = 'BulbBeam')
    BulbDepth = lp.PStates(name = 'BulbDepth')
    
    
    CBulbMid = lp.PStates(name = 'CBulbMid')
    
    

    
    
    #TODO: fix this with class to construct rules graph!
    
    A_mid = A_mid == ia(10.,20.) #issue: order of assignment "matters"
    BulbBeam = BulbBeam == ia(5.,10.)
    CBulbMid = CBulbMid == ia(.5,1.)
    
    print 'problem:'
    BulbDepth = BulbDepth == ia(0.1,40.) #makes it blow down
    #BulbDepth = BulbDepth == ia(0.05,0.2)  
    print 'BulbDepth = BulbDepth == ia(0.,40.) #makes it blow down'
    ##
    ##----
    ##
    #"""
    # add CBuldMid to st
    #"""
    #initialize, if not already done
    st, vars_ = CBulbMid.construct(CBulbMid, 
                                   lp.States(lp.State({})), 
                                   {})
    fundict = CBulbMid.compile(CBulbMid,vars_)
    for el in fundict:
        st = fundict[el](st)
    #"""
    
    #add A_mid to st
    #initialize, if not already done
    #st, vars_ = A_mid.construct(A_mid, 
    #                             lp.States(lp.State({})), 
    #                               {})
    st, vars_ = A_mid.construct(A_mid, st, {})
    fundict = A_mid.compile(A_mid,vars_)
    for el in fundict:
        st = fundict[el](st)
        
    
    # add BulbBeam to st
    BulbBeam = BulbBeam == ia(5.,10.)
    #BulbBeam.verbose = True
    st, vars_ = BulbBeam.construct(BulbBeam, st, {})
    
     
    #"""
    fundict = BulbBeam.compile(BulbBeam,vars_)
    for el in fundict:
        st = fundict[el](st)
    #"""
     
    """
    # add BulbDepth to st  (not necessary as this gets picked up in the graph)
    st, vars_ = BulbDepth.construct(BulbDepth, st, {})
    fundict = BulbDepth.compile(BulbDepth,vars_)
    for el in fundict:
        st = fundict[el](st)
    #"""
    
    ##
    ##-----
    ##
    """
    CBulbMid = lp.PStates(name = 'CBulbMid')
    #CBulbMid = A_mid/(BulbBeam*BulbDepth)
    CBulbMid = CBulbMid == A_mid/(BulbBeam*BulbDepth)
    CBulbMid.name = 'CBulbMid'
    print CBulbMid
    #"""
    ##-----
    ##
    #"""
    #CBulbMid = lp.PStates(name = 'CBulbMid')
    """"""
    CBulbMid = A_mid/(BulbBeam*BulbDepth)
    """"""
    #CBulbMid = CBulbMid == BulbBeam*BulbDepth
    """"""
    #CBulbMid = CBulbMid == BulbBeam*BulbDepth
    """"""
    #CBulbMid = CBulbMid == BulbDepth*BulbBeam
    
    #CBulbMid.name = 'CBulbMid'
    #print CBulbMid
    #"""
    
    print '\n\n state before:'
    print st
    
    
    
    CBulbMid.verbose = True
    #"""
    st, vars_ = CBulbMid.construct(CBulbMid, st, {})
    fundict = CBulbMid.compile(CBulbMid,vars_)
    #"""
    #"""
    for el in fundict:
        st = fundict[el](st)
        st = fundict[el](st)
        st = fundict[el](st)
    #"""
    
    
    
    print 'this is it:'
    """
    ia(10.0,20.0) /( ia(0.5,1.0) * ia(5.0,10.0)  )
    """
    
    #bgp = RulesGraphProcessor()
    """
    bgp.add_one_rule(A_mid,'A_mid')
    bgp.add_one_rule(BulbBeam,'BulbBeam')
    bgp.add_one_rule(CBulbMid,'CBulbMid')
    bgp.compute_rules_graph()
    bgp.add_one_rule(BulbDepth,'BulbDepth') #should not be needed!
    
    bgp.compute_rules_graph()
    bgp.AC_revise()
    
    #"""
    #CBulbMid = lp.PStates(name = 'CBulbMid')
    #CBulbMid = A_mid/(BulbBeam*BulbDepth)
    #CBulbMid.name = 'CBulbMid'
    #print CBulbMid
    #bgp.add_one_rule(CBulbMid,'CBulbMid') #should not be needed!
    
    #bgp.compute_rules_graph()
    #"""
    #CBulbMid = A_mid/(BulbBeam*BulbDepth)
    #bgp.add_one_rule(CBulbMid,'CBulbMid')
    
    #bgp.compute_rules_graph()
    #"""
    
    print '\n\n state after:'
    print st
    
    
    print 'problem!!!:'
    print 'here:  BulbDepth = BulbDepth == ia(0.,40.) #makes it all blow down?'
    print 'When it is input above'
    ##
        
    """
    TODO:
        1.)  instantiation of ia value 
                with assignment -> very simple computational graph
        2.)  use returned vars_ 
                to build the association map
                for efficient rule propagation
    """
    #"""
#    #d = d.equal(ia(3.,4.))
#    d = d == ia(3.,4.)
#    d.name = 'd'
#    st, vars_d1_ = d.construct(d, st, {})
#    funclist_d1_ = d.compile(d,vars_d1_)
#    
#    treenode = d
#    vars_ = vars_d1_
#    funclist = vars_d1_
    
    #"""
    #st = funclist_d1_[0](st)
    #st = funclist_d1_[1](st)
    
    #for el in funclist_d1_:
    #    st = funclist_d1_[el](st)
    """
    d = d == ia(1.,2.)
    d.name = 'd'
    st, vars_d1_ = d.construct(d, st, {})
    funclist_d1_ = d.compile(d,vars_d1_)
    
    for el in funclist_d1_:
        st = funclist_d1_[el](st)
    #"""
    
    
    
    """ #commented this just to focus on bulb rules above
    c = lp.PStates(name='c')
    d = lp.PStates(name='d')
    d = d == ia(3.,4.)
    d.name = 'd'
    #st, vars_d1_ = d.construct(d, st, {})
    #funclist_d1_ = d.compile(d,vars_d1_)
    
    #treenode = d
    #vars_ = vars_d1_
    #funclist = vars_d1_
    
    
    
    rga = RulesGraphProcessor()
    rga.add_one_rule(d, 
                     'set d')
    
    #rule_node = d
    #self = rga
    #rule = self.nodes[0]
    #state = self.env
    #st = self.env
    
    
    rga.compute_rules_graph()
    #"""
    
    
    print '\n\n\n Do Not Use This File'
    print 'this file is of historical value only\n\n\n'