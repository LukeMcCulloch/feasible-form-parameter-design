#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 23:07:05 2017

@author: luke
"""

import copy
import sqKanren as lp
import sobol
ia = lp.ia
#
#from extended_interval_arithmetic import ia 
#use import from lp instead to match isinstance

GraphScaleNode = lp.GraphScaleNode
RulesGraphProcessor = lp.RulesGraphProcessor



if __name__ == """__main__""":
    A_mid = lp.PStates(name='A_mid')
    A_lateral = lp.PStates(name='A_lateral')
  
    
    #use the principle of the minimum square enclosing box (cartesian b/c ia is cartesian)
    BulbBeam    = lp.PStates(name = 'BulbBeam')     #Bulb half beam
    BulbDepth   = lp.PStates(name = 'BulbDepth')    #Bulb depth
    BulbLength  = lp.PStates(name = 'BulbLength')   #Bulb max length (min square enclosing box length)
    
    
    BulbVolume  = lp.PStates(name = 'BulbVolume') 
        
    
    CBulbMid = lp.PStates(name = 'CBulbMid') #Bulb midsection coefficient

    CBulbBlock     = lp.PStates(name = 'CBulbBlock')
    CBulbPrismatic = lp.PStates(name = 'CBulbPrismatic') 
    

    #TODO: fix this with class to construct rules graph!
    
    A_mid = A_mid == ia(10.,20.) #issue: order of assignment "matters"
    BulbBeam = BulbBeam == ia(5.,10.)
    BulbDepth = BulbDepth == ia(-10.1,10.)
    
    
    rgp = RulesGraphProcessor()
    #"""
    
    #rgp.add_one_rule(BulbDepth,'BulbDepth')
    rgp.add_one_rule(A_mid,'A_mid')
    rgp.add_one_rule(BulbBeam,'BulbBeam')
    rgp.compute_fresh_rules_graph()
    
    
    ##
    ## TEST add one rule:
    ##
    rgp.add_one_rule(CBulbMid,'CBulbMid')
    """
    self = rgp
    rule_node = CBulbMid
    rulename = rule_node.name
    
    st = self.env
        
    st, vars_ = rule_node.construct(rule_node, st, {})
    dflist = rule_node.compile(rule_node,vars_)
    
    self.env = st
    #"""
    
    CBulbMid = CBulbMid == A_mid/(BulbBeam*BulbDepth)
    
    
    
    ##
    ## TEST compute_fresh_rules_graph:
    ##
    self = rgp
    
    #rule = self.current_rules.pop()
    #self.compute_one_rule(rule) 
    
    #AC 1
    print 'starting flists AC_revise'
    self.get_rules_set()
    print 'GOT rules set = ',self.rules_set
    
    print 'reseting var_sets'
    self.var_sets = set([])
    
    """
    rule = self.rules_set.pop()
    fundict = rule.fundict #actually a dict
    vars_ = rule.vars
    save_state = copy.deepcopy(self.env)
    for key in fundict:
        self.env = self.compute_one_relation(fundict[key], 
                                             self.env)
    self.set_updates(self.env.states, save_state.states, vars_)
    
    #problem:  self.var_sets is not empty!
    state_new, state_old = self.setup_diff_debug(self.env.states, 
                                          save_state.states)
    
    el = rule.vars.keys()[0]
    q_new = state_new.get_by_str(el)
    q_old = state_old.get_by_str(el)
    #"""