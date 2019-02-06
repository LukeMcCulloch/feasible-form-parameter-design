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

print 'issue : did not normalize coefficients to ia(0.,1.)'


def ini_coeffs(rgp,Coefficients):#clist, rgp):
    """initialize coefficients to ia(0.,1.)
    """
    for co in Coefficients:
        co   = co == ia(0.,1.)
    
        rgp.add_one_rule(co,co.name)
    
    rgp.compute_fresh_rules_graph()
    return rgp

def do_first(sobol_seq=None,
             ellist=None):
    
    if sobol_seq is None:
        sobol_seq = sobol.sobolSeq([1,1],[1,1])
    if ellist is None:
        ellist = [BulbDepth,
                  CBulbMid,
                  CBulbBlock]
    ith = -1
    #x=.5
    redo = []
    #last_valid = None
    for i in range(len(ellist)):
        x = sobol_seq.next()
        var = ellist[i]
        mk_name, mk_val = rgp.get_name_and_value(var)
        print mk_name, mk_val
        val = mk_val[ith].getpoint(x)
        value = ia(val,val)
        var = var == value
        #-----------------------------------------
        ret_state = copy.copy(rgp.env)
        rgp.add_one_rule(var,var.name)
        rgp.compute_fresh_rules_graph()
        #-----------------------------------------
        
        if len(rgp.env.states)==0:
            #rgp.nodes.pop()
            #rgp.env = ret_state
            rgp.reset(ret_state,var)
            #rgp.env.states = ret_state.states
            redo.append(var)
            #last_valid = ret_state
        else:
            if rgp._verbose: print 'done ', mk_name,'=>', value
    return rgp,redo, None, sobol_seq


def test_rest(sobol_seq=None,
             ellist=None):
    
    if sobol_seq is None:
        sobol_seq = sobol.sobolSeq([1,1],[1,1])
    if ellist is None:
        ellist = [BulbBeam,
              BulbLength,
              BulbVolume,
              CBulbPrismatic
              ]
    
    ith = -1
    #x=.5
    redo = []
    #last_valid = None
    for i in range(len(ellist)):
        x = sobol_seq.next()
        #x=.5
        var = ellist[i]
        mk_name, mk_val = rgp.get_name_and_value(var)
        print mk_name, mk_val
        val = mk_val[ith].getpoint(x)
        value = ia(val,val)
        var = var == value # ia(5.1,6.9) #
        #-----------------------------------------
        ret_state = copy.copy(rgp.env)
        rgp.add_one_rule(var,var.name)
        rgp.compute_fresh_rules_graph()
        #-----------------------------------------
        
        if len(rgp.env.states)==0:
            #rgp.nodes.pop()
            #rgp.env = ret_state
            rgp.reset(ret_state,var)
            
            redo.append(var)
            #last_valid = ret_state
        else:
            if rgp._verbose: print 'done ', mk_name,'=>', value
    return rgp,redo, None, sobol_seq

def do_final_test(sobol_seq=None,
                  ellist = None):
    if ellist is None:
        ellist = [BulbDepth,
                  CBulbMid,
                  CBulbBlock,
                  BulbBeam,
                  BulbLength,
                  BulbVolume,
                  A_mid,
                  #A_lateral,
                  #A_BBwl,
                  #CBulbCtrPln,
                  #CBulbWtrPln,
                  CBulbPrismatic
                  ]
        
    if sobol_seq is None: sobol_seq = sobol.sobolSeq([1,1],[1,1])
    ith = -1
    #x=.5
    redo = []
    #last_valid = None
    for i in range(len(ellist)):
        x = sobol_seq.next()
        var = ellist[i]
        mk_name, mk_val = rgp.get_name_and_value(var)
        print mk_name, mk_val
        val = mk_val[ith].getpoint(x)
        print val
        value = ia(val,val)
        var = var == value
        #-----------------------------------------
        ret_state = copy.copy(rgp.env)
        rgp.add_one_rule(var,var.name)
        rgp.compute_fresh_rules_graph()
        #-----------------------------------------
        
        if len(rgp.env.states)==0:
            #rgp.nodes.pop()
            #rgp.env = ret_state
            rgp.reset(ret_state,var)
            redo.append(var)
            #last_valid = ret_state
        else:
            if rgp._verbose: print 'done ', mk_name,'=>', value
            
    return rgp, redo, None, sobol_seq
    


def checker():
    print '\n CBulbMid'
    print gv(A_mid) / (gv(BulbBeam)*gv(BulbDepth))
    print gv(CBulbMid)
    
    #print '\n CBulbCtrPln'
    #print gv(A_lateral)/(gv(BulbLength)*gv(BulbDepth))
    #print gv(CBulbCtrPln)
    
    #print '\n CBulbWtrPln'
    #print gv(A_BBwl)/(gv(BulbLength)*gv(BulbBeam))
    #print gv(CBulbWtrPln)
    
    print '\n CBulbBlock'
    print gv(BulbVolume)/(gv(BulbLength)*gv(BulbBeam)*gv(BulbDepth))
    print gv(CBulbBlock)
    
    print '\n CBulbPrismatic'
    print gv(BulbVolume)/(gv(BulbLength)*gv(A_mid))
    print gv(CBulbPrismatic)
    return


def gv(var):
    return rgp.get_name_and_value(var)[1][0]



def check_rule(rule,
               rgp,
               sol_map=None):
    """
    rgp.vars : mk_vars
    how can we map from a rule
    to the envrionrment's mk_values
    to check the rule?
    
    this thing will not work if you -make- a new version
    of the same rule using the top level stuff
    
    
        CBulbMid = CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
        
    if you -redo- the rule so as to toss it in the 
    checker, it will error out on the new 
    
    -vXX- numbered connective
    
    """
    if sol_map is None:
        sol_map = {}
    print 'doing rule: ',rule
    try:
        print 'name = ',rule.name
    except:
        print 'rule ',rule,'has no name'
    sol_map[rule.name] = rgp.env(rgp.varsmap[rule.name])
    if len(rule.args)>0:
        for arg in rule.args:
            if lp.isa(arg,lp.PStates):
                sol_map = check_rule(arg,rgp,sol_map)
    return sol_map



'''
    *******************************************
    *******************************************-
    Checkers
'''
def check_node(node,rgp,sol_map=None):
    if sol_map is None:
        sol_map = {}
        
    for key in node.vars:
        var = node.vars[key]
        sol_map[var.name] = rgp.env(var)
    
    return sol_map


#flist = rgp.varsmap['BulbDepth'].flist
#flist = rgp.varsmap['CBulbMid'].flist
def check_flist(flist,rgp):
    mmap = {}
    for i, node in enumerate(flist):
        #mmap[i] = {}
        #for rule in node:
        mmap = check_node(node,rgp,mmap)
    return mmap


def get_map(var):
    return rgp.varsmap[var.name].flist

def check_var(var):
    flist = get_map(var)
    mmap = check_flist(flist, rgp)
    return mmap


'''
    *******************************************
    *******************************************
    Done with Checkers
'''

if __name__ == """__main__""":
    
    start = True
    after_ship = True
    check_interval_splitting = True
    start_ship_to_bulb = True
    final_test = True
    
    if start:
    
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
        #BulbDepth = BulbDepth == ia(-10.1,10.)
        
        
        rgp = RulesGraphProcessor()
        #"""
        
        #rgp.add_one_rule(BulbDepth,'BulbDepth')
        rgp.add_one_rule(A_mid,'A_mid')
        rgp.add_one_rule(BulbBeam,'BulbBeam')
        rgp.add_one_rule(CBulbMid)
        
        #
        #*************************************
        #
        #TLM critical test:
        #BulbVolume = BulbVolume == ia(100.,500000.)
        #rgp.add_one_rule(BulbVolume) #not only should this not be needed
        #but what is happening when it isn't there??
        #
        #*************************************
        #
        
        
        #rgp.add_one_rule(BulbDepth,'BulbDepth') #should not be needed!
        
        rgp.compute_fresh_rules_graph()
        #rgp.compute_rules_graph()
        
        #rgp.AC_revise()
        #rgp.env
        
        ##
        ##************************************* end bulb
        ##
        
        ##
        ##************************************* ship
        ##
        ship_beam = lp.PStates(name='ship_beam') 
        # note: could use bare_hull var names instead. 
        #       e.g. lp.PStates(name=self.init_ship_beam.name)
        ship_depth = lp.PStates(name='ship_depth')
        ship_Lpp = lp.PStates(name='ship_Lpp')
        #
        #quantities of m**2
        ship_Amshp = lp.PStates(name='ship_Amshp')
        ship_Acp = lp.PStates(name='ship_Acp')
        ship_Awp = lp.PStates(name='ship_Awp')
        #
        #quantities of m**3
        ship_Vol = lp.PStates(name='ship_Vol')
        
        
        #-----------------------------------------
        # set the ship values in the bulb environement:
        ship_beam   = ship_beam == ia(17.4663142374, 17.4663142374)
        
        ship_depth  = ship_depth == ia(16.2051841085, 16.2051841085)
        
        ship_Lpp    = ship_Lpp == ia(111.099919763, 111.099919763)
        
        ship_Amshp  = ship_Amshp == ia(261.639572047, 261.639572047)
        
        ship_Acp    = ship_Acp == ia(1656.36308186, 1656.36308186)
        
        ship_Awp    = ship_Awp == ia(1736.75296874, 1736.75296874)
        
        ship_Vol    = ship_Vol == ia(27043.7825521, 27043.7825521)
        #-----------------------------------------
        
        
        
        #-----------------------------------------
        rgp.add_one_rule(ship_beam,'ship_beam')
        rgp.add_one_rule(ship_depth,'ship_depth')
        rgp.add_one_rule(ship_Lpp,'ship_Lpp')
        rgp.add_one_rule(ship_Amshp,'ship_Amshp')
        rgp.add_one_rule(ship_Acp,'ship_Acp')
        rgp.add_one_rule(ship_Awp,'ship_Awp')
        rgp.add_one_rule(ship_Vol,'ship_Vol')
        #-----------------------------------------
        ##
        ##************************************* end ship
        ##
        rgp.compute_fresh_rules_graph()
        #rgp.compute_rules_graph()
        
    if after_ship:
        
        ##
        ##************************************* bulb rules
        ##
        """-----------------------------------------------
        Rule:  Midbulb_Area < max_Beam * max_Depth
        CBulbMid -> [0.,1.]
        CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
        """
        CBulbMid = CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
        
        
        
        """-----------------------------------------------
        Rule: z-y_area < max_length * max_Depth
        """
        #CBulbCtrPln = CBulbCtrPln == A_lateral/(BulbLength*BulbDepth)
        
        """-----------------------------------------------
        Rule: wL_area < max_length * max_Depth
        """
        #CBulbWtrPln = CBulbWtrPln == A_BBwl/(BulbLength*BulbDepth)
        #CBulbWtrPln = CBulbWtrPln == A_BBwl/(BulbLength*BulbBeam)
        
        
        
        #2 more rules!    
        """-----------------------------------------------
        Rule: Bulb_vol < max_length * max_Depth * max *BulbBeam
        """
        CBulbBlock = CBulbBlock == BulbVolume/(BulbLength*BulbDepth*BulbBeam)
        #
        """-----------------------------------------------
        Rule: Bulb_vol < max_length * mid_bulb_area
        """
        CBulbPrismatic = CBulbPrismatic == BulbVolume/(BulbLength*A_mid)    
        #
    
    
        
        rgp.add_one_rule(CBulbMid,'CBulbMid')
        #sm = check_rule(CBulbMid,rgp)
        
        #rgp.add_one_rule(CBulbCtrPln,'CBulbCtrPln')
        #rgp.add_one_rule(CBulbWtrPln,'CBulbWtrPln')
        
        rgp.add_one_rule(CBulbBlock)
        rgp.add_one_rule(CBulbPrismatic)
    
    if check_interval_splitting:
        
        #
        # Check for interval splitting:
        #
        CBulbMid = CBulbMid == ia(-1.,1.)
        rgp.add_one_rule(CBulbMid,'CBulbMid')
        
        rgp.compute_fresh_rules_graph()
        #rgp.compute_rules_graph()
        #"""
        
        #BulbLength  = BulbLength == ia(10.,15.)
        BulbLength  = BulbLength == ia(1.,15.)
        #CBulbCtrPln = CBulbCtrPln == ia(.5,1.)
        #CBulbCtrPln = CBulbCtrPln == ia(0.,1.)
        rgp.add_one_rule(BulbLength)
        #rgp.add_one_rule(CBulbCtrPln,'CBulbCtrPln')
        
        
        #still good--------------------------------------------------
        
        rgp.compute_fresh_rules_graph()
        #rgp.compute_rules_graph()
        
        #BulbDepth = BulbDepth == ia(1.,10.)
        BulbDepth = BulbDepth == ia(5.,10.)
        rgp.add_one_rule(BulbDepth)
        
        
        
        
        #single space left, but still good---------------------------
        
        
        
        CBulbBlock  = CBulbBlock == ia(0.,1.)
        rgp.add_one_rule(CBulbBlock)
        
        CBulbPrismatic = CBulbPrismatic == ia(0.,1.)
        rgp.add_one_rule(CBulbPrismatic)
        
        
        rgp.compute_fresh_rules_graph()
        
        
        #single space left, but still good ------------------------
        
        
        ##
        ##************************************* end bulb rules
        ##
    
    if start_ship_to_bulb:
        
        ##
        ##************************************* Ship to Bulb Coefficients
        ##
        #
        #linear
        #-----------------------------------------
        Cbb    = lp.PStates(name='Cbb')
        Clpr   = lp.PStates(name='Clpr')
        Czb    = lp.PStates(name='Czb')
        #
        #nonlinear
        #-----------------------------------------
        Cabt   = lp.PStates(name='Cabt')
        Cabl   = lp.PStates(name='Cabl')
        Cvpr   = lp.PStates(name='Cvpr')
        
        #still goood-------------------------------------------------
        
        
        
        #
        ##
        ##************************************* end Ship to Bulb Coefficients
        ##
        
        
        ##
        ##************************************* nonlinear relations
        ##
        #Kracht nonlinear relations:
        Cabt    = Cabt ==  A_mid/ship_Amshp
        Cabl    = Cabl ==  A_lateral/ship_Acp #departure from Kracht
        Cvpr    = Cvpr ==  BulbVolume/ship_Vol
        #
        #"""
        rgp.add_one_rule(Cabt)
        rgp.add_one_rule(Cabl)
        rgp.add_one_rule(Cvpr)
        #"""
        #
        rgp.compute_fresh_rules_graph()
        #    
        #    
        #still goood!------------------------------------
        
    #    #
    #    ##
    #    ##************************************* end nonlinear relations
    #    ##
    #    
    #    ##
    #    ##************************************* linear relations
    #    ##
    #    #Kracht linear relations:
        Cbb     = Cbb  ==  BulbBeam/ship_beam
        Clpr    = Clpr ==  BulbLength/ship_Lpp
        Czb     = Czb  ==  BulbLength/ship_depth
        #
        #"""
        rgp.add_one_rule( Cbb)
        rgp.add_one_rule(Clpr)
        rgp.add_one_rule( Czb)
        #"""
        #
        rgp.compute_fresh_rules_graph()
    #    #
    #    ##
    #    ##************************************* end linear relations
    #    ##
    #    
    #    
    #    print '\n\n state after:'
    #    #print st
    #    
    #    print rgp.env
        
        
    
    
    if final_test:
        self = rgp
        sbs = sobol.sobolSeq([1,1],[1,1])
        redo=None
        
        
        
        ellist = [CBulbMid]
        rgp, redo, last_valid, sbs = do_first(sbs,
                                              ellist)
        rgp, redo, last_valid, sbs = do_first(sbs,
                                              redo)
        
        
        ellist = [CBulbBlock]
        rgp, redo, last_valid, sbs = do_first(sbs,
                                              ellist)
        rgp, redo, last_valid, sbs = do_first(sbs,
                                              redo)
        
        
        
        ellist = [BulbDepth]
        rgp, redo, last_valid, sbs = do_first(sbs,
                                              ellist)
        rgp, redo, last_valid, sbs = do_first(sbs,
                                              redo)
        
        
        ellist = [A_mid]
        rgp, redo, last_valid, sbs = do_first(sbs,
                                              ellist)
        rgp, redo, last_valid, sbs = do_first(sbs,
                                              redo)
        
        
        
        ellist = [BulbDepth,
                  CBulbMid,
                  CBulbBlock]
        rgp, redo, last_valid, sbs = test_rest(sbs,
                                               ellist)
        rgp, redo, last_valid, sbs = do_first(sbs,
                                              redo)
        #rgp.compute_rules_graph()
        
        
        
        ellist = [BulbBeam,
                  BulbLength,
                  BulbVolume,
                  CBulbPrismatic
                  ]
        
        rgp, redo, last_valid, sbs = do_final_test(sbs,
                                                   ellist)
        
        #rgp.compute_rules_graph()
        
        ellist = [BulbDepth,
                  CBulbMid,
                  CBulbBlock,
                  BulbBeam,
                  BulbLength,
                  BulbVolume,
                  A_mid,
                  CBulbPrismatic
                  ]
        
        
        rgp, redo, last_valid, sbs = do_final_test(sbs,
                                                   ellist)
        rgp, redo, last_valid, sbs = do_final_test(sbs,
                                                   redo)
        rgp, redo, last_valid, sbs = do_final_test(sbs,
                                                   redo)
        checker()