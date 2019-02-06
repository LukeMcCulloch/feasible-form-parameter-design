#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 18:09:31 2017

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


def ini_coeffs(clist, rgp):
    """initialize coefficients to ia(0.,1.)
    """
    for co in clist:
        co   = co == ia(0.,1.)
    
        rgp.add_one_rule(co,co.name)
    
    rgp.compute_rules_graph()
    return rgp

class BulbGenerator(object):
    """
        Rules for the Feasible Design
        of a Bulbous Bow
    """
    def __init__(self,
                 ship_beam=None,
                 ship_depth=None,
                 ship_Lpp=None,
                 ship_Amshp=None,
                 ship_Acp=None,
                 ship_Awp=None,
                 ship_Vol=None,
                 bare_hull=None):
        
        self.init_ship_beam  = ship_beam
        self.init_ship_depth = ship_depth
        self.init_ship_Lpp   = ship_Lpp
        self.init_ship_Amshp = ship_Amshp
        self.init_ship_Acp   = ship_Acp
        self.init_ship_Awp   = ship_Awp
        self.init_ship_Vol   = ship_Vol
        self.bare_hull       = bare_hull # sqKanren States() processed by opt_simple_hull

        self.rgp = RulesGraphProcessor(verbose=False)
        self.sobol_seq = sobol.sobolSeq([1,1],[1,1])
        
        
    
    def get_value(self, name):
        """
            put in this level of indirection
            because the hull is technically a list of sates
            it's not so wise to [0] everywhere!
        """
        if isinstance(name, lp.Variable):
            val = self.bare_hull(name)[0]
        return val
    
    
    
    def get_parameter(self, name):
        """Not used so far
        """
        param = self.rgp(name)
        return param
    
    
    def get_from_rgp(self, name):
        """
        -its odd that this works,
        yet dropping the treednode straight into the env
         does nothing
        yet construct returns the right answers
         supposedly without adding a new varialb to env
         
         
        cur_el = self.rgp.env.bind(self.BulbBeam.name)
        >>> self.rgp.env(cur_el)
        [ia(5.0, 8.0)]
        """
        mk_name, mk_value = self.rgp.get_name_and_value(name)
        return mk_name, mk_value
    
    
    
    def tree_search(self,
                    sobol_seq = None):
        """States tracks it's self history
        This may cause objects never to be 
        garbage collected?
        Anyway thats how design_sq was working in opt_simple_hull
        """
        if sobol_seq is None:sobol_seq=self.sobol_seq
        #if initial_state is None: initial_state = self.rgp
        
        #design = initial_state
        
        k1 = set(self.Primary)
        k2 = set(self.Coefficients)
        k3 = set(self.Areas)
        
        self.wrapped_search(k1,sobol_seq)
        self.wrapped_search(k2,sobol_seq)
        self.wrapped_search(k3,sobol_seq)
        
        
        #self.search_loop(k1,sobol_seq)
        #self.search_loop(k2,sobol_seq)
        #self.search_loop(k3,sobol_seq)
        return
    
    def wrapped_search(self,
                       klist,
                       sobol_seq,
                       maxinner = 20):
        loop = True
        inner = 0
        redolist = klist
        while loop and inner<maxinner:
            redolist = self.search_loop(redolist,
                                        sobol_seq)
            if len(redolist) == 0:
                loop=False
                print '*'
            inner +=1
        return
    
    
        
    def search_loop(self, 
                    searchlist,
                    generator,
                    maxinner = 20):
        """
            searchlist = k1
            generator = self.sobol_seq  
            
            key = k1.pop()
            key = k2.pop()
            key = k3.pop()
            
        why in the world are there two BulbLength vars ??!!
        """
        redolist = []
        for key in searchlist:
            loop = True
            inner = 0
            while loop and inner<maxinner:
                #c += 1
                x = generator.next() 
                lastnode = self.rgp.env.get_latest_child()
                checker, env = self.set_this(key,x)
                if checker:
                    nextnode = lp.States(env.states,
                                         parent = lastnode)
                    lastnode.children.append(nextnode)
                    self.rgp.env.states = nextnode.states
                    break
                else:
                    print 'now get the states one level up'
                    self.rgp.env.states = env.states
                    inner +=1
            if not checker:
                redolist.append(key)
        return set(redolist)
    
    

    
    def set_this(self, var, x, ith=-1):
        """
            var = key
            ith=-1
            
            ?Do this the old fashioned way 
            to go back one state in the event
            that an equality rule returns the null state
            
            or, -since that would require adding 
                  def equal()
                  to the code,
                -Which would then mean future users
                  building similar classes would have to
                  be smart about internals...
            
            ?how about copy.deepcopy(self.rgp.env)
            try the equality
            accept new env only if it works.
            
        """
        mk_name, mk_val = self.get_from_rgp(var) #
        val = mk_val[ith].getpoint(x)
        value = ia(val,val)
        var = var == value
        #-----------------------------------------
        ret_state = copy.deepcopy(self.rgp.env)
        self.rgp.add_one_rule(var,var.name)
        self.rgp.compute_rules_graph()
        #-----------------------------------------
        #
        #now check if ok or not
        #(checking here instead of 
          #RulesGraphProcessor's AC_revise function )
        if len(self.rgp.env.states)==0:
            #self.rgp.env = copy.deepcopy(ret_state)
            return False, ret_state
        else:
            if self.rgp._verbose: print 'done ', mk_name,'=>', value
            return True,  self.rgp.env
    
    
    
    
    def DOIT(self, rgp =None):
        if rgp is None:
            rgp = self.rgp
        #
        #linear
        #-----------------------------------------
        self.Cbb    = lp.PStates(name='Cbb')
        self.Clpr   = lp.PStates(name='Clpr')
        self.Czb    = lp.PStates(name='Czb')
        #
        #nonlinear
        #-----------------------------------------
        self.Cabt   = lp.PStates(name='Cabt')
        self.Cabl   = lp.PStates(name='Cabl')
        self.Cvpr   = lp.PStates(name='Cvpr')
        
        
    #
        #--------------------
        # quantities of m**1
        self.ship_beam = lp.PStates(name='ship_beam') 
        # note: could use bare_hull var names instead. 
        #       e.g. lp.PStates(name=self.init_ship_beam.name)
        self.ship_depth = lp.PStates(name='ship_depth')
        self.ship_Lpp = lp.PStates(name='ship_Lpp')
        #
        #quantities of m**2
        self.ship_Amshp = lp.PStates(name='ship_Amshp')
        self.ship_Acp = lp.PStates(name='ship_Acp')
        self.ship_Awp = lp.PStates(name='ship_Awp')
        #
        #quantities of m**3
        self.ship_Vol = lp.PStates(name='ship_Vol')
        #-----------------------------------------
        #
        #existing ship:
        #-----------------------------------------
        # quantities of m**1
        ship_beam   = self.ship_beam 
        ship_depth  = self.ship_depth
        ship_Lpp    = self.ship_Lpp 
        #
        #quantities of m**2
        ship_Amshp  = self.ship_Amshp 
        ship_Acp    = self.ship_Acp 
        ship_Awp    = self.ship_Awp 
        #
        #quantities of m**3
        ship_Vol    = self.ship_Vol 
        #-----------------------------------------
        #
        #-----------------------------------------
        # set the ship values in the bulb environement:
        self.ship_beam   = self.ship_beam == self.get_value(
                                    self.init_ship_beam)
        self.ship_depth  = self.ship_depth == self.get_value(
                                    self.init_ship_depth)
        self.ship_Lpp    = self.ship_Lpp == self.get_value(
                                    self.init_ship_Lpp)
        self.ship_Amshp  = self.ship_Amshp == self.get_value(
                                    self.init_ship_Amshp)
        self.ship_Acp    = self.ship_Acp == self.get_value(
                                    self.init_ship_Acp)
        self.ship_Awp    = self.ship_Awp == self.get_value(
                                    self.init_ship_Awp)
        self.ship_Vol    = self.ship_Vol == self.get_value(
                                    self.init_ship_Vol)
        #-----------------------------------------
        
        #-----------------------------------------
        rgp.add_one_rule(self.ship_beam,'ship_beam')
        rgp.add_one_rule(self.ship_depth,'ship_depth')
        rgp.add_one_rule(self.ship_Lpp,'ship_Lpp')
        rgp.add_one_rule(self.ship_Amshp,'ship_Amshp')
        rgp.add_one_rule(self.ship_Acp,'ship_Acp')
        rgp.add_one_rule(self.ship_Awp,'ship_Awp')
        rgp.add_one_rule(self.ship_Vol,'ship_Vol')
        #-----------------------------------------
        rgp.compute_rules_graph()
    
    
    
    #def initialize_bulb_parameters(self, rgp=None):
    
        #bulb areas
        self.A_mid      = lp.PStates(name='A_mid')
        self.A_lateral  = lp.PStates(name='A_lateral')
        self.A_flat     = lp.PStates(name='A_flat')
        self.A_BBwl     = lp.PStates(name='A_BBwl')
        
        
        #use the principle of the minimum square enclosing box (cartesian b/c ia is cartesian)
        self.BulbBeam    = lp.PStates(name = 'BulbBeam')     #Bulb half beam
        self.BulbDepth   = lp.PStates(name = 'BulbDepth')    #Bulb depth
        self.BulbLength  = lp.PStates(name = 'BulbLength')   #Bulb max length (min square enclosing box length)
        
        self.BulbVolume  = lp.PStates(name = 'BulbVolume') 
        
        self.CBulbMid    = lp.PStates(name = 'CBulbMid') #Bulb midsection coefficient
        self.CBulbCtrPln = lp.PStates(name = 'CBulbCtrPln') #Bulb centerplane profile area coefficient
        self.CBulbWtrPln = lp.PStates(name = 'CBulbWtrPln') #Bulb waterplane area coefficient
        
        self.CBulbBlock     = lp.PStates(name = 'CBulbBlock') #Bulb block coefficient 
        self.CBulbPrismatic = lp.PStates(name = 'CBulbPrismatic') #Bulb prismatic coefficient 
        
   
        """-----------------------------------------------
        Rule:  Midbulb_Area < max_Beam * max_Depth
        CBulbMid -> [0.,1.]
        CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
        """
        self.CBulbMid = self.CBulbMid ==  self.A_mid/(self.BulbBeam*self.BulbDepth)
        #
        """-----------------------------------------------
        Rule: z-y_area < max_length * max_Depth
        """
        self.CBulbCtrPln = self.CBulbCtrPln == self.A_lateral/(self.BulbLength*self.BulbDepth)
        #
        """-----------------------------------------------
        Rule: wL_area < max_length * max half-Beam
        """
        self.CBulbWtrPln = self.CBulbWtrPln == self.A_BBwl/(self.BulbLength*self.BulbBeam)
        #
        """-----------------------------------------------
        Rule: Bulb_vol < max_length * max_Depth * max *BulbBeam
        """
        self.CBulbBlock = self.CBulbBlock == self.BulbVolume/(self.BulbLength*self.BulbDepth*self.BulbBeam)
        #
        """-----------------------------------------------
        Rule: Bulb_vol < max_length * mid_bulb_area
        """
        self.CBulbPrismatic = self.CBulbPrismatic == self.BulbVolume/(self.BulbLength*self.A_mid)
        #
        rgp.add_one_rule(   self.CBulbMid, 'CBulbMid'   )
        rgp.add_one_rule(self.CBulbCtrPln, 'CBulbCtrPln')
        rgp.add_one_rule(self.CBulbWtrPln, 'CBulbWtrPln')
        rgp.add_one_rule(self.CBulbBlock, 'CBulbBlock')
        rgp.add_one_rule(self.CBulbPrismatic, 'CBulbPrismatic')
        #
        rgp.compute_rules_graph()
        #
    
    #def linear_parameters(self, rgp=None):
        """Kracht Design of Bulbous Bows, 1978
        """
        
        #Kracht linear relations:
        self.Cbb     = self.Cbb  ==  self.BulbBeam/self.ship_beam
        self.Clpr    = self.Clpr ==  self.BulbLength/self.ship_Lpp
        self.Czb     = self.Czb  ==  self.BulbLength/self.ship_depth
        #
        rgp.add_one_rule( self.Cbb, 'Cbb')
        rgp.add_one_rule(self.Clpr, 'Clpr')
        rgp.add_one_rule( self.Czb, 'Czb')
        #
        rgp.compute_rules_graph()
        #
    
    
    
    
    #def nonlinear_parameters(self, rgp=None):
        """Kracht Design of Bulbous Bows, 1978
        """
        
        #Kracht nonlinear relations:
        self.Cabt    = self.Cabt ==  self.A_mid/self.ship_Amshp
        self.Cabl    = self.Cabl ==  self.A_lateral/self.ship_Acp #departure from Kracht
        self.Cvpr    = self.Cvpr ==  self.BulbVolume/self.ship_Vol
        #
        self.rgp.add_one_rule(self.Cabt)
        self.rgp.add_one_rule(self.Cabl)
        self.rgp.add_one_rule(self.Cvpr)
        #
        rgp.compute_rules_graph()
        #
    
    #def experiment_bulb_parameters(self, rgp=None):
        
        self.A_mid = self.A_mid == ia(10.,15.)#20.) #issue: order of assignment "matters"
        self.BulbBeam = self.BulbBeam == ia(5.,10.)
        self.CBulbMid = self.CBulbMid == ia(.5,1.)
        self.CBulbWtrPln = self.CBulbWtrPln == ia(.5,1.)
        #BulbDepth = BulbDepth == ia(-10.1,10.)
    
        #bgp.add_one_rule(BulbDepth,'BulbDepth')
        rgp.add_one_rule(self.A_mid,'A_mid')
        rgp.add_one_rule(self.BulbBeam,'BulbBeam')
        rgp.add_one_rule(self.CBulbMid,'CBulbMid')
        rgp.add_one_rule(self.CBulbWtrPln,'CBulbWtrPln')
        #self.bgp.add_one_rule(self.BulbDepth,'BulbDepth') 
        
        self.BulbLength     = self.BulbLength == ia(10.,15.)
        self.BulbDepth      = self.BulbDepth == ia(5.,10.)
        self.CBulbCtrPln    = self.CBulbCtrPln == ia(.5,1.)
        rgp.add_one_rule(self.BulbLength,'BulbLength')
        rgp.add_one_rule(self.BulbDepth,'BulbDepth')
        rgp.add_one_rule(self.CBulbCtrPln,'CBulbCtrPln')
        
        rgp.compute_rules_graph()
        
        
        return

    
    
    def initialize_lists(self):
        
        self.Primary = [self.CBulbBlock,
                        self.BulbLength,
                        self.BulbDepth,
                        self.BulbBeam,
                        self.BulbVolume]
        
        self.Areas = [self.A_mid,
                      self.A_lateral,
                      self.A_flat,
                      self.A_BBwl]
        
        
        
        self.Coefficients = [self.Cbb,
                             self.Clpr,
                             self.Czb,
                             self.Cabt,
                             self.Cabl,
                             self.Cvpr,
                             self.CBulbMid,
                             self.CBulbCtrPln,
                             self.CBulbWtrPln,
                             self.CBulbBlock,
                             self.CBulbPrismatic]
        return 
    
    
    
    
    
    def print_bulb_state(self):
        self.print_list(self.Primary)
        self.print_list(self.Coefficients)
        self.print_list(self.Areas)
        return
        
    def print_list(self, list_):
        for key in list_:
            print self.get_from_rgp(key)
        return
    
    
    
    

def generate_a_thin_bulb_design():
    return
    
def setup_dummy_bare_hull(bbobj, rgp=None):
    if rgp is None:
        rgp = bbobj.rgp
    #-----------------------------------------
    # quantities of m**1
    bbobj.ship_beam = lp.PStates(name='ship_beam') 
    # note: could use bare_hull var names instead. 
    #       e.g. lp.PStates(name=self.init_ship_beam.name)
    bbobj.ship_depth = lp.PStates(name='ship_depth')
    bbobj.ship_Lpp = lp.PStates(name='ship_Lpp')
    #
    #quantities of m**2
    bbobj.ship_Amshp = lp.PStates(name='ship_Amshp')
    bbobj.ship_Acp = lp.PStates(name='ship_Acp')
    bbobj.ship_Awp = lp.PStates(name='ship_Awp')
    #
    #quantities of m**3
    bbobj.ship_Vol = lp.PStates(name='ship_Vol')
    #-----------------------------------------
    #
    #existing ship:
    #-----------------------------------------
    # quantities of m**1
    ship_beam   = bbobj.ship_beam 
    ship_depth  = bbobj.ship_depth
    ship_Lpp    = bbobj.ship_Lpp 
    #
    #quantities of m**2
    ship_Amshp  = bbobj.ship_Amshp 
    ship_Acp    = bbobj.ship_Acp 
    ship_Awp    = bbobj.ship_Awp 
    #
    #quantities of m**3
    ship_Vol    = bbobj.ship_Vol 
    #-----------------------------------------
    #
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
    rgp.compute_rules_graph()
    return bbobj, rgp



    
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

def old_test():
    print 'Kraft Bulbous Bow Parameters'
    print'\n Linear Parameters:'
    
    
    
    A_BT = lp.PStates(name='A_BT') #cross section area at fwd pp
    A_BL = lp.PStates(name='A_BL') #area in longitudinal plane
    Bb = lp.PStates(name='Bb')
    #c = lp.PStates(name='c')
    #d = lp.PStates(name='d')
    
    
    C_BB = lp.PStates(name='C_BB') #breadth parameter
    A_Bb = lp.PStates(name='A_Bb') #
    
    
    A_Bmid = lp.PStates(name='A_Bmid')
    A_Llateral = lp.PStates(name='A_Llateral')
    A_Lflat = lp.PStates(name='A_Lflat')
    
    A_mid = lp.PStates(name='A_mid')
    A_lateral = lp.PStates(name='A_lateral')
    A_flat = lp.PStates(name='A_flat')
    A_BBwl = lp.PStates(name='A_BBwl')
    
    
    #use the principle of the minimum square enclosing box (cartesian b/c ia is cartesian)
    BulbBeam    = lp.PStates(name = 'BulbBeam')     #Bulb half beam
    BulbDepth   = lp.PStates(name = 'BulbDepth')    #Bulb depth
    BulbLength  = lp.PStates(name = 'BulbLength')   #Bulb max length (min square enclosing box length)
    
    
    CBulbMid = lp.PStates(name = 'CBulbMid') #Bulb midsection coefficient
    CBulbCtrPln = lp.PStates(name = 'CBulbCtrPln') #Bulb centerplane profile area coefficient
    CBulbWtrPln = lp.PStates(name = 'CBulbWtrPln') #Bulb waterplane area coefficient
    
    

    
    
    #TODO: fix this with class to construct rules graph!
    
    A_mid = A_mid == ia(10.,20.) #issue: order of assignment "matters"
    BulbBeam = BulbBeam == ia(5.,10.)
    CBulbMid = CBulbMid == ia(.5,1.)
    CBulbWtrPln = CBulbWtrPln == ia(.5,1.)
    #BulbDepth = BulbDepth == ia(-10.1,10.)
    
    
    bgp = RulesGraphProcessor()
    #"""
    
    #bgp.add_one_rule(BulbDepth,'BulbDepth')
    bgp.add_one_rule(A_mid,'A_mid')
    bgp.add_one_rule(BulbBeam,'BulbBeam')
    bgp.add_one_rule(CBulbMid,'CBulbMid')
    bgp.add_one_rule(CBulbWtrPln,'CBulbWtrPln')
    bgp.compute_rules_graph()
    #bgp.add_one_rule(BulbDepth,'BulbDepth') #should not be needed!
    
    bgp.compute_rules_graph()
    #bgp.AC_revise()
    #bgp.env
    
    
    
    
    
    """-----------------------------------------------
    Rule:  Midbulb_Area < max_Beam * max_Depth
    CBulbMid -> [0.,1.]
    CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
    """
    CBulbMid = CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
    
    """-----------------------------------------------
    Rule: z-y_area < max_length * max_Depth
    """
    CBulbCtrPln = CBulbCtrPln == A_lateral/(BulbLength*BulbDepth)
    
    """-----------------------------------------------
    Rule: wL_area < max_length * max_Depth
    """
    #CBulbWtrPln = CBulbWtrPln == A_BBwl/(BulbLength*BulbDepth)
    CBulbWtrPln = CBulbWtrPln == A_BBwl/(BulbLength*BulbBeam)

    
    bgp.add_one_rule(CBulbMid,'CBulbMid')
    bgp.add_one_rule(CBulbCtrPln,'CBulbCtrPln')
    bgp.add_one_rule(CBulbWtrPln,'CBulbWtrPln')
    
    bgp.compute_rules_graph()
    #"""
    
    BulbLength  = BulbLength == ia(10.,15.)
    CBulbCtrPln = CBulbCtrPln == ia(.5,1.)
    bgp.add_one_rule(BulbLength,'BulbLength')
    bgp.add_one_rule(CBulbCtrPln,'CBulbCtrPln')
    
    bgp.compute_rules_graph()
    
    BulbDepth = BulbDepth == ia(5.,10.)
    bgp.add_one_rule(BulbDepth,'BulbDepth')
    bgp.compute_rules_graph()
    
    
    print '\n\n state after:'
    #print st
    
    print bgp.env
    return

if __name__ == "__main__":
    #old_test()
    #self = BulbGenerator()
    #self.rgp._verbose = True
    #print self.rgp.env
    
    #type(self.rgp.nodes[0].vars['ship_beam'])
    
    low_down_dirty_testing = False
    x1 = lp.PStates(name='x1')
    x2 = lp.PStates(name='x2')
    x3 = lp.PStates(name='x3')
    #x4 = lp.PStates(name='x4')
    
    x1 = x1 == ia(1.,10.)
    x2 = x2 == ia(1.,10.)
    x3 = x3 == ia(1.,10.)
    
    env = lp.States(lp.State())
    
    # add x1 rule to env
    env, vars_ = x1.construct(x1, env, {})
    fundict = x1.compile(x1,vars_)
    for el in fundict:
        env = fundict[el](env)
         
    # add x2 rule to env
    env, vars_ = x2.construct(x2, env, vars_)
    fundict = x2.compile(x2,vars_)
    for el in fundict:
        env = fundict[el](env)
        
    # add x3 rule to env
    env, vars_ = x3.construct(x3, env, vars_)
    fundict = x3.compile(x3,vars_)
    for el in fundict:
        env = fundict[el](env)
#    
#    
#    
    rgp = lp.RulesGraphProcessor()
#
    rgp.add_one_rule(x1,x1.name)
    rgp.add_one_rule(x2,x2.name)
    rgp.add_one_rule(x3,x3.name)
#    
    rgp.compute_fresh_rules_graph()
#    
    
    
    x3 = x3 == x1*x2#*x4
    env, vars_ = x3.construct(x3, env, vars_)
    fundict = x3.compile(x3,vars_)
    for el in fundict:
        env = fundict[el](env)
        
    rgp.add_one_rule(x3,x3.name)
    rgp.compute_fresh_rules_graph()
    
    
    x1 = x1 == ia(2.,2.5)
    env, vars_ = x1.construct(x1, env, vars_)
    fundict = x1.compile(x1,vars_)
    for el in fundict:
        env = fundict[el](env)
        
    
    rgp.add_one_rule(x1,x1.name)
    rgp.compute_fresh_rules_graph()
    print rgp.env
    
    # a = ia(2.,2.5)
    # b = ia(2.,5.)
    # c = ia(2.,10.)
    
    
    if low_down_dirty_testing:
        
        self = rgp
        rule = self.current_rules.pop()
        self.compute_one_rule(rule)
        #AC
        self.get_rules_set()
        self.var_sets = set([])  
        
        
        
        
        
        #
        # DO AC REVISE Step By Step!
        #
        """
        #AC step 1
        rule = self.rules_set.pop()
        fundict = rule.fundict
        vars_ = rule.vars
        save_state = copy.deepcopy(self.env)
        for key in fundict:
            self.env = self.compute_one_relation(fundict[key], 
                                                 self.env)
        print self.env.states
        print save_state.states
        states = self.env.states
        states_old = save_state.states
        state = states[0]
        statei = states_old[0]
        for el in vars_:
            print el,':'
            print state(vars_[el]) , statei(vars_[el]) ,\
            state(vars_[el]) == statei(vars_[el])
            
        
        #AC step 2
        rule = self.rules_set.pop()
        fundict = rule.fundict
        vars_ = rule.vars
        save_state = copy.deepcopy(self.env)
        for key in fundict:
            self.env = self.compute_one_relation(fundict[key], 
                                                 self.env)
        print self.env.states
        print save_state.states
        states = self.env.states
        states_old = save_state.states
        state = states[0]
        statei = states_old[0]
        for el in vars_:
            print el,':'
            print state(vars_[el]) , statei(vars_[el]) ,\
            state(vars_[el]) == statei(vars_[el])
            
        
        #AC step 2
        rule = self.rules_set.pop()
        fundict = rule.fundict
        vars_ = rule.vars
        save_state = copy.deepcopy(self.env)
        for key in fundict:
            self.env = self.compute_one_relation(fundict[key], 
                                                 self.env)
        print self.env.states
        print save_state.states
        states = self.env.states
        states_old = save_state.states
        state = states[0]
        statei = states_old[0]
        for el in vars_:
            print el,':'
            print state(vars_[el]) , statei(vars_[el]) ,\
            state(vars_[el]) == statei(vars_[el])
        #"""