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


def ini_coeffs(bbobj):#clist, rgp):
    """initialize coefficients to ia(0.,1.)
    """
    for co in bbobj.Coefficients:
        co   = co == ia(0.,1.)
    
        bbobj.rgp.add_one_rule(co,co.name)
    
    bbobj.rgp.compute_fresh_rules_graph()
    return bbobj

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
        
        if bare_hull is None:
            self = setup_dummy_bare_hull(self)
            self.rgp.compute_fresh_rules_graph()
        else:
            self.initialize_ship_parameters_and_values()
        #        self.rgp = self.initialize_bulb_parameters(self.rgp)
        #        self.rgp = self.coupling_constants(self.rgp)
        #        
        #        self.rgp = self.linear_parameters(self.rgp)
        #        self.rgp = self.nonlinear_parameters(self.rgp)
        #        
        #        
        #        self.initialize_lists()
        #        self.rgp = ini_coeffs(self.Coefficients, self.rgp)
        
        
        #self.rgp = self.experiment_bulb_parameters(self.rgp)
        #self.rgp = self.basic_bulb_rules(self.rgp)
        
        self.initialize_bulb_parameters()
        self.coupling_constants()
        
        self.linear_parameters()
        self.nonlinear_parameters()
        
        
        self.initialize_lists()
        self = ini_coeffs(self)#.Coefficients, self.rgp)
        self.experiment_bulb_parameters()
        self.basic_bulb_rules()
        
        #        self.clist = [self.Cbb,
        #                      self.Clpr,
        #                      self.Czb,
        #                      self.Cabt,
        #                      self.Cabl,
        #                      self.Cvpr,
        #                      self.CBulbMid,
        #                      self.CBulbCtrPln,
        #                      self.CBulbWtrPln,
        #                      self.CBulbBlock,
        #                      self.CBulbPrismatic]
        
        
        #BulbObject.rgp = ini_coeffs(BulbObject.clist, BulbObject.rgp)
        
    
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
        
        k1 = self.Primary
        k2 = self.Coefficients
        #k3 = self.Areas
        
        self.wrapped_search(k1,sobol_seq)
        self.wrapped_search(k2,sobol_seq)
        #self.wrapped_search(k3,sobol_seq)
        
        
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
            if len(redolist) == 0: #now this is redundant
                loop=False
                print '*'
            else:
                print 'redo list contents: ',redolist
            inner +=1
        return
    
    
        
    def search_loop_old(self, 
                    inlist,
                    generator=None,
                    maxinner = 20):
        """
            searchlist = k1
            generator = self.sobol_seq  
            
            key = k1.pop()
            key = k2.pop()
            key = k3.pop()
            
        why in the world are there two BulbLength vars ??!!
        """
        if generator is None: generator = self.sobol_seq
        #redolist = []
        searchlist = copy.copy(inlist)
        for i in range(len(searchlist)):
            print i
            key = searchlist.pop(0)
            loop = True
            inner = 0
            
            while loop:# and inner<maxinner:
                #c += 1
                x = generator.next() 
                x=.5
                #lastnode = self.rgp.env.get_latest_child()
                checker, env = self.set_this(key,x)
                if checker:
                    #nextnode = lp.States(env.states,
                    #                     parent = lastnode)
                    #lastnode.children.append(nextnode)
                    self.rgp.env.states = env.states#nextnode.states
                    break
                else:
                    print 'now get the states one level up'
                    self.rgp.env.states = env.states
                    inner +=1
                    print inner
            if not checker:
                searchlist.append(key)
        if len(searchlist)==0:
            return  searchlist
        else:
            return self.search_loop(inlist, generator)
        
    
    def search_loop(self, 
                    inlist,
                    generator=None,
                    maxinner = 20,
                    ith = -1):
        if generator is None: generator = self.sobol_seq
        redo = []
        for i in range(len(inlist)):
            var = inlist[i]
            #x = generator.next() 
            x = .5
            checker = self.set_this(var,x)
            #self.rgp.env.states = env.states
            #self.rgp.env = env
            if not checker:
                redo.append(var)
        if len(redo)==0:
            return redo
        else:
            print 'redo',redo
            return self.search_loop(redo,
                                    generator)
    
    

    
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
            
            ?how about copy.copy(self.rgp.env)
            try the equality
            accept new env only if it works.
            
        """
        mk_name, mk_val = self.get_from_rgp(var) #
        val = mk_val[ith].getpoint(x)
        value = ia(val,val)
        var = var == value
        #-----------------------------------------
        ret_state = copy.copy(self.rgp.env)
        self.rgp.add_one_rule(var,var.name)
        self.rgp.compute_fresh_rules_graph()
        #-----------------------------------------
        #
        #now check if ok or not
        #(checking here instead of 
          #RulesGraphProcessor's AC_revise function )
        if len(self.rgp.env.states)==0:
            #self.rgp.env = copy.copy(ret_state)
            #
            self.rgp.env = ret_state
            return False #, ret_state
        else:
            if self.rgp._verbose: print 'done ', mk_name,'=>', value
            return True #,  self.rgp.env
    
    
    
    
    def coupling_constants(self, rgp =None):
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
        return rgp
    
    
    def initialize_ship_parameters_and_values(self):#, rgp=None):
        #if rgp is None:
        #    rgp = self.rgp
        #-----------------------------------------
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
        self.rgp.add_one_rule(self.ship_beam,self.ship_beam.name)
        self.rgp.add_one_rule(self.ship_depth,self.ship_depth.name)
        self.rgp.add_one_rule(self.ship_Lpp,self.ship_Lpp.name)
        self.rgp.add_one_rule(self.ship_Amshp,self.ship_Amshp.name)
        self.rgp.add_one_rule(self.ship_Acp,self.ship_Acp.name)
        self.rgp.add_one_rule(self.ship_Awp,self.ship_Awp.name)
        self.rgp.add_one_rule(self.ship_Vol,self.ship_Vol.name)
        #-----------------------------------------
        self.rgp.compute_fresh_rules_graph()
        return #rgp
    
    
    
    def initialize_bulb_parameters(self):#, rgp=None):
        #if rgp is None:
        #    rgp = self.rgp
        
        #bulb areas
        self.A_mid      = lp.PStates(name='A_mid')
        self.A_lateral  = lp.PStates(name='A_lateral')
        #self.A_flat     = lp.PStates(name='A_flat') #A_BBwl instead
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
        
        return #rgp
    
    
    def initialize_lists(self):
        
        self.Primary = [#self.CBulbBlock,
                        self.BulbLength,
                        self.BulbDepth,
                        self.BulbBeam,
                        self.BulbVolume]
        
        self.Areas = [self.A_mid,
                      self.A_lateral,
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
    
    
    
    
    def basic_bulb_rules(self):#, rgp=None):
        #if rgp is None:
        #    rgp = self.rgp
        
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
        self.rgp.add_one_rule(   self.CBulbMid, self.CBulbMid.name  )
        self.rgp.add_one_rule(self.CBulbCtrPln, 'CBulbCtrPln')
        self.rgp.add_one_rule(self.CBulbWtrPln, 'CBulbWtrPln')
        self.rgp.add_one_rule(self.CBulbBlock, 'CBulbBlock')
        self.rgp.add_one_rule(self.CBulbPrismatic, 'CBulbPrismatic')
        #
        self.rgp.compute_fresh_rules_graph()
        #
        #self.rgp = rgp
        return #rgp
    
    
    def linear_parameters(self):#, rgp=None):
        """Kracht Design of Bulbous Bows, 1978
        """
        #if rgp is None:
        #    rgp = self.rgp
        
        #Kracht linear relations:
        self.Cbb     = self.Cbb  ==  self.BulbBeam/self.ship_beam
        self.Clpr    = self.Clpr ==  self.BulbLength/self.ship_Lpp
        self.Czb     = self.Czb  ==  self.BulbLength/self.ship_depth
        #
        self.rgp.add_one_rule( self.Cbb, 'Cbb')
        self.rgp.add_one_rule(self.Clpr, 'Clpr')
        self.rgp.add_one_rule( self.Czb, 'Czb')
        #
        self.rgp.compute_fresh_rules_graph()
        #
        #self.rgp = rgp
        return #rgp
    
    
    
    
    def nonlinear_parameters(self):#, rgp=None):
        """Kracht Design of Bulbous Bows, 1978
        """
        #if rgp is None:
        #    rgp = self.rgp
        
        
        #Kracht nonlinear relations:
        self.Cabt    = self.Cabt ==  self.A_mid/self.ship_Amshp
        self.Cabl    = self.Cabl ==  self.A_lateral/self.ship_Acp #departure from Kracht
        self.Cvpr    = self.Cvpr ==  self.BulbVolume/self.ship_Vol
        #
        self.rgp.add_one_rule( self.Cbb, 'Cbb')
        self.rgp.add_one_rule(self.Clpr, 'Clpr')
        self.rgp.add_one_rule( self.Czb, 'Czb')
        #
        self.rgp.compute_fresh_rules_graph()
        #
        #self.rgp = rgp
        return #rgp
    
    
    
    def experiment_bulb_parameters(self):#, rgp=None):
        #if rgp is None:
        #    rgp = self.rgp
        
        
        
        self.A_mid = self.A_mid == ia(10.,15.)#20.) #issue: order of assignment "matters"
        self.BulbBeam = self.BulbBeam == ia(5.,10.)
        self.CBulbMid = self.CBulbMid == ia(.5,1.)
        self.CBulbWtrPln = self.CBulbWtrPln == ia(.5,1.)
        #BulbDepth = BulbDepth == ia(-10.1,10.)
        #now add them
        #rgp.add_one_rule(BulbDepth,'BulbDepth')
        self.rgp.add_one_rule(self.A_mid,'A_mid')
        self.rgp.add_one_rule(self.BulbBeam,'BulbBeam')
        self.rgp.add_one_rule(self.CBulbMid,'CBulbMid')
        self.rgp.add_one_rule(self.CBulbWtrPln,'CBulbWtrPln')
        #self.rgp.add_one_rule(self.BulbDepth,'BulbDepth') 
        
        self.BulbLength     = self.BulbLength == ia(10.,15.)
        self.BulbDepth      = self.BulbDepth == ia(5.,10.)
        self.CBulbCtrPln    = self.CBulbCtrPln == ia(.5,1.)
        #now add them
        self.rgp.add_one_rule(self.BulbLength,'BulbLength')
        self.rgp.add_one_rule(self.BulbDepth,'BulbDepth')
        self.rgp.add_one_rule(self.CBulbCtrPln,'CBulbCtrPln')
        
        self.rgp.compute_fresh_rules_graph()
        
        
        return #rgp
    

def generate_a_thin_bulb_design():
    return
    
def setup_dummy_bare_hull(bbobj):#, rgp=None):
    #if rgp is None:
    #    rgp = bbobj.rgp
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
    #    ship_beam   = bbobj.ship_beam 
    #    ship_depth  = bbobj.ship_depth
    #    ship_Lpp    = bbobj.ship_Lpp 
    #    #
    #    #quantities of m**2
    #    ship_Amshp  = bbobj.ship_Amshp 
    #    ship_Acp    = bbobj.ship_Acp 
    #    ship_Awp    = bbobj.ship_Awp 
    #    #
    #    #quantities of m**3
    #    ship_Vol    = bbobj.ship_Vol 
    #-----------------------------------------
    #
    #-----------------------------------------
    # set the ship values in the bulb environement:
    bbobj.ship_beam   = bbobj.ship_beam == ia(17.4663142374, 17.4663142374)
    
    bbobj.ship_depth  = bbobj.ship_depth == ia(16.2051841085, 16.2051841085)
    
    bbobj.ship_Lpp    = bbobj.ship_Lpp == ia(111.099919763, 111.099919763)
    
    bbobj.ship_Amshp  = bbobj.ship_Amshp == ia(261.639572047, 261.639572047)
    
    bbobj.ship_Acp    = bbobj.ship_Acp == ia(1656.36308186, 1656.36308186)
    
    bbobj.ship_Awp    = bbobj.ship_Awp == ia(1736.75296874, 1736.75296874)
    
    bbobj.ship_Vol    = bbobj.ship_Vol == ia(27043.7825521, 27043.7825521)
    #-----------------------------------------
    
    #-----------------------------------------
    bbobj.rgp.add_one_rule(bbobj.ship_beam,'ship_beam')
    bbobj.rgp.add_one_rule(bbobj.ship_depth,'ship_depth')
    bbobj.rgp.add_one_rule(bbobj.ship_Lpp,'ship_Lpp')
    bbobj.rgp.add_one_rule(bbobj.ship_Amshp,'ship_Amshp')
    bbobj.rgp.add_one_rule(bbobj.ship_Acp,'ship_Acp')
    bbobj.rgp.add_one_rule(bbobj.ship_Awp,'ship_Awp')
    bbobj.rgp.add_one_rule(bbobj.ship_Vol,'ship_Vol')
    #-----------------------------------------
    #rgp.compute_fresh_rules_graph()
    return bbobj#, rgp



    
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
    
 
    
    A_mid = lp.PStates(name='A_mid')

    
    #use the principle of the minimum square enclosing box (cartesian b/c ia is cartesian)
    BulbBeam    = lp.PStates(name = 'BulbBeam')     #Bulb half beam
    BulbDepth   = lp.PStates(name = 'BulbDepth')    #Bulb depth

    
    CBulbMid = lp.PStates(name = 'CBulbMid') #Bulb midsection coefficient


    
    
    #TODO: fix this with class to construct rules graph!
    
    A_mid = A_mid == ia(10.,20.) #issue: order of assignment "matters"
    BulbBeam = BulbBeam == ia(5.,10.)
    CBulbMid = CBulbMid == ia(.5,1.)
    
    
    rgp = RulesGraphProcessor()
    
    rgp.add_one_rule(A_mid,A_mid.name)
    rgp.add_one_rule(BulbBeam,BulbBeam.name)
    rgp.add_one_rule(CBulbMid,CBulbMid.name)
    
    
    rgp.compute_fresh_rules_graph()
    
    
    """-----------------------------------------------
    Rule:  Midbulb_Area < max_Beam * max_Depth
    CBulbMid -> [0.,1.]
    CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
    """
    CBulbMid = CBulbMid ==  A_mid/(BulbBeam*BulbDepth)
    
    #rgp.add_one_rule(BulbDepth,BulbDepth.name)
    rgp.add_one_rule(CBulbMid,CBulbMid.name)
    
    rgp.compute_fresh_rules_graph()
    #"""
    
    #rgp.compute_fresh_rules_graph()
    
    #BulbDepth = BulbDepth == ia(5.,10.)
    #rgp.add_one_rule(BulbDepth,BulbDepth.name)
    #rgp.compute_fresh_rules_graph()
    
    
    
    
    
    print '\n\n state after:'
    #print st
    
    print rgp.env
    
    
    
    ellist = [BulbDepth,
              CBulbMid,
              BulbBeam,
              A_mid]
    
    ith = -1
    x=.5
    redo = []
    for i in range(len(ellist)):
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
        #rgp.compute_rules_graph()
        #-----------------------------------------
        
        if len(rgp.env.states)==0:
            rgp.env = ret_state
            redo.append(var)
        else:
            if rgp._verbose: print 'done ', mk_name,'=>', value
        
    
    return



def small_test():
    print 'Kraft Bulbous Bow Parameters'
    print'\n Linear Parameters:'
    
    
    
    
    ##
    ##************************************* bulb
    ##
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
    rgp.add_one_rule(CBulbMid,'CBulbMid')
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
    #rgp.add_one_rule(CBulbCtrPln,'CBulbCtrPln')
    #rgp.add_one_rule(CBulbWtrPln,'CBulbWtrPln')
    
    
    #
    # Check for interval splitting:
    #
    CBulbMid = CBulbMid == ia(-1.,1.)
    rgp.add_one_rule(CBulbMid,'CBulbMid')
    
    rgp.add_one_rule(CBulbBlock, 'CBulbBlock')
    rgp.add_one_rule(CBulbPrismatic, 'CBulbPrismatic')
    
    rgp.compute_fresh_rules_graph()
    #rgp.compute_rules_graph()
    #"""
    
    #BulbLength  = BulbLength == ia(10.,15.)
    BulbLength  = BulbLength == ia(1.,15.)
    #CBulbCtrPln = CBulbCtrPln == ia(.5,1.)
    #CBulbCtrPln = CBulbCtrPln == ia(0.,1.)
    rgp.add_one_rule(BulbLength,'BulbLength')
    #rgp.add_one_rule(CBulbCtrPln,'CBulbCtrPln')
    
    
    #still good--------------------------------------------------
    
    rgp.compute_fresh_rules_graph()
    #rgp.compute_rules_graph()
    
    #BulbDepth = BulbDepth == ia(1.,10.)
    BulbDepth = BulbDepth == ia(5.,10.)
    rgp.add_one_rule(BulbDepth,'BulbDepth')
    
    
    #single space left, but still good---------------------------
    
    
    
    CBulbBlock  = CBulbBlock == ia(0.,1.)
    CBulbPrismatic = CBulbPrismatic == ia(0.,1.)
    rgp.add_one_rule(BulbLength,'CBulbBlock')
    rgp.add_one_rule(CBulbPrismatic,'CBulbPrismatic')
    
    rgp.compute_fresh_rules_graph()
    #rgp.compute_rules_graph()
    
    
    #single space left, but still good ------------------------
    
    
    ##
    ##************************************* end bulb rules
    ##
    
    
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
    #    Cabt    = Cabt ==  A_mid/ship_Amshp
    #    Cabl    = Cabl ==  A_lateral/ship_Acp #departure from Kracht
    #    Cvpr    = Cvpr ==  BulbVolume/ship_Vol
    #    #
    #    #"""
    #    rgp.add_one_rule( Cbb, 'Cbb')
    #    rgp.add_one_rule(Clpr, 'Clpr')
    #    rgp.add_one_rule( Czb, 'Czb')
    #    #"""
    #    #
    #    rgp.compute_fresh_rules_graph()
    #    #rgp.compute_rules_graph()
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
#    Cbb     = Cbb  ==  BulbBeam/ship_beam
#    Clpr    = Clpr ==  BulbLength/ship_Lpp
#    Czb     = Czb  ==  BulbLength/ship_depth
#    #
#    """
#    rgp.add_one_rule( Cbb, 'Cbb')
#    rgp.add_one_rule(Clpr, 'Clpr')
#    rgp.add_one_rule( Czb, 'Czb')
#    #"""
#    #
#    rgp.compute_fresh_rules_graph()
#    #rgp.compute_rules_graph()
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
    
    
    
    
    ellist = [BulbDepth,
              CBulbMid,
              CBulbBlock,
              BulbBeam,
              BulbLength,
              BulbVolume,
              #A_mid,
              #A_lateral,
              #A_BBwl,
              #CBulbCtrPln,
              #CBulbWtrPln,
              CBulbPrismatic]
    
    sobol_seq = sobol.sobolSeq([1,1],[1,1])
    ith = -1
    #x=.5
    redo = []
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
        #rgp.compute_rules_graph()
        #-----------------------------------------
        
        if len(rgp.env.states)==0:
            rgp.env = ret_state
            #rgp.env.states = ret_state.states
            redo.append(var)
        else:
            if rgp._verbose: print 'done ', mk_name,'=>', value
        
    
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
    
    checker()
    return redo, rgp

def older_test():
    print 'Kraft Bulbous Bow Parameters'
    print'\n Linear Parameters:'
    
    
    
    
    ##
    ##************************************* bulb
    ##
    A_mid = lp.PStates(name='A_mid')
    A_lateral = lp.PStates(name='A_lateral')
    #A_flat = lp.PStates(name='A_flat')  #BBwl instead
    A_BBwl = lp.PStates(name='A_BBwl')
    
    
    #use the principle of the minimum square enclosing box (cartesian b/c ia is cartesian)
    BulbBeam    = lp.PStates(name = 'BulbBeam')     #Bulb half beam
    BulbDepth   = lp.PStates(name = 'BulbDepth')    #Bulb depth
    BulbLength  = lp.PStates(name = 'BulbLength')   #Bulb max length (min square enclosing box length)
    
    
    BulbVolume  = lp.PStates(name = 'BulbVolume') 
        
    
    CBulbMid = lp.PStates(name = 'CBulbMid') #Bulb midsection coefficient
    CBulbCtrPln = lp.PStates(name = 'CBulbCtrPln') #Bulb centerplane profile area coefficient
    CBulbWtrPln = lp.PStates(name = 'CBulbWtrPln') #Bulb waterplane area coefficient
    
    CBulbBlock     = lp.PStates(name = 'CBulbBlock')
    CBulbPrismatic = lp.PStates(name = 'CBulbPrismatic') 
    
    
    
    
    #TODO: fix this with class to construct rules graph!
    
    A_mid = A_mid == ia(10.,20.) #issue: order of assignment "matters"
    BulbBeam = BulbBeam == ia(5.,10.)
    CBulbMid = CBulbMid == ia(.5,1.)
    CBulbWtrPln = CBulbWtrPln == ia(.5,1.)
    #BulbDepth = BulbDepth == ia(-10.1,10.)
    
    
    rgp = RulesGraphProcessor()
    #"""
    
    #rgp.add_one_rule(BulbDepth,'BulbDepth')
    rgp.add_one_rule(A_mid,'A_mid')
    rgp.add_one_rule(BulbBeam,'BulbBeam')
    rgp.add_one_rule(CBulbMid,'CBulbMid')
    rgp.add_one_rule(CBulbWtrPln,'CBulbWtrPln')
    #rgp.add_one_rule(BulbDepth,'BulbDepth') #should not older_testbe needed!
    
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
    CBulbCtrPln = CBulbCtrPln == A_lateral/(BulbLength*BulbDepth)
    
    """-----------------------------------------------
    Rule: wL_area < max_length * max_Depth
    """
    #CBulbWtrPln = CBulbWtrPln == A_BBwl/(BulbLength*BulbDepth)
    CBulbWtrPln = CBulbWtrPln == A_BBwl/(BulbLength*BulbBeam)
    
    
    
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
    rgp.add_one_rule(CBulbCtrPln,'CBulbCtrPln')
    rgp.add_one_rule(CBulbWtrPln,'CBulbWtrPln')
    
    
    rgp.add_one_rule(CBulbBlock, 'CBulbBlock')
    rgp.add_one_rule(CBulbPrismatic, 'CBulbPrismatic')
    
    rgp.compute_fresh_rules_graph()
    #rgp.compute_rules_graph()
    #"""
    
    #BulbLength  = BulbLength == ia(10.,15.)
    BulbLength  = BulbLength == ia(1.,15.)
    #CBulbCtrPln = CBulbCtrPln == ia(.5,1.)
    CBulbCtrPln = CBulbCtrPln == ia(0.,1.)
    rgp.add_one_rule(BulbLength,'BulbLength')
    rgp.add_one_rule(CBulbCtrPln,'CBulbCtrPln')
    
    rgp.compute_fresh_rules_graph()
    #rgp.compute_rules_graph()
    
    #BulbDepth = BulbDepth == ia(1.,10.)
    BulbDepth = BulbDepth == ia(5.,10.)
    rgp.add_one_rule(BulbDepth,'BulbDepth')
    
    
    
    CBulbBlock  = CBulbBlock == ia(0.,1.)
    CBulbPrismatic = CBulbPrismatic == ia(0.,1.)
    rgp.add_one_rule(BulbLength,'CBulbBlock')
    rgp.add_one_rule(CBulbPrismatic,'CBulbPrismatic')
    
    rgp.compute_fresh_rules_graph()
    #rgp.compute_rules_graph()
    
    
    ##
    ##************************************* end bulb rules
    ##
    
    
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
    """
    rgp.add_one_rule( Cbb, 'Cbb')
    rgp.add_one_rule(Clpr, 'Clpr')
    rgp.add_one_rule( Czb, 'Czb')
    #"""
    #
    rgp.compute_fresh_rules_graph()
    #rgp.compute_rules_graph()
    #
    ##
    ##************************************* end nonlinear relations
    ##
    
    ##
    ##************************************* linear relations
    ##
    #Kracht linear relations:
    Cbb     = Cbb  ==  BulbBeam/ship_beam
    Clpr    = Clpr ==  BulbLength/ship_Lpp
    Czb     = Czb  ==  BulbLength/ship_depth
    #
    """
    rgp.add_one_rule( Cbb, 'Cbb')
    rgp.add_one_rule(Clpr, 'Clpr')
    rgp.add_one_rule( Czb, 'Czb')
    #"""
    #
    rgp.compute_fresh_rules_graph()
    #rgp.compute_rules_graph()
    #
    ##
    ##************************************* end linear relations
    ##
    
    
    print '\n\n state after:'
    #print st
    
    print rgp.env
    
    
    
    
    ellist = [BulbDepth,
              CBulbMid,
              CBulbBlock,
              BulbBeam,
              BulbLength,
              BulbVolume,
              #A_mid,
              #A_lateral,
              #A_BBwl,
              CBulbCtrPln,
              CBulbWtrPln,
              CBulbPrismatic]
    
    sobol_seq = sobol.sobolSeq([1,1],[1,1])
    ith = -1
    #x=.5
    redo = []
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
        #rgp.compute_rules_graph()
        #-----------------------------------------
        
        if len(rgp.env.states)==0:
            rgp.env = ret_state
            #rgp.env.states = ret_state.states
            redo.append(var)
        else:
            if rgp._verbose: print 'done ', mk_name,'=>', value
        
    
    def checker():
        print '\n CBulbMid'
        print gv(A_mid) / (gv(BulbBeam)*gv(BulbDepth))
        print gv(CBulbMid)
        
        print '\n CBulbCtrPln'
        print gv(A_lateral)/(gv(BulbLength)*gv(BulbDepth))
        print gv(CBulbCtrPln)
        
        print '\n CBulbWtrPln'
        print gv(A_BBwl)/(gv(BulbLength)*gv(BulbBeam))
        print gv(CBulbWtrPln)
        
        print '\n CBulbBlock'
        print gv(BulbVolume)/(gv(BulbLength)*gv(BulbBeam)*gv(BulbDepth))
        print gv(CBulbBlock)
        
        print '\n CBulbPrismatic'
        print gv(BulbVolume)/(gv(BulbLength)*gv(A_mid))
        print gv(CBulbPrismatic)
        return
    
    
    def gv(var):
        return rgp.get_name_and_value(var)[1][0]
    
    checker()
    return redo, rgp



if __name__ == "__main__":
    pass
    redo, rgp = older_test()
    #redo1, rgp1= small_test()
    
    """
    self = BulbGenerator()
    self.rgp._verbose = True
    #"""
    #print self.rgp.env
    
    #type(self.rgp.nodes[0].vars['ship_beam'])
    
    
    """
    
    self.CBulbMid = self.CBulbMid ==  self.A_mid/(
            self.BulbBeam*self.BulbDepth)
    self.rgp.add_one_rule(self.CBulbMid,
                          self.CBulbMid.name)
    self.rgp.compute_fresh_rules_graph()
    #"""