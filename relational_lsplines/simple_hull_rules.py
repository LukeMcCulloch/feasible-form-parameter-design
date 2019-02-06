#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 13:38:58 2017

@author: luke

contains: Simplified Rules System for a bare hull

    Simplifications:
        -This level will encompass only SAC, DWL, CPKeel
        - midship curves follow directly and are not used to 
        propogate back into the hull curves
        - BFC and SFC follow directly and are not used
        to propogate design information back into hull curves
        
        There is no specification of Fairness Curves
        These are derived from the hull curves
        which are specified here.
        
    Why this is good:
        e.g.
        Before I had Xc and loc_bfc = F(Abfc,Amsh) staking out
            at what station BFC should live.
        This is a double specification, but is too vauge to make into a rule
        I would say it gave rise to badly conditioned systems of constarints
        in practice.
        
    Other simplifications:
        Lfos -> Lfdwl - lfsac
        Lfokeel -> lfcpk - lfsac
        
    
    Heuristics:
        - SAC_entrace_area ~ SAC_run_area
        - LCG of SAC portions allowed to be such that
        square fwd/aft shapes would make sense
    Heuristic Todo:
        -make LCG CPK, DWL, SAC kind of follow each other
        
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

"""
NOTES:
    -breadths and beams are full beams.  
    -Areas and area ratios use breadth as such
    -Currently assuming draft is deepest at the mid section
    with no dsmax

TODO: 
    -Done: accessor on states which returns a map from states to valsues
    for that particular parameter.  Actually this lives in states.
    
additions:
    Dec 26, 2017: added Clcg, a coefficient on LCG
    relating it to lwl in hopes of eliminating some 
    of the more troublesome SAC shapes that get solved for...
    
"""
#
#class rule(object):
#    """
#    thinking variadically
#        rule = (op, var, var, var, var ...)
#        http://norvig.com/lispy.html
#    """
#    def __init__(self, rule = None):
#        self.rule = rule
#        self.rule_dict = {}
#        
#
#    def process_rule(self, rule):
#        self.rule_dict = {}
#        self.head = rule[0]
#        self.tail = rule[1:]
#        oper = rule[0]
#        variables = rule[1:]
#        
#        return
#    
    

#def list_setter(design,x,designlist):
#    """
#        if using coefficients list
#        be sure and initialize to [0.,1.]
#        to set sesible vars here!
#    """
#    #print designlist
#    for el in designlist:
#        mvar = design(el)[0]
#        #print el, mvar
#        if isinstance(mvar, ia):
#            val = mvar.getpoint(x)
#            #print 'int', val
#            val = ia(val,val)
#            design.__setattr__(el.name,val)
#        else:
#            print 'mavar is a variale'
#    return
#        
#def set_lots(self):
#    """Fooling around
#    with setting methods
#    for lists of parameters
#    """
#    x=1.
#    list_setter(self, x, self.Coefficients)
#    list_setter(self, x, self.Primary)
#    list_setter(self, x, self.BowCurveCoeff)
#    list_setter(self, x, self.SternCurveCoeff)
#    #list_setter(self, x, self.Areas)
#    #list_setter(self, x, self.list_SAC)
#    return
#    
#
#
#class DList(list):
#    def __init__(self):
#        pass
#    def getpoint(self, pt):
#        return self[0].getpoint(pt)
#        
#        
#class DesignTree(lp.States):
#    """Holds the tree of design states
#    """
#    #def __init__(self, states, parent=None, children=None):
#    #list.__init__(self, states)
#    #        self.parent = parent
#    #        if children is None:
#    #            self._children = []
#    #        else:
#    #            self.children = children
#        
#    def get_latest_child(self):
#        lc = len(self.children)
#        #nst = len(self)
#        if lc == 0:
#            return lp.States(self.states)#, nst
#            #return lp.States(self[:])#, nst
#        else:
#            return self.get_child()
#    
#    def __call__(self, x):
#        """we could inherit this call from States
#        """
#        d = []
#        #for st in self: #inherit from list
#        for st in self.states:
#            d.append(st(x))
#        return d
#    
#    def __len__(self):
#        return len(self.states)
#    
#            
#class CompoundCurve(object):
#    def __init__(self):
#        pass        

class C3PartCurve(object):
    """
    design_states:
        the state
        to be used with each call to a compound rule
    """
    def __init__(self,
                 length, l1,l2,l3,
                 total_area, a1,a2,a3,
                 Xc, x1,x2,x3,
                 max_height):
        self.length         = length
        self.entrance_len   = l1
        self.mid_len        = l2
        self.run_len        = l3
        self.dl1 = lp.Variable('dl1')
        self.dl2 = lp.Variable('dl2')
        #self.dl3 = lp.Variable('dl3')
        #self.dl4 = lp.Variable('dl4')
        self.total_area     = total_area
        self.entrance_area  = a1
        self.mid_area       = a2
        self.run_area       = a3
        self.da1 = lp.Variable('da1') #area dummies
        self.da2 = lp.Variable('da2')
        self.da3 = lp.Variable('da3') #a_fwd ~ a_aft
        self.da4 = lp.Variable('da4') #a_fwd ~ a_aft
        self.da5 = lp.Variable('da5') #a_fwd ~ a_aft
        self.da6 = lp.Variable('da6') #a_fwd ~ a_aft
        self.Xc             = Xc
        self.fwd_Xc         = x1
        self.mid_Xc         = x2
        self.run_Xc         = x3
        self.dx1 = lp.Variable('dx1') #length dummies
        self.dx2 = lp.Variable('dx2')
        #
        self.sc1 = lp.Variable('sc1') #area*length dummies
        self.sc2 = lp.Variable('sc2')
        self.sc3 = lp.Variable('sc3')
        self.sct = lp.Variable('sct')
        #
        self.dx3 = lp.Variable('dx3')
        self.dx4 = lp.Variable('dx4')
        self.dx5 = lp.Variable('dx5')
        self.dx6 = lp.Variable('dx6')
        #
        # Harries Poly Centroid Approx
        #
        self.xstart = lp.Variable('xstart')
        self.xend = lp.Variable('xend')
        self.xrange = lp.Variable('xrange')
        self.ystart = lp.Variable('ystart')
        self.yend = lp.Variable('yend')
        self.yrange = lp.Variable('yrange')
        
        self.max_height     = max_height
    
    
    
    def C3_length_heuristic(self, design_states, niter=3):
        """
        lengths of entrance and run hull protions
        must not take up the whole hull.
        
        NOTE:  this is where the midbody flat section comes from.
            i.e. SAC_mid_len is implicitly defined here
            when we add C3_length below
        """
        lwl             = self.length
        entrance_len    = self.entrance_len
        mid_len         = self.mid_len
        run_len         = self.run_len
        #loc_bfc             = self.loc_bfc
        #loc_sfc             = self.loc_sfc
        s       = design_states.states
        statesi = copy.copy(s)
        
        vars = [lwl,entrance_len,mid_len,run_len]
        
        states = (design_states * (ia(.29,.7),lwl,entrance_len))
        #states = (states * (ia(.25,.4),lwl,mid_len))
        states = (states * (ia(.29,.7),lwl,run_len))
        
        
        design_states.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = design_states.clean_states(states.states,[])
        return states

        
        
    def C3_length(self, design_states, niter=3):
        """ generic 3 part curve function: called with any section
        to ensure that the __quantities__ of the 3 sections
        always equal the total quantity of the vessel
        """
        total_length    = self.length
        entrance_len   = self.entrance_len
        mid_len        = self.mid_len
        run_len        = self.run_len
        c1 = self.dl1
        c2 = self.dl2
        
        vars = [total_length, entrance_len, 
                mid_len, run_len,
                c1,
                c2]

        
        s = design_states.states
        statesi = copy.copy(s)
        dlist = []

        states = (design_states + (entrance_len, mid_len,c1))
        for i in range(niter):
            states = (states - (total_length,c1,run_len))
            states = (states + (entrance_len, mid_len,c1))
            #
            #defacto in this version!
            #states = (states == (SAC_mid_len,lfsac) )
            #
        
        design_states.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars,
                         dlist=dlist)
        states = design_states.clean_states(states.states,dlist)
        return states
        
    
    
    
    def C3_area_heuristic(self, design_states, niter=3):
        """
            entrance_area   = ia(.3,.7)*total_area
            run_area        = ia(.3,.7)*total_area
            
            entrance_area ~ run_area  :: TODO: do I really need that heuristic?
        """
        total_area      = self.total_area
        entrance_area   = self.entrance_area
        mid_area        = self.mid_area
        run_area        = self.run_area
        c1 = self.da3
        c2 = self.da4
        c3 = self.da5
        c4 = self.da6
        
        s       = design_states.states
        statesi = copy.copy(s)
        
        vars = [total_area, entrance_area, mid_area,run_area,
                c1,c2,c3,c4]
        
        
        states = (design_states * (ia(.3,.7),total_area,entrance_area))
        for i in range(niter):
            states = (states * (ia(.3,.7),total_area,entrance_area))
            states = (states * (ia(.3,.7),total_area,run_area))
            #*****************************************************
            # STRONG HEURISTIC here:
            # entrance_area ~ run_area
            #
            states = (states * (ia(1.25,1.25),entrance_area,c1))
            states = (states * (ia(1.25,1.25),run_area,c2))
            states = (states * (ia(.75,.75),  entrance_area,c3))
            states = (states * (ia(.75,.75),  run_area,c4))
            states = (states <= (entrance_area, c2))
            states = (states <= (run_area, c1))
            states = (states >= (entrance_area, c4))
            states = (states >= (run_area, c3))
            
        design_states.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = design_states.clean_states(states.states,[])
        return states
        
    
        
    def C3_area(self, design_states, niter=3):
        """ 
        Rule:
            
            entrance_area + mid_area = total_area - run_area .
        
            self.da1 : is the facilitating variable
        """
        total_area      = self.total_area
        entrance_area   = self.entrance_area
        mid_area        = self.mid_area
        run_area        = self.run_area
        c1 = self.da1
        
        vars = [total_area, entrance_area, 
                mid_area, run_area,
                c1]
        
        s = design_states.states
        statesi = copy.copy(s)
        dlist = []

        states = (design_states + (entrance_area, mid_area,c1))
        for i in range(niter):
            states = (states - (total_area,c1,run_area))
            states = (states + (entrance_area, mid_area,c1))
            
        
        design_states.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars,
                         dlist=dlist)
        states = design_states.clean_states(states.states,dlist)
        return states
        
        
        
        
    def C3_mid_section_rule(self, 
                     design_states, niter=3):
        """Rule:
            states = (states * (max_height,mid_len,mid_area))  
        """
        mid_area        = self.mid_area
        mid_len         = self.mid_len
        max_height      = self.max_height
        
        s = design_states.states
        statesi = copy.copy(s)
        vars = [mid_area,
                mid_len,
                max_height]
        dlist = []
        
        states = (design_states * (max_height,mid_len,mid_area)) 
        
        design_states.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = design_states.clean_states(states.states,dlist)
        return states
        
        
        
    def C3_centroid(self, 
                     design_states, 
                     niter=3):
        """Rules:
            
            mid_Xc = entrance_len + ia(.5,.5)*mid_len
            
            fwd_Xc = ia(.64,.68)*entrance_len
            fwd_Xc < entrance_len
            
            run_Xc = b3 + b5
                   = (mid_len + entrance_len) + ia(.32,.36)*b4
                   = (mid_len + entrance_len) + ia(.32,.36)*(lwl - b3)
                   
                   could shorten to:
                   = (mid_len + entrance_len) + ia(.32,.36)*run_len
                       
            run_Xc > mid_len + entrance_len
            
        """
        lwl             = self.length
        LCG             = self.Xc
        total_area      = self.total_area
        
        entrance_area   = self.entrance_area
        mid_area        = self.mid_area
        run_area        = self.run_area
        
        entrance_len    = self.entrance_len
        mid_len         = self.mid_len
        run_len         = self.run_len
        
        fwd_Xc          = self.fwd_Xc
        mid_Xc          = self.mid_Xc
        run_Xc          = self.run_Xc
        
        max_height      = self.max_height
        
        c1 = self.dx1       #c1 = sc1 + sc2
        c2 = self.dx2       #c2 = c1 + sc3 = sc1+sc2+sc3
        
        sc1 = self.sc1      #sc1 = entrance_area * fwd_Xc
        sc2 = self.sc2      #sc2 = mid_area * mid_Xc
        sc3 = self.sc3      #sc3 = run_area * run_Xc
        sct = self.sct      #sc4 = total_area * LCG
        
        b1  = self.dx3
        b2  = self.dx4
        b3  = self.dx5
        b5  = self.dx6
        
        
        s = design_states.states
        statesi = copy.copy(s)
        
        vars = [lwl,
                LCG,
                total_area,
                entrance_area,
                mid_area,
                run_area,
                entrance_len,
                mid_len,
                run_len,
                fwd_Xc,
                mid_Xc,
                run_Xc,
                max_height,
                c1,c2,
                sc1,sc2,sc3,sct,
                b1,b2,b3,b5]
        
        #        c1 = lp.Variable('c1')      #c1 = sc1 + sc2
        #        self.set_val(c1,None)
        #        c2 = lp.Variable('c2')      #c2 = c1 + sc3 = sc1+sc2+sc3 
        #self.set_val(c2,None)
        #c3 = lp.Variable('c3')
        #self.set_val(c3,None)
        
        #        sc1 = lp.Variable('sc1')    #sc1 = entrance_area * fwd_Xc
        #        self.set_val(sc1,None)
        #        sc2 = lp.Variable('sc2')    #sc2 = mid_area * mid_Xc
        #        self.set_val(sc2,None)
        #        sc3 = lp.Variable('sc3')    #sc3 = run_area * run_Xc
        #        self.set_val(sc3,None)
        #        sct = lp.Variable('sct')    #sct = total_area * LCG
        #        self.set_val(sct,None)
        
        #        a1 = lp.Variable('a1')
        #        self.set_val(a1,None)
        #        a2 = lp.Variable('a2')
        #        self.set_val(a2,None)
        
        
        #        b1 = lp.Variable('b1')  #b1 = 0.5*mid_len
        #        self.set_val(b1,None)
        #        b2 = lp.Variable('b2')  #b2 = b1+entrance_len
        #        self.set_val(b2,None)
        #        b3 = lp.Variable('b3')
        #        self.set_val(b3,None)   #b3 = mid_len + entrance_len
        #        #b4 = lp.Variable('b4')
        #        #self.set_val(b4,None)   #b4 = (lwl - b3)  #replace with run_len
        #        b5 = lp.Variable('b5')
        #        self.set_val(b5,None)   #b5 = ia(.32,.36)*b4
        
        #f1 = lp.Variable('f1')
        #self.set_val(f1,None)
        
        #        dlist = [c1,c2,c3,
        #                 sc1,sc2,sc3,sct,
        #                 a1,a2,
        #                 b1,b2,b3,b4,b5,
        #                 f1]
        dlist = []

        
        states = (design_states <= (run_Xc, lwl))
        states = (states <= (mid_Xc, run_Xc))
        states = (states <= (fwd_Xc, mid_Xc))
        states = (states * (ia(.5,.5), mid_len, b1) ) 
        for i in range(niter):
            
            #
            # rule: mid_Xc = entrance_len + ia(.5,.5)*mid_len
            #states = (states * (ia(.5,.5), mid_len, b1) )
            states = (states + (b1,entrance_len,b2) )
            states = (states == (mid_Xc,b2) )
            #
            #
            # fwd_Xc = ia(.64,.68)*entrance_len
            states = (states * (ia(.30,.68), entrance_len, fwd_Xc))
            #
            # entrance Xc < entrance_len
            states = (states <= (fwd_Xc,entrance_len) )        
            #
            #
            # run Xc must be in the length of the run somewhere:
            states = (states + (mid_len, entrance_len, b3) )
            states = (states >= (run_Xc,b3))
            #
            # rule:
            #  run_Xc = entrance_len + mid_len + ia(.32,.34)*run_len
            #
            #subrule1:  (lwl - b3) => b4
            #states = (states - (lwl, b3, b4) ) #b4 == run_len
            #subrule2: (ia(.3,.37) * b4 ) => b5
            #states = (states * (ia(.32,.36),b4,b5))  #TLM Feb 4
            states = (states * (ia(.32,.55),run_len,b5))#TLM Feb 4 change to run_len directly
            states = (states + (b3,b5,run_Xc))
            #
            #
            #
            #        # Sum ( A*Xc )local = (total_area * LCG)
            #        # SAC - area- == some real total_areaume
            #"""
            states = (states * (total_area,LCG,sct) ) #total total_area * LCG = sct
            #states = (states == (c2,sct))             #equiv of area barycenters #shortcut!!  eliminate TLM FEB 4
            #
            states = (states * (entrance_area, fwd_Xc, sc1))
            states = (states * (mid_area,      mid_Xc, sc2) )
            states = (states * (run_area,      run_Xc, sc3) )
            states = (states + (sc1,sc2,c1) )
            #states = (states + (c1,sc3,c2) ) #sum of (local_areas * center of areaXc_local)
            states = (states + (c1,sc3,sct) ) #shortcut!!  TLM FEB 4
            
            states = (states * (ia(.5,.5), mid_len, b1) )
            #states = (states * (ia(.1,.3), run_len, run_Xc))
        
            
        design_states.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = design_states.clean_states(states.states,dlist)
        return states
        
        
        
    def C3_harries_centroid_rule(self, 
                     design_states, 
                     niter=3):
        """Rules:
            
            TOBD: polynomial box rule for centroids
            
            NOTE:  DO THIS FOR LCG
            NOT for individual centroids
            
            not withouth n-ary constraints!
            
            issue:  Y values in the formula
            
        """
        lwl             = self.length
        LCG             = self.Xc
        total_area      = self.total_area
        
        entrance_area   = self.entrance_area
        mid_area        = self.mid_area
        run_area        = self.run_area
        
        entrance_len    = self.entrance_len
        mid_len         = self.mid_len
        run_len         = self.run_len
        
        fwd_Xc          = self.fwd_Xc
        mid_Xc          = self.mid_Xc
        run_Xc          = self.run_Xc
        
        max_height      = self.max_height
 
        xi              = self.xstart
        xe              = self.xend
        xm              = self.xrange
        yi              = self.ystart
        ye              = self.yend
        ym              = self.yrange
        
        
        
        
        s = design_states.states
        statesi = copy.copy(s)
        
        vars = [lwl,
                LCG,
                total_area,
                entrance_area,
                mid_area,
                run_area,
                entrance_len,
                mid_len,
                run_len,
                fwd_Xc,
                mid_Xc,
                run_Xc,
                max_height,
                xi,xe,xm,
                yi,ye,ym]
        
        dlist = []

        #xb = x1 = x2 = 0.
        #xe = x5 = x4 = design_states(lwl)
        
        

        
        states = (design_states == (xi, ia(0.,0.) ))
        states = (states == (xe, lwl))
        states = (states == (xm, lwl))
        
        for i in range(niter):
            pass

            
        design_states.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = design_states.clean_states(states.states,dlist)
        return states
        
        

 
            
class HullDesignNode(lp.States):
    """
        Key Insight: This class should intialize the design space
        and pass that initialized state off to a design tree
            -desisions made in that tree act on [state] that is fed
                back into this class
            -But new states are sent from here to the design tree
            -Which can do what it wants with them
            
        Two ways to go:
            -thin optimization:  reify down to one design point
                This is what is currently done.
            -topological optimization: compute on the design space
                holistically.  Keep the entire space, just subdivide into
                regions of:
                    -feasabiliyt
                    -infeasibility
                    -not yet determined
            
        
        Improve by:
            -each fundemental curve of the bare hull
            should be its own particular design rules class
            -the bare hull design class should integrate them
            
            
        currently used variables:
            
            
        self.hull_list = [self.vol,    # displacement
                          self.lwl,    # waterline length
                          self.bwl,    # max breadth
                          self.draft,  # Max draft
                          self.LCG]    # Longitudinal Center of Buoyancy
                          
        self.Coefficients = [self.Cb,       # block 
                             self.Cmidshp,  # midship 
                             self.Cp,       # prismatic
                             self.Cwp,      # water plane
                             self.Ccp,      # center plane
                             self.Clcg]     # Coeff long. center of gravity
        
        self.Areas = [self.Amsh,  # midship
                      self.Awp,   # water plane
                      self.Acp]   # center plane
                      
        self.Primary = [self.Cb,  #see above
                        self.lwl,
                        self.vol,
                        self.draft,
                        self.bwl,
                        self.LCG]
        
        self.BowCurveCoeff = [self.Abfc,  # area bow fairness curve
                              self.Cbfc,  # Cbfc = Abfc/(bbfc*dbfc)
                              self.bbfc,  # breadth at bow fairness curve (bfc)
                              self.dbfc]  # depth at bfc
        
        self.SternCurveCoeff = [self.Asfc,  # stern fairness curve (sfc)
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
    """
    SMALL = 1.e-3
    def __init__(self, strd_coeff=True, parent=None, verbose=False):
        self.nstates = 0
        self.propagate = True
        self._verbose = verbose
        self.valid_designs = True
        self._states = None
        super(HullDesignNode, self).__init__( states = [self._setup_()] )
        self.nstates = len(self.states)
        self.keys = [el for el in self.states[0].values.keys() 
                    if isinstance(el, lp.Variable)]
        self.make_lists()
        self.parent = parent
        self.SAC_rules = C3PartCurve(self.lwl, 
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
                                     self.SAC_run_Xc,
                                     self.Amsh)
        
        if strd_coeff:
            self.ini_coeff()
            self.ini_lengths()
            self.ini_areas()
            #self.ini_locations() #fairness curve locations
            self.ini_flats()
            #self.compute_stern_fairness_section_ini()
            

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
                          
        #        self.Coefficients = [#self.Cb,
        #                             self.Cmidshp,
        #                             self.Cp,
        #                             self.Cwp,
        #                             self.Ccp,
        #                             self.Cbfc,
        #                             self.Csfc]
                          
        self.Coefficients = [self.Cb,
                             self.Cmidshp,
                             self.Cp,
                             self.Cwp,
                             self.Ccp,
                             self.Clcg]
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
        
                         
        self.alists=self.Coefficients +\
                    self.Areas+\
                    self.Primary+\
                    self.list_SAC
        
        self.keys = [el for el in self.keys if el not in self.BowCurveCoeff]
        self.keys = [el for el in self.keys if el not in self.SternCurveCoeff]
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
                                          self.ini_lengths,
                                          self.compute_LCG_coefficient])
                                          
        self._bwl     = lp.Variable('bwl',
                                    arcs=[self.compute_Cb,
                                          self.compute_Cwp,
                                          self.compute_midship_coefficient,
                                          self.compute_FOWL,
                                          self.compute_FOS])#,
                                          #self.derive_bow_fairness_section_loc_rule_fromDWL,
                                          #self.derive_stern_fairness_section_loc_rule_fromDWL])
        self._draft   = lp.Variable('draft',
                                    arcs=[self.compute_Cb,
                                          self.compute_Ccp,
                                          #self.compute_drafts,
                                          self.compute_FOK,
                                          self.compute_FOS,
                                          self.compute_midship_coefficient])#,
                                          #self.derive_bow_fairness_section_loc_rule_fromCPKeel,
                                          #self.derive_stern_fairness_section_loc_rule_fromCPKeel])
        
        self._dsmax   = lp.Variable('dsmax',
                                    arcs=[self.compute_midship_coefficient])#,
                                          #self.compute_drafts])
        
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
                                          self.sac_XC_ini,
                                          self.compute_LCG_coefficient])  
        self._Clcg   = lp.Variable('Clcg',
                                    arcs=[self.compute_LCG,
                                          self.compute_SAC_section_properties,
                                          self.sac_XC_ini,
                                          self.compute_LCG_coefficient]) 
        # 3 part curve use only:
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
                                          self.compute_FOSAC])#,
                                          #self.derive_bow_fairness_section_loc_rule_fromSAC,
                                          #self.derive_stern_fairness_section_loc_rule_fromSAC])
        
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
        
        #        self._lfos   = lp.Variable('lfos',
        #                                    arcs=[self.compute_FOS,
        #                                         self.compute_flat_relations])
        
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
                                                     self.ini_areas])#,
                                                     #self.derive_bow_fairness_section_loc_rule_fromSAC,
                                                     #self.derive_bow_fairness_section_loc_rule_fromDWL,
                                                     #self.derive_bow_fairness_section_loc_rule_fromCPKeel])
        self._SAC_mid_len = lp.Variable('SAC_mid_len',
                                             arcs = [self.compute_SAC_section_properties,
                                                     self.ini_lengths,
                                                     self.ini_areas,
                                                     self.compute_flat_relations])
        self._SAC_run_len = lp.Variable('SAC_run_len',
                                             arcs = [self.compute_SAC_section_properties,
                                                     self.ini_lengths,
                                                     self.ini_areas])#,
                                                     #self.derive_stern_fairness_section_loc_rule_fromSAC,
                                                     #self.derive_stern_fairness_section_loc_rule_fromDWL,
                                                     #self.derive_stern_fairness_section_loc_rule_fromCPKeel])
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
                              #self._dsmax    : None,
                              self._vol      : None,
                              self._LCG      : None,
                              self._Clcg     : None,
                              self._Cb       : None,
                              self._lwl      : None,
                              self._bwl      : None,
                              self._Awp      : None,
                              self._Amsh     : None,
                              self._Cwp      : None,
                              self._Cmidshp  : None,
                              self._Cp       : None,
                              #self._lfos     : None,
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
    
    def AC_revise(self, print_=False, maxit = 8):
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
    
    def equal(self, c1,c2, propagate=None):
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
        if propagate is None:
            propagate = self.propagate
        
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
        if propagate:
            self.AC_revise()
        self.delete_redundant_states()
        if self.nstates == 0:
            self.valid_designs = False
        return
    
    def clean_states(self, states, special=None):
        """TODO: Update for many states...
        Launguage is not quite up to par yet
        making do, clean extranious vars thatcrop up
        during rule computations
        
        special : temporary variables that should(?) be deleted
        """
        #return states
        for state in states:
            snm = [el.name for el in special]
            dkeys = []
            for key in state.values:
                #
                # Using this kills the propogation!
                # which will be noticable in greatly increased 
                # iteration (AC_revise always going to max)
                #
                #if not isinstance(key, lp.Variable):
                #    dkeys.append(key)
                #
                # End of Propogation Killer!
                #
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
        #l1 = len(states) #new
        #l2 = len(states_old) #old
        #iters = np.
        #if l1==l2:
       # c = 0
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
    
    #    @property
    #    def dsmax(self):
    #        """midship draft dso
    #        TODO: rename
    #        """
    #        return self._dsmax
    #    @dsmax.setter
    #    def dsmax(self, dsmax):
    #        self.equal(self.dsmax,dsmax)
        
        
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
    def Clcg(self):
        return self._Clcg
    @Clcg.setter
    def Clcg(self, Clcg):
        print 'setting Clcg'
        self.equal(self.Clcg,Clcg)
    
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
        """Area Water Plane"""
        return self._Awp
    @Awp.setter
    def Awp(self, Awp):
        self.equal(self.Awp,Awp)
        
        
    @property
    def Cwp(self):
        """Water Plane Coefficient"""
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
        """Center Plane Coefficient"""
        return self._Ccp
    @Ccp.setter
    def Ccp(self, Ccp):
        self.equal(self.Ccp,Ccp)
        
    @property
    def Amsh(self):
        """Misdship Area"""
        return self._Amsh
    @Amsh.setter
    def Amsh(self, Amsh):
        self.equal(self.Amsh,Amsh)

        
    @property
    def Cmidshp(self):
        """Midship Coefficient"""
        return self._Cmidshp
    @Cmidshp.setter
    def Cmidshp(self, Cmidshp):
        self.equal(self.Cmidshp,Cmidshp)

        
    @property
    def Cp(self):
        """Prismatic Coefficient"""
        return self._Cp
    @Cp.setter
    def Cp(self, Cp):
        self.equal(self.Cp,Cp)
    
    @property
    def lfwl(self):
        """Length of Flat of Water Line
        lfwl-lfsac = len flat of side
        """
        return self._lfwl
    @lfwl.setter
    def lfwl(self, lfwl):
        self.equal(self.lfwl,lfwl)
        
    #    @property
    #    def lfos(self):
    #        """Length of Flat of Side
    #        lfwl-lfsac = len flat of side"""
    #        return self._lfos
    #    @lfos.setter
    #    def lfos(self, lfos):
    #        self.equal(self.lfos,lfos)
    
    @property
    def lfsac(self):
        """Length of Flat of SAC"""
        return self._lfsac
    @lfsac.setter
    def lfsac(self, lfsac):
        self.equal(self.lfsac,lfsac)
        
    @property
    def lfcp(self):
        """Length of Flat of Center Plane
        lfcp - lfsac = flat of bottom extension
        """
        return self._lfcp
    @lfcp.setter
    def lfcp(self, lfcp):
        self.equal(self.lfcp,lfcp)
        
    #    @property
    #    def Afos(self):
    #        """not used"""
    #        return self._Afos
    #    @Afos.setter
    #    def Afos(self, Afos):
    #        self.equal(self.Afos,Afos)
        
    ##
    ## Bow Fairness curve
    ##
    @property
    def loc_bfc(self):
        """location of bow fairness curve"""
        return self._loc_bfc
    @loc_bfc.setter
    def loc_bfc(self, loc_bfc):
        self.equal(self._loc_bfc,loc_bfc)
    @property
    def bbfc(self):
        """mean at bow fairness curve
        """
        return self._bbfc
    @bbfc.setter
    def bbfc(self, bbfc):
        self.equal(self.bbfc,bbfc)
    @property
    def dbfc(self):
        """draft at bow fairness curve"""
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
    ##********************************************************************
    ##
    """
    def rule_parser(self, sexprs ):
        ops = [op]
        for el in sexprs:
            if op
    """
    ##
    ##********************************************************************
    ##
    def compute_FOSAC(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        #self.states = self.flat_of_sac_constraint()
        self.states = self.SAC_rules.C3_mid_section_rule(self)
        return
    def compute_LCG(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.LCG_constraints()
        return
    def compute_LCG_coefficient(self, **args):
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.states = self.CLCG_constraints()
        return
    
    #    def compute_drafts(self, **args):
    #        """
    #            dsmax == draft
    #        """
    #        if args:
    #            for arg in args:
    #                key = getattr(self, arg) 
    #                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
    #        self.states = self.relate_drafts()
    #        return
        
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
        """
        Not used anymore 
        
        define lfos == lfwl - lfsac
        define lfob == lfcpk - lfsac
        """
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        #self.states = self.flat_of_side_constraint()
        return
        
    def compute_FOB(self, **args):
        """
        Not used anymore 
        
        define lfos == lfwl - lfsac
        define lfob == lfcpk - lfsac
        """
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        #self.states = self.flat_of_bottom_constraint()
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
        self.states = self.SAC_rules.C3_length(self)
        self.states = self.SAC_rules.C3_length_heuristic(self)
        self.states = self.SAC_rules.C3_area(self)
        self.states = self.SAC_rules.C3_area_heuristic(self)
        self.states = self.SAC_rules.C3_centroid(self)
        
        #        self.states = self.SAC_rules.C3_curve(self.vol,
        #                                              section = self.SAC_entrance_area,
        #                                              sections = [  self.SAC_entrance_area,
        #                                                            self.SAC_mid_area,
        #                                                            self.SAC_run_area])
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
        #        self.states = self.SAC_run_area_consistency()
        #        #self.states = self.SAC_entrance_len_consistency()
        #        self.states = self.constrain_SAC_Xc() #SAC_run_area
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
        #self.SAC_length_ini()
        self.states = self.SAC_rules.C3_length(self)
        self.states = self.SAC_rules.C3_length_heuristic(self)
        return
    
    def ini_flats(self):
        self.states = self.initialize_flats()
        
    def ini_areas(self):
        """
            make sure all areas are >= 0.
            
            how abou limits on how small a SAC section can be?
        """
        #self.states = self.SAC_run_area_ini()
        self.states = self.SAC_rules.C3_area(self)
        self.states = self.SAC_rules.C3_area_heuristic(self)
        return
    
    def sac_XC_ini(self):
        self.states = self.SAC_rules.C3_mid_section_rule(self)
        self.states = self.SAC_rules.C3_centroid(self)
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

        #states = (self * (lwl,br,c1) * (c1,draft,c2) / (vol,c2,Cb) ) 
        states = (self * (lwl,br,c1) / (vol,c2,Cb)  * (c1,draft,c2) ) #heuristic:  compute dummy termed terms last! This gets draft
        states = (states * (lwl,br,c1) * (c1,draft,c2) / (vol,c2,Cb) ) # but it's good to mix things up - this gets Cb
        #states = (states * (lwl,br,c1) * (c1,draft,c2) / (vol,c2,Cb) ) 
        states = (states * (lwl,br,c1) / (vol,c2,Cb)  * (c1,draft,c2) )
        
        
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
        draft   : simplify life a bit and just use a single draft
        Cmidshp : Midship Coefficient
        
                Amsh/(Bwl*draft) = Cmidshp
        """
        bwl = self.bwl
        #dsmax = self.dsmax
        draft = self.draft
        Amsh = self.Amsh
        Cmidshp = self.Cmidshp
        s       = self.states
        statesi = copy.copy(s)
        
        c1 = lp.Variable('c1')
        self.set_val(c1,None)
        
        #states = (self == (draft,dsmax))
        states = (self * (bwl,draft,c1) / (Amsh,c1,Cmidshp))
        states = (states * (bwl,draft,c1) / (Amsh,c1,Cmidshp))
                          
        vars = [bwl,Amsh,Cmidshp,draft]
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[c1])
        return states
        
        
        
    def prismatic_coefficient(self):
        """
        Cb / Cmidshp = Cp
        
        test:
        cb = ia(0.95, 0.95)
        cmidshp = ia(0.989583333336, 0.989583333337)
        cp = ia(0.96, 0.96)
        """
        Cb      = self.Cb
        Cmidshp = self.Cmidshp
        Cp      = self.Cp
        
        s       = self.states
        statesi = copy.copy(s)
        
        states = (self / (Cb,Cmidshp,Cp) )
        states = (states * (Cmidshp,Cp,Cb) )
        
        vars = [Cb,Cmidshp,Cp]
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[])
        return states
    
    def LCG_constraints(self):
        """
        0.0 < LCG < lwl
        """
        lwl     = self.lwl
        LCG     = self.LCG
        
        s       = self.states
        statesi = copy.copy(s)
        states = (self <= (LCG,lwl))
        states = (states >= (LCG,ia(0.,0.)))
        
        vars = [lwl,LCG]
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[])
        return states
    
    
    
    def CLCG_constraints(self):
        """
        Clcg*lwl = LCG
        """
        #print 'CLCG_constraints'
        lwl     = self.lwl
        LCG     = self.LCG
        Clcg     = self.Clcg
        
        s       = self.states
        statesi = copy.copy(s)
        states = (self * (Clcg,lwl,LCG))
        
        vars = [lwl,LCG,Clcg]
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[])
        return states
    
    
    
        
        
        
    ##
    ## ------------------------------------------------- Curve Flat Interactions
    ##
    def initialize_flats(self, niter=3):
        """
        Definition:  Lfos is just the piece of the DWL which is both flat
        and not stationwise coincident with Lfsac
        SO I AM DELETING IT as it is described by the difference  lwl-lfsac
            Lfos : length of flat of side
            Lfsac + Lfos <= LfwL
            LfwL <= lwl
            Lfsac + Lfob <= Lfcp where LFlatOfoBottom was never defined!
        """
        lwl     = self.lwl
        lfwl    = self.lfwl
        lfsac   = self.lfsac
        lfcp    = self.lfcp
        
        s       = self.states
        statesi = copy.copy(s)
        vars = [lwl,lfwl,lfsac,lfcp]

        
        states = (self   <= (lfwl,  lwl )) #Feb 3 change!
        states = (states <= (lfcp,  lwl )) #Feb 3 change!
        states = (states <= (lfsac, lfcp)) #Feb 3 change!
        states = (states <= (lfsac, lfwl)) #Feb 3 change!
        
        states = (states >= (lfwl,  ia(0.,0.))) #Feb 3 change!
        states = (states >= (lfcp,  ia(0.,0.))) #Feb 3 change!
        states = (states >= (lfsac, ia(0.,0.))) #Feb 3 change!
            
        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[])
        return states
        
        
    def relate_flats(self, niter=3):
        """
        Definition:  Lfos is just the piece of the DWL which is both flat
        and not stationwise coincident with Lfsac
        SO I AM DELETING IT 
            as it is described by the difference  lfos = lwl-lfsac
            
        Definitions:
            Lfos: length of flat of side
            Lfsac + Lfos <= LfwL so lfos is getting deleted
            because actually Lfsac + Lfos == LfwL
            LfwL <= lwl
            Lfsac + Lfob <= Lfcp where LFlatOfoBottom was never defined!
            and actually Lfsac+Lfob == Lfcp
        """
        lwl         = self.lwl
        lfwl        = self.lfwl
        lfsac       = self.lfsac
        SAC_mid_len = self.SAC_mid_len
        lfcp        = self.lfcp
        s = self.states
        statesi = copy.copy(s)
        
        c = lp.Variable('c')
        self.set_val(c,None)
        c1 = lp.Variable('c1')
        self.set_val(c1,None)
        c2 = lp.Variable('c2')
        self.set_val(c2,None)
        
        vars = [lwl,lfwl,lfsac,lfcp,SAC_mid_len]

        
        states = (self <=   (lfwl, lwl)) #Feb 3 2017 change!
        states = (states <= (lfcp, lwl)) #Feb 3 change!
        for i in range(niter):
            states = (states <= (lfsac, lfcp)) #Feb 3 change!
            states = (states <= (lfsac, lfwl)) #Feb 3 change!
            states = (states == (lfsac, SAC_mid_len))
            #        for i in range(niter):
            #            states = (states <= (lfsac, lfwl))
            #            states = (states <= (lfsac, lfcp))

        self.set_updates(states = states.states,
                         states_old = statesi,
                         vars=vars)
        states = self.clean_states(states.states,[c,c1,c2])
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
    
    
    #    def flat_of_side_constraint(self, niter=6):
    #        """
    #            flat of side, FOS, is the distance:
    #                from : end of flat of sac
    #                  to : farthest edge of flat of WLarea
    #        
    #        lfos*bwl*draft < ( vol - lfsac*Amsh )
    #        <--> ?same? as:
    #        lfos < ( lwl - lfsac )
    #        
    #        other relations, for use with AFOS:
    #            AFOS*bwl < Length(FOS)*draft*bwl
    #            AFOS < Length(FOS)*draft
    #            Length(FOS)*draft > AFOS        
    #        """
    #        lwl     = self.lwl
    #        bwl     = self.bwl
    #        draft   = self.draft
    #        #lfos    = self.lfos
    #        vol     = self.vol
    #        lfsac   = self.lfsac
    #        Amsh    = self.Amsh
    #        s       = self.states
    #        statesi = copy.copy(s)
    #        
    #        #vars = [bwl,draft,lfos,vol,lfsac,Amsh,lwl]
    #        vars = [bwl,draft,vol,lfsac,Amsh,lwl]
    #  
    #  
    #        c1 = lp.Variable('c1')
    #        self.set_val(c1,None)
    #        c2 = lp.Variable('c2')
    #        self.set_val(c2,None)
    #        c3 = lp.Variable('c3')
    #        self.set_val(c3,None)
    #        c4 = lp.Variable('c4')
    #        self.set_val(c4,None)
    #        c5 = lp.Variable('c5')
    #        self.set_val(c5,None)
    #        
    #        states = ( (self * (bwl,draft,c1))  * (lfsac,Amsh,c2) )
    #        for i in range(niter):
    #            """#needed?"""
    #            # lfos*bwl*draft < ( vol - lfsac*Amsh )
    #            #        g = lp.Goal.both(
    #            #                lp.Goal.both(
    #            #                    lp.Goal.both(lp.Goal.both(lp.Goal.mulo(bwl,draft,c1),
    #            #                                              lp.Goal.mulo(lfsac,Amsh,c2)),
    #            #                                 lp.Goal.both(lp.Goal.subo(vol,c2,c3),
    #            #                                              lp.Goal.divo(c3,c1,c4))),
    #            #                    lp.Goal.both(lp.Goal.lto(lfos,c4),
    #            #                                 lp.Goal.both(lp.Goal.divo(c3,c1,c4),
    #            #                                              lp.Goal.subo(vol,c2,c3)))
    #            #                            ),
    #            #                 lp.Goal.both(lp.Goal.mulo(lfsac,Amsh,c2),
    #            #                              lp.Goal.mulo(bwl,draft,c1))
    #            #                              )
    #            #        s = g(s)[0]  #ok but needed?
    #    
    #     
    #            """#succinct but doesn't back prop"""
    #            """# lfos*bwl*draft < ( vol - lfsac*Amsh )"""            
    #            states = ( ( (states - (vol,c2,c3))  / (c3,c1,c4) ) <= (lfos,c4))
    #    
    #            
    #            #very important relations below 
    #            #(should be subsummed into relate_flats but relate flats is failing...):
    #            #or just re-write flats with these.
    #            #        """# lfsac  < (lwl - lfos)"""
    #            #        g = lp.Goal.both(lp.Goal.subo(lwl,lfos,c5),
    #            #                         lp.Goal.lto(lfsac,c5) )
    #            #        s = g(s)[0]
    #            #        
    #            
    #            """# lfos  < (lwl - lfsac)"""
    #            #        g = lp.Goal.both(lp.Goal.subo(lwl,lfsac,c5),
    #            #                         lp.Goal.lto(lfos,c5) )
    #            states = ( (states * (bwl,draft,c1))  * (lfsac,Amsh,c2) )
    #        
    #        self.set_updates(states = states.states,
    #                         states_old = statesi,
    #                         vars=vars)
    #        states = self.clean_states(states.states,[c1,c2,c3,c4,c5])
    #        return states

        
        
    
        
if __name__ == "__main__":
    def mysetter(ukanrenobj, myclass):
        myclass.states = ukanrenobj
        return 
    self = HullDesignNode(strd_coeff=True,verbose=False)
    
    #"""
    self.Amsh       = ia(10.,20.)
    self.bwl        = ia(5.,10.)
    self.Cmidshp    = ia(0.5,1.)
    
    print self(self.Amsh)
    #[ia(10.0, 20.0)]
    print self(self.bwl)
    #[ia(5.0, 10.0)]
    print self(self.Cmidshp)
    #[ia(0.5, 1.0)]
    print self(self.draft)
    #[ia(1.0, 8.0)]
    
    self.AC_revise()
    
    #"""
    
    self.draft = ia(-2.,4.)
    self.vol = ia(100.,200.)
    self.lwl = ia(8.,10.)
    self.bwl = ia(2.,4.) 
    
    self.Cp = ia(-1.,1.)
    #self.draft = ia(0.,20.)
    #self.Amsh = ia(400.,400.)
    
    #self.states = (self >= (self.SAC_mid_area, ia(0.,0.))).states

    c1 = lp.Variable('c1')
    c2 = lp.Variable('c2')

    st = self.block_coefficient()
    #self.compute_Cb()
    self.states = self.block_coefficient()
    
    def printCb():
        print 'self(self.Cb) = {}'.format(self(self.Cb))
        print 'self(self.vol) = {}'.format(self(self.vol))
        print 'self(self.lwl) = {}'.format(self(self.lwl ))
        print 'self(self.bwl) = {}'.format(self(self.bwl))
        print 'self(self.draft) = {}'.format(self(self.draft))
        return
        
        
    self = HullDesignNode(strd_coeff=True,verbose=False)
    self.draft = ia(-2.,4.)
    self.vol = ia(100.,200.)
    self.lwl = ia(8.,10.)
    self.bwl = ia(2.,4.) 
    
    printCb()
    print '\nPick a thin interval'
    #self.Cb = ia(.95,.95) #check out self.Cb = ia(.95,.95) later!
    self.Cb = ia(.96,.96)  
    printCb()
    print '\nPick another thin interval'
    #self.vol = ia(100.,100.)
    #self._verbose = True
    #self.vol = ia(125.,125.)
    self.vol = ia(153.6,153.6)
    printCb()