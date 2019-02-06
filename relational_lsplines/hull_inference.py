# -*- coding: utf-8 -*-
"""
Created on Sat Jul 23 18:44:37 2016

@author: lukemcculloch

help(lp.run_br)
"""



import numpy as np
import copy
from automatic_differentiation import ad
from interval_arithmetic import ia

import uKanren as lp


#def ops(tris,state,ops):
#    for pa, op in zip(tris,ops):
#        st=lp.run_br(pa[0],pa[1],pa[2], state, op)
#        
#    
#    return

def waterplane_coefficient(lwl,bwl,Awp,Cwp,s):
    """
    Awp     : wateplane area
    lwl     : waterline length
    bwl     : waterline breadth
    Cwp     : waterplane coefficient
    
            Awp/(Lwl*Bwl) = Cwp
    """
    c1 = lp.Variable('c1')

    state = copy.copy(s)
    state.values[c1]=None
    
    goal = lp.Goal.both(lp.Goal.mulo(lwl,bwl,c1),
                        lp.Goal.divo(Awp,c1,Cwp)
                        )
    
    state = goal(goal(state)[0])[0]
    del(state.values[c1])
    dkeys = []
    for key in state.values:
        if not isinstance(key, lp.Variable):
            dkeys.append(key)
    for key in dkeys:
        del(state.values[key])
    return state

def block_coefficient(lpp,br,draft,vol,Cb,s):
    """
    lpp     : length between perpendiculars or wL length
    br      : some breadth (e.g. breadth at the waterline)
    draft   : draft
    vol     : volumetric displacement/water weight
    
            Vol/(Lwl*Bwl*Draft) = Cb
    """
    c1 = lp.Variable('c1')
    c2 = lp.Variable('c2')

    state = copy.copy(s)
    state.values[c1]=None
    state.values[c2]=None
    
    goal = lp.Goal.both(lp.Goal.both(lp.Goal.mulo(lpp,br,c1),
                                     lp.Goal.mulo(c1,draft,c2)),
                        lp.Goal.divo(vol,c2,Cb)
                        )
#    def go(s):
#        s = lp.Goal.mulo(lpp,br,c1)(s)[0]
#        s = lp.Goal.mulo(c1,draft,c2)(s)[0]
#        s = lp.Goal.divo(vol,c2,Cb)(s)[0]
#        return s
#    s = go(state)
#    s=go(s)
#    s=go(s)
#    
#    def narrow(s): #?doesn't narrow
#        s = lp.Goal.divo(c1,br,lpp)(s)[0]
#        s = lp.Goal.divo(c2,draft,c1)(s)[0]
#        s = lp.Goal.mulo(Cb,c2,vol)(s)[0]
#        return s
#    s = go(state)
#    s=go(s)
#    s=go(s)
    #state = s
                                         
    #
    # search procedure:
    #
    #state = goal(state)
    state = goal(goal(goal(state)[0])[0])[0] #forward chaining!  maybe decorate or otherwise move outside

    del(state.values[c1])
    del(state.values[c2])
    dkeys = []
    for key in state.values:
        if not isinstance(key, lp.Variable):
            dkeys.append(key)
    for key in dkeys:
        del(state.values[key])
    return state #s

def midship_coefficient(bwl,dsmax,Amsh,Cmidshp,s):
    """
    Amsh    : Midship area
    Dsmax   : Draft at station of max section area
    Cmidshp    : Midship Coefficient
    
            Amsh/(Bwl*Dsmax) = Cmidshp
    """
    c1 = lp.Variable('c1')

    state = copy.copy(s)
    state.values[c1]=None
    
    goal = lp.Goal.both(lp.Goal.mulo(dsmax,bwl,c1),
                        lp.Goal.divo(Amsh,c1,Cmidshp)
                        )
    state = goal(goal(state)[0])[0]
    del(state.values[c1])
    dkeys = []
    for key in state.values:
        if not isinstance(key, lp.Variable):
            dkeys.append(key)
    for key in dkeys:
        del(state.values[key])
    return state
    
def prismatic_coefficient(lwl,bwl,draft,vol,Cb,
                            dsmax,Amsh,Cmidshp,
                            Cp,s):
    """
    Cb / Cmidshp = Cp
    """
    
    print Cb,       s.value_of(Cb)
    print Cmidshp,     s.value_of(Cmidshp)
    goal = lp.Goal.divo(Cb,Cmidshp,Cp)
    s = goal(s)[0]
    s = block_coefficient(lwl,bwl,draft,vol,Cb,s)
    s = midship_coefficient(bwl,dsmax,Amsh,Cmidshp,s)
    
    return s


class Hull(object):
    def __init__(self, state=None):
        self.graph = {}
        self.state = state or self._setup_()
#        if state:  NO - > how will you grab the correct pointers to variables without passing them in?
#            self.draft   = 
#            self.dsmax   =
        
        
    def _setup_(self):
        self.draft   = lp.Variable('draft')
        self.dsmax   = lp.Variable('dsmax')
        
        self.vol     = lp.Variable('vol')
        self.Cb      = lp.Variable('Cb')
        
        self.lwl     = lp.Variable('lwl')
        self.bwl     = lp.Variable('bwl')
        
        self.Awp     = lp.Variable('Awp') #waterplane area
        self.Cwp     = lp.Variable('Cwp')
        
        self.Amsh    = lp.Variable('Amsh') #Area midship - area of largest midship section
        self.Cmidshp    = lp.Variable('Cmidshp') #midship coeff
        
        self.Cp      = lp.Variable('Cp') 
        
        s  = lp.State(values={self.draft    : None,
                              self.dsmax    : None,
                              self.vol      : None,
                              self.Cb       : None,
                              self.lwl      : None,
                              self.bwl      : None,
                              self.Awp      : None,
                              self.Amsh     : None,
                              self.Cwp      : None,
                              self.Cmidshp  : None,
                              self.Cp       : None})
        return s
    
    
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
    
        state = copy.copy(s)
        state.values[c1]=None
        state.values[c2]=None
        
        goal = lp.Goal.both(lp.Goal.both(lp.Goal.mulo(lpp,br,c1),
                                         lp.Goal.mulo(c1,draft,c2)),
                            lp.Goal.divo(vol,c2,Cb)
                            )

        state = goal(goal(goal(state)[0])[0])[0] #forward chaining!  maybe decorate or otherwise move outside
    
        del(state.values[c1])
        del(state.values[c2])
        dkeys = []
        for key in state.values:
            if not isinstance(key, lp.Variable):
                dkeys.append(key)
        for key in dkeys:
            del(state.values[key])
        return state    
    
    def waterplane_coefficient(self):
        """
        Awp     : wateplane area
        lwl     : waterline length
        bwl     : waterline breadth
        Cwp     : waterplane coefficient
        
                Awp/(Lwl*Bwl) = Cwp
        """
        lwl = self.lwl
        bwl = self.bwl
        Awp = self.Awp
        Cwp = self.Cwp
        s = self.state
        
        c1 = lp.Variable('c1')
        state = copy.copy(s)
        state.values[c1]=None
        
        goal = lp.Goal.both(lp.Goal.mulo(lwl,bwl,c1),
                            lp.Goal.divo(Awp,c1,Cwp)
                            )
        
        state = goal(goal(state)[0])[0]
        del(state.values[c1])
        dkeys = []
        for key in state.values:
            if not isinstance(key, lp.Variable):
                dkeys.append(key)
        for key in dkeys:
            del(state.values[key])
        return state
        
        
    def midship_coefficient(self):
        """
        Amsh    : Midship area
        Dsmax   : Draft at station of max section area
        Cmidshp    : Midship Coefficient
        
                Amsh/(Bwl*Dsmax) = Cmidshp
        """
        bwl = self.bwl
        dsmax = self.dsmax
        Amsh = self.Amsh
        Cmidshp = self.Cmidshp
        s = self.state
        
        c1 = lp.Variable('c1')
        state = copy.copy(s)
        state.values[c1]=None
        goal = lp.Goal.both(lp.Goal.mulo(dsmax,bwl,c1),
                            lp.Goal.divo(Amsh,c1,Cmidshp)
                            )
        state = goal(goal(state)[0])[0]
        del(state.values[c1])
        dkeys = []
        for key in state.values:
            if not isinstance(key, lp.Variable):
                dkeys.append(key)
        for key in dkeys:
            del(state.values[key])
        return state
        
    def prismatic_coefficient(self):
        """
        Cb / Cmidshp = Cp
        """
        bwl = self.bwl
        draft = self.draft
        vol = self.vol
        Cb = self.Cb
        dsmax = self.dsmax
        Amsh = self.Amsh
        Cmidshp = self.Cmidshp
        Cp = self.Cp
        s = self.state
        
        print Cb,       s.value_of(Cb)
        print Cmidshp,  s.value_of(Cmidshp)
        goal = lp.Goal.divo(Cb,Cmidshp,Cp)
        s = goal(s)[0]
        s = block_coefficient(lwl,bwl,draft,vol,Cb,s)
        s = midship_coefficient(bwl,dsmax,Amsh,Cmidshp,s)
        
        return s 
        
if __name__=='__main__':

    draft   = lp.Variable('draft')
    dsmax   = lp.Variable('dsmax')
    
    vol     = lp.Variable('vol')
    Cb      = lp.Variable('Cb')
    
    lwl     = lp.Variable('lwl')
    bwl     = lp.Variable('bwl')
    
    Awp     = lp.Variable('Awp') #waterplane area
    Cwp     = lp.Variable('Cwp')
    
    Amsh    = lp.Variable('Amsh') #Area midship - area of largest midship section
    Cmidshp    = lp.Variable('Cmidshp') #midship coeff
    
    Cp      = lp.Variable('Cp') 
    
    s  = lp.State(values={draft : None,
                          dsmax : None,
                          vol   : None,
                          Cb    : None,
                          lwl   : None,
                          bwl   : None,
                          Awp   : None,
                          Amsh  : None,
                          Cwp   : None,
                          Cmidshp  : None,
                          Cp    : None})
                         
    #s = lp.Goal.eq(lpp,  ia(100.,120.)  )(s)[0]
    #s = lp.Goal.eq(lpp, 120.)(s)[0]
    #s = lp.Goal.eq(br, ia(30.,50.)      )(s)[0]#(s[0])#(s)#
    #s = lp.Goal.eq(br, ia(1.,218.75) )(s)[0]
    #s = lp.Goal.eq(br, 30.)(s)[0]
    
    s = lp.Goal.eq(vol,ia(2000.,300000.))(s)[0]
    s = lp.Goal.eq(Cb, ia(.5,.7)    )(s)[0]
    
    s = lp.Goal.eq(lwl,  ia(120.,120.) )(s)[0]
    s = lp.Goal.eq(bwl,  ia(20.,20.)   )(s)[0]
    s = lp.Goal.eq(draft,ia(20.,20.)   )(s)[0]
    
    s = lp.Goal.eq(Cwp, ia(.2,.4)     )(s)[0]
    #s = Block_Coefficient_long(lpp,br,draft,vol,Cb,s)
    #
    
    s = lp.Goal.eq(Cmidshp,  ia(.7,.99) )(s)[0]
    s = lp.Goal.eq(dsmax,  ia(25.,35.) )(s)[0]
    
    #s = lp.Goal.eq(Cp,  ia(.01,.99) )(s)[0]
    
    s = block_coefficient(lwl,bwl,draft,vol,Cb,s)
    s = waterplane_coefficient(lwl,bwl,Awp,Cwp,s)
    s = midship_coefficient(bwl,dsmax,Amsh,Cmidshp,s)
    #s = prismatic_coefficient(lwl,bwl,draft,vol,Cb,
    #                          dsmax,Amsh,Cmidshp,Cb,s)
    goal = lp.Goal.divo(Cb,Cmidshp,Cp)
    s = goal(s)[0]
    #for el in s.values:
    #    print el, s.values[el]
    print lwl,      s.value_of(lwl)
    print bwl,      s.value_of(bwl)
    print draft,    s.value_of(draft)
    print dsmax,    s.value_of(dsmax)
    print vol,      s.value_of(vol)
    print Cb,       s.value_of(Cb)
    print Awp,      s.value_of(Awp)
    print Cwp,      s.value_of(Cwp)
    print Amsh,     s.value_of(Amsh)
    print Cmidshp,     s.value_of(Cmidshp)
    print Cp,       s.value_of(Cp)
    
#    x = lp.Variable('x')
#    y = lp.Variable('y')
#    s = lp.State(values={x:None,y:None})
#    g = lp.Goal.both(lp.Goal.eq(y,ia(0.,1.)),
#                          lp.Goal.eq(y,x))
#    s = g(s)
    
    a = ia(0.,100.)
    b = ia(50.,150.)
    c = a&b
    
    
    
    self = Hull()
    s = self.block_coefficient()
    