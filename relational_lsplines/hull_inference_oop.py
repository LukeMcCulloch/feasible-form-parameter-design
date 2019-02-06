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
        self._draft   = lp.Variable('draft')
        self._dsmax   = lp.Variable('dsmax')
        
        self._vol     = lp.Variable('vol')
        self._Cb      = lp.Variable('Cb')
        
        self._lwl     = lp.Variable('lwl')
        self._bwl     = lp.Variable('bwl')
        
        self._Awp     = lp.Variable('Awp') #waterplane area
        self._Cwp     = lp.Variable('Cwp')
        
        self._Amsh    = lp.Variable('Amsh') #Area midship - area of largest midship section
        self._Cmidshp    = lp.Variable('Cmidshp') #midship coeff
        
        self._Cp      = lp.Variable('Cp') 
        
        s  = lp.State(values={self._draft    : None,
                              self._dsmax    : None,
                              self._vol      : None,
                              self._Cb       : None,
                              self._lwl      : None,
                              self._bwl      : None,
                              self._Awp      : None,
                              self._Amsh     : None,
                              self._Cwp      : None,
                              self._Cmidshp  : None,
                              self._Cp       : None})
        return s
    
    @property
    def lwl(self):
        return self._lwl
    @lwl.setter
    def lwl(self, lwl):
        goal = lp.Goal.eq(self.lwl, lwl)
        self.state = goal(self.state)[0]


    @property
    def bwl(self):
        return self._bwl
    @bwl.setter
    def bwl(self, bwl):
        goal = lp.Goal.eq(self.bwl, bwl)
        self.state = goal(self.state)[0]
        
        
    @property
    def draft(self):
        return self._draft
    @draft.setter
    def draft(self, draft):
        goal = lp.Goal.eq(self.draft, draft)
        self.state = goal(self.state)[0]
    
    @property
    def dsmax(self):
        return self._dsmax
    @dsmax.setter
    def dsmax(self, dsmax):
        goal = lp.Goal.eq(self.dsmax, dsmax)
        self.state = goal(self.state)[0]
        
    @property
    def vol(self):
        return self._vol
    @vol.setter
    def vol(self, vol):
        goal = lp.Goal.eq(self.vol, vol)
        self.state = goal(self.state)[0]
        
        
    @property
    def Cb(self):
        return self._Cb
    @Cb.setter
    def Cb(self, Cb):
        goal = lp.Goal.eq(self.Cb, Cb)
        self.state = goal(self.state)[0]
        
        
    @property
    def Awp(self):
        return self._Awp
    @Awp.setter
    def Awp(self, Awp):
        goal = lp.Goal.eq(self.Awp, Awp)
        self.state = goal(self.state)[0]
        
        
    @property
    def Cwp(self):
        return self._Cwp
    @Cwp.setter
    def Cwp(self, Cwp):
        goal = lp.Goal.eq(self.Cwp, Cwp)
        self.state = goal(self.state)[0]
        
        
    @property
    def Amsh(self):
        return self._Amsh
    @Amsh.setter
    def Amsh(self, Amsh):
        goal = lp.Goal.eq(self.Amsh, Amsh)
        self.state = goal(self.state)[0]
        
        
    @property
    def Cmidshp(self):
        return self._Cmidshp
    @Cmidshp.setter
    def Cmidshp(self, Cmidshp):
        goal = lp.Goal.eq(self.Cmidshp, Cmidshp)
        self.state = goal(self.state)[0]
        
        
    @property
    def Cp(self):
        return self._Cp
    @Cp.setter
    def Cp(self, Cp):
        goal = lp.Goal.eq(self.Cp, Cp)
        self.state = goal(self.state)[0]
    
    def compute_Cb(self, **args):
        """
        """
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.block_coefficient()
        return
    
    def compute_Cwp(self, **args):
        """
        """
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.waterplane_coefficient()
        return
        
    def compute_midship_coefficient(self, **args):
        """
        """
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.midship_coefficient()
        return
        
    def compute_prismatic_coefficient(self, **args):
        """
        """
        if args:
            for arg in args:
                key = getattr(self, arg) 
                self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
        self.state = self.prismatic_coefficient()
        return
        

    def set_generic(self, func, **args):
        def wrap(**args):
            if args:
                for arg in args:
                    key = getattr(self, arg) 
                    self.state = lp.Goal.eq(self.__dict__['_'+str(key)] , args[arg])(self.state)[0]
            self.state = func()
            return
        return wrap
        
    
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
        Cb = self.Cb
        Cmidshp = self.Cmidshp
        Cp = self.Cp
        s = self.state
        
        #c1 = lp.Variable('c1')
        #c2 = lp.Variable('c2')
        #state = copy.copy(s)
        #state.values[c1]=None
        
        #self.compute_Cb()
        #self.compute_midship_coefficient()
        

        goal = lp.Goal.divo(Cb,Cmidshp,Cp)
        state = goal(s)[0]
        
        #s = self.compute_Cb()
        #s = self.compute_midship_coefficient()
        
        #del(state.values[c1])
        #del(state.values[c2])
        dkeys = []
        for key in state.values:
            if not isinstance(key, lp.Variable):
                dkeys.append(key)
        for key in dkeys:
            del(state.values[key])
        
        return state
        
if __name__=='__main__':

    
    a = ia(0.,100.)
    b = ia(50.,150.)
    c = a&b
    
    
    
    self = Hull()
    #goal = lp.Goal.eq(self.lwl, ia(2000.,300000.))
    #self.state = goal(self.state)[0]
    self.vol =ia(2000.,300000.)
    self.lwl = ia(120.,120.)
    self.bwl = ia(20.,20.) 
    #self.Cb = ia(.5,.7) 
    self.draft = ia(20.,20.) 
    self.Cwp = ia(.2,.4) 
    self.Cmidshp =  ia(.7,.99) 
    self.dsmax =  ia(25.,35.)
    print self.state
    
    
    self.compute_prismatic_coefficient()
    print self.state
    
    self.compute_Cb()
    print self.state
    
    #self.compute_Cb(**{'Cb':ia(.5,.55)})
    self.compute_Cb()
    print self.state
    
    self.compute_Cwp()
    print self.state
    
    self.compute_midship_coefficient()
    print self.state
    
    self.compute_prismatic_coefficient()
    print self.state
    
    
def two(func):
    def fwrap():
        return '2', func()
    return fwrap

@two
def one():
    return '1!'
    

