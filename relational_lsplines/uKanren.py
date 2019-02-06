# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 11:17:25 2016

@author: lukemcculloch? No I am a mere(copy)cat:
From 
  https://codon.com/hello-declarative-world
To
  https://github.com/logpy/logpy

uKanren
    -minus addo, mulo, .. as at least addo kept recursing forever
    -exteded with binary mathematical relations and inerval unification
"""
#
#
#
from collections import OrderedDict
import itertools  #.permutations to have func with combi_explosivity
from inspect import isfunction
import operator #python native opps
import copy
from interval_arithmetic import ia
#
#
#
class Variable(object):
    def __init__(self, name='', nbrs=None, arcs=None):
        self.name = name
        self.nbrs = nbrs
        self.arcs = arcs

    def __repr__(self):
        return str(self.name)

    def inspect(self):
        print self.name
        return

#class Number(object):
#    def __init__(self, name='',value=None):
#        self.name = name
#        self.value = value
#
#    def __repr__(self):
#        return str(self.value)
#
#    def inspect(self):
#        print self.value
#        return

class States(object):
    
    def __init__(self, states=None):
        if isinstance(states, State):
            self.states = [states]
        elif isinstance(states, list):
            self.states = states
        elif isinstance(states, tuple):
            self.states = list(states)
            
    def __repr__(self):
        return 'States({})'.format(self.states)
    
    def __str__(self):
        return 'States({})'.format(self.states)
        
    def __add__(self, v):
        """addo(v1,v2,v3)(state)
        """
        states = []
        for s in self.states:
            sn = Goal.addo(v[0],v[1],v[2])(s)
            for sni  in sn:
                states.append(sni)
        return States(states)
        
    def __sub__(self, v):
        """subo(v1,v2,v3)(state)
        """
        states = []
        for s in self.states:
            sn = Goal.subo(v[0],v[1],v[2])(s)
            for sni  in sn:
                states.append(sni)
        return States(states)
        
    def both(self, a,b):
        states = []
        for s in self.states:
            sn = Goal.both(a,b)(s)
            for sni  in sn:
                states.append(sni)
        sf = []
        for s in states:
            sn = Goal.both(a,b)(s)
            for sni  in sn:
                sf.append(sni)
        return States(sf)
        
    
    #    def goaleval(goal):
    #        if callable(goal):  # goal is already a function like eq(x, 1)
    #            return goal
    #        if isinstance(goal, States):
    #            return goal.states
    #    
    #    def conj(self, stuff):
    #        """recursive both??
    #        """
    #        for goal in stuff:
    #            
    #        states = []
    #        for s in self.states:
    #            sn = Goal.both(a,b)(s)
    #            for sni  in sn:
    #                states.append(sni)
    #        sf = []
    #        for s in states:
    #            sn = Goal.both(a,b)(s)
    #            for sni  in sn:
    #                sf.append(sni)
    #        return States(sf)
        
    def both_combi_explosion(self, vars):
        """list(itertools.permutations([1,2,3]))
            [(1, 2, 3), (1, 3, 2), (2, 1, 3), 
            (2, 3, 1), (3, 1, 2), (3, 2, 1)]
            etc.. stick to two goals max...
        """
        states = []
        for s in self.states:
            sn = Goal.both(a,b)(s)
            for sni  in sn:
                states.append(sni)
        sf = []
        for s in states:
            sn = Goal.both(a,b)(s)
            for sni  in sn:
                sf.append(sni)
        return States(sf)
        
class State(object):
    
    def __init__(self, values=None):
        self.values = values or {}#OrderedDict()
        #self.vars = []
        
    def __repr__(self):
        st = ''
        for val in self.values:
            st+='  {}=>{}  \n'.format(val,self.values[val])
        return '{\n'+st+'}'
        
    def __call__(self, x):
        return self.value_of(x)
    
    def __add__(self, v):
        """s + (x,y,c)
        """
        return Goal.addo(v[0],v[1],v[2])(self)
    
    def __sub__(self, v):
        """s + (x,y,c)
        """
        return Goal.subo(v[0],v[1],v[2])(self)
    
    def bind(self, name):
        """bind a name to state with no assignment"""
        values = copy.copy(self.values)
        var = Variable(name)
        values[var] = None
        return (State(values), var)
    
    def assign(self, var, value):
        values = copy.copy(self.values)
        values[var] = value
        return State(values)
    
    def assign_bind(self,var,value):
        s,var = self.bind(var)
        return s.assign(var, value)
    
    def fresh(self, vars, func=None):
        values = copy.copy(self.values)
        for var in vars:
            var = Variable(var)
            values[var] = None
        if func is not None:
            return func(self)
        else:
            return (State(values),vars)
        
    def value_of(self, key):
        if key in self.values:
            value = self.values[key]
            if value is None:
                return key
            else:
                return self.value_of(value)
        elif isinstance(key, tuple):
            return (self.value_of(key[0]), 
                    self.value_of(key[1]))
        else:
            return key
            
    def unify(self, a, b):
        ivars = a,b
        a = self.value_of(a)
        b = self.value_of(b)
        if a == b:
            return self
        elif isinstance(a, Variable):
            return self.assign(a, b)
        elif isinstance(b, Variable):
            return self.assign(b, a)
        elif isinstance(a, ia) and isinstance(b, ia):
            values = copy.copy(self.values)
            values[ivars[0]] = a & b
            values[ivars[1]] = values[ivars[0]] & b
            #values[ivars[0]] = values[ivars[0]] & values[ivars[1]]
            return State(values)
        elif isinstance(a, ia) and isinstance(b, Variable):
            values = copy.copy(self.values)
            values[ivars[1]] = a#ia(a,a)
            return State(values)
        elif isinstance(b, ia) and isinstance(a, Variable):
            values = copy.copy(self.values)
            values[ivars[0]] = b#ia(b,b)
            return State(values)
        elif isinstance(a, tuple) and isinstance(b, tuple):
            unify_left = self.unify(a[0], b[0])
            if unify_left:
                unify_right = unify_left.unify(a[1], b[1])
                return unify_right
            else:
                return None
        else:
            return None


    def lessthanify(self, a, b, option=None):
        """a <= b
        for intervals (a,b)
        
        is intersection ncessary?
        """
        ivars = a,b
        a = self.value_of(a)
        b = self.value_of(b)
        if a == b:
            return self
        elif isinstance(a, Variable) and isinstance(b, ia):
            values = copy.copy(self.values)
            values[ivars[0]] = ia(min(0.,b.inf),b.sup)
            return State(values)
        elif isinstance(b, Variable) and isinstance(a, ia):
            values = copy.copy(self.values)
            values[ivars[1]] = ia(min(0.,a.inf) ,a.sup) #assumes + only!
            return State(values)
        elif isinstance(a, ia) and isinstance(b, ia):
            values = copy.copy(self.values)
            c = a&b
            if c.isempty and a<b:
                values[ivars[0]] = a
                values[ivars[1]] = b
            else:
                values[ivars[0]] = ia(min(a.inf,c.inf),min(c.sup,a.sup))
                values[ivars[1]] = ia(max(b.inf,c.inf),max(b.sup,c.sup))
            return State(values)
        else:
            return self
    
    def greaterthanify(self, a, b, option=None):
        """ a >= b
        for intervals (a,b)
        
        is intersection ncessary?
        """
        ivars = a,b
        a = self.value_of(a)
        b = self.value_of(b)
        if a == b:
            return self
        elif isinstance(a, Variable):
            return self
        elif isinstance(b, Variable):
            return self
        elif isinstance(a, ia) and isinstance(b, ia):
            values = copy.copy(self.values)
            c = a&b
            values[ivars[0]] = ia(max(a.inf,c.inf),max(a.sup,c.sup))
            values[ivars[1]] = ia(min(b.inf,c.inf),min(c.sup,b.sup))
            return State(values)

            

class EarlyGoalError(Exception):
    """ A Goal has been constructed prematurely"""
    
def _reify(o, s):
    return o 
    
def reify(e, s):
    """ Replace variables of expression with substitution

    >>> from logpy.unification import reify, var
    >>> x, y = var(), var()
    >>> e = (1, x, (3, y))
    >>> s = {x: 2, y: 4}
    >>> reify(e, s)
    (1, 2, (3, 4))

    >>> e = {1: x, 3: (y, 5)}
    >>> reify(e, s)
    {1: 2, 3: (4, 5)}
    
    Not yet used
    """
    if isinstance(e,Variable):
        return reify(s[e], s) if e in s else e
    return _reify(e, s)
    

    
def binop(op, revop=None):
    """ Transform binary operator into goal

    >>> from logpy.arith import binop
    >>> import operator
    >>> add = binop(operator.add, operator.sub)

    >>> from logpy import var, run
    >>> x = var('x')
    >>> next(add(1, 2, x)({}))
    {~x: 3}
    """

    def goal(x, y, z):
        if not isinstance(x, Variable) and not isinstance(y, Variable):
            return Goal.eq(op(x, y), z)
        if not isinstance(y, Variable) and not isinstance(z, Variable) and revop:
            return Goal.eq(x, revop(z, y))
        if not isinstance(x, Variable) and not isinstance(z, Variable) and revop:
            return Goal.eq(y, revop(z, x))
        else:
            return Goal.eq(None,None)#(z, revop(x,y))#NULL#State() #raise EarlyGoalError()#null #

    goal.__name__ = op.__name__
    return goal
    
add = binop(operator.add, operator.sub)
add.__doc__ = """ x + y == z """
mul = binop(operator.mul, operator.truediv)
mul.__doc__ = """ x * y == z """
mod = binop(operator.mod)
mod.__doc__ = """ x % y == z """

def null(s):
    return Goal.eq(None,None)(s) #[NULL]

def run_ck_relation(vars, state, goalfun):
    """not used
    """
    #valsi = [state(var) for var in vars]
    statei = copy.copy(state)
    
    #(was?) ternary op, now list comprehension
    chg_lst = [var for var in vars if not state(var) == statei(var)]

    return

def run_unary_relation(x,y, state, goal):#, option=None):
    x1 = state.value_of(x)
    y1 = state.value_of(y)
    #state = goal(x1,y)(state)[0]
    #state = goal(x,y1)(state)
    state = goal(x1,y1)(state) #equiv?why
    return state
    
def run_binary_relation(x,y,z, state, goal):#, option=None):
    x1 = state.value_of(x)
    y1 = state.value_of(y)
    z1 = state.value_of(z)
    state = goal(x1,y1,z)(state)[0]
    state = goal(x1,y,z1)(state)[0]
    state = goal(x,y1,z1)(state)[0]
    state = goal(x1,y1,z1)(state)
    return state
    
def run_br(x,y,z, state, goal):
    return run_binary_relation(x,y,z, state, goal)
    
run_binary_relation.__doc__ = """ x1 o x2 <=> x3 with state and goal"""
run_br.__doc__ = """ x1 o x2 <=> x3 with state and goal"""

class Goal(object):
    
    def __init__(self, func):
        self.func = func

    def pursue(self, state):
        return self.func(state)
        
    @staticmethod
    def null(a,b):
        def pursue(state):
            return [state]
        return
        
    @staticmethod
    def bind(names, func):
        def pursue(state):
            vars = []
            for name in names:
                (state, var) = state.bind(name)
                vars.append(var)
            #state.vars = vars
            goal = func(vars)
            return goal(state)
            
        return pursue
    
    @staticmethod
    def eq(a, b):
        #print 'equality ',a,'?:=',b
        def pursue(state):
            new_state = state.unify(a, b)
            if new_state:
                return [new_state]
            else:
                return []
    
        return pursue
        
    @staticmethod
    def lt_ia_(x, y):
        """ Actually implemented in lessthanunify
        x : x < y """
        a = x&y
        return ia(x.inf,min(a.sup,x.sup))
        
    @staticmethod
    def lto(x, y):
        """  Less than Constrant
        x : x < y """
        def pursue(state):
            new_state = state.lessthanify(x, y)
            if new_state:
                return [new_state]
            else:
                return []
    
        return pursue
        
    
    @staticmethod
    def gto(x, y):
        """  Less than Constrant
        x | x > y """
        def pursue(state):
            new_state = state.greaterthanify(x, y)
            if new_state:
                return [new_state]
            else:
                return []
    
        return pursue
    

        
    @staticmethod
    def gt(x, y):
        """ x > y """
        if not isinstance(x, Variable) and not isinstance(y, Variable):
            return Goal.eq(x > y, True)
        else:
            raise Goal.EarlyGoalError()
    
    @staticmethod
    def lt(x, y):
        """ x < y """
        if not isinstance(x, Variable) and not isinstance(y, Variable):
            return Goal.eq(x < y, True)
        else:
            raise Goal.EarlyGoalError()
    
    
        
        
    @staticmethod
    def either(a, b):
        def pursue(state):
            #print state
            return a(state) + b(state)
    
        return pursue
    
    @staticmethod
    def both(a, b):
        def pursue(state):
            states = []
            for a_state in a(state):
                states += b(a_state)
            #for b_state in b(state):
            #    states += a(b_state)
            return states
    
        return pursue
    
        
    @staticmethod        
    def append(x, y, z):
        def x_not_null(vars):
            [first, x_rest, z_rest] = vars
            return Goal.both(
                Goal.both(
                    Goal.eq(x, (first, x_rest)),
                    Goal.eq(z, (first, z_rest))
                ),
                Goal.append(x_rest, y, z_rest)
            )
    
        return Goal.either(
            Goal.both(
                Goal.eq(x, NULL),
                Goal.eq(y, z)
            ),
            Goal.bind(['first', 'x_rest', 'z_rest'], x_not_null)
        )
        
    @staticmethod
    def add(x,y,z):
        return add(x,y,z)
        
    @staticmethod
    def addo(x,y,z):
        def pursue(state):
            s = run_binary_relation(x,y,z,state,Goal.add)
            return s
        return pursue
        
    @staticmethod
    def sub(x, y, z):
        """ x - y == z """
        return add(y, z, x)
        
    @staticmethod
    def subo(x,y,z):
        def pursue(state):
            s = run_binary_relation(x,y,z,state,Goal.sub)
            return s
        return pursue
        
        
    @staticmethod
    def mul(x,y,z):
        return mul(x,y,z)
        
    @staticmethod
    def mulo(x,y,z):
        def pursue(state):
            s = run_binary_relation(x,y,z,state,Goal.mul)
            return s
        return pursue
    

    
    @staticmethod
    def mulo_le(x,y,z):
        """may not be advisable to do this way"""
        def pursue(state):
            s = run_binary_relation(x,y,z,state,Goal.mul)#,option='<=')
            
            return s
        return pursue
        
        
        
    @staticmethod
    def div(x, y, z):
        """ x - y == z """
        return mul(z, y, x)
        
    @staticmethod
    def divo(x,y,z):
        def pursue(state):
            s = run_binary_relation(x,y,z,state,Goal.div)
            return s
        return pursue
        
    
    @staticmethod        
    def add_ruby(x, y, z):
        def x_not_null(vars):
            [sx, sz] = vars
            print 'sx ',sx
            print 'sz ',sz
            return Goal.both(
                Goal.both(
                    Goal.eq(x, (sx, INC)),
                    Goal.eq(z, (sz, INC))
                    ),
                Goal.add(sx, y, sz)
            )
    
        return Goal.either(
            Goal.both(
                Goal.eq(x, ZERO),
                Goal.eq(y, z)
                ),
            Goal.bind(['sx', 'sz'], x_not_null)
        )
    
    @staticmethod        
    def add_old(x, y, z):
        def x_not_null(vars):
            [x_rest, z_rest] = vars
            return Goal.both(
                Goal.both(
                    Goal.eq(x, (INC, x_rest)),
                    Goal.eq(z, (INC, z_rest))
                    ),
                Goal.append(x_rest, y, z_rest)
            )
    
        return Goal.either(
            Goal.both(
                Goal.eq(x, ZERO),
                Goal.eq(y, z)
                ),
            Goal.bind(['x_rest', 'z_rest'], x_not_null)
        )
        
    
    @staticmethod
    def HalfAdder(A,B):
        #S = int((A and not B) or (not A and B))
        #C = int(A and B)
        return (S, C)
    
    @staticmethod
    def FullAdder(A,B,C):
        #AB = int((A and not B) or (not A and B))
        #S = int((C and not AB) or (not C and AB))
        #CAB = int(C and AB)
        #C1 = int((A and B) or (CAB))
        return (S,C1)
        
def from_list(list):
    if list:
        return (list[0], from_list(list[1:]))
    else:
        return NULL

def to_list(tuple):
    if tuple == NULL:
        return []
    else:
        return [tuple[0]] + to_list(tuple[1])


NULL = object()

#
#def with_vars(vars):
#            [|smaller_x, smaller_z] = vars
#            return Goal.both(
#                            Goal.both(
#                              Goal.equal(x, Pair.new(INC, smaller_x)),
#                              Goal.equal(z, Pair.new(INC, smaller_z))
#                            ),
#                            add(smaller_x, y, smaller_z)
#                          )
#        goal = Goal.bind(['x'],with_vars)
#
#def notadd(x, y, z):
#    def locg():
#        Goal.either(
#            Goal.both(
#                  Goal.equal(x, ZERO),
#                  Goal.equal(y, z)
#                ),
#            Goal.with_variables { |smaller_x, smaller_z|
#              Goal.both(
#                    Goal.both(
#                      Goal.equal(x, Pair.new(INC, smaller_x)),
#                      Goal.equal(z, Pair.new(INC, smaller_z))
#                    ),
#                    add(smaller_x, y, smaller_z)
#          )
#        }
#  )

def search(goallist, state):
    n = len(goallist)
    rlist = [False for el in goallist]
    goal_hash = {}
    for i, goal in enumerate(goallist):
        goal_hash[i] = (goal,{'done':False})
    todo=[]
    for i, goal in enumerate(goallist):
        try:
            state = goal(state)
            goal_hash[i][1]['done'] = True
            rlist[i] = True
        except:
            todo.append(i)

def to_peano(num):
    if num==0:
        return ZERO
    else:
        return (INC, to_peano(num-1))
        
def from_peano(peano):
    if peano==ZERO:
        return 0
    else:
        return from_peano(peano[1]) + 1      

        

ZERO = Variable(name='0')
INC = Variable(name=':+')



if __name__ == '__main__':
    world = State()
    
    t1  = False# True# 
    t2  = False# True# 
    t3  = False# True# 
    t4  = False# True# 
    t5  = False# True# 
    t6  = False# True# 
    t7  = False# True#  
    t8  = False# True# 
    t9  = False# True# 
    t10 = False# True# 
    i1  = False# True# 
    i2  = False# True# 
    t11 = False #True# False# 
    i3  = True# False# 
    
    
    if t1:
        
        def with_vars(vars):
            [x, y] = vars
            return Goal.eq(x, 5)
        
        goal = Goal.bind(['x', 'y'], func = with_vars)
        w1 = goal(world)
        print w1
        
    
    if t2:
        print 'either'
        def with_vars(vars):
            [x, y] = vars
            return Goal.either(Goal.eq(x, 5), 
                               Goal.eq(x, 6))
        
        goal = Goal.bind(['x', 'y'], with_vars)
        w2 = goal(world)
        print w2
        
    if t3:
        print 'both'
        def with_vars(vars):
            [x, y] = vars
            return Goal.both(Goal.eq(x, 'wigwam'), 
                             Goal.eq(y, 'wotsit'))
        
        goal = Goal.bind(['x', 'y'], with_vars)
        w3 = goal(world)
        print w3
        
    if t3:
        print 'either'
        def with_vars(vars):
            [x, y] = vars
            return Goal.either(Goal.eq(x, 'wigwam'), 
                             Goal.eq(y, 'wotsit'))
        
        goal = Goal.bind(['x', 'y'], with_vars)
        w3 = goal(world)
        print w3
        
    if t4:
        print 'Pairs'#works
        """
        a=(3,x)
        b=(y,Pair(5,y))
        self = State()
        
        x = Variable(name='x')
        y = Variable(name='y')
        goal = with_vars([x,y])
        w4 = goal(world)
        print w4
        
        
        
        def with_vars(vars):
            [x, y] = vars
            return Goal.eq(( from_list(['l','o']),x ), 
                           ( y,(5,y) )
                           )
        goal = Goal.bind(['x','y'],with_vars)
        w47 = goal(world)
        print w47
        """
        def with_vars(vars):
            [x, y] = vars
            return Goal.eq((3,x), 
                           (y,(5,y)))
        goal = Goal.bind(['x', 'y'], with_vars)
        
        
        goal = Goal.bind(['x','y'],with_vars)
        w4 = goal(world)
        print w4
        
        def with_vars(vars):
            [x, y] = vars
            return Goal.eq(
                            from_list([1,2,3,4]) , 
                            (1,(2,(x)) )
                           )
        goal = Goal.bind(['x', 'y'], with_vars)
        w4c = goal(world)
        print to_list(w4c[0].value_of(w4c[0].values.keys()[0]))
        
    if t5:
        print 'Lists'#works
        #print from_list([1,2,3,4])
        #print to_list(from_list([1,2,3,4]))
        def with_vars(vars):
            [x, y, z] = vars
            return Goal.eq(from_list([x,2,z]), 
                              from_list([1,y,3]))
        goal = Goal.bind(['x','y','z'],with_vars)
        w5 = goal(world)
        print w5
        
    if t6:
        print 'append'#works
        """
        x = Variable(name='x')
        goal = with_vars([x])
        g6 = goal(world)
        for key in g6[0].values.keys():
            try:
                print key, to_list(g6[0].value_of(key))
            except:
                pass
        print g6[0].value_of(x)
        """
        def with_vars(vars):
            [x] = vars
            return Goal.append(from_list(['h','e']), 
                               from_list(['l','l','o']),
                               x)
        goal = Goal.bind(['x'],with_vars)
        w6 = goal(world)
        #print to_list(w6[0].value_of(w6[0].values.keys()[5]))
        #print to_list(w6[0].value_of(w6[0].values.keys()[1]))
        for key in w6[0].values.keys():
            try:
                print key, to_list(w6[0].value_of(key))
            except:
                pass
            
    x = Variable()
    s = State({x:None})
    g = Goal.append(from_list(['h','e']), 
                               from_list(['l','l','o']),
                               x)
    s = g(s)
    ls = to_list(s[0].value_of(x))
    print 'x:',ls
    
    if t7:
        print 't7'
        print 'reverse append'#works
        """
        x = Variable(name='x')
        y = from_list(['l','o'])
        z = from_list(['h','e','l','l','o'])
        
        first = x
        x_rest = y
        z_rest = z
        
        x1=Goal.eq(x, (first, x_rest))
        x2=Goal.eq(z, (first, z_rest))
        
        def x_not_null(vars):
            [first, x_rest, z_rest] = vars
            return Goal.both(
                Goal.both(
                    Goal.eq(x, (first, x_rest)),
                    Goal.eq(z, (first, z_rest))
                ),
                Goal.append(y, x_rest, z_rest)
            )
        
        goal = Goal.either(Goal.both(
                        Goal.eq(x, NULL),
                        Goal.eq(y, z)
                    ),
                    Goal.bind(['first', 'x_rest', 'z_rest'], x_not_null)
                )
        
        goal=Goal.eq(x, NULL)
        pp = goal(world) #empty object
        
        goal = Goal.eq(y, z)
        a=y
        b=z
        new_state = world.unify(a, b)
        pp = goal(world)
        
        
        goal = Goal.bind(['first', 'x_rest', 'z_rest'], x_not_null)
        goal = Goal.bind(['x','y','z'], x_not_null)
        
        pp = goal(State())
            
        goal = with_vars([x])
        #or
        g7 = goal(world)
        
        
        for key in g7[0].values.keys():
            try:
                print key, to_list(g7[0].value_of(key))
            except:
                pass
        #or
        goal = Goal.append(x, 
                           from_list(['l','o']),
                           from_list(['h','e','l','l','o'])
                           )
        w7 = goal(world)
        """
        def with_vars(vars):
            [x] = vars
            return Goal.append(x, 
                               from_list(['l','o']),
                               from_list(['h','e','l','l','o'])
                               )
        goal = Goal.bind(['x'],with_vars)
        
        
        w7 = goal(world)
        #print w7
        the_key = -1
        for key in w7[0].values.keys():
            try:
                #print key, to_list(w7[0].value_of(key))
                if 'x' == key.name:
                    the_key=key
            except:
                pass
        print the_key,'<=>',to_list(w7[0].value_of(the_key))
        
        x = Variable()
        s = State({x:None})
        g = Goal.append(x, from_list(['l','o']),
                           from_list(['h','e','l','l','o'])
                           )
        s = g(s)
        print 'x:',to_list(s[0].value_of(x))
        
    if t8:
        print 'append all'# works
        """
        x = Variable(name='x')
        y = Variable(name='y')
        Goal.append(x, 
                   y,
                   from_list(['h','e','l','l','o'])
                   )
        """
        def with_vars(vars):
            [x,y] = vars
            return Goal.append(x, 
                               y,
                               from_list(['h','e','l','l','o'])
                               )
        goal = Goal.bind(['x','y'],with_vars)
        w8 = goal(world)
        nms = ['x','y']
        for st in w8:
            for key in st.values:
                if key.name in nms:
                    print key,'<=>',to_list(st.value_of(key)) 
        ### same as:
        x = Variable('x')
        y = Variable('y')
        nms = ['x','y']
        s = State({x:None,y:None})
        g = Goal.append(x, y,
                           from_list(['h','e','l','l','o'])
                           )
        s = g(s)
        for st in s:
            print to_list(st.value_of(x)), to_list(st.value_of(y))
        ### same as:
        s = State()
        names = ['x','y']
        func = with_vars
        vars = []
        for name in names:
            (s, var) = s.bind(name)
            vars.append(var)
        g = func(vars)#(s)
        s1 = g(s)
        nms = ['x','y']
        for st in w8:
            for key in st.values:
                if key.name in nms:
                    print key,'<=>',to_list(st.value_of(key))
    
    if t9:
        print 'numbers'# 
        print to_peano(3)
        print from_peano(to_peano(3))
        print from_peano(
                (INC,
                 (INC,
                  (INC,ZERO))))
                  
        five=to_peano(5)
        three=to_peano(3)
        def with_vars(vars):
            [x] = vars
            return Goal.add(   five,
                               three,
                               x)
        goal = Goal.bind(['x'],with_vars)
        w9 = goal(world)
        x=w9[0].values.keys()[0]
        from_peano(w9[0](x))
    if t10:
        
        def with_vars(vars):
            [x, y] = vars
            return add(0,x,5)
            #return Goal.add(Goal.eq(y,5),x,3)
        
        goal = Goal.bind(['x', 'y'], func = with_vars)
        w10 = goal(world)
        #print w10
        
        the_key = -1
        """
        for key in w10[0].values.keys():
            try:
                if 'x' == key.name:
                    the_key=key
            except:
                pass
        print w10[0].value_of(the_key)
        #"""
        print 'Addition'
        x = Variable('x')
        g=add(0,-5,x)
        s=State(values={x:None})
        r=g(s)
        print '0 - 5 = (x={})'.format(r[0].value_of(x))
        
        x = Variable('x')
        g=add(0,x,5)
        s=State(values={x:None})
        r=g(s)
        print '0 + (x={}) = 5'.format(r[0].value_of(x))
        
        x = Variable('x')
        sv = OrderedDict()
        sv[x]=None
        s=State(values=OrderedDict())
        g=add(x,0,-5)
        r=g(s)
        print '(x={}) + 0 = -5'.format(r[0].value_of(x))
        
        
        x = Variable('x')
        y = Variable('y')
        s=State(values={x:None,y:None})
        g = Goal.eq(y,-3)
        s = g(s)
        #g=add(x,y,5)
        const = 5
        r=run_binary_relation(x,y,const, s[0], Goal.sub)
        print '(x={}) - (y={}) = {}'.format(r[0].value_of(x), 
                                            r[0].value_of(y),
                                            const)
        
        
        x = Variable('x')
        sv = OrderedDict()
        sv[x]=None
        s=State(values=sv)
        g=add(x,0,-5)
        a=g(s)
        #r=run(sv, g)
        print '(x={}) + 0 = -5'.format(r[0].value_of(x)) 
    
    if i1:
        print '\ni1'
        print '\ntesting forward unification, y<=>ia(0.,1.)'
        x = Variable('x')
        y = Variable('y')
        s = State(values={x:None,y:None})
        g = Goal.eq(y,ia(0.,1.))
        s = g(s)
        print s
        print 'testing reverse unification, ia(0.,1.)<=>y'
        s = State(values={x:None,y:None})
        g = Goal.eq(ia(0.,1.),y)
        s = g(s)
        print s
        print 'testing compound addition'
        const = 5
        r=run_binary_relation(x,y,const, s[0], Goal.add)
        print '(x={}) + (y={}) = {}'.format(r[0].value_of(x), 
                                            r[0].value_of(y),
                                            const)
        print 'testing compound subtraction'
        const = 5
        r=run_binary_relation(x,y,const, s[0], Goal.sub)
        print '(x={}) - (y={}) = {}'.format(r[0].value_of(x), 
                                            r[0].value_of(y),
                                            const)
        print '\n\nswitch x and y\n'
        print 'testing compound addition'
        const = 5
        r=run_binary_relation(y,x,const, s[0], Goal.add)
        print '(y={}) + (x={}) = {}'.format(r[0].value_of(y), 
                                            r[0].value_of(x),
                                            const)
        print 'testing compound subtraction'
        const = 5
        r=run_binary_relation(y,x,const, s[0], Goal.sub)
        print '(y={}) - (x={}) = {}'.format(r[0].value_of(y), 
                                            r[0].value_of(x),
                                            const)
                                            
        print 'testing compound multiplication'
        x = Variable('x')
        y = Variable('y')
        s = State(values={x:None,y:None})
        #g = Goal.eq(y,3)
        g = Goal.eq(ia(1.,2.),y)
        s = g(s)
        const = 5
        print s
        r=run_binary_relation(x,y,const, s[0], Goal.mul)
        print '(x={}) * (y={}) = {}'.format(r[0].value_of(x), 
                                            r[0].value_of(y),
                                            const)
        print 'testing compound division'
        const = 5
        r=run_binary_relation(x,y,const, s[0], Goal.div)
        print '(x={}) / (y={}) = {}'.format(r[0].value_of(x), 
                                            r[0].value_of(y),
                                            const)
        print '\n\nswitch x and y\n'
        print 'testing compound multiplication'
        x = Variable('x')
        y = Variable('y')
        s = State(values={x:None,y:None})
        #g = Goal.eq(y,3)
        g = Goal.eq(ia(1.,2.),y)
        s = g(s)
        const = 5
        print s
        r=run_binary_relation(y,x,const, s[0], Goal.mul)
        print '(y={}) * (x={}) = {}'.format(r[0].value_of(y), 
                                            r[0].value_of(x),
                                            const)
        print 'testing compound division'
        const = 5
        r=run_binary_relation(y,x,const, s[0], Goal.div)
        print '(y={}) / (x={}) = {}'.format(r[0].value_of(y), 
                                            r[0].value_of(x),
                                            const)
        
        print 'testing equality of intervals'
        a = Variable('a')
        b = Variable('b')
        si = State()
        si = si.assign(a,None)
        si = si.assign(b,None)
        #si = State({a:None, b:None})
        si = Goal.eq(a, ia(0.,2.))(si)[0]
        si = Goal.eq(b, ia(1.,3.))(si)[0]
        print 'interval a: ',si.value_of(a)
        print 'interval b: ',si.value_of(b)
        sg = Goal.eq(a,b)(si)
        print 'unified a <?> b:',sg[0].value_of(a),'<=',sg[0].value_of(a)==sg[0].value_of(b),'=>',sg[0].value_of(b)

    if i2:
        print '\ni2'
        print '\ntesting forward unification, y<=>ia(0.,1.)'
        x = Variable('x')
        y = Variable('y')
        s = State(values={x:None,y:None})
        g = Goal.eq(y,ia(0.,1.))
        s = g(s)
        print s
        print 'testing reverse unification, ia(0.,1.)<=>y'
        s = State(values={x:None,y:None})
        g = Goal.eq(ia(0.,1.),y)
        s = g(s)
        print s
        print 'testing compound addition'
        const = 5
        r=run_binary_relation(x,y,const, s[0], Goal.add)
        print '(x={}) + (y={}) = {}'.format(r[0].value_of(x), 
                                            r[0].value_of(y),
                                            const)
        
        def with_vars(vars):
            [x, y, z] = vars
            return Goal.both( Goal.both(Goal.eq(x, .9), 
                                        Goal.eq(y, x)),
                             Goal.eq(y, z))
        goal = Goal.bind(['x', 'y','z'], with_vars)
        w3 = goal(world) # x<=>1 y<=>1 z<=>1
        
        ##scratch work 1
        #x = Variable('x')
        #y = Variable('y')
        #z = Variable('z')
        #s = State(values={x:None,
        #                  y:None,
        #                  z:None})
        #s=Goal.eq(y,3)(s)[0]
        #g=Goal.addo(x,y,4)
        #s=g(s)[0]
        
        ##scratch work 2
        #g1=Goal.eq(y,3) 
        #g2 = Goal.addo(x,y,4)
        #states = []
        #for st in g1(s):
        #    states += g2(st)
        
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        s = State(values={x:None,
                          y:None,
                          z:None})
        const = 4.
        
        
        g = Goal.both(  Goal.both(Goal.eq(y,3),
                                              Goal.addo(x,y,4)),
                        Goal.both(Goal.divo(49,z,7),
                                              Goal.mulo(z,y,21) ))
        s=g(s)
        
            
        
        
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        s = State(values={x:None,
                          y:None,
                          z:None})
        const = 4.
        #g=Goal.eq(y,3)
        #s=g(s)[0]
        g = Goal.either(Goal.addo(x,y,4),
                        Goal.eq(y,3)
                        )
        s=g(s)
        states = []
        for st in s:
            states.append(g(st))  #gets it right, but only in 1 state


        print'Hardo Problem'
        
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        s = State(values={x:None,
                          y:None,
                          z:None})
        const = 4.
        
        g = Goal.either(Goal.either(Goal.eq(x,7.),
                                    Goal.eq(5.,x)),
                        Goal.either(Goal.eq(x,7.),
                                    Goal.eq(5.,x))
                        )
        s = g(s)
        
        
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        s = State(values={x:None,
                          y:None,
                          z:None})
        const = 4.
        
        g = Goal.either(Goal.either(Goal.eq(y,3.),
                                    Goal.addo(x,y,4)),
                        Goal.either(Goal.addo(x,y,4),
                                    Goal.eq(y,3.))
                        )
        s = g(s)
        
        
        
        
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        s = State(values={x:None,
                          y:None,
                          z:None})
        const = 4.
        
        g = Goal.both(Goal.both(Goal.eq(y,3.),
                                    Goal.addo(x,y,4)),
                        Goal.both(Goal.addo(x,y,4),
                                    Goal.eq(y,3.))
                        )
        s = g(s)
        
        
        
        #the above is not needed?
        
        #just do this?
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        s = State(values={x:None,
                          y:None,
                          z:None})
        const = 4.
        
        g = Goal.both(Goal.addo(x,y,4),
                      Goal.eq(y,3.))
        s = g(s)
        #nope, the order matters, making for ugly code down the line
        
        #instead do this:
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        s = State(values={x:None,
                          y:None,
                          z:None})
        S = States(s)
        S1 = S.both(Goal.addo(x,y,4),
                    Goal.eq(y,3.))
        
        #or this:
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        s = State(values={x:None,
                          y:None,
                          z:None})
        S = States(s)
        S1 = (S.both(Goal.addo(x,y,4),
                    Goal.eq(y,3.)) + (z,y,x))
        
        #or what about this:
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        s = State(values={x:None,
                          y:None,
                          z:None})
        S = States(s)
        S1 = (S.both(Goal.addo(z,y,x),
                    Goal.eq(y,3.)) + (x,y,4))
        #ah, darn it, the last function cannot push info
        # back up to the first two
                    
        
        #or what about this:
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        s = State(values={x:None,
                          y:None,
                          z:None})
        S = States(s)
        S1 = (S.both( S + (x,y,4),
                       S.both(Goal.addo(z,y,x),
                            Goal.eq(y,3.)) ))  #Fails.
        # back up to the first two
                    
        
        #or what about this:
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        s = State(values={x:None,
                          y:None,
                          z:None})
        S = States(s)
        S1 = (S.both(S.both(Goal.addo(z,y,x),
                            Goal.eq(y,3.)),
                     (S + (x,y,4)))
                     ) #Fails
                    
        
        #or what about this:
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        s = State(values={x:None,
                          y:None,
                          z:None})
        S = States(s)
        S1 = (S.both(S.both(Goal.addo(z,y,x),
                            Goal.eq(y,3.)),
                     (Goal.addo(x,y,4)))
                     ) #Fails.  Cannot S.both(States, Goal)
                     
        #or what about this:
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        s = State(values={x:None,
                          y:None,
                          z:None})
        S = States(s)
        S1 = S.both(Goal.both(Goal.addo(z,y,x),
                              Goal.addo(x,y,4)),
                    Goal.both(Goal.addo(z,y,x),Goal.eq(y,3.)))#Success!
                     
        #or what about this:
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        s = State(values={x:None,
                          y:None,
                          z:None})
        S = States(s)
        S1 = S.both(Goal.both(Goal.addo(z,y,x),
                              Goal.addo(x,y,4)),
                    Goal.eq(y,3.))#Unify's ony x and y
                     
        #or what about this:
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        s = State(values={x:None,
                          y:None,
                          z:None})
        S = States(s)
        S1 = S.both(Goal.both(Goal.addo(z,y,x),
                              Goal.eq(y,3.)),
                    Goal.both(Goal.addo(z,y,x),Goal.addo(x,y,4)))#Success!
                    
        
        #or what about this:
#        x = Variable('x')
#        y = Variable('y')
#        z = Variable('z')
#        s = State(values={x:None,
#                          y:None,
#                          z:None})
#        S = States(s)
#        S1 = (S.conj(S.both(Goal.addo(z,y,x),
#                            Goal.eq(y,3.)),
#                     (Goal.addo(x,y,4)))
#                     ) #Fails.  Cannot S.both(States, Goal)
        
        
                    
                    
                    
                      
        
        
    #        try:
    #            g = Goal.either(Goal.either(Goal.eq(y,3),
    #                                      Goal.addo(x,y,4)),
    #                            Goal.either(Goal.addo(x,y,4),
    #                                      Goal.eq(y,3)))
    #            s = g(s)
    #        except:
    #            g = Goal.both(Goal.eq(y,3),
    #                          Goal.addo(x,y,4))
    #            s = g(s)
                                
    #        g = Goal.either( Goal.addo(x,y,4),
    #                            Goal.eq(y,3) )
    #        s = g(s)
            
    #        g2 = Goal.both(Goal.divo(49,z,7),Goal.mulo(z,y,21) )
    #        states=[]
    #        for els in s:
    #            try:
    #                states.append(g2(els))
    #            except:
    #                pass
                                                   
    #        g = Goal.both(  Goal.either( Goal.both(Goal.addo(x,y,4),
    #                                               Goal.eq(y,3)),
    #                                     Goal.both(Goal.eq(y,3),
    #                                              Goal.addo(x,y,4))), 
    #                        Goal.both(Goal.divo(49,z,7),
    #                                              Goal.mulo(z,y,21) ))
    #        s=g(s)
    if t11:
        x=Variable('x')
        s=State({x:None})
        g = Goal.addo(x,3,8)
        s = g(s)
        
        five=to_peano(5)
        three=to_peano(3)
        eight = to_peano(8)
        def with_vars(vars):
            [x,y] = vars
            return Goal.add_ruby(x,y,eight)
            
        goal = Goal.bind(['x','y'],with_vars)
        s = goal(world)
        
        
        """recursion limit
        x = Variable('x')
        y = Variable('y')
        nms = ['x','y']
        s = State({x:None,y:None})
        g = Goal.add_ruby(x,y,eight)
        s = g(s)
        for st in s:
            print from_peano(st.value_of(x)), from_peano(st.value_of(y))
            print (st.value_of(x)), (st.value_of(y))
        #"""
#        x = Variable('x')
#        y = Variable('y')
#        nms = ['x','y']
#        s = State({x:None,y:None})
#        g = Goal.FullAdder(x,y,eight)
#        s = g(s)
#        for st in s:
#            print from_peano(st.value_of(x)), from_peano(st.value_of(y))
#            print (st.value_of(x)), (st.value_of(y))
        
        
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        
        const = 4.
        
        goallist = [Goal.eq(x,7.),Goal.addo(x,y,4)]
        try:
            s = State(values={x:None,
                              y:None,
                              z:None})
            s = Goal.addo(x,y,4)(s)
            ns=[]
            for st in s:
                ns.append(Goal.eq(x,7.)(st))
        except:
            print 'HERE'
            s = State(values={x:None,
                              y:None,
                              z:None})
            s = Goal.eq(x,7.)(s)
            ns=[]
            for st in s:
                ns.append(Goal.addo(x,y,4)(st))
                
#        try:
#            s = State(values={x:None,
#                              y:None,
#                              z:None})
#            s = Goal.eq(x,7.)(s)
#            ns=[]
#            for st in s:
#                ns.append(Goal.addo(x,y,4)(st))
#        except:
#            s = State(values={x:None,
#                              y:None,
#                              z:None})
#            s = Goal.eq(x,7.)(s)
#            ns=[]
#            for st in s:
#                ns.append(Goal.addo(x,y,4)(st))
                
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        
        const = 4.
        s = State(values={x:None,
                              y:None,
                              z:None})
        st = Goal.eq(x,y)(s)[0]
        sf=Goal.eq(x,7.)(st)
    
    if i3:
        print '\ni3'
        print '\ntesting forward unification, y<=>ia(0.,1.)'
        x = Variable('x')
        y = Variable('y')
        c = Variable('c')
        s = State(values={x:None,y:None,c:None})
        g = Goal.eq(y,ia(0.5,1.))
        s = g(s)[0]
        g = Goal.eq(x,ia(1.0,1.5))
        s = g(s)[0]
        g = Goal.eq(c,ia(0.0,1.9))
        s = g(s)[0]
        g = Goal.addo(x,y,c)
        s = g(s)[0]
        
        a = ia(1.0,1.5)
        b = ia(0.5,1.0)
        c1 = a+b #should trim contraint when done in s!
        print '\nnaive c = {}'.format(c1)
        print 'minimal relational c = {}\n'.format(s.value_of(c))
        
        q = Variable('q')
        w = Variable('w')
        e = Variable('e')
        st = State(values={q:None,w:None,e:None})
        g = Goal.eq(q, ia(.1,1.5))
        st = g(st)[0]
        g = Goal.eq(w, ia(1.0,1.9))
        st = g(st)[0]
        g = Goal.eq(q,w) #unifies via intersection
        st = g(st)[0]
        
        
        
        
        x = Variable('x')
        y = Variable('y')
        c = Variable('c')
        d = Variable('d')
        s = State(values={x:None,y:None,c:None,d:None})
        g = Goal.eq(y,ia(0.5,1.))
        s = g(s)[0]
        g = Goal.eq(x,ia(1.0,1.5))
        s = g(s)[0]
        g = Goal.eq(c,ia(0.0,1.9))
        s = g(s)[0]
        s1 = (s - (y,x,c))#g(s)[0]
        a = ia(1.0,1.5)
        b = ia(0.5,1.0)
        c1 = b-a #should trim contraint when done in s!
        #c1 & ia(0.,1.9) = (0.,0.)
        print '\nnaive c = {}'.format(c1)
        print 'minimal relational c = {}\n'.format(s1[0].value_of(c))
        
        s2 = (s + (x,y,c))#g(s)[0]
        a = ia(1.0,1.5)
        b = ia(0.5,1.0)
        c1 = a+b #should trim contraint when done in s!
        print '\nnaive c = {}'.format(c1)
        print 'minimal relational c = {}\n'.format(s2[0].value_of(c))
        
        q = Variable('q')
        w = Variable('w')
        e = Variable('e')
        st = State(values={q:None,w:None,e:None})
        g = Goal.eq(q, ia(.1,1.5))
        st = g(st)[0]
        g = Goal.eq(w, ia(1.0,1.9))
        st = g(st)[0]
        g = Goal.eq(q,w) #unifies via intersection
        st = g(st)[0]
        
        s3 = ((s + (y,x,c))[0] + (x,y,d))#g(s)[0]
        a = ia(1.0,1.5)
        b = ia(0.5,1.0)
        c1 = b+a #should trim contraint when done in s!
        #c1 & ia(0.,1.9) = (0.,0.)
        print '\nnaive c = {}'.format(c1)
        print 'minimal relational c = {}\n'.format(s3[0].value_of(c))
        print 'minimal relational d = {}\n'.format(s3[0].value_of(d))
        
        
        
        st = States(s)
        st1 = ( (st + (x,y,c) ) + (x,y,d) )
        
        st2 = (st + (x,y,c))
        st2 = (st2 + (x,y,d))