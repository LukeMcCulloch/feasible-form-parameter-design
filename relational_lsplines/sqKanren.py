# -*- coding: utf-8 -*-
"""
Extended uKanren
sq :: working multi-state interval constraint logic

uKanren with extended interval arithmetic

Created on Sun Jul 17 11:17:25 2016

@author: lukemcculloch? well, for parts of this, yes I am.

-Genisis of the idea, and some basics of miniKanren,
  including computing with multivalued returns the miniKanren way:
      https://codon.com/hello-declarative-world
-Some nuts and bolts for hooking up relational operators 
   in Python:
      https://github.com/logpy/logpy

-Parts I created:
       *rules graph processor and associated classes
          these enable the relational programming to 
          act like python on the front end
          (genesis of the idea goes to Dr. Kearfott)
          -as to what the 'parser/accumulator' should accumulate? 
          - Bruce Christensen's reverse mode AD paper. _cite_here__
       *use of intervals in the minikanren paradigm
          genesis of the idea was from fwd bckwd contractors
          a paper by ____cite_here__.  
       


uKanren
    -minus addo, mulo, .. as at least addo kept recursing forever
    -exteded with binary mathematical relations and inerval unification
    
    
idea: 
    
    build a parse tree of the higher order constraint
    and then run it fwd - bckwd according to HC4 in 
    Benhamou's 'revising hull and box consistency'
    
    parse tree:
        http://interactivepython.org/runestone/static/pythonds/Trees/ParseTree.html

Notes:
    -There is currently a layer of seperation between parser 
    variables, which are the PStates class
    and
    miniKanren variables, which are of the variable class
    Why not make them the same?
    
        *it would be cool to see the graph of latest rule applied
        whenever you wanted.
        MORE IMPORTANTLY
        *the user could look up the value of her variables
        directly in the databases using her PStates objects.
        *an alternative might be to overload a getter in
        the classes which use PStates directly - so as 
        to have those classes go directly to the env(name) value
        from the PState variables.
        AND
        *It would be 'easier' to use the system (directly, again)
        in the functional way - where we 
        pass the env through the rules filter style.
        BUT
        *maybe its better to do this with 
        a bit of quasi-fancy footwork in the stuff you
        have now.  -Separation of concerns keeps the PStates
        and Kanren classes simpler.  People are going 
        to get lost looking at rules graph processor etc.
        
        
    TODO:
        -Write a Gauss Siedel solver in the States langague
        -Give it some facility to self vary it's AC_revision
        from 0 iterations to (?open ended?) when the solver gets stuck
        
    Clean house:
        need tail recurions?  Trampoline style?
        https://eli.thegreenplace.net/2017/on-recursion-continuations-and-trampolines/
        
        see also monad -pattern encapsulation cited over in 
            hull_from_simple_design_space.
"""
#
#
#
from collections import OrderedDict, Iterable
import itertools  #.permutations to have func with combi_explosivity
from inspect import isfunction
import operator as op#python native opps
#from math import log as mlog#needed for log base n -> just use the log rule dummy
import numpy as np #just using this for logs
import copy
from interval_arithmetic import ia as dmyia
from extended_interval_arithmetic import ia
from functools import wraps #decorator inside a class
#
from functools import partial
from kanren_utils import *
#
#easy decorators with args:
#from funcy import decorator

#
#import math
#
#
#
#Symbol = str          # A Scheme Symbol is implemented as a Python str
#List   = list         # A Scheme List is implemented as a Python list
#Number = (int, float) # A Scheme Number is implemented as a Python int or float
#
#def tokenize(chars):
#    "Convert a string of characters into a list of tokens."
#    return chars.replace('(', ' ( ').replace(')', ' ) ').split()
#
#
#def read_from_tokens(tokens):
#    "Read an expression from a sequence of tokens."
#    if len(tokens) == 0:
#        raise SyntaxError('unexpected EOF while reading')
#    token = tokens.pop(0)
#    if '(' == token:
#        L = []
#        while tokens[0] != ')':
#            L.append(read_from_tokens(tokens))
#        tokens.pop(0) # pop off ')'
#        return L
#    elif ')' == token:
#        raise SyntaxError('unexpected )')
#    else:
#        return atom(token)
#
#def atom(token):
#    "Numbers become numbers; every other token is a symbol."
#    try: return int(token)
#    except ValueError:
#        try: return float(token)
#        except ValueError:
#            return Symbol(token)
#        
#def parse(program):
#    "Read a Scheme expression from a string."
#    return read_from_tokens(tokenize(program))


#class Env(dict):
#    "An environment: a dict of {'var':val} pairs, with an outer Env."
#    def __init__(self, parms=(), args=(), outer=None):
#        self.update(zip(parms, args))
#        self.outer = outer
#    def find(self, var):
#        "Find the innermost Env where var appears."
#        return self if (var in self) else self.outer.find(var) 


#global_env = standard_env()
#
#
#

#def rule_generator(rule, state):
"""
rule = str version of a rule 

where substitutions.... happen?
where dummies are ... turned set equal to LHS
and then incorporated in the next level of relations...?

say rule:=
    x+y+z = A

which has to be delt with like:
    (x+y) = d1
    d1 + z = A

so drop in syntatical form:  rule = (x+y = _0_) + z = A

parse (x + y = _0_)

rule_generator(rule = 'x+y=_0_')

set nth+1 dummy state.d(N+1) <= _0_

state = (state + (x,y,_0_))

or return def r(): _add_ (x,y,_0_)

def r(): _0_ + Z = A


Edit, DONE, See Rules GraphProcessor
"""
   # return 






def exp(a,b):
    meta = 'lambda x: {}*{}*x**2'.format(a,b)
    return meta
meta = exp(2,3)
func = eval(meta)
#print func(2.) #the nefarious 24. has been found.

#class BinaryOp(object):
#    def __init__(self,op,a,b):
#        self.op = op
#        self.a  = a
#        self.b  = b
#
#class TernaryRule(object):
#    def __init__(self,op,a,b,c):
#        self.op = op
#        self.a  = a
#        self.b  = b
#        self.c  = c
#        self.r1 = BinaryRule(self.op, self.a,self.b)
#
#
#class Rule(object):
#    def __init__(self, args):
#        self.args = args
#        return      
    

        
class DList(list):
    """NOT USED:x
        Overload List to probide a custom accessor
        in Class Hull
    """
    def getpoint(self, pt):
        return self[0].getpoint(pt)
        
        
class Variable(object):
    """sqKanren logic variable class
    These are the objects which will be combined in relations
    -the geometric logic programming rules used
    to enforce feasible form parametric design.
    
    When passed into a sqKanren environment,
    they map to other data, typically 
    either interval values
    or simply the name of the variable itself,
    if data is lacking.
    
    Notes:
    --------------
        Currently undefined variables map to themselves.  
        And undefined strings map to themselves.
        Perhaps this is ....
        ok, but not terribly _interesting_ as long as one is
        aware of it.... so it seems.
        
    Dev:
    --------------
        had I know what PStates was going to be
        maybe I would have built that in here.
        But maybe it is good to have more separation of
        concerns.  I don't know.
    """
    def __init__(self, name='', nbrs=None, 
                 arcs=None, dummy=False,
                 flist=None):
        #print 'making new var : ',name
        self.name = name
        self.nbrs = nbrs
        self.dummy = dummy
        if name == '_':
            self.dummy = True
        #micro scale (ternary unifiable) rules:
        if arcs is None:
            self.arcs = []
        else:
            self.arcs = arcs
        
        #macro scale rules:
        if flist is None:
            self.flist = []
        else:
            self.flist = flist

    def __repr__(self):
        return str(self.name)

    def inspect(self):
        print self.name
        return
    
class _(object):
    pass
_0_ = _()

class Node(object):
    def __init__(self,State=None):
        self.values = State
        self.children = None
        pass
    
    #@children.setter
    #def children(self, children):
    #    pass

    @property
    def values(self):
        return self._values
    @values.setter
    def values(self, values):
        pass
#
#def bind(self, name):
#    """bind a name to state with no assignment"""
#    values = copy.copy(self.values)
#    var = Variable(name)
#    values[var] = None
#    return (State(values), var)
#
#def assign(self, var, value):
#    values = copy.copy(self.values)
#    values[var] = value
#    return State(values)
#
#def assign_bind(self,var,value):
#    s,var = self.bind(var)
#    return s.assign(var, value)
#
#def fresh(self, vars, func=None):
#    values = copy.copy(self.values)
#    for var in vars:
#        var = Variable(var)
#        values[var] = None
#    if func is not None:
#        #TLM March 2017:
#        nself = copy.copy(self)
#        nself.values = values
#        return func(nself)
#        #return func(self)
#    else:
#        return (State(values),vars)

#def with_vars(vars):
#    [x, y] = vars
#    return Goal.eq(x, 5)

#goal = Goal.bind(['x', 'y'], func = with_vars)
#w1 = goal(world)
#print w1

#def procedure(parms, rule, env):
#    "A user-defined Scheme procedure."
#    newstate = Goal.fresh(parms,env)
#    def with_vars(parms):
#        return
#    Goal.bind(['x', 'y'], func = with_vars)
#    #
#    return lambda *args: eval(rule, State(parms, args, env))

#
#******************************************************************
#
#playing around with Norvig's lispy:
isa = isinstance
Symbol = str          # A Scheme Symbol is implemented as a Python str
List   = list         # A Scheme List is implemented as a Python list
Number = (int, float, ia, dmyia) # A Scheme Number is implemented as a Python int or float
# these definitions are used.
# nothing else of lispy is used.
#
#******************************************************************
#
def tokenize(chars):
    "Convert a string of characters into a list of tokens."
    return chars.replace('(', ' ( ').replace(')', ' ) ').split()


def parse(program):
    "Read a Scheme expression from a string."
    return read_from_tokens(tokenize(program))

def read_from_tokens(tokens):
    "Read an expression from a sequence of tokens."
    if len(tokens) == 0:
        raise SyntaxError('unexpected EOF while reading')
    token = tokens.pop(0)
    if '(' == token:
        L = []
        while tokens[0] != ')':
            L.append(read_from_tokens(tokens))
        tokens.pop(0) # pop off ')'
        return L
    elif ')' == token:
        raise SyntaxError('unexpected )')
    else:
        return atom(token)

def atom(token):
    "Numbers become numbers; every other token is a symbol."
    try: return int(token)
    except ValueError:
        try: return float(token)
        except ValueError:
            return Symbol(token)

#
#******************************************************************
#
#
#def recurse(cls,el):
#    """el := (x,'==',[(y,'==',[c,ia(0.3,2.2)]),ia(1.,2.)])"""
#    if not isa(el,(Variable,ia,Number)):
#        return cls.eval(el)
#    else:
#        return cls, el  
#
#@decorator
#def process_args(func):
#    """ """
#    v = func._args
#    if isa(v,list):
#            for i,el in enumerate(v):
#                cls, v[i] = recurse(el)
#    return func(v)



##
##*****************************************************
##
# args not arg1, arg2 
# because binary operations
# can always be ordered this way...
class PStates(object):
    """
    
    Function
    --------------
        syntax aware class
        which uses Python itself to construct computational graphs
        of expressions input by user.
        It is like a first stage of a reverse mode AD tool,
        but for other things than differentiation.
        
    To add a command:
    --------------
        
        1.) add a representation of it to the PStates class 
        so that it may be recognized when parsed
        
        2.)  Optional:  add a representation of it to 
        the States.eval method
        what does that do:  just syntatic sugar for use outside the compiler
        
        3.) Add the operator to the Goal.eval method
        this dispatches the rule to the appropriate goal


    Notes:
    --------------
        
        Once a rule is compiled, it would be best to 
        sanitize the graphs built up in the 
        PStates themselves.
        
        If a PStates var holding a graph is 
        used at the base of another rule, 
        then the old graph will be tied into the new one.
        
        This will not result in any incorrect behavior,
        but possibly inefficient behavior as this interdependency
        is handled in general with AC revise and the 
        (construct then compile) process.
        
    TODO:  wipe the args tree when a rule has been compiled 
    to the rule graph processor
    
    
    Credit:
    ---------------
    Baker Kearfott's suggestion 
    (on the day of my PhD {prospectus:'design review'}) 
    to use operator overloading to
    build up a reverse mode AD tool was instrumental in the design
    of this class.

    For the specifics of how to do it,    
    Bruce Christianson's paper: 'automatic hessians by reverse accumulation'
    helped me decide how to structure the computational graph building
    via returned/stored stuff as you see in PStates.
    
    Some nice tutorials on reverse mode AD helped me more quickly
    check 'this kind of thinking worked for that'.  
    Building confidence to move to other use cases.  ;)  
        e.g. 
        https://justindomke.wordpress.com/2009/03/24/a-simple-explanation-of-reverse-mode-automatic-differentiation/
        https://rufflewind.com/2016-12-30/reverse-mode-automatic-differentiation
        
        better?:
            https://stats.stackexchange.com/questions/224140/step-by-step-example-of-reverse-mode-automatic-differentiation
            http://www.columbia.edu/~ahd2125/post/2015/12/5/
    
    As for the 'rules compiler' - I winged it after reading a good bit of
    Sussman and Abelson's SICP and Norvig's AIP
    and after having developed a 'parser helper' 
                                    in the Goal class below.
    Especially discussions with Will Byrd 
                                    were the inspiration for 
    'taking it further'.
    
    
    
    
    TRIVIAL EXAMPLES
        --------------
        demonstrating what construct and compile do:
        (normally this is hidden away in the 
        RulesGraphProcessor function 'add_one_rule',
        where it is used as the pre-compile stage.
        so we use it explicity in these examples)
        
        Finally, note that RulesGraphProcessor fuctnion
            add_graph_scale_node(dflist,vars_,rulename)
            
            is the missing ingredient in the examples below
            which updates each logic variable (objects of class Variable)
            with the rules in which it is used.
            This is the final peice of the puzzle for efficent AC revision
        --------------
        
        Example 1: 
            
            (from test_RGP_add_one_rule.py)
            #ADD A RULE BY HAND, step one, construct:
            
            CBulbMid = lp.PStates(name = 'CBulbMid') 
            rgp = RulesGraphProcessor()
            
            #rgp.add_one_rule(CBulbMid,'CBulbMid')
            
            rule_node = CBulbMid   #note a node is a PStates variable
                                   #     simple computations are overloaded
                                   #     and the class builds up its
                                   #     own 'computation tree' in response
                                   #     to usage.
                                   #     -almost exactly like a reverse mode 
                                   #     AutoDiff tool.
                                   #     -but all that happens automatically
                                   #     in PStates itself, before construct
                                   #     would ever be called.
                                   #     -nevertheless one can succesfully
                                   #     call construct on a 'bare PStates node'
                                   #     and see what happens.
                                   #     That is what this example shows:
                                            
            
            
            st = rgp.env
                
            st, vars_ = rule_node.construct(rule_node, st, {})  # CBulbMid alias 'rule_node'
                                                                # is just to impress
                                                                # upon the user that 
                                                                # CBulbMid is a node in a 
                                                                # computational graph
                                                                
            
            vars_                               # vars_  is a map from PStates.names to Variables.
            >>> {'CBulbMid': CBulbMid}          # That is to say,
                                                # from node.names of the  
                                                # computational graph nodes 
                                                # involved in the rule
                                                # to logic Variables,
                                                # aka Variables class objects
                                                # which are the variables of the States environment
                                                # (RulesGraphProcessor.env - States object)
                                                # RulesGraphProcessor == rgp in this code snippet
                                                
                                                
            st                                  # st is an instance of the States class of this file
            >>>                                 # it is also known as the RGP.env (environment)
            States([{                           # the construct function 
              CBulbMid:CBulbMid                 # adds the Variables involved in the rule
            }])                                 # to the environment, as seen here.
    
    
        #aside: 
        #  -engineering parlance would call States the design space, which may be disjoint
        #  -functional programmers building things like lambda calculus inturpreters
        #  call this sort of thing the environment.
            
    
    
        
        Example 2:
        (assume the above does not exist, but it doesn't really matter)
        
        
            A_mid = lp.PStates(name='A_mid')
            A_mid = A_mid == ia(10.,20.) 
            
            >>> print A_mid
            'A_mid'==
               'ia(10.0,20.0)'
               
            
            #ADD A RULE BY HAND, step one, construct:
            
            rgp = RulesGraphProcessor()
            
            st = rgp.env #must use a consistent environment!
                         #this is one reason to automate.
                
            st, vars_ = A_mid.construct(A_mid, st, {})
            
            >>> vars_
            {'A_mid': A_mid, 'ia(10.0,20.0)': ia(10.0, 20.0)}
            
            >>>st 
            States([{ 
              A_mid:A_mid 
            }])
            
    
    
            st, vars_ = rule_node.construct(A_mid, st, {})
            >>> vars_
            {'A_mid': A_mid, 'ia(10.0,20.0)': ia(10.0, 20.0)}
            
            #jumping ahead, by passing the node A_mid into the compiler,
            together with the mapped rules, the following happens:
            
            dflist = rule_node.compile(A_mid,vars_)
            >>> dflist
            {0: <function sqKanren.pursue>} # a function which assigns
                                            # A_mid the value of ia(10.,20.)
                                            
            Though pointless here, since it is a simple fixed assignment,
            this example shows that relational expressions are 
            converted into rules.  
            These rules are actually __continuations__ which,
            when passed an environment, i.e. a States object,
            return a narrowed States object.
            Thus we are constructing lists of functions to execute 
            in any sequence.
            That will always act to narrow the design space.
            
            tl,dr:
                This seemingly baroque construction
                is going to help us compile to rules lists.
                Rules lists can be run over and over in any order.
                Rules lists propogate information in all connected directions.
                That's what we need.
                Thus we don't care about running a rule right now.
            
        
        Example 3:
        --------------
        
        
        #initial code:
        
            A_mid = lp.PStates(name='A_mid')
            BulbBeam    = lp.PStates(name = 'BulbBeam')    
            BulbDepth   = lp.PStates(name = 'BulbDepth')
            CBulbMid = lp.PStates(name = 'CBulbMid')
        
        
        
            BulbBeam = BulbBeam == ia(5.,10.)
            BulbDepth = BulbDepth == ia(-10.1,10.)
            CBulbMid = CBulbMid == ia(0.,1.)
            
            rgp = RulesGraphProcessor()
            st = rgp.env
        
            
            st, vars_ = BulbBeam.construct(BulbBeam, st, {})
            
            
        # what does that give us?:
            
        >>> print BulbBeam
        'BulbBeam'==
           'ia(5.0,10.0)' 
        
        >>> st
        >>> 
        States([{
          BulbBeam:BulbBeam  
        }])
        
        >>> vars_
        >>> {'BulbBeam': BulbBeam, 'ia(5.0,10.0)': ia(5.0, 10.0)}
        
        
        #back to the code:
        
        
        
            st, vars_ = BulbBeam.construct(BulbDepth, st, 
                                                existing=vars_)
            
            

        # what does that give us?:
        
        
        >>> st
        >>> 
        States([{
          BulbDepth:BulbDepth  
          BulbBeam:BulbBeam  
        }])
        
        >>> vars_
        >>>  
        {'BulbBeam': BulbBeam,
         'BulbDepth': BulbDepth,
         'ia(-10.1,10.0)': ia(-10.1, 10.0),
         'ia(5.0,10.0)': ia(5.0, 10.0)}
        
        #back to the code:
        
            st, vars_ = BulbBeam.construct(CBulbMid, st, 
                                                existing=vars_)
            
            
        # what does that give us?:
        
        >>> st
        >>> 
        States([{
          BulbDepth:BulbDepth  
          BulbBeam:BulbBeam  
          CBulbMid:CBulbMid  
        }])
        
        >>> vars_
        >>> 
        {'BulbBeam': BulbBeam,
         'BulbDepth': BulbDepth,
         'CBulbMid': CBulbMid,
         'ia(-10.1,10.0)': ia(-10.1, 10.0),
         'ia(0.0,1.0)': ia(0.0, 1.0),
         'ia(5.0,10.0)': ia(5.0, 10.0)}
        
        #back to the code:
        
        
            # now compile these vars with values and the environment 
            #to get a dictionary of logical, relational, continuations
            # which together make up the entire rule.
            # for these first simple rules, the rule dict
            # only contains one relation.
        
            # side notes:
            # self doesn't do anything in the compiler
            # (s)could make this thing a classmethod.
        
            
            code_list = []
            dflist = BulbBeam.compile(BulbBeam,vars_) 
            code_list.append(dflist)
            
        >>> dflist
         {0: <function sqKanren.pursue>}
            
            dflist = BulbBeam.compile(BulbDepth,vars_)
            code_list.append(dflist)
            dflist = BulbBeam.compile(CBulbMid,vars_)
            code_list.append(dflist)
        
        
                    
            #critical:  compile CBulbMid (as we just did above) 
            #before this next rule:
            
            CBulbMid = CBulbMid == A_mid/(BulbBeam*BulbDepth)
            
            
        >>> print CBulbMid
         'CBulbMid'==
            'v3'/
               'A_mid' 
                 
               'v2'*
                  'BulbBeam'==
                     'ia(5.0,10.0)' 
                 
                  'BulbDepth'==
                     'ia(-10.1,10.0)' 
                  
            #ooo, that's interesting.  It's a real relation.
            # furthermore, python + PStates have seen fit to 
            # construct dummy variable connectives between binary compositions
            # (binary compositions <=?=> ternary relation  => not sure of my words here)
            # construct as before.
            # Compile as before.
            
            st, vars_ = CBulbMid.construct(CBulbMid, st, 
                                                existing=vars_)
            dflist = BulbBeam.compile(CBulbMid,vars_)
            code_list.append(dflist)
            
        
        #we now we have a list of rules (relational 'filter functions') 
        #                       composed of smaller relations that all work together 
        #                       so that the code of an entire relation is always called 'with itself'
        # and state, st.
        #
        # so let's try to filter the state through this!
        
        >> st
        >> 
        States([{
          v3:v3  
          BulbBeam:BulbBeam  
          A_mid:A_mid  
          BulbDepth:BulbDepth  
          CBulbMid:CBulbMid  
          v2:v2  
        }])
        
        for rule in code_list:
            for key in rule:
                st = rule[key](st)
                
        >>> st
        >>>  
        States([{
          v3:ia(0.0,1.0)  
          BulbBeam:ia(5.0,10.0)  
          BulbDepth:ia(-10.1,10.0)  
          A_mid:A_mid  
          CBulbMid:ia(0.0,1.0)  
          v2:ia(-101.0,100.0)  
        }])
    
        ## hmm... one more pass?:
    
        for rule in code_list:
            for key in rule:
                st = rule[key](st)
            
        >>> st
        >>> States([{
          v3:ia(0.0,1.0)  
          v2:ia(-101.0,100.0)  
          A_mid:ia(-101.0,100.0)  
          BulbDepth:ia(-10.1,10.0)  
          CBulbMid:ia(0.0,1.0)  
          BulbBeam:ia(5.0,10.0)  
        }])
    
        #it seems that's as far as we go.
        #  It is the purpose of the Rules Graph Processor
        # to automate this kind of thing (AC_revise), and to 
        # make it more efficient, together with the compiler, which adds links 
        # from mkVariable to it's own flist of rules it's involved in
    """
    count = 0
    def __init__(self, 
                 args=[],
                 op=None, 
                 opstr=None, 
                 name=None,
                 value=None,
                 verbose=False):
        """ 
        """
        if op is None: 
            op=PStates.__init__
            opstr = 'init'
        self.verbose = verbose
        self.value = value
        self.op = op
        self.opstr = opstr
        self.args = args
        self.nodelist = [] #list of this vars Big Rules {}
        if name is None:
            PStates.count += 1
            self.name = 'v'+str(PStates.count)
            self.discount = True
        else:
            self.name = name
            self.discount = False
        
            
        
    def __str__(self, level=0):
        """Fancy printer...
        walks the computational graph, printing nodes and their 
        associated ops
        """
        ret = "   "*level+repr(self.name)+""+self.op.__doc__+"\n"
        for el in self.args:
            if isinstance(el,Number):
                dc = PStates(value=el)
                dc.name = str(el)
                #ret += "   "*level+str(el)
                ret += dc.__str__(level+1)
            else:
                ret += el.__str__(level+1)
        return ret
    
    def __repr__(self):
        return 'PS({})'.format(self.name)
    
#    @property
#    def value(self):
#        return self._value
#    @value.setter
#    def value(self, v):
#        #print 'self = ',self
#        #print 'v    = ',v
#        self.value = v
#        return #self #.__eq__(v)
    
    def equal(self, v):
        """=="""
        return PStates(op=PStates.__eq__,
                       opstr = '==',
                       args = [self,v],
                       name = self.name)#,
                       #value = v)
    
    
    def __eq__(self, v):
        """=="""
        #print 'v',v
        return PStates(op=PStates.__eq__,
                       opstr = '==',
                       args = [v],
                       name = self.name)
        
    def __le__(self, v):
        """<="""
        return PStates(op=PStates.__le__,
                       opstr = '<=',
                       #args = [v], #modified __le__ to make this a reality
                       args = [self,v], #1.)-args would be better if it were just v here-
                       name = self.name)#2.)pass the name to avoid dummy name issues!
                                        #and make this an 'official setter' part of the language
        
    def __ge__(self, v):
        """>="""
        return PStates(op=PStates.__ge__,
                       opstr = '>=',
                       #args = [v], 
                       args = [self,v], 
                       name = self.name)
        
    def __add__(self, v):
        """+"""
        return PStates(op=PStates.__add__,
                       opstr = '+',
                       args = [self,v])
    def __sub__(self, v):
        """-"""
        return PStates(op=PStates.__sub__,
                       opstr = '-',
                       args = [self,v])
    def __mul__(self, v):
        """*"""
        return PStates(op=PStates.__mul__,
                       opstr = '*',
                       args = [self,v])
    def __div__(self, v):
        """/"""
        return PStates(op=PStates.__div__,
                       opstr = '/',
                       args = [self,v])
        
        
    
    #    def query(self, treenode, st, existing=None):
    #        #TODO: if not exist then make Var with name
    #        """Takes a computation tree varaible
    #        and returns it's Logic variable and value
    #        in the environment
    #        
    #        """
    #        if existing is None:
    #            existing = {}
    #        #if treenode.name == 'unknown_var':
    #        #if isinstance(treenode, Number):
    #        #    cur_el = treenode
    #        #    existing[treenode.name] = cur_el
    #        
    #        #scary:?
    #        if treenode.name not in existing:
    #            cur_el = st.bind(treenode.name)
    #            existing[treenode.name] = cur_el
    #        
    #        #if hasattr(treenode, 'args'):
    #        #    for el in treenode.args:
    #        #        self.construct(el, st, existing)
    #        return st, existing
        
    
    def construct(self, treenode, st, existing=None,
                  tracked_vars=None):
        #TODO: tracked_vars is no longer used.  Use existing instead.
        #       remove tracked_vars
        #TODO: if not exist then make Var with name
        """Takes a computation tree built by overloading
        and returns 2 things:
        (1) vars: a set of sqKanren.Variable class (logic variables) for it
        (2) environment: sqKanren.States class: the current environment
        that the computation tree works with
        
        
        inputs
        --------------
        treenode    : PStates object, possibly after building some
                        computational graph by using it in an expression
                        with other PStates objects or interval constants
                        
        st          : States object, an environment in which the 
                        results of assignments and other computations 
                        filters through.
                        
        existing    : if one already has started composing a dictionary of
                        vars_ in use, then pass it in here. The map goes like
                        existing = {PStates.name : Variable}
                        where Variable is the Variable class from sqKanren.py
                        
        tracked_vars : Nothing.  Testing for safe deletion.
                        
        
        returns
        --------------
        st          : States object, updated with new data constructed here
        
        existing    : updated map of PStates object names to Variables.
        
        
        TODO 
        --------------
        add binder for values?  (not really needed)
        
        
        side notes
        --------------
            self doesn't do anything in this construct function
            -could make this a classmethod.
        
        """
#        trackthem = False
#        if tracked_vars is not None:
#            trackthem = True
        if existing is None:
            existing = {}
        #if treenode.name == 'unknown_var':
        if isinstance(treenode, Number):
            cur_el = treenode
            existing[treenode.name] = cur_el
#            if trackthem: 
#                if cur_el not in tracked_vars:
#                    tracked_vars.append(cur_el)
        elif not treenode.name in existing:
            cur_el = st.bind(treenode.name)
            existing[treenode.name] = cur_el
#            if trackthem: 
#                if cur_el not in tracked_vars:
#                    print 'adding ',cur_el
#                    tracked_vars.append(cur_el)
#                else:
#                    print 'not adding ',cur_el
#                    print existing
        if hasattr(treenode, 'args'):
            for el in treenode.args:
                self.construct(el, st, existing)
        return st, existing
    
    
        
    def compile(self, treenode, vars_, funclist = None, i=0, st=None):
        """construct the space with 
        logic vars (use construct, above)
        before calling Goal.eval to 
        construct the rules themselves.
        
        inputs
        --------------
        treenode,   : A node in a computational graph.
                        If a user is inputing this node, it
                        should be the top node of the 
                        'PStates style' rule expression.
                        
        vars_       : dict map from {PStates.names : Variables}
                        aka from the name of a computational graph node 
                        to a sqKanren logical variable
                        
        funclist    : dictionary of rules, possibly empty to start
        
        i           : is immediately overwritten, slated for removal
        
        st          : States object (aka environement aka rgp.env)
        
        returns
        --------------
        funclist : list of rules dictionaries.  Each rule dictionary
                    contains all the individual (atomistic in fact - as far as relations go) 
                    rules used in composing the overall rule expression
                    each atomistic rule is a function closure 
                    from the Goal class
                    which awaits a States object.  Upon recieving one,
                    it will filter the State and return 
                    a new States, narrowed where possible by the
                    interval relational rule.
        
        
        TODO:  fix so NUMBERS dont require NAMES(?)
        
        side notes
        --------------
            self doesn't do anything in this compile function
            -could make this a classmethod.
        """
        #TODO: Write and use a new ``States'' where states are
        # passed in, filtered, and returned
        # make rules lazy (again)
        #TODO: delete the links of the graph after compiling
        #so that the next time a PStates var is used in a relation 
        # (at PStates level)
        #the relation doesn't get a second 'shot' of another rule!
        # this is one possible source of inefficiency.
        # though perhaps a long shot at fixing the 'hang-ups'
        # which happen from time to time.
        
        #idea:
        nconectives = 0 #track num of times a rule's functions
        # needs to be called to fully propogate the info
        
        #idea:
        # rules are saved with fixed ia(_,_) 
        # these are then creeping back into the env after having been
        # reduced out?
        
        if funclist is not None:
            i=len(funclist)
            #supernode name here - name for the rule at top of tree?
            #probably better to reverse this into the
            #mk_var arcs at "GraphScaleNode time"
        else:
            funclist = {}
            i=len(funclist)
        #print '\n\ntreenode data...'
        #print 'treenode.name = ',treenode.name
        #print 'treenode.args = ',treenode.args
        #print 'treenode = ',treenode
        if len(treenode.args)>0: #not PEP8!
            if len(treenode.args)==1:
                op = treenode.opstr
                operands = (vars_[treenode.args[0].name],)
                if self.verbose: 
                    print treenode.args[0].name
                    print treenode
                    print 'ops (1):-',op, operands
            if len(treenode.args)==2: #not PEP8!.. TODO: write tests and fix such as this
                op = treenode.opstr
                operands = (vars_[treenode.args[0].name],
                            vars_[treenode.args[1].name])
    #            operands = (vars_[treenode.args[0]],#.name],
    #                        vars_[treenode.args[1]])#.name])
                if self.verbose: print 'ops (2):-',op, operands
            result = vars_[treenode.name]
            if self.verbose: print 'ops (F):=>','(',op, operands, result,')'
            func = Goal.eval(operands, op, result) #turn the rule into a function
            #if func is not None:
            assert(func is not None)
            funclist[i] = func #put the function on the 'list'
            # alternative is to run the rule gen system every time?
            #
            for el in treenode.args:
                #if not isinstance(el, Number):
                #if  isinstance(el, PStates) or isinstance(el, Variable):
                if hasattr(el, 'args'):
                    funclist = self.compile(el, vars_, funclist, i+1)
        
        for var in vars_:
            mk_var = vars_[var]
            if isinstance(mk_var, Variable):
                #mk_var.flist.append(funclist) #code list onto mkvar
                for func in funclist:
                    fi = funclist[func]
                    mk_var.arcs.append(fi)
                    
                    
        """
        
        self.current_rules.append( 
                        GraphScaleNode(dflist, 
                                       vars_, 
                                       rulename))
        """
        
        #these are not gauranteed to remain constant...
#        for var in vars_:
#            print 'appending funclist to ',var,'?'
#            if isa(var,str):
#                pass
#            else:
#                var.nodelist.append(GraphScaleNode(funclist))
        
        """
        the million dollar question is (was),
        
        are the db elements??
        
        I mean:
        for var in vars_:
            for func in funclist:
                el = vars_[var]
            
        That is, the el themselves?
        (turtles = turtles ==  all( the, way, down) )
        """
        return funclist
    
    
    


##
##*****************************************************
##        
"""
ps = PStates(name = 'ps')
a = PStates(name='a')
b = PStates(name='b')
c = PStates(name='c')
d = PStates(name='d')
ps1 = (ps == (a,b) )

c =  ( ps + (a,b,c) )

ps1 = ( ps == (d,(ps + (a,b,c)) ) )

ps1.name = 'ps1'

ps2 = ( ps == (a,(ps - (b,c))))
"""
    
    
class States(object):
    """
    for use in/with...
        -design_trees as design node
        -HullCLP types
        -lists
        
        Mathmatical functions MUST be written so that:
            
            a * b <=> c
            
        NEVER expect or write anything like:
            
            a <=> b * c
        
        BECAUSE 
            the prefix notation currently used will
            not account for it!
            
        Note
            If you use the language built on top of this though,
            i.e. mathematical composition of PStates 
            and rules compilation,
            Then any statement form is fair game.
    """
    def __init__(self,states=None,
                 parent=None, children=None):
        self.parent = parent 
        self.dummy_count = 0
        if children is None:
            self.children = []
        else:
            self.children = children
        if isinstance(states, State):
            self.states = [states]
        elif isinstance(states, list):
            self.states = states
        elif isinstance(states, tuple):
            self.states = list(states)
            
    def bind_if_not_extant(self, var):
        for state in self.states:
            state.bind_if_not_extant(var)
        return
    
    def eval(self, eval_tuple ):
        """States.eval
        convert a symbol to a operation
        usage:
            st.eval((x,st.__eq__,(x,ia(5.,6.))))
            st.eval((x,'==',(x,ia(5.,6.))))
            
            st.eval( (x,'==',[(y,'==',[c,ia(0.3,2.2)]),ia(1.,2.)]) )
            -sets c = ia(.3,2.2) in these states
            -sets y = ia(1.,2.) in these states
            -returns x
        
        notes:
            This is mostly syntatic sugar
            - not used for the rules compiler
            - is used as a sort of hybrid miniKanren 'run' device
            - this was a transitional piece of code for me
            - see examples for usage.  maybe not totally useless.(redundant?)
        """
        a, rule, v = eval_tuple
        if not isa(rule, str):
            #temp = rule(v)
            return rule(v), a #temp(a)
        if rule == '==':
            rule = self.__eq__
        elif rule == '+':
            rule = self.__add__
        elif rule == '-':
            rule = self.__sub__
        elif rule == '*':
            rule = self.__mul__
        elif rule == '/':
            rule = self.__div__
        elif rule == '<=':
            rule = self.__le__
        elif rule == '>=':
            rule = self.__ge__
        elif rule == 'split':
            rule = self.split
        #temp = rule(v)
        return rule(v), a #temp(a)
        

    def __repr__(self):
        return 'States({})'.format(self.states)

    def __str__(self):
        return 'States({})'.format(self.states)

    def __call__(self, x):
        d = []
        for st in self.states:
            d.append(st(x))
        return d

    def get_latest_child(self):
        lc = len(self.children)
        #nst = len(self)
        if lc == 0:
            return self#States(self.states)#, nst
            #return lp.States(self[:])#, nst
        else:
            return self.children[0].get_latest_child()
    
    
    def update(self, dict_obj):
        """
            update to take dict_obj 
            or State class
        """
        for st in self.states:
            st.values.update(dict_obj)
        return
    
    
    
    def bind(self, vars_):
        if isa(vars_, list):
            for var in vars_:
                assert isa(var, str), 'string binding from list not supported'
                #var = self.bind_str(var)
                #assert isa(vars_, Variable), 'vars_ {} is not string'.format(vars_)
                #self.update({var:None})
                self.bind(var)
        elif isa(vars_, Variable):
            self.update({vars_:None})
        else:
            var_ = self.bind_str(vars_)
            return var_
        return
    
    
    def bind_str(self, var):
        assert isa(var, str), 'vars_ {} is not string'.format(var)
        instates = True
        for state in self.states:
            names = [ el.name for el in state.values.keys()]
            if var not in names:
                instates = False
        if not instates: #True:#
            #not found so make lp.Variable out of var
            var_ = Variable(var)
            self.bind(var_)
        else:
            #found by string name, do not duplicate var
            for el in state.values:
                if var == el.name:
                    var_ = el
                    break
        return var_
    
    def get_by_str(self, var):
        assert isa(var, str), 'vars_ {} is not string'.format(var)
        found = False
        var_ = None
        for state in self.states:
            if not found:
                for el in state.values:
                    if var == el.name:
                        var_ = el
                        found = True
                        break
        
        return var_
        
    
    
    def get_parent(self):
        return self.parent
    def get_child(self):
        return self.children

    def eq(self, v):
        """#print 'equality ',a,'?:=',b"""
        states = []
        for s in self.states:
            sn = Goal.eq(v[0],v[1])(s)
            for sni in sn:
                if sni is not None:
                    states.append(sni)
        return States(states)
        
        
    #    def diff(self,v):
    #        #        states = []
    #        #        for s in self.states:
    #        #            sn = Goal.diff(v[0],v[1],v[2])(s)
    #        #            for sni in sn:
    #        #                if sni is not None:
    #        #                    states.append(sni)
    #        if not isinstance(v[0], Variable) and not isinstance(v[1], Variable):
    #            dl = abs(v[0]-v[1])
    #            tol = v[2]
    #            if dl>tol:
    #                return True
    #            else:
    #                return False
    #        else:
    #            return False
    #        #return States(states)
    def recurse(self, el):
        """
        usage:
            
        def el := (x,'==',[(y,'==',[c,ia(0.3,2.2)]),ia(1.,2.)])
        
        """
        if not isa(el,(Variable,ia,Number)):
            return self.eval(el)
        else:
            return self, el
    def process_args(self, v):
        """
        evaluate a list of arguments,
        return modified state and a desired list of new arguments
        args in use:
        v = [(x,'==',[x,ia(2.,4.)]),
             (y,'==',[y,ia(-2.,2.)]),
             (z,'==',[z,ia(-10.,10.)])]
        
        if this were not called at the front of every evaluation,
        each function would do what process_args does. 
        Namely, 
            evaluate the arguments of 'v':
                
                if not isa(v[0], (Variable,ia,Number)):
                    self,v[0] = self.eval(v[0])
                if not isa(v[1], (Variable,ia,Number)):
                    self,v[1] = self.eval(v[1])
                if not isa(v[2], (Variable,ia,Number)):
                    self,v[2] = self.eval(v[2])
        """
        if isa(v,list):
            for i,el in enumerate(v):
                self, v[i] = self.recurse(el)
        return self,v
    
    def split(self, v):
        print 'spliting...'
        print 'v0 = ',v[0]
        print 'v1 = ',v[1]
        print 'v2 = ',v[2]
        states = []
        for s in self.states:
            s1 = Goal.eq(v[0],v[1])(s)
            s2 = Goal.eq(v[0],v[2])(s)
            states += s1
            states += s2
        return States(states)
    
    #@process_args
    def __eq__(self, v):
        """Not a comparison operator
        instead tries to make v[0] = v[1]
        
        hmm: how does this work with PStates parser?
        there, self.args = [v]
        i.e. no self in args of PStates
        """
        """#print 'equality ',a,'?:=',b"""
        self,v = self.process_args(v)
        states = []
        for s in self.states:
            sn = Goal.eq(v[0],v[1])(s)
            for sni in sn:
                if sni is not None:
                    states.append(sni)
        return States(states)
    
    
    def __le__(self, v):
        """
        currently built to requst7 = (st6 == (d,c)) 
        self.args = [self,v]
        in the PStates parser
        
        
        usage:
            
        say we have
            y == ia(2.,4.)
            z = z
        then apply: st = (st <= (z,y))
        result:
            z == ia(0.,4.)
        """
        """#print 'inequality ',a,'?:<=',b"""
        self,v = self.process_args(v)
        states = []
        for s in self.states:
            sn = Goal.lto(v[0],v[1])(s)
            for sni in sn:
                if sni is not None:
                    states.append(sni)
        return States(states)
    
    def __ge__(self, v):
        """#print 'inequality ',a,'?:>=',b"""
        self,v = self.process_args(v)
        states = []
        for s in self.states:
            sn = Goal.gto(v[0],v[1])(s)
            for sni in sn:
                if sni is not None:
                    states.append(sni)
        return States(states)

    #@staticmethod
    def both_states(self, a,b):
        s=[]
        for s1 in self.states:
            st1 = a(s1)
            for st2 in st1:
                s.append(b(st2))
        return States(s)

    def both(self, a,b):
        """not needed
        """
        s=[]
        #for s1 in self.states:
        st1 = a(self.states)
        #for st2 in st1:
        s += b(st1.states).states
        return States(s)


    def __add__(self, v):
        """addo(v1,v2,v3)(state)
        """
        #stuff replaced by func process_args:
        #        if not isa(v[0], (Variable,ia,Number)):
        #            self,v[0] = self.eval(v[0])
        #        if not isa(v[1], (Variable,ia,Number)):
        #            self,v[1] = self.eval(v[1])
        #        if not isa(v[2], (Variable,ia,Number)):
        #            self,v[2] = self.eval(v[2])
        self,v = self.process_args(v)
        states = []
        for s in self.states:
            sn = Goal.addo(v[0],v[1],v[2])(s)
            for sni  in sn:
                if sni is not None:
                    states.append(sni)
        return States(states)

    def __sub__(self, v):
        """subo(v1,v2,v3)(state)
        """
        self,v = self.process_args(v)
        states = []
        for s in self.states:
            sn = Goal.subo(v[0],v[1],v[2])(s)
            for sni  in sn:
                if sni is not None:
                    states.append(sni)
        return States(states)
    
    def __mul__(self, v):
        """mulo(v1,v2,v3)(state)
        """
        self,v = self.process_args(v)
        states = []
        for s in self.states:
            sn = Goal.mulo(v[0],v[1],v[2])(s)
            for sni  in sn:
                if sni is not None:
                    states.append(sni)
        return States(states)

    def __div__(self, v):
        """extended division
        on relationl intervals
        this way doesn't work
        """
        self,v = self.process_args(v)
        states_ = []
        for s in self.states:
            sn = Goal.extended_divo(v[0],v[1],v[2])(s)
            for sni  in sn:
                if sni is not None:
                    states_.append(sni)
        return States(states_)
        
    
    def __pow__(self, v):
        """order must make some sense here
        or we must interpret the input
        """
        self,v = self.process_args(v)
        states = []
        for s in self.states:
            sn = Goal.powo(v[0],v[1],v[2])(s)
            for sni  in sn:
                if sni is not None:
                    states.append(sni)
        return States(states)
    

    def __div__NOGOOD(self, v):
        """extended division
        on relationl intervals  .works.
        .and. works with:
        state = goal(x1,y1,z)(state)
        state = goal(x1,y,z1)(state)
        state = goal(x,y1,z1)(state)
        state = goal(x1,y1,z1)(state)
        """
        states = []
        for s in self.states:
            v0 = s(v[0])
            v1 = s(v[1])
            v2 = s(v[2])
            sn = Goal.extended_division(v0,v1,v2)(s)
            #sn = Goal.extended_division(v[0],v[1],v[2])(s)
            for sni  in sn:
                if sni is not None:
                    states.append(sni)
        return States(states)

    def __add__NOGOOD(self, v):
        """addo(v1,v2,v3)(state)
        """
        states = []
        for s in self.states:
            v0 = s(v[0])
            v1 = s(v[1])
            v2 = s(v[2])
            sn = Goal.add(v0,v1,v2)(s)
            for sni  in sn:
                if sni is not None:
                    states.append(sni)
        return States(states)

    def __sub__NOGOOD(self, v):
        """subo(v1,v2,v3)(state)
        """
        states = []
        for s in self.states:
            v0 = s(v[0])
            v1 = s(v[1])
            v2 = s(v[2])
            sn = Goal.sub(v0,v1,v2)(s)
            for sni  in sn:
                if sni is not None:
                    states.append(sni)
        return States(states)

    def __mul__NOGOOD(self, v):
        """mulo(v1,v2,v3)(state)
        """
        states = []
        for s in self.states:
            v0 = s(v[0])
            v1 = s(v[1])
            v2 = s(v[2])
            sn = Goal.mul(v0,v1,v2)(s)
            for sni  in sn:
                if sni is not None:
                    states.append(sni)
        return States(states)

    @staticmethod
    def both_old(self, a,b):
        """
        Problem:

        These:
        st6 = st.both((st + (x,y,c) ),
                      (st + (x,y,d)))

        are already evaluated
        """
        def pursue(state):
            states = []
            for a_state in a(state):
                states += b(a_state)
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


class StateTree(States):
    def __init__(self, State=None):
        if isinstance(states, State):
            self.states = [states]
        else:
            print 'Error, this State Tree must be initialized'
            print 'with a single state'
            
class State(object):

    def __init__(self, values=None):
        self.values = values or {}#OrderedDict()
        self.valid_design = True
        self.excludes = ['c1','c2','c3','c4','c5']
        #self.vars = []

    def __repr__(self):
        st = ''
        for val in self.values:
            if val.name not in self.excludes:
                if isinstance(val, ia):
                    pass
                else:
                    #st+='  {}=>{}  \n'.format(val,self.values[val])
                    st+='  {}:{}  \n'.format(val,self.value_of(val))
        return '{\n'+st+'}'

    def __call__(self, x):
        return self.value_of(x)
    
    def bind_if_not_extant(self, var):
        if var in self.values:
            pass
        else:
            self.update({var:None})
        return
    
    def update(self, dict_obj):
        if isinstance(dict_obj, dict):
            self.values.update(dict_obj)
        else:
            for el in dict_obj:
                self.values
        return
    
    #    def split(self, var,x):
    #        return
    
    def update_or_instate(self):
        print 'update_or_instate:'
        print 'think harder about this!'
        return

    def __add__(self, other):
        """s + (x,y,c)
        """
        return self.state + other#Goal.addo(v[0],v[1],v[2])(self)

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
    
    
    ## Oct 4, 2017, just get things done.
    def get_by_str(self, var):
        assert isa(var, str), 'vars_ {} is not string'.format(var)
        #found = False
        var_ = None
        for el in self.values:
            if var == el.name:
                var_ = el
                break
        
        return var_
    
    #def exist_ckvar(self, var):
        

    def fresh(self, vars, func=None):
        values = copy.copy(self.values)
        for var in vars:
            var = Variable(var)
            values[var] = None
        if func is not None:
            #TLM March 2017:
            nself = copy.copy(self)
            nself.values = values
            return func(nself)
            #return func(self)
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

    def unify(self, a, b): #moved simple var lookups downwards
        """The Heart of a Logic Programming System
        """
        ivars = a,b
        a = self.value_of(a)
        b = self.value_of(b)
        #print 'ivar = ',ivars
        #print 'a = ',a
        #print 'b = ',b
        if a == b:
            return self
        #
        elif isinstance(a, ia) and isinstance(b, ia):
            values = copy.copy(self.values)
            test = a&b
            if test.isempty:
                #print 'empty intersection!'
                #print 'a={}, b={}, test={}, type={}'.format(a,b,test, type(test) )
                return None
            #else:
            #assert(ivars[0] in values.keys())
            values[ivars[0]] = test
            values[ivars[1]] = values[ivars[0]] & b
            values[ivars[0]] = values[ivars[0]] & values[ivars[1]]
            #print 'success!, values = ',values[ivars[0]], values[ivars[1]]
            return State(values)
        elif isinstance(a, ia) and isinstance(b, Variable):
            values = copy.copy(self.values)
            values[ivars[1]] = a
            return State(values)
        elif isinstance(b, ia) and isinstance(a, Variable):
            values = copy.copy(self.values)
            values[ivars[0]] = b
            return State(values)
        elif isinstance(a, Variable): #simple varlookup
            return self.assign(a, b)
        elif isinstance(b, Variable): #simple varlookup
            return self.assign(b, a)
        elif isinstance(a, tuple) and isinstance(b, tuple):
            unify_left = self.unify(a[0], b[0])
            if unify_left:
                unify_right = unify_left.unify(a[1], b[1])
                return unify_right
            else:
                return None
        elif isinstance(a, list) and isinstance(b, list):
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
        
        #print ivars
        #print a
        #print b
        if a == b:
            #print '1'
            return self
        elif isinstance(a, Variable) and isinstance(b, ia):
            #print '2'
            values = copy.copy(self.values)
            values[ivars[0]] = ia(min(0.,b.inf),b.sup)
            return State(values)
        elif isinstance(b, Variable) and isinstance(a, ia):
            #print '3'
            values = copy.copy(self.values)
            values[ivars[1]] = ia(min(0.,a.inf) ,a.sup) #assumes + only!
            return State(values)
        elif isinstance(a, ia) and isinstance(b, ia):
            #print '4'
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
            #print '5'
            return self

    def greaterthanify(self, a, b, option=None):
        """ a >= b
        for intervals (a,b)

        is intersection ncessary?
        """
        ivars = a,b
        a = self.value_of(a)
        b = self.value_of(b)
        
        
        #print ivars
        #print a
        #print b
        if a == b:
            #print '1'
            return self
        elif isinstance(a, Variable) and isinstance(b, ia):
            #print '2'
            #values = copy.copy(self.values)
            #values[ivars[0]] = ia(b.inf,max(0.,b.sup))
            return self#State(values)
        elif isinstance(b, Variable) and isinstance(a, ia):
            #print '3'
            #values = copy.copy(self.values)
            #values[ivars[1]] = ia(min(0.,a.inf) ,a.sup) 
            return self#State(values)
        elif isinstance(a, ia) and isinstance(b, ia):
            #print '4'
            values = copy.copy(self.values)
            c = a&b
            if c.isempty and a>b:
                values[ivars[0]] = a
                values[ivars[1]] = b
            else:
                values[ivars[0]] = ia(max(a.inf,c.inf),max(a.sup,c.sup))
                values[ivars[1]] = ia(min(b.inf,c.inf),min(c.sup,b.sup))
            return State(values)
        else:
            #print '5'
            return self

    def greaterthanify_old(self, a, b, option=None):
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
    """
    if isinstance(e,Variable):
        return reify(s[e], s) if e in s else e
    return _reify(e, s)



def binop(op, revop=None):
    """ Transform binary operator into goal
    """

    def goal(x, y, z):
        if not isinstance(x, Variable) and not isinstance(y, Variable):
            #print 'binop 1: {}*{} = {}'.format(x,y,z)
            return Goal.eq(op(x, y), z)
        
        if not isinstance(y, Variable) and not isinstance(z, Variable) and revop:
            #print 'binop 2'
            return Goal.eq(x, revop(z, y))
        
        if not isinstance(x, Variable) and not isinstance(z, Variable) and revop:
            #print 'binop 3'
            return Goal.eq(y, revop(z, x))
        
        else:
            #print 'binop 4th - no answer'
            return Goal.eq(None,None)#(z, revop(x,y))#NULL#State() #raise EarlyGoalError()#null #

    goal.__name__ = op.__name__+'binop'
    return goal

add = binop(op.add, op.sub)
add.__doc__ = """ x + y == z """
#mul = binop(op.mul, op.truediv)
mul = binop(op.mul, op.div)
mul.__doc__ = """ x * y == z """
mod = binop(op.mod)
mod.__doc__ = """ x % y == z """
#divex = binop(op.div, op.mul)
#divex.__doc__ = """[x]/[y]=[z]"""
pow_ = binop(ia.__pow__, ia.__pow__)
pow_.__doc__ = """  x**y == z  """


split = ia.split_operator



def null(s):
    return Goal.eq(None,None)(s) #[NULL]



def take(cgoal, states):
    
    return 
    
    
def run(n,x, states):
    """ Run a logic program.  Pass a variable through 
    a goal unified down to state to get a result
    the state can be the result of a relational program
    
    x = Variable('x')
    y = Variable('y')
    c = Variable('c')
    d = Variable('d')
    s = State(values={x:None,y:None,c:None,d:None})
    st = States(s)
    
    
    st = ((st == (x,ia(1.,2.))) == (run(0, x, st  == (x,ia(1.,2.))),y))
    
    
    st = (st == (run(0, x, st  == (x,ia(1.,2.))),y))
    st = (st == (run(0, x, st  == (x,ia(1.,2.))),x))
    
    Not Really Used, but something like this
        would be used to make a miniKanren style machine
    """
    
    return states(x)[n]



def run_ck_relation(vars, state, goalfun):
    """not used
    """
    #valsi = [state(var) for var in vars]
    statei = copy.copy(state)

    #(was?) ternary op, now list comprehension
    chg_lst = [var for var in vars if not state(var) == statei(var)]

    return chg_lst

#def run_unary_relation(x,y, state, goal):#, option=None):
#    x1 = state.value_of(x)
#    y1 = state.value_of(y)
#    #state = goal(x1,y)(state)[0]
#    #state = goal(x,y1)(state)
#    state = goal(x1,y1)(state) #equiv?why
#    return state

def run_binary_relation(x,y,z, state, goal):#, option=None):
    x1 = state.value_of(x)
    y1 = state.value_of(y)
    z1 = state.value_of(z)
    #print 'run bin relation:'
    #print x,y,z
    #print 'states are:'
    #print x1,y1,z1

    states = goal(x1,y1,z1)(state)
    ns = []
    while None in states: states.remove(None)
    for s in states:
        #print 's = ',s
        ns += goal(x1,y1,z)(s)
    if len(ns)>0:
        states = copy.copy(ns)
    ns=[]
    while None in states: states.remove(None)
    for s in states:    
        #print 'states len = ',len(states)
        #print 's = ',s
        ns += goal(x1,y,z1)(s)
    if len(ns)>0:
        states = copy.copy(ns)
    ns=[]
    while None in states: states.remove(None)
    for s in states:
        #print 's = ',s
        ns += goal(x,y1,z1)(s)
    if len(ns)>0:
        states = copy.copy(ns)
    while None in states: states.remove(None)
    #ns=[]
    #for s in states:
    #    ns +=goal(x1,y1,z1)(state)#[0]
    #states = copy.copy(ns)
    return states


def run_binary_relation__(x,y,z, state, goal):#, option=None):
    x1 = state.value_of(x)
    y1 = state.value_of(y)
    z1 = state.value_of(z)
    #if isinstance(z1, Variable):
    states = goal(x1,y1,z)(state)
    #if isinstance(y1, Variable):
    ex = 0
    for i, s in enumerate(states):
        test = goal(x1,y,z1)(s)
        if len(test)==1:
            states[i+ex] = test
        else:
            states[i+ex] = test.pop(0)
            for el in test:
                states.insert()
                ex+=1

    #if isinstance(x1, Variable):
    for s in states:
        states += goal(x,y1,z1)(s)
    for state in states:
        states += goal(x1,y1,z1)(state)#[0]
    return states
#def run_all_binary_relations(x,y,z, state, goal):#, option=None):
#    x1 = state.value_of(x)
#    y1 = state.value_of(y)
#    z1 = state.value_of(z)
#    print x1
#    print y1
#    print z1
#    states = goal(x1,y1,z)(state)
#    #s1 = []
#    #for s in states:
#    #    s1.append(goal(x1,y,z1)(state))
#    #state = goal(x,y1,z1)(state)
#    return s1#state

def run_br(x,y,z, state, goal):
    return run_binary_relation(x,y,z, state, goal)

run_binary_relation.__doc__ = """ x1 o x2 <=> x3 with state and goal"""
run_br.__doc__ = """ x1 o x2 <=> x3 with state and goal"""



class Goal(object):
    """
    To add an operator to the language 
    you will need to add it to eval here.
    
    
    
    
    Notes (continuation implemenation aplogetics)
    ----------
    At some point somebody is going to say I'm using 
    broken features of Python.  They are making a mistake, 
    at least as far as the implementation of Goals is concerned.
    
    
    Take mulo, a goal for relatoinal multiplication implemented in this
    Goals class as follow:
        
        
    @staticmethod
    def mulo(x,y,z):
        def pursue(state):
            s = run_binary_relation(x,y,z,state,Goal.mul)
            return s
        return pursue
    
    
    Often smart folk want to talk about how Python scoping is broken.
                (lexical, I believe, sometimes, but maybe dynamical at others)

    
    -One cannot mutate x,y, or z, in a function _like_ mulo
    -Good thing we wouldn't attempt such a thing by accident.  
                (anywehere, right?.. right??)
    -instead the idea with miniKanren is to mutate... nothing at all!
                (logic programming is quite functional, after all)
    -so what is this mulo an instance of?
                (closure?  continuation?)
                
                summary of the below, its a closure.  Norvig is my answer.
                It doesn't really do control flow.
                It doesn't tell the program how to 'continue' computation.
                e.g. we don't feed a funtion to the continuation to stear
                the continued computation.  We just close the function.
                
    -so why is this okay?
        *we do not mutate x,y, or z.  we just feed them in 
        at the creation of the closure.
    
    Details (contunuations and closures and how this is not a thunk)
    ----------     
    continuation:  
                    -this may not be a continuation 
                        because mulo really does return.  It just happens
                        to return a function.
                    -A real continuation calls a function, instead of
                        just returning it.
                    -We do not call the final function, we've
                        just composed it and returned it.
                    Oh, a snag?: we will be calling it soon enough.
                        Which brings us to closures:
    closure:        
                    'closures are simply the current lexical stack'
                        That is, we capture data from the environment
                        (x,y,z)
                        
                    Norvig PAIP, page 280
                    closures: specify the something as a function, and
                    call the function at a later time.
                    
                    PAIP, page 92:  a closure is a function coupled with 
                    free lexical variables.  (Ours are the states to be
                    passed in later, to close the closure.)
    
    *https://stackoverflow.com/questions/25953545/function-closure-versus-continuation-in-general-and-sml?lq=1
    
    continuation is a more abstract concept
    
    *https://stackoverflow.com/questions/14019341/whats-the-difference-between-a-continuation-and-a-callback/14022348
    
    A continuation is actually a reification 
    of the control state of the program: 
    a snapshot of the state of the program 
    at a certain point in time. 
    The fact that it can be called 
    like a normal function is irrelevant. 
    Continuations are not actually functions.
     
    Callbacks on the other hand are actually functions. 
    That's the real difference between continuations and callbacks. 
    Nevertheless JS doesn't support first-class continuations. 
    Only first-class functions. 
    Hence continuations written in CPS in JS are simply functions. 
    
    
    
    Here we are:
        
        
        A continuation actually represents the instructions remaining 
        to complete a computation, or the remainder of the computation 
        from this point in time. 
        You can think of a continuation as a hole that needs to be filled in.
        
    
    Continuations are good in Python:
    *https://www.ps.uni-saarland.de/~duchier/python/continuations.html
    
    Is it a thunk?
    No, this not really a kind of a thunk:
        https://en.wikipedia.org/wiki/Thunk
        if only because mulo returns a funtion
        which takes states as an argument.  It is not totally encapulated
        it cannot be memoized.  These are my heuristic tells for a thunk
        at the moment, given my un-closured understanding... 
    If, instead of explicitly passing states
        to and return states from functions,
        we instead mutated state internally,
        Then we could use thunks just like we now use
        closures.  -I don't think this is the spirit of the thing.
    """
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
        #print 'goal equality ',a,'?:=',b
        def pursue(state):
            st = []
            if isinstance(a, list):
                for ai in a:
                    #print 'unify def:',state.unify
                    #print 'unify {}<=>{}'.format(ai,b)
                    st.append( state.unify(ai, b) )
                    #print 'good'
                if len(st)>0:
                    return st
                else:
                    return []
            elif isinstance(b, list):
                for bi in b:
                    st.append( state.unify(a, bi) )
                if len(st)>0:
                    return st
                else:
                    return []
            else:
                if state is not None:
                    new_state = state.unify(a, b)
                    if new_state:
                        return [new_state]
                    else:
                        return []
                else:
                    return []
        return pursue

    @staticmethod
    def lt_ia_(x, y):
        """ Actually implemented in lessthanunify
        x | x < y """
        a = x&y
        return ia(x.inf,min(a.sup,x.sup))

    @staticmethod
    def lto(a, b):
        """  Less than Constrant
        x | x < y """
        #print 'goal lto ',a,'?:=',b
        def pursue(state):
            st = []
            if isinstance(a, list):
                for ai in a:
                    st.append( state.lessthanify(ai, b) )
                if len(st)>0:
                    return st
                else:
                    return []
            elif isinstance(b, list):
                for bi in b:
                    st.append( state.lessthanify(a, bi) )
                if len(st)>0:
                    return st
                else:
                    return []
            else:
                new_state = state.lessthanify(a, b)
                if new_state:
                    return [new_state]
                else:
                    return []

        return pursue


    @staticmethod
    def gto(a, b):
        """  Less than Constrant
        x | x > y """
        #print 'goal gto ',a,'?:=',b
        def pursue(state):
            st = []
            if isinstance(a, list):
                for ai in a:
                    st.append( state.greaterthanify(ai, b) )
                if len(st)>0:
                    return st
                else:
                    return []
            elif isinstance(b, list):
                for bi in b:
                    st.append( state.greaterthanify(a, bi) )
                if len(st)>0:
                    return st
                else:
                    return []
            else:
                new_state = state.greaterthanify(a, b)
                if new_state:
                    return [new_state]
                else:
                    return []
        return pursue



    @staticmethod
    def diff(x, y, tol):
        """ x > y """
        if not isinstance(x, Variable) and not isinstance(y, Variable):
            diff = abs(x-y)
            return Goal.eq(diff > tol, True)
        else:
            raise Goal.EarlyGoalError()
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
            
    #    @staticmethod
    #    def eval(vals):
    #        """thinking about higher order funcs...
    #        >>> g = getattr(Goal, 'add')(x,y,z)(State)
    #        """
    #        getattr(self, vals)
    
    @staticmethod
    def freify(vals,dummy):
        """NOT USED
         thinking about higher order funcs...
        """
        dvar = Variable(vals[dummy])
        return
        
    @staticmethod
    def reify_eval(func, vals,index):
        """NOT USED
         thinking about higher order funcs...
         reify a rule into a spot
        rule:-
         new_state = (state * (_0_, x, y))
         new_state = (state * (reify_eval puts a number here, x, y))
        """
        dummy = Variable('-')
        vals[index] = dummy
        def eval_and_reify(state):
            return reify(dummy, func(vals)(state))
        return  eval_and_reify
            
    @staticmethod
    def funcify(func, vals):
        """NOT USED
        thinking about higher order funcs...
        parse the rule composition...
        
        rule:-
         new_state = (state * (state * (_0_, x, y)), C, D )
        """
        #head = vals[0]
        #tail = vals[1:]
        for var in vals:
            if isinstance(var, Iterable):
                Goal.funcify(var)
            if isinstance(var, _):
                dummy_index = vals.index(var)
                val = Goal.reify_eval(func, vals, dummy_index)
                vals[dummy_index] = val
            if callable(var):
                #return self.eval(head, tail)
                print 'no functions should be in the body!'
        #    if not isinstance(var, Variable):
                
        #a=isinstance(vals)
        return func(vals),dummy_index #not working yet - state will give dummy's value...

    #    @staticmethod
    #    def split(var,x):
    #        return split(var,x)

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
    def div(x, y, z):
        """ x / y == z """
        return mul(z, y, x)

    @staticmethod
    def divo(x,y,z):
        def pursue(state):
            s = run_binary_relation(x,y,z,state,Goal.div)
            return s
        return pursue


    @staticmethod
    def extended_division(x, y, z):
        """ x / y == z """
        return mul(z, y, x)#divex(x, y, z)
        #return eia.__div__(x, y, z) #cannot work unless x is extended ia

    @staticmethod
    def extended_divo(x,y,z):
        def pursue(state):
            #s = run_binary_relation(x,y,z,state,divex)
            s = run_binary_relation(x,y,z,state,Goal.extended_division)
            return s
        return pursue


    @staticmethod
    def pow_(x,y,z):
        """
            pow is not symmetric in its three operations
            so we parse them here.
            
            this will be a common operation as 
            relational arithemitic is extended to
            a larger class of mathemeatical operators
        """
        #print 'x = ',x
        #print 'y = ',y
        #print 'z = ',z
        #print 'checking pow type....'
        if not isinstance(x, Variable) and not isinstance(y, Variable):
            #print 'pow 1'
            #a = pow_(x,y,z)#, invert=True)
            #print 'ans = ',a
            return pow_(x,y,z)
        if not isinstance(y, Variable) and not isinstance(z, Variable):
            #print 'pow',2
            p = 1./y
            return pow_(z,p,x)
        if not isinstance(x, Variable) and not isinstance(z, Variable):
            #print 'pow',3
            p = np.log(z)/np.log(x)
            return Goal.eq(y,p)  #pow_(x,p,z)
        else:
            #print 'issue with pow'
            return pow_(x,y,z)
            
    @staticmethod
    def powo(x,y,z):
        def pursue(state):
            s = run_binary_relation(x,y,z,state,Goal.pow_)
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
    
    
        return States(states)
    
    
    """--------------------------------------------"""
    @staticmethod
    def eval(operands,op,result ):
        """Goal.eval
        usage:with PStates variables which self assemble 
        computational expression trees
        
        This takes a rule and turns it into a (Goal) 
        - a function which represents that rule
        
        question:  what about [var == ia()]?
        these need to take the latest info,
        always
        """
        #print 'op = ',op
        if op == '==':
            rule = Goal.__eq__
            func = rule( (operands[0],result) )
        elif op == '+':
            rule = Goal.__add__
            func = rule( (operands[0],operands[1],result) )
        elif op == '-':
            rule = Goal.__sub__
            func = rule( (operands[0],operands[1],result) )
        elif op == '*':
            rule = Goal.__mul__
            func = rule( (operands[0],operands[1],result) )
        elif op == '/':
            rule = Goal.__div__
            func = rule( (operands[0],operands[1],result) )
        elif op == '<=':
            rule = Goal.__le__
            #func = rule( (operands[0],operands[1],result) )
            func = rule( (result,operands[0]) ) #Jan 9 2018: flipped order here to make this work
        elif op == '>=':
            rule = Goal.__ge__
            #func = rule( (operands[0],operands[1],result) )
            func = rule( (result,operands[0]) ) #Jan 9 2018: flipped order here to make this work
        else:
            assert(op == 'init')
            func = None
        func.name = [op,operands,result]
        return func
    

    @staticmethod
    def __eq__(v):
        """Not a comparison operator
        instead tries to make v[0] = v[1]
        """
        """#print 'equality ',a,'?:=',b"""
        #states,v = states.process_args(v)
        def pursue(states):
            states_ = []
            for s in states.states:
                sn = Goal.eq(v[0],v[1])(s)
                for sni in sn:
                    if sni is not None:
                        states_.append(sni)
            return States(states_)
        return pursue
    
    @staticmethod
    def __le__(v):
        """#print 'inequality ',a,'?:<=',b"""
        #print 'inequality ',v[0],'?:<=',v[1]
        #sttes,v = states.process_args(v)
        def pursue(states):
            states_ = []
            for s in states.states:
                sn = Goal.lto(v[0],v[1])(s)
                for sni in sn:
                    if sni is not None:
                        states_.append(sni)
            return States(states_)
        return pursue
    
    @staticmethod
    def __ge__(v):
        """#print 'inequality ',a,'?:<=',b"""
        #print 'inequality ',v[0],'?:>=',v[1]
        #states,v = self.process_args(v)
        def pursue(states):
            states_ = []
            for s in states.states:
                sn = Goal.gto(v[0],v[1])(s)
                for sni in sn:
                    if sni is not None:
                        states_.append(sni)
            return States(states_)
        return pursue


    @staticmethod
    def __add__(v):
        """addo(v1,v2,v3)(state)
        """
        #states,v = self.process_args(v)
        def pursue(states):
            states_ = []
            for s in states.states:
                sn = Goal.addo(v[0],v[1],v[2])(s)
                for sni  in sn:
                    if sni is not None:
                        states_.append(sni)
            return States(states_)
        return pursue

    @staticmethod
    def __sub__(v):
        """subo(v1,v2,v3)(state)
        """
        #states,v = self.process_args(v)
        def pursue(states):
            states_ = []
            for s in states.states:
                sn = Goal.subo(v[0],v[1],v[2])(s)
                for sni  in sn:
                    if sni is not None:
                        states_.append(sni)
            return States(states_)
        return pursue
    
    @staticmethod
    def __mul__(v):
        """mulo(v1,v2,v3)(state)
        """
        #states,v = Goal.process_args(v)
        def pursue(states):
            states_ = []
            for s in states.states:
                sn = Goal.mulo(v[0],v[1],v[2])(s)
                for sni  in sn:
                    if sni is not None:
                        states_.append(sni)
            return States(states_)
        return pursue
    
    @staticmethod
    def __div__(v):
        """extended division
        on relationl intervals
        this way doesn't work
        """
        #states,v = Goal.process_args(v)
        def pursue(states):
            states_ = []
            for s in states.states:
                sn = Goal.extended_divo(v[0],v[1],v[2])(s)
                for sni  in sn:
                    if sni is not None:
                        states_.append(sni)
            return States(states_)
        return pursue
        
    @staticmethod
    def __pow__(v):
        """order must make some sense here
        or we must interpret the input
        """
        #states,v = Goal.process_args(v)
        def pursue(states):
            states_ = []
            for s in states.states:
                sn = Goal.powo(v[0],v[1],v[2])(s)
                for sni  in sn:
                    if sni is not None:
                        states_.append(sni)
            return States(states_)
        return pursue
    """--------------------------------------------"""

    @staticmethod
    def HalfAdder(A,B):
        S = int((A and not B) or (not A and B))
        C = int(A and B)
        return (S, C)

    @staticmethod
    def FullAdder(A,B,C):
        AB = int((A and not B) or (not A and B))
        S = int((C and not AB) or (not C and AB))
        CAB = int(C and AB)
        C1 = int((A and B) or (CAB))
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

        

class HigherGoal(Goal):
    pass

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

##
##*****************************************************
##



class GraphScaleNode(object):
    """A compound relational rule
    compound:  a rule expression 
                composed of at least one ternary logical relational rule
    
    Parameters
    ----------
        Takes a compiled rule: fundict
        and its list of variables: varlist
    """
    def __init__(self, fundict, varlist, name=None):
        self.name = name
        self.fundict = fundict
        self.vars = varlist
        #for var in varlist:
        #    if isa(var,PStates)
        #        var.nodelist.append()
        
    
    
    def __repr__(self):
        return 'GraphScaleNode(\n name: {}\n, rules: {}\n, vars: {})'.format(
                self.name,
                self.fundict,
                self.vars)
        

class RulesGraphProcessor(object):
    """
       Let this class handle the state 
       to ensure consistent assignment of vars 
       and rules composed from the vars
    """
    def __init__(self, nodes=None, env=None,
                 verbose=True):#funclist, varlist):
        if nodes is None:
            self.nodes = []
        else:
            self.nodes = nodes
        
        self.current_rules = copy.copy(self.nodes)
        #self.current_rules = copy.deepcopy(self.nodes)
        self.arc_sets = set([]) #active arcs
        self.var_sets = set([]) #active variables
        self.rules_set = set([])
        self.vars = []#set([])
        self.varsmap = {}
        #for func in funclist:
        #    self.add(func, varlist)
        if env is None: env = States(State())
        self.env = env
        self._verbose = verbose
        
    
    def __call__(self,var):
        return self.get_value(var)
    
    def reset(self, new_env,var):
        self.arc_sets = set([]) #active arcs
        self.var_sets = set([]) #active variables
        self.rules_set = set([])
        self.env = new_env
        self.nodes.pop()
        #self.current_rules.pop() 
        #list is self emptying 
        #(and does so during compute_fresh_rules_graph)
        self.varsmap[var.name].flist.pop()
        #for var in self.vars:
        #    self.env.update({var:new_env(var)})
        return
    
    def check_intervals(self):
        invtlist = []
        for el in self.vars:
            mknm, mkval = self.get_name_and_value(el.name)
            ckv = mkval[0]
            if ckv.inf>ckv.sup:
                invtlist.append(el)
        return invtlist
        
    
    def get_name_and_value(self, var):
        if isa(var,PStates):
            #print 'searching PStates'
            mk_name = self.env.get_by_str(var.name)
            return mk_name, self.env(mk_name)
        elif isa(var,Variable):
            #print 'searching Variable ',var
            mk_name = self.env.get_by_str(var.name)
            #print var is mk_name
            #return mk_name, self.env(mk_name.name)
            return var, self.env(var.name)
        elif isa(var, str):
            #print 'searching String'
            mk_name = self.env.get_by_str(var)
            return mk_name, self.env(mk_name)
        else:
            #print 'var is type: ',type(var)
            return (None, [None])
    
    
    def add_graph_scale_node(self, dflist,vars_,rulename):
        """when adding a rule,
        also add the rule to each 
        mk_var's flist (arc like 'list')
        """
        gsn = GraphScaleNode(dflist,
                             vars_,
                             rulename)
        
        self.current_rules.append(gsn)
        
        self.nodes.append(gsn)
        
        
        for var in vars_:
            mk_var = vars_[var]
            #if isinstance(mk_var, Variable):
            if isa(mk_var, Variable):
                mk_var.flist.append(gsn)
                if mk_var not in self.vars:
                    #mk_var = vars_[var]
                    self.vars.append(mk_var)
                    self.varsmap[var] = mk_var
                
#        for var in vars_:
#            mkvar = vars_[var]
#            if isa(mkvar, Variable):
#                if mkvar not in self.vars:
#                    self.vars.append(vars_[var])
        
        return
        
    def testfunc(self,node):
        print node.name
        return node
        
    def walk_apply(self,node,apply_func):
        """
        Parameters
        ----------
        node        : a node in a computational graph
        apply_func  : function to apply recursively 
                      from root to leaf
        Returns
        ----------
        node        : the starting node
        """
        node = apply_func(node)
        if node.args:
            for nd in node.args:
                self.walk_apply(nd,apply_func)
        else:
            return node
    
    def delete_one_link(self,node):
        """ Recursively cleanse the coputation graph
        held by PStates after the graph has been 
        compiled.
        
        usage
        ----------
            walk_apply_in_reverse(top_node,
                                  apply_func = delete_one_line)
            
        """
        #print '\n deleting link from ',node.name
        #print 'to ',node.args
        node.args = []
        return node
        
    def walk_apply_in_reverse(self,
                              node,apply_func):
        """Walk a computational graph
        Delete links from leaf to root
        -so as not to define the same rule twice!
        
        Parameters
        ----------
        node        : a node in a computational graph
        apply_func  : function to apply recursively 
                      from leaf to root
        Returns
        ----------
        node        : the starting node
        
        >>> self.walk_apply_in_reverse(v1,
                                       delete_one_link)
        """
        if isa(node,Number):
            return
        #elif not node.args:
        elif not hasattr(node, 'args'):
            return
        else:
            for nd in node.args:
                self.walk_apply_in_reverse(nd,apply_func)
            node = apply_func(node)
        return # nothing returned
        
    
    def add_one_rule(self, 
                     rule_node,
                     rulename=None):
        if rulename is None: rulename = rule_node.name
        
        st = self.env
        """-----------------------------------------------
        #construct the environment, {PStates:mkVars} 
        # intermediate data structures
        # using the pre-graphed/parsed rule_node:
        #"""
        st, vars_ = rule_node.construct(rule_node, st, {},
                                        self.vars)
        
        """-----------------------------------------------
        #compile a list of relational rule dictionaries:
        #"""
        dflist = rule_node.compile(rule_node,vars_)
        
        
        """-----------------------------------------------
        #store the updated envronment (it may know about
        #                   new variables, intervals, etc.)
        #"""
        self.env = st
        
        
        """-----------------------------------------------
        # add the rule to each 
        #  mk_var's flist 
        #    (arc-like 'list' which connects from Varaibles
        #                         to relational compound rule
        #                         aka a GraphScaleNode)
        #"""
        if len(dflist)>0:
            self.add_graph_scale_node(dflist,vars_,rulename)
         
        
        """-----------------------------------------------
        #clean the graph (rule_node links, recursive, leaf to root):
        #"""
        self.walk_apply_in_reverse(node=rule_node,
                                apply_func=self.delete_one_link)
        
        return
    
    
    def compute_rules_graph(self):
        """
            TODO: edit this to just compute 
            a short list of new rules.
            
            AC_revise through any additional
            rules that come up
        """
        for rule in self.nodes:
            self.compute_one_rule(rule) #hmmm:
            self.AC_revise() #if adding a rule updates var_sets, 
            #then only AC_revise would be needed.
        return
    
    
    def compute_fresh_rules_graph(self, maxiter=None):
        """
            TODO: edit this to just compute 
            a short list of new rules.
            
            AC_revise through any additional
            rules that come up
        """
        #for rule in self.nodes:
        if maxiter is None:
            while len(self.current_rules) >0:
                rule = self.current_rules.pop()
                self.compute_one_rule(rule) #hmmm:
                self.AC_revise() #if adding a rule updates var_sets, 
                #then only AC_revise would be needed.
        else:
            iter_ = 0
            while len(self.current_rules) >0 and iter_ < maxiter:
                rule = self.current_rules.pop()
                self.compute_one_rule(rule) #hmmm:
                self.AC_revise() #if adding a rule updates var_sets, 
                #then only AC_revise would be needed.
                iter_ += 1
        return
    
    
    
    def compute_one_rule(self, rule, state=None):
        """
            outcome:  rule => state  & var_sets is updated
        """
        if state is None:
            state = self.env
        fundict = rule.fundict #actually a dict
        vars_ = rule.vars
        save_state = copy.deepcopy(state)
        for el in fundict:
            state = self.compute_one_relation(fundict[el], 
                                                  state)
        #if 'unknown_var' in vars_: del(vars_['unknown_var'])
        self.set_updates(state.states, save_state.states, vars_)
        self.env = state
        #print '\nOne rule computation',vars_
        #print 'now var_sets = ',self.var_sets
        #print 'One rule done, env = ',self.env
        #print ''
        
        return
        
    def compute_one_relation(self, rule, state):
        new_state = rule(state)
        return new_state
        
    
        
    def set_updates_old(self, states, states_old, vars_, tol=.001):
        """TODO: make sure you compare the same state from
        old to new.  That way you really do track real
        changes.
        
        -Put a unique identifier on each state 
        
        inputs:  States.states, states_old.states
        (dont be fooled!)
        """
        for state, statei in zip(states, states_old):
            chg_lst = [vars_[el] for el in vars_ if (not state(vars_[el]) == statei(vars_[el]) )]
            #chg_lst = [var for var in vars if not self.diff( ( state(var), statei(var),tol))]
            
            var_set = self.var_sets.union(set(chg_lst)) 
            self.var_sets = var_set
        #print 'set updates returning:'
        #print 'var_sets:',self.var_sets
        #print 'end set updates'
        #if len(self.var_sets)>0:
        #    for el in vars_:
        #        print '{} == {}'.format(state(vars_[el]) , statei(vars_[el]))
        return 
    
    
    def get_arc_sets(self):
        """Find the constraints that have inputs which changed
        
        is it okay to ignore ia(_,^)?  I think the answer is yes
        but it has been a while since I was down in this code...
        """
        #lac = len(self.arc_sets)
        #lvs = len(self.var_sets)
        #for c, var_set in enumerate(self.var_sets):
            #if c>=lac:
            #    self.arc_sets.append(set())
        for var in self.var_sets:
            if isinstance(var, ia):
                pass
            else:
                #self.arc_sets[c] = self.arc_sets[c].union(set(var.arcs))
                self.arc_sets = self.arc_sets.union(set(var.arcs))
    #            else:
    #                for var in var_set:
    #                    self.arc_sets[c] = self.arc_sets[c].union(set(var.arcs))
        return
    
    
    
    
        
    
    
    def get_rules_set(self):
        self.rules_set = set([])
        for var in self.var_sets:
            #            if isinstance(var, ia):
            #                print 'ia = ',var
            #                pass
            #            else:
            mk_var = self.get_name_and_value(var)[0]
            #print 'var = ',var,' mk_var = ', mk_var
            if mk_var is not None:
                self.rules_set = self.rules_set.union(set(mk_var.flist)) 
        return
    
    
    def setup_diff_debug(self, states_new, states_old):
        for state, statei in zip(states_new, states_old):
            pass
        return  state, statei
    
    
    def diff(self, state_new, state_old,el,tol):
        q_new = state_new.get_by_str(el)
        q_old = state_old.get_by_str(el)
        #if svalnew.dist(svalold) < tol:
        nresult = state_new(q_new)
        oresult = state_old(q_old)
        if isa(nresult, Variable) and isa(oresult, Variable):
            return False
        else:
            if state_new(q_new) == state_old(q_old):
                return False
            else:
                return True
        
    def set_updates(self, states, states_old, vars_, tol=.001):
        """
        inputs:  States.states, states_old.states
        
        outcome:  var_sets is updated
        """
        var_set = set([])
        for state, statei in zip(states, states_old):
            #chg_lst = [vars_[el] for el in vars_ if (not state(vars_[el]) == statei(vars_[el]) )]
            chg_lst = [vars_[el] for el in vars_ if self.diff( state, statei, el, tol )]
            
            var_set = self.var_sets.union(set(chg_lst)) 
        self.var_sets = var_set
        #print 'set updates returning:'
        #print 'var_sets:',self.var_sets
        #print 'end set updates'
        return 
    
    
    #def AC_revise_flists(self, print_=False, maxit = 8):
    def AC_revise(self, print_=False, maxit = 8):
        """Arc Consistency for the Hull Parameters
        
        issues:  
            -arcs don't know about rules
                *so they don't know about var_sets to compare on
                *so we are not sure computation is getting 
                -really converged- and 
                certainly no idea if it is done so in efficient
                manner
            -the solution to this was to have the 
            parser-compiler system know about a 'rule'
            as the collection of relations that come 
            in with a PStates var when 'setting a rule'
            (DONE)
            -hence now we 'get_rules_sets'
            instead of 'get_arc_sets'
        """
        #print 'starting flists AC_revise'
        self.get_rules_set()
        #print 'GOT rules set = ',self.rules_set
        
        #print 'reseting var_sets'
        self.var_sets = set([])  
        
        
        
        count=0
        while len(self.rules_set)>0:
            for rule in self.rules_set:
                fundict = rule.fundict #actually a dict
                vars_ = rule.vars
                save_state = copy.deepcopy(self.env)
                for key in fundict:
                    self.env = self.compute_one_relation(fundict[key], 
                                                         self.env)
                self.set_updates(self.env.states, save_state.states, vars_)
                #set updates: => var_sets
                    
            #self.var_sets = set([])
            #save_state = copy.deepcopy(self.env)
            self.rules_set = set([])
            self.get_rules_set() #replace get arc_sets 
                                 #because 
                                 # rules are several arc connected together
                                 # and this 'solver' uses them together.
            #
            self.var_sets = set([])
            
            count += 1
            #if print_ or self._verbose: print count
            if print_: print count
            
            
            """#this is a huge tiny line of code!
            # performance vs exactitude
            #"""
            #if count<maxit:
            #    break
            """# but interval width comparisons
            # may be even more of a factor.
            #"""
            
        return
    
    
    
    
    #def AC_revise(self, print_=False, maxit = 8):
    def AC_revise_old(self, print_=False, maxit = 8):
        """Arc Consistency for the Hull Parameters
        
        issues:  
            -arcs don't know about rules
                *so they don't know about var_sets to compare on
                *so we are not sure computation is getting 
                -really converged- and 
                certainly no idea if it is done so in efficient
                manner
            -?
        """
        self.get_arc_sets()
        
        #        if self._verbose:
        #            #print 'arcs of interest:\n'
        #            for el in self.arc_sets:
        #                print el
        
        #
        #self.var_sets = set([])  #with this gone,
        #
        # set_updates 
        # really gets going, whether or not it is efficient
        # is another matter
        #
        
        #for c, arc_set in enumerate(self.arc_sets):
        count=0
        while len(self.arc_sets)>0:
            #for arc in self.arc_sets:
            while len(self.arc_sets)>0:
                arc = self.arc_sets.pop()
                #print len(self.arc_sets),self.var_sets
                #if print_ or self._verbose: print '\ndoing arc', arc
                save_state = copy.deepcopy(self.env)
                self.env = arc(self.env)
                self.get_arc_sets()
                # compute one rule?
                #
                self.set_updates(self.env.states, save_state.states, 
                                 { key: value for (key, value) in enumerate(self.var_sets) }
                                  )
                #self.var_sets = set([]) 
                #if (self.env.ismepty) #TODO: add this? - nah, check where you really need it
            self.arc_sets = set([])
            self.get_arc_sets()
            self.var_sets = set([])
            #            self.arc_sets[c] = set([])
            #            self.get_arc_sets()
            #            self.var_sets[c] = set([])
            count += 1
            if print_ or self._verbose: print 'ac:',count
            
            #if count<maxit:
            #    break
            
        return
    
    def print_state(self):
        for var in self.vars:
            print var, ' = ', self.env(var)



##
##*****************************************************
##



def wisdom_byrd():
    print 'define my-append:'
    print 'takes 2 args'
    print ''
    print 'in general when checking lists'
    print 'the first test should be to'
    print 'check if the list is empty'
    print 'do this before checking list components, always!'
    print 'hangout 1 is here:'
    print 'https://www.youtube.com/watch?v=a5p8DPbaokE'
    print ''
    print 'do we check 1st or 2nd for null-ness first?'
    print 'probably better to check the first list'
    print ''
    print 'no recurse'
    print 'since we are checking the 1st element then we need to contract it'
    print 'to eventually terminate '
    print ''
    print 'to write a recursive definition'
    print 'play make believe'
    print 'https://www.youtube.com/watch?v=a5p8DPbaokE'
    print '1:40:00'
    return
    
def wisdom_presentation():
    print 'find important stuff here:'
    print 'hull rules : simple_hull_rules'
    return

if __name__ == '__main__':
    world = State()

    t1  = False# True#
    t2  = False# True#
    t3  = False# True#
    t4  = False# True#
    t5  = False# True#
    t6  = True# False# 
    t7  = True# False# True#
    t8  = True# False# True#
    t9  = True# False# True#
    t10 = True# False# True#
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
    print 'x:',to_list(s[0].value_of(x))

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

    #if i3:
    if False:
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
        #g = Goal.addo(x,y,c) #does not work with binop on lists
        #s = g(s)[0]

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
        #s1 = (s - (y,x,c))#g(s)[0] #does not work with binop on lists!
        a = ia(1.0,1.5)
        b = ia(0.5,1.0)
        c1 = b-a #should trim contraint when done in s!
        #c1 & ia(0.,1.9) = (0.,0.)
        print '\nnaive c = {}'.format(c1)
        print 'minimal relational c = {}\n'.format(s1[0].value_of(c))

        #s2 = (s + (x,y,c))#g(s)[0]#does not work with binop on lists!
        a = ia(1.0,1.5)
        b = ia(0.5,1.0)
        c1 = a+b #should trim contraint when done in s!
        #print '\nnaive c = {}'.format(c1)
        #print 'minimal relational c = {}\n'.format(s2[0].value_of(c))

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

        #s3 = ((s + (y,x,c))[0] + (x,y,d))#g(s)[0]
        a = ia(1.0,1.5)
        b = ia(0.5,1.0)
        c1 = b+a #should trim contraint when done in s!
        #c1 & ia(0.,1.9) = (0.,0.)
        print '\nnaive c = {}'.format(c1)
        #print 'minimal relational c = {}\n'.format(s3[0].value_of(c))
        #print 'minimal relational d = {}\n'.format(s3[0].value_of(d))


    if True: #States():

        x = Variable('x')
        y = Variable('y')
        c = Variable('c')
        d = Variable('d')
        s = State(values={x:None,y:None,c:None,d:None})
        g = Goal.eq(x,ia(1.0,1.5))
        s = g(s)[0]
        g = Goal.eq(y,ia(0.5,1.))
        s = g(s)[0]
        g = Goal.eq(c,ia(-1.0,1.9))
        stest = g(s)[0]
        stest = States(stest)

        x = Variable('x')
        y = Variable('y')
        c = Variable('c')
        d = Variable('d')
        s = State(values={x:None,y:None,c:None,d:None})
        st = States(s)
        st = (st == (x,ia(1.0,1.5) ))
        st = (st == (y,ia(0.5,1.) ))
        st = (st == (c,ia(-1.0,1.9)) )
        
        print 'integrate inequality'
        st1_test1 = ( (st + (x,y,c) ) <= (y,d))
        st1_test2 = ( (st - (x,y,c) ) <= (y,d))

        st0 = (st / (x,y,d))
        st0 = (st / (x,d,y))
        st0 = (st / (d,x,y))

        #f0 = (st + (x,y,c) )
        st1 = (st + (x,y,c) )#(st.states)  c=ia(1.5,1.9)
        a = ia(1.0,1.5)
        b = ia(0.5,1.0)
        c1 = a+b #should trim constraint when done in s!
        print '\nnaive c = {}'.format(c1)
        print 'minimal relational c = {}\n'.format(st1.states[0].value_of(c))


        #st1 = ( (st + (x,y,c) )(st.states) + (x,y,d)(st.states) )
        st1 = ( (st + (x,y,c) ) + (x,y,d))
        #f1 = (st + (x,y,c))
        #f2 = (st + (x,y,d))
        #st1 = st.both( f1, f2 )#(st.states)
        #f2 = (f1 + (x,y,d))


        #st2 = (st + (x,y,c))(st.states)
        #st2 = (st2 + (x,y,d))(st2.states)
        #st2 = (st + (x,y,c))
        #st2 = (st2 + (x,y,d))
        st2 = (st + (x,y,c) + (x,y,d))

        print ''
        print 'extended division:'
        print ' ia(0.5,1.0)/ia(-1.0,1.9) ={}\n\n'.format(
                ia(0.5,1.0)/ia(-1.0,1.9) )

        st3 = (st * (c,y,d)) #d = ia(-1.0, 1.9)
        print 'st3'
        print st3
        st4 = (st / (y,c,d)) #splits in two!  d=ia(0.263157894737,1e+25), d= ia(-1e+25,-0.5)
        print 'st4'
        print st4

        st5 = (st / (d,y,c) )

        #        st6 = States.both((st + (x,y,c) ),
        #                          (st + (x,y,d)))

        x = Variable('x')
        y = Variable('y')
        c = Variable('c')
        d = Variable('d')
        s = State(values={x:None,y:None,c:None,d:None})
        st6 = States(s)
        #st7 = st6.eq((d,c))
        #st7 = st7.eq( (ia(.9,1.5), x) )
        #st7 = st7.eq( (ia(.5,1.2), y) )

        st7 = (st6 == (d,c))
        #st7 = (st6 == (ia(.9,1.5), x) )
        st7 = (st7 == (ia(.9,1.5), x) )
        st7 = (st7 == (ia(.5,1.2), y) )
        st7 = (st7 == (ia(-1.,1.9),c) )
        print '\nst7\n',st7
        st8 = (st7 / (x,y,c) )
        
        
#        st7p = (st6 == (ia(.9,1.5), x) )
#        st7p = (st7p == (ia(.5,1.2), z) )
#        st7p = (st7p == (ia(-1.,1.9),c) )        
#        st8p = (st7p / (x,y,z))
        
        vol = Variable('vol')
        c2 = Variable('c2')
        Cb = Variable('Cb')
        st = State(values={vol:None,c2:None,Cb:None})
        st = States(st)
        st = (st == ( ia(0.0, 1.0), Cb) )
        st = (st == ( ia(2000.0, 300000.0),vol ) )
        #st / (Cb,c2,vol)
        
        st / (vol,c2,Cb) #3
        st * (Cb,c2,vol) #3
        
        st / (vol,Cb,c2) #2
        st * (c2,Cb,vol) #2
        
        gg = ia(2000.0, 300000.0) / ia(0.0, 1.0)
        st == (c2,gg[0])

        self = st7
        """
        s=self.states[0]
        v=(x,y,c)
        x,y,z=v[0],v[1],v[2]
        state = self.states[0]
        goal = Goal.add
        x1 = state.value_of(x)
        y1 = state.value_of(y)
        z1 = state.value_of(z)
        #"""
        print 'st9 is a list of states resulting from interval split'
        st9 = (st7 / (y,c,d)  )
        st9 = (st9 / (y,c,d)  )#now this
        #st9 = (st9 == (d,c)) #or this is needed
        print '\nst9\n',st9


        c1 = Variable('c1')
        draft = Variable('draft')
        c2 = Variable('c2')
        s = State(values={c1:None,draft:None,c2:None})
        ht = States(s)
        ht = (ht == ( c1, ia(2400.0,2400.0) ) )
        ht = (ht == ( c2, ia(0.0, 300000.0)  ) )
        ht = (ht * (c1,draft,c2) )  #draft=ia(0.0,125.0)

        c1 = Variable('c1')
        draft = Variable('draft')
        c2 = Variable('c2')
        s = State(values={c1:None,draft:None,c2:None})
        h1 = States(s)
        h1 = (h1 == ( c1, ia(2400.0,2400.0) ) )
        h1 = (h1 == ( c2, ia(0.0, 300000.0)  ) )
        h1 = (h1 == ( draft, ia(-20.,20.)  ) )
        h1 = (h1 * (c1,draft,c2) ) #
        """
        v=(c1,draft,c2)
        s = h1.states[0]
        state = s
        x,y,z = v[0],v[1],v[2]
        goal = Goal.mul
        def func_to_see_extended_divo_issues():
            self = st

            sq = run_all_binary_relations(ia(0.5,1.0),
                                        ia(-1.0,1.9),
                                        d,
                                        self.states[0],Goal.extended_division)

            g = Goal.extended_division(ia(0.5,1.0),
                                    ia(-1.0,1.9),
                                    d)
            g1 = g(s)
            v=(y,c,d)
            states = []
            for s in self.states:
                sn = Goal.extended_divo(v[0],v[1],v[2])(s)
                for sni  in sn:
                    states.append(sni)
            return
        #"""
            
            
        # what about:
        """
            if we could do 
            #stL = (st6 == (x,y),(d,c))   # Did not work 'right'!
            stL = (st6 == ((x,y),(d,c)) )  # Does work right!  (tuple assignment style)
            stL = (stL == (d,2))          # d=2, x=2
            
            #stL = (stL == (x,y),(d,c))
        """
        #stL = (st6 == (x,y),(d,c))
        #stL = (stL == (x,c))
        
        
        x = Variable('x')
        y = Variable('y')
        z = Variable('z')
        a = Variable('a')
        b = Variable('b')
        c = Variable('c')
        d = Variable('d')
        s = State(values={x:None,y:None,z:None,
                          a:None,b:None,c:None,d:None})
        sto = States(s)
        
        stL = (sto == (a,1))
        stL = (stL == (b,1))
        stL = (stL == (x,3))
        
        #stL = (stL + (x,y,z) ) #have to do this twice to get y=-3
        #stL = (stL + (z,a,b) ) #sets z=0
        
        stL = ((stL + (z,a,b)  )  + (x,y,z) )  #assigns y=-3 in on go
        stL = ((stL + (x,y,z)  )  + (z,a,b) )  #requires two passes
        
        ###stL = ( (stL + (x,y,z)  )  + (z,a,b) )  # don't know how to do higher orderconstraints in one pass
        
        
        x = Variable('x')
        y = Variable('y')
        c = Variable('c')
        d = Variable('d')
        s = State(values={x:None,y:None,c:None,d:None})
        st10 = States(s)
        ss = ( (st10 == (x,ia(.1,10.))) == (ia(5.,7.),y))
        ss1 = (ss >= (x,y)) 
        

        def evali(a,b):
            if isinstance(a, float):
                print a, 'is float'
            elif isinstance(a, Variable):
                print a, 'is Variable'
            elif callable(a):
                print a,'is callable'
                return a(b)
                
    
    #import kanren
    #from kanren import run, var, eq
    #from kanren import lall
    #x = kanren.var()
    #x = Variable()
    #kanren.run(1, x, eq(x, 1))
    x = Variable('x')
    y = Variable('y')
    s = State({x:None,y:None})
    st = States(s)
    st = (st == (x,ia(1.,2.)))
    st = (st == (y,ia(3.,4.)))
    
    def r1(x,y,st):
        st = (st + (x,y,ia(6.,6.)) )
        return st
        
    st = r1(x,y,st)
    
    """
    test = (st == (y,ia(2.,3.)))
    len(test.states)  
    >>> 0
    #"""