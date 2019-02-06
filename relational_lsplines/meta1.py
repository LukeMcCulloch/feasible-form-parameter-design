#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 16:53:55 2017

@author: luke
"""
#import numpy as np
from extended_interval_arithmetic import ia
#from hull_inference_ob_graph import Hull as hullclp

#import uKanren as lp #original
#import eKanren as lp#nicer! NOT backwards compatible
import sqKanren as lp

import math
import operator as op
import numpy as np

def parall():
    return
def madd(mlist):
    return

Symbol = str          # A Scheme Symbol is implemented as a Python str
List   = list         # A Scheme List is implemented as a Python list
Number = (int, float) # A Scheme Number is implemented as a Python int or float

def tokenize(chars):
    "Convert a string of characters into a list of tokens."
    return chars.replace('(', ' ( ').replace(')', ' ) ').split()


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
        
def parse(program):
    "Read a Scheme expression from a string."
    return read_from_tokens(tokenize(program))


#Env = dict          # An environment is a mapping of {variable: value}
class Env(dict):
    "An environment: a dict of {'var':val} pairs, with an outer Env."
    def __init__(self, parms=(), args=(), outer=None):
        self.update(zip(parms, args))
        self.outer = outer
    def find(self, var):
        "Find the innermost Env where var appears."
        return self if (var in self) else self.outer.find(var)



def standard_env():
    "An environment with some Scheme standard procedures."
    env = Env()
    #env.update(vars(math)) # sin, cos, sqrt, pi, ...
    env.update({
        '+':op.add, '-':op.sub, '*':op.mul, '/':op.div, 
        '>':op.gt, '<':op.lt, '>=':op.ge, '<=':op.le, '=':op.eq, 
        'abs':     abs,
        'append':  op.add,  
        'apply':   apply,
        'begin':   lambda *x: x[-1],
        'car':     lambda x: x[0],
        'cdr':     lambda x: x[1:], 
        'cons':    lambda x,y: [x] + y,
        'eq?':     op.is_, 
        'equal?':  op.eq, 
        'length':  len, 
        'list':    lambda *x: list(x), 
        'list?':   lambda x: isinstance(x,list), 
        'map':     map,
        'max':     max,
        'min':     min,
        'not':     op.not_,
        'null?':   lambda x: x == [], 
        'number?': lambda x: isinstance(x, Number),   
        'procedure?': callable,
        'round':   round,
        'symbol?': lambda x: isinstance(x, Symbol),
    })
    return env

global_env = standard_env()


#def eval(x, env=global_env):
#    "Evaluate an expression in an environment."
#    if isinstance(x, Symbol):      # variable reference
#        return env[x]
#    elif not isinstance(x, List):  # constant literal
#        return x                
#    elif x[0] == 'if':             # conditional
#        (_, test, conseq, alt) = x
#        exp = (conseq if eval(test, env) else alt)
#        return eval(exp, env)
#    elif x[0] == 'define':         # definition
#        (_, var, exp) = x
#        env[var] = eval(exp, env)
#    else:                          # procedure call
#        proc = eval(x[0], env)
#        args = [eval(arg, env) for arg in x[1:]]
#        return proc(*args)
def eval(x, env=global_env):
    "Evaluate an expression in an environment."
    if isinstance(x, Symbol):      # variable reference
        return env.find(x)[x]
    elif not isinstance(x, List):  # constant literal
        return x                
    elif x[0] == 'quote':          # quotation
        (_, exp) = x
        return exp
    elif x[0] == 'if':             # conditional
        (_, test, conseq, alt) = x
        exp = (conseq if eval(test, env) else alt)
        return eval(exp, env)
    elif x[0] == 'define':         # definition
        (_, var, exp) = x
        env[var] = eval(exp, env)
    elif x[0] == 'set!':           # assignment
        (_, var, exp) = x
        env.find(var)[var] = eval(exp, env)
    elif x[0] == 'lambda':         # procedure
        (_, parms, body) = x
        return Procedure(parms, body, env)
    else:                          # procedure call
        proc = eval(x[0], env)
        args = [eval(arg, env) for arg in x[1:]]
        return proc(*args)

def repl(prompt='lis.py> '):
    "A prompt-read-eval-print loop."
    while True:
        val = eval(parse(raw_input(prompt)))
        if val is not None: 
            print(schemestr(val))

def schemestr(exp):
    "Convert a Python object back into a Scheme-readable string."
    if isinstance(exp, List):
        return '(' + ' '.join(map(schemestr, exp)) + ')' 
    else:
        return str(exp)
    
    
    
class Procedure(object):
    "A user-defined Scheme procedure."
    def __init__(self, parms, body, env):
        self.parms, self.body, self.env = parms, body, env
    def __call__(self, *args): 
        return eval(self.body, Env(self.parms, args, self.env))




if __name__ == "__main__":
    
    x = lp.Variable('x')
    y = lp.Variable('y')
    z = lp.Variable('z')
    a = lp.Variable('a')
    b = lp.Variable('b')
    c = lp.Variable('c')
    d = lp.Variable('d')
    s = lp.State(values={x:None,y:None,z:None,
                      a:None,b:None,c:None,d:None})
    st = lp.States(s)
    
    h = parse('(states + (a,b,c) )')