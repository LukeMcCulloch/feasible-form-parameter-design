#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 18:55:00 2017

@author: luke
"""
import numpy as np
import math
import operator as op
from extended_interval_arithmetic import ia
import sqKanren as lp

Symbol = str          # A Scheme Symbol is implemented as a Python str
List   = list         # A Scheme List is implemented as a Python list
Number = (int, float) # A Scheme Number is implemented as a Python int or float

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
        
#Env = dict          # An environment is a mapping of {variable: value}
class Env(dict):
    "An environment: a dict of {'var':val} pairs, with an outer Env."
    def __init__(self, parms=(), args=(), outer=None):
        self.update(zip(parms, args))
        self.outer = outer
    def find(self, var):
        "return the innermost Env where var appears."
        return self if (var in self) else self.outer.find(var)
    def __call__(self, var):
        return self[var]
    
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
        '=='    :  lp.Goal.eq,
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
        'pi'     : np.pi
    })
    return env

global_env = standard_env()


def eval(x, env=global_env):
    "Evaluate an expression in an environment."
    if isinstance(x, Symbol):      # variable reference
        if x in env:
            return env[x]
        else:
            return x
    elif not isinstance(x, List):  # constant literal
        print 'constant literal ',x
        return x                
    elif x[0] == 'if':             # conditional
        (_, test, conseq, alt) = x
        exp = (conseq if eval(test, env) else alt)
        return eval(exp, env)
    elif x[0] == 'define':         # definition
        (_, var, exp) = x
        env[var] = eval(exp, env)
    elif isinstance(x[0],lp.States):
        proc = eval(x[1], env)
    else:                          # procedure call
        proc = eval(x[0], env)
        args = [eval(arg, env) for arg in x[1:]]
        return proc(*args)
    
def eparse(it):
    return eval(parse(it))

if __name__ == '__main__':
        
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
    
    h = parse('(st + (a,b,c) )')
    
    eval(parse('(define r 10)'))
    
    eval(parse('(define x y)'))
    print eval(parse('r'))
    
    print eval(parse(('(if (> (* 11 11) 120) (* 7 6) oops)')))
    print eval(parse(('(if True (* 7 6) oops)')))
    print eval(parse(('(if (equal? 120 120) (* 7 6) oops)')))
    #eval(parse('(define circle-area (lambda (r) (* pi (* r r)))'))
    #print eval(parse((circle-area 10)))
    
    print eval(parse('x'))
    print eval(parse('y'))