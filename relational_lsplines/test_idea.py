#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 13:41:20 2017

@author: luke

    Testing examples from 
    
    Interval Constraint Logic Programming
    
    by Frederic Benhamou

"""
#import numpy as np
from extended_interval_arithmetic import ia
#from hull_inference_ob_graph import Hull as hullclp

#import uKanren as lp #original
#import eKanren as lp#nicer! NOT backwards compatible
import sqKanren as lp


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
    
    #st1 = (st == (y,ia(1.,3.)))
    #st1 = (st == (ia(-100.,100.),x))
    #st1 = (st1 * (x,ia(2.,2.),y))
    
    #st1 = (st1 ** (x,2.,y))
    #st1 = (st1 ** (y,.5,x))
    #st1 = (st1 ** (y,-1.,x))
    
    # finding exponenets must have positive y
    #st1 = (st == (y,ia(1.,3.)))
    #st1 = (st1 ** (ia(1.73205080757,1.73205080757) , x, y ) )
    
    
    st1 = (st == (y,ia(-1.,3.)))
    #st1 = (st1 ** (y,-2.,x))
    st1 = (st1 ** (y,-1.,x))
    
    
    
    
    print 'Query: y in {},y = x**2'.format(st1(y))
    print 'answer:'
    print 'y = {}'.format(st1(y)) 
    print 'x = {}'.format(st1(x)) 
    
    print 'NOTE:  this does not support either '
    print 'integer reasoning of CLP(BNR) Prolog'
    print 'nor Boolean reasoning of CLP(BNR)'
    print 'at the moment'
    
    
    
    #E=0.
    #st2 = (st + (x,1.,a))
    #st2 = (st2 ** (x,19.,b))
    #st2 = (st2 * (E,b,c))
    #st2 = (st2 + (a,c,d))
    
    
    x = lp.Variable('x')
    y = lp.Variable('y')
    c = lp.Variable('c')
    d = lp.Variable('d')
    s = lp.State(values={x:None,y:None,c:None,d:None})
    st6 = lp.States(s)
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
    
    print 'st9 is a list of states resulting from interval split'
    st9 = (st7 / (y,c,d)  )
    
    
    
    
    ##
    ##
    ##
    lp.reify(x,{x:'Hello Reify'})
    
    print 'reification allows us to compute things and return values'
    st10 = (st6  * ( lp.reify(x,{x:ia(1.,2.)}), y, ia(5.,10) )   ) 
    
    dummy = lp.Variable()
    #lp.reify(dummy,{dummy:st10(x)/st1})
    st10 = (st10  * (x,y,dummy)  )
    
    print '\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n'
    """
    ###
    #############
    ###
     Extended Division
    ###
    #############
    ###
    """
    X = ia(2.,4.)
    Y = ia(-2.,2.)
    Z = ia(-10.,10)
    def printst():
        print 'x = ', state1(x)
        print 'y = ', state1(y)
        print 'z = ', state1(z)
        return
        
        
    x = lp.Variable('x')
    y = lp.Variable('y')
    c = lp.Variable('z')
    s = lp.State(values={x:None,y:None,z:None})
    state1 = lp.States(s)
    
    print '\n\n----------------------------------------------------'
    print '----------------------------------------------------'
    print '\n'
    print 'Start the example for the talk'
    print '\n'
    #st7 = (st6 == (ia(.9,1.5), x) )
    state1 = (state1 == (ia(2.,4.), x) )
    state1 = (state1 == (ia(-2.,2.), y) )
    state1 = (state1 == (ia(-10.,10),z) )
    
    print '\n\n----------------------------------------------------'
    print '----------------------------------------------------'
    print 'initial x y z values:\n'
    printst()
    print '\n y contains 0!'
    print ''
    #print '\nstate1\n',state1
    print 'lets compute  '
    print '(state1 / (x,y,z))'
    print '\n In words, x over y is z.'
    print ''
    state1 = (state1 / (x,y,z) )
    print '----------------------------------------------------'
    print '----------------------------------------------------'
    
    