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
import copy


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
    state1 = (state1 == (ia(-10.,10.),z) )
    
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
    
    print '\n\n----------------------------------------------------'
    print '----------------------------------------------------'
    print 'result:\n'
    printst()
    
    print '\n###\n#############\n###'
    print '\nAll states are immediately non-singular'
    print 'what happened? \n'
    
    print 'To see what happened,\n'
    print ''
    print '-------------------------'
    print 'FIRST:  Look at x/y:'
    print ''
    print 'x/y -> ia(2.,4.)/ia(-2.,2.) = ?\n'
    print 'splitting the interval'
    print  '[ia(-4e+25, -1.0), ia(1.0, 4e+25)]'
    print ''
    print 'z = z & (x/y)'
    print '  = ia(-10. , 10.) &  [ia(-4e+25, -1.0), ia(1.0, 4e+25)]'
    print '  = [ia(-10., -1.) , ia(1., 10.)]'
    
    print ''
    print '\nBUT the relational interpreter can '
    print ' use this for further narrowing...'
    print 'y = y & x /z\n'
    print ''#'\nSO'
    print ''
    print '\n\n-------------------------'
    print 'Second:'
    print '\nSo what is (x/z)?'
    print '(x/z) = ia(2.0, 4.0) * [ia(-10., -1.) , ia(1., 10.)] '
    print '(x/z) = [ia(-4.0, -0.2) , ia(0.2, 4.0) ] \n\n'
    
    """Get there quick with
    X / (Z & (X/Y)[0])
    
    X / (Z & (X/Y)[1])
    
    """
                  
    print 'y = ia(-2., 2.) & [ia(-4.0, -0.2) , ia(0.2, 4.0) ]  '
    
    """Get there quick with
    Y & ( X / (Z & (X/Y)[0]))
    
    Y &  ( X / (Z & (X/Y)[1]))
    """
    print 'y = [ ia(-2.0, -0.2) , ia(0.2, 2.0) ]  '
    print ''
    print ''
    print 'That does it so far.  we now have 2 disjoint intervals for z and y'
    print ''
    print 'What happens if we run the original relation again on this state?'
    print ''
    print 'computing \n state1 = (state1 / (x,y,z) ) \n again...\n'
    state1 = (state1 / (x,y,z) )
    printst()
    print '\nwhat happens if we apply it yet again?\n'
    print 'computing \n state1 = (state1 / (x,y,z) ) \n again...\n'
    state1 = (state1 / (x,y,z) )
    printst()
    print '\n\n Ah, finally there was no change.  We have reached'
    print 'the maximum we can contract with the given information'
    
    
    
    #print st
    #print st.eval((x,st.__eq__,(x,ia(5.,6.))))
    
    sv = copy.copy(st)
    #st  = (st == (x,ia(5.,6.)) )
    
    #set y, return Varaible x:
    a,b = st.eval( [x, '==', (y, ia(1.0, 2.0))]  )
    v = ( (x,'==',(y,ia(1.,2.))) , ia(2.,4.) ) #package -not for anything in particular?
    
    #
    si1 = (st == (ia(2.,4.), x) )
    
    #
    si2 = (st == (x,ia(2.,4.)) )
    
    #set both:
    si3 = (st == [ (x,'==',[y,ia(1.,2.)]), ia(2.,4.)] 
          )
    
    #st = copy.copy(sv)
    
    sp = (st + [ (x,'==',[(y,'==',[c,ia(0.3,2.2)]),ia(1.,2.)]), y, ia(2.,4.)] 
          )
    #v = [ (x,'==',[y,ia(1.,2.)]), y, ia(2.,4.)]
    #st,b = st.eval(v[0])
    #sp = (st + (x,y,ia(2.,4.))
    #        )
    
    #st = copy.copy(sv)
    
    sp2 = (st == (y,ia(1.,2.)) + (x,y,ia(2.,4.)) )
    sp3 = (st + (x,y,ia(2.,4.)) == (y,ia(1.,2.)) )
    sp4 = (st == (y,ia(1.,2.)) )
    sp4 = (sp4 + (x,y,ia(2.,4.)) )
    
    """TODO?: make pre-class """
    sm = (st - [ (x,'==',[y,ia(1.,2.)]), y, ia(2.,4.)] 
          )
    """
        The pre-class implements all the methods of the
        States class: ==, +, -, *, / **
        but instead of doing the rules,
        it stores the computational graph for them
    """
    
    #smm = st + (x,y,z,a) ## x+y+z = a
    _ = lp.Variable('_')
    #st.update({_:None})
    #st.bind(_)
    st.bind_if_not_extant(_)
    #st.states[0].update({_:None})
    #st.states[0] = st.states[0].bind('_')
    #smm = (st + ( ( _,'+',[x,y,_] ), z, a) )
    
        
    #    import inspect
    #    
    #    def fi(st):
    #        spec = inspect.getargs(st)
    #        return spec,st
    #    
    #    spec,sc = fi([y,ia(1.,2.)]), y, ia(2.,4.))
            
    ps = lp.PStates(name = 'ps')
    a = lp.PStates(name='a')
    b = lp.PStates(name='b')
    c = lp.PStates(name='c')
    d = lp.PStates(name='d')
    #ps1 = (ps == (a,b) ) #does not work this way anymore...
    ps1 = a == b
    
    #e =  ( ps + (a,b,c) )
    e = a+b+c
    e.name = 'e'
    
    f = a+b+c*d
    f.name = 'f'
    
    
    ps1 = ( ps == (d,(ps + (a,b,c)) ) )
    #mK language : in some state, ==  (d(+ (+ (a,b,_) c) )
    
    a = a == ( ia(4.,5.) ) #todo->DONE : make this work?
    
    d = a+b+c #careful, this will impact f, above!
    d.name = 'd'
    
    ps1.name = 'ps1'
    
    ps2 = ( ps == (a,(ps - (b,c))))
    
    
    """
    st =  States(State())
    tree = d
    cur_el = st.bind(tree.name)
    
    st, vars_ = d.construct(d, st, {})
    """
    print 'construct state for expression tree d:'
    st, vars_ = d.construct(d, lp.States(lp.State({})), {})
    #ss,vv = d.construct(d, st, vars_)
    #st, vars_ = a.construct(a+b+c, lp.States(lp.State({})), {})
    funclist = d.compile(d,vars_)
    #a = a == ( ia(4.,5.) )
    
    a1 = a == ( ia(4.,5.) )
    a2 = a.equal( ia(4.,5.) )
    #b = b == ()
    #st, vars_ = d.construct()
    
    print a1
    print f
    print 'NOTE THAT a and f represent different rules!'
    print 'maybe you do not have to make f smart about a'
    
    
    