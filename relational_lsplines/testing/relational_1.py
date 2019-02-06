#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May 19 13:19:08 2018

@author: luke

Current rule operations in production
    ==      :: two atoms should have the same value    
    +       :: addition
    -       :: subtraction
    *       :: multiplication
    /       :: division
    <=      :: less than or equal to
"""

import relational_lsplines as rlspline 
#from relational_lsplines import simple_hull_rules_language as srl

from relational_lsplines.simple_hull_rules_language import HullGeometryGenerator


lp = rlspline.lp #relational logic programming
ia = rlspline.ia #interval analysis


hullrulesnet = HullGeometryGenerator(rtype='gauss',
                                     verbose=True) 



"""Make some design space variables:
"""
lwl = lp.PStates(name='lwl')
bwl = lp.PStates(name='bwl')
draft = lp.PStates(name='draft')
vol = lp.PStates(name='vol')
disp = lp.PStates(name='disp')

Cb = lp.PStates(name='Cb')
Cp = lp.PStates(name='Cp')


#x = lp.PStates(name='x')
#y = lp.PStates(name='y')
#z = lp.PStates(name='z')
#a = lp.PStates(name='a')
#b = lp.PStates(name='b')

#
#***********************************
# Basic Dimension Relations
"""
#rule: rename vol to disp
"""
disp = disp == vol

"""-----------------------------------------------
#rule: block coefficient
"""
Cb = Cb == vol/(lwl*bwl*draft)

print Cb



#
#***********************************
# 

hullrulesnet.set_rules(Cb,disp)


hullrulesnet.rgp.compute_fresh_rules_graph()




Cb = Cb == ia(-1.,1.)


lwl = lwl == ia(80.,120.)

vol= vol == ia(101000.,202000.)

draft = draft == ia(-20.,40.)


hullrulesnet.set_rules(Cb, 
                      lwl,
                      vol,
                      draft)
print hullrulesnet.rgp.env



#
#*****************************************************************************
#
Cb = Cb == ia(.7,.85)
hullrulesnet.set_rules(Cb)
print hullrulesnet.rgp.env
print '------------------'
#
#*****************************************************************************
#
draft = draft == ia(20.,40.)
hullrulesnet.set_rules(draft)
print hullrulesnet.rgp.env
print '------------------'