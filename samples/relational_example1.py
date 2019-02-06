#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 00:18:46 2018

@author: luke mcculloch


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

print "\nI think it's fun to look at the computation tree"
print "self assembled in these rules.\n"
print 'For example, the computation tree for this rule:\n'
print '   Cb = Cb == vol/(lwl*bwl*draft)'
print "\nlooks like this:\n"
print Cb
print 'As you can see, the contents are a tree.'
print "Cb is what I call the 'top node' of this tree"
print "It's also the node which is kind of at the "
print "recieving end of the rule."



#
#***********************************
# 
print "\nAdd the above rules to the Database like so: \n"
print 'hullrulesnet.set_rules(Cb)'
hullrulesnet.set_rules(Cb,disp)

print 'There is nothing stoping us from constructing several rules and'
print 'Adding them to the rules network in one go.'
print 'Except one critical feature!'
print 'You cannot safely add a rule featuring the _same_ PStates variables'
print 'in two different rules at the same time.'
print ''
print "This is because each variable involved in a rule "
print "keeps track of it's children in the rules graph."
print "It holds them in an attribute called 'args'" 
print ""
print "It could be the case that a particular PStates variable"
print "already has it's args list occupied"
print "and thus if we were to use it in another rule,"
print "this would overwrite the args list without ever setting "
print "the first rule!"

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