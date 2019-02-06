#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 00:18:46 2018

@author: luke
"""
print '--------------------------------------------------------------\n'
print ''
print '     This file corresponds to the Jupyter notebook'
print '         quick_hull.ipynb'
print '     Run that notebook - or view it without running(?)'
print '         to see the full documentation'
print ''
print 'The idea is to show the quickest way'
print 'to use relational_lsplines'
print 'to generate some Feasible Form Parameter Design Ship Hull Shapes'
print ''
print '--------------------------------------------------------------\n'
    
import relational_lsplines as rlspline 

#
#*******************************************************************************
# search tool for bare hull design space is encapsulated here
from relational_lsplines.simple_hull_rules_language import HullGeometryGenerator

#*******************************************************************************
# class which wraps "everything" and funcitons as the real API
# (we could get to the search tool above through this)
from relational_lsplines.opt_simple_hull import ShipDesigner

lp = rlspline.lp #logic programming with intervals (sqKanren)
ia = rlspline.ia #extended interval arithmetic (sometimes divide by zero and still contract the space!)



def main():
    """make a hull
    """
    hullrulesnet = HullGeometryGenerator(rtype='gauss',
                                         verbose=True) 
    
    
    
    print '\n--------------------------------------------------------------'
    print ' Defining Variables for the design space '
    print '\n--------------------------------------------------------------'
    lwl = lp.PStates(name='lwl')
    bwl = lp.PStates(name='bwl')
    draft = lp.PStates(name='draft')
    vol = lp.PStates(name='vol')
    disp = lp.PStates(name='disp')
    
    Cb = lp.PStates(name='Cb')
    Cp = lp.PStates(name='Cp')
    
    Awp = lp.PStates(name='Awp')
    Cwp = lp.PStates(name='Cwp')
    
    Acp = lp.PStates(name='Acp')
    Ccp = lp.PStates(name='Ccp')
    
    Amsh = lp.PStates(name='Amsh')
    Cmidshp = lp.PStates(name='Cmidshp')
    
    LCG = lp.PStates(name='LCG')
    Clcg = lp.PStates(name='Clcg')
    
    Cbl = lp.PStates(name='Cbl') #beam to length ratio
    Cdl = lp.PStates(name='Cdl') #displacement to length ratio
    
    
    print '\n--------------------------------------------------------------'
    print ' Defining Rules for the design space '
    print '\n--------------------------------------------------------------'
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
    
    """-----------------------------------------------
    #rule: prismatic coefficient
    """
    Cp = Cp == Cb/Cmidshp
    
    
    """-----------------------------------------------
    #rule: waterplane coefficient
    #"""
    Cwp = Cwp == Awp/(lwl*bwl)
    
    
    """-----------------------------------------------
    #rule: Centerplane Coefficient
    #"""
    Ccp = Ccp == Acp/(lwl*draft)
    
    
    """-----------------------------------------------
    #rule: midship coefficient
    #"""
    Cmidshp = Cmidshp == Amsh/(bwl*draft)
    
    
    """-----------------------------------------------
    #rule: LCG Coefficient
    #"""
    LCG = LCG == Clcg*lwl
    
    
    
    
    #
    #***********************************
    # Flats
    lfwl = lp.PStates(name='lfwl')      #flat water line
    lfcp = lp.PStates(name='lfcp')      #flat center plane
    lfsac = lp.PStates(name='lfsac')    #flat of SAC
    
    
    """-----------------------------------------------
    #rule:  flat of WL FOWL
    #"""
    lfwl = lfwl <= Cp*Awp/bwl
    
    
    """-----------------------------------------------
    #rule: flat of CL Keel FOCP
    #"""
    lfcp = lfcp <= Cp*Acp/draft
    
    """-----------------------------------------------
    #rule:
        states = (states * (max_height,mid_len,mid_area))  
        mid_area <= Amsh*lfsac
        lfsac = mid_len
        max_height = Amsh
    #"""
    #lfsac = lfsac <= Cp*vol*ia(.2,.4)/Amsh
    lfsac = lfsac <= Cp*vol*ia(.5,.75)/Amsh
        
    #
    #***********************************
    # Add flat rules all to the Database:
    hullrulesnet.set_rules(lfwl,
                           lfcp,
                           lfsac)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    
    
    
    
    #
    #***********************************
    # Add the above rules to the Database
    # like so:
    hullrulesnet.set_rules(disp,
                           Cb,
                           Cp,
                           Cwp,
                           Ccp,
                           Cmidshp,
                           LCG)
    #***********************************
    #
    hullrulesnet.rgp.compute_fresh_rules_graph()
    """-----------------------------------------------
    #rule: beam/Length ratio (to be added)
    #"""
    Cbl = Cbl == bwl/lwl
    
    """-----------------------------------------------
    #rule: displacement/Length ratio (to be added)
    #"""
    Cdl = Cdl == vol/lwl
    #
    #***********************************
    # compile to the relational rules base again:
    hullrulesnet.set_rules(Cbl,Cdl)
    #
    """
    #rule: Cbl interval is set to some fractional range
    """
    Cbl = Cbl == ia(.2,.4)
    #
    #***********************************
    # compile to the relational rules base again:
    hullrulesnet.set_rules(Cbl)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    
    """-----------------------------------------------
    #rule: RULE: flat of SAC <= lfwl
    """
    lfsac = lfsac <= lfwl#*ia(.5,.7) - bad to modify this rule.  why? 
    # because we already have a rule for the fraction relating these two
    # it's [> lfwl = lfwl <= lfsac*ia(1.0,1.2) <] this guy
    #
    #***********************************
    # Add rule to the Database:
    hullrulesnet.set_rules(lfsac)
    
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    """-----------------------------------------------
    #rule: flat of SAC <= lfcp
    #"""
    lfsac = lfsac <= lfcp
    
    
    #
    #***********************************
    # Add rule to the Database:
    hullrulesnet.set_rules(lfsac)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    
    
    #
    #***********************************
    # more Flat rules
    
    #rule:  flat of WL FOWL
    lfwl = lfwl <= lfsac*ia(2.,3.5)
    
    
    #rule: flat of CL Keel FOCP
    lfcp = lfcp <= lfsac*ia(1.1,1.5)
        
        
    #
    #***********************************
    # Add flat rules all to the Database:
    hullrulesnet.set_rules(lfwl,
                           lfcp)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    
    
    lfcp = lfcp <= lfwl
    hullrulesnet.set_rules(lfcp)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    #
    #***********************************
    # Set the design space intervals themselves
    # not done here.  Using something else to do the same
    
    
    #
    #***********************************
    # rule: instantiate Cp with tight limits
    #Cp = Cp == ia(.7,.88)
    hullrulesnet.set_rules(Cp)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    print 'A few more rules  - '
    
    
    
        
    print '\n--------------------------------------------------------------'
    print 'Now the Design Space is Pretty Well Built Up '
    print 'Lets make some specific choices about it'
    print '--------------------------------------------------------------\n'
    
    lwl = lwl == ia(80.,130.) #how about a design with some ballpark near 100 "meter" hull?
    #lwl is the allowable range of the design space length along the waterline
    # this is the only measure of "global length" we will use in these simplistic examples.
    
    
    draft = draft == ia(12.,21.) #design space for draft - an interval ranging from 12 to 21 length units.
    
    bwl = bwl == ia(22.,34.) #design space for breadth 
    #(NOTE!: full width required expected here for the form parameter code - to come after design space narrowing.)
    
    vol = vol == ia(1000.,25000.) #volume of the ship hull 
    #(not worrying about density of fluid this second - yes our weightless vessel has a waterline.... I am sorry!)
    
    LCG = LCG == ia(30.,70.) #longitudinal center of volume #better to specify this closely with an 
    # interval coefficient.  See Clcg.
    
    Clcg = Clcg == ia(.45,.49) #ratio of LCG to lwl - see rule above.  This interval allows the LCG to be in a 
    # narrow range around midship.
    
    #"""
    hullrulesnet.print_hull_state()
    hullrulesnet.set_rules(lwl,
                          draft,
                          bwl,
                          vol,
                          LCG,
                          Clcg)
    
    
    
    Cb = Cb == ia(0.3,0.65)  # block coeficient allowable range of the design space. (maybe a good osv is not to fat?)
    
    Cmidshp = Cmidshp == ia(0.94,.99) #full midsection desired
    
    
    Cwp = Cwp == ia(.81,.95) #water plane coefficient
    
    Ccp = Ccp == ia(.76,.85) #centerplane coefficient
    
    hullrulesnet.set_rules(Cb,
                           Cmidshp,
                           Cwp,
                           Ccp)
    
    ## See rules above for the definitions of these rules.  Here we mostly just assign some 
    ## reasonable interval ranges for our design space.
    ##
    
        
    print '\n--------------------------------------------------------------'
    print 'Start Generation of Hull with Bulbous Bow '
    print 'representation type: THB-spline'
    print '--------------------------------------------------------------\n'
    SD = ShipDesigner(design_space=hullrulesnet)
    SD.bare_hull_random_design_selection_new() #filter the space
    SD.make_bare_hull() #FPD
    
        
    print '\n--------------------------------------------------------------'
    print 'Generate Bulbous Bow'
    print '--------------------------------------------------------------\n'
    SD.bulbous_bow_random_design_selection() #filter the space
    SD.make_bulbous_bow() #FPD
        
    print '\n--------------------------------------------------------------'
    print 'Combine Hull as Truncated Hierarchical B-spline Surface'
    print '--------------------------------------------------------------\n'
    #SD.make_THBspline_complex_hull() #THB-splines
    print 'Oops - I am depending on caching the projection matrices it seems'
    print '(in this file)'
    print 'TODO: generate the projection matrices and memoize them!'
    return SD

def plot_stuff(SD):
    SD.hull.plot_primary_curves()
    SD.plot_bare_hull_results()
    SD.hull.plot_hull_transverse_curves_and_longitudinal_curves()
    #SD.plot_THB_ship_hull()
    return 

if __name__ == '__main__':
    ShipDesigner = main()
    #plot_stuff()
    #ShipDesigner.plot...