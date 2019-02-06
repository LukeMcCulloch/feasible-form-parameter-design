#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 21:42:32 2017

@author: luke


TODO: -How about generic hooks to store bare hull and bbow?
"""
#import numpy as np

import relational_lsplines as rlspline # everything (short of AD direct surface solvers) for the PHD
import relational_lsplines.opt_simple_hull as opt # design space search and goemtry creation tools
#****s*************************************************************************
# Done with primary import
ShipDesigner = opt.ShipDesigner
DesignSpecification = opt.DesignSpecification  # trivial Design space class
shrl = rlspline.simple_hull_rules_language     # api to the lp rules graph processor and setup function
#
#******************************************************************************
#
#hooks to use the 'real' language for bare hull specification:
HullGeometryGenerator = shrl.HullGeometryGenerator #rules graph processor api
#
#******************************************************************************
#
ia = opt.ia #interval analysis
lp = opt.lp #logic classes
#
spline = rlspline.spline #Bspline curves
IntervalLagrangeSpline = rlspline.IntervalLagrangeSpline #I.A. Bspline curve solver
Lagrangian = rlspline.Lagrangian #curve constraints




def make_hull_rules_net(hullrulesnet):
    """Primary Bare Hull Rules
    -----------------------------------
    
    This shows how to set your own relationships
    for the design space.
    Unfortunately you still need to 
    set the right names in order to have the intended effects!
    
    *There is flexibility in the relationships specified here.
       -Make your own rules!
       -Set Coefficients or other paremters to particular interval ranges
    
    *These relationships are automatically tunred into
     a relational constraint graph 
       -this enforces feasibility in the parameters
       -which later become form parameters in the curve generation
       
    *Drawback:
        -Still no automatic generation of the optimization
         problems which then turn those rules into B-splines
    
    hullrulesnet = HullGeometryGenerator()
    
    >>>  Awp  =  [ia(2125.0, 4290.0)]
    >>>  bwl  =  [ia(25.0, 33.3333333333)]
    >>>  lfwl  =  [ia(0.0, 171.6)]
    
    
    hullrulesnet.bwl = hullrulesnet.bwl == ia(25.,30.)
    hullrulesnet.rgp.add_one_rule(hullrulesnet.bwl)
    hullrulesnet.rgp.AC_revise()
    
    
    >>>  Awp  =  [ia(2125.0, 3861.0)]
    >>>  bwl  =  [ia(25.0, 30.0)]
    >>>  lfwl  =  [ia(0.0, 154.44)]
    
    hullrulesnet.rgp.print_state()
    
    
    ----------------------------------
    
    hullrulesnet.lfwl = hullrulesnet.lfwl == ia(30.,30.)
    hullrulesnet.rgp.add_one_rule(hullrulesnet.bwl)
    hullrulesnet.rgp.AC_revise()
    
    
    hullrulesnet.bwl = hullrulesnet.bwl == ia(25.,30.)
    hullrulesnet.rgp.add_one_rule(hullrulesnet.bwl)
    hullrulesnet.rgp.AC_revise()
    
    
    >>>  Awp  =  [ia(2125.0, 3861.0)]
    >>>  bwl  =  [ia(25.0, 30.0)]
    >>>  lfwl  =  [ia(30.0, 30.0)]
    ----------------------------------
    
    """
    
    #
    #***********************************
    # Instantiate the design space variables
    lwl = lp.PStates(name='lwl')
    bwl = lp.PStates(name='bwl')
    draft = lp.PStates(name='draft')
    vol = lp.PStates(name='vol')
    
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
    
    
    #
    #***********************************
    # Basic Dimension Relations
    
    #RULE: block coefficient
    Cb = Cb == vol/(lwl*bwl*draft)
    
    #RULE: prismatic coefficient
    Cp = Cp == Cb/Cmidshp
    
    #RULE: waterplane coefficient
    Cwp = Cwp == Awp/(lwl*bwl)
    
    #RULE: Centerplane Coefficient
    Ccp = Ccp == Acp/(lwl*draft)
    
    #RULE: midship coefficient
    Cmidshp = Cmidshp == Amsh/(bwl*draft)
    
    #RULE: LCG Coefficient
    LCG = LCG == Clcg*lwl
    
    #-------------------------------
    #TODO: add these new rules
    #RULE: beam/Length ratio
    Cbl = Cbl == bwl/lwl
    #RULE: displacement/Length ratio
    Cdl = Cdl == vol/lwl
    #-------------------------------
    
    #
    #***********************************
    # Add the above rules to the Database
    # like so:
    hullrulesnet.set_rules(Cb,
                           Cp,
                           Cwp,
                           Ccp,
                           Cmidshp,
                           LCG)
    #
    #***********************************
    # compile some new rules:
    hullrulesnet.set_rules(Cbl,Cdl)
    
    hullrulesnet.rgp.compute_fresh_rules_graph()
    #
    #***********************************
    # Flats
    lfwl = lp.PStates(name='lfwl')      #flat water line
    lfcp = lp.PStates(name='lfcp')      #flat center plane
    lfsac = lp.PStates(name='lfsac')    #flat of SAC
    
    #
    #***********************************
    # Flat rules
    
    #rule:  flat of WL FOWL
    lfwl = lfwl <= Cp*Awp/bwl
    
    
    #rule: flat of CL Keel FOCP
    lfcp = lfcp <= Cp*Acp/draft
    
    #RULE: flat of SAC <= Cprismatic*Displacement/Amsh
    lfsac = lfsac <= Cp*vol/Amsh
        
        
    #
    #***********************************
    # Add flat rules all to the Database:
    hullrulesnet.set_rules(lfwl,
                           lfcp,
                           lfsac)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    #
    #***********************************
    # More Flat rules
    # RULE: flat of SAC <= lfwl
    lfsac = lfsac <= lfwl#*ia(.2,.5) - bad to modify this rule.  why? 
    # because we already have a rule for the fraction relating these two
    # it's [> lfwl = lfwl <= lfsac*ia(1.0,1.2) <] this guy
    #
    #***********************************
    # Add rule to the Database:
    hullrulesnet.set_rules(lfsac)
    
    #RULE: flat of SAC <= lfcp
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
    lfwl = lfwl <= lfsac*ia(1.0,1.2)
    
    
    #rule: flat of CL Keel FOCP
    lfcp = lfcp <= lfsac*ia(1.0,1.2)
        
        
    #
    #***********************************
    # Add flat rules all to the Database:
    hullrulesnet.set_rules(lfwl,
                           lfcp)
    
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    
    # NOTE!
    # design space intervals can be set directly:
    
    #
    #***********************************
    # rule: instantiate Cp with tight limits
    Cp = Cp == ia(.7,.83)
    hullrulesnet.set_rules(Cp)
    
    
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    """--------------------------------
    #rule: alias 'vol' with 'disp'
    """
    disp = lp.PStates(name='disp')
    disp = disp == vol
    hullrulesnet.set_rules(disp)
    
    
    return hullrulesnet



def set_design_space():
    """
    USER INPUT HERE:
    
        *Choose your design space input parameters here
        
    
    TODO: add design space hook for the bbow?
    -at present internal design space settings
    define the bulb rules in the langauge (same as the hull above)
    -rules are also hidden
    -All of these bulb setting are found within GeometryGenerator
    """
    return DesignSpecification(
                            lwl = (80.,130.),
                            draft = (12.,21.),
                            bwl = (20.,30.), #full width of vessel
                            vol=(4000.,25000.),
                            LCG = (30.,70.),
                            Clcg = (.48,.515), #location Coeff for the LCG
                            Cb = (0.5,0.82),
                            Cwp = (.75,.93),#water plane coefficient
                            Cmidshp = (0.9,.99), #midship coefficient
                            Ccp = (.75,.93) #centerplane coefficient
                            )



def postprocess_example(SD):
    """Using the ADILS SolverPostProcessor class
    to check Lspline data for unmet parameters
    """
    spp = rlspline.ADILS.SolverPostProcessor
    try:
        print '\n--------------------------------------------------------------'
        print 'SAC: checking Lspline Optimization:'
        SACpp = SD.hull.Lsplines.SAC.postprocess
    except:
        print '\n--------------------------------------------------------------'
        print 'SAC failed to solve'
        Lspline = SD.hull.Lsplines.SAC
        SACpp = spp(Lspline)
    try:
        print '\n--------------------------------------------------------------'
        print 'Stern fairning curve: checking Lspline Optimization:'
        sfcpp = SD.hull.Lsplines.sternfairing.postprocess
    except:
        print '\n--------------------------------------------------------------'
        print 'stern fairness curve failed to solve'
        sfcpp = spp(SD.hull.Lsplines.sternfairing)
    #
    # post process of reparameterized CPkeel longitudinal section
    # for when/if it is solved for using Lsplines:
    #
    # removed long_keel: doing knot insert delete to match old CProfile split portions
    #    try:
    #        longkeel = spp(SD.hull.Lsplines.long_keel)
    #        print '\n Checking Longitudinal Keel Split Curve:'
    #        print 'split from the CProfile curve and reparameterized'
    #        print longkeel.badparameters
    #    except:
    #        print 'longitudinal keel split from CProfile failed to solve'
    #        sfcpp = spp(SD.hull.Lsplines.sternfairing)
    return SACpp,sfcpp

def main():

    """
    setup the design space
    #
    # use the rules graph processor to add simple Python-like
    # rules making on top of the sqKanren logical states and goals classes
    # these rules expressions will be compiled 
    # into a list of logical functions,
    # with mapped interdependencies
    # suitable for efficient A.C. (arc consistency) revision
    #
    """
    #hullrulesnet = make_hull_rules_net(HullGeometryGenerator('sobol'))
    hullrulesnet = make_hull_rules_net(HullGeometryGenerator(rtype='gauss',
                                                             verbose=True) )
    DS = set_design_space()
    #
    # now add the rules to the hull rules net:
    #
    hullrulesnet = shrl.combine_space_with_rules(DS, hullrulesnet)
    #
    # API to the FPD tooling used to generate the specific hull
    # curves using the design space, once we narrow it down
    SD = ShipDesigner(design_space = hullrulesnet)
    #
    #SD = ShipDesigner(DS) #old stuff
    #SD.design_space.hdp.print_coeff() #old stuff
    #SD.design_space.hdp.print_primary() #old stuff
    return SD




if __name__ == '__main__':
    

    
    
    SD = main()
    hullrulesnet = SD.design_space
    #hullrulesnet.set_random_generator('sobol') #think harder meta-setattr
    #
    SD.demorun() 
    #------------------------------------
    #demorun does the following:
        #1.)
        #SD.bare_hull_random_design_selection_new()
        #2.)
        #SD.make_bare_hull()
        #3.)
        #SD.bulbous_bow_random_design_selection()
        #4.)
        #SD.make_bulbous_bow()
    #------------------------------------
    SD.hull.plot_primary_curves()
    SD.plot_bare_hull_results()
    
        
        
    SD.hull.SAC.plot3DmultiList(
            SD.hull.lcurvenet[:1],
            [SD.hull.CProfile])
    
    # Add the Bulbous bow to the bare hull as a THB surface:
    #------------------------------------
    #"""
    SD.make_THBspline_complex_hull()
    SD.plot_THB_ship_hull()
    #"""