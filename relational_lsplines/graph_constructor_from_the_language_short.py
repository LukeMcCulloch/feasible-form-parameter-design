#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 21:13:35 2017

@author: luke

TODO:
    
    THIS IS NOT RIGHT YET!
"""
from __future__ import print_function


#from design_tree import HullDesignNode as hullclp
#from simple_hull_rules import HullDesignNode as hullclp

#from simple_hull_rules_language import HullGeometryGenerator as hullclp
from simple_hull_rules_language import *
import relational_lsplines as rlspline 

#from design_tree import DesignTree
#import inspect
import pydot
import networkx as nx


def make_hull_rules_net(hullrulesnet):
    """Primary Bare Hull Rules
    -----------------------------------
    
    -minimal requirements for usage
    -this kind of thing works 'on the fly'
        -what do you mean?
        -how do you demonstrate this?
    
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
    
    
    TODO:
    ----------------------------------
        -rule to constrain the total length of the
        flat of side
        
        -how about regulating where the bow and stern
        fairness curves can be?
        
        -how about the fore and aft extents of the flat spot?
        
        DevNotes: somehow DWL and CPKeel (but especially DWL)
        work 'better' than SAC rules, despite having fewer points
    """
    
    #
    #***********************************
    # Instantiate the design space variables
    lwl = lp.PStates(name='lwl')
    bwl = lp.PStates(name='bwl')
    draft = lp.PStates(name='draft')
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
    
    Clb = lp.PStates(name='Clb') #length to beam ratio
    Cdl = lp.PStates(name='Cdl') #displacement to length ratio
    
    
    Cld = lp.PStates(name='Cld') #length to depth ratio (not really used)
    
    
    
    #
    #***********************************
    # Flats
    lfwl = lp.PStates(name='lfwl')      #flat water line
    lfcp = lp.PStates(name='lfcp')      #flat center plane
    lfsac = lp.PStates(name='lfsac')    #flat of SAC
    
    
    #***********************************
    #***********************************
    #***********************************
    #
    # 3 Section SAC Curve Variables
    #
    # 3 lengths, areas, centers
    #
    #***********************************
    # rules:  3-part variables
    SAC_entrance_area = lp.PStates(name='SAC_entrance_area')
    SAC_mid_area = lp.PStates(name='SAC_mid_area')
    SAC_run_area = lp.PStates(name='SAC_run_area')
    
    SAC_entrance_len = lp.PStates(name='SAC_entrance_len')
    #lfsac = lp.PStates(name='lfsac')
    SAC_run_len = lp.PStates(name='SAC_run_len')
    
    SAC_fwd_Xc = lp.PStates(name='SAC_fwd_Xc')
    SAC_mid_Xc = lp.PStates(name='SAC_mid_Xc')
    SAC_run_Xc = lp.PStates(name='SAC_run_Xc')
    
    
    
    
    
    """-----------------------------------------------
    #rule: block coefficient
    """
    Cb = Cb == disp/(lwl*bwl*draft)
    hullrulesnet.set_rules(Cb)
    
    
    """-----------------------------------------------
    #rule: prismatic coefficient
    """
    Cp = Cp == Cb/Cmidshp
    hullrulesnet.set_rules(Cp)
    
    
    """-----------------------------------------------
    #rule: waterplane coefficient
    #"""
    Cwp = Cwp == Awp/(lwl*bwl)
    hullrulesnet.set_rules(Cwp)
    
    
    """-----------------------------------------------
    #rule: Centerplane Coefficient
    #"""
    Ccp = Ccp == Acp/(lwl*draft)
    hullrulesnet.set_rules(Ccp)
    
    
    """-----------------------------------------------
    #rule: midship coefficient
    #"""
    Cmidshp = Cmidshp == Amsh/(bwl*draft)
    hullrulesnet.set_rules(Cmidshp)
    
    
    """-----------------------------------------------
    #rule: LCG Coefficient
    #"""
    #LCG = LCG == Clcg*lwl
    Clcg = Clcg == LCG/lwl
    hullrulesnet.set_rules(Clcg)
    
    
    
    
    return hullrulesnet
    




def set_design_space():
    """
    USER INPUT HERE:
    
        *Choose your design space input parameters here
        
    
    TODO: add design space rules input hooks for the bbow?
    
    
    use 
    SD.export('designspace')
    SD.export('rhino')
    to save the design
    and the simple curves for ComplexHull maker
    to turn into a 5 part surface for Rhino.
    
    
    1.) narrow the mid-flat zone
        to make more volume go to towards the ends
    """
    return DesignSpecification(
                            #lwl = (80.,130.), #making tubs
                            lwl = (110.,160.),
                            draft = (15.,26.),
                            bwl = (25.,34),#35.), #full width of vessel
                            disp=(1000.,35000.),
                            #
                            # no more tubs, so we have to adapt LCG (no bearing on tubiness though)
                            LCG = (55.,85.), #square DWL aft may make 
                            #
                            # for more disp than you realize
                            Clcg = (.48,.52), #location Coeff for the LCG
                            Cb = (.5,.75), #high prismatic, low block.
                            #
                            Cmidshp = (0.9,.95), #midship coefficient
                            #
                            Cwp = (.89,.94),
                            Ccp = (.78,.87)
                            )
"""Cwp should be > Ccp
since the only difference (besides the 
antisymmetric appearance)
is that Awp is -really square- aft
while Acp is not square fwd.

(At other opposite ends they are similar.)


Cwp should be > Ccp
because Ccp has curveature fore and aft.
Cwp is a box aft and curvy fwd.  ergo Cwp has more area.
"""





def combine_space_with_rules(DS,hullrulesnet):
    """TODO: make ia accept tuples of len 2
    """
    for key in DS.__dict__:
        new_value = DS.__dict__[key]
        var,val = hullrulesnet.rgp.get_name_and_value(key)
        #
        # We could set this as if using '=' on the calss attribute:
        #setattr(hullrulesnet,var.name,ia(new_value[0],new_value[1])) 
        #
        # Or set this using the logic language directly
        #
        # still have to get this class attribute on the fly though, so:
        this_var = getattr(hullrulesnet,key) 
        #
        # here is the rule
        this_var = this_var == ia(new_value[0],new_value[1]) 
        #
        # add it to the database
        hullrulesnet.rgp.add_one_rule(this_var, 
                                      this_var.name) 
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    return hullrulesnet




def pydot_graph_of_rules(rlist,
                         rgp):
    """
    #rlist = SD.design_space.get_node_list()
    rlist = newlist
    rgp=SD.design_space.rgp
    
    
    rlist = ls3
    nlist = ls2
    
    """
    ## map 
    ##       function.name : function
    #nodes = rgp.nodes
    methods = {}
    for elem in rgp.nodes:#rlist:
        #parameter,val = rgp.get_name_and_value(elem)
        parameter = elem.name
        #for arc in parameter.arcs:
        #for gs_node in parameter.flist: #they are connected at random, so you can't do this.
        #    if gs_node.name in methods.keys():
        #        pass
        #    else:
        #        methods[gs_node.name]=gs_node
        methods[parameter] = elem
    #
    ## graph 
    ##       initialize
    graph = pydot.Dot(graph_type='graph',
                      overlap='scale',
                      splines='true')
    #
    ## map 
    ##       graphed_function.name : pydot.Node(shape = 'box', style = 'filled)
    mrules = {}
    for method in methods:
        rule = pydot.Node(method, 
                          shape='box',
                          style='filled')
        mrules[method] = rule
        graph.add_node(rule)
    #
    #
    #
    for parameter in rlist:
        pnode = pydot.Node(parameter.name)
        graph.add_node(pnode)
        
        
        parameter,val = rgp.get_name_and_value(parameter)
        #for method in parameter.arcs:
        #print parameter
        for method in parameter.flist:
            graph.add_edge(pydot.Edge(pnode, mrules[method.name]))
    return graph
    
    
   
def plot_some_rules(graph,name):
    #
    #
    #
    subfolder = name.split('.')[0]
    #
    #----------------------------------------------- Plotting
    #
    graph.write_pdf("graphs/"+subfolder+"/_circlegraph_"+name,
                    prog='circo')
    graph.write_pdf("graphs/"+subfolder+"/_neato_"+name,
                    prog='neato')
    return
    
    
    
    
def graph_andplot_rules(rlist,
                        rgp,
                        name):
    graph = pydot_graph_of_rules(rlist,rgp)
    plot_some_rules(graph,name)
    return 
    
    
def plot_netx_of_pydot_graph(pydot_graph):
    """
    input:      pydot graph object
    
    retult:     Plot a networkx graph made from a pydot graph
    """
    netx_primary = nx.nx_pydot.from_pydot(pydot_graph)
    nx.draw(netx_primary)
    return
    


if __name__ == '__main__':
    #hdp = hullclp()
    
    
    
    #print '\n--------------------------------------------------------------'
    #print 'Start Graph of Hull Relations '
    #print '--------------------------------------------------------------\n'
    spp = rlspline.ADILS.SolverPostProcessor
    
    hullrulesnet = make_hull_rules_net(HullGeometryGenerator(rtype='gauss',
                                                             verbose=True) ) 
    
#    DS = set_design_space()
#    
#    hullrulesnet = combine_space_with_rules(DS,
#                                            hullrulesnet)
#        
#    
#    SD = ShipDesigner(design_space=hullrulesnet)
#    #    
#    
#    RuleList = SD.design_space.get_node_list()
#    
#    
#    
#    graph_andplot_rules(rlist = RuleList,
#                        rgp = hullrulesnet.rgp,
#                        name='Big.pdf')
    
    
    RuleList = hullrulesnet.get_allnode_list()

    graph_andplot_rules(rlist = RuleList,
                        rgp = hullrulesnet.rgp,
                        name='Big.pdf')

    #
    #
    #
    #
    #
    #dotgraph = pydot_graph_of_rules(hdp.alists)
    #plot_netx_of_pydot_graph(pydot_graph=dotgraph)