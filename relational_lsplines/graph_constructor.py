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
from simple_hull_rules import HullDesignNode as hullclp
#from design_tree import DesignTree
#import inspect
import pydot
import networkx as nx



def pydot_graph_of_rules(rlist):
    
    ## map 
    ##       function.name : function
    methods = {}
    for parameter in rlist:
        for arc in parameter.arcs:
            if arc.__name__ in methods.keys():
                pass
            else:
                methods[arc.__name__]=arc
    #
    #
    #
    ## graph 
    ##       initialize
    graph = pydot.Dot(graph_type='graph',
                      overlap='scale',
                      splines='true')
    #
    #
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
    ## GRAPH 
    ##       node(variable.name) : pydot.Node(default)
    for parameter in rlist:
        pnode = pydot.Node(parameter.name)
        graph.add_node(pnode)
    
        for method in parameter.arcs:
            graph.add_edge(pydot.Edge(pnode, mrules[method.__name__]))
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
    
    
    
    
def graph_andplot_rules(rlist,name):
    graph = pydot_graph_of_rules(rlist)
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
    
    hdp = hullclp()
    #
    #
    #
    graph_andplot_rules(hdp.Primary, name='primary.pdf')
    graph_andplot_rules(hdp.Coefficients, name='coefficients.pdf')
    graph_andplot_rules(hdp.alists, name='all.pdf')
    #
    #
    #
    dotgraph = pydot_graph_of_rules(hdp.alists)
    plot_netx_of_pydot_graph(pydot_graph=dotgraph)