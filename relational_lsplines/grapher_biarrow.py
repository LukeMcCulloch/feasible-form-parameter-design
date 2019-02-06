#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 28 10:57:41 2017

@author: 
    http://www.programcreek.com/python/example/5579/pydot.Dot
    
"""

# -*- encoding: utf-8 -*-
"""
pydot graph example
@author: Francisco Portillo
@url: https://gist.github.com/fportillo/a499f9f625c6169524b8
"""
import pydot

# Create the graph
#graph = pydot.Dot(graph_type='digraph')
graph = pydot.Dot(graph_type='graph')


# Create Parameters
disp    = pydot.Node("disp")
bwl     = pydot.Node("bwl")
lwl     = pydot.Node("lwl")
draft   = pydot.Node("draft")

Awp   = pydot.Node("Awp")
Amsh   = pydot.Node("Amsh")

#Rules (which are also parameters ;)
Cp      = pydot.Node("Cp")
Cb      = pydot.Node("Cb")
Cwp      = pydot.Node("Cwp")
Ccp      = pydot.Node("Ccp")
Cmidship      = pydot.Node("Cmidship")



graph.add_node(Cb)
graph.add_node(Cp)
graph.add_node(Cwp)
graph.add_node(Ccp)
graph.add_node(Cmidship)
graph.add_node(disp)
graph.add_node(bwl)
graph.add_node(lwl)
graph.add_node(draft)
graph.add_node(Awp)


#
# Block Coefficient
graph.add_edge(pydot.Edge(lwl, Cb,dir="both"))
graph.add_edge(pydot.Edge(bwl, Cb,dir="both"))
graph.add_edge(pydot.Edge(draft, Cb,dir="both"))
graph.add_edge(pydot.Edge(disp, Cb,dir="both"))

#
# waterplane coefficient
graph.add_edge(pydot.Edge(lwl, Cwp,dir="both"))
graph.add_edge(pydot.Edge(bwl, Cwp,dir="both"))
graph.add_edge(pydot.Edge(Awp, Cwp,dir="both"))

#
# Centerplane Coefficient
graph.add_edge(pydot.Edge(lwl, Ccp,dir="both"))
graph.add_edge(pydot.Edge(draft, Ccp,dir="both"))
graph.add_edge(pydot.Edge(Awp, Ccp,dir="both"))

#
#
# Cmidship
graph.add_edge(pydot.Edge(Amsh, Cmidship,dir="both"))
graph.add_edge(pydot.Edge(bwl, Cmidship,dir="both"))
graph.add_edge(pydot.Edge(draft, Cmidship,dir="both"))


graph.write_png("graphs/BlockCoefficient.png")

