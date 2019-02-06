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
graph = pydot.Dot(graph_type='graph',
                  splines='true')#,
                #dimen="10,10")



# Create Parameters
disp    = pydot.Node("disp")
bwl     = pydot.Node("bwl")
lwl     = pydot.Node("lwl")
draft   = pydot.Node("draft")

Awp   = pydot.Node("Awp")
Acp   = pydot.Node("Acp")
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
graph.add_node(Acp)
graph.add_node(Amsh)


#
# Block Coefficient
Block_Coefficient      = pydot.Node(
                    "Block Coefficient \n Disp/(Lwl*Bwl*Draft) = Cb ",
                    shape='box')
graph.add_node(Block_Coefficient)
graph.add_edge(pydot.Edge(lwl, Block_Coefficient,
                          decorate=1))
graph.add_edge(pydot.Edge(bwl, Block_Coefficient))
graph.add_edge(pydot.Edge(draft, Block_Coefficient))
graph.add_edge(pydot.Edge(disp, Block_Coefficient))
graph.add_edge(pydot.Edge(Cb, Block_Coefficient))

#
# waterplane coefficient
Waterplane_Coefficient      = pydot.Node(
                    "Waterplane Coefficient \n Awp/(lwl*bwl) = Cwp ",
                    shape='box')
graph.add_node(Waterplane_Coefficient)
graph.add_edge(pydot.Edge(lwl, Waterplane_Coefficient))
graph.add_edge(pydot.Edge(bwl, Waterplane_Coefficient))
graph.add_edge(pydot.Edge(Awp, Waterplane_Coefficient))
graph.add_edge(pydot.Edge(Cwp, Waterplane_Coefficient))

#
# Centerplane Coefficient
Centerplane_Coefficient      = pydot.Node(
                    "Centerplane Coefficient \n Acp/(Lwl*draft) = Ccp ",
                    shape='box')
graph.add_node(Centerplane_Coefficient)
graph.add_edge(pydot.Edge(lwl, Centerplane_Coefficient,
               splines=True))
graph.add_edge(pydot.Edge(draft, Centerplane_Coefficient,
               splines=True))
graph.add_edge(pydot.Edge(Acp, Centerplane_Coefficient,
               splines=True))
graph.add_edge(pydot.Edge(Ccp, Centerplane_Coefficient,
               splines=True))

#
#
# Cmidship
midship_Coefficient      = pydot.Node(
                    "midship Coefficient \n Amsh/(Bwl*Dsmax) = Cmidshp ",
                    shape='box')
graph.add_node(midship_Coefficient)
graph.add_edge(pydot.Edge(Amsh, midship_Coefficient))
graph.add_edge(pydot.Edge(bwl, midship_Coefficient))
graph.add_edge(pydot.Edge(draft, midship_Coefficient))
graph.add_edge(pydot.Edge(Cmidship, midship_Coefficient))

#
#
# Prismatic Coefficient
Prismatic_Coefficient      = pydot.Node(
                    "Prismatic Coefficient \n Cb / Cmidshp = Cp ",
                    shape='box')
graph.add_node(Prismatic_Coefficient)
graph.add_edge(pydot.Edge(Cb, Prismatic_Coefficient,
               splines='spline'))
graph.add_edge(pydot.Edge(Cmidship, Prismatic_Coefficient))
graph.add_edge(pydot.Edge(Cp, Prismatic_Coefficient))
#graph.add_edge(pydot.Edge(Cp, Cp))

#graph.write_png("graphs/BlockCoefficient.png",
#                prog='neato')

graph.write_png("graphs/BlockCoefficient.png",
                prog='circo')


graph.write_png("graphs/dot.png",
                prog='twopi')