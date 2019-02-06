#!/usr/bin/env python2
#
#To set a better module encoding 
# add the following comment 
# to the first or second line of the Python module 
# using the Unicode literal:
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 15:52:55 2018

@author: luke
"""

""" #Bizarrely, this does not work right?:-
#import relational_lsplines as rlspline 
##from relational_lsplines import simple_hull_rules_language as srl
##from relational_lsplines.simple_hull_rules_language import HullGeometryGenerator
#
##import opt_simple_hull as opt
#
#ShipDesigner = rlspline.ShipDesigner
#
#
#make_hull_rules_net         = rlspline.simple_hull_rules_language.make_hull_rules_net
#set_design_space            = rlspline.simple_hull_rules_language.set_design_space
#combine_space_with_rules    = rlspline.simple_hull_rules_language.combine_space_with_rules

#hullrulesnet = make_hull_rules_net(
#                rlspline.HullGeometryGenerator(rtype='gauss',
#                                               verbose=True))
#
#DS = set_design_space()
#
#hullrulesnet = combine_space_with_rules(DS,
#                                        hullrulesnet)
#"""
## 
##################################################################
## jinja2 latexing...
#http://eosrei.net/articles/2015/11/latex-templates-python-and-jinja2-generate-pdfs
#import jinja2
#import os
#from jinja2 import Template
##   
##################################################################
## 
import FileTools #my stuff for file handling with slow ease.
##   
##################################################################
## ship designer setup

## but this does work right!:-
import opt_simple_hull as opt
#
ShipDesigner = opt.ShipDesigner#Design Helper class
# which abstracts out the various processes of designing a hull
# 'the random way'
#
DesignSpecification = opt.DesignSpecification #tiny class
# which maps from tuples to intervals
# the idea is to instantiate the design space with something familiar (tuple)
#
# personally, I prefer to use the tiny language itself as that's what
# captures the design space intent.
# And that's where the 'real machinary' is.
#
#
#
lp = opt.lp #everything in sqKanren.py
np = lp.np #numpy
#
ia = lp.ia #use ia import from lp instead to match isinstance
#
RulesGraphProcessor = lp.RulesGraphProcessor
#
import random
import sobol #actually using python's random.random at the moment.
# -this was only because I percieved sobol not to be 'random enough'
# -but after some experience with random.random, things seem similar.
import copy 
##
##
##


import simple_hull_rules_language as shrl



hullrulesnet = shrl.make_hull_rules_net(
                shrl.HullGeometryGenerator(rtype='gauss',
                                               verbose=True))

DS = shrl.set_design_space()

hullrulesnet = shrl.combine_space_with_rules(DS,
                                        hullrulesnet)

##  End ship designer setup
##################################################################
## Jinja2 setup for latex:


#latex_jinja_env = jinja2.Environment(
#                block_start_string= '\BLOCK{',
#                block_end_string = '}',
#                variable_start_string = '\VAR{',
#                variable_end_string = '}',
#                comment_start_string = '\#{',
#                comment_end_string = '}',
#                line_statement_prefix = '%%',
#                line_comment_prefix = '%#',
#                trim_blocks = True,
#                autoescape = False,
#                loader = jinja2.FileSystemLoader(os.path.abspath('.'))
#                )
#
#pre = '/home/luke/Documents/computational_naval_architecture/projects/relational_hull_design'
#directory = pre+'/relational_lsplines/relational_lsplines'
#filename =   'TemplateTable.tex'

#template = latex_jinja_env.get_template(filename)
#ans = template.render()


def get_template():
    pre = '/home/luke/Documents/computational_naval_architecture/projects/relational_hull_design'
    directory = pre+'/relational_lsplines/relational_lsplines'
    filename =   'TemplateTable.tex'
    testlines = FileTools.GetLines(directory,filename)
    for line in testlines:
        print line
    return testlines

## 
##################################################################
## 

if __name__ == "__main__":

    
    
    
    
    SD = ShipDesigner(design_space=hullrulesnet)
    HyperParameters = SD.hp
    
    importthis = 'designspace_wavy_March1'
    importthis = 'designspace_March3'
    
    SD.Import(option='designspace')
    
    
    tl = SD.design_space.get_thick_list()
    if tl:
        for el in tl:
            print SD.design_space.rgp.get_name_and_value(el)
            
            
        tl = SD.design_space.get_node_list()
        for el in tl:
            print SD.design_space.rgp.get_name_and_value(el)
    
    
    


    #from string import Template
    lorem = "Lorem ipsum dolor sit amet {GIBBERISH}, consectetur adipiscing elit {DRIVEL}. Expectoque quid ad id, quod quaerebam, respondeas."
    #loremtpl = Template("Lorem ipsum dolor sit amet $GIBBERISH, consectetur adipiscing elit $DRIVEL. Expectoque quid ad id, quod quaerebam, respondeas.")
    d = dict(GIBBERISH='FOOBAR', DRIVEL = 'RAXOOP')
    print lorem.format(**d)
    #print loremtpl.substitute(d)
    
    
    
    ## 
    ##################################################################
    ## map kanren to latex
    tl = SD.design_space.get_node_list()
    k2l = {}
    #l2k = {}
    for el in tl:
        k2l[el] = el.name
        #l2k[el.name] = None
    l2k={'Acp'  : r'\acp',#done in code
         'Amsh' : r'\am', #done in code
         'Awp'  : r'\awp',#done in code
         'Cb'   : r'\cb',
         'Ccp'  : r'\ccp',
         #'Ccp'  : r'\CPK',
         'Cdl'  : r'\dlr', #Disp to Length Ratio, approx page 93, BasicRules.tex 
         'Clb'  : r'\Clb', #length to beam ratio
         'Clcg' : r'\Clcg',
         'Cmidshp'  : r'\cmdshp',
         'Cp'   : r'\cp',
         'Cwp'  : r'\cwp',
         'LCG'  : r'\lcb',
         'bwl'  : r'\beam',
         'disp' : r'\volDisp',
         'draft': r'\draft',
         'lfcp' : r'\lenFCPK',
         'lfsac': r'\lfsac',
         'lfwl' : r'\lenFWL', #\lenFSAC
         'lwl'  : r'\lwl',
         'vol'  : r'\volDisp'} #done in code
    
    il2k = {}
    for key in l2k:
        v = l2k[key]
        il2k[v] = key
    
    
    #for el in k2l:
    #    var, val = SD.design_space.rgp.get_name_and_value(el)
    #    print l2k[var.name]
        
    lines = get_template()
    #
    #  line = lines[25]
    #
    for line in lines:
        tokens = line.split()
        if not tokens:
            continue #skips to top of loop
        ident = tokens[0]
        if '$' in ident:
            keyv = ident[1:-1]
            try:
                #print keyv, il2k[keyv]
                var, val = SD.design_space.rgp.get_name_and_value(il2k[keyv])
                #print var, val
                d = dict(val=val[0],
                         act = None,
                         dif = 0.0)
                print line.format(**d)
            except:
                pass
        else:
            print line