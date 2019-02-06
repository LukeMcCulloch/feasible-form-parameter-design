

try:
    import  rhinoscriptsyntax as rs
    #from THBLSplineHullcomplex import *
    #import pickle as pickle#C pickle with no subclassing
    #import cPickle as pickle#C pickle with no subclassing
    tryRhino = True
    from Rhino.Geometry import Point3d
    from Rhino.Geometry import NurbsCurve,NurbsSurface
except:
    tryRhino = False
    print 'RhinoScript only availible in Rhino... not compatible (directly) with Other Pythons :('

"""
#
# usage:
#   First
#       start Rhino and type
#           "StartAtomEditorListener"
#           in the command line there
#
#   Then
#       packages/RhinoPython/save and run in Rhino
#
# That's it!  (Assuming the curves are where they should be and paths are good)
# TODO:  Update this doc to explain what curve files are needed and where
#           any required files must "live"
#
#
#
#
# to send the file to Rhino for execution press the ctrl + alt + r keys.
#
#
#
# see also
#   https://atom.io/packages/rhino-python
"""

the_filename = 'transverse_curves.txt'
with open(the_filename, 'r') as f:
    tlist = []
    ci = []
    for line in f:
        if not 'new' in line:
            ci.append([float(el) for el in line.split(',')])
        else:
            tlist.append(ci)
            ci = []



the_filename = 'transverse_knots.txt'
with open(the_filename, 'r') as f:
    tknots = []
    ci = []
    for line in f:
        if not 'new' in line:
            ci.append([float(el) for el in line.split()])
        else:
            tknots.append(ci[0])
            ci = []

the_filename = 'longitudinal_curves.txt'
with open(the_filename, 'r') as f:
    llist = []
    ci = []
    for line in f:
        if not 'new' in line:
            ci.append([float(el) for el in line.split(',')])
        else:
            llist.append(ci)
            ci = []

the_filename = 'longitudinal_knots.txt'
with open(the_filename, 'r') as f:
    lknots = []
    ci = []
    for line in f:
        if not 'new' in line:
            ci.append([float(el) for el in line.split()])
        else:
            lknots.append(ci[0])
            ci = []

if tryRhino:
    for verts, knots in zip(tlist,tknots):
        rs.AddNurbsCurve(verts, knots,3)

    for verts, knots in zip(llist,lknots):
        rs.AddNurbsCurve(verts, knots,3)



#TODO:  rs.AddNurbsSurface  !!
