import  rhinoscriptsyntax as rs
import pickle as pickle#C pickle with no subclassing

from Rhino.Geometry import Point3d
from Rhino.Geometry import NurbsCurve,NurbsSurface
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


rs.AddNurbsCurve
    Returns:
        String
            -The identifier of the new object if successful.
        Null
            -If not successful, or on error.

rs.ProjectCurveToSurface

"""

def doTcurves(key):
    """
    """
    the_filename = 'complex_Tcurvesurface#data_curves.txt'
    a,b = the_filename.split('#')
    the_filename = a+key+b
    with open(the_filename, 'r') as f:
        tlist = []
        ci = []
        for line in f:
            if not 'new' in line:
                ci.append([float(el) for el in line.split(',')])
            else:
                tlist.append(ci)
                ci = []

    the_filename = 'complex_Tcurvesurface#data_knots.txt'
    a,b = the_filename.split('#')
    the_filename = a+key+b
    with open(the_filename, 'r') as f:
        tknots = []
        ci = []
        for line in f:
            if not 'new' in line:
                ci.append([float(el) for el in line.split()])
            else:
                tknots.append(ci[0])
                ci = []

    for verts, knots in zip(tlist,tknots):
        rs.AddNurbsCurve(verts, knots,3)
    return


def doLcurves(key):
    """
    """
    the_filename = 'complex_Lcurvesurface#data_curves.txt'
    a,b = the_filename.split('#')
    the_filename = a+key+b
    with open(the_filename, 'r') as f:
        llist = []
        ci = []
        for line in f:
            if not 'new' in line:
                ci.append([float(el) for el in line.split(',')])
            else:
                llist.append(ci)
                ci = []


    the_filename = 'complex_Lcurvesurface#data_knots.txt'
    a,b = the_filename.split('#')
    the_filename = a+key+b
    with open(the_filename, 'r') as f:
        lknots = []
        ci = []
        for line in f:
            if not 'new' in line:
                ci.append([float(el) for el in line.split()])
            else:
                lknots.append(ci[0])
                ci = []

    for verts, knots in zip(llist,lknots):
        rs.AddNurbsCurve(verts, knots,3)
    return



def doHcurves(key):
    """Hull Curves
    """
    the_filename = 'complex_Hullcurvesurface#data_curves.txt'
    a,b = the_filename.split('#')
    the_filename = a+key+b
    with open(the_filename, 'r') as f:
        tlist = []
        ci = []
        for line in f:
            if not 'new' in line:
                ci.append([float(el) for el in line.split(',')])
            else:
                tlist.append(ci)
                ci = []

    the_filename = 'complex_Hullcurvesurface#data_knots.txt'
    a,b = the_filename.split('#')
    the_filename = a+key+b
    with open(the_filename, 'r') as f:
        tknots = []
        ci = []
        for line in f:
            if not 'new' in line:
                ci.append([float(el) for el in line.split()])
            else:
                tknots.append(ci[0])
                ci = []

    for verts, knots in zip(tlist,tknots):
        rs.AddNurbsCurve(verts, knots,3)
    return


if __name__ == """__main__""":


    keys = ['0','1','2','3','4']
    for key in keys:
        doTcurves(key)
        doLcurves(key)
        doHcurves(key)
