#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May  5 21:25:52 2018

@author: luke

rhino interface classes
-no numpy allowed in here!
-no pickle allo‰wed in here!
-open text files, read in the points.... seriously.
"""

#"""
try:
    import  rhinoscriptsyntax as rs
    from Rhino.Geometry import Point3d
    from Rhino.Geometry import NurbsCurve,NurbsSurface
    tryRhino = True
    print 'Usage in Rhino:'
    print 'type Command:  StartAtomEditorListener'
    print 'Open this file in Atom (assuming Atom has installed the Rhinoscript plugin)'
    print 'run with: ctrl + alt + r'
except:
    tryRhino = False
    print 'RhinoScript only availible in Rhino... not compatible (directly) with Other Pythons :('
#"""



def CurvetoRhino(object):
    def __init__(self):
        pass



class SurfacetoRhino(object):
    def __init__(self):
        self.point_count = None
        self.control_vertices = []
        self.uknots = []
        self.vknots = []
        self.up = None
        self.vp = None

        self.dataswitch = {0:self.set_point_count,
                           1:self.set_control_vertices,
                           2:self.set_uknots,
                           3:self.set_vknots,
                           4:self.set_udegree,
                           5:self.set_vdegree}

    def set_point_count(self, line):
        self.point_count = [int(line[0]), int(line[1])]
        return

    def set_control_vertices(self, line):
        if tryRhino:
            pt = [float(el) for el in line]
            #pt = rs.AddPoint(pt[0],pt[1],pt[2])
            pt = rs.coerce3dpoint(pt)
            self.control_vertices.append(pt)
        else:
            self.control_vertices.append([float(el) for el in line])
        self.weights = [1. for el in self.control_vertices]
        return

    def set_uknots(self, line):
        self.uknots = [float(el) for el in line]
        return

    def set_vknots(self, line):
        self.vknots = [float(el) for el in line]
        return

    def set_udegree(self, line):
        self.up = int(line[0])
        return

    def set_vdegree(self, line):
        self.vp = int(line[0])
        return

    def GetRhinoSurface(self,
                        the_filename='surface_data.txt'
        ):

        dispatch = 0
        with open(the_filename, 'r') as f:
            for line in f:
                line = line.lower().split()
                if line:
                    if 'new' and 'surface' in line:
                        pass
                    if 'end' and 'surface' in line:
                        pass
                    elif 'point' and 'count' in line:
                        dispatch = 0
                        pass
                    elif 'vertices' in line:
                        dispatch = 1
                        pass
                    elif 'uknot' in line:
                        dispatch = 2
                        pass
                    elif 'vknot' in line:
                        dispatch = 3
                        pass
                    elif 'degree' and 'u' in line:
                        dispatch = 4
                        pass
                    elif 'degree' and 'v' in line:
                        dispatch = 5
                        pass
                    else:
                        self.dataswitch[dispatch](line)


        return

    def print_properties(self):
        print 'point count:'
        print self.point_count
        print 'len of control points list'
        print len(self.control_vertices)
        print 'type of control points'
        print type(self.control_vertices[0][0])
        print 'uknots : \n',self.uknots
        print 'vknots : \n',self.vknots
        print 'udegree : \n',self.up
        print 'vdegree : \n',self.vp
        return


class ComplexHulltoRhino(object):
    def __init__(self):

        the_filename = 'surface_data.txt'
        keys = [0,1,2,3,4]
        self.maker = {}
        for key in keys:
            self.maker[key] = SurfacetoRhino()
            self.maker[key].GetRhinoSurface(
                the_filename='surface'+str(key)+'data.txt'
                )
            s1 = self.maker[key]
            rs.AddNurbsSurface(point_count = s1.point_count,
                                points = s1.control_vertices,
                                knots_u = s1.uknots,
                                knots_v = s1.vknots,
                                degree = [s1.up,s1.vp],
                                weights=None)
        pass


class SingleHulltoRhino(object):
    def __init__(self):

        the_filename = 'longsinglesurf.txt'
        keys = [0]
        self.maker = {}
        for key in keys:
            self.maker[key] = SurfacetoRhino()
            self.maker[key].GetRhinoSurface(
                the_filename=the_filename
                )
            s1 = self.maker[key]
            rs.AddNurbsSurface(point_count = s1.point_count,
                                points = s1.control_vertices,
                                knots_u = s1.uknots,
                                knots_v = s1.vknots,
                                degree = [s1.up,s1.vp],
                                weights=None)
        pass


if __name__ == """__main__""":
    s1 = SurfacetoRhino()
    the_filename='surface_data.txt'
    self = s1
    s1.GetRhinoSurface()

    if tryRhino:
        #Rhino.Command("_-Intersect", 0)
        self = ComplexHulltoRhino()
        #self = SingleHulltoRhino()
        #        rs.AddNurbsSurface(point_count = s1.point_count,
        #                            points = s1.control_vertices,
        #                            knots_u = s1.uknots,
        #                            knots_v = s1.vknots,
        #                            degree = [s1.up,s1.vp],
        #                            weights=None)
