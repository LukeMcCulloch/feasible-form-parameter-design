# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 16:39:03 2015

@author: lukemcculloch

Contains:
    -FormParameterDict Class
    -Formparmeter Class
    
    Helper functions to 
    properly setup the constraints for 
    Form Parameter Curve design in 2 dimensions
    with some (essential, but limited) capability
    in 3D
    
    
    
TODO:
    bring back surface constraints, e.g. gaussian curvature
    Automatic Differentation of Surfaces and the surface solver
"""
import numpy               as np
#import AF
from  interval_arithmetic import ia
import copy
from area import area, area_to_y_axis, area_to_any_x_axis
from moments import Xc, Yc,Xc_Yaxis,Yc_Yaxis
from vertex import x_vertex, y_vertex, z_vertex, \
                    relative_x_vertex, relative_y_vertex
from points import x_point, y_point, z_point
from tangent import tangent,                 \
                    vertical_tangent,        \
                    vertical_XoverZ_tangent, \
                    tangent_XoverZ,          \
                    harries_tangent
from curvature import curvature, curvature_XoverZ
from fairness import E1, E2, E3, ArcLength, ArcLengthApprox
from attractor import xAttractor, yAttractor

pi = np.pi


def generalized_aattractor(function, method, vertices, weight, value, attraction):
    #result = ((function(method, vertices, value) - attraction)**2)*weight
    result = function(method, vertices, value)*attraction
    return result

class FormParameterDict(object):
    
    def __init__(self, curve,parent=None,child=None):
        self.valid_kinds = ['equality', 'max', 'min', 'LS','fix']
        self.nec        = 0
        self.nic        = 0
        self.nc         = self.nec + self.nic
        self.curve      = curve
        self.FormParam  = {}
        self.FixParam   = {}
        self.THBchange  = {}
        self.weights    = []
        #*************************************
        # hierarchical structure for multigrid
        self.parent     = parent
        self.child      = child
        return
    
    
    def reset(self):#, curve):
        """essentially this is the __init__ method
        """
        self.__init__(self.curve)
        #        self.valid_kinds = ['equality', 'max', 'min', 'LS']
        #        self.nec        = 0
        #        self.nic        = 0
        #        self.nc         = self.nec + self.nic
        #        self.curve      = curve
        #        self.FormParam  = {}
        #        self.weights    = []
        return
    
    def new_prolongation_level(self, new_curve):
        """Reform the same FPD constraints at a new
        prolonged (finer) dyadic level
        (wavelet dyadic nonuniform refinement)
        
        notes
        ----------
            formerly called prolongate
            1 level of 
            wavelet dyadic nonuniform refinement
        """
        NewFPD = FormParameterDict(new_curve,parent=self)
        for i in self.FormParam:
            c = self.FormParam[i]
            c.callback(NewFPD) #implicit setting of the constraints
        
        parent = new_curve.parent
        #hasfixity = False
        for i in self.FixParam:
            #hasfixity = True
            fp = self.FixParam[i]
            index = fp.index_of_vertex
            #            ga = parent.greville_abscissa(index)
            #            child_real_basis = new_curve.active_basis_at_u(ga)
            #            # = parent
            #            child_real_basis[-1]+=1
            #            child_real_basis = range(*child_real_basis)
            #            #print 'child_real_basis = ',child_real_basis
            #            cindex = child_real_basis[new_curve.k/2]
            #            #print 'parent index = ',index 
            #            #print 'child index  = ',cindex
            
            cindex = new_curve.given_ParentIndex_get_ChildIndex(index)
            
            if fp.track_prolongation:
                fp.callback(NewFPD, newindex = cindex) 
            else:
                fp.callback(NewFPD, newindex = index) 
        
        #if hasfixity:
            #print 'updating the fixity vertex from {} to {}'.format(index, cindex)
            #NewFPD.FixParam[key].index_of_vertex = cindex
            
        self.child = NewFPD
        return NewFPD
    
        
    def add_E1(self, kind = 'LS', location=None, value = 0., weight = 1.0 ):
        #if FPDprolong is None:
        #    constraint_rank = len(self.FormParam)
        #    method = Formparmeter(self.curve, kind, location, value, weight)
        #else:
        #    constraint_rank = len(FPDprolong.FormParam)
        #    method = Formparmeter(FPDprolong.curve, kind, location, value, weight)
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_E1()
        #
        def cbk(FPDprolong):
            return FPDprolong.add_E1(kind, location, value, weight)
        method.callback = cbk
        #
        self.FormParam[constraint_rank] = method
        return
        
    def add_E2(self, kind = 'LS', location=None, value = 0., weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_E2()
        #
        def cbk(FPDprolong):
            return FPDprolong.add_E2(kind, location, value, weight)
        method.callback = cbk
        #
        self.FormParam[constraint_rank] = method
        return
    
    def add_E3(self, kind = 'LS', location=None, value = 0., weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_E3()
        #
        def cbk(FPDprolong):
            return FPDprolong.add_E3(kind, location, value, weight)
        method.callback = cbk
        #
        self.FormParam[constraint_rank] = method
        return

    def add_ArcLength(self, kind = 'LS', location=None, value = 0., weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_ArcLength()
        #
        def cbk(FPDprolong):
            return FPDprolong.add_ArcLength(kind, location, value, weight)
        method.callback = cbk
        #
        self.FormParam[constraint_rank] = method
        return
    
    def add_ArcLengthApprox(self, kind = 'LS', location=None, value = 0., weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_ArcLengthApprox()
        #
        def cbk(FPDprolong):
            return FPDprolong.add_ArcLengthApprox(kind, location, value, weight)
        method.callback = cbk
        #
        self.FormParam[constraint_rank] = method
        return
    
    def add_xAttractor(self, kind = 'LS', location=None, value = None, weight = 1.0, index = None ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_xAttractor(index)
        method.callback = self.add_xAttractor
        self.FormParam[constraint_rank] = method
        return
    
    def add_yAttractor(self, kind = 'LS', location=None, value = None, weight = 1.0, index = None ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_yAttractor(index)
        method.callback = self.add_yAttractor
        self.FormParam[constraint_rank] = method
        return
        
        
    """
    -c
    - o
    -  n
    -   s
    -    t
    -     r
    -      a
    -       i
    -        n
    -         t
    -          s
    """
    #"""
    def add_relative_xVertexConstraint(self, kind = 'equality', location=None, 
                                       value = None, weight = 1.0, 
                                       index = None, index2 = None, 
                                       seperation = 0. ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_relative_x_vertexConstraint(index, index2, seperation)
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_relative_xVertexConstraint(kind, location, 
                                                             value, weight,
                                                             index, index2,
                                                             seperation)
        method.callback = cbk
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
        
    def add_relative_yVertexConstraint(self, kind = 'equality', location=None, 
                                       value = None, weight = 1.0, 
                                       index = None, index2 = None, 
                                       seperation = 0. ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_relative_y_vertexConstraint(index, index2, seperation)
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_relative_yVertexConstraint(kind, location, 
                                                             value, weight,
                                                             index, index2,
                                                             seperation)
        method.callback = cbk
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
    #"""
    #
    #**************************************************************************
    #
    def add_xFixity(self, kind = 'fix', location=None, 
                    value = None, weight = 1.0, index = None,
                                        track_prolongation=True):
        """track_prolongation=False if you do not want the index to 
        update during prolongation of a THB lagrangian to a dyadicly 
        refined level
        """
        constraint_rank = len(self.FixParam)
        method = Formparmeter(self.curve, kind='fix', 
                              value=value,mask=0)
        method.is_x_vertexFixity(index)
        #********************************************
        def cbk(FPDprolong, newindex):
            return FPDprolong.add_xFixity(kind, location, 
                                    value, weight, newindex,
                                    track_prolongation)
        method.callback = cbk
        method.track_prolongation = track_prolongation
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FixParam[constraint_rank] = method
        return
    
    def add_yFixity(self, kind = 'fix', location=None, 
                                value=None, weight = 1.0, index=None,
                                                track_prolongation=True):
        """track_prolongation=False if you do not want the index to 
        update during prolongation of a THB lagrangian to a dyadicly 
        refined level
        """
        constraint_rank = len(self.FixParam)
        method = Formparmeter(self.curve, kind='fix', 
                              value=value, mask=0)
        method.is_y_vertexFixity(index)
        #********************************************
        def cbk(FPDprolong, newindex):
            return FPDprolong.add_yFixity(kind, location, 
                                    value, weight, newindex,
                                    track_prolongation)
        method.callback = cbk
        method.track_prolongation = track_prolongation
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FixParam[constraint_rank] = method
        return
    
    
    def add_zFixity(self, kind = 'fix', location=None, 
                                 value = None, weight = 1.0, index = None,
                                                     track_prolongation=True):
        """track_prolongation=False if you do not want the index to 
        update during prolongation of a THB lagrangian to a dyadicly 
        refined level
        """
        constraint_rank = len(self.FixParam)
        method = Formparmeter(self.curve, kind='fix', 
                              value=value, mask=0)
        method.is_z_vertexFixity(index)
        self.THBchange[constraint_rank] = method
        #********************************************
        def cbk(FPDprolong, newindex):
            return FPDprolong.add_zFixity(kind, location, 
                                    value, weight, newindex,
                                    track_prolongation)
        method.callback = cbk
        method.track_prolongation = track_prolongation
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FixParam[constraint_rank] = method
        return
    #
    #**************************************************************************
    #
    #"""
    def add_xVertexConstraint(self, kind = 'equality', location=None, 
                              value = None, weight = 1.0, index = None ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_x_vertexConstraint(index)
        #********************************************
        def cbk(FPDprolong,newkind):
            return FPDprolong.add_xVertexConstraint(
                                kind, location, value, weight,index)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
        
    def add_yVertexConstraint(self, kind = 'equality', location=None, 
                              value = None, weight = 1.0, index = None ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_y_vertexConstraint(index)
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_yVertexConstraint(
                                kind, location, value, weight,index)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
        
    def add_zVertexConstraint(self, kind = 'equality', location=None, 
                              value = None, weight = 1.0, index = None ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_z_vertexConstraint(index)
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_zVertexConstraint(
                                kind, location, value, weight,index)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
        
    def add_xPointConstraint(self, 
                             kind = 'equality', 
                             location=None, 
                             value = None, 
                             weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_x_pointConstraint()
        # this is how you make the constraints adapt!
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_xPointConstraint(
                            kind, location, value, weight)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
        
    def add_yPointConstraint(self, 
                             kind = 'equality', 
                             location=None, 
                             value = None, 
                             weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_y_pointConstraint()
        #*************
        def cbk(FPDprolong):
            return FPDprolong.add_yPointConstraint(
                            kind, location, value, weight)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #***********
        self.FormParam[constraint_rank] = method
        return 
        
    def add_zPointConstraint(self, kind = 'equality', 
                    location=None, value = None, weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_z_pointConstraint()
        #***************
        def cbk(FPDprolong):
            return FPDprolong.add_zPointConstraint(
                            kind, location, value, weight)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #*******************
        self.FormParam[constraint_rank] = method
        return 
    """
    # still using 'special' for this...
    def add_functionConstraint(self, kind = 'equality', location=None, fvp = None, weight = 1.0 ):
        """"""fvp :: [function,value] 
                : a function, value pair
        """"""
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value = fvp, weight)
        method.is_functionConstraint()
        method.callback = self.add_functionConstraint
        self.FormParam[constraint_rank] = method
        return 
    #"""
    def add_AngleConstraint(self, kind = 'equality', 
                location=None, value = None, weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_AngleConstraint()
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_AngleConstraint(
                                kind, location, value, weight)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
    
    def add_HarriesTangentConstraint(self, 
                                     kind       = 'equality', 
                                     location   = None, 
                                     value      = None, 
                                     weight     = 1.0, 
                                     index      = None,
                                     rise_index = 1,
                                     run_index = 0):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_HarriesTangentConstraint(index,
                                           rise_index = rise_index,
                                           run_index = run_index)
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_HarriesTangentConstraint(
                        kind, location, value, weight, index)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
        
    def add_verticalAngleConstraint(self, kind = 'equality', 
                    location=None, value = None, weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_Vertical_AngleConstraint()
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_verticalAngleConstraint(
                                kind, location, value, weight)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
    
    #tangent_XoverZ
    def add_AngleXoverZConstraint(self, 
                                  kind = 'equality', 
                                  location=None, 
                                  value = None, 
                                  weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_AngleXoverZConstraint()
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_AngleXoverZConstraint(
                                kind, location, value, weight)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
        
    def add_verticalXoverZAngleConstraint(self, kind = 'equality', 
                                          location=None, value = None, 
                                          weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_VerticalXoverZAngleConstraint()
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_verticalXoverZAngleConstraint(
                                kind, location, value, weight)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
    
    def add_CurvatureConstraint(self, kind = 'equality', location=None, value = None, weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_CurvatureConstraint()
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_CurvatureConstraint(
                                kind, location, value, weight)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
    
    #curvature_XoverZ
    
    def add_CurvatureXoverZConstraint(self, kind = 'equality', 
                                      location=None, value = None, 
                                      weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_CurvatureXoverZConstraint()
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_CurvatureXoverZConstraint(
                                kind, location, value, weight)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
    
    def add_AreaConstraint(self, kind = 'equality', location=None, 
                           value = None, weight = 1.0, axis = 0.):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_AreaConstraint()#y_axis = axis)
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_AreaConstraint(
                                kind, location, value, weight,axis)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
    
    def add_AreaToAnyXAxisConstraint(self, kind = 'equality', location=None, 
                                     value = None, weight = 1.0, x_axis_loc = 0.):
        """ Meant to be used only for 
        area relativeto an -upper- x-axis
        """        
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_Area_to_anyX_Constraint(x_axis_loc)
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_AreaToAnyXAxisConstraint(
                                kind, location, value, weight,x_axis_loc)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
    
    def add_Area_to_y_Constraint(self, kind = 'equality', location=None, 
                                 value = None, weight = 1.0, axis = 0.):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_Area_to_y_Constraint()#x_axis = axis)
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_Area_to_y_Constraint(
                                kind, location, value, weight,axis)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
        
    def add_XcConstraint(self, kind = 'equality', location=None, 
                         value = None, weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_XcConstraint()
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_XcConstraint(
                                kind, location, value, weight)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
        
    def add_YcConstraint(self, kind = 'equality', location=None, 
                         value = None, weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_YcConstraint()
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_YcConstraint(
                                kind, location, value, weight)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return
    
    def add_Xc_to_y_Constraint(self, kind = 'equality', location=None, 
                               value = None, weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_Xc_to_yConstraint()
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_Xc_to_y_Constraint(
                                kind, location, value, weight)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return 
        
    def add_Yc_to_yConstraint(self, kind = 'equality', location=None, 
                              value = None, weight = 1.0 ):
        constraint_rank = len(self.FormParam)
        method = Formparmeter(self.curve, kind, location, value, weight)
        method.is_Yc_to_yConstraint()
        #********************************************
        def cbk(FPDprolong):
            return FPDprolong.add_Yc_to_yConstraint(
                                kind, location, value, weight)
        method.callback = cbk
        self.THBchange[constraint_rank] = method
        #********************************************
        self.FormParam[constraint_rank] = method
        return
        
    
        
    
class Formparmeter(object):
    """
        store class (kind) instance of form parameter
        store which function in the class we will be using
    """
    def __init__(self, 
                 curve, kind='equality', 
                 location=None, value = None, 
                 weight=1.0,mask=1):#function = None):
        self.curve                  = curve
        self.kind                   = kind
        self.computed_value         = None
        self.loc                    = location
        self.Lagrange               = 1.0
        self.slack                  = 1.0
        self.interval_Lagrange      = None
        self.has_contractor         = False #I believe we only contract on real control points
        self.value                  = value
        self.weight                 = weight
        self.mask                   = mask
        
    #
    #**************************************************************************
    #
    def is_x_vertexFixity(self, index):
        self.type                   = 'xfix'
        self.index_of_vertex        = index
        self.computeUpdate          = None
        self.object                 = None
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = 'TBD'#x_point_constraint.contractor
        return
    def is_y_vertexFixity(self, index):
        self.type                   = 'yfix'
        self.index_of_vertex        = index
        self.computeUpdate          = None
        self.object                 = None
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = 'TBD'#x_point_constraint.contractor
        return
    def is_z_vertexFixity(self, index):
        self.type                   = 'zfix'
        self.index_of_vertex        = index
        self.computeUpdate          = None
        self.object                 = None
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = 'TBD'#x_point_constraint.contractor
        return
    #"""
    #
    #**************************************************************************
    #
    def is_relative_x_vertexConstraint(self, index, index2, seperation):
        x_vertex_constraint = relative_x_vertex(self.curve, 
                                                index, index2, seperation)
        self.type                   = 'relative_x_vertex'
        self.computeUpdate          = x_vertex_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = seperation
        self.given                  = seperation
        self.has_contractor         = False
        self.contractor             = 'TBD'#x_point_constraint.contractor
        return
        
    def is_relative_y_vertexConstraint(self, index, index2, seperation):
        y_vertex_constraint = relative_y_vertex(self.curve, 
                                                index, index2, seperation)
        self.type                   = 'relative_y_vertex'
        self.computeUpdate          = y_vertex_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = seperation
        self.given                  = seperation
        self.has_contractor         = False
        self.contractor             = 'TBD'#x_point_constraint.contractor
        return
    #"""
    def is_x_vertexConstraint(self, index):
        x_vertex_constraint = x_vertex(self.curve, index)
        self.type                   = 'x_vertex'
        self.computeUpdate          = x_vertex_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = 'TBD'#x_point_constraint.contractor
        return
        
    def is_y_vertexConstraint(self, index):
        y_vertex_constraint = y_vertex(self.curve, index)
        self.type                   = 'y_vertex'
        self.computeUpdate          = y_vertex_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = 'TBD'#y_point_constraint.contractor
        return
        
    def is_z_vertexConstraint(self, index):
        z_vertex_constraint = z_vertex(self.curve, index)
        self.type                   = 'z_vertex'
        self.computeUpdate          = z_vertex_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = 'TBD'#z_point_constraint.contractor
        return
        
    def is_xAttractor(self, index):
        x_attractor                 = xAttractor(self.curve, index)
        self.type                   = 'xAttract'
        self.computeUpdate          = x_attractor
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = 'NA'#x_attractor.contractor
        return

    def is_yAttractor(self, index):
        y_attractor                 = yAttractor(self.curve, index)
        self.type                   = 'yAttract'
        self.computeUpdate          = y_attractor
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = 'NA'#y_attractor.contractor
        return
        
    def is_x_pointConstraint(self):
        x_point_constraint = x_point(self.curve, self.loc)
        self.type                   = 'x_point'
        self.computeUpdate          = x_point_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = 'TBD'#x_point_constraint.contractor
        return
        
    def is_y_pointConstraint(self):
        y_point_constraint = y_point(self.curve, self.loc)
        self.type                   = 'y_point'
        self.computeUpdate          = y_point_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = 'TBD'#y_point_constraint.contractor
        return
        
    def is_z_pointConstraint(self):
        z_point_constraint = z_point(self.curve, self.loc)
        self.type                   = 'z_point'
        self.computeUpdate          = z_point_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = 'TBD'#z_point_constraint.contractor
        return
    """
    def is_functionConstraint(self):
        gfunc                       = self.value[0](self.curve, self.loc)
        self.type                   = 'gfunc'
        self.computeUpdate          = gfunc
        self.object                 = self.computeUpdate
        self.pass_value             = self.value[1]
        self.given                  = self.value[1]
        self.has_contractor         = False
        self.contractor             = 'TBD'#x_point_constraint.contractor
        return
    #""" 
    def is_AngleConstraint(self):
        tangent_constraint          = tangent(self.curve, self.loc)
        value_radians               = self.value*pi/180.
        self.type                   = 'angle'
        self.computeUpdate          = tangent_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = value_radians
        self.given                  = value_radians
        self.has_contractor         = True
        self.contractor             = tangent_constraint.contractor
        return
    
    def is_HarriesTangentConstraint(self, index,
                                          run_index = 0,
                                          rise_index = 1):
        value_radians               = self.value*pi/180.
        tangent_constraint          = harries_tangent(curve = self.curve, 
                                                      s = self.loc,
                                                      alpha = value_radians,
                                                      index = index,
                                                      run_index = run_index,
                                                      rise_index = rise_index)
        self.type                   = 'angle'
        self.computeUpdate          = tangent_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = value_radians
        self.given                  = value_radians
        self.has_contractor         = False
        self.contractor             = None
        return
        
    def is_Vertical_AngleConstraint(self):
        vertical_tangent_constraint = vertical_tangent(self.curve, self.loc)
        value_radians               = self.value*pi/180.
        self.type                   = 'angle_vertical'
        self.computeUpdate          = vertical_tangent_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = value_radians
        self.given                  = value_radians
        self.has_contractor         = True
        self.contractor             = vertical_tangent_constraint.contractor
        return
    
    #tangent_XoverZ
    def is_AngleXoverZConstraint(self):
        tangent_constraint          = tangent_XoverZ(self.curve, self.loc)
        value_radians               = self.value*pi/180.
        self.type                   = 'angle_XoverZ'
        self.computeUpdate          = tangent_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = value_radians
        self.given                  = value_radians
        self.has_contractor         = False
        self.contractor             = None
        return
    
    def is_VerticalXoverZAngleConstraint(self):
        vertical_tangent_constraint = vertical_XoverZ_tangent(self.curve, self.loc)
        value_radians               = self.value*pi/180.
        self.type                   = 'angle_XoverZ_vertical'
        self.computeUpdate          = vertical_tangent_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = value_radians
        self.given                  = value_radians
        self.has_contractor         = False
        self.contractor             = None#vertical_tangent_constraint.contractor
        return
        
    def is_CurvatureConstraint(self):
        curvature_constraint = curvature(self.curve, self.loc)
        value_radians               = self.value*pi/180.
        self.type                   = 'curvature'
        self.computeUpdate          = curvature_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = value_radians
        self.given                  = value_radians
        self.has_contractor         = True
        self.contractor             = curvature_constraint.contractor
        return
        
    def is_CurvatureXoverZConstraint(self):
        curvature_constraint = curvature_XoverZ(self.curve, self.loc)
        value_radians               = self.value*pi/180.
        self.type                   = 'curvature'
        self.computeUpdate          = curvature_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = value_radians
        self.given                  = value_radians
        self.has_contractor         = False
        self.contractor             = None
        return
    def is_AreaConstraint(self):
        area_constraint = area(self.curve)
        self.type                   = 'area'
        self.computeUpdate          = area_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = True
        self.contractor             = area_constraint.contractor
        #self.y_axis                 = y_axis
        return
        
    def is_Area_to_y_Constraint(self):
        yarea_constraint = area_to_y_axis(self.curve)
        self.type                   = 'area_y_axis_method'
        self.computeUpdate          = yarea_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = True
        self.contractor             = yarea_constraint.contractor
        #self.x_axis                 = x_axis
        return
        
    def is_Area_to_anyX_Constraint(self, x_axis_loc = 0.):
        area_constraint = area_to_any_x_axis(self.curve, x_axis_loc)
        self.type                   = 'area'
        self.computeUpdate          = area_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = True
        self.contractor             = area_constraint.contractor
        self.x_axis                 = x_axis_loc
        return
        
    def is_XcConstraint(self):
        Xc_constraint = Xc(self.curve)
        self.type                   = 'Xc'
        self.computeUpdate          = Xc_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = Xc_constraint.contractor
        return
        
    def is_YcConstraint(self):
        Yc_constraint = Yc(self.curve)
        self.type                   = 'Yc'
        self.computeUpdate          = Yc_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = Yc_constraint.contractor
        return
    
    def is_Xc_to_yConstraint(self):
        Xc_constraint = Xc_Yaxis(self.curve)
        self.type                   = 'Xc'
        self.computeUpdate          = Xc_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = Xc_constraint.contractor
        return
        
    def is_Yc_to_yConstraint(self):
        Yc_constraint = Yc_Yaxis(self.curve)
        self.type                   = 'Yc'
        self.computeUpdate          = Yc_constraint
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = Yc_constraint.contractor
        return
        
    def is_E1(self):
        fairness_function = E1(self.curve)
        self.type                   = 'E1'
        self.computeUpdate          = fairness_function
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = fairness_function.contractor
        return
        
    def is_E2(self):
        fairness_function = E2(self.curve)
        self.type                   = 'E2'
        self.computeUpdate          = fairness_function
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = fairness_function.contractor
        return
    
    def is_E3(self):
        fairness_function = E3(self.curve)
        self.type                   = 'E3'
        self.computeUpdate          = fairness_function
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = fairness_function.contractor
        return
    
    def is_ArcLength(self):
        arc_length_function         = ArcLength(self.curve)
        self.type                   = 'AL'
        self.computeUpdate          = arc_length_function
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = arc_length_function.contractor
        return
        
    def is_ArcLengthApprox(self):
        arc_length_function         = ArcLengthApprox(self.curve)
        self.type                   = 'AL'
        self.computeUpdate          = arc_length_function
        self.object                 = self.computeUpdate
        self.pass_value             = self.value
        self.given                  = self.value
        self.has_contractor         = False
        self.contractor             = arc_length_function.contractor
        return
        
    
        
    

        
        
#def vector_AND_(xi, delta):  #|TLM1> #|000>
#    return_values = copy.deepcopy(xi)
#    for i in range(len(xi)):
#        #element of xi
#        a = xi[i].min.value
#        b = xi[i].real.value
#        c = xi[i].max.value
#        x_simple = AF.fuzzyNumber(a,b,c)
#        
#        #element of delta
#        d = delta[i].min.value
#        e = delta[i].real.value
#        f = delta[i].max.value
#        dz_simple = AF.fuzzyNumber(d,e,f)
#        
#        new_values = x_simple & dz_simple
#        
#        
#        return_values[i].min.value = new_values.min
#        return_values[i].max.value = new_values.max
#        #return_values[i].real.value = new_values.real #krawczyk fails with this
#        return_values[i].real.value = new_values.midpoint()
#        
#    return return_values
    

if __name__ == '__main__':
    pass