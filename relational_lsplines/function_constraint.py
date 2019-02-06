# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 21:33:32 2015

@author: lukemcculloch
"""

import numpy as np
from utilities import vector_AND_



class x_point(object):
    """
        Assumes Area will be used
    """
    def __init__(self, curve, s):
        self.curve = curve
        self.loc    = s
        self.compute_basis()
        self.method_dict = {'equality':self.compute,
                            'max':self.compute_max,
                            'min':self.compute_min,
                            'LS':self.compute_LS}
                            
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : vertices
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
        
    def plot(self, 
             constraint,
             canvas = None, 
             Lspline = None, 
             color = 'green',
             legendloc = None,
             cn = None):
        if canvas == None:print 'function constraint has no canvas'
        
        if Lspline == None:
            print 'Error in function constraint plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        xi,xe,yi,ye = curve.extreme_C0()
        ymin = min(yi,ye)
        
        V = curve.CurvePoint(constraint.loc)
        
        name = r'$h_{} := $'.format(cn) + \
            ' '+ r'$\mathrm{function}$'  + \
            ' - ' +r'$\mathrm{}$'.format(constraint.pass_value) +\
            ' = ' +r'$\mathrm{}$'.format(np.round(V[0]-constraint.pass_value))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (V[0],V[1]+3.),
                            xycoords='data')
        else:
            canvas.annotate(name, 
                            xy = (legendloc[0],legendloc[1]),
                            xycoords='data')
        return canvas
        
    def compute_basis(self):
        localBasis = np.zeros((self.curve.n),float)
        span = self.curve.FindSpan(self.loc)
        self.curve.BasisFuns(span,self.loc,localBasis[span-self.curve.p:span+1])
        self.localBasis = localBasis
        return
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        xpts = vertices[0]
        ypts = vertices[1]
        qx = np.dot(xpts,self.localBasis) 
        return qx
        
    def compute_LS(self, *args, **kwargs):
        function = args[1]
        pt = self.compute(*args, **kwargs)
        return (pt - desired_pt)**2
        
    def compute_max(self, *args, **kwargs):
        function = args[1]
        pt = self.compute(*args, **kwargs)
        return (pt-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        function = args[1]
        pt = self.compute(*args, **kwargs)
        return (pt-minvalue)