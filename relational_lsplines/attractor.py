"""
    Generalized attractor class
    instantiate with object to attract    
"""
import numpy as np
from utilities import vector_AND_



class xAttractor(object):
    def __init__(self, curve, index):
        self.curve = curve
        self.index = index
        self.method_dict = {'equality':self.compute,
                            'LS':self.compute_LS}
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : vertices
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        xpts = vertices[0]
        return xpts[self.index]
        
    def compute_LS(self, *args, **kwargs):
        xpt = self.compute(*args, **kwargs)
        return (xpt-1000.)**2
        

    
class yAttractor(object):
    def __init__(self, curve, index):
        self.curve = curve
        self.index = index
        self.method_dict = {'equality':self.compute,
                            'LS':self.compute_LS}
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : vertices
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        ypts = vertices[1]
        return ypts[self.index]
        
    def compute_LS(self, *args, **kwargs):
        ypt = self.compute(*args, **kwargs)
        return (ypt-1000.)**2
        
