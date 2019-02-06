"""
    Generalized attractor class
    instantiate with object to attract    
"""
import numpy as np
from utilities import vector_AND_




class x_vertex(object):
    def __init__(self, curve, index):
        self.curve = curve
        self.index = index
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
             cn = None,
             font_size = 14):
        if canvas == None:print 'x_vertex plot has no canvas'
        
        if Lspline == None:
            print 'Error in x_vertex plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        xi,xe,yi,ye = curve.extreme_C0()
        ymin = min(yi,ye)
        
        V = []
        V.append(curve.xpts[self.index].value)#curve.CurvePoint(constraint.loc)
        V.append(curve.ypts[self.index].value)
        
        name = r'$h_{%d} := $' % cn + \
            ' '+ r'${x vertex}$'  + \
            ' - ' +r'${}$'.format(constraint.pass_value) +\
            ' = ' +r'${}$'.format(np.round(V[0]-constraint.pass_value))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (V[0],V[1]+3.),
                            xycoords='data')
        else:
            """
            canvas.annotate(name, 
                            xy = (legendloc[0],legendloc[1]),
                            xycoords='data')
            #"""
        canvas.plot(V[0],
                    V[1], 
                    marker = 'o', 
                    color = 'black', 
                    alpha = .4,
                    label=name)
        return canvas
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        xpts = vertices[0]
        return xpts[self.index]
        
    def compute_LS(self, *args, **kwargs):
        desired_loc = args[1]
        xpt = self.compute(*args, **kwargs)
        return (xpt - desired_loc)**2
        
    def compute_max(self, *args, **kwargs):
        maxvalue = args[1]
        value = self.compute(*args, **kwargs)
        return (value-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        minvalue = args[1]
        value = self.compute(*args, **kwargs)
        return (value-minvalue)
    
    
class y_vertex(object):
    def __init__(self, curve, index):
        self.curve = curve
        self.index = index
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
             cn = None,
             font_size = 14):
        if canvas == None:print 'y_vertex plot has no canvas'
        
        if Lspline == None:
            print 'Error in y_vertex plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        xi,xe,yi,ye = curve.extreme_C0()
        ymin = min(yi,ye)
        
        V = []
        V.append(curve.xpts[self.index].value)#curve.CurvePoint(constraint.loc)
        V.append(curve.ypts[self.index].value)
        
        name = r'$h_{%d} := $' % cn + \
            ' '+ r'${y vertex}$'  + \
            ' - ' +r'${}$'.format(constraint.pass_value) +\
            ' = ' +r'${}$'.format(np.round(V[1]-constraint.pass_value))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (V[0],V[1]+3.),
                            xycoords='data')
        else:
            """
            canvas.annotate(name, 
                            xy = (legendloc[0],legendloc[1]),
                            xycoords='data')
            #"""
        canvas.plot(V[0],
                    V[1], 
                    marker = 'o', 
                    color = 'black', 
                    alpha = .4,
                    label=name)
        return canvas
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        ypts = vertices[1]
        return ypts[self.index]
        
    def compute_LS(self, *args, **kwargs):
        desired_loc = args[1]
        ypt = self.compute(*args, **kwargs)
        return (ypt - desired_loc)**2
        
    def compute_max(self, *args, **kwargs):
        maxvalue = args[1]
        value = self.compute(*args, **kwargs)
        return (value-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        minvalue = args[1]
        value = self.compute(*args, **kwargs)
        return (value-minvalue)
    
     
        
    
    
class z_vertex(object):
    def __init__(self, curve, index):
        self.curve = curve
        self.index = index
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
             cn = None,
             font_size = 14):
        if canvas == None:print 'y_vertex plot has no canvas'
        
        if Lspline == None:
            print 'Error in y_vertex plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts, curve.zpts]
        xi,xe,yi,ye = curve.extreme_C0()
        ymin = min(yi,ye)
        
        V = []
        V.append(curve.xpts[self.index].value)#curve.CurvePoint(constraint.loc)
        V.append(curve.ypts[self.index].value)
        V.append(curve.zpts[self.index].value)
        
        name = r'$h_{%d} := $' % cn + \
            ' '+ r'${y vertex}$'  + \
            ' - ' +r'${}$'.format(constraint.pass_value) +\
            ' = ' +r'${}$'.format(np.round(V[1]-constraint.pass_value))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (V[0],V[1]+3.),
                            xycoords='data')
        else:
            """
            canvas.annotate(name, 
                            xy = (legendloc[0],legendloc[1]),
                            xycoords='data')
            #"""
        canvas.plot(V[0],
                    V[1], 
                    marker = 'o', 
                    color = 'black', 
                    alpha = .4,
                    label=name)
        return canvas
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        zpts = vertices[2]
        return zpts[self.index]
        
    def compute_LS(self, *args, **kwargs):
        desired_loc = args[1]
        zpt = self.compute(*args, **kwargs)
        return (zpt - desired_loc)**2
        
    def compute_max(self, *args, **kwargs):
        maxvalue = args[1]
        value = self.compute(*args, **kwargs)
        return (value-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        minvalue = args[1]
        value = self.compute(*args, **kwargs)
        return (value-minvalue)
    
    


class relative_x_vertex(object):
    def __init__(self, curve, index, index2, seperation):
        self.curve = curve
        self.index = index
        self.index2 = index2
        self.seperation = seperation
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
             cn = None,
             font_size = 14):
        if canvas == None:print 'x_vertex plot has no canvas'
        
        if Lspline == None:
            print 'Error in x_vertex plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        xi,xe,yi,ye = curve.extreme_C0()
        ymin = min(yi,ye)
        
        V = []
        V.append(curve.xpts[self.index].value)#curve.CurvePoint(constraint.loc)
        V.append(curve.ypts[self.index].value)
        
        name = r'$h_{%d} := $' % cn + \
            ' '+ r'${x vertex}$'  + \
            ' - ' +r'${}$'.format(constraint.pass_value) +\
            ' = ' +r'${}$'.format(np.round(V[0]-constraint.pass_value))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (V[0],V[1]+3.),
                            xycoords='data')
        else:
            canvas.annotate(name, 
                            xy = (legendloc[0],legendloc[1]),
                            xycoords='data')
        return canvas
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        xpts = vertices[0]
        return xpts[self.index]-xpts[self.index2]
        
    def compute_LS(self, *args, **kwargs):
        desired_loc = args[1]
        xpt = self.compute(*args, **kwargs)
        return (xpt - desired_loc)**2
        
    def compute_max(self, *args, **kwargs):
        maxvalue = args[1]
        value = self.compute(*args, **kwargs)
        return (value-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        minvalue = args[1]
        value = self.compute(*args, **kwargs)
        return (value-minvalue)*(-1.0)
    
    
class relative_y_vertex(object):
    def __init__(self, curve, index, index2, seperation):
        self.curve = curve
        self.index = index
        self.index2 = index2
        self.seperation = seperation
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
             cn = None,
             font_size = 14):
        if canvas == None:print 'y_vertex plot has no canvas'
        
        if Lspline == None:
            print 'Error in y_vertex plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        xi,xe,yi,ye = curve.extreme_C0()
        ymin = min(yi,ye)
        
        V = []
        V.append(curve.xpts[self.index].value)#curve.CurvePoint(constraint.loc)
        V.append(curve.ypts[self.index].value)
        
        name = r'$h_{%d} := $' % cn + \
            ' '+ r'${y vertex}$'  + \
            ' - ' +r'${}$'.format(constraint.pass_value) +\
            ' = ' +r'${}$'.format(np.round(V[1]-constraint.pass_value))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (V[0],V[1]+3.),
                            xycoords='data')
        else:
            canvas.annotate(name, 
                            xy = (legendloc[0],legendloc[1]),
                            xycoords='data')
        return canvas
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        ypts = vertices[1]
        return ypts[self.index]-ypts[self.index2]
        
    def compute_LS(self, *args, **kwargs):
        desired_loc = args[1]
        ypt = self.compute(*args, **kwargs)
        return (ypt - desired_loc)**2
        
    def compute_max(self, *args, **kwargs):
        maxvalue = args[1]
        value = self.compute(*args, **kwargs)
        return (value-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        minvalue = args[1]
        value = self.compute(*args, **kwargs)
        return (value-minvalue)*(-1.0)
    