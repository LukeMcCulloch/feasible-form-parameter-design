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
             cn = None,
             font_size = 14):
        if canvas == None:print 'x_point plot has no canvas'
        
        if Lspline == None:
            print 'Error in x_point plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            #vertices = [curve.xpts, curve.ypts]
            vertices = curve.ADvertex
        xi,xe,yi,ye = curve.extreme_C0()
        ymin = min(yi,ye)
        
        V = curve.CurvePoint(constraint.loc)
        
        name = r'$h_{%d} := $' % cn+ \
            ' '+ r'$q_{x}$' + \
            ' - ' +r'${}$'.format(constraint.pass_value) +\
            ' = ' +r'${}$'.format(np.round(abs(V[0]-constraint.pass_value)))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = V +[0.,3.5],
                            xycoords='data')
        else:
            """
            canvas.annotate(name, 
                            xy = (legendloc[0],
                                  legendloc[1]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            pstring = r'$q_{x}$' + \
            ' = ' +r'${}$'.format(constraint.pass_value)
            canvas.annotate(pstring, 
                            xy = V +[0.,3.5],
                            xycoords='data',
                            fontsize=font_size)
            #"""
            canvas.plot(V[0],
                        V[1], 
                        marker = 'o', 
                        color = 'black', 
                        alpha = .4,
                        label=name)
       
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
        #ypts = vertices[1]
        qx = np.dot(xpts,self.localBasis) 
        return qx
        
    def compute_LS(self, *args, **kwargs):
        desired_pt = args[1]
        pt = self.compute(*args, **kwargs)
        return (pt - desired_pt)**2
        
    def compute_max(self, *args, **kwargs):
        maxvalue = args[1]
        pt = self.compute(*args, **kwargs)
        return (pt-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        minvalue = args[1]
        pt = self.compute(*args, **kwargs)
        return (pt-minvalue)
        
        

class y_point(object):
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
             cn = None,
             font_size = 14):
        if canvas == None:print 'y_point plot has no canvas'
        
        if Lspline == None:
            print 'Error in y_point plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            #vertices = [curve.xpts,curve.ypts]
            vertices = curve.ADvertices
        else:
            curve = Lspline.curve
            #vertices = [curve.xpts, curve.ypts]
            vertices = curve.ADvertices
        xi,xe,yi,ye = curve.extreme_C0()
        ymin = min(yi,ye)
        
        V = curve.CurvePoint(constraint.loc)
        
        name = r'$h_{%d} := $' % cn  + \
            ' '+ r'$q_{y}$' + \
            ' - ' +r'${}$'.format(constraint.pass_value) +\
            ' = ' +r'${}$'.format(np.round(abs(V[1]-constraint.pass_value)))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = V+[.0,2.5],
                            xycoords='data')
        else:
            """
            canvas.annotate(name, 
                            xy = (legendloc[0],
                                  legendloc[1]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            pstring = r'$q_{y}$' + \
            ' = ' +r'${}$'.format(constraint.pass_value)
            canvas.annotate(pstring, 
                            xy = V +[0.,2.5],
                            xycoords='data',
                            fontsize=font_size)
            #"""
            canvas.plot(V[0],
                        V[1], 
                        marker = 'o', 
                        color = 'black', 
                        alpha = .4,
                        label=name)
        return canvas
        
    def compute_basis(self):
        localBasis = np.zeros((self.curve.n),float)
        span = self.curve.FindSpan(self.loc)
        self.curve.BasisFuns(span,self.loc,localBasis[span-self.curve.p:span+1])
        self.localBasis = localBasis
        return
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        #xpts = vertices[0]
        ypts = vertices[1]
        qy = np.dot(ypts,self.localBasis) 
        return qy
        
    def compute_LS(self, *args, **kwargs):
        desired_pt = args[1]
        pt = self.compute(*args, **kwargs)
        return (pt - desired_pt)**2
        
    def compute_max(self, *args, **kwargs):
        maxvalue = args[1]
        pt = self.compute(*args, **kwargs)
        return (pt-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        minvalue = args[1]
        pt = self.compute(*args, **kwargs)
        return (pt-minvalue)
    
class z_point(object):
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
             cn = None,
             font_size = 14):
        if canvas == None:print 'y_point plot has no canvas'
        
        if Lspline == None:
            print 'Error in y_point plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = curve.ADvertex
        else:
            curve = Lspline.curve
            vertices = curve.ADvertex
        xi,xe,yi,ye = curve.extreme_C0()
        ymin = min(yi,ye)
        
        V = curve.CurvePoint(constraint.loc)
        
        name = r'$h_{%d} := $' % cn  + \
            ' '+ r'$q_{y}$' + \
            ' - ' +r'${}$'.format(constraint.pass_value) +\
            ' = ' +r'${}$'.format(np.round(abs(V[1]-constraint.pass_value)))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = V+[.0,2.5],
                            xycoords='data')
        else:
            """
            canvas.annotate(name, 
                            xy = (legendloc[0],
                                  legendloc[1]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            pstring = r'$q_{y}$' + \
            ' = ' +r'${}$'.format(constraint.pass_value)
            canvas.annotate(pstring, 
                            xy = V +[0.,2.5],
                            xycoords='data',
                            fontsize=font_size)
            #"""
            canvas.plot(V[0],
                        V[1], 
                        marker = 'o', 
                        color = 'black', 
                        alpha = .4,
                        label=name)
        return canvas
        
    def compute_basis(self):
        localBasis = np.zeros((self.curve.n),float)
        span = self.curve.FindSpan(self.loc)
        self.curve.BasisFuns(span,self.loc,localBasis[span-self.curve.p:span+1])
        self.localBasis = localBasis
        return
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        #xpts = vertices[0]
        zpts = vertices[2]
        qz = np.dot(zpts,self.localBasis) 
        return qz
        
    def compute_LS(self, *args, **kwargs):
        desired_pt = args[1]
        pt = self.compute(*args, **kwargs)
        return (pt - desired_pt)**2
        
    def compute_max(self, *args, **kwargs):
        maxvalue = args[1]
        pt = self.compute(*args, **kwargs)
        return (pt-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        minvalue = args[1]
        pt = self.compute(*args, **kwargs)
        return (pt-minvalue)