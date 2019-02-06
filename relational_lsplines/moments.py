import numpy as np
import copy
from utilities import vector_AND_



class Xc(object):
    """X-centroid constraint
    
        memoize the curve's integrated basis product matrix
        used for the x-centroid computation
    
    Assumptions:
    --------------
        Area will be used
    """
    def __init__(self, curve):
        self.curve = curve
        self.curve.MomentMatrices()
        self.AM2 = curve.AM2
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
        if canvas == None:print 'Xc plot has no canvas'
        
        if Lspline == None:
            print 'Error in Xc plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        xi,xe,yi,ye = curve.extreme_C0()
        ymin = min(yi,ye)
        
        curve.compute_moments()
        Xc = curve.My/curve.area.value
        Yc = curve.Mx/curve.area.value
        name = r'$h_{} := $'.format(str(cn)) + \
            ' '+ r'$x_{c}$'  + \
            ' - ' +r'${}$'.format(constraint.pass_value) +\
            ' = ' +r'${}$'.format(np.round(abs(curve.Xc.value-constraint.pass_value),decimals=2))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (Xc-1.5,Yc- 2.5),
                            xycoords='data')
        else:
            if Xc <=.5*legendloc[0]:
                xloc = legendloc[2]
            else:
                xloc = legendloc[0]
            """
            canvas.annotate(name, 
                            xy = (xloc,legendloc[1]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            mstring = r'$x_{c}$'  + \
            ' = ' +r'${}$'.format(constraint.pass_value)
            canvas.annotate(mstring, 
                            xy = (Xc+.5,Yc-.5),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            canvas.plot(Xc,Yc, 
                        marker = 'o', 
                        color = 'black', 
                        alpha = .4,
                        label=name)
        return canvas
        
    def compute(self, *args, **kwargs):
        vertices    = args[0]
        xpts        = vertices[0]
        ypts        = vertices[1]
        n           = self.curve.n

        temp_array0=[]
        temp_array1=[]
        for l in range(n):
            temp_array0.append(list())
            temp_array1.append(list())
            for i in range(n):
                temp_array0[l].append( np.dot(xpts,np.transpose(self.curve.AM2[l,i,:])) )
                temp_array1[l].append( np.dot(ypts,np.transpose(self.curve.AM2[l,i,:])) )
    
        temp_array0=np.asarray(temp_array0)
        temp_array1=np.asarray(temp_array1)
                
        temp5  = np.dot(xpts,np.dot(ypts,temp_array0))
        temp6  = np.dot(xpts,np.dot(xpts,temp_array1))
        store4 = (ypts[-1]*xpts[-1]*xpts[-1]) - (ypts[0]*xpts[0]*xpts[0])
        self.curve.My   = (store4 + (temp5-temp6))*(1./3.) #sign swapped to '+' <midnight June14 - this was essential I believe.
        
#        try:
#            self.curve.Xc  = self.curve.My/self.curve.area
#        except:
        dotprdct   = np.dot(ypts,np.dot(self.curve.AM1,xpts))
        smple_bdry  = (xpts[-1]*ypts[-1] - xpts[0]*ypts[0])
        self.curve.area = (dotprdct + smple_bdry)*.5
        self.curve.Xc  = self.curve.My/self.curve.area
        return self.curve.Xc
        
    def compute_LS(self, *args, **kwargs):
        desired_Xc = args[1]
        #        if len(args)>2:
        #            print 'activating decorator'
        #            decorator = args[2]
        #            Xc = decorator(self.compute, *args, **kwargs)
        #        else:
        Xc = self.compute(*args, **kwargs)
        return ((Xc-desired_Xc)**2).sqrt()
        
    def compute_max(self, *args, **kwargs):
        max_Xc = args[1]
        Xc = self.compute(*args, **kwargs)
        return (Xc-max_Xc)*(-1.)
        
    def compute_min(self, *args, **kwargs):
        min_Xc = args[1]
        Xc = self.compute(*args, **kwargs)
        return (Xc - min_Xc)

    def contractor(self, *args, **kwargs):
        vertices    = copy.deepcopy(args[0])
        nrange = len(vertices[0])
        xpts = []
        ypts = []
        for i in range(nrange):
            xpts.append(vertices[0][i].value)
            ypts.append(vertices[1][i].value)
        constraint  = copy.deepcopy(args[1])
        
        
        print 'Xc contractor Not Implemented Yet'
        for i in range(nrange):
            vertices[0][i].value = xpts[i]
            vertices[1][i].value = ypts[i]
        return vertices



class Yc(object):
    """Y-centroid constraint
    
        *memoize the curve's integrated basis product matrix
        used for the x-centroid computation
        
    Assumptions:
    --------------
        *Area will be used
    """
    def __init__(self, curve):
        self.curve = curve
        self.curve.MomentMatrices()
        self.AM2 = curve.AM2
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
        if canvas == None:print 'Yc plot has no canvas'
        
        if Lspline == None:
            print 'Error in Yc plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        xi,xe,yi,ye = curve.extreme_C0()
        ymin = min(yi,ye)
        
        curve.compute_moments()
        Xc = curve.My/curve.area.value
        Yc = curve.Mx/curve.area.value
        name = r'$h_{} := $'.format(str(cn)) + \
            ' '+ r'$y_{c}$'  + \
            ' - ' +r'${}$'.format(abs(constraint.pass_value)) +\
            ' = ' +r'${}$'.format(np.round(abs(curve.Yc.value-constraint.pass_value),decimals=2))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (Xc-1.5,Yc- 2.5),
                            xycoords='data')
        else:
            if False:
                xloc = legendloc[2]
            else:
                xloc = legendloc[0]
            """
            canvas.annotate(name, 
                            xy = (xloc,legendloc[1]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            mstring = r'$y_{c}$'  + \
            ' = ' +r'${}$'.format(abs(constraint.pass_value))
            canvas.annotate(mstring, 
                            xy = (Xc+.5,Yc-.5),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            canvas.plot(Xc,Yc, 
                        marker = 'o', 
                        color = 'black', 
                        alpha = .4,
                        label=name)
        return canvas
        
    def compute(self, *args, **kwargs):
        vertices    = args[0]
        xpts        = vertices[0]
        ypts        = vertices[1]
        n           = self.curve.n

        temp_array0=[]
        temp_array1=[]
        for l in range(n):
            temp_array0.append(list())
            temp_array1.append(list())
            for i in range(n):
                temp_array0[l].append( np.dot(xpts,np.transpose(self.AM2[l,i,:])) )
                temp_array1[l].append( np.dot(ypts,np.transpose(self.AM2[l,i,:])) )
    
        temp_array0=np.asarray(temp_array0)
        temp_array1=np.asarray(temp_array1)
                
        temp3  = np.dot(ypts,np.dot(ypts,temp_array0))
        temp4  = np.dot(ypts,np.dot(xpts,temp_array1))
        store3 = (xpts[-1]*ypts[-1]*ypts[-1])-(xpts[0]*ypts[0]*ypts[0])
        self.curve.Mx  = (store3+((temp3-temp4)*2.))*(1./6)
#        try:
#            self.curve.Yc  = self.curve.Mx/self.curve.area
#        except:
        dotprdct   = np.dot(ypts,np.dot(self.curve.AM1,xpts))
        smple_bdry  = (xpts[-1]*ypts[-1] - xpts[0]*ypts[0])
        self.curve.area = (dotprdct + smple_bdry)*.5
        self.curve.Yc  = self.curve.Mx/self.curve.area
        return self.curve.Yc
        
    def compute_LS(self, *args, **kwargs):
        desired_Yc = args[1]
        Yc = self.compute(*args, **kwargs)
        return ((Yc-desired_Yc)**2).sqrt()
        
    def compute_max(self, *args, **kwargs):
        max_Yc = args[1]
        Yc = self.compute(*args, **kwargs)
        return (Yc-max_Yc)*(-1.)
        
    def compute_min(self, *args, **kwargs):
        min_Yc = args[1]
        Yc = self.compute(*args, **kwargs)
        return (Yc - min_Yc)

    def contractor(self, *args, **kwargs):
        vertices    = copy.deepcopy(args[0])
        nrange = len(vertices[0])
        xpts = []
        ypts = []
        for i in range(nrange):
            xpts.append(vertices[0][i].value)
            ypts.append(vertices[1][i].value)
        constraint  = copy.deepcopy(args[1])
        
        
        print 'Yc contractor Not Implemented Yet'
        for i in range(nrange):
            vertices[0][i].value = xpts[i]
            vertices[1][i].value = ypts[i]
        return vertices
        
        
##-----------------------------------------------------------------------------
## -------------------- Moments of Area with the Y Axis -----------------------
##-----------------------------------------------------------------------------
class Xc_Yaxis(object):
    """
        Xcenter of the Area from the
        curve to the Y axis
    """
    def __init__(self, curve):
        self.curve = curve
        self.curve.MomentMatrices()
        self.AM2 = curve.AM2
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
        if canvas == None:print 'Xc plot has no canvas'
        
        if Lspline == None:
            print 'Error in Xc plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        xi,xe,yi,ye = curve.extreme_C0()
        ymin = min(yi,ye)
        
        curve.compute_moments_to_y()
        Xc = curve.My_to_y/curve.area_to_y.value
        Yc = curve.Mx_to_y/curve.area_to_y.value
        name = r'$h_{} := $'.format(cn) + \
            ' '+ r'$x_{c}$'  + \
            ' - ' +r'${}$'.format(constraint.pass_value) +\
            ' = ' +r'${}$'.format(np.round(abs(curve.Xc_to_y.value-constraint.pass_value),decimals=2))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (Xc-1.5,Yc- 2.5),
                            xycoords='data')
        else:
            """
            canvas.annotate(name, 
                            xy = (legendloc[0],
                                  legendloc[1]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            canvas.annotate(str(cn), 
                            xy = (Xc,Yc-.5),
                            xycoords='data',
                            fontsize=font_size)
            #tacked on Janurary 31 2016:
            canvas.plot(Xc,Yc, 
                        marker = 'o', 
                        color = 'black', 
                        alpha = .4,
                        label=name)
        return canvas
        
    def compute(self, *args, **kwargs):
        vertices    = args[0]
        xpts        = vertices[0]
        ypts        = vertices[1]
        n           = self.curve.n

        temp_array0=[]
        temp_array1=[]
        for l in range(n):
            temp_array0.append(list())
            temp_array1.append(list())
            for i in range(n):
                temp_array0[l].append( np.dot(xpts,np.transpose(self.curve.AM2[l,i,:])) )
                temp_array1[l].append( np.dot(ypts,np.transpose(self.curve.AM2[l,i,:])) )
    
        temp_array0=np.asarray(temp_array0)
        temp_array1=np.asarray(temp_array1)
                
        temp5  = np.dot(xpts,np.dot(ypts,temp_array0))
        temp6  = np.dot(xpts,np.dot(xpts,temp_array1))
        store4 = -( (ypts[-1]*xpts[-1]*xpts[-1]) - (ypts[0]*xpts[0]*xpts[0]) )
        self.curve.My_to_y   = -(store4 + (temp5-temp6))*(1./3.) 
        
#        try:
#            self.curve.Xc  = self.curve.My/self.curve.area
#        except:
        dotprdct   = np.dot(ypts,np.dot(self.curve.AM1,xpts))
        smple_bdry  = -(xpts[-1]*ypts[-1] - xpts[0]*ypts[0])
        self.curve.area_to_y = -(dotprdct + smple_bdry)*.5
        self.curve.Xc_to_y  = self.curve.My_to_y/self.curve.area_to_y
        return self.curve.Xc_to_y
        
    def compute_LS(self, *args, **kwargs):
        desired_Xc = args[1]
        #        if len(args)>2:
        #            print 'activating decorator'
        #            decorator = args[2]
        #            Xc = decorator(self.compute, *args, **kwargs)
        #        else:
        Xc = self.compute(*args, **kwargs)
        return ((Xc-desired_Xc)**2).sqrt()
        
    def compute_max(self, *args, **kwargs):
        max_Xc = args[1]
        Xc = self.compute(*args, **kwargs)
        return (Xc-max_Xc)*(-1.)
        
    def compute_min(self, *args, **kwargs):
        min_Xc = args[1]
        Xc = self.compute(*args, **kwargs)
        return (Xc - min_Xc)

    def contractor(self, *args, **kwargs):
        vertices    = copy.deepcopy(args[0])
        nrange = len(vertices[0])
        xpts = []
        ypts = []
        for i in range(nrange):
            xpts.append(vertices[0][i].value)
            ypts.append(vertices[1][i].value)
        constraint  = copy.deepcopy(args[1])
        
        
        print 'Xc contractor Not Implemented Yet'
        for i in range(nrange):
            vertices[0][i].value = xpts[i]
            vertices[1][i].value = ypts[i]
        return vertices


class Yc_Yaxis(object):
    """
        Assumes Area will be used
    """
    def __init__(self, curve):
        self.curve = curve
        self.curve.MomentMatrices()
        self.AM2 = curve.AM2
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
        if canvas == None:print 'Yc plot has no canvas'
        
        if Lspline == None:
            print 'Error in Yc plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        xi,xe,yi,ye = curve.extreme_C0()
        ymin = min(yi,ye)
        
        curve.compute_moments_to_y()
        Xc = curve.My_to_y/curve.area_to_y.value
        Yc = curve.Mx_to_y/curve.area_to_y.value
        name = r'$h_{} := $'.format(cn) + \
            ' '+ r'$y_{c}$'  + \
            ' - ' +r'${}$'.format(abs(constraint.pass_value)) +\
            ' = ' +r'${}$'.format(np.round(abs(curve.Yc_to_y.value-constraint.pass_value),decimals=2))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (Xc-1.5,Yc- 2.5),
                            xycoords='data')
        else:
            if Xc <=.5*legendloc[0]:
                xloc = legendloc[2]
            else:
                xloc = legendloc[0]
            """
            canvas.annotate(name, 
                            xy = (xloc,legendloc[1]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            mstring = r'$y_{c}$'  + \
            ' = ' +r'${}$'.format(abs(constraint.pass_value)) 
            canvas.annotate(mstring, 
                            xy = (Xc,Yc-.5),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            canvas.plot(Xc,Yc, 
                        marker = 'o', 
                        color = 'black', 
                        alpha = .4,
                        label=name)
        return canvas
        
    def compute(self, *args, **kwargs):
        vertices    = args[0]
        xpts        = vertices[0]
        ypts        = vertices[1]
        n           = self.curve.n

        temp_array0=[]
        temp_array1=[]
        for l in range(n):
            temp_array0.append(list())
            temp_array1.append(list())
            for i in range(n):
                temp_array0[l].append( np.dot(xpts,np.transpose(self.AM2[l,i,:])) )
                temp_array1[l].append( np.dot(ypts,np.transpose(self.AM2[l,i,:])) )
    
        temp_array0=np.asarray(temp_array0)
        temp_array1=np.asarray(temp_array1)
                
        temp3  = np.dot(ypts,np.dot(ypts,temp_array0))
        temp4  = np.dot(ypts,np.dot(xpts,temp_array1))
        store3 = -( (xpts[-1]*ypts[-1]*ypts[-1])-(xpts[0]*ypts[0]*ypts[0]) )
        self.curve.Mx_to_y  = -(store3+((temp3-temp4)*2.))*(1./6)
#        try:
#            self.curve.Yc  = self.curve.Mx/self.curve.area
#        except:
        dotprdct   = np.dot(ypts,np.dot(self.curve.AM1,xpts))
        smple_bdry  = -(xpts[-1]*ypts[-1] - xpts[0]*ypts[0])
        self.curve.area_to_y = -(dotprdct + smple_bdry)*.5
        self.curve.Yc_to_y  = self.curve.Mx_to_y/self.curve.area_to_y
        return self.curve.Yc_to_y
        
    def compute_LS(self, *args, **kwargs):
        desired_Yc = args[1]
        Yc = self.compute(*args, **kwargs)
        return ((Yc-desired_Yc)**2).sqrt()
        
    def compute_max(self, *args, **kwargs):
        max_Yc = args[1]
        Yc = self.compute(*args, **kwargs)
        return (Yc-max_Yc)*(-1.)
        
    def compute_min(self, *args, **kwargs):
        min_Yc = args[1]
        Yc = self.compute(*args, **kwargs)
        return (Yc - min_Yc)

    def contractor(self, *args, **kwargs):
        vertices    = copy.deepcopy(args[0])
        nrange = len(vertices[0])
        xpts = []
        ypts = []
        for i in range(nrange):
            xpts.append(vertices[0][i].value)
            ypts.append(vertices[1][i].value)
        constraint  = copy.deepcopy(args[1])
        
        
        print 'Yc contractor Not Implemented Yet'
        for i in range(nrange):
            vertices[0][i].value = xpts[i]
            vertices[1][i].value = ypts[i]
        return vertices